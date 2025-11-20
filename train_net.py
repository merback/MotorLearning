# SPDX-License-Identifier: GPL-3.0

import numpy as np
import torch
import torch.nn as nn
import pandas as pd
import sys
from time import time

DEVICE = "cpu"

# Function to read all CSV files and create a data matrix
def create_data_matrix(folder_path, file_name, num_samples):
    data = []
    if IS_AIRGAP:
        n=462
    else:
        n=6354
    for i in range(num_samples):
        d = pd.read_csv(folder_path + '/' + file_name + str(i+1) + '.csv', header=None) 
        x = np.reshape(d.values, n)
        data.append(x)
    return np.array(data).T


def normalize_data(data, ranges):

    normalized = np.zeros_like(data)
    for i, (min_val, max_val) in enumerate(ranges):
        normalized[i, :] = (data[i, :] - min_val) / (max_val - min_val)
    return normalized

def perform_pod(data_matrix, K, num_modes):
     # 1. Weighted correlation matrix in snapshot space
    C = data_matrix.T @ K @ data_matrix  # (n_snapshots, n_snapshots)

    # 2. Solve eigenproblem
    eigvals, eigvecs = np.linalg.eig(C)

    # 3. Sort eigenvalues descending
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]

    # 4. Project back to coefficient space
    modes = data_matrix @ eigvecs  # (n_coeffs, n_snapshots)

    # 5. Normalize each mode
    for i in range(num_modes):
        phi = modes[:, i]
        norm = np.sqrt(eigvals[i])
        if norm > 1e-12:
            modes[:, i] = phi / norm
        else:
            modes[:, i] = 0.0

    return modes[:,:num_modes], eigvals[:num_modes]

# Project data onto the top k POD modes
def project_onto_pod_modes(data_matrix, U, K):
    return U.T @ K @ data_matrix

def loss_function_pod(y_pred, y_true_red, K_red):
    # y_pred: predicted data in pod space
    # y_true: true data in full space
    y_pred = y_pred * feature_stds_t + feature_means_t

    diff = y_pred - y_true_red
    numerator = torch.sum((diff@K_red) * diff, dim=1)
    denominator = torch.sum((y_true_red@K_red) * y_true_red, dim=1)
    errors = numerator / (denominator + 1e-12)
    return torch.mean(errors)  

def calc_error(model,x,Y,K):
    model.eval()
    with torch.no_grad():
        y_test_norm = model(x)
        y_e = (y_test_norm * feature_stds_t) + feature_means_t
        y_e =  U @ y_e.T 
        errors = relative_error(Y,y_e,K)
        errors = np.array(errors)
        error_max = np.max(errors)
        error_mean = np.mean(errors)
        error_std = np.std(errors)
    model.train()
    return error_max, error_mean, error_std, errors

def analyze_errors(model, epoch):
    # Calculate training errors
    max_err_train, mean_err_train, std_err_train, train_errors = calc_error(model,x_train,Y_train,K)

    # Calculate validation errors
    max_err_val, mean_err_val, std_err_val, val_errors = calc_error(model,x_val,Y_val,K)

    print(f"\n=== Error Analysis - Epoch {epoch} ===")
    print(f"Training Error - Max: {max_err_train:.6f}, Mean: {mean_err_train:.6f}, Std: {std_err_train:.6f}")
    print(f"Validation Error - Max: {max_err_val:.6f}, Mean: {mean_err_val:.6f}, Std: {std_err_val:.6f}")

    return max_err_train, mean_err_train, std_err_train, max_err_val, mean_err_val, std_err_val

def relative_error(original, reconstructed, K):
     # original, reconstructed: (n_features, n_samples)  
    # K: (n_features, n_features)
    
    diff = reconstructed - original  # (n_features, n_samples)
    
    # Compute numerator for all samples at once: (x-y)^T @ K @ (x-y) for each sample
    diff_K = diff.T @ K  # (n_samples, n_features)
    numerator = torch.sum(diff_K * diff.T, dim=1)  # (n_samples,)
    
    # Compute denominator for all samples at once: y^T @ K @ y for each sample
    original_K = original.T @ K  # (n_samples, n_features)  
    denominator = torch.sum(original_K * original.T, dim=1)  # (n_samples,)
    
    # Compute relative error for each sample (square root of relative error squared)
    errors = torch.sqrt(numerator / denominator)  # (n_samples,)
    
    # Convert to list
    return [errors[i] for i in range(errors.shape[0])]

def train(model, optimizer):
    model.train()
    optimizer.zero_grad()
    outputs = model(x_train)
    outputs = outputs.to(DEVICE)
    loss = loss_function_pod(outputs, y_train_red, K_red)
    loss.backward()
    optimizer.step()
    epoch_loss = loss.item()

    return epoch_loss

def define_model(features, input_dim, modes):
    layers = []
    in_features = input_dim
    for h in features:
        layers.append(nn.Linear(in_features, h))
        layers.append(nn.ReLU())
        in_features = h
    layers.append(nn.Linear(in_features, modes))
    return nn.Sequential(*layers)

def run(in_feature, modes):

    hidden_dim = [336,336,336] if not IS_AIRGAP else [190,110,180]
    n_layers = 3
    lr = 0.00634457140869262 if not IS_AIRGAP else 1e-2
    wdecay = 1.23e-6 if not IS_AIRGAP else 1e-6
    gamma = 0.973 if not IS_AIRGAP else 0.9776622845249979
    step_size = 200
    epochs = 11000
    
    model = define_model(hidden_dim, in_feature,modes);# [170, 170, 210]); #5, [150, 80, 100, 180, 130]);# [200, 200, 200]) #[140,170,180]# Define the model with 6 layers and 100 units each
    model = model.to(DEVICE).double()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wdecay)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=step_size, gamma=gamma)

    for epoch in range(epochs):
        epoch_loss = train(model, optimizer)

        if epoch_loss < EARLY_STOPPING_THRESHOLD:
            print(f'\n Early stopping at epoch {epoch}: Loss {epoch_loss:.4e} below threshold {EARLY_STOPPING_THRESHOLD:.4e}')
            break
        scheduler.step()
       # if (epoch + 1) % 10 == 0:
       #     print(f'Epoch {epoch+1}/{num_epochs}, Loss: {epoch_loss:.8f}')
        if (epoch + 1) % 1000 == 0:
            analyze_errors(model, epoch + 1)


    return model


def run_training():
       
    start_time = time()
    torch.manual_seed(100)

    global U, K, K_red, x_train, y_train_red, Y_train, x_val, Y_val, feature_means_t, feature_stds_t
    
    if IS_AIRGAP:
        domain = 'gap'
        modes =14
        in_feature=4

    else:
        domain = 'full'
        modes = 90 
        in_feature=14 

    print(f'Starting training process for {domain} model')

    # Define paths and file names
    train_folder_path = f'sol_{domain}_training'
    train_file_name = f'sol_{domain}_training_'
    val_folder_path = f'sol_{domain}_validation'
    val_file_name = f'sol_{domain}_validation_'

    samples = 1024 # Number of samples
    samples_val = 128  # Number of validation samples


    # Read data
    data_train = create_data_matrix(train_folder_path, train_file_name, samples)
    data_val =create_data_matrix(val_folder_path, val_file_name, samples_val)

    input_train = np.array(pd.read_csv(f'samples_{domain}_training.csv').values)
    input_val = np.array(pd.read_csv(f'samples_{domain}_validation.csv').values)

    if IS_AIRGAP:
        ranges = [
            (1.5e-3, 12e-3), # x1
            (7e-3, 23e-3),   # x2
            (5e-3, 15e-3),   # x3
            (0, 20)          # x4
        ]
    else:
        ranges = [
        (2e-3, 12e-3),  # x0
        (15e-3, 25e-3),    # x1
        (5e-3, 15e-3),    # x2
        (0, 20),           # x3 (only used if dimensions == 4)
        (1e-3, 3e-3),     # x4
        (2e-3, 4e-3),     # x5
        (18, 22),     # x6
        (0.1,0.6),     # x7
        (9, 11),           # x8
        (9, 11),           # x9
        (9, 11),           # x10
        (-5, 5),     # x11
        (-5, 5),     # x12
        (-5, 5)      # x13
        ]

    input_train_norm = normalize_data(input_train[:samples,1:].T, ranges[:in_feature]).T
    input_val_norm = normalize_data(input_val[:samples_val,1:].T, ranges[:in_feature]).T

    matrix = pd.read_csv(f'MatrixK_{domain}.csv', header=None)
    K_mat = matrix.to_numpy().astype(np.float64)
    # Perform POD
    U, eigenvalues = perform_pod(data_train, K_mat, modes)


    if np.any(np.abs(np.imag(U)) > 1e-12):
        print("Warning: POD modes contain imaginary parts, which may lead to unexpected behavior.")

    U = np.real(U)  # Ensure U is real-valued
    
    data_train_red = project_onto_pod_modes(data_train, U, K_mat).T

    # Transform data
    feature_means = np.mean(data_train_red, axis=0, keepdims=True)
    feature_stds = np.std(data_train_red, axis=0, keepdims=True)
    # Normalize the data matrix
    data_train_norm = (data_train_red - feature_means) / feature_stds 

    x_train = torch.tensor(input_train_norm, dtype=torch.float64, device=DEVICE)
    x_val = torch.tensor(input_val_norm, dtype=torch.float64, device=DEVICE)
    y_train_red = torch.tensor(data_train_red, dtype=torch.float64, device=DEVICE)
    Y_train = torch.tensor(data_train, dtype=torch.float64, device=DEVICE)
    Y_val = torch.tensor(data_val, dtype=torch.float64, device=DEVICE)

    # Read matrix A for error computation and create a tensor
    K = torch.tensor(matrix.values, dtype=torch.float64, requires_grad=False, device=DEVICE)
    K_red = torch.tensor(np.dot(U.T, np.dot(matrix.values, U)), dtype=torch.float64, requires_grad=False, device=DEVICE)

    feature_stds_t = torch.tensor(feature_stds, dtype=torch.float32, device=DEVICE)
    feature_means_t = torch.tensor(feature_means, dtype=torch.float32, device=DEVICE)
   
    U = torch.tensor(U, dtype=torch.float64, requires_grad=False, device=DEVICE)

    model = run(in_feature, modes)

    end_time = time()-start_time
    print(f'Training completed in {end_time:.2f} seconds')

    # Final comprehensive error analysis
    with torch.no_grad():
        final_results = analyze_errors(model, "Final")

    torch.save(model.state_dict(), f"model_{domain}.pt")

    return

def main(is_airgap):
    global IS_AIRGAP, EARLY_STOPPING_THRESHOLD
    IS_AIRGAP = is_airgap
    if IS_AIRGAP:
        EARLY_STOPPING_THRESHOLD = 3.13e-5
    else:
        EARLY_STOPPING_THRESHOLD = 2.2e-4 
    run_training()

# Run the main function
if __name__ == '__main__':
    # Check if an argument is provided; default to "true" if not
    arg = sys.argv[1].lower() if len(sys.argv) > 1 else "true"
    # Map common string values to a boolean
    false_vals = {"false", "0", "f", "no", "n"}
    true_vals = {"true", "1", "t", "yes", "y"}
    if arg in false_vals:
        is_airgap = False
    elif arg in true_vals:
        is_airgap = True
    else:
        # fallback: interpret unknown as True but warn
        print(f"Warning: unrecognized argument '{arg}', defaulting to airgap=True")
        is_airgap = True
    main(is_airgap)
