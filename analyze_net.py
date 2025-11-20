# SPDX-License-Identifier: GPL-3.0

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import os
import sys
import time

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


def calc_error(model,x,Y,K):
    model.eval()
    with torch.no_grad():
        y_test_norm = model(x)
        y_e = (y_test_norm * feature_stds_t) + feature_means_t
        y_e =  U_t @ y_e.T 
        errors = relative_error(Y,y_e,K)
        errors = np.array(errors)
    model.train()
    return errors

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

def define_model(features, input_dim, modes):
    layers = []
    in_features = input_dim
    for h in features:
        layers.append(nn.Linear(in_features, h))
        layers.append(nn.ReLU())
        in_features = h
    layers.append(nn.Linear(in_features, modes))
    return nn.Sequential(*layers)

def analyze_net():
    global U_t, K, x_train, Y_train, x_val, Y_val, feature_means_t, feature_stds_t
    
    if PAPER:
        opt = "_paper"
    else:
        opt = ""
    
    if IS_AIRGAP:
        domain = 'gap'
        modes =14
        in_feature=4
        hidden_dim = [190, 110, 180]

    else:
        domain = 'full'
        modes = 90 
        in_feature=14 
        hidden_dim = [336, 336, 336]

    # Define paths and file names
    train_folder_path = f'sol_{domain}_training{opt}'
    train_file_name = f'sol_{domain}_training_'
    val_folder_path = f'sol_{domain}_validation{opt}'
    val_file_name = f'sol_{domain}_validation_'
    test_folder_path = f'sol_{domain}_testing{opt}'
    test_file_name = f'sol_{domain}_testing_'

    samples = 1024 # Number of samples
    samples_val = 128  # Number of validation samples
    samples_test = 128 # Number of test samples
    
    print(f'Analyze the neural network for {"airgap" if IS_AIRGAP else "full"} model')

    # Read data
    data_train = create_data_matrix(train_folder_path, train_file_name, samples);
    data_val = create_data_matrix(val_folder_path, val_file_name, samples_val);
    data_test = create_data_matrix(test_folder_path, test_file_name, samples_test);
   
    input_train = np.array(pd.read_csv(f'samples_{domain}_training.csv').values)
    input_val = np.array(pd.read_csv(f'samples_{domain}_validation.csv').values)
    input_test = np.array(pd.read_csv(f'samples_{domain}_testing.csv').values)

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
    input_test_norm = normalize_data(input_test[:samples_test,1:].T, ranges[:in_feature]).T

  
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

    # Convert to Torch Tensor
    x_train = torch.tensor(input_train_norm, dtype=torch.float32, device=DEVICE)
    x_val = torch.tensor(input_val_norm, dtype=torch.float32, device=DEVICE)  
    x_test = torch.tensor(input_test_norm, dtype=torch.float32, device=DEVICE)
    Y_train = torch.tensor(data_train, dtype=torch.float64, device=DEVICE)
    Y_val = torch.tensor(data_val, dtype=torch.float64, device=DEVICE)
    Y_test = torch.tensor(data_test, dtype=torch.float64, device=DEVICE)  

    # Read matrix A for error computation and create a tensor
    K = torch.tensor(matrix.values, dtype=torch.float64, requires_grad=False, device=DEVICE)
    
    feature_stds_t = torch.tensor(feature_stds, dtype=torch.float64, device=DEVICE)
    feature_means_t = torch.tensor(feature_means, dtype=torch.float64, device=DEVICE)

    U_t = torch.tensor(U, dtype=torch.float64, requires_grad=False, device=DEVICE)

    # --- Load and evaluate the model ---
    model = define_model(hidden_dim, in_feature,modes);
    model.load_state_dict(torch.load(f"model_{domain}{opt}.pt", map_location=DEVICE))
    model.eval()

    errors_train = calc_error(model, x_train, Y_train, K)
    errors_val = calc_error(model, x_val, Y_val, K)
    errors_test = calc_error(model, x_test, Y_test, K)
   
    errors = {'train': errors_train, 'test': errors_test, 'val': errors_val}

    # Print statistics for each error set
    for key in ['train', 'test', 'val']:
        arr = np.array(errors[key])
        print(f"{key.capitalize()} error: mean={arr.mean():.4e}, std={arr.std():.4e}, max={arr.max():.4e}")
        

    times = []
    # Compute average prediction times
    for i in range(samples_test):
        start_time = time.time()
        z = model(torch.tensor(input_val_norm[i,:], dtype=torch.float32, device=DEVICE))
        z = z.detach().cpu().numpy()
        z =  z *feature_stds + feature_means
        x = U @ z.T
        end_time = time.time()
        elapsed = end_time - start_time
        times.append(elapsed)
    
    print(f"Average prediction time: {np.mean(times):.6f} seconds per sample.")

    # Save training, testing and validation predictions
    datasets = [
        ('training', input_train_norm, samples),
        ('validation',   input_val_norm, samples_val),
        ('testing',  input_test_norm, samples_test),
    ]

    times_summary = {}
    for key, input_norm, n_samples in datasets:
        folder = f'Predictions{opt}_{domain}_{key}'
        os.makedirs(folder, exist_ok=True)
        times = []
        for i in range(n_samples):
            inp = torch.tensor(input_norm[i, :], dtype=torch.float32, device=DEVICE)
            with torch.no_grad():
                z = model(inp).detach().cpu().numpy().reshape(-1)
            # denormalize and reconstruct to full space
            z = z * feature_stds + feature_means
            x = U @ z.T
            # save prediction
            out_path = os.path.join(folder, f'{domain}_{key}_{i+1}.csv')
            pd.DataFrame(x).to_csv(out_path, index=False, header=False)
        print(f"Saved {n_samples} of {key} predictions to '{folder}'.")


def main(is_airgap, paper):
    global IS_AIRGAP, PAPER
    IS_AIRGAP = is_airgap
    PAPER = paper
    analyze_net()

# Run the main function
if __name__ == '__main__':
    # Check if an argument is provided; default to "true" if not
    arg = sys.argv[1].lower() if len(sys.argv) > 1 else "true"
    arg_paper = sys.argv[2].lower() if len(sys.argv) > 2 else "true"
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
    if arg_paper in false_vals:
        paper = False
    elif arg_paper in true_vals:
        paper = True
    else:
        # fallback: interpret unknown as True but warn
        print(f"Warning: unrecognized argument '{arg_paper}', defaulting to paper=True")
        paper = True
    main(is_airgap, paper)
