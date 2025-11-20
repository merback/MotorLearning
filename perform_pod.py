# SPDX-License-Identifier: GPL-3.0

import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt



# Function to read all CSV files and create a data matrix
def create_data_matrix(folder_path, file_name, num_samples):
    data = []
    if IS_AIRGAP:
        n= 462
    else:
        n= 6354
    for i in range(num_samples):
        d = pd.read_csv(folder_path + '/' + file_name + str(i+1) + '.csv', header=None) 
        x = np.reshape(d.values, n)
        data.append(x)
    return np.array(data).T


# Perform POD 
def perform_pod(data_matrix, K, max_modes):
   
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
    for i in range(max_modes):
        phi = modes[:, i]
        norm = np.sqrt(eigvals[i])
        if norm > 1e-12:
            modes[:, i] = phi / norm
        else:
            modes[:, i] = 0.0

    return modes, eigvals


# Function to determine the number of eigenvalues needed for a given percentage of information
def determine_eigenvalues_needed(eigenvalues, percentage):
    cumulative_sum = np.cumsum(eigenvalues)
    total_sum = np.sum(eigenvalues)
    threshold = percentage * total_sum
    num_eigenvalues = np.searchsorted(cumulative_sum, threshold) + 1
    return num_eigenvalues


# Project data onto the top k POD modes
def project_onto_pod_modes(data_matrix, U, K, k):
    # Phi: (features, k), data_matrix: (features, samples)
    return U[:, :k].T @ K @ data_matrix

# Reconstruct data from the top k POD modes
def reconstruct_from_pod_modes(projections, U, k):
    return U[:, :k] @ projections

def compute_relative_error(original, reconstructed, K):
    errors = []
    for i in range(original.shape[1]):
        x = reconstructed[:, i]
        y = original[:, i]
        numerator = (x.T - y.T) @ K @ (x - y)
        denominator = y.T @ K @ y
        error = np.sqrt(numerator / denominator)
        errors.append(error)
    
    errors = np.array(errors)
    mean_error = np.mean(errors)
    std_error = np.std(errors)
    max_error = np.max(errors)
    
    return mean_error, std_error, max_error


# Main function
def main(is_airgap=True):
    global IS_AIRGAP
    IS_AIRGAP = is_airgap

    if IS_AIRGAP:
        domain = 'gap'
        threshold = 0.001
    else:
        domain = 'full'
        threshold = 0.006

    print(f'Running POD analysis for {"airgap" if IS_AIRGAP else "full"} model with threshold {threshold*100}%')

    # Define paths and file names
    train_folder_path = f'sol_{domain}_training'
    train_file_name = f'sol_{domain}_training_'
    val_folder_path = f'sol_{domain}_validation'
    val_file_name = f'sol_{domain}_validation_'

    samples_train = 1024
    samples_val = 128
    
    # Create data matrix
    data_matrix = create_data_matrix(train_folder_path, train_file_name, samples_train)
    print(f'Size of data matrix: {data_matrix.shape[0]} features, {data_matrix.shape[1]} samples')

    K =  pd.read_csv(f'MatrixK_{domain}.csv', header=None).values 
    
    max_modes = 100

    # Perform POD
    u, eigenvalues = perform_pod(data_matrix, K, max_modes)

    # Project the data onto the top k POD modes and analyze reconstruction error for different k
    mean_errors = []
    max_errors = []
    std_errors = []

    val_matrix = create_data_matrix(val_folder_path, val_file_name,samples_val)
    
    min_modes = 1
    for k in range(min_modes, max_modes + 1):
        projections = project_onto_pod_modes(val_matrix, u, K, k)
        reconstructed = reconstruct_from_pod_modes(projections, u, k)

        # Compute relative error
        mean_error, std_error, max_error = compute_relative_error(val_matrix, reconstructed, K)
        mean_errors.append(mean_error)
        max_errors.append(max_error)
        std_errors.append(std_error)

    # Plot mean reconstruction error vs number of modes
    plt.figure(figsize=(8, 5))
    plt.semilogy(range(min_modes, max_modes + 1), np.array(mean_errors), label='Mean Reconstruction Error')
    plt.axvline(x=14, color='red', linestyle='--', label='14 Modes') 
    plt.xlabel('Number of POD Modes')
    plt.ylabel('Mean Relative Reconstruction Error')
    plt.title('Mean Reconstruction Error vs Number of POD Modes')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Get number of modes for threshhold.
    indices = np.where(np.array(mean_errors) < threshold)[0]
    if len(indices) > 0:
        k = min_modes +indices[0]
   
    elif max_modes < min(data_matrix.shape[0], data_matrix.shape[1]):
        raise RuntimeError(f"No number of modes below max_modes ({max_modes}) achieves the desired mean error threshold, raise the maximum number of modes.")
    else:
        raise RuntimeError("No number of modes achieves the desired mean error threshold.")
    
    print(f'Number of POD modes needed to achieve mean relative error below {threshold*100}%: {k}')
    validation_projections = project_onto_pod_modes(val_matrix, u, K, k)
    validation_reconstructed = reconstruct_from_pod_modes(validation_projections, u, k)
    # Compute and print relative error
    mean_error, std_error, max_error= compute_relative_error(val_matrix, validation_reconstructed,K)
    print(f'Max relative reconstruction error for the validation set using top {k} POD modes: {max_error* 100:.4f}%')
    print(f'Mean relative reconstruction error for the validation set using top {k} POD modes: {mean_error* 100:.4f}%')
    print(f'Std relative reconstruction error for the validation set using top {k} POD modes: {std_error* 100:.4f}%')


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
