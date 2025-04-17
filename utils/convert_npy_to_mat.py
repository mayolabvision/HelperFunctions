import numpy as np
import os
from scipy.io import savemat

def convert_npy_to_mat(directory):
    if not os.path.isdir(directory):
        print(f"Directory does not exist: {directory}")
        return  # Exit the function gracefully

    for file in os.listdir(directory):
        if file.endswith('.npy') and file != 'ops.npy' and file != 'pc_features.npy':
            npy_path = os.path.join(directory, file)
            mat_path = os.path.join(directory, file.replace('.npy', '.mat'))

            if os.path.exists(mat_path):
                print(f"Skipped (already exists): {file} -> {os.path.basename(mat_path)}")
                continue

            try:
                data = np.load(npy_path, allow_pickle=True)
                var_name = file.replace('.npy', '')

                try:
                    savemat(mat_path, {var_name: data})
                    print(f"Converted: {file} -> {os.path.basename(mat_path)}")
                except Exception as e:
                    print(f"Initial save failed for {file}: {e}")
                    try:
                        # Try converting to float64
                        data_converted = data.astype(np.float64)
                        savemat(mat_path, {var_name: data_converted})
                        print(f"Converted with float64 fallback: {file} -> {os.path.basename(mat_path)}")
                    except Exception as e2:
                        print(f"Failed to convert {file} even with float64 fallback: {e2}")

            except Exception as e:
                print(f"Failed to load {file}: {e}")
