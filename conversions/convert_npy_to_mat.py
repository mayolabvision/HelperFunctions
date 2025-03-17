import numpy as np
import os
from scipy.io import savemat

def convert_npy_to_mat(directory):
    for file in os.listdir(directory):
        if file.endswith('.npy'):
            npy_path = os.path.join(directory, file)
            try:
                if file != 'ops.npy':
                    data = np.load(npy_path, allow_pickle=True)

                    mat_path = os.path.join(directory, file.replace('.npy', '.mat'))
                    savemat(mat_path, {file.replace('.npy', ''): data})
                    print(f"Converted: {file} -> {file.replace('.npy', '.mat')}")

            except Exception as e:
                print(f"Failed to convert {file}: {e}")
