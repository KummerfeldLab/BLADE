import h5py
import os

def inspect_hdf5_structure(hdf5_path):
    """
    Save the structure of an HDF5 file (groups, datasets, and their shapes/dtypes) to a text file.
    """
    if not os.path.exists(hdf5_path):
        raise FileNotFoundError(f"File not found: {hdf5_path}")

    # Save output in the hdf5_utils folder in the repo
    repo_dir = os.path.dirname(os.path.abspath(__file__))
    base_name = os.path.basename(hdf5_path)
    output_path = os.path.join(repo_dir, f"{base_name}_structure.txt")
    lines = []

    def collect_structure(name, obj):
        if isinstance(obj, h5py.Group):
            lines.append(f"Group: {name}")
        elif isinstance(obj, h5py.Dataset):
            lines.append(f"Dataset: {name}, shape: {obj.shape}, dtype: {obj.dtype}")

    with h5py.File(hdf5_path, 'r') as f:
        f.visititems(collect_structure)

    with open(output_path, 'w') as out:
        for line in lines:
            out.write(line + '\n')


    print(f"Structure saved to {output_path}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python inspect_hdf5.py <path_to_hdf5_file>")
        sys.exit(1)
    hdf5_path = sys.argv[1]
    inspect_hdf5_structure(hdf5_path)
