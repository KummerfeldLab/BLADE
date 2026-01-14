hdf5_path = "/users/2/wan00232/data_folder/Tissue_data/spaceranger_V11D13-034-D1/filtered_feature_bc_matrix.h5"
f = h5py.File(hdf5_path, 'r')

import h5py
import os
import shutil

def print_hdf5_details(hdf5_path):
	"""
	Recursively print all groups, datasets, shapes, dtypes, and attributes in an HDF5 file.
	"""
	def print_attrs(name, obj):
		print(f"{'Group' if isinstance(obj, h5py.Group) else 'Dataset'}: {name}")
		if isinstance(obj, h5py.Dataset):
			print(f"  shape: {obj.shape}, dtype: {obj.dtype}")
		if obj.attrs:
			print(f"  attributes:")
			for key, val in obj.attrs.items():
				print(f"    {key}: {val}")

	with h5py.File(hdf5_path, 'r') as f:
		f.visititems(print_attrs)

def copy_hdf5_file(src_path, dst_path):
	"""
	Copy an HDF5 file byte-for-byte to a new location.
	"""
	shutil.copy2(src_path, dst_path)
	print(f"Copied {src_path} to {dst_path}")

# Example usage:
print_hdf5_details(hdf5_path)
# copy_hdf5_file("/path/to/source.h5", "/path/to/destination.h5")