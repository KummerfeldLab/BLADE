def load_parquet(parquet_path, nrows=5):
	print(f"\nPreview of {parquet_path}:")
	df = pd.read_parquet(parquet_path)
	print(df.head(nrows))

import os
import pandas as pd
import h5py

def list_output_files(output_dir):
	print(f"\nFiles in {output_dir}:")
	for root, dirs, files in os.walk(output_dir):
		for file in files:
			print(os.path.relpath(os.path.join(root, file), output_dir))

def load_h5_matrix(h5_path):
	print(f"\nLoading HDF5 matrix: {h5_path}")
	with h5py.File(h5_path, 'r') as f:
		print("Groups:", list(f.keys()))
		# Print shape info if present
		if 'matrix' in f:
			matrix = f['matrix']
			for k in matrix.keys():
				print(f"  matrix/{k}: {matrix[k].shape if hasattr(matrix[k], 'shape') else 'N/A'}")

def load_csv(csv_path, nrows=5):
	print(f"\nPreview of {csv_path}:")
	df = pd.read_csv(csv_path, nrows=nrows)
	print(df.head())

def inspect_visium_hd_output(output_dir):
	list_output_files(output_dir)

	# Try to load common files
	h5_files = [f for f in os.listdir(output_dir) if f.endswith('.h5')]
	for h5_file in h5_files:
		load_h5_matrix(os.path.join(output_dir, h5_file))

	csv_files = [f for f in os.listdir(output_dir) if f.endswith('.csv')]
	for csv_file in csv_files:
		load_csv(os.path.join(output_dir, csv_file))

	parquet_files = [f for f in os.listdir(output_dir) if f.endswith('.parquet')]
	for parquet_file in parquet_files:
		load_parquet(os.path.join(output_dir, parquet_file))

if __name__ == "__main__":
	# Example usage: change this path to your output directory
	output_dir = "/users/2/wan00232/data_folder/Visium_HD"
	inspect_visium_hd_output(output_dir)



import pyarrow.parquet as pq
import pandas as pd

# Read the Parquet file into an Arrow Table
table = pq.read_table('/users/2/wan00232/data_folder/Visium_HD/binned_outputs/square_002um/spatial/tissue_positions.parquet')

table['cell_id'][1:10]


/users/2/wan00232/data_folder/Visium_HD/segmented_outputs

import h5py
import numpy as np

file_path = '/users/2/wan00232/data_folder/Visium_HD/binned_outputs/filtered_feature_cell_matrix.h5'
with h5py.File(file_path, 'r') as f:
    # File operations go here
    pass