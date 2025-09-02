import pandas as pd

def analyze_tissue_positions(csv_path):
    df = pd.read_csv(csv_path, header=0)
    print(f"Total rows: {len(df)}")
    in_tissue_count = (df["in_tissue"] == 1).sum()
    print(f"Rows with in_tissue == 1: {in_tissue_count}")

if __name__ == "__main__":
    analyze_tissue_positions("hdf5_utils/tissue_positions.csv")
