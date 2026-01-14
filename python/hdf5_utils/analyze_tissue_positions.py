import polars as pl


def analyze_tissue_positions(csv_path):
    df = pl.read_csv(csv_path)
    print(f"Total rows: {df.height}")
    in_tissue_count = df.filter(pl.col("in_tissue") == 1).height
    print(f"Rows with in_tissue == 1: {in_tissue_count}")


if __name__ == "__main__":
    analyze_tissue_positions("hdf5_utils/tissue_positions.csv")
