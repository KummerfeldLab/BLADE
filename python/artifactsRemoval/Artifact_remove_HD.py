import os
import copy
import warnings
import numpy as np
import polars as pl
import h5py
from scipy.sparse import csc_matrix

from artifactsRemoval.polars_neighbors import NeighborIndex


class TissueHD:
    """Visium HD reader that loads matrices from 10x-style HDF5 files."""

    def __init__(
        self,
        data_dir=None,
        h5_path=None,
        positions_path=None,
        neighbor_mode="k",
        k_neighbors=8,
        neighbor_radius=None,
    ):
        self.data_dir = data_dir
        self.h5_path = h5_path or self._auto_find_h5(data_dir)
        if self.h5_path is None:
            raise FileNotFoundError(
                "HDF5 file not provided and automatic discovery failed. "
                "Pass `h5_path` or set `data_dir` to a folder containing "
                "filtered_feature_cell_matrix.h5 or filtered_feature_bc_matrix.h5."
            )

        (
            self.tissue_matrix,
            self.barcode_list,
            self.feature_list,
        ) = self._read_h5_matrix(self.h5_path)

        pos_path = positions_path or self._auto_find_positions(data_dir, self.h5_path)
        self.tissue_position = self._read_positions(pos_path)

        self.dict_barcode_to_coor = self._barcode_to_coor()
        self.dict_coor_to_barcode = self._coor_to_barcode()
        self.dict_barcode_to_column = self._barcode_to_column()
        self.dict_feature_to_row = self._feature_to_row()

        self._neighbor_mode = neighbor_mode
        self._k_neighbors = int(k_neighbors) if k_neighbors else 8
        self._user_radius = float(neighbor_radius) if neighbor_radius else None
        self._neighbor_index = NeighborIndex(
            self.tissue_position,
            default_mode=self._neighbor_mode,
            radius=self._user_radius,
            k=self._k_neighbors,
        )

    def _auto_find_h5(self, data_dir):
        if not data_dir:
            return None
        candidates = [
            os.path.join(data_dir, "filtered_feature_cell_matrix.h5"),
            os.path.join(data_dir, "filtered_feature_bc_matrix.h5"),
        ]
        for cand in candidates:
            if os.path.exists(cand):
                return cand
        return None

    def _auto_find_positions(self, data_dir, h5_path):
        if data_dir:
            root = data_dir
        elif h5_path:
            root = os.path.dirname(h5_path)
        else:
            return None
        candidates = [
            os.path.join(root, "spatial", "tissue_positions.parquet"),
            os.path.join(root, "spatial", "tissue_positions.csv"),
        ]
        for cand in candidates:
            if os.path.exists(cand):
                return cand
        return None

    def _read_h5_matrix(self, h5_path):
        with h5py.File(h5_path, "r") as f:
            grp = f["matrix"]
            data = grp["data"][:]
            indices = grp["indices"][:]
            indptr = grp["indptr"][:]
            shape = tuple(grp["shape"][:].tolist())
            barcodes = grp["barcodes"][:]

            feat_grp = grp["features"]
            feature_id = feat_grp["id"][:]
            feature_name = feat_grp["name"][:]
            feature_type = feat_grp["feature_type"][:] if "feature_type" in feat_grp else None

        matrix = csc_matrix((data, indices, indptr), shape=shape)

        def _decode(arr):
            return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr]

        barcode_df = pl.DataFrame({"barcode": _decode(barcodes)})
        feat_cols = {"feature_id": _decode(feature_id), "feature_name": _decode(feature_name)}
        if feature_type is not None:
            feat_cols["feature_type"] = _decode(feature_type)
        feature_df = pl.DataFrame(feat_cols)
        return matrix, barcode_df, feature_df

    def _read_positions(self, positions_path):
        if positions_path is None:
            raise FileNotFoundError(
                "Positions file not provided and automatic discovery failed. "
                "Provide `positions_path` pointing to tissue_positions.parquet or CSV."
            )

        if positions_path.endswith(".parquet"):
            try:
                df = pl.read_parquet(positions_path)
            except Exception as exc:
                raise ImportError(
                    "Reading Parquet positions requires the 'pyarrow' package. "
                    "Install pyarrow or provide a CSV positions file."
                ) from exc
        else:
            df = pl.read_csv(positions_path)

        id_col = None
        for cand in ["barcode", "cell_id", "bin_id", "barcode_id"]:
            if cand in df.columns:
                id_col = cand
                break
        if id_col is None:
            raise ValueError("Positions file must contain a barcode/cell identifier column.")

        coord_pairs = [
            ("x_mtx", "y_mtx"),
            ("x", "y"),
            ("pxl_col_in_fullres", "pxl_row_in_fullres"),
        ]
        x_col = y_col = None
        for x_name, y_name in coord_pairs:
            if x_name in df.columns and y_name in df.columns:
                x_col, y_col = x_name, y_name
                break
        if x_col is None:
            raise ValueError("Positions file must contain coordinate columns (e.g. x_mtx/y_mtx).")

        if "in_tissue" not in df.columns:
            df = df.with_columns(pl.lit(1).alias("in_tissue"))

        out = df.select(
            [
                pl.col(id_col).cast(pl.Utf8).alias("barcode"),
                pl.col("in_tissue").cast(pl.Int64).alias("in_tissue"),
                pl.col(x_col).cast(pl.Float64).alias("x_mtx"),
                pl.col(y_col).cast(pl.Float64).alias("y_mtx"),
                pl.col(x_col).cast(pl.Float64).alias("x_coor"),
                pl.col(y_col).cast(pl.Float64).alias("y_coor"),
            ]
        )
        return out

    def _barcode_to_coor(self):
        bars = self.tissue_position["barcode"].to_list()
        xs = self.tissue_position["x_mtx"].to_list()
        ys = self.tissue_position["y_mtx"].to_list()
        return {b: [x, y] for b, x, y in zip(bars, xs, ys)}

    def _coor_to_barcode(self):
        mapping = {}
        for row in self.tissue_position.iter_rows(named=True):
            key = f"{row['x_mtx']} {row['y_mtx']}"
            mapping[key] = row["barcode"]
        return mapping

    def _barcode_to_column(self):
        bars = self.barcode_list["barcode"].to_list()
        return {b: idx for idx, b in enumerate(bars)}

    def _feature_to_row(self):
        genes = self.feature_list.get_column("feature_name").to_list()
        return {g: idx for idx, g in enumerate(genes)}

    def neighbors_by_barcode(self, barcode, radius=None, k=None, include_center=False):
        return self._neighbor_index.neighbors_by_barcode(
            barcode, radius=radius, k=k, include_center=include_center
        )


class ArtifactDetectHD(TissueHD):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_sum(self):
        col_sums = np.asarray(self.tissue_matrix.sum(axis=0)).ravel()
        bars = []
        xs = []
        ys = []
        sums = []
        for row in self.tissue_position.iter_rows(named=True):
            b = row["barcode"]
            idx = self.dict_barcode_to_column.get(b)
            if idx is None:
                continue
            bars.append(b)
            xs.append(row["x_mtx"])
            ys.append(row["y_mtx"])
            sums.append(col_sums[idx])
        return pl.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": sums})

    def seqdepth(self):
        df = pl.DataFrame(
            self.tissue_matrix.toarray().transpose(),
            schema=self.feature_list.get_column("feature_name").to_list(),
        )
        df = df.with_columns(pl.Series("barcode", self.barcode_list["barcode"].to_list()))
        return df

    def get_neighbour_1(self, barcode):
        bars = self.neighbors_by_barcode(barcode)
        return [self.dict_barcode_to_column.get(b) for b in bars]

    def get_neighbour_1_bar(self, barcode):
        return self.neighbors_by_barcode(barcode)

    def get_border(self, expected_neighbors=8, strategy="extreme", border_width=0.0, include_counts=True):
        mask_df = getattr(self, "spot_inclusion_condition", self.tissue_position.select(["barcode", "in_tissue"]))
        covered_barcodes = set(mask_df.filter(pl.col("in_tissue") == 1)["barcode"].to_list())
        covered = self.tissue_position.filter(pl.col("barcode").is_in(covered_barcodes))
        if covered.height == 0:
            return pl.DataFrame({"barcode": []})

        bars = []
        xs = []
        ys = []
        gene_count = []

        if strategy == "extreme":
            x_vals = covered["x_mtx"].to_numpy()
            y_vals = covered["y_mtx"].to_numpy()
            xmin, xmax = np.min(x_vals), np.max(x_vals)
            ymin, ymax = np.min(y_vals), np.max(y_vals)
            bw = float(border_width or 0.0)
            for row in covered.iter_rows(named=True):
                x = float(row["x_mtx"])
                y = float(row["y_mtx"])
                if not (
                    abs(x - xmin) <= bw
                    or abs(x - xmax) <= bw
                    or abs(y - ymin) <= bw
                    or abs(y - ymax) <= bw
                ):
                    continue
                b = row["barcode"]
                idx = self.dict_barcode_to_column.get(b)
                if idx is None:
                    continue
                bars.append(b)
                xs.append(x)
                ys.append(y)
                if include_counts:
                    gene_count.append(self.tissue_matrix[:, idx].sum())
        else:
            rad = self._user_radius if self._user_radius is not None else self._neighbor_index.default_radius
            for row in covered.iter_rows(named=True):
                b = row["barcode"]
                ngh = self.neighbors_by_barcode(b, radius=rad)
                if len(ngh) < expected_neighbors:
                    idx = self.dict_barcode_to_column.get(b)
                    if idx is None:
                        continue
                    bars.append(b)
                    xs.append(row["x_mtx"])
                    ys.append(row["y_mtx"])
                    if include_counts:
                        gene_count.append(self.tissue_matrix[:, idx].sum())

        if include_counts:
            return pl.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": gene_count})
        return pl.DataFrame({"barcode": bars})

    def get_edge(self, include_counts=True, radius=None):
        mask_df = getattr(self, "spot_inclusion_condition", self.tissue_position.select(["barcode", "in_tissue"]))
        zero_barcodes = set(mask_df.filter(pl.col("in_tissue") == 0)["barcode"].to_list())
        covered_barcodes = set(mask_df.filter(pl.col("in_tissue") == 1)["barcode"].to_list())

        neighbors = set()
        for zb in zero_barcodes:
            ngh = self.neighbors_by_barcode(zb, radius=radius)
            neighbors.update(ngh)
        neighbors.intersection_update(covered_barcodes)
        neighbors.discard(None)

        if not include_counts:
            return pl.DataFrame({"barcode": list(neighbors)})

        bars = []
        xs = []
        ys = []
        gene_count = []
        for b in neighbors:
            coor = self.dict_barcode_to_coor.get(b)
            idx = self.dict_barcode_to_column.get(b)
            if coor is None or idx is None:
                continue
            bars.append(b)
            xs.append(coor[0])
            ys.append(coor[1])
            gene_count.append(self.tissue_matrix[:, idx].sum())
        return pl.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": gene_count})


class Artifact_remove_HD(ArtifactDetectHD):
    def __init__(self, *args, compute_depth_on_init=False, max_depth_limit=4, **kwargs):
        super().__init__(*args, **kwargs)
        self.spot_inclusion_condition = self._init_inclusion_condition()
        self.tissue_depth_df = None
        self.tissue_depth = 0
        if compute_depth_on_init:
            self.tissue_depth_df = self.fn_spot_dis_to_edge()
        self.remove_procedure = []
        self.removed_spots = {}
        self.remaining_depth = self.tissue_depth
        self.max_depth_limit = int(max_depth_limit)

    def _init_inclusion_condition(self):
        return self.tissue_position.select(["barcode", "in_tissue"]).clone()

    def fn_spot_dis_to_edge(self):
        bars = self.tissue_position["barcode"].to_list()
        in_tissue = self.tissue_position["in_tissue"].to_list()
        depth_map = {b: 0 for b in bars}

        mask_df = self.spot_inclusion_condition
        zero_barcodes = set(mask_df.filter(pl.col("in_tissue") == 0)["barcode"].to_list())
        covered_barcodes = set(mask_df.filter(pl.col("in_tissue") == 1)["barcode"].to_list())

        frontier = set()
        for zb in zero_barcodes:
            frontier.update(self.get_neighbour_1_bar(zb))

        border_barcodes = self.get_border(include_counts=False)["barcode"].to_list()
        frontier.update(border_barcodes)
        frontier.intersection_update(covered_barcodes)
        frontier.discard(None)

        for b in frontier:
            depth_map[b] = 1
        visited = set(frontier)
        level = 2
        remaining = covered_barcodes - visited

        while remaining:
            next_frontier = set()
            for b in frontier:
                next_frontier.update(self.get_neighbour_1_bar(b))
            next_frontier.intersection_update(covered_barcodes)
            next_frontier.discard(None)
            next_frontier -= visited
            if not next_frontier:
                break
            for b in next_frontier:
                depth_map[b] = level
            visited.update(next_frontier)
            frontier = next_frontier
            remaining = covered_barcodes - visited
            level += 1

        self.tissue_depth = level
        self.tissue_depth_df = pl.DataFrame(
            {"barcode": bars, "in_tissue": in_tissue, "spot_depth": [depth_map[b] for b in bars]}
        )
        return self.tissue_depth_df

    def _set_inclusion_zero(self, barcodes):
        if not barcodes:
            return
        mask = pl.col("barcode").is_in(barcodes)
        self.spot_inclusion_condition = self.spot_inclusion_condition.with_columns(
            pl.when(mask).then(0).otherwise(pl.col("in_tissue")).alias("in_tissue")
        )

    def remove_border(self, distance=None, strategy="extreme", expected_neighbors=8):
        self.remove_procedure.append("border")
        if strategy == "extreme":
            df_border = self.get_border(strategy="extreme", border_width=(distance or 0.0))
        else:
            df_border = self.get_border(strategy="degree", expected_neighbors=expected_neighbors)
        self._set_inclusion_zero(df_border["barcode"].to_list())
        msg = f"We removed {df_border.height} border spots"
        if strategy == "extreme":
            msg += f" within border width {distance or 0.0}"
        else:
            msg += f" (degree < {expected_neighbors})"
        print(msg)

    def remove_edge(self, distance=1, strategy="extreme", expected_neighbors=8, seed_mode="zeros", border_width=None):
        self.remove_procedure.append("edge")
        if distance <= 0:
            print("Distance must be >= 1; nothing removed.")
            return
        if distance > self.max_depth_limit:
            print(f"Requested distance {distance} exceeds max limit {self.max_depth_limit}; clamping.")
            distance = self.max_depth_limit

        mask_df = self.spot_inclusion_condition
        covered_barcodes = set(mask_df.filter(pl.col("in_tissue") == 1)["barcode"].to_list())

        frontier = set()
        if seed_mode in ("zeros", "both"):
            zero_barcodes = mask_df.filter(pl.col("in_tissue") == 0)["barcode"].to_list()
            margin = self._user_radius if self._user_radius is not None else (self._neighbor_index.default_radius or 0.0)
            band = margin * 2.0 if margin > 0 else None
            if band is not None and covered_barcodes:
                cov_pos = self.tissue_position.filter(pl.col("barcode").is_in(covered_barcodes))
                xmin = cov_pos["x_mtx"].min()
                xmax = cov_pos["x_mtx"].max()
                ymin = cov_pos["y_mtx"].min()
                ymax = cov_pos["y_mtx"].max()
            else:
                xmin = ymin = xmax = ymax = None

            for b in zero_barcodes:
                if band is not None and xmin is not None:
                    coor = self.dict_barcode_to_coor.get(b)
                    if coor is None:
                        continue
                    if not (
                        xmin - band <= coor[0] <= xmax + band
                        and ymin - band <= coor[1] <= ymax + band
                    ):
                        continue
                frontier.update(self.get_neighbour_1_bar(b))

        if seed_mode in ("border", "both"):
            if strategy == "extreme":
                bw = border_width if border_width is not None else (
                    self._user_radius if self._user_radius is not None else (self._neighbor_index.default_radius or 0.0)
                )
                border_df = self.get_border(strategy="extreme", border_width=bw, include_counts=False)
            else:
                border_df = self.get_border(
                    strategy="degree", expected_neighbors=expected_neighbors, include_counts=False
                )
            frontier.update(border_df["barcode"].to_list())

        frontier.intersection_update(covered_barcodes)
        frontier.discard(None)

        visited = set(frontier)
        current_frontier = frontier
        for _ in range(1, distance):
            next_frontier = set()
            for b in current_frontier:
                next_frontier.update(self.get_neighbour_1_bar(b))
            next_frontier.intersection_update(covered_barcodes)
            next_frontier.discard(None)
            next_frontier -= visited
            if not next_frontier:
                break
            visited.update(next_frontier)
            current_frontier = next_frontier

        mask = self.spot_inclusion_condition["barcode"].is_in(visited)
        removed = int(self.spot_inclusion_condition.filter(mask)["in_tissue"].sum())
        self._set_inclusion_zero(list(visited))
        print(f"We removed {removed} edge spots within {distance} layer(s)")

    def remove_malfunction(self, use_clusters=True, min_cluster_size=15, radius=None, radius_scale=1.0):
        if use_clusters:
            return self.remove_malfunction_clustered(
                min_cluster_size=min_cluster_size, radius=radius, radius_scale=radius_scale
            )

        self.remove_procedure.append("malfunction")
        df_outlier = Artifact_remove_HD.outlier(self.get_sum())
        self._set_inclusion_zero(df_outlier["barcode"].to_list())
        print(f"We removed {df_outlier.height} outlier spots (non-clustered)")

    def remove_malfunction_clustered(self, min_cluster_size=15, radius=None, radius_scale=1.0):
        self.remove_procedure.append("malfunction_clustered")
        df_out = Artifact_remove_HD.outlier(self.get_sum())
        include_mask = set(self.spot_inclusion_condition.filter(pl.col("in_tissue") == 1)["barcode"].to_list())
        outliers = [b for b in df_out["barcode"].to_list() if b in include_mask]
        if not outliers:
            print("No outlier spots among currently included spots.")
            return

        coords = []
        kept_barcodes = []
        for b in outliers:
            coor = self.dict_barcode_to_coor.get(b)
            if coor is None:
                continue
            kept_barcodes.append(b)
            coords.append(coor)
        if not coords:
            print("No coordinates found for outlier spots.")
            return

        pts = np.asarray(coords, dtype=float)
        base_r = radius
        if base_r is None:
            base_r = self._user_radius if self._user_radius is not None else self._neighbor_index.default_radius
        if base_r is None:
            base_r = 1.0
        rad = float(base_r) * float(radius_scale)

        n = len(kept_barcodes)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        # naive clustering with radius search using numpy distances
        for i in range(n):
            dists = np.linalg.norm(pts - pts[i], axis=1)
            inds = np.where(dists <= rad)[0].tolist()
            for j in inds:
                if j != i:
                    union(i, j)

        clusters = {}
        for i in range(n):
            clusters.setdefault(find(i), []).append(kept_barcodes[i])

        removed = 0
        for members in clusters.values():
            if len(members) >= int(min_cluster_size):
                mask = self.spot_inclusion_condition["barcode"].is_in(members)
                removed += int(self.spot_inclusion_condition.filter(mask)["in_tissue"].sum())
                self._set_inclusion_zero(members)

        print(
            f"Clustered malfunction removal: considered {n} outliers, removed {removed} spots "
            f"in clusters >= {min_cluster_size} (radius={rad:.3f})."
        )

    def simple_cleanser(self):
        self.remove_procedure = []
        self.removed_spots = {}
        self.spot_inclusion_condition = self._init_inclusion_condition()
        print("Back to the initial Tissue Sample")

    def review_removing(self):
        for proc in self.remove_procedure:
            print(proc)

    def save(self, out_dir):
        df = self.tissue_position.clone()
        df = df.drop("in_tissue")
        df = df.join(self.spot_inclusion_condition, on="barcode", how="left")
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "tissue_positions.csv")
        df.write_csv(out_path)
        print(f"Saved the tissue position file to {out_path}")
        return df

    def save_hdf5(self, out_path, filter_in_tissue=True):
        if not out_path:
            raise ValueError("out_path must be provided")

        out_dir = os.path.dirname(out_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        matrix = self.tissue_matrix
        if not hasattr(matrix, "tocsc"):
            raise ValueError("tissue_matrix must be a scipy sparse matrix")
        matrix = matrix.tocsc()
        total_spots = matrix.shape[1]

        barcodes_series = pl.Series(self.barcode_list["barcode"].to_list())
        if filter_in_tissue and hasattr(self, "spot_inclusion_condition"):
            inclusion = (
                self.spot_inclusion_condition
                .join(barcodes_series.to_frame("barcode"), on="barcode", how="right")
                .select(pl.col("in_tissue").fill_null(0).cast(pl.Int64))
            )
            keep_mask = inclusion.to_numpy().astype(bool).ravel()
            matrix = matrix[:, keep_mask]
            matrix = matrix.tocsc()
            barcodes_series = barcodes_series.filter(pl.Series(keep_mask))
        retained_spots = matrix.shape[1]
        barcodes = barcodes_series.to_numpy().astype("S")

        features = self.feature_list
        feature_id = features.get_column("feature_id").to_numpy().astype("S")
        feature_name = features.get_column("feature_name").to_numpy().astype("S")
        if "feature_type" in features.columns:
            feature_type = features.get_column("feature_type").to_numpy().astype("S")
        else:
            feature_type = np.array([b"Gene Expression"] * len(features), dtype="S")
        genome = np.array([b"unknown"] * len(features), dtype="S")
        all_tag_keys = np.array([b"feature_type"], dtype="S")

        with h5py.File(out_path, "w") as f:
            grp = f.create_group("matrix")
            grp.create_dataset("barcodes", data=barcodes)
            grp.create_dataset("data", data=matrix.data)
            grp.create_dataset("indices", data=matrix.indices)
            grp.create_dataset("indptr", data=matrix.indptr)
            grp.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int32))

            feat_grp = grp.create_group("features")
            feat_grp.create_dataset("id", data=feature_id)
            feat_grp.create_dataset("name", data=feature_name)
            feat_grp.create_dataset("feature_type", data=feature_type)
            feat_grp.create_dataset("genome", data=genome)
            feat_grp.create_dataset("_all_tag_keys", data=all_tag_keys)

        print(
            f"Saved Visium HD HDF5 to {out_path} (retained {retained_spots} of {total_spots} spots)."
        )

    @staticmethod
    def outlier(df, out_rate=0.05):
        tol = 0.01
        max_iter = 200
        step_scale = 0.25
        min_step = 0.001
        max_step = 0.05
        min_prob = 1e-4
        center_gap = 0.005

        counts = df["gene_count"].to_numpy()
        total = len(counts)
        if total == 0:
            return df.head(0)

        probs = [0.02, 0.98]
        data_outlier = df
        out_rate_epr = 0.0

        for _ in range(max_iter):
            quart_1 = np.quantile(counts, probs[0])
            quart_2 = np.quantile(counts, probs[1])
            iqr = np.quantile(counts, 0.75) - np.quantile(counts, 0.25)
            try:
                import scipy.stats as st

                skew = st.skew(counts)
            except Exception:
                skew = 0.0

            lower = quart_1 - 1.5 * np.exp(-4 * skew) * iqr
            upper = quart_2 + 1.5 * np.exp(3 * skew) * iqr
            mask = ~((counts > lower) & (counts < upper))
            data_outlier = df.filter(pl.Series(mask))
            out_rate_epr = len(data_outlier) / total

            if abs(out_rate_epr - out_rate) <= tol:
                break

            delta = float(np.clip((out_rate_epr - out_rate) * step_scale, -max_step, max_step))
            step = max(abs(delta), min_step)

            if delta >= 0:
                new_lower = max(probs[0] - step, min_prob)
                new_upper = min(probs[1] + step, 1.0 - min_prob)
            else:
                new_lower = min(probs[0] + step, 0.5 - center_gap)
                new_upper = max(probs[1] - step, 0.5 + center_gap)

            if new_lower >= new_upper - min_prob:
                break

            if new_lower == probs[0] and new_upper == probs[1]:
                break

            probs = [new_lower, new_upper]
        else:
            warnings.warn(
                f"Outlier detection stopped after {max_iter} iterations (target rate={out_rate:.3f}, "
                f"observed={out_rate_epr:.3f}). Results may deviate from the target rate.",
                RuntimeWarning,
            )

        return data_outlier


class ArtifactRemoveHD(Artifact_remove_HD):
    pass
