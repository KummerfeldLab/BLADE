import os
import copy
import warnings
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree


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

        self._build_kdtree()

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

        barcode_df = pd.DataFrame(_decode(barcodes))
        feat_cols = {0: _decode(feature_id), 1: _decode(feature_name)}
        if feature_type is not None:
            feat_cols[2] = _decode(feature_type)
        feature_df = pd.DataFrame(feat_cols)
        return matrix, barcode_df, feature_df

    def _read_positions(self, positions_path):
        if positions_path is None:
            raise FileNotFoundError(
                "Positions file not provided and automatic discovery failed. "
                "Provide `positions_path` pointing to tissue_positions.parquet or CSV."
            )

        if positions_path.endswith(".parquet"):
            try:
                import pyarrow  # noqa: F401
                df = pd.read_parquet(positions_path, engine="pyarrow")
            except ImportError as exc:
                raise ImportError(
                    "Reading Parquet positions requires the 'pyarrow' package. "
                    "Install pyarrow or provide a CSV positions file."
                ) from exc
        else:
            df = pd.read_csv(positions_path)

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
            df["in_tissue"] = 1

        out = pd.DataFrame(
            {
                "barcode": df[id_col].astype(str).values,
                "in_tissue": df["in_tissue"].astype(int).values,
                "x_mtx": df[x_col].astype(float).values,
                "y_mtx": df[y_col].astype(float).values,
                "x_coor": df[x_col].astype(float).values,
                "y_coor": df[y_col].astype(float).values,
            }
        )
        return out

    def _barcode_to_coor(self):
        return {row["barcode"]: [row["x_mtx"], row["y_mtx"]] for _, row in self.tissue_position.iterrows()}

    def _coor_to_barcode(self):
        mapping = {}
        for _, row in self.tissue_position.iterrows():
            key = f"{row['x_mtx']} {row['y_mtx']}"
            mapping[key] = row["barcode"]
        return mapping

    def _barcode_to_column(self):
        mapping = {}
        for idx, row in self.barcode_list.iterrows():
            mapping[row[0]] = idx
        return mapping

    def _feature_to_row(self):
        mapping = {}
        for idx, row in self.feature_list.iterrows():
            mapping[row[1]] = idx
        return mapping

    def _build_kdtree(self):
        pos_df = self.tissue_position[["barcode", "x_mtx", "y_mtx"]].copy()
        self._pos_df = pos_df.reset_index(drop=True)
        self._coords = self._pos_df[["x_mtx", "y_mtx"]].to_numpy(dtype=float)
        if len(self._coords) == 0:
            self._tree = None
            self._pos_index = {}
            self._grid_radius = None
            return
        self._tree = cKDTree(self._coords)
        self._pos_index = {b: i for i, b in enumerate(self._pos_df["barcode"].values)}

        dists, _ = self._tree.query(self._coords, k=2)
        nn = dists[:, 1]
        nn = nn[np.isfinite(nn) & (nn > 0)]
        self._grid_radius = float(np.median(nn) * 1.1) if nn.size else None

    def neighbors_by_barcode(self, barcode, radius=None, k=None, include_center=False):
        if self._tree is None:
            return []
        idx = self._pos_index.get(barcode)
        if idx is None:
            return []

        mode = self._neighbor_mode
        if radius is not None:
            mode = "radius"
        elif k is not None:
            mode = "k"

        if mode == "radius":
            rad = radius
            if rad is None:
                rad = self._user_radius if self._user_radius is not None else self._grid_radius
            inds = self._tree.query_ball_point(self._coords[idx], r=rad)
            if not include_center:
                inds = [i for i in inds if i != idx]
        else:
            k_eff = k if k is not None else self._k_neighbors
            k_query = k_eff + (0 if include_center else 1)
            dists, inds = self._tree.query(self._coords[idx], k=k_query)
            if np.isscalar(inds):
                inds = [int(inds)]
                dists = [float(dists)]
            if not include_center:
                pairs = [(i, d) for i, d in zip(inds, dists) if i != idx]
                inds = [i for i, _ in pairs][:k_eff]

        return self._pos_df.iloc[inds]["barcode"].tolist()


class ArtifactDetectHD(TissueHD):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_sum(self):
        col_sums = np.asarray(self.tissue_matrix.sum(axis=0)).ravel()
        bars = []
        xs = []
        ys = []
        sums = []
        for _, row in self.tissue_position.iterrows():
            b = row["barcode"]
            idx = self.dict_barcode_to_column.get(b)
            if idx is None:
                continue
            bars.append(b)
            xs.append(row["x_mtx"])
            ys.append(row["y_mtx"])
            sums.append(col_sums[idx])
        return pd.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": sums})

    def seqdepth(self):
        df = pd.DataFrame(
            self.tissue_matrix.toarray().transpose(),
            index=self.barcode_list[0].values,
            columns=self.feature_list[1].values,
        )
        return df.sum(axis=1)

    def get_neighbour_1(self, barcode):
        bars = self.neighbors_by_barcode(barcode)
        return [self.dict_barcode_to_column.get(b) for b in bars]

    def get_neighbour_1_bar(self, barcode):
        return self.neighbors_by_barcode(barcode)

    def get_border(self, expected_neighbors=8, strategy="extreme", border_width=0.0, include_counts=True):
        mask_df = getattr(self, "spot_inclusion_condition", self.tissue_position[["barcode", "in_tissue"]])
        covered_barcodes = set(mask_df.loc[mask_df["in_tissue"] == 1, "barcode"])
        covered = self.tissue_position[self.tissue_position["barcode"].isin(covered_barcodes)]
        if len(covered) == 0:
            return pd.DataFrame({"barcode": []})

        bars = []
        xs = []
        ys = []
        gene_count = []

        if strategy == "extreme":
            x_vals = covered["x_mtx"].to_numpy(dtype=float)
            y_vals = covered["y_mtx"].to_numpy(dtype=float)
            xmin, xmax = np.min(x_vals), np.max(x_vals)
            ymin, ymax = np.min(y_vals), np.max(y_vals)
            bw = float(border_width or 0.0)
            for _, row in covered.iterrows():
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
            rad = self._user_radius if self._user_radius is not None else self._grid_radius
            for _, row in covered.iterrows():
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
            return pd.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": gene_count})
        return pd.DataFrame({"barcode": bars})

    def get_edge(self, include_counts=True, radius=None):
        mask_df = getattr(self, "spot_inclusion_condition", self.tissue_position[["barcode", "in_tissue"]])
        zero_barcodes = set(mask_df.loc[mask_df["in_tissue"] == 0, "barcode"])
        covered_barcodes = set(mask_df.loc[mask_df["in_tissue"] == 1, "barcode"])

        neighbors = set()
        for zb in zero_barcodes:
            ngh = self.neighbors_by_barcode(zb, radius=radius) if radius is not None else self.get_neighbour_1_bar(zb)
            neighbors.update(ngh)
        neighbors.intersection_update(covered_barcodes)
        neighbors.discard(None)

        if not include_counts:
            return pd.DataFrame({"barcode": list(neighbors)})

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
        return pd.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": gene_count})


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
        return copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])

    def fn_spot_dis_to_edge(self):
        df = copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])
        df["spot_depth"] = 0

        mask_df = self.spot_inclusion_condition
        zero_barcodes = set(mask_df.loc[mask_df["in_tissue"] == 0, "barcode"])
        covered_barcodes = set(mask_df.loc[mask_df["in_tissue"] == 1, "barcode"])

        frontier = set()
        for zb in zero_barcodes:
            frontier.update(self.get_neighbour_1_bar(zb))

        border_barcodes = self.get_border(include_counts=False)["barcode"].tolist()
        frontier.update(border_barcodes)
        frontier.intersection_update(covered_barcodes)
        frontier.discard(None)

        for b in frontier:
            df.loc[df.barcode == b, "spot_depth"] = 1
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
                df.loc[df.barcode == b, "spot_depth"] = level
            visited.update(next_frontier)
            frontier = next_frontier
            remaining = covered_barcodes - visited
            level += 1

        self.tissue_depth = level
        self.tissue_depth_df = df
        return df

    def remove_border(self, distance=None, strategy="extreme", expected_neighbors=8):
        self.remove_procedure.append("border")
        if strategy == "extreme":
            df_border = self.get_border(strategy="extreme", border_width=(distance or 0.0))
        else:
            df_border = self.get_border(strategy="degree", expected_neighbors=expected_neighbors)
        for b in df_border["barcode"]:
            self.spot_inclusion_condition.loc[self.spot_inclusion_condition["barcode"] == b, "in_tissue"] = 0
        msg = f"We removed {len(df_border['barcode'])} border spots"
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
        covered_barcodes = set(mask_df.loc[mask_df["in_tissue"] == 1, "barcode"])

        frontier = set()
        if seed_mode in ("zeros", "both"):
            zero_barcodes = mask_df.loc[mask_df["in_tissue"] == 0, "barcode"].tolist()
            margin = (self._user_radius if self._user_radius is not None else self._grid_radius) or 0.0
            band = margin * 2.0 if margin > 0 else None
            if band is not None and covered_barcodes:
                cov_pos = self.tissue_position[self.tissue_position["barcode"].isin(covered_barcodes)]
                xmin, xmax = cov_pos["x_mtx"].min(), cov_pos["x_mtx"].max()
                ymin, ymax = cov_pos["y_mtx"].min(), cov_pos["y_mtx"].max()
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
                bw = border_width if border_width is not None else ((self._user_radius or self._grid_radius) or 0.0)
                border_df = self.get_border(strategy="extreme", border_width=bw, include_counts=False)
            else:
                border_df = self.get_border(strategy="degree", expected_neighbors=expected_neighbors, include_counts=False)
            frontier.update(border_df["barcode"].tolist())

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

        mask = self.spot_inclusion_condition["barcode"].isin(visited)
        removed = int(self.spot_inclusion_condition.loc[mask, "in_tissue"].sum())
        self.spot_inclusion_condition.loc[mask, "in_tissue"] = 0
        print(f"We removed {removed} edge spots within {distance} layer(s)")

    def remove_malfunction(self, use_clusters=True, min_cluster_size=15, radius=None, radius_scale=1.0):
        """
        Remove malfunction (outlier) spots.

        By default, performs cluster-aware removal and only removes clusters with size
        >= min_cluster_size. Set use_clusters=False to fall back to the original
        non-spatial removal.
        """
        if use_clusters:
            return self.remove_malfunction_clustered(
                min_cluster_size=min_cluster_size, radius=radius, radius_scale=radius_scale
            )

        self.remove_procedure.append("malfunction")
        df_outlier = Artifact_remove_HD.outlier(self.get_sum())
        for b in df_outlier["barcode"]:
            self.spot_inclusion_condition.loc[self.spot_inclusion_condition["barcode"] == b, "in_tissue"] = 0
        print(f"We removed {len(df_outlier['barcode'])} outlier spots (non-clustered)")

    def remove_malfunction_clustered(self, min_cluster_size=15, radius=None, radius_scale=1.0):
        self.remove_procedure.append("malfunction_clustered")
        df_out = Artifact_remove_HD.outlier(self.get_sum())
        include_mask = set(self.spot_inclusion_condition.loc[self.spot_inclusion_condition["in_tissue"] == 1, "barcode"])
        outliers = [b for b in df_out["barcode"] if b in include_mask]
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
        tree = cKDTree(pts)
        base_r = radius
        if base_r is None:
            base_r = self._user_radius if self._user_radius is not None else self._grid_radius
        if base_r is None:
            dists, _ = tree.query(pts, k=2)
            nn = dists[:, 1]
            nn = nn[np.isfinite(nn) & (nn > 0)]
            base_r = float(np.median(nn)) if nn.size else 1.0
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

        for i in range(n):
            inds = tree.query_ball_point(pts[i], r=rad)
            for j in inds:
                if j != i:
                    union(i, j)

        clusters = {}
        for i in range(n):
            clusters.setdefault(find(i), []).append(kept_barcodes[i])

        removed = 0
        for members in clusters.values():
            if len(members) >= int(min_cluster_size):
                mask = self.spot_inclusion_condition["barcode"].isin(members)
                removed += int(self.spot_inclusion_condition.loc[mask, "in_tissue"].sum())
                self.spot_inclusion_condition.loc[mask, "in_tissue"] = 0

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
        df = copy.deepcopy(self.tissue_position)
        df["in_tissue"] = self.spot_inclusion_condition["in_tissue"].values
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, "tissue_positions.csv")
        df.to_csv(out_path, index=False)
        print(f"Saved the tissue position file to {out_path}")
        return df

    def save_hdf5(self, out_path, filter_in_tissue=True):
        """Save the HD matrix in 10x-style HDF5 format, optionally filtering to in-tissue spots.

        Parameters
        ----------
        out_path : str
            Destination path for the output .h5 file.
        filter_in_tissue : bool, optional
            If True (default), keep only spots currently marked as in-tissue.
        """
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

        barcodes_series = self.barcode_list.iloc[:, 0].astype(str).reset_index(drop=True)
        if filter_in_tissue and hasattr(self, "spot_inclusion_condition"):
            inclusion = (
                self.spot_inclusion_condition.set_index("barcode")["in_tissue"]
                .reindex(barcodes_series)
                .fillna(0)
                .astype(bool)
            )
            keep_mask = inclusion.to_numpy(dtype=bool)
            matrix = matrix[:, keep_mask]
            matrix = matrix.tocsc()
            barcodes_series = barcodes_series[keep_mask]
        retained_spots = matrix.shape[1]
        barcodes = barcodes_series.to_numpy(dtype="S")

        features = self.feature_list
        feature_id = features.iloc[:, 0].astype(str).to_numpy(dtype="S")
        feature_name = features.iloc[:, 1].astype(str).to_numpy(dtype="S")
        if features.shape[1] > 2:
            feature_type = features.iloc[:, 2].astype(str).to_numpy(dtype="S")
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
            return df.iloc[0:0]

        probs = [0.02, 0.98]
        data_outlier = df
        out_rate_epr = 0.0

        for iteration in range(max_iter):
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
            data_outlier = df[mask]
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
