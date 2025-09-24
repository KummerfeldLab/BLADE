
import os
import copy
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree


class TissueHD:
    """
    HD tissue object that reads expression from a Visium/Visium HD HDF5 matrix and
    spatial positions from CSV/Parquet. Provides fast neighbor queries via cKDTree.

    Expected HDF5 layout (10x-style):
      - group `matrix/`
        - datasets: barcodes, data, indices, indptr, shape
        - group `features/` with datasets: id, name, feature_type, genome, _all_tag_keys

    Position file expectations (flexible):
      - must contain one ID column among: barcode, cell_id, bin_id
      - must contain two coordinate columns among: (x_mtx,y_mtx), (x,y), (pxl_col_in_fullres,pxl_row_in_fullres)
      - optional `in_tissue` column; assumed 1 if missing
    """

    def __init__(self, data_dir=None, h5_path=None, positions_path=None,
                 neighbor_mode='k', k_neighbors=8, neighbor_radius=None):
        self.data_dir = data_dir
        self.h5_path = h5_path or self._auto_find_h5(data_dir)
        if self.h5_path is None:
            raise FileNotFoundError("HDF5 file not provided and not found in data_dir")

        # Load counts, barcodes, features
        (
            self.tissue_matrix,
            self.barcode_list,
            self.feature_list,
        ) = self._read_h5_matrix(self.h5_path)

        # Load and normalize positions
        self.tissue_position = self._read_positions(positions_path or self._auto_find_positions(data_dir))

        # Dictionaries for fast lookup
        self.dict_barcode_to_coor = self._barcode_to_coor()
        self.dict_coor_to_barcode = self._coor_to_barcode()
        self.dict_barcode_to_column = self._barcode_to_column()
        self.dict_feature_to_row = self._feature_to_row()

        # Neighbor configuration
        self._neighbor_mode = neighbor_mode  # 'k' or 'radius'
        self._k_neighbors = int(k_neighbors) if k_neighbors is not None else 8
        self._user_radius = float(neighbor_radius) if neighbor_radius is not None else None

        # Build KDTree for neighbor queries on available barcodes
        self._build_kdtree()

    # ----------------- IO helpers -----------------
    def _auto_find_h5(self, data_dir):
        if not data_dir:
            return None
        candidates = [
            os.path.join(data_dir, 'filtered_feature_cell_matrix.h5'),  # Visium HD binned output
            os.path.join(data_dir, 'filtered_feature_bc_matrix.h5'),    # Visium classic
        ]
        for p in candidates:
            if os.path.exists(p):
                return p
        return None

    def _auto_find_positions(self, data_dir):
        if not data_dir:
            return None
        candidates = [
            # Common Visium HD binned outputs
            os.path.join(data_dir, 'spatial', 'tissue_positions.parquet'),
            os.path.join(data_dir, 'spatial', 'tissue_positions.csv'),
            # Some pipelines write positions alongside matrix
            os.path.join(os.path.dirname(self.h5_path), 'spatial', 'tissue_positions.parquet'),
            os.path.join(os.path.dirname(self.h5_path), 'spatial', 'tissue_positions.csv'),
        ]
        for p in candidates:
            if p and os.path.exists(p):
                return p
        return None

    def _read_h5_matrix(self, h5_path):
        with h5py.File(h5_path, 'r') as f:
            grp = f['matrix']
            data = grp['data'][:]
            indices = grp['indices'][:]
            indptr = grp['indptr'][:]
            shape = tuple(grp['shape'][:].tolist())
            barcodes = grp['barcodes'][:]

            feat_grp = grp['features']
            feature_id = feat_grp['id'][:]
            feature_name = feat_grp['name'][:]
            feature_type = feat_grp['feature_type'][:] if 'feature_type' in feat_grp else None

        # Build CSC matrix (10x stores columns = barcodes)
        matrix = csc_matrix((data, indices, indptr), shape=shape)

        # Build DataFrames similar to non-HD code for compatibility
        barcode_df = pd.DataFrame([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in barcodes])

        feat_cols = {
            0: [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in feature_id],
            1: [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in feature_name],
        }
        if feature_type is not None:
            feat_cols[2] = [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in feature_type]
        feature_df = pd.DataFrame(feat_cols)

        return matrix, barcode_df, feature_df

    def _read_positions(self, positions_path):
        if positions_path is None:
            raise FileNotFoundError("Positions file not provided and not found automatically. Please pass positions_path.")

        # Parquet-only: explicitly require pyarrow to avoid silent engine errors
        if positions_path.endswith('.parquet'):
            try:
                import pyarrow  # noqa: F401
                df = pd.read_parquet(positions_path, engine='pyarrow')
            except ImportError as e:
                raise ImportError(
                    "Reading Parquet positions requires the 'pyarrow' package. "
                    "Please install pyarrow or convert your positions to CSV."
                ) from e
        else:
            df = pd.read_csv(positions_path)

        # Normalize columns
        id_col = None
        for cand in ['barcode', 'cell_id', 'bin_id', 'barcode_id']:
            if cand in df.columns:
                id_col = cand
                break
        if id_col is None:
            raise ValueError("Positions file must contain an ID column: barcode/cell_id/bin_id/barcode_id")

        # Coordinates
        x_col = None
        y_col = None
        for x_cand, y_cand in [
            ('x_mtx', 'y_mtx'),
            ('x', 'y'),
            ('pxl_col_in_fullres', 'pxl_row_in_fullres'),
        ]:
            if x_cand in df.columns and y_cand in df.columns:
                x_col, y_col = x_cand, y_cand
                break
        if x_col is None or y_col is None:
            # Fallback if columns look like generic indices
            raise ValueError("Positions file must contain coordinate columns, e.g., x_mtx/y_mtx or x/y")

        if 'in_tissue' not in df.columns:
            df['in_tissue'] = 1

        out = pd.DataFrame({
            'barcode': df[id_col].astype(str).values,
            'in_tissue': df['in_tissue'].astype(int).values,
            'x_mtx': df[x_col].values,
            'y_mtx': df[y_col].values,
            'x_coor': df[x_col].values,
            'y_coor': df[y_col].values,
        })
        return out

    # ----------------- dictionaries -----------------
    def _barcode_to_coor(self):
        return {row['barcode']: [row['x_mtx'], row['y_mtx']] for _, row in self.tissue_position.iterrows()}

    def _coor_to_barcode(self):
        d = {}
        for _, row in self.tissue_position.iterrows():
            key = f"{row['x_mtx']} {row['y_mtx']} "
            d[key] = row['barcode']
        return d

    def _barcode_to_column(self):
        d = {}
        # barcode_list is a one-column DataFrame like non-HD version
        for index, row in self.barcode_list.iterrows():
            d[row[0]] = index
        return d

    def _feature_to_row(self):
        d = {}
        # feature_list col 1 is gene name per non-HD convention
        for index, row in self.feature_list.iterrows():
            d[row[1]] = index
        return d

    # ----------------- KDTree neighbors -----------------
    def _build_kdtree(self):
        # Index ALL positions (including off-tissue) to support geometric neighbor queries
        pos_df = self.tissue_position[['barcode', 'x_mtx', 'y_mtx']].copy()
        self._pos_df = pos_df.reset_index(drop=True)
        self._coords = self._pos_df[['x_mtx', 'y_mtx']].to_numpy(dtype=float)
        if len(self._coords) == 0:
            self._tree = None
            self._pos_index = {}
            self._grid_radius = None
            return
        self._tree = cKDTree(self._coords)
        self._pos_index = {b: i for i, b in enumerate(self._pos_df['barcode'].values)}
        # Estimate neighbor radius from the first non-zero neighbor distance (median across points)
        dists, _ = self._tree.query(self._coords, k=2)
        # dists[:, 0] is zero (self), dists[:, 1] is nearest neighbor
        nn = dists[:, 1]
        nn = nn[np.isfinite(nn) & (nn > 0)]
        self._grid_radius = float(np.median(nn) * 1.1) if nn.size else None

    def neighbors_by_barcode(self, barcode, radius=None, k=None, include_center=False):
        if self._tree is None:
            return []
        idx = self._pos_index.get(barcode)
        if idx is None:
            return []
        # Determine neighbor query mode
        mode = self._neighbor_mode
        # Allow explicit overrides
        if k is not None:
            mode = 'k'
        if radius is not None:
            mode = 'radius'

        if mode == 'radius':
            if radius is None:
                radius = self._user_radius if self._user_radius is not None else self._grid_radius
            inds = self._tree.query_ball_point(self._coords[idx], r=radius)
            if not include_center:
                inds = [i for i in inds if i != idx]
        else:
            k_eff = k if k is not None else self._k_neighbors
            kq = k_eff + (0 if include_center else 1)
            dists, inds = self._tree.query(self._coords[idx], k=kq)
            if np.isscalar(inds):  # handle degenerate case
                inds = [int(inds)]
                dists = [float(dists)]
            if not include_center:
                # remove self (distance 0)
                paired = [(i, d) for i, d in zip(inds, dists) if i != idx]
                inds = [i for i, _ in paired][: k_eff]
        return self._pos_df.iloc[inds]['barcode'].tolist()


class ArtifactDetectHD(TissueHD):
    def __init__(self, data_dir=None, h5_path=None, positions_path=None, **kwargs):
        # Forward optional neighbor args (neighbor_mode, k_neighbors, neighbor_radius)
        super().__init__(data_dir=data_dir, h5_path=h5_path, positions_path=positions_path, **kwargs)

    # Keep public API similar to non-HD version
    def get_sum(self):
        # Vectorized column sums from sparse matrix
        col_sums = np.asarray(self.tissue_matrix.sum(axis=0)).ravel()
        # Align sums to tissue_position order via barcode mapping
        sums = []
        xs = []
        ys = []
        bars = []
        for _, row in self.tissue_position.iterrows():
            bar = row['barcode']
            idx = self.dict_barcode_to_column.get(bar)
            if idx is None:
                continue
            bars.append(bar)
            xs.append(row['x_mtx'])
            ys.append(row['y_mtx'])
            sums.append(col_sums[idx])
        return pd.DataFrame({"barcode": bars, "x_mtx": xs, "y_mtx": ys, "gene_count": sums})

    def seqdepth(self):
        df = pd.DataFrame(
            self.tissue_matrix.toarray().transpose(),
            index=self.barcode_list[0].values,
            columns=self.feature_list[1].values,
        )
        return df.sum(axis=1)

    # --- Neighbor API (compat with old names) ---
    def get_neighbour_1(self, barcode):
        bars = self.neighbors_by_barcode(barcode)
        return [self.dict_barcode_to_column.get(b) for b in bars]

    def get_neighbour_1_bar(self, barcode):
        return self.neighbors_by_barcode(barcode)

    # --- Graph helpers ---
    def fn_adjmatrix(self):
        n = self.tissue_position.shape[0]
        adjM = np.zeros((n, n), dtype=int)
        # position index maps barcode -> row index in tissue_position
        pos_index = {b: i for i, b in enumerate(self.tissue_position['barcode'].values)}
        for b in self.tissue_position['barcode'].values:
            i = pos_index[b]
            for nb in self.get_neighbour_1_bar(b):
                j = pos_index.get(nb)
                if j is not None:
                    adjM[i, j] = 1
        return adjM

    def BFS(self, point, full_list):
        adjmatrix = self.fn_adjmatrix()
        visited = []
        queue = []
        visited.append(point)
        queue.append(point)
        while queue:
            m = queue.pop(0)
            temp_ngh = np.where(adjmatrix[m, :] == 1)[0].tolist()
            temp_sub_ngh = list(set(temp_ngh) & set(full_list))
            for neighbour in temp_sub_ngh:
                if neighbour not in visited:
                    visited.append(neighbour)
                    queue.append(neighbour)
        return visited

    def GET_cluster(self, barcode_list):
        point_list = [self.dict_barcode_to_column.get(i) for i in barcode_list]
        clusters = []
        visited = []
        for point in point_list:
            if point in visited:
                continue
            cluster_temp = self.BFS(point, point_list)
            visited = visited + cluster_temp
            clusters.append(cluster_temp)
        return clusters

    # --- Region helpers (HD adaptations) ---
    def get_border(self, expected_neighbors=8, strategy='extreme', border_width=0.0, include_counts=True):
        """
        HD border definition options:
        - strategy='extreme': border are covered spots near the min/max coordinate edges
          within `border_width` (same units as coordinates). If border_width=0, it means
          exactly min/max edges.
        - strategy='degree': covered spots with fewer than `expected_neighbors` neighbors
          within KDTree grid neighborhood. Use this if extremes are unreliable.
        """
        covered_barcodes = set(self.spot_inclusion_condition.loc[self.spot_inclusion_condition['in_tissue'] == 1, 'barcode'])
        covered = self.tissue_position[self.tissue_position['barcode'].isin(covered_barcodes)]

        if len(covered) == 0:
            return pd.DataFrame({"barcode": [], "x_mtx": [], "y_mtx": [], "gene_count": []})

        bars = []
        x_mtx = []
        y_mtx = []
        gene_count = []

        if strategy == 'extreme':
            x_vals = covered['x_mtx'].to_numpy(dtype=float)
            y_vals = covered['y_mtx'].to_numpy(dtype=float)
            xmin, xmax = float(np.min(x_vals)), float(np.max(x_vals))
            ymin, ymax = float(np.min(y_vals)), float(np.max(y_vals))
            bw = float(border_width or 0.0)
            for _, row in covered.iterrows():
                x = float(row['x_mtx']); y = float(row['y_mtx'])
                is_border = (abs(x - xmin) <= bw) or (abs(x - xmax) <= bw) or (abs(y - ymin) <= bw) or (abs(y - ymax) <= bw)
                if not is_border:
                    continue
                b = row['barcode']
                i = self.dict_barcode_to_column.get(b)
                if i is None:
                    continue
                bars.append(b)
                x_mtx.append(row['x_mtx'])
                y_mtx.append(row['y_mtx'])
                if include_counts:
                    gene_count.append(self.tissue_matrix[:, i].sum())
        else:
            # Use radius-based neighbor counting for degree strategy
            rad = self._user_radius if self._user_radius is not None else self._grid_radius
            for _, row in covered.iterrows():
                b = row['barcode']
                ngh = self.neighbors_by_barcode(b, radius=rad)
                if len(ngh) < expected_neighbors:
                    i = self.dict_barcode_to_column.get(b)
                    if i is None:
                        continue
                    bars.append(b)
                    x_mtx.append(row['x_mtx'])
                    y_mtx.append(row['y_mtx'])
                    if include_counts:
                        gene_count.append(self.tissue_matrix[:, i].sum())
        if include_counts:
            return pd.DataFrame({"barcode": bars, "x_mtx": x_mtx, "y_mtx": y_mtx, "gene_count": gene_count})
        else:
            return pd.DataFrame({"barcode": bars})

    def get_edge(self, include_counts=True, radius=None):
        """
        Edge = covered spots that neighbor at least one off-tissue spot (current mask),
        robust for internal/irregular edges.

        - include_counts: whether to compute per-spot total counts.
        - radius: if provided, use radius-based neighbor query; otherwise use configured mode.
        """
        mask_df = self.spot_inclusion_condition
        zero_barcodes = set(mask_df.loc[mask_df['in_tissue'] == 0, 'barcode'])
        covered_barcodes = set(mask_df.loc[mask_df['in_tissue'] == 1, 'barcode'])

        neighbors = set()
        for zb in zero_barcodes:
            ngh = self.neighbors_by_barcode(zb, radius=radius) if radius is not None else self.get_neighbour_1_bar(zb)
            neighbors.update(ngh)
        neighbors.intersection_update(covered_barcodes)
        neighbors.discard(None)

        if not include_counts:
            return pd.DataFrame({"barcode": list(neighbors)})

        bars = []
        x_mtx = []
        y_mtx = []
        gene_count = []
        for b in neighbors:
            coor = self.dict_barcode_to_coor.get(b)
            i = self.dict_barcode_to_column.get(b)
            if i is None or coor is None:
                continue
            bars.append(b)
            x_mtx.append(coor[0])
            y_mtx.append(coor[1])
            gene_count.append(self.tissue_matrix[:, i].sum())
        return pd.DataFrame({"barcode": bars, "x_mtx": x_mtx, "y_mtx": y_mtx, "gene_count": gene_count})

    def get_inner(self, exclude_df):
        tissue_position = self.tissue_position.loc[self.tissue_position['in_tissue'] == 1]
        gene_count = []
        x_mtx = []
        y_mtx = []
        barcode = []
        for _, row in tissue_position.iterrows():
            bar = row['barcode']
            i = self.dict_barcode_to_column.get(bar)
            if i is None:
                continue
            x_mtx.append(row['x_mtx'])
            y_mtx.append(row['y_mtx'])
            barcode.append(bar)
            gene_count.append(self.tissue_matrix[:, i].sum())
        df = pd.DataFrame({"barcode": barcode, "x_mtx": x_mtx, "y_mtx": y_mtx, "gene_count": gene_count})
        return df[~df.barcode.isin(exclude_df.barcode)]


class Artifact_remove_HD(ArtifactDetectHD):
    def __init__(self, data_dir=None, h5_path=None, positions_path=None, compute_depth_on_init=False, max_depth_limit=4, **kwargs):
        super().__init__(data_dir=data_dir, h5_path=h5_path, positions_path=positions_path, **kwargs)
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        self.tissue_depth = 0
        self.tissue_depth_df = None
        if compute_depth_on_init:
            self.tissue_depth_df = self.fn_spot_dis_to_edge()
        self.remove_procedure = []
        self.removed_spots = {}
        self.remaining_depth = self.tissue_depth
        self.max_depth_limit = int(max_depth_limit)

    def fn_spot_inclusion_condition(self):
        return copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])

    def fn_spot_dis_to_edge(self):
        df = copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])
        df['spot_depth'] = [0] * (df.shape[0])

        # Edge spots: neighbors of off-tissue and border spots
        tissue_position = self.tissue_position
        zero_idx = tissue_position[tissue_position.in_tissue == 0].index.to_list()
        covered_idx = tissue_position[tissue_position.in_tissue == 1].index.to_list()
        neighbors = set()
        for idx in zero_idx:
            barcode = self.tissue_position.at[idx, 'barcode']
            neighbors.update(self.get_neighbour_1_bar(barcode))
        covered_tissue = set(self.tissue_position.loc[covered_idx, 'barcode'])

        # Border spots (HD definition)
        df_border_spot = self.get_border()
        neighbors.update(df_border_spot.barcode)

        neighbors.intersection_update(covered_tissue)
        neighbors.discard(None)
        for barcode in neighbors:
            df.loc[df.barcode == barcode, "spot_depth"] = 1
        remained = covered_tissue.difference(neighbors)

        # Grow layers inward by KDTree neighborhoods
        level = 2
        while len(remained) != 0:
            ngh_iter = copy.deepcopy(neighbors)
            for barcode in list(ngh_iter):
                neighbors.update(self.get_neighbour_1_bar(barcode))
            neighbors.intersection_update(covered_tissue)
            neighbors.discard(None)
            added = neighbors.difference(ngh_iter)
            for barcode in added:
                df.loc[df.barcode == barcode, "spot_depth"] = level
            level += 1
            remained = covered_tissue.difference(neighbors)

        self.tissue_depth = level
        print(f"The depth of tissue is {level}")
        return df

    # --- Removal operations ---
    def remove_edge(self, distance=1, strategy='extreme', expected_neighbors=8, seed_mode='zeros', border_width=None):
        """
        Remove spots within `distance` neighbor layers from the edge, without computing the full depth map.
        Capped by `self.max_depth_limit` to avoid deep expansions on huge datasets.

        strategy: for border detection, 'extreme' (coordinate edges with width) or 'degree'.
        expected_neighbors: used only when strategy='degree'.
        seed_mode: 'zeros' (default), 'border', or 'both'. 'zeros' strictly matches edge definition and captures internal/irregular edges.
        """
        self.remove_procedure.append("edge")

        if distance <= 0:
            print("Distance must be >= 1; nothing removed.")
            return

        if distance > self.max_depth_limit:
            print(f"Requested distance {distance} exceeds max limit {self.max_depth_limit}; clamping.")
            distance = self.max_depth_limit

        # Seeds: neighbors of off-tissue spots (optional) + border spots
        mask_df = self.spot_inclusion_condition
        covered_tissue = set(mask_df.loc[mask_df['in_tissue'] == 1, 'barcode'])

        frontier = set()

        # Optionally seed from zeros (current mask)
        if seed_mode in ('zeros', 'both'):
            zero_idx = mask_df[mask_df['in_tissue'] == 0].index.to_list()
            margin = (self._user_radius if self._user_radius is not None else self._grid_radius) or 0.0
            band = margin * 2.0 if margin > 0 else None
            if band is not None and len(covered_tissue) > 0:
                cov_pos = self.tissue_position[self.tissue_position['barcode'].isin(covered_tissue)]
                xmin = float(np.min(cov_pos['x_mtx'])); xmax = float(np.max(cov_pos['x_mtx']))
                ymin = float(np.min(cov_pos['y_mtx'])); ymax = float(np.max(cov_pos['y_mtx']))
                zpos = self.tissue_position.loc[zero_idx, ['barcode', 'x_mtx', 'y_mtx']]
                zmask = (zpos['x_mtx'] >= xmin - band) & (zpos['x_mtx'] <= xmax + band) & \
                        (zpos['y_mtx'] >= ymin - band) & (zpos['y_mtx'] <= ymax + band)
                zero_barcodes = zpos.loc[zmask, 'barcode'].tolist()
            else:
                zero_barcodes = self.tissue_position.loc[zero_idx, 'barcode'].tolist()
            for barcode in zero_barcodes:
                frontier.update(self.get_neighbour_1_bar(barcode))

        # Optionally seed from border barcodes
        if seed_mode in ('border', 'both'):
            if strategy == 'extreme':
                bw = border_width if border_width is not None else ((self._user_radius or self._grid_radius) or 0.0)
                df_border_spot = self.get_border(strategy='extreme', border_width=bw, include_counts=False)
            else:
                df_border_spot = self.get_border(strategy='degree', expected_neighbors=expected_neighbors, include_counts=False)
            frontier.update(df_border_spot.barcode)

        # Keep only covered
        frontier.intersection_update(covered_tissue)
        frontier.discard(None)

        # Expand up to (distance-1) more layers
        visited = set(frontier)
        for _ in range(1, distance):
            next_frontier = set()
            for b in list(frontier):
                next_frontier.update(self.get_neighbour_1_bar(b))
            next_frontier.intersection_update(covered_tissue)
            next_frontier.discard(None)
            next_frontier.difference_update(visited)
            if not next_frontier:
                break
            visited.update(next_frontier)
            frontier = next_frontier

        # Remove visited from inclusion_condition
        mask = self.spot_inclusion_condition['barcode'].isin(visited)
        before = int(self.spot_inclusion_condition.loc[mask, 'in_tissue'].sum())
        self.spot_inclusion_condition.loc[mask, 'in_tissue'] = 0
        removed_spot_n = before  # spots that were previously in_tissue
        print(f"We removed {removed_spot_n} edge spots within {distance} layer(s)")

    def remove_border(self, distance=None, strategy='extreme', expected_neighbors=8):
        self.remove_procedure.append("border")
        if strategy == 'extreme':
            df_border_spot = self.get_border(strategy='extreme', border_width=(distance or 0.0))
        else:
            df_border_spot = self.get_border(strategy='degree', expected_neighbors=expected_neighbors)
        for barcode in df_border_spot.barcode:
            self.spot_inclusion_condition.loc[self.spot_inclusion_condition.barcode == barcode, 'in_tissue'] = 0
        msg = f"We removed {len(df_border_spot.barcode)} border spots"
        if strategy == 'extreme':
            msg += f" within border width {distance or 0.0}"
        else:
            msg += f" (degree< {expected_neighbors})"
        print(msg)

    def remove_malfunction(self, use_clusters=True, min_cluster_size=15, radius=None, radius_scale=1.0):
        """
        Remove malfunction (outlier) spots.

        By default, performs cluster-aware removal (Option C) and only removes clusters
        of outliers with size >= min_cluster_size (default 15). Set use_clusters=False
        to use the original, purely distributional removal without spatial clustering.
        """
        if use_clusters:
            return self.remove_malfunction_clustered(min_cluster_size=min_cluster_size,
                                                     radius=radius,
                                                     radius_scale=radius_scale)
        # Fallback: original non-clustered behavior
        self.remove_procedure.append("malfunction")
        df_outlier_spot = Artifact_remove_HD.outlier(self.get_sum())
        for barcode in df_outlier_spot.barcode:
            self.spot_inclusion_condition.loc[self.spot_inclusion_condition.barcode == barcode, 'in_tissue'] = 0
        n_outlier = len(df_outlier_spot.barcode)
        print(f"We removed {n_outlier} outlier spots (non-clustered)")

    def remove_malfunction_clustered(self, min_cluster_size=15, radius=None, radius_scale=1.0):
        """
        Cluster-aware malfunction removal (Option C):
        - Detect outliers via count distribution
        - Build KDTree over outlier coordinates only
        - Find connected components using radius-based adjacency
        - Remove only clusters with size >= min_cluster_size

        Args:
            min_cluster_size: minimum number of outlier spots in a cluster to remove (default 15)
            radius: neighbor radius to define connectivity; defaults to grid/user radius
            radius_scale: multiply the base radius by this factor (e.g., 1.0–1.2)
        """
        self.remove_procedure.append("malfunction_clustered")

        # 1) Outlier candidates (barcodes) intersected with currently included spots
        df_out = Artifact_remove_HD.outlier(self.get_sum())
        include_mask = set(self.spot_inclusion_condition.loc[self.spot_inclusion_condition['in_tissue'] == 1, 'barcode'])
        out_barcodes = [b for b in df_out.barcode if b in include_mask]
        if len(out_barcodes) == 0:
            print("No outlier spots among currently included spots.")
            return

        # 2) Coordinates array for subset KDTree
        coords = []
        valid_barcodes = []
        for b in out_barcodes:
            coor = self.dict_barcode_to_coor.get(b)
            if coor is None:
                continue
            valid_barcodes.append(b)
            coords.append(coor)
        if len(coords) == 0:
            print("No coordinates found for outlier spots.")
            return

        import numpy as np
        from scipy.spatial import cKDTree
        pts = np.asarray(coords, dtype=float)
        tree = cKDTree(pts)

        # 3) Determine radius
        base_r = radius
        if base_r is None:
            base_r = self._user_radius if self._user_radius is not None else self._grid_radius
        if base_r is None:
            # Fallback: estimate from subset 1-NN
            dists, _ = tree.query(pts, k=2)
            nn = dists[:, 1]
            nn = nn[np.isfinite(nn) & (nn > 0)]
            base_r = float(np.median(nn)) if nn.size else 1.0
        rad = float(base_r) * float(radius_scale)

        # 4) Union-Find over neighbors within radius
        n = len(valid_barcodes)
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

        # Query neighbors for each point and union
        for i in range(n):
            inds = tree.query_ball_point(pts[i], r=rad)
            for j in inds:
                if j != i:
                    union(i, j)

        # 5) Gather clusters
        clusters = {}
        for i in range(n):
            r = find(i)
            clusters.setdefault(r, []).append(valid_barcodes[i])

        # 6) Remove clusters meeting size threshold
        removed = 0
        for members in clusters.values():
            if len(members) >= int(min_cluster_size):
                mask = self.spot_inclusion_condition['barcode'].isin(members)
                removed += int(self.spot_inclusion_condition.loc[mask, 'in_tissue'].sum())
                self.spot_inclusion_condition.loc[mask, 'in_tissue'] = 0

        print(f"Clustered malfunction removal: considered {n} outliers, removed {removed} spots in clusters >= {min_cluster_size} (radius={rad:.3f}).")

    # --- Utility / persistence ---
    def simple_cleanser(self):
        self.remove_procedure = []
        self.removed_spots = {}
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        print("Back to the initial Tissue Sample")

    def review_removing(self):
        [print(i) for i in self.remove_procedure]

    def save(self, out_dir):
        df = copy.deepcopy(self.tissue_position)
        df.in_tissue = self.spot_inclusion_condition.in_tissue
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, 'tissue_positions.csv')
        df.to_csv(out_path, index=False)
        print(f"Saved the tissue position file to {out_path}")
        return df

    @staticmethod
    def outlier(df, out_rate=0.05):
        out_rate_epr = 0
        probs = [.02, .98]
        data_outlier = df
        # Iterate to calibrate whiskers to target rate
        while abs(out_rate_epr - (out_rate)) > 0.01:
            if (out_rate_epr - (out_rate) > 0):
                probs = [probs[0] - 0.0005, probs[1] + 0.0005]
            if (out_rate_epr - (out_rate) < 0):
                probs = [probs[0] + 0.0005, probs[1] - 0.0005]

            quart_1 = np.quantile(df.gene_count, probs[0])
            quart_2 = np.quantile(df.gene_count, probs[1])
            IQR = np.quantile(df.gene_count, 0.75) - np.quantile(df.gene_count, 0.25)

            # Robust skew-aware whiskers
            # If scipy.stats.skew is needed, import locally to avoid heavy global deps
            try:
                import scipy.stats as st
                medcp = st.skew(df.gene_count)
            except Exception:
                medcp = 0.0

            Lower = quart_1 - 1.5 * np.exp(-4 * medcp) * IQR
            Upper = quart_2 + 1.5 * np.exp(3 * medcp) * IQR
            temp_i = np.logical_not((df.gene_count > Lower) & (df.gene_count < Upper))

            data_outlier = df[temp_i]
            out_rate_epr = len(data_outlier) / len(df.gene_count)

        return data_outlier

    def save_hdf5(self, out_path):
        """
        Save the filtered matrix (keeping spots with in_tissue==1), barcodes, and features
        in a 10x-style HDF5 file at `out_path`.
        """
        # Build inclusion map from current state
        include_map = dict(self.spot_inclusion_condition[['barcode', 'in_tissue']].values)

        keep_indices = []
        keep_barcodes = []
        for idx, row in self.barcode_list.iterrows():
            b = str(row[0])
            if include_map.get(b, 0) == 1:
                keep_indices.append(idx)
                keep_barcodes.append(b)

        # Filter matrix columns
        matrix = self.tissue_matrix.tocsc() if hasattr(self.tissue_matrix, 'tocsc') else self.tissue_matrix
        if keep_indices:
            filt_matrix = matrix[:, keep_indices]
        else:
            # Empty matrix with same number of rows
            from scipy.sparse import csc_matrix as _csc
            filt_matrix = _csc((matrix.shape[0], 0))

        # Prepare datasets
        barcodes_arr = np.array([x.encode() if not isinstance(x, (bytes, np.bytes_)) else x for x in keep_barcodes])

        # Features
        feat_id = self.feature_list.iloc[:, 0].astype(str).apply(lambda x: x.encode()).values
        feat_name = self.feature_list.iloc[:, 1].astype(str).apply(lambda x: x.encode()).values
        if self.feature_list.shape[1] > 2:
            feat_type = self.feature_list.iloc[:, 2].astype(str).apply(lambda x: x.encode()).values
        else:
            feat_type = np.array([b'Gene Expression'] * len(self.feature_list))

        # Try to carry genome and _all_tag_keys from source HDF5 if present
        genome = None
        all_tag_keys = None
        try:
            with h5py.File(self.h5_path, 'r') as f:
                if 'features' in f['matrix'] and 'genome' in f['matrix']['features']:
                    genome = f['matrix']['features']['genome'][:]
                if 'features' in f['matrix'] and '_all_tag_keys' in f['matrix']['features']:
                    all_tag_keys = f['matrix']['features']['_all_tag_keys'][:]
        except Exception:
            pass
        if genome is None:
            genome = np.array([b'unknown'] * len(self.feature_list))
        if all_tag_keys is None:
            all_tag_keys = np.array([b'feature_type'])

        # Write HDF5
        os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
        with h5py.File(out_path, 'w') as f:
            grp = f.create_group('matrix')
            grp.create_dataset('barcodes', data=barcodes_arr)
            grp.create_dataset('data', data=filt_matrix.data)
            grp.create_dataset('indices', data=filt_matrix.indices)
            grp.create_dataset('indptr', data=filt_matrix.indptr)
            grp.create_dataset('shape', data=np.array(filt_matrix.shape, dtype=np.int32))

            feat_grp = grp.create_group('features')
            feat_grp.create_dataset('id', data=feat_id)
            feat_grp.create_dataset('name', data=feat_name)
            feat_grp.create_dataset('feature_type', data=feat_type)
            feat_grp.create_dataset('genome', data=genome)
            feat_grp.create_dataset('_all_tag_keys', data=all_tag_keys)

        print(f"Saved filtered matrix in HDF5 format to {out_path}")


# Backwards-compatible alias with the original stub class name
class ArtifactRemoveHD(Artifact_remove_HD):
    pass
