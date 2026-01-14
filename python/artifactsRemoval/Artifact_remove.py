import copy
import numpy as np
import scipy
import scipy.io
import scipy.stats
import polars as pl
import h5py
from scipy.sparse import csc_matrix

from artifactsRemoval.polars_neighbors import NeighborIndex


class Tissue_obj:
    def __init__(self, dir):
        self.dir = dir
        self.tissue_matrix = self.read_matrix()
        self.tissue_position = self.read_tissue_position()
        self.feature_list = self.read_features()
        self.barcode_list = self.read_barcodes()
        self.dict_barcode_to_coor = self.fn_barcode_to_coor()
        self.dict_coor_to_barcode = self.fn_coor_to_barcode()
        self.dict_barcode_to_column = self.fn_barcode_to_column()
        self.dict_feature_to_row = self.fn_feature_to_row()
        self._neighbor_index = NeighborIndex(self.tissue_position, default_mode="radius", k=6)

    def read_matrix(self):
        dataM = scipy.io.mmread(self.dir + "/filtered_feature_bc_matrix/matrix.mtx.gz")
        tissue_matrix = dataM.tocsr()
        return tissue_matrix

    def read_barcodes(self):
        return pl.read_csv(
            self.dir + "/filtered_feature_bc_matrix/barcodes.tsv.gz",
            separator="\t",
            has_header=False,
            new_columns=["barcode"],
        )

    def read_features(self):
        df = pl.read_csv(
            self.dir + "/filtered_feature_bc_matrix/features.tsv.gz",
            separator="\t",
            has_header=False,
        )
        rename = {}
        cols = df.columns
        names = ["feature_id", "feature_name", "feature_type"]
        for idx, name in enumerate(names):
            if idx < len(cols):
                rename[cols[idx]] = name
        if rename:
            df = df.rename(rename)
        return df

    def read_tissue_position(self):
        tissue_position = pl.read_csv(self.dir + "/spatial/tissue_positions.csv")
        if tissue_position.width < 6:
            raise ValueError("tissue_positions.csv must have at least 6 columns")
        new_cols = ["barcode", "in_tissue", "x_mtx", "y_mtx", "x_coor", "y_coor"]
        rename = {tissue_position.columns[i]: new_cols[i] for i in range(6)}
        tissue_position = tissue_position.rename(rename)
        return tissue_position

    def fn_barcode_to_coor(self):
        bars = self.tissue_position["barcode"].to_list()
        xs = self.tissue_position["x_mtx"].to_list()
        ys = self.tissue_position["y_mtx"].to_list()
        return {b: [x, y] for b, x, y in zip(bars, xs, ys)}

    def fn_coor_to_barcode(self):
        bars = self.tissue_position["barcode"].to_list()
        xs = self.tissue_position["x_mtx"].to_list()
        ys = self.tissue_position["y_mtx"].to_list()
        mapping = {}
        for b, x, y in zip(bars, xs, ys):
            mapping[f"{x} {y}"] = b
        return mapping

    def fn_barcode_to_column(self):
        bars = self.barcode_list["barcode"].to_list()
        return {b: idx for idx, b in enumerate(bars)}

    def fn_feature_to_row(self):
        if "feature_name" in self.feature_list.columns:
            genes = self.feature_list["feature_name"].to_list()
        else:
            genes = self.feature_list.select(pl.all()).to_numpy()[:, 1].tolist()
        return {g: idx for idx, g in enumerate(genes)}

    def neighbors_by_barcode(self, barcode, radius=None, k=None, include_center=False):
        return self._neighbor_index.neighbors_by_barcode(
            barcode, radius=radius, k=k, include_center=include_center
        )


class Artifact_detect:
    def __init__(self, param):
        if type(param) is str:
            self.Tissue = Tissue_obj(param)
        elif type(param) is Tissue_obj:
            self.Tissue = param
        else:
            ValueError("Either a directory or an instantiated Tissue_obj must be provided")

    def get_sum(self):
        gene_count = []
        x_mtx = []
        y_mtx = []
        barcode = []
        for row in self.Tissue.tissue_position.iter_rows(named=True):
            bar = self.Tissue.dict_coor_to_barcode.get(f"{row['x_mtx']} {row['y_mtx']}")
            i = self.Tissue.dict_barcode_to_column.get(bar)
            if i is None:
                continue
            x_mtx.append(row["x_mtx"])
            y_mtx.append(row["y_mtx"])
            barcode.append(bar)
            gene_count.append(self.Tissue.tissue_matrix[:, i].sum())

        d = {"barcode": barcode, "x_mtx": x_mtx, "y_mtx": y_mtx, "gene_count": gene_count}
        return pl.DataFrame(d)

    def get_neighbour_1(self, barcode):
        bars = self.Tissue.neighbors_by_barcode(barcode)
        return [self.Tissue.dict_barcode_to_column.get(b) for b in bars]

    def get_neighbour_1_bar(self, barcode):
        return self.Tissue.neighbors_by_barcode(barcode)

    def nh_variable(self, gene_sub_list):
        df = self.Tissue.barcode_list.clone()
        for gene in gene_sub_list:
            values = self.gen_table_by_gene(gene).get_column("gene_value")
            df = df.with_columns(values.alias(f"{gene}_nh"))
        return df

    def fn_adjmatrix(self):
        adjM = np.zeros(
            (self.Tissue.tissue_position.height, self.Tissue.tissue_position.height), dtype=int
        )
        for row in self.Tissue.barcode_list.iter_rows(named=True):
            temp1 = self.get_neighbour_1(row["barcode"])
            for i in temp1:
                if i is None:
                    continue
                adjM[self.Tissue.dict_barcode_to_column.get(row["barcode"]), i] = 1
        return adjM

    def BFS(self, point, full_list):
        adjmatrix = self.fn_adjmatrix()
        visited = []
        queue = []
        visited.append(point)
        queue.append(point)
        while queue:
            m = queue.pop(0)
            print(m, end=" ")
            temp_ngh = np.where(adjmatrix[m, :] == 1)[0].tolist()
            temp_sub_ngh = list(set(temp_ngh) & set(full_list))
            for neighbour in temp_sub_ngh:
                if neighbour not in visited:
                    visited.append(neighbour)
                    queue.append(neighbour)

        return visited

    def GET_cluster(self, barcode_list):
        point_list = [self.Tissue.dict_barcode_to_column.get(i) for i in barcode_list]
        clusters = []
        visited = []
        for point in point_list:
            if point in visited:
                continue
            cluster_temp = self.BFS(point, point_list)
            visited = visited + cluster_temp
            clusters.append(cluster_temp)

        return clusters

    def gen_table_by_gene(self, gene):
        gen_row = self.Tissue.dict_feature_to_row.get(gene)
        gen_table = []
        for row in self.Tissue.barcode_list.iter_rows(named=True):
            if self.Tissue.dict_barcode_to_coor.get(row["barcode"]) is None:
                gen_table.append(0)
                continue
            temp1 = self.get_neighbour_1(row["barcode"])
            value1 = []
            for i in temp1:
                if i is None:
                    value1.append(0)
                    continue
                value1.append(self.Tissue.tissue_matrix[gen_row, i].tolist())

            gen_table.append(sum(value1) / 6)
        return pl.DataFrame({"gene_value": gen_table})

    @staticmethod
    def outlier(df, out_rate=0.05):
        out_rate_epr = 0
        probs = [0.02, 0.98]
        data_outlier = df
        while abs(out_rate_epr - out_rate) > 0.01:
            if out_rate_epr - out_rate > 0:
                probs = [probs[0] - 0.0005, probs[1] + 0.0005]

            if out_rate_epr - out_rate < 0:
                probs = [probs[0] + 0.0005, probs[1] - 0.0005]

            counts = df["gene_count"].to_numpy()
            quart_1 = np.quantile(counts, probs[0])
            quart_2 = np.quantile(counts, probs[1])

            iqr = np.quantile(counts, 0.75) - np.quantile(counts, 0.25)
            medcp = scipy.stats.skew(counts)

            lower = quart_1 - 1.5 * np.exp(-4 * medcp) * iqr
            upper = quart_2 + 1.5 * np.exp(3 * medcp) * iqr
            temp_i = ~((counts > lower) & (counts < upper))

            data_outlier = df.filter(pl.Series(temp_i))
            out_rate_epr = len(data_outlier) / len(counts)

        return data_outlier

    def get_border(self):
        gene_count = []
        x_mtx = []
        y_mtx = []
        barcode = []
        for row in self.Tissue.tissue_position.iter_rows(named=True):
            if row["in_tissue"] == 0:
                continue
            if not (
                (row["x_mtx"] == 0)
                or (row["x_mtx"] == 77)
                or (row["y_mtx"] == 126)
                or (row["y_mtx"] == 127)
                or (row["y_mtx"] == 1)
                or (row["y_mtx"] == 0)
            ):
                continue
            bar = self.Tissue.dict_coor_to_barcode.get(f"{row['x_mtx']} {row['y_mtx']}")
            i = self.Tissue.dict_barcode_to_column.get(bar)
            if i is None:
                continue
            x_mtx.append(row["x_mtx"])
            y_mtx.append(row["y_mtx"])
            barcode.append(bar)
            gene_count.append(self.Tissue.tissue_matrix[:, i].sum())
        d = {"barcode": barcode, "x_mtx": x_mtx, "y_mtx": y_mtx, "gene_count": gene_count}
        return pl.DataFrame(d)

    def get_edge(self):
        tissue_position = self.Tissue.tissue_position
        zero_idx = tissue_position.filter(pl.col("in_tissue") == 0).get_column("barcode").to_list()
        covered_idx = tissue_position.filter(pl.col("in_tissue") == 1).get_column("barcode").to_list()

        neighbors = set()
        for barcode in zero_idx:
            neighbors.update(self.get_neighbour_1_bar(barcode))
        covered_tissue = set(covered_idx)

        neighbors.intersection_update(covered_tissue)

        gene_count = []
        x_mtx = []
        y_mtx = []
        bar = []
        for barcode in neighbors:
            coor = self.Tissue.dict_barcode_to_coor.get(barcode)
            i = self.Tissue.dict_barcode_to_column.get(barcode)
            if i is None or coor is None:
                continue
            x_mtx.append(coor[0])
            y_mtx.append(coor[1])
            bar.append(barcode)
            gene_count.append(self.Tissue.tissue_matrix[:, i].sum())

        d = {"barcode": bar, "x_mtx": x_mtx, "y_mtx": y_mtx, "gene_count": gene_count}
        return pl.DataFrame(d)


class Artifact_remove(Artifact_detect):
    def __init__(self, dir):
        super().__init__(dir)
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        self.tissue_depth = 0
        self.tissue_depth_df = self.fn_spot_dis_to_edge()
        self.remove_procedure = []
        self.removed_spots = {}
        self.remaining_depth = self.tissue_depth

    def fn_spot_dis_to_edge(self):
        bars = self.Tissue.tissue_position["barcode"].to_list()
        in_tissue = self.Tissue.tissue_position["in_tissue"].to_list()
        depth_map = {b: 0 for b in bars}

        tissue_position = self.Tissue.tissue_position
        zero_idx = tissue_position.filter(pl.col("in_tissue") == 0).get_column("barcode").to_list()
        covered_idx = tissue_position.filter(pl.col("in_tissue") == 1).get_column("barcode").to_list()
        neighbors = set()
        for barcode in zero_idx:
            neighbors.update(self.get_neighbour_1_bar(barcode))
        covered_tissue = set(covered_idx)
        df_border_spot = self.get_border()
        neighbors.update(df_border_spot["barcode"].to_list())

        neighbors.intersection_update(covered_tissue)
        neighbors.discard(None)
        for barcode in neighbors:
            depth_map[barcode] = 1
        remained = covered_tissue.difference(neighbors)

        level = 2
        while len(remained) != 0:
            ngh_iter = copy.deepcopy(neighbors)
            for barcode in ngh_iter:
                neighbors.update(self.get_neighbour_1_bar(barcode))
            neighbors.intersection_update(covered_tissue)
            neighbors.discard(None)
            added = neighbors.difference(ngh_iter)
            for barcode in added:
                depth_map[barcode] = level
            level = level + 1
            remained = covered_tissue.difference(neighbors)
        self.tissue_depth = level
        print(f"The depth of tissue is {level}")

        return pl.DataFrame(
            {
                "barcode": bars,
                "in_tissue": in_tissue,
                "spot_depth": [depth_map[b] for b in bars],
            }
        )

    def fn_spot_inclusion_condition(self):
        return self.Tissue.tissue_position.select(["barcode", "in_tissue"]).clone()

    def _set_inclusion_zero(self, barcodes):
        if not barcodes:
            return
        mask = pl.col("barcode").is_in(barcodes)
        self.spot_inclusion_condition = self.spot_inclusion_condition.with_columns(
            pl.when(mask).then(0).otherwise(pl.col("in_tissue")).alias("in_tissue")
        )

    def remove_edge(self, distance=1):
        self.remove_procedure.append("edge")
        depth_df = self.tissue_depth_df
        mask_df = self.spot_inclusion_condition.rename({"in_tissue": "current_in_tissue"})
        joined = depth_df.join(mask_df, on="barcode", how="left")
        removed_spot_n = joined.filter(
            (pl.col("spot_depth") <= distance) & (pl.col("current_in_tissue") == 1)
        ).height
        barcodes = depth_df.filter(pl.col("spot_depth") <= distance)["barcode"].to_list()
        self._set_inclusion_zero(barcodes)
        print(
            f"We removed {removed_spot_n} edge spots with distance to edge {distance} and less"
        )

    def remove_border(self):
        self.remove_procedure.append("border")
        df_border_spot = self.get_border()
        self._set_inclusion_zero(df_border_spot["barcode"].to_list())
        print(f"We removed {df_border_spot.height} border spots")

    def remove_malfunction(self):
        self.remove_procedure.append("malfunction")
        df_outlier_spot = Artifact_detect.outlier(self.get_sum())
        self._set_inclusion_zero(df_outlier_spot["barcode"].to_list())
        n_outlier = df_outlier_spot.height
        print(f"We removed {n_outlier} outlier spots")

    def simple_cleanser(self):
        self.remove_procedure = []
        self.removed_spots = {}
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        print("Back to the initial Tissue Sample")

    def review_removing(self):
        [print(i) for i in self.remove_procedure]

    def save(self, dir):
        df = self.Tissue.tissue_position.clone()
        df = df.drop("in_tissue")
        df = df.join(self.spot_inclusion_condition, on="barcode", how="left")
        df.write_csv(dir + "/tissue_positions.csv")
        print(f"Saved the tissue position file to {dir}/tissue_positions.csv")
        return df

    def save_hdf5(self, out_path):
        matrix = self.Tissue.tissue_matrix
        if not hasattr(matrix, "tocsc"):
            raise ValueError("tissue_matrix must be a scipy sparse matrix")
        matrix = matrix.tocsc()

        barcodes = self.Tissue.barcode_list.get_column("barcode").to_numpy().astype("S")
        features = self.Tissue.feature_list

        feature_id = features.get_column("feature_id").to_numpy().astype("S")
        feature_name = features.get_column("feature_name").to_numpy().astype("S")
        if "feature_type" in features.columns:
            feature_type = features.get_column("feature_type").to_numpy().astype("S")
        else:
            feature_type = np.array([b"Gene Expression"] * len(features))
        genome = np.array([b"GRCh38"] * len(features))
        all_tag_keys = np.array([b"feature_type"])

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

        print(f"Saved filtered matrix in Visium HDF5 format to {out_path}")
