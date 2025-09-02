
import pandas as pd
import numpy as np
import scipy
import scipy.io
import scipy.stats
import copy
import h5py
from scipy.sparse import csc_matrix

# Import the efficient edge detection function
from utils.spatial_utils import find_edge_spots




class Tissue_obj:

    def __init__(self, dir):
        self.dir = dir
        self.tissue_matrix = self.read_matrix()
        self.tissue_position = self.read_tissue_position()
        self.feature_list = self.read_features()
        self.barcode_list = self.read_barcodes()
        ########### Here are useful dictionaries
        self.dict_barcode_to_coor = self.fn_barcode_to_coor()
        self.dict_coor_to_barcode = self.fn_coor_to_barcode()
        self.dict_barcode_to_column = self.fn_barcode_to_column()
        self.dict_feature_to_row = self.fn_feature_to_row()
        #self.dict_gene_name_to_id = self.fn_gene_name_to_id()
        ########### 
    
    #def read_gene_list(self):
    #   return pd.read_csv(self.dir + 'data/' + self.gene_list, names=['gene_name', 'gene_id'], header=0)

    def read_matrix(self):
        dataM = scipy.io.mmread(self.dir + "/filtered_feature_bc_matrix/matrix.mtx.gz") # here need to revise
        tissue_matrix = dataM.tocsr()
        return tissue_matrix

    def read_barcodes(self):
        barcode_list = pd.read_csv(self.dir + '/filtered_feature_bc_matrix/barcodes.tsv.gz',sep = '\t',compression='gzip', header=None)
        return barcode_list

    def read_features(self):
        feature_list = pd.read_csv(self.dir + '/filtered_feature_bc_matrix/features.tsv.gz',sep = '\t',compression='gzip', header=None)
        return feature_list


    def read_tissue_position(self):
        tissue_position = pd.read_csv(self.dir + "/spatial/tissue_positions.csv")
        tissue_position.columns = ['barcode', 'in_tissue', 'x_mtx', 'y_mtx', 'x_coor', 'y_coor']
        #tissue_position = tissue_position.loc[tissue_position['in_tissue'] == 1]
        #tissue_position = self.tissue_position.loc[self.tissue_position['barcode'].isin(self.barcode_list)]
        return tissue_position


    def fn_barcode_to_coor(self):
        dict = {}
        for index, row in self.tissue_position.iterrows():
            dict[row['barcode']] = [row['x_mtx'], row['y_mtx']]
        return dict


# Move Artifact_detect to module level
class Artifact_detect(Tissue_obj):
    def __init__(self, dir):
        super().__init__(dir)
        # Prepare coordinates and tissue mask for spatial edge detection
        self.coords_xy = self.tissue_position[['x_mtx', 'y_mtx']].values.astype(float)
        self.in_tissue = self.tissue_position['in_tissue'].astype(bool).values

    def get_sum(self):
        # Efficiently compute gene count for each spot
        gene_count = np.array(self.tissue_matrix.sum(axis=0)).flatten()
        df = pd.DataFrame({
            "barcode": self.barcode_list[0],
            "x_mtx": self.tissue_position['x_mtx'],
            "y_mtx": self.tissue_position['y_mtx'],
            "gene_count": gene_count
        })
        return df

    def get_edge(self, radius_factor=1.35):
        # Use spatial_utils KD-tree based edge detection
        edge_mask = find_edge_spots(self.coords_xy, self.in_tissue, radius_factor=radius_factor)
        edge_barcodes = self.tissue_position.loc[edge_mask, 'barcode'].values
        edge_x = self.tissue_position.loc[edge_mask, 'x_mtx'].values
        edge_y = self.tissue_position.loc[edge_mask, 'y_mtx'].values
        gene_count = np.array(self.tissue_matrix.sum(axis=0)).flatten()[edge_mask]
        df = pd.DataFrame({
            "barcode": edge_barcodes,
            "x_mtx": edge_x,
            "y_mtx": edge_y,
            "gene_count": gene_count
        })
        return df


    def GET_cluster(self, barcode_list):
        # Note point_list is a list with position list
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



#----- Generate Gene tables ------
        """ 
        Goal: For each 'tissue', for each 'gene', gene rate a table, such that, one row represent the expression of 
        'gene' in the 'center' spot and surronding 12 spots.
        """
    def gen_table_by_gene(self, gene):
        gen_row = self.dict_feature_to_row.get(gene)
        gen_table = []
        for index, row in self.barcode_list.iterrows():
            if self.dict_barcode_to_coor.get(row[0]) == None:
                gen_table.append(0)
                continue
            temp1 = self.access_neighbour_1(gene,row[0])
            value1 = []
            for i in temp1:
                if i == None: 
                    value1.append(0) # Note, here originally should be None
                    continue
                else: 
                    value1.append(self.tissue_matrix[gen_row,i].tolist())

            gen_table.append( sum(value1)/6)
        df = pd.DataFrame(gen_table)
        return df
    
    

#--------------

    def outlier(df, out_rate = 0.05):
        out_rate_epr = 0
        probs = [.02,.98]
        data_outlier = df
        while(abs(out_rate_epr -(out_rate)) > 0.01 ):  

            probs_temp = probs
            if (out_rate_epr -(out_rate) > 0):
                probs = [probs[0]-0.0005, probs[1]+0.0005]

            if (out_rate_epr -(out_rate) < 0):
                probs = [probs[0]+0.0005, probs[1]-0.0005]


            quart_1 = np.quantile(df.gene_count, probs[0])
            quart_2 = np.quantile(df.gene_count, probs[1])

            IQR = np.quantile(df.gene_count, 0.75) - np.quantile(df.gene_count, 0.25)
            medcp = scipy.stats.skew(df.gene_count)

            Lower = quart_1 - 1.5* np.exp(-4*medcp)*IQR
            Upper = quart_2 + 1.5* np.exp(3*medcp)*IQR 
            temp_i = np.logical_not((df.gene_count > Lower) & (df.gene_count < Upper))

            data_outlier = df[temp_i]
            out_rate_epr = len(data_outlier)/len(df.gene_count)


        return(data_outlier)





    def get_border(self): 
        # we now change border test to level 2
        gene_count = []
        x_mtx =[]
        y_mtx = []
        barcode =[]
        for index, row in self.tissue_position.iterrows():
            if row.iloc[1] == 0:
                continue
            # the following conditioning controls the deepth of testing
            if not ((row.iloc[2] == 0) | (row.iloc[2]==77) | (row.iloc[3] == 126) | (row.iloc[3] == 127) |(row.iloc[3] == 1) | (row.iloc[3]==0) ):
                continue    
            bar = self.dict_coor_to_barcode.get(str(row.iloc[2]) + ' ' + str(row.iloc[3]) + ' ')
            i = self.dict_barcode_to_column.get(bar)
            if i == None:
                continue
            x_mtx.append(row.iloc[2])
            y_mtx.append(row.iloc[3])
            barcode.append(bar)
            gene_count.append(self.tissue_matrix[:,i].sum())
        d = {"barcode":barcode,
            "x_mtx":x_mtx,
                          "y_mtx": y_mtx,
                          "gene_count": gene_count
        }
        df = pd.DataFrame(d)
        return df

    def get_edge(self): 
        # change edge test to level 2
        tissue_position = self.tissue_position
        zero_idx = tissue_position[tissue_position.in_tissue == 0].index.to_list()
        covered_idx = tissue_position[tissue_position.in_tissue == 1].index.to_list()

        neighbors = set()
        
        for index in zero_idx:
                barcode = self.tissue_position.at[index,'barcode']
                neighbors.update(self.get_neighbour_1_bar(barcode))
        covered_tissue = set(self.tissue_position.loc[covered_idx,'barcode'])

        neighbors.intersection_update(covered_tissue)
        
        ngh_iter = neighbors.copy()
        
        #for barcode in ngh_iter:
        #    neighbors.update(self.get_neighbour_1(barcode))
        #neighbors.intersection_update(covered_tissue)
        ###
        # update the zero and covered_idx
        ###
        
        
        
        gene_count = []
        x_mtx =[]
        y_mtx = []
        bar=[]
        for barcode in neighbors:
            coor = self.dict_barcode_to_coor.get(barcode)
            i = self.dict_barcode_to_column.get(barcode)
            if i == None:
                continue
            x_mtx.append(coor[0])
            y_mtx.append(coor[1])
            bar.append(barcode)
            gene_count.append(self.tissue_matrix[:,i].sum())

        d = {"barcode":bar,
            "x_mtx":x_mtx,
                          "y_mtx": y_mtx,
                          "gene_count": gene_count
        }
        df = pd.DataFrame(d)
        return df
      
    
    def get_inner(self, exclude_df):
        # this is just temporary, better to use to edge distance for further analysis
        tissue_position = self.tissue_position.loc[self.tissue_position['in_tissue'] == 1]
        gene_count = []
        x_mtx =[]
        y_mtx = []
        barcode = []
        for index, row in tissue_position.iterrows():
            bar = self.dict_coor_to_barcode.get(str(row[2]) + ' ' + str(row[3]) + ' ')
            i = self.dict_barcode_to_column.get(bar)
            if i == None:
                continue
            x_mtx.append(row[2])
            y_mtx.append(row[3])
            barcode.append(bar)
            gene_count.append(self.tissue_matrix[:,i].sum())
        d = {"barcode":barcode,
            "x_mtx":x_mtx,
                          "y_mtx": y_mtx,
                          "gene_count": gene_count
        }
        df = pd.DataFrame(d)
        #print(df)
        
        df_return = df[~df.barcode.isin(exclude_df.barcode)]
        return df_return
    

########
# This class contain functions removing function
#
#
#
########


# Ensure Artifact_detect is defined before Artifact_remove

# ...existing code for Artifact_detect...

class Artifact_remove(Artifact_detect):

    def __init__(self, dir):
        super().__init__(dir)
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        self.tissue_depth = 0
        self.tissue_depth_df = self.fn_spot_dis_to_edge()
        self.remove_procedure = []
        self.removed_spots = {}

        self.remaining_depth = self.tissue_depth
        

    def fn_spot_dis_to_edge(self, radius_factor=1.35):
        """
        Assigns each tissue spot a depth (layer) from the edge, using spatial neighbor search.
        Returns a DataFrame with columns: barcode, in_tissue, spot_depth
        """
        import numpy as np
        from scipy.spatial import cKDTree
        import copy

        df = copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])
        df['spot_depth'] = 0

        coords_xy = self.coords_xy
        in_tissue = self.in_tissue
        N = coords_xy.shape[0]

        # Build KD-tree for all spots
        tree = cKDTree(coords_xy)
        d1 = tree.query(coords_xy, k=2)[0][:, 1]
        r = radius_factor * np.median(d1)

        # Build adjacency list for all spots (tissue only)
        pairs = list(tree.query_pairs(r=r))
        adjacency = [[] for _ in range(N)]
        for i, j in pairs:
            adjacency[i].append(j)
            adjacency[j].append(i)

        # Only consider tissue spots
        tissue_idx = np.where(in_tissue)[0]
        unassigned = set(tissue_idx)

        # Find edge spots (layer 1)
        edge_mask = find_edge_spots(coords_xy, in_tissue, radius_factor=radius_factor)
        edge_idx = set(np.where(edge_mask)[0])
        df.loc[edge_idx, 'spot_depth'] = 1
        current_layer = edge_idx.copy()
        assigned = edge_idx.copy()
        unassigned -= assigned
        depth = 1

        # Iteratively find next layers
        while unassigned:
            depth += 1
            next_layer = set()
            for idx in current_layer:
                for ngh in adjacency[idx]:
                    if ngh in unassigned:
                        next_layer.add(ngh)
            if not next_layer:
                break
            df.loc[list(next_layer), 'spot_depth'] = depth
            assigned |= next_layer
            unassigned -= next_layer
            current_layer = next_layer

        self.tissue_depth = depth
        print(f"The depth of tissue is {depth}")
        return df

    def fn_spot_inclusion_condition(self):
        df = copy.deepcopy(self.tissue_position[["barcode", "in_tissue"]])

        return df




    def remove_edge(self, distance = 1):
        self.remove_procedure.append("edge")
        # return warn "all spots are removed"
        removed_spot_n = sum( (self.tissue_depth_df.spot_depth <= distance) & (self.spot_inclusion_condition.in_tissue.astype('bool')))
        self.spot_inclusion_condition.loc[self.tissue_depth_df.spot_depth <= distance, 'in_tissue'] = 0
        print(f"We removed {removed_spot_n} edge spots with distance to edge {distance} and less")
        return 

    def remove_border(self):
        self.remove_procedure.append("border")
        df_border_spot = self.get_border()
        for barcode in df_border_spot.barcode:
            self.spot_inclusion_condition.loc[ self.spot_inclusion_condition.barcode == barcode, 'in_tissue'] = 0
        print(f"We removed {len(df_border_spot.barcode)} border spots")
        return

    def remove_malfunction(self):
        # Get all malfunction points, add to remove_procedure list
        self.remove_procedure.append("malfunction")
        df_outlier_spot = self.outlier(self.get_sum())
        for barcode in df_outlier_spot.barcode:
            self.spot_inclusion_condition.loc[self.spot_inclusion_condition.barcode == barcode, 'in_tissue'] = 0
        n_outlier = len(df_outlier_spot.barcode)
        print(f"We removed {n_outlier} outlier spots")
        return
    



#----------- 
# Functions for check object current conditions (in our case check what spots removed, 
# procedure done for removing)

    def simple_cleanser(self): 
        self.remove_procedure = []
        self.removed_spots = {}
        self.spot_inclusion_condition = self.fn_spot_inclusion_condition()
        print("Back to the initial Tissue Sample")
        return 
    
    def review_removing(self):
        [print(i) for i in self.remove_procedure]
        return 
    
    def save(self, dir):
        # Save the position list dataframe 
        df = copy.deepcopy(self.tissue_position)

        df.in_tissue = self.spot_inclusion_condition.in_tissue
        df.to_csv(dir + '/tissue_positions.csv', index = False)

        print(f"Saved the tissue position file to {dir}/tissue_positions.csv")

        return df

    def save_hdf5(self, out_path):
        """
        Save the filtered matrix, barcodes, and features in Visium HDF5 format.
        out_path: path to the output .h5 file
        """

        # Ensure matrix is in CSC format for Visium compatibility
        matrix = self.tissue_matrix
        if not hasattr(matrix, 'tocsc'):
            raise ValueError("tissue_matrix must be a scipy sparse matrix")
        matrix = matrix.tocsc()

        # Prepare barcodes and features
        barcodes = self.barcode_list.iloc[:,0].values.astype('S')
        features = self.feature_list

        # Prepare features subfields
        feature_id = features.iloc[:,0].values.astype('S')
        feature_name = features.iloc[:,1].values.astype('S')
        feature_type = features.iloc[:,2].values.astype('S') if features.shape[1] > 2 else np.array([b'Gene Expression']*len(features))
        genome = np.array([b'GRCh38']*len(features))  # or set appropriately
        all_tag_keys = np.array([b'feature_type'])

        with h5py.File(out_path, 'w') as f:
            grp = f.create_group('matrix')
            # Main matrix datasets
            grp.create_dataset('barcodes', data=barcodes)
            grp.create_dataset('data', data=matrix.data)
            grp.create_dataset('indices', data=matrix.indices)
            grp.create_dataset('indptr', data=matrix.indptr)
            grp.create_dataset('shape', data=np.array(matrix.shape, dtype=np.int32))

            # Features group
            feat_grp = grp.create_group('features')
            feat_grp.create_dataset('id', data=feature_id)
            feat_grp.create_dataset('name', data=feature_name)
            feat_grp.create_dataset('feature_type', data=feature_type)
            feat_grp.create_dataset('genome', data=genome)
            feat_grp.create_dataset('_all_tag_keys', data=all_tag_keys)

        print(f"Saved filtered matrix in Visium HDF5 format to {out_path}")