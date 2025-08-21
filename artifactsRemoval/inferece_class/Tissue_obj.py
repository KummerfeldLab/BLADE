import pandas as pd
import numpy as np
import scipy 






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
        #need to read from the tissue position file
        for index, row in self.tissue_position.iterrows():
            dict[row['barcode']] = [row['x_mtx'], row['y_mtx']]
        
        return dict
            

    def fn_coor_to_barcode(self):
        #need to read from the tissue position file
        dict = {}
        for index, row in self.tissue_position.iterrows():
            temp_coor = [row['x_mtx'], row['y_mtx']]
            dict[''.join(str(x) + ' ' for x in temp_coor)] = row['barcode']

        return dict
    

    def fn_barcode_to_column(self):
        dict = {}
        for index, row in self.barcode_list.iterrows():
            dict[row[0]] = index

        return dict
    
    def fn_feature_to_row(self):
        dict = {}
        for index, row in self.feature_list.iterrows():
            dict[row[1]] = index
        return dict
    