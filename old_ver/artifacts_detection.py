import pandas as pd
import numpy as np
import scipy 
import os
import gzip
import copy


class Artifact_detect: 
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
        tissue_position.columns = ['barcode', 'is_covered', 'x_mtx', 'y_mtx', 'x_coor', 'y_coor']
        #tissue_position = tissue_position.loc[tissue_position['is_covered'] == 1]
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
    
###########
    def detection(df1 , df2):
        if len(df1) == 0 | len(df2) == 0:
            return "Missing Border or Edge"

        pvalue = scipy.stats.ttest_ind(df1.gene_count, df2.gene_count,equal_var = False)

        return pvalue 


   

###########

    def get_neighbour_1(self, barcode):
        # Note this function is a transformation of access_neighbour_1 function, however, we do not need gene input
        # at here
        """ Test partial independence of central spot to neighbours. 
        Does s screeen off C from *
                 *  *  *  
               *  s1 s2  *
              * s6  C  s3 *
               *  s5 s4  *
                 *  *  * 

            Note: this function would also return the center point for convinence.
        """


        # firstly get the coor of barcode from self.dict_barcode_to_coor
        center = self.dict_barcode_to_coor.get(barcode)
        # secondly get the neighbour one by one
            # check if in the coor_to_barcode dictionary, if not, enter None
            # if yes, find the right column number and access 
        s1 = [ center[0] - 1 , center[1] - 1]
        s2 = [ center[0] + 1 , center[1] - 1]
        s3 = [ center[0] + 2 , center[1]    ]
        s4 = [ center[0] + 1 , center[1] + 1]
        s5 = [ center[0] - 1 , center[1] + 1]
        s6 = [ center[0] + 1 , center[1] - 1]
        # then, change coor to barcode, then to column.
        bar1 = self.dict_coor_to_barcode.get( ''.join(str(x) + ' ' for x in s1) )
        bar2 = self.dict_coor_to_barcode.get( ''.join(str(x) + ' ' for x in s2) )
        bar3 = self.dict_coor_to_barcode.get( ''.join(str(x) + ' ' for x in s3) )
        bar4 = self.dict_coor_to_barcode.get( ''.join(str(x) + ' ' for x in s4) )
        bar5 = self.dict_coor_to_barcode.get( ''.join(str(x) + ' ' for x in s5) )
        bar6 = self.dict_coor_to_barcode.get( ''.join(str(x) + ' ' for x in s6) )

        

        # finnally add everything to list
        return {bar1,bar2,bar3,bar4,bar5,bar6}

###########



    def get_border(self): 
        # we now change border test to level 2
        gene_count = []
        x_mtx =[]
        y_mtx = []
        barcode =[]
        for index, row in self.tissue_position.iterrows():
            if row[1] == 0:
                continue
            # the following conditioning controls the deepth of testing
            if not ((row[2] == 0) | (row[2]==77) | (row[3] == 126) | (row[3] == 127) |(row[3] == 1) | (row[3]==0) ):
                continue    
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
        return df

    def get_edge(self): 
        # change edge test to level 2
        tissue_position = self.tissue_position
        zero_idx = tissue_position[tissue_position.is_covered == 0].index.to_list()
        covered_idx = tissue_position[tissue_position.is_covered == 1].index.to_list()

        neighbors = set()
        
        for index in zero_idx:
                barcode = self.tissue_position.at[index,'barcode']
                neighbors.update(self.get_neighbour_1(barcode))
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
        tissue_position = self.tissue_position.loc[self.tissue_position['is_covered'] == 1]
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
    

