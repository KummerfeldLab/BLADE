import pandas as pd
import numpy as np
import scipy 
from scipy import io
import os
import gzip


class Malfunction: 
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
        #self.adjmatrix = self.fn_adjmatrix()

    
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
        tissue_position = tissue_position.loc[tissue_position['is_covered'] == 1]
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

    
#############################
    def get_sum(self): 
        # we now change border test to level 2
        gene_count = []
        x_mtx =[]
        y_mtx = []
        barcode =[]
        for index, row in self.tissue_position.iterrows():
 
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

#############################

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

        
        col1 = self.dict_barcode_to_column.get(bar1)
        col2 = self.dict_barcode_to_column.get(bar2)
        col3 = self.dict_barcode_to_column.get(bar3)
        col4 = self.dict_barcode_to_column.get(bar4)
        col5 = self.dict_barcode_to_column.get(bar5)
        col6 = self.dict_barcode_to_column.get(bar6)

        # finnally add everything to list
        return [ col1,col2,col3,col4,col5,col6]

    

#----- Saving gene tables -------


    def nh_variable(self, gene_sub_list):
        df = pd.DataFrame(self.barcode_list)
        for i in gene_sub_list:
            df[i+'_nh'] = self.gen_table_by_gene(i)
        return df

    def fn_adjmatrix(self):
        adjM = np.zeros((self.tissue_position.shape[0], self.tissue_position.shape[0]),dtype=int)

        for index, row in self.barcode_list.iterrows():

            temp1 = self.get_neighbour_1(row[0])
            for i in temp1:
                if i == None:
                    continue
                adjM[index, i] = 1

        return adjM


#----- Work on adjMatrix -------
    def BFS(self, point, full_list):
        #cluster_temp = []
        adjmatrix = self.fn_adjmatrix()
        visited = []
        queue = [] 
        visited.append(point)
        queue.append(point)
        while queue:          
            m = queue.pop(0) 
            print (m, end = " ") 
            temp_ngh = np.where(adjmatrix[m,:] == 1)[0].tolist()
            temp_sub_ngh = list( set(temp_ngh) & set(full_list) )
            for neighbour in temp_sub_ngh:
                if neighbour not in visited:
                    visited.append(neighbour)
                    queue.append(neighbour)

        return visited

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

