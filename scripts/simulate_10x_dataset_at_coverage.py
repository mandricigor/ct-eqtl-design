
import sys
from scipy.sparse import csr_matrix, save_npz
import pandas
from math import log10
import numpy as np 
from collections import Counter
from scanpy import api
import pandas as pd
from scipy.sparse import csr_matrix, vstack
import time
import pickle


exname = sys.argv[1]
number_individuals = int(sys.argv[2])
singlets_per_individual = int(sys.argv[3])
multiplets_per_individual = int(sys.argv[4])
reads_per_singlet = int(sys.argv[5])
reads_per_doublet = int(sys.argv[6])





data = api.read_h5ad("IGOR_120-individuals_3000-cells.h5ad")

matrix_ = data.X
observations_ = data.obs
observations_.index = observations_.index.astype(int)


genes = list(data.var["gene"])
individuals = list(observations_["ind_cov"].unique())
if number_individuals < len(individuals):
    individuals = np.random.choice(list(observations_["ind_cov"].unique()), size=number_individuals, replace=False)


good_rows = []
doublets = {}
for individual in individuals:
    all_cells_ = observations_[observations_["ind_cov"] == individual].sample(singlets_per_individual + 2 * multiplets_per_individual).index
    doublet_cells = np.random.choice(all_cells_, size=2*multiplets_per_individual, replace=False)
    doublets_in = doublet_cells[:multiplets_per_individual]
    doublets_out = doublet_cells[multiplets_per_individual:]
    for x1, y1 in zip(doublets_in, doublets_out):
        doublets[x1] = y1
    good_rows.extend(list(set(all_cells_) - set(doublets_out)))
good_rows = sorted(good_rows)



observations = observations_.loc[good_rows]





matrices = []
iii = 0
for i in observations.index:
    cellName = i
    row = matrix_[i, ]
    if cellName in doublets:
    	ill_reads = reads_per_doublet
    	rowdata = row.todense() + matrix_[int(doublets[cellName])].todense()
    	rowdata = np.squeeze(np.asarray(rowdata))
    	sumrowdata = np.sum(rowdata)
    
    	row_probs = [1.0 / sumrowdata for x in range(int(sumrowdata))]
    	if sum(row_probs) > 1:
            diff = (sum(row_probs) - 1)
            for ii in range(len(row_probs)):
                if row_probs[ii] > diff:
                    row_probs[ii] -= diff
                    break
    	ill_reads = reads_per_doublet
    	illumina = np.random.multinomial(ill_reads, row_probs)
    	illumina = list(map(lambda x: int(x > 0), illumina))

    	actual_umi_data = np.zeros(len(rowdata), dtype=matrix_.dtype)
    	pointer = 0
    	for kk in range(len(rowdata)):
            umisum = 0
            for jj in range(int(rowdata[kk])):
                umisum += illumina[pointer]
                pointer += 1
            actual_umi_data[kk] = umisum

    	b_row2 = csr_matrix(actual_umi_data)
    	matrices.append(b_row2)

    else:
        rowdata = row.data
        sumrowdata = sum(rowdata)
    
        row_probs = [1.0 / sumrowdata for x in range(int(sumrowdata))]
        if sum(row_probs) > 1:
            diff = (sum(row_probs) - 1)
            for ii in range(len(row_probs)):
                if row_probs[ii] > diff:
                    row_probs[ii] -= diff
                    break
        ill_reads = reads_per_singlet
        illumina = np.random.multinomial(ill_reads, row_probs)
        illumina = list(map(lambda x: int(x > 0), illumina))

        actual_umi_data = np.zeros(len(rowdata), dtype=matrix_.dtype)
        pointer = 0
        for kk in range(len(rowdata)):
            umisum = 0
            for jj in range(int(rowdata[kk])):
                umisum += illumina[pointer]
                pointer += 1
            actual_umi_data[kk] = umisum
    
        b_row = row.tocoo()
        b_row2 = csr_matrix((actual_umi_data, (b_row.row, b_row.col)), shape=row.shape)
        matrices.append(b_row2)

    iii += 1
    if iii % 100 == 0:
        print (iii, "CELLS PROCESSED")




downsampled_sparse = vstack(matrices)
observations.index = range(len(observations.index))


save_npz("dense_CT_june/downsampling_ALLTYPES_%s_%s_%s_%s_%s_%s.MATRIX" % (exname, number_individuals, singlets_per_individual, multiplets_per_individual, reads_per_singlet, reads_per_doublet), downsampled_sparse.tocoo()) 
observations.to_csv("dense_CT_june/downsampling_ALLTYPES_%s_%s_%s_%s_%s_%s.OBSERVATIONS.csv" % (exname, number_individuals, singlets_per_individual, multiplets_per_individual, reads_per_singlet, reads_per_doublet))
genes_df = pd.DataFrame({"genes": genes})
genes_df.to_csv("dense_CT_june/downsampling_ALLTYPES_%s_%s_%s_%s_%s_%s.GENES.csv" % (exname, number_individuals, singlets_per_individual, multiplets_per_individual, reads_per_singlet, reads_per_doublet))


