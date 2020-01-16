
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
prevalence = float(sys.argv[7])


demux_perturbation = 1 - np.array([1972,3108,3707,4041,4275,4410,4485,4563,4599,4628,4661,4678,4689,4708,4722,4730,4739,4747,4749,4800]) / 4800.0
cov_milestones = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40]
dmx_perturb = 0
cov_round = reads_per_singlet / 1000 
for i in range(len(cov_milestones)):
    if i == 0:
        if cov_round <= cov_milestones[i]:
            dmx_perturb = demux_perturbation[i]
    elif i > 0 and i <= len(cov_milestones) - 1:
        if cov_milestones[i - 1] < cov_round <= cov_milestones[i]:
            dmx_perturb = demux_perturbation[i]



data = api.read_h5ad("GROUND_TRUTH_CT_DMX/IGOR_120-individuals_3000-cells.h5ad")

matrix_ = data.X
observations_ = data.obs
observations_.index = observations_.index.astype(int)


genes = list(data.var["gene"])
individuals = list(observations_["ind_cov"].unique())

print (Counter(observations_["ind_cov"]))
# permute the individuals labels due to multiplexing performance at difference levels of coverage

ind_minus = {}
for x in individuals:
    ind_minus[x] = [y for y in individuals if y != x]
if dmx_perturb > 0:
    random_cells = np.random.choice(observations_.index, size = int(dmx_perturb * observations_.shape[0]) + 1, replace = False)
    print (len(random_cells))
    for rc in random_cells:
        ni = np.random.choice(ind_minus[observations_.loc[rc]["ind_cov"]])
        observations_.at[rc, "ind_cov"] = ni

print (Counter(observations_["ind_cov"]))
        
# these are individuals who have at least 550 CD4 T cells
good_individuals = ["1019_1019","1054_1054","1056_1056","1127_1127","1132_1132","1196_1196","1203_1203","1219_1219","1221","1243_1243","1250_1250","1262_1262","1279_1279","1338_1338","1340_1340","1404_1404","1419_1419","1420_1420","1479_1479","1492_1492","1496_1496","1510_1510","1514_1514","1545_1545","1558_1558","1596_1596","1597_1597","1602_1602","1615_1615","1621_1621","1623_1623","1667_1667","1716_1716","1726_1726","1730_1730","1743_1743","1767_1767","1768_1768","1775_1775","1791_1791","1827_1827","1848_1848","1891","1892_1892","1958_1958","900033200_900033200","900034200_900034200","900216200_900216200","900759200_900759200","900805200_900805200","901347200_901347200","901457200_901457200","901560200_901560200","902289200_902289200","902299200_902299200","902991200_902991200","903398200_903398200","903648200_903648200","904194200_904194200","904326200_904326200","904344200_904344200","904463200_904463200","904464200_904464200","904477200_904477200","IGTB1290","IGTB143","IGTB1506","IGTB1540","IGTB1542","IGTB1650","IGTB195","IGTB256","IGTB469","IGTB498","IGTB508","IGTB514","IGTB645","IGTB670","IGTB884","IGTB986"]


if number_individuals < len(individuals):
    individuals = np.random.choice(list(observations_["ind_cov"].unique()), size=number_individuals, replace=False)


good_rows = []
doublets = {}
for individual in individuals:
    totalcells = singlets_per_individual + 2 * multiplets_per_individual
    cd4tcells = int(prevalence * totalcells)
    restcells = totalcells - cd4tcells
    cd4_cells_ = observations_[(observations_["ind_cov"] == individual) & (observations_["ct_cov"] == "CD4 T cells")].sample(cd4tcells, replace = True).index
    rest_cells_ = observations_[(observations_["ind_cov"] == individual) & (observations_["ct_cov"] != "CD4 T cells")].sample(restcells, replace = True).index
    print (individual, totalcells, len(cd4_cells_), len(rest_cells_))
    all_cells_ = list(cd4_cells_) + list(rest_cells_)
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


save_npz("dense_CT_november_prevalence/downsampling_ALLTYPES_%s_%s_%s_%s_%s_%s_%s.MATRIX" % (exname, number_individuals, singlets_per_individual, multiplets_per_individual, reads_per_singlet, reads_per_doublet, prevalence), downsampled_sparse.tocoo()) 
observations.to_csv("dense_CT_november_prevalence/downsampling_ALLTYPES_%s_%s_%s_%s_%s_%s_%s.OBSERVATIONS.csv" % (exname, number_individuals, singlets_per_individual, multiplets_per_individual, reads_per_singlet, reads_per_doublet, prevalence))
genes_df = pd.DataFrame({"genes": genes})
genes_df.to_csv("dense_CT_november_prevalence/downsampling_ALLTYPES_%s_%s_%s_%s_%s_%s_%s.GENES.csv" % (exname, number_individuals, singlets_per_individual, multiplets_per_individual, reads_per_singlet, reads_per_doublet, prevalence))


