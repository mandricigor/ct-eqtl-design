

import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


data_csv = sys.argv[1]
ground_csv = sys.argv[2]
#out_csv = sys.argv[3]

data = pd.read_csv(data_csv, header=0, index_col=0)
ground = pd.read_csv(ground_csv, header=0, index_col=0)


newcols = []
for col in data.columns:
    if col.startswith("X"):
        newcols.append(col[1:])
    else:
        newcols.append(col)
data.columns = newcols


newcols = []
for col in ground.columns:
    if col.startswith("X"):
        newcols.append(col[1:])
    else:
        newcols.append(col)
ground.columns = newcols




genes = ground.index


ground = ground[data.columns]


r2s = []

for gene in genes:
    gene1 = gene
    gene2 = gene
    if gene1 == "RP11-442N24__B.1":
        gene1 = "RP11-442N24--B.1"
    if gene1 == "RP11-99J16__A.2":
        gene1 = "RP11-99J16--A.2"
    if gene1 == "RP11-59D5__B.2":
        gene1 = "RP11-59D5--B.2"
    if gene1 == "RP11-445L13__B.3":
        gene1 = "RP11-445L13--B.3"
    if gene1 == "RP11-544L8__B.4":
        gene1 = "RP11-544L8--B.4"
    if gene1 == "XXyac-YX65C7_A.2":
        gene1 = "XXyac-YX65C7-A.2"
    if gene1 == "XXyac-YX65C7_A.3":
        gene1 = "XXyac-YX65C7-A.3"
    if gene1 == "RP11-524D16__A.3":
        gene1 = "RP11-524D16--A.3"
    if gene1 == "RP11-453F18__B.1":
        gene1 = "RP11-453F18--B.1"
    if gene1 == "Metazoa_SRP":
        gene1 = "Metazoa-SRP"
    if gene1 == "Y_RNA":
        gene1 = "Y-RNA"
    if gene1 == "XX-DJ76P10__A.2":
        gene1 = "XX-DJ76P10--A.2"
    if gene1 == "5S_rRNA":
        gene1 = "5S-rRNA"
    if gene1 == "Y_RNA-1":
        gene1 = "Y-RNA-1"
    if gene1 == "Y_RNA-2":
        gene1 = "Y-RNA-2"
    if gene1 == "RP11-1157N2__B.2":
        gene1 = "RP11-1157N2--B.2"
    if gene1 == "RP1-213J1P__B.1":
        gene1 = "RP1-213J1P--B.1"
    if gene1 == "RP1-213J1P__B.2":
        gene1 = "RP1-213J1P--B.2"
    if gene1 == "RP4-633O19__A.1":
        gene1 = "RP4-633O19--A.1"
    if gene1 == "RP4-754E20__A.5":
        gene1 = "RP4-754E20--A.5"
    if gene1 == "CTA-280A3__B.2":
        gene1 = "CTA-280A3--B.2"
    r2s.append(pearsonr(data.loc[gene1], ground.loc[gene2])[0] ** 2)


pd_cor = pd.DataFrame({"r2": r2s}, index=genes)

#pd_cor.to_csv(out_csv, na_rep="NA")

#print (data_csv, np.nanmean(r2s))

celltype = data_csv.split("_")[3]
money = 35000 #int(data_csv.split("/")[0].split("_")[2][4:])
#money = int(data_csv.split("/")[0].split("_")[2][4:])
prevalence = ".".join(data_csv.split("_")[10].split(".")[:2])
nr_individuals = int(data_csv.split("_")[5])
nr_cells = 50 * (int((int(data_csv.split("_")[6]) + int(data_csv.split("_")[7])) / 50.0) + 1)
coverage = 500 * (int(int(data_csv.split("_")[8]) / 500.0))
with open("CT_DMX_correlations_prevalence/%s_%s_%s_%s_%s_%s.csv" % (celltype, money, nr_individuals, nr_cells, coverage, prevalence), "w") as f:
    f.write("%s %s %s %s %s %s %s %s\n" % (celltype, money, nr_individuals, nr_cells, coverage, "%.2f" % np.nanmean(r2s), prevalence, int(round(np.nanmean(r2s) * nr_individuals))))


