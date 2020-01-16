
import sys
import os
import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd

inputdir = "/u/scratch/i/imandric/UCLA_SINGLE_CELL/dense_CT_november"


#files = list(set([x.split(".")[0] for x in os.listdir(inputdir) if "igor35000" in x]))

# remove this

prefix = sys.argv[1]


files = [prefix]


for x in files:
    expdict = {} # we take only CD4 T cells
    matdata = np.load(inputdir + "/" + x + ".MATRIX.npz")
    mat = csr_matrix((matdata["data"], (matdata["row"], matdata["col"])), shape=matdata["shape"])
    #print (x, mat.shape)
    infctfile = "CT_DMX_igor35000/downsampling_INFERRED_CELL_TYPES_%s.csv" % "_".join(x.split("_")[5:])
    with open(infctfile) as f:
        lines = f.readlines()
    infcelltypes = list(map(lambda x: " ".join(x.strip().split()[1:]), lines[1:]))
    obsfile = inputdir + "/" + x + ".OBSERVATIONS.csv"
    obs = pd.read_csv(obsfile)
    obs = list(obs["ind_cov"])
    for xx, y, z in zip(range(len(infcelltypes)), infcelltypes, obs):
        if y == "CD4 T cells":
            if z not in expdict:
                expdict[z] = []
            expdict[z].append(xx)

    genefile = inputdir + "/" + x + ".GENES.csv"
    genes = pd.read_csv(genefile)
    genes = genes["genes"]
    print (len(infcelltypes), mat.shape, len(obs))
    total = 0
    expdict_var = {}
    for x, y in expdict.items():
        y = mat[y, :]
        suka = y.todense()
        suka = np.asarray(suka)
        suka = 1 + 10000 * (suka.T / suka.sum(axis=1)).T
        suka = np.log(suka)
        expdict_var[x] = suka.var(axis = 0)

    expdict_var = pd.DataFrame(expdict_var, index = genes)
    expdict_var.to_csv("EQTL_INPUT/VARIANCE/%s.csv" % prefix)


