
import sys
import os
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import spearmanr
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

    cogenefile = "co-expression-eqtl-genes.txt"
    with open(cogenefile) as f:
        lines = f.readlines()
    cogenes = list(map(lambda x: x.strip(), lines))

    cogene_orderdict = {}
    for i, g in enumerate(genes):
        if g in cogenes:
            cogene_orderdict[g] = i

    cogene_index = []
    for i in range(len(cogenes)):
        for j in range(i + 1, len(cogenes)):
            cogene_index.append("%s?%s" % (cogenes[i], cogenes[j]))

    total = 0
    expdict_var = {}
    for x, y in expdict.items():
        y = mat[y, :]
        suka = y.todense()
        suka = np.asarray(suka)
        suka = 1 + 10000 * (suka.T / suka.sum(axis=1)).T
        suka = np.log(suka)
        coex = []
        for i in range(len(cogenes)):
            for j in range(i + 1, len(cogenes)):
                u, v = cogene_orderdict[cogenes[i]], cogene_orderdict[cogenes[j]]
                coex.append(spearmanr(suka[:, u], suka[:, v])[0])
        expdict_var[x] = coex


    expdict_var = pd.DataFrame(expdict_var, index = cogene_index)
    expdict_var.to_csv("EQTL_INPUT/CO_EXPRESSION/%s.csv" % prefix)


