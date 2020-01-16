

import sys
from collections import Counter
import pandas as pd

#checkfile = "CT_DMX_igor15000/downsampling_INFERRED_CELL_TYPES_48_2427_72_3554_6468.csv"
#obsfile = "dense_CT_november/downsampling_ALLTYPES_CT_DMX_igor15000_48_2427_72_3554_6468.OBSERVATIONS.csv"

checkfile = sys.argv[1]
obsfile = sys.argv[2]
groundfile = "GROUND_TRUTH_CT_DMX/IGOR_120-individuals_3000-cells_INFERRED_CELL_TYPES.csv"
groundrow = "GROUND_TRUTH_CT_DMX/IGOR_120-individuals_3000-cells.OBSERVATIONS.csv"


with open(checkfile) as f:
    a = f.readlines()
    a = list(map(lambda x: x.strip().split(), a[1:]))
    inferredCT = list(map(lambda x: " ".join(x[1:]), a))

obs = pd.read_csv(obsfile)
row = list(obs["row"])

with open(groundfile) as f:
    a = f.readlines()
    a = list(map(lambda x: x.strip().split(), a[1:]))
    groundCT = list(map(lambda x: " ".join(x[1:]), a))

obs2 = pd.read_csv(groundrow)
row2 = list(obs2["row"])

groundCTdict = {}
for x, y in zip(row2, groundCT):
    groundCTdict[x] = y

correct = 0
oll = 0
for i in range(len(row)):
    oll += 1
    if inferredCT[i] == groundCTdict[row[i]]:
        correct += 1

correctCT = correct * 1.0 / oll

#CT_DMX_igor15000/downsampling_INFERRED_CELL_TYPES_48_2427_72_3554_6468.cs

cost = int(checkfile.split("_")[2].split("/")[0][4:])
nrindiv = int(checkfile.split("_")[6])
nrcells = int(checkfile.split("_")[7]) + int(checkfile.split("_")[8])
nrcells = int(50 * round(float(nrcells) / 50))


print (cost, nrindiv, nrcells, "%.1f" % (100 - 100 * correctCT))



