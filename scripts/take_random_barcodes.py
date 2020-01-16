
import sys
import numpy as np

nr_barcodes = int(sys.argv[1])

with open("YE_8_30_barcodes-1-1-0-0-0-0-0-0-0-0-0-1.all.txt") as f:
    a = f.readlines()

a = list(map(lambda x: x.strip().split(), a))

adict = {}

for u, v in a:
    if v not in adict:
        adict[v] = []
    adict[v].append(u)

take = []
for u, v in adict.items():
    take.extend(np.random.choice(v, size = nr_barcodes, replace = False))


with open("YE_8_30_barcodes-1-1-0-0-0-0-0-0-0-0-0-1.sample-%s.txt" % nr_barcodes, "w") as f:
    for x in take:
        f.write("%s\n" % x)




