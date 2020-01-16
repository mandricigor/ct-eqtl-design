

import sys

# file1 - ground truth
# file2 - inferred

file1 = sys.argv[1]
file2 = sys.argv[2]

with open(file1) as f:
    a = f.readlines()

a = map(lambda x: x.strip().split(), a)
a = list(a)

with open(file2) as f:
    b = f.readlines()

b = map(lambda x: x.strip().split(), b[1:])
b = list(b)

adict = {}
for x, y in a:
    adict[x] = y

bdict = {}
for x in b:
    bdict[x[0]] = x[6]

count = 0

for key in set(adict.keys()) & set(bdict.keys()):
    if adict[key] == bdict[key]:
        count += 1

print count



