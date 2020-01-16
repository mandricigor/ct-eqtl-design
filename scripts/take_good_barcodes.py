
import sys
import pysam

inputfile = sys.argv[1]
outputfile = sys.argv[2]
barcodes = sys.argv[3]

with open(barcodes) as f: 
    a = f.readlines() 
a = list(map(lambda x: x.strip().split()[0], a))


#samfile = pysam.AlignmentFile("immvarYE_0831_1.splitted_1.bam", "rb")
samfile = pysam.AlignmentFile(inputfile, "rb")
sam = samfile.fetch()

samdict = {}
for barcode in a:
    samdict[barcode] = []

for uu in sam:
    try:
        barcode = uu.get_tag("CB")
        if barcode in samdict:
            samdict[barcode].append(uu)
    except Exception as e:
        pass

header = samfile.header

with pysam.AlignmentFile(outputfile, "wb", header=header) as outf: 
    for u, v in samdict.items(): 
        for w in v:
            outf.write(w) 


