
library('MatrixEQTL')
library('vcfR')

# Number of samples
#n = 100;

# Number of variables
#ngs = 2000;

# Common signal in all variables (population stratification)
#pop = 0.2 * rnorm(n);

args = commandArgs(trailingOnly=TRUE)


fil <- args[1]
vcf <- args[2]
coords <- args[3]
out <- args[4]

vcf <- read.vcfR(vcf, verbose = FALSE)
vcf_gt <- extract.gt(vcf, as.numeric=T)


samples <- colnames(vcf_gt)
pos <- getPOS(vcf)
chr <- getCHROM(vcf)
snpid <- getID(vcf)


snips <- data.frame(snpid, chr, pos)
genes <- read.csv(coords, header=T)

trimX <- function(a) {ifelse(substr(a, 1, 1) == "X", substr(a, 2, nchar(a)), a)}

x <- read.csv(fil, row.names=1, header=T)



colnames(x) <- trimX(colnames(x))

common <- intersect(colnames(x), samples)





x <- x[common]



y <- as.matrix(x)
vcf_gt <- vcf_gt[,common]


# data matrices
snps.mat = vcf_gt #matrix(rnorm(n*ngs), ncol = ngs) + pop;
gene.mat = y #matrix(rnorm(n*ngs), ncol = ngs) + pop + snps.mat*((1:ngs)/ngs)^9/2;



# data objects for Matrix eQTL engine
snps1 = SlicedData$new( snps.mat );
gene1 = SlicedData$new( gene.mat );
cvrt1 = SlicedData$new( );
rm(snps.mat, gene.mat)

# Slice data in blocks of 500 variables
snps1$ResliceCombined(500);
gene1$ResliceCombined(500);

# name of temporary output file
filename = tempfile();

# Perform analysis recording information for
# a histogram


#meh = Matrix_eQTL_engine(
#  snps = snps1,
#  gene = gene1,
#  cvrt = cvrt1,
#  output_file_name = filename,
#  pvOutputThreshold = 1e-100,
#  useModel = modelLINEAR,
#  errorCovariance = numeric(),
#  verbose = TRUE,
#  pvalue.hist = 100);


meh <- Matrix_eQTL_main( snps1, 
                  gene1, 
                  cvrt = cvrt1, 
                  output_file_name = filename, 
                  pvOutputThreshold = 0,
                  useModel = modelLINEAR, 
                  errorCovariance = numeric(), 
                  verbose = TRUE, 
                  output_file_name.cis = out, 
                  pvOutputThreshold.cis = 0.05,
                  snpspos = snips, 
                  genepos = genes,
                  cisDist = 1e5,
                  pvalue.hist = FALSE,
                  min.pv.by.genesnp = FALSE,
                  noFDRsaveMemory = FALSE)

#write.table(meh$all$eqtls, file="eqtls.csv", quote=F, sep="\t")
unlink( filename );


