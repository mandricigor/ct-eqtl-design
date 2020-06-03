
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
fil_mean <- args[2]
vcf <- args[3]
coords <- args[4]
out <- args[5]

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


z <- read.csv(fil_mean, row.names=1, header=T)
colnames(z) <- trimX(colnames(z))
z <- z[common]





y <- as.matrix(x)
vcf_gt <- vcf_gt[,common]


# data matrices
snps.mat = vcf_gt #matrix(rnorm(n*ngs), ncol = ngs) + pop;
gene.mat = y #matrix(rnorm(n*ngs), ncol = ngs) + pop + snps.mat*((1:ngs)/ngs)^9/2;

covs <- read.csv("covariates.txt", row.names = 1)
covs <- covs[common, ]

# data objects for Matrix eQTL engine
snps1 = SlicedData$new( snps.mat );
#gene1 = SlicedData$new( gene.mat );
#cvrt1 = SlicedData$new( t(covs) );
#rm(snps.mat, gene.mat)

# Slice data in blocks of 500 variables
snps1$ResliceCombined(500);

# name of temporary output file
#filename = tempfile();

# Perform analysis recording information for
# a histogram



results <- data.frame(matrix(ncol = 6, nrow = 0))
xxx <- c("snps", "gene", "statistic", "pvalue", "FDR", "beta")
colnames(results) <- xxx


for (gen in rownames(x)) {
    try({
    covar <- t(covs)
    covar <- as.matrix(rbind(covar, z[gen,colnames(covar)]))
    cvrt1 = SlicedData$new( covar );
    gsuka <- t(as.matrix(gene.mat[gen,]))
    rownames(gsuka) <- c(gen)
    gene1 = SlicedData$new( gsuka );
    meh <- Matrix_eQTL_main( snps1, 
                  gene1, 
                  cvrt = cvrt1, 
                  output_file_name = "/dev/null", 
                  pvOutputThreshold = 0,
                  useModel = modelLINEAR, 
                  errorCovariance = numeric(), 
                  verbose = TRUE, 
                  output_file_name.cis = out, 
                  pvOutputThreshold.cis = 0.05,
                  snpspos = snips, 
                  genepos = as.data.frame(genes[genes$geneid == gen,]),
                  cisDist = 1e5,
                  pvalue.hist = FALSE,
                  min.pv.by.genesnp = FALSE,
                  noFDRsaveMemory = FALSE)
    results <- rbind(results, meh$cis$eqtls)},
    silent = T)
}

results <- results[c("snps", "gene", "beta", "statistic", "pvalue", "FDR")]
colnames(results) <- c("SNP", "gene", "beta", "t-stat", "p-value", "FDR")
results$FDR <- p.adjust(results$FDR, method = "fdr")
write.table(results, file=out, quote=F, sep = "\t", row.names = F)



