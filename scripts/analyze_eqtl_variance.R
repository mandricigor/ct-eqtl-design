
suppressWarnings(suppressMessages(library(tidyverse)))

args <- commandArgs(trailingOnly=TRUE)
#groundDir <- "GROUND_TRUTH_VAR_EQTL/IGOR_120-individuals_3000-cells_cd4_cells.csv/"
groundDir <- "GROUND_TRUTH_MEAN_EQTL/IGOR_120-individuals_3000-cells_cd4_cells.csv/"
dataDir <- args[1]

groundName <- tail(strsplit(groundDir, "/")[[1]], n = 1)
dataName <- tail(strsplit(dataDir, "/")[[1]], n = 1)




#ground <- read.table("data/super_10x_ground/ground_b_cells.csv/ground_b_cells.csv_chr1_output.txt", header=T)
#data <- read.table("data/15000/downsampling_b_cells_56_1232_17_14872_27067.csv/downsampling_b_cells_56_1232_17_14872_27067.csv_chr1_output.txt", header=T)



ground_egenes_all <- rep(0, 22)
ground_eqtls_all <- rep(0, 22)
data_egenes_all <- rep(0, 22)
data_eqtls_all <- rep(0, 22)
common_egenes_all <- rep(0, 22)
common_eqtls_all <- rep(0, 22)


for (chrom in 1:22) {    
    ground <- read.table(sprintf("%s/chromosome.%s", groundDir, chrom), header=T)
    ground <- ground[ground$FDR < 0.05,]
    ground <- ground %>% unite("eqtl", c("SNP", "gene"), remove=F)

    ground_egenes <- unique(ground$gene)
    ground_eqtls <- unique(ground$eqtl)

    data <- read.table(sprintf("%s/chromosome.%s", dataDir, chrom), header=T)
    data <- data[data$FDR < 0.05,]
    data <- data %>% unite("eqtl", c("SNP", "gene"), remove=F)

    data_egenes <- unique(data$gene)
    data_eqtls <- unique(data$eqtl)

    ground_egenes_all[chrom] <- length(ground_egenes)
    ground_eqtls_all[chrom] <- length(ground_eqtls)
    data_egenes_all[chrom] <- length(data_egenes)
    data_eqtls_all[chrom] <- length(data_eqtls)
    common_egenes_all[chrom] <- length(intersect(ground_egenes, data_egenes))
    common_eqtls_all[chrom] <- length(intersect(ground_eqtls, data_eqtls))
}


ground_egenes_all <- sum(ground_egenes_all)
ground_eqtls_all <- sum(ground_eqtls_all)
data_egenes_all <- sum(data_egenes_all)
data_eqtls_all <- sum(data_eqtls_all)
common_egenes_all <- sum(common_egenes_all)
common_eqtls_all <- sum(common_eqtls_all)

#ground_eqtls_all
#ground_egenes_all
#data_eqtls_all
#data_egenes_all

eqtl_sens <- common_eqtls_all / ground_eqtls_all
eqtl_ppv <- common_eqtls_all / data_eqtls_all
egenes_sens <- common_egenes_all / ground_egenes_all
egenes_ppv <- common_egenes_all / data_egenes_all

print(sprintf("%s %.3f %.3f %.3f %.3f", dataName, eqtl_sens, eqtl_ppv, egenes_sens, egenes_ppv))

