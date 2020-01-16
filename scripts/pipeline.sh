


expname=$1
nr_individuals=$2
nr_singlets=$3
nr_multiplets=$4
cov_singlets=$5
cov_multiplets=$6

python3.7 super_igor_down_10x_super_multiplex_CT_DMX.py ${expname} ${nr_individuals} ${nr_singlets} ${nr_multiplets} ${cov_singlets} ${cov_multiplets}

Rscript classify_pbmc_cells_by_chunks.november.R dense_CT_november/downsampling_ALLTYPES_${expname}_${nr_individuals}_${nr_singlets}_${nr_multiplets}_${cov_singlets}_${cov_multiplets} ${expname}/downsampling_CELLTYPE_${nr_individuals}_${nr_singlets}_${nr_multiplets}_${cov_singlets}_${cov_multiplets}.csv

