sudo apt-get install libhdf5-dev
micromamba install -c conda-forge r-ggrastr
micromamba install -n R4.4.3 -c conda-forge r-units r-sf ggrastr -y


BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'lme4', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'HDF5Array','ggrastr'))
remotes::install_github("bnprks/BPCells/r")
devtools::install_github('cole-trapnell-lab/monocle3')