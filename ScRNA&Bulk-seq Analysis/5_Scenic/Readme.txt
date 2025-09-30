wget https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
wget -c https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
wget -c https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt

micromamba create -n pyscenic python=3.10 -y
micromamba activate pyscenic

micromamba install -c conda-forge pyarrow
micromamba install -y -c anaconda cytoolz
pip install scanpy -i https://pypi.tuna.tsinghua.edu.cn/simple

pip install pyscenic  -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install umap -i https://pypi.tuna.tsinghua.edu.cn/simple

pip install distributed==2023.12.1（如果报错TypeError: Must supply at least one delayed object）


install.package("BiocManager")
BiocManager::install(c("AUCell", "RcisTarget","GENIE3","zoo", "mixtools", "rbokeh","DT", "NMF", "pheatmap", "R2HTML", "Rtsne","doMC", "doRNG","scRNAseq"))
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC")



1、先用Step5_1. scRNA_for_scenic.ipynb生成for.pyscenic.csv
2、修改Step5_2. scenic_mouse.bash中的文件目录
3、nohup bash scenic_mouse.bash  1>pySCENIC.log 2>&1 &
4、Step5_3. R_Vis.ipynb可视化
