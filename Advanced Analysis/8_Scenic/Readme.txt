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



1、先用gen_for_scenic.ipynb生成for.pyscenic.csv
2、修改scenic_mouse.bash中的文件目录
3、nohup bash scenic_mouse.bash  1>pySCENIC.log 2>&1 &


分析的时候需要安装R包：（飞书文档）
https://g2qg57unnc.feishu.cn/docx/HqaGdgMChoMc5Sxhl3oc6P0Gnvd