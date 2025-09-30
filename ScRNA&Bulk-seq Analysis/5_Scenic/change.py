import os,sys
import loompy as lp
import numpy as np
import scanpy as sc

x=sc.read_csv("for.pyscenic.csv")
row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs)