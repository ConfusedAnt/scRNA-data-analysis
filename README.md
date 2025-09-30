# scRNA-data-analysis
Fundamental Analysis &amp; Advanced Analysis

The single-cell RNA-seq and bulk RNA-seq datasets analyzed during the current study are available in the [GSE276682, GSE292681 and GSE193876]. 
For SCENIC analysis, the required files—including mouse feather files, motif2tf, and tf_lists—were obtained from https://resources.aertslab.org/cistarget/



# scRNA-data-analysis

This repository contains the code for the paper: [**Clec4e Potentially Mediates Macrophage Inflammatory Responses in LPS-Induced Lung Injury via Syk–NF-κB Signaling**]. 

## Introduction

Acute lung injury (ALI) induced by lipopolysaccharide (LPS) involves complex remodeling of the lung microenvironment. Using single-cell RNA sequencing, we mapped the cellular landscape of LPS-induced lung injury and identified substantial shifts in cellular composition compared with healthy lungs. Refined clustering and phenotype-linked subset analysis, together with CellChat-based cell–cell communication mapping, revealed that macrophage subsets, particularly Macrophage_c2, c9, and c10, act as primary recipients of intercellular signals, highlighting their central role in coordinating immune responses. Pseudotime analysis revealed that Macrophage_c2 exhibits three distinct differentiation states, reflecting its functional plasticity during lung injury. SCENIC and bulk RNA-seq integration identified Clec4e as a key functional gene specifically enriched in Macrophage_c2. In vitro validation via qPCR and ELISA supported a potential role for Clec4e in mediating macrophage inflammatory responses, potentially through activation of the Syk–NF-κB signaling pathway. Collectively, this study provides a comprehensive single-cell atlas of LPS-induced ALI, delineates macrophage-centered intercellular communication and dynamic polarization, and identifies Clec4e as a promising target for further mechanistic and therapeutic investigation. 

### Graphical abstract
Integration of scRNA-seq and bulk RNA-seq revealed pathogenic immune subpopulations and regulatory genes in LPS-induced ALI. Key findings were supported by qPCR and ELISA, suggesting that Clec4e may act through the Syk–NF-κB pathway to mediate macrophage inflammatory responses.

![FEAOF](./docs/Architecture.png)



### Raw and Split Data

download the Data from:[Google Drive](https://drive.google.com/file/d/1l4zSZVukms-DZDgp-mJ9Zj3ZbXIndTXd/view?usp=sharing) and put it in the following path:

```
---data
  ---raw
  # Wash Smiles
    ---- herg_raw_0528.csv
    ---- herg_raw_0528_clean.csv
    ---- 1_eval_set_herg_60_ini.csv
    ---- 1_eval_set_herg_70_ini.csv

  ---processed
  # Split Datasets
    ---- Train_Val.csv
    ---- Test_1.csv
    ---- Test_2.csv
```
### All Features Data

download the Data from:[Google Drive](https://drive.google.com/file/d/1l4zSZVukms-DZDgp-mJ9Zj3ZbXIndTXd/view?usp=sharing) and put it in the following path:

```
---data
  ---IF_Data
    # Docking Result
    ---- Docking_IF_1.csv
    ---- Docking_IF_2.csv
    ---- Docking_IF_3.csv
    ---- 2_eval_set_60_interaction.csv
    ---- 2_eval_set_70_interaction.csv

    # Select Interaction Fingerprint
    ---- Train_Val_IF.csv
    ---- Test_1_IF.csv
    ---- Test_2_IF.csv
    ---- 5_eval_set_60_IF_filter.csv
    ---- 5_eval_set_70_IF_filter.csv

    # Merge All Features
    ---- Train_All_Features.pkl
    ---- Val_All_Features.pkl
    ---- Test_1_All_Features.pkl
    ---- Test_2_All_Features.pkl
    ---- Set60_All_Features_modified.pkl
    ---- Set70_All_Features_modified.pkl
```

### Trained Models

- download the cpkts in the following link: [Google Drive](https://drive.google.com/file/d/1qWF9zhevi33jVZwWYAG1djefvxiZvm5M/view?usp=drive_link)
- unzip the cpkts to the following path:
```
--- trained_models
    --- DL
        --- CNN-SMILES.pkl
        --- Transformer-TOKENS.pkl
        --- GCN-GRAPH.pkl
        --- MPNN-GRAPH.pkl
    --- ML
        --- GBM_2FP_11Des_1441IF.pkl
        --- RF_2FP_11Des_1441IF.pkl
        --- SVM_2FP_11Des_1441IF.pkl
    --- FEAOF
        --- 0FC_FEAOF.pkl
        --- 1FC_FEAOF.pkl
        --- 3FC_FEAOF.pkl
        --- 5FC_FEAOF.pkl
        --- Protein_0FC_FEAOF.pkl
        --- Protein_3FC_FEAOF.pkl
```

### Retrain and test models
```bash
python Train_FEAOF.py
python Train_DL.py
Train_Eval_ML.ipynb
Eval_FEAOF.ipynb
Eval_DL.ipynb
Eval_FEAOF_valu_set_github.ipynb
```

### Model Performance
![Performance](./docs/Performance.png)



