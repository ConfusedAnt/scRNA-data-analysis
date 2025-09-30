# scRNA-data-analysis

This repository contains the code for the paper: [**Clec4e Potentially Mediates Macrophage Inflammatory Responses in LPS-Induced Lung Injury via Syk–NF-κB Signaling**]. 

## Introduction

Acute lung injury (ALI) induced by lipopolysaccharide (LPS) involves complex remodeling of the lung microenvironment. Using single-cell RNA sequencing, we mapped the cellular landscape of LPS-induced lung injury and identified substantial shifts in cellular composition compared with healthy lungs. Refined clustering and phenotype-linked subset analysis, together with CellChat-based cell–cell communication mapping, revealed that macrophage subsets, particularly Macrophage_c2, c9, and c10, act as primary recipients of intercellular signals, highlighting their central role in coordinating immune responses. Pseudotime analysis revealed that Macrophage_c2 exhibits three distinct differentiation states, reflecting its functional plasticity during lung injury. SCENIC and bulk RNA-seq integration identified Clec4e as a key functional gene specifically enriched in Macrophage_c2. In vitro validation via qPCR and ELISA supported a potential role for Clec4e in mediating macrophage inflammatory responses, potentially through activation of the Syk–NF-κB signaling pathway. Collectively, this study provides a comprehensive single-cell atlas of LPS-induced ALI, delineates macrophage-centered intercellular communication and dynamic polarization, and identifies Clec4e as a promising target for further mechanistic and therapeutic investigation. 

### Graphical abstract
Integration of scRNA-seq and bulk RNA-seq revealed pathogenic immune subpopulations and regulatory genes in LPS-induced ALI. Key findings were supported by qPCR and ELISA, suggesting that Clec4e may act through the Syk–NF-κB pathway to mediate macrophage inflammatory responses.

![scRNA-data-analysis](./docs/Graphical_abstract.png)

### Installation
```bash
# Install micromamba
curl -L https://micro.mamba.pm/install.sh | bash
source ~/.bashrc

# Install Processing Packages 
micromamba create -n R4.4.3 jupyter r-base=4.4.3 r-irkernel r-seurat=4.4.0  -c conda-forge -c bioconda -c r -y
micromamba activate R4.4.3
jupyter kernelspec list

# Install of Analysis Packages in the Readme.txt file under each subdirectory of the ScRNA&Bulk-seq Analysis directory
```

### ScRNA-seq Processing

The single-cell RNA-seq datasets processed during the current study are available in the [GSE276682, GSE292681]. The data is organized under the `1.10X_files/matrix/` directory.  
Each dataset (e.g., `GSM8504078`) includes the following files:
```
1.10X_files/
└── matrix/
└── GSM8504078/
├── barcodes.tsv.gz # Cell barcodes
├── features.tsv.gz # Gene/feature annotations
└── matrix.mtx.gz # Sparse expression matrix
```

```
---- Step1.Read_10X_data.ipynb
---- Step2.process.ipynb
---- Step3.Cell_frequency.ipynb
```
### ScRNA&Bulk-seq Analysis

The single-cell RNA-seq and bulk RNA-seq datasets analyzed during the current study are available in the [GSE276682, GSE292681 and GSE193876]

```
---- 1_ bulk_pheno_data
  ---- Step1. bulk_pheno.ipynb
---- 2_Scissor_pheno
  ---- step2. Scissor_pheno.ipynb
---- 3_cellchat
  ---- step3. Run_cellchat.ipynb
  ---- run_cell_chat.r
---- 4_monocle3
  ---- Step1. bulk_pheno.ipynb
---- 5_Scenic
  ---- Step5_1. scRNA_for_scenic.ipynb
  ---- Step5_2. scenic_mouse.bash
  ---- Step5_3. R_Vis.ipynb
  ---- Readme.txt
---- 6_Bulk_RNAseq
  ---- Step6. Bulk_RNAseq_pro.ipynb
```



