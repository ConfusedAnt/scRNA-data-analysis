#remotes::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(Seurat)
# library(SeuratData)
library(qs)
# AvailableData()
# InstallData("pbmc3k")

seurat.data = qread("Combined_Selected_Celltypes.qs")

## 1.2 构建cellchat对象
#pbmc3k里的seurat_annotations有一些NA注释，过滤掉
data.input = seurat.data@assays$RNA@data
meta.data =  seurat.data@meta.data
meta.data = meta.data[!is.na(meta.data$celltype),]
data.input = data.input[,row.names(meta.data)]

#设置因子水平
meta.data$celltype = factor(meta.data$celltype,
                                    #   levels = c("Epithelial_c0", "Epithelial_c1", "Epithelial_c2", "Epithelial_c3", "Epithelial_c4", 
                                    #                 "Epithelial_c5", "Epithelial_c6", "Epithelial_c7", "Epithelial_c8", "Epithelial_c9")
                                                    )

### 1.3 Create a CellChat object
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "celltype")


###可在cellchat对象的meta插槽中添加表型信息
# 添加meta.data信息
cellchat <- addMeta(cellchat, meta = meta.data)

# 设置默认的labels
# levels(cellchat@idents) # show factor levels of the cell labels
# cellchat <- setIdent(cellchat, ident.use = "new.labels") 
# groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

### 1.4 加载CellChat受配体数据库
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)


# Fix "H2-BI" and "H2-Ea-ps" typo in CellChatDB.mouse
CellChatDB$interaction["H2-BI_KIR3DL1","interaction_name"] = "H2-BL_KIR3DL1"
CellChatDB$interaction["H2-BI_KIR3DL1","ligand"] = "H2-Bl"
CellChatDB$interaction["H2-BI_KIR3DL1","interaction_name_2"] = "H2-bl - Kir3dl1"
rownames(CellChatDB$interaction)[rownames(CellChatDB$interaction) == "H2-BI_KIR3DL1"] = "H2-BL_KIR3DL1"

CellChatDB$interaction["H2-EA-PS_CD4","interaction_name"] = "H2-EA_CD4"
CellChatDB$interaction["H2-EA-PS_CD4","ligand"] = "H2-Ea"
CellChatDB$interaction["H2-EA-PS_CD4","interaction_name_2"] = "H2-ea - Cd4"
rownames(CellChatDB$interaction)[rownames(CellChatDB$interaction) == "H2-EA-PS_CD4"] = "H2-EA_CD4"

CellChatDB.use <- CellChatDB # simply use the default CellChatDB
cellchat@DB <- CellChatDB.use

### 1.5 对表达数据进行预处理，用于细胞间通讯分析
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 40) # do parallel


# cellchat <- updateCellChat(cellchat) 如果报错Error in data.use[RsubunitsV, ] : subscript out of bounds)，执行这行代码

print(paste("identifyOverExpressedGenes"))

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

qsave(cellchat, file = "CellChat_Res_0.qs")

print(paste("CellChat_Res_0.qs saved"))

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,population.size = F)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)

df.pathway = subsetCommunication(cellchat,slot.name = "netP")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

qsave(cellchat, file = "CellChat_Res_1.qs")