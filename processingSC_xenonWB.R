#https://satijalab.org/seurat/articles/integration_introduction.html
#https://satijalab.org/seurat/articles/hashing_vignette.html

library(Seurat)
library(patchwork)
library(tidyverse)
library(openxlsx)

options(future.globals.maxSize = 4000 * 1024^5)

##################################  DATA PREPARATION  ############################################
projectname = 'xenonZY2'

# Read in the data
data1 <- Read10X(data.dir = paste0('data_files/3245_atmo'))
data2 <- Read10X(data.dir = paste0('data_files/3406_atmo'))
data3 <- Read10X(data.dir = paste0('data_files/3246_xenon'))
data4 <- Read10X(data.dir = paste0('data_files/3247_xenon'))

sobj1 <- CreateSeuratObject(counts = data1, project = '3245_atmo')
sobj2 <- CreateSeuratObject(counts = data2, project = '3406_atmo')
sobj3 <- CreateSeuratObject(counts = data3, project = '3246_xenon')
sobj4 <- CreateSeuratObject(counts = data4, project = '3247_xenon')

############################  MERGING ###########################################
# Merge files
scrna_merged <- merge(sobj1, y = c(sobj2, sobj3, sobj4), add.cell.ids = c('3245_atmo', '3406_atmo', '3246_xenon', '3247_xenon'), project = projectname)
table(scrna_merged$orig.ident)

######################### MITOCHONDRIAL REGRESSION ############################################

# Store mitochondrial gene statatistics in your Seurat object
scrna_merged[['percent_mt']] <- PercentageFeatureSet(scrna_merged, pattern = '^mt-')

# Basic QC plot to set cutoffs
VlnPlot(scrna_merged, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)

# Filter data to remove unwanted cells
scrna_final <- subset(scrna_merged, nFeature_RNA > 200 & nCount_RNA < 25000 & percent_mt < 20)

Idents(scrna_final) <- 'orig.ident'
scrna_final$orig.ident<- factor(x = scrna_final$orig.ident, levels = c('3245_atmo', '3406_atmo', '3246_xenon', '3247_xenon'))

sample_counts <- table(scrna_final$orig.ident)
sample_counts
write.xlsx(sample_counts, file = paste0('results/sampleCounts2_', projectname, '.xlsx'))

dev.off()
pdf(file = paste0('plots/QC/vln_QC_', projectname, '.pdf'), pointsize = 10)
VlnPlot(scrna_final, features = c('nFeature_RNA', 'nCount_RNA', 'percent_mt'), ncol = 3)
dev.off()

pdf(file = paste0('plots/QC/scatter_QC_', projectname, '.pdf'), pointsize = 10)
scat1 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'percent_mt')
scat2 <- FeatureScatter(scrna_final, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
scat1 + scat2
dev.off()

########################### NORMALIZATION AND SCALING #######################################
# Normalize RNA expression data and scale to variable features
scrna_final <- NormalizeData(scrna_final)
scrna_final <- FindVariableFeatures(scrna_final, selection.method = 'vst')
scrna_final <- ScaleData(scrna_final, features = VariableFeatures(scrna_final), vars.to.regress = 'percent_mt')


#################################### DIMENSIONAL REDUCTIONS##############################
# Run principal component analysis
scrna_final <- RunPCA(scrna_final, features = VariableFeatures(object = scrna_final))
scrna_final <- JackStraw(object = scrna_final, num.replicate = 50, prop.freq=0.025, dims = 50)
scrna_final <- ScoreJackStraw(scrna_final, dims = 1:50)


# Visualize PCA to ensure merged samples are comparable
Idents(scrna_final) <- 'orig.ident'
pdf(file = paste0('plots/QC/pca_QC_', projectname, '.pdf'), pointsize = 10)
DimPlot(scrna_final, reduction = 'pca')
dev.off()

# Visualize component strengths to decide how many to use
pdf(file = paste0('plots/QC/components_QC_', projectname, '.pdf'), pointsize = 10)
JackStrawPlot <- JackStrawPlot(object = scrna_final, dims = 1:50, xmax = 0.05) + guides(col = guide_legend(ncol = 1)) + theme(legend.text = element_text(size = 6), legend.key.size = unit(0.02, "cm"))
ElbowPlot <- ElbowPlot(object = scrna_final, ndims = 50)
JackStrawPlot + ElbowPlot
dev.off()


# Cluster cells according to elbow plot dimension choice
scrna_final <- FindNeighbors(scrna_final, dims = 1:26)
scrna_final <- FindClusters(scrna_final, resolution = 0.8)

# Run UMAP reduction to visualize clusters
scrna_final <- RunUMAP(scrna_final, dims = 1:32)

# Plot UMAP
pdf(file = paste0('plots/umap_', projectname, '.pdf'), width = 12, height = 9)
DimPlot(scrna_final, reduction = "umap", split.by = 'orig.ident', ncol = 2)
dev.off()

# Save final object
saveRDS(scrna_final, file = paste0('seurat_objects/SeuratObject_', projectname, '.rds'))

