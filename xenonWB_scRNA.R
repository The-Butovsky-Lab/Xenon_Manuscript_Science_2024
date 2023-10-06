library(tidyverse)
library(Seurat)
library(Signac)
library(enrichR)
library(openxlsx)
library(patchwork)
library(data.table)
library(dittoSeq)
library(ggplot2)
library(SCP)
library(RColorBrewer)
library(EnhancedVolcano)
library(tidyverse)
library(ggsignif)
library(ggpubr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(limma)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(data.table)
library(stringr)
library(fgsea)
library(enrichplot)


projectname = 'xenon_KK'


########################  DIFFERENTIAL FEATURE PROCESSING ##########################################
# Load Seurat data set
dataset <- readRDS('seurat_objects/SeuratObject_260.8xenon.rds')

######################## METADATA MANIPULATION AND CELLTYPE ASSIGNMENT ##########################
# Cluster mapping
Idents(dataset_filtered) <- 'seurat_clusters'
levels(dataset_filtered)
new_mappings <- c('atmosphere', 'atmosphere', 'xenon', 'xenon')
new_mappings <- c('M0','MGnD','M0','Ribosome-Intermediate','HSP-Intermediate','Cytokine-Intermediate',
                  'Interferon','Cycling-S','M0','Cycling-G2M','MGnD','Macrophage','MGnD-Antigen')
new_mappings <- c('M0','MGnD','M0','Intermediate','Intermediate','Intermediate',
                  'Interferon','Cycling','M0','Cycling','MGnD','Macrophage','MGnD')
new_mappings <- c('M0','MGnD','M0','Lineage1','Lineage2','Cytokine-Intermediate',
                  'Lineage2','Cycling-S','M0','Cycling-G2M','MGnD','Macrophage','MGnD-Antigen')
new_mappings <- c('1-M0','6-MGnD','1-M0','3-Ribosome-Intermediate','4-HSP-Intermediate','2-Cytokine-Intermediate',
                  '5-Interferon','Cycling-S','1-M0','Cycling-G2M','6-MGnD','Macrophage','7-MGnD-Antigen')
new_mappings <- c('M0','MGnD','M0','Pre-MGnD_Ribosome','Pre-MGnD','Pre-MGnD_Cytokine',
                  'Pre-MGnD','Cycling','M0','Cycling','MGnD','Macrophage','MGnD')


names(new_mappings) <- levels(dataset_filtered)
dataset_filtered <- RenameIdents(dataset_filtered, new_mappings)

dataset_filtered$SubCellType_ordered <- Idents(dataset_filtered)
dataset_filtered$Treatment <- Idents(dataset_filtered)
dataset_filtered$SubCellType <- Idents(dataset_filtered)
dataset_filtered$CellType <- Idents(dataset_filtered)
dataset_filtered$Lineage <- Idents(dataset_filtered)
dataset_filtered$SubCellType_Intermediate <- Idents(dataset_filtered)
dataset_filtered$SubCellType_PreMGnD <- Idents(dataset_filtered)

#Filter out Macrophages
dataset_filtered_no_macrophage_no_cycling <- subset(dataset_filtered, idents= c('M0','MGnD','Pre-MGnD_Ribosome','Pre-MGnD','Pre-MGnD_Cytokine'))

########################  QC (Ribosomal, mitochondiral cell cycle) ##########################################
CellCycleScoring(dataset, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> data
data[[]]
as_tibble(data[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

as_tibble(data[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))

Cell_cycle_per_cluster <- data@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  dplyr::count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster") +
  scale_fill_manual(values=c('blue', 'forestgreen', 'red')) +
  theme_minimal()
dev.off()
pdf(file = paste0('plots/cell_cycle_per_cluster', projectname, '.pdf'), width = 12, height = 9)
Cell_cycle_per_cluster
dev.off()

Cell_cycle_cluster <- ClassDimPlot(
  srt = data, group.by = "seurat_clusters", stat.by = "Phase",
  reduction = "UMAP", theme_use = "theme_blank"
)

dev.off()
pdf(file = paste0('plots/cell_cycle_per_cluster_vizualised_', projectname, '.pdf'), width = 12, height = 9)
Cell_cycle_cluster
dev.off()

dataset <- RunCellQC(srt = dataset)
QC_plot <- ClassDimPlot(srt = dataset, group.by = "CellQC", reduction = "UMAP")
dev.off()
pdf(file = paste0('plots/QC_plot_', projectname, '.pdf'), width = 12, height = 9)
QC_plot
dev.off()

ClassStatPlot <-  ClassStatPlot(srt = dataset, stat.by = "CellQC", group.by = "seurat_clusters", label = TRUE)
dev.off()
pdf(file = paste0('plots/QC_barplot_per_cluster_', projectname, '.pdf'), width = 12, height = 9)
ClassStatPlot
dev.off()

#Filtering taking only confidently assigned cells
Idents(dataset) = 'seurat_clusters'
dataset_filtered <- subset(dataset, idents = c('0','1','2','3','4','5','6','7','8','9','11','12','13'))
Idents(dataset_filtered) = 'CellQC'
dataset_filtered <- subset(dataset_filtered, idents = 'Pass')


####### TOTAL COUNTS ########
Idents(dataset_filtered) <- 'Treatment'
dataset_filtered_xenon <- subset(dataset_filtered, ident = 'xenon')
dataset_filtered_air <- subset(dataset_filtered, ident = 'atmosphere')

# Gene expression
Idents(dataset_filtered_xenon) <- 'SubCellType'
exprTable <- AverageExpression(dataset_filtered_xenon)
exprTable <- exprTable[['RNA']] %>% as.data.frame() 
exprTable <- exprTable %>% select(order(colnames(exprTable))) %>% select(everything())
write.xlsx(exprTable, paste0('results/geneExpression_', 'dataset_filtered_xenon', 'SubCellType','.xlsx'), rowNames = T, overwrite = T)

# Gene expression
Idents(dataset_filtered) <- 'seurat_clusters'
exprTable <- AverageExpression(dataset_filtered)
exprTable <- exprTable[['RNA']] %>% as.data.frame() 
exprTable <- exprTable %>% select(order(colnames(exprTable))) %>% select(gene, everything())
write.xlsx(exprTable, paste0('results/geneExpression_', 'dataset_filtered_no_macrophage', 'SubCellType','.xlsx'), rowNames = T, overwrite = T)



################################## DIFFERENTIAL EXPRESSION #######################################
Idents(dataset_filtered) = 'seurat_clusters'
Idents(dataset_filtered) = 'orig.ident'
Idents(dataset_filtered) = 'Treatment'
Idents(dataset_filtered) = 'SubCellType'
Idents(dataset_filtered_no_macrophage) = 'SubCellType'

dataset_filtered_intermediate <-  subset(dataset_filtered, idents = c('Ribosome-Intermediate','HSP-Intermediate','Cytokine-Intermediate',
                                         'Interferon'))

dataset_filtered_subcelltype_PreMGnD_markers <- FindAllMarkers(dataset_filtered_no_macrophage_no_cycling, only.pos = T, min.pct = 0.05, logfc.threshold = 0.25)
TOP_dataset_filtered_subcelltype_markers <- dataset_filtered_subcelltype_PreMGnD_markers %>% group_by(cluster) %>% top_n(n=50, wt = avg_log2FC)

CytokineIntermediate_xenon_vs_control <- FindMarkers(dataset_filtered, ident.1 = 'xenon', ident.2 = 'atmosphere', group.by = 'Treatment', 
                                                     subset.ident = 'Cytokine-Intermediate', 
                                                     min.pct = 0.05, logfc.threshold = 0.05, min.diff.pct = 0.0) %>%rownames_to_column('gene')
heatmap_markers <- cluster_markers %>% slice_min(order_by = p_val, n = 100)
write.xlsx(dataset_filtered_intermediate_no_macrophage_markers, paste0('results/dataset_filtered_intermediate_no_macrophage_markers','.xlsx'))
#All microglia
All_microglia_xenon_vs_control <- FindMarkers(dataset_filtered_no_macrophage, ident.1 = 'xenon', ident.2 = 'atmosphere', group.by = 'Treatment', 
                                             only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0) %>% rownames_to_column('gene')

write.xlsx(TOP_dataset_filtered_subcelltype_markers, paste0('results/TOP_dataset_filtered_subcelltype_markers.xlsx'))

####### DEGs PER SUBCELLTYPE #########
Idents(dataset_filtered_no_macrophage_no_cycling) <- 'SubCellType_PreMGnD'

SubCellTypes_no_macrophage <- levels(dataset_filtered_no_macrophage_no_cycling)
SubCellTypes_no_macrophage
for (i in SubCellTypes_no_macrophage) {
  Subcelltype_xenon_vs_control <- FindMarkers(dataset_filtered_no_macrophage_no_cycling, ident.1 = 'xenon', ident.2 = 'atmosphere', group.by = 'Treatment', 
                                              subset.ident = i, only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0) %>% rownames_to_column('gene')
  assign(paste0(i,'_xenon_vs_control_volcano'), Subcelltype_xenon_vs_control)
}

rm(df_total)
df_total <- data.frame()
for (i in dataset_filtered_no_macrophage_no_cycling) {
  count <- nrow(get(paste0(i,'_xenon_vs_control_padj02')))
  df <- data.frame(count)
  df_total <- rbind(df_total, df)
}

df_total <- cbind(df_total,SubCellTypes_no_macrophage)

pointSize = 20 
lineWidth = 1

df_total$SubCellTypes_no_macrophage <- factor(df_total$SubCellTypes_no_macrophage, levels=c("M0",
          "Ribosome-Intermediate",
          "HSP-Intermediate",
          "Cytokine-Intermediate",
          "Interferon",
          "Cycling-S",
          "Cycling-G2M",
          "MGnD",
          "MGnD-Antigen"))
  
degs_bar <- ggplot(df_total,aes(x = SubCellTypes_no_macrophage, y = count, fill = SubCellTypes_no_macrophage))+
            geom_col() + scale_fill_manual(values = c("#A6CEE3",  "#B2DF8A", "#33A02C",  "#FDBF6F", "#FF7F00" ,"#FB9A99", "#E31A1C","#1F78B4","#CAB2D6"))+
        theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"))+
  labs(y= 'Count of DEGs (Pval < 0.01)', x=NULL) +
  ggtitle('Number of DEGs per cluster (Xenon vs Air)')+
  scale_y_continuous( expand = c(0,0.1,0,0)) 
degs_bar

ggsave(
  paste0('plots/barplot_DEGs_per_SubCellType_no_macrophage_CUSTOM_pval01', '.pdf'),
  plot = degs_bar,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 400,
  height = 300,
  units = c("mm"),
  dpi = 600)




################################### VISUALIZATION ##############################################
Idents(dataset_filtered) = 'seurat_clusters'
Idents(dataset_filtered_no_macrophage) = 'SubCellType_PreMGnD'

#level active ident 
levels(dataset_filtered_no_macrophage)
dataset_filtered_no_macrophage@active.ident <- factor(dataset_filtered_no_macrophage@active.ident,
                                                      levels=c( "M0",
                                                                "Ribosome-Intermediate",
                                                                "HSP-Intermediate",
                                                                "Cytokine-Intermediate",
                                                                "Interferon",
                                                                "Cycling-S",
                                                                "Cycling-G2M",
                                                                "MGnD",
                                                                "MGnD-Antigen"))
levels(dataset_filtered_no_macrophage)

######## COLORS ##########
#Set up color scheme
cols <- brewer.pal(5, "Paired")
cols
cols_SubCellType <- c("#A6CEE3", "#B2DF8A", "#33A02C","#E31A1C","#FF7F00", "#FB9A99","#FDBF6F","#1F78B4","#6A3D9A")

levels(dataset_filtered_no_macrophage)
Idents(dataset_filtered_no_macrophage) <- 'SubCellType_PreMGnD' 
cols_CellType <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6","#CAB2D6","#CAB2D6")
cols_SubCellType <- c("#A6CEE3", #M0
                      "#1F78B4", #MGnD
                      "#B2DF8A", #Ribosome
                      "#33A02C", #HSP
                      "#FDBF6F", #Cytokine
                      "#FF7F00", #Interferon
                      "#FB9A99", #Cycling-S
                      "#E31A1C", #Cycling-G2M
                      "#6A3D9A") #MGnD-Antigen

cols_SubCellType_intermediate <- c("#A6CEE3", #M0
                      "#1F78B4", #MGnD
                      "#B2DF8A", #Ribosome
                      "#33A02C", #Interferon
                      "#FDBF6F") #Cytokine


######## GENERAL DIMPLOT FEATUREPLOT ########
pdf(file = paste0('plots/umap_SubCellType_PreMGnD', 'no_macrophage_no_cycling_split_treatment', '.pdf'), width = 18, height = 9)
DimPlot(dataset_filtered_no_macrophage_no_cycling, label = F, label.box = F, cols = cols_SubCellType_intermediate, pt.size = 1, split.by = 'Treatment')
dev.off()

pdf(file = paste0('plots/umap_Featureplot_','Mrc1','Clec12a', '.pdf'), width = 10, height = 5)
FeaturePlot(dataset_filtered,features = c('Mrc1','Clec12a'), pt.size = 1, order =T, min.cutoff = 0.5)
dev.off()

######## MODULE PLOT ########
#INF
INF_module = list(dataset_filtered_subcelltype_PreMGnD_markers %>% filter(cluster == 'Pre-MGnD') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset_filtered_no_macrophage_no_cycling,
                             features = INF_module,
                             name="INF_module")
INF_module_featureplot <- FeaturePlot(dataset_filtered_no_macrophage,
                                      features = "INF_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Interferon_Module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
INF_module_featureplot
dev.off()


#M0
M0_module = list(c('P2ry12', 'Fosb', 'Ifngr1', 'Maf', 'Tmem119', 'Jund','Jun', 'Egr1'))
M0_module <- AddModuleScore(dataset_filtered_no_macrophage,
                            features = M0_module,
                            name="M0_module")
M0_module_featureplot <- FeaturePlot(M0_module,
                                     features = "M0_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_M0_Module_', 'dataset_filtered_no_macrophage','.pdf'), width = 9, height = 9)
M0_module_featureplot
dev.off()


#MGnD
MGnD_module = list(c('Clec7a','Axl','Igf1','Lpl','Spp1','Cst7'))
MGnD_module <- AddModuleScore(dataset_filtered_no_macrophage,
                              features = MGnD_module,
                              name="MGnD_module")
MGnD_module_featureplot <- FeaturePlot(MGnD_module,
                                       features = "MGnD_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_MGnD_module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
MGnD_module_featureplot
dev.off()

#S-Phase
S_phase_module = list(c('Mcm5',
                     'Mcm6',
                     'Lig1',
                     'Mcm4',
                     'Mcm7',
                     'Stmn1',
                     'Mcm3'))
S_phase_module <- AddModuleScore(dataset_filtered_no_macrophage,
                              features = S_phase_module,
                              name="S_phase_module")
S_phase_module_featureplot <- FeaturePlot(S_phase_module,
                                       features = "S_phase_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_S_phase_module_featureplot_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
S_phase_module_featureplot
dev.off()

#Ribosome
Ribosome_module = list(Sub_celltype_markers_filtered_top_markers %>% filter(cluster == 'Ribosome-Intermediate') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset_filtered_no_macrophage,
                                 features = Ribosome_module,
                                 name="Ribosome_module")
Ribosome_module_featureplot <- FeaturePlot(dataset_filtered_no_macrophage,
                                          features = "Ribosome_module1", label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Ribosome_module_featureplot_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
Ribosome_module_featureplot
dev.off()

#G2M phase
G2M_module = list(c('Top2a', 'Mki67' ))
G2M_module <- AddModuleScore(dataset_filtered_no_macrophage,
                                  features = G2M_module,
                                  name="G2M_module")
G2M_module_featureplot <- FeaturePlot(G2M_module,
                                           features = "G2M_module1",label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_G2M_module_featureplot_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
G2M_module_featureplot
dev.off()

#AntigenMGnD
Antigen_module = list(c('H2-Aa','H2-Eb1', 'H2-Ab1'))
Antigen_module <- AddModuleScore(dataset_filtered_no_macrophage,
                             features = Antigen_module,
                             name="Antigen_module")
Antigen_module_featureplot <- FeaturePlot(Antigen_module,
                                      features = "Antigen_module1",label = TRUE, repel = TRUE) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/Featureplot_Antigen_module_','dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
Antigen_module_featureplot
dev.off()

#Macrophages
Macrophage_module = list(c('Mrc1','Clec12a', 'Ms4a7'))
Macrophage_module <- AddModuleScore(dataset_filtered,
                                 features = Macrophage_module,
                                 name="Macrophage_module")
Macrophage_module_featureplot <- FeaturePlot(Macrophage_module,
                                          features = "Macrophage_module1", label = F, pt.size = 1, cols = c('grey', 'red')) 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Macrophage_module_featureplot_', projectname, '.pdf'), width = 12, height = 9)
Macrophage_module_featureplot
dev.off()

#HSP-module
HSP_module = list(Sub_celltype_markers_filtered_top_markers %>% filter(cluster == 'HSP-Intermediate') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset_filtered_no_macrophage,
                                    features = HSP_module,
                                    name="HSP_module")
HSP_module_featureplot <- FeaturePlot(HSP_module,features = "HSP_module1", label = T,  cols = c('grey', 'red')) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_HSP_module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
HSP_module_featureplot
dev.off()


#Cytokine-module
Cytokine_module = list(Sub_celltype_markers_filtered_top_markers %>% filter(cluster == 'Cytokine-Intermediate') %>% pull(gene))
dataset_filtered_no_macrophage <- AddModuleScore(dataset_filtered_no_macrophage,
                             features = Cytokine_module,
                             name="Cytokine_module")
Cytokine_module_featureplot <- FeaturePlot(Cytokine_module,features = "Cytokine_module1", label = T,  cols = c('grey', 'red')) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pdf(file = paste0('plots/umap_Featureplot_Cytokine_module_', 'dataset_filtered_no_macrophage', '.pdf'), width = 9, height = 9)
Cytokine_module_featureplot
dev.off()



######## VIOLIN PLOT ############
dittoPlot(dataset, var = 'Apoe', group.by = "CellType", split.by = 'Treatment', plots = c('vlnplot', 'jitter'), jitter.size = 0.5, split.ncol = 1)

Module = list(c(Sub_celltype_markers_filtered %>% filter(cluster == 'MGnD-Antigen' & p_val_adj < 0.05) %>% pull(gene)))
dataset_filtered_no_macrophage <- AddModuleScore(dataset_filtered_no_macrophage,
                                     features = Module,
                                     name="MGnDAntigen")

Vln_plot <- ExpVlnPlot(
  srt = dataset_filtered_no_macrophage, group.by = "Treatment", 
  features = c("MGnDAntigen1"),
  comparisons = list(
    c("xenon", "atmosphere")
  ),
  multiplegroup_comparisons = FALSE
)

Vln_plot

dev.off()
pdf(file = paste0('plots/Violin_plot_', 'Intermediate_combined','_Xenon_vs_Air', '.pdf'), width = 12, height = 9)
Vln_plot
dev.off()

######## CUSTOM ENRICHMENT VIOLIN PLOT #########
Module_score <- dataset_filtered_no_macrophage@meta.data$INF_module1
Treatment <- dataset_filtered_no_macrophage@meta.data$Treatment

vln_data <- map2_dfr(Module_score, Treatment,~ tibble(Module_score = .x, Treatment = .y))

lineWidth = 0.5
pointSize = 15

Vln_plot_custom <- ggplot(vln_data, aes(x = Treatment, y = Module_score, fill = Treatment,group = Treatment)) + 
  geom_violin() +
  stat_summary(fun = "median",geom = "crossbar",width = 0.5,colour = "black")+
  geom_point(aes(x = Treatment), position = position_jitterdodge(jitter.width = 1, jitter.height=0, dodge.width=0.9), size = 0.5, alpha = 0.1) +
  theme(text = element_text(size = pointSize, colour = "black"),
    rect = element_blank(),
    line = element_line(size = lineWidth, colour = "black"),
    plot.title  = element_text(color="black", size=20),
    axis.title  = element_text(size = pointSize, colour = "black"),
    axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = pointSize , colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = pointSize, colour = "black"),
    legend.text = element_text(size = pointSize , colour = "black"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    axis.line = element_line(size = lineWidth, colour = "black"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + 
  labs(title = 'Pre-MGnD Interferon Module') +
  stat_compare_means(comparison =  list(c('atmosphere', 'xenon')), p.adjust.method = "bonferroni",  method = 'wilcox.test', label = "p.signif") +  
  stat_compare_means(label.x = 1.2, label.y = 2)+
  scale_y_continuous(expand = c(0, 0, .05, 0))+
  #scale_fill_manual(values = c(  "#FDBF6F",'#b37627')) #Cytokine Colors
  #scale_fill_manual(values = c("#FF7F00",'#ab5500')) #HSP Colors
  scale_fill_manual(values = c("#33A02C",'#10360e')) #HSP Colors
  
Vln_plot_custom


dev.off()
pdf(file = paste0('plots/Violin_plot_', 'Custom_Pre_MGnD_Module_enrichment','_Xenon_vs_Air', '.jpeg'), width = 5, height = 9)
Vln_plot_custom
dev.off()

######## BARPLOT ##########
dev.off()
Idents(dataset_filtered_no_macrophage_no_cycling) <- 'Treatment'
barplot <- dittoBarPlot(dataset_filtered_no_macrophage_no_cycling, 
                        var = 'SubCellType_PreMGnD',  
                        group.by = "Treatment",
                        color.panel = cols_SubCellType_intermediate,
                        var.labels.reorder = c(1,2,3,5,4))

pdf(file = paste0('plots/barplot_SubCellType_PreMGnD_no_macrophage_', '.pdf'), width = 5, height = 9)
barplot
dev.off()

barplot_data <- as.data.frame(barplot$data)
#barplot_data
#write.xlsx(barplot_data,paste0('results/barplot_counts_subbcelltype_per_origident.xlsx'))
gg_bar <- ggplot(barplot_data, aes(x=grouping, y=percent, fill = label)) +
  geom_bar(position = "fill", stat = "identity",width=0.9) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0.01,0,0)) +
  geom_text(aes(label = paste0(round(percent*100,digits = 2),"%")), 
            position = position_stack(vjust = 0.5), size = 5)+
  scale_fill_manual(values = c("#A6CEE3",  "#B2DF8A", "#33A02C",  "#FDBF6F", "#FF7F00" ,"#FB9A99", "#E31A1C","#1F78B4","#CAB2D6"))+
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20, face="bold.italic"),
        axis.title  = element_text(size = pointSize, colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 45, hjust = 1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.line = element_line(size = lineWidth, colour = "black"))

ggsave(
  paste0('plots/barplot_SubCellType_no_macrophage_CUSTOM', '.pdf'),
  plot = gg_bar,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 200,
  height = 300,
  units = c("mm"),
  dpi = 600)

######## DONUT PLOT #######
barplot_data_air <- barplot_data %>%filter(grouping == 'atmosphere')
barplot_data_xenon <- barplot_data %>%filter(grouping == 'xenon')

# Compute the cumulative percentages (top of each rectangle)
barplot_data_air$ymax = cumsum(barplot_data_air$percent)
barplot_data_xenon$ymax = cumsum(barplot_data_xenon$percent)

# Compute the bottom of each rectangle
barplot_data_air$ymin = c(0, head(barplot_data_air$ymax, n=-1))
barplot_data_xenon$ymin = c(0, head(barplot_data_xenon$ymax, n=-1))
barplot_data_air

# Compute label position
barplot_data_air$labelPosition <- (barplot_data_air$ymax + barplot_data_air$ymin) / 2

# Compute a good label
barplot_data_air$label_good <- paste0(barplot_data_air$label, "\n value: ", barplot_data_air$count)

ggplot(barplot_data_air, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=label)) +
  geom_rect() +
  geom_text(x=2, aes(y=labelPosition, label=label_good, color=label), size=6) + # x here controls label position (inner / outer)
  scale_fill_brewer(palette=3) +
  scale_color_brewer(palette=3) +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.position = "none")


percentage <- bar$data
write.xlsx(percentage, 'results/percentage_subclusters_no_macrophage.xlsx')

dev.off()
pdf(file = paste0('plots/barplot_scaled_seurat_clusters_filtered_', projectname, '.pdf'), width = 5, height = 9)
dittoBarPlot(dataset_filtered, 
             var = 'Lineage',  
             group.by = "Treatment")
dev.off()


######## HEATMAP #########
DoHeatmap(dataset_filtered,
          features = TOP_dataset_filtered_subcelltype_markers$gene,
          raster = TRUE)



ht <- ExpHeatmap(
  srt = dataset_filtered, features = TOP_dataset_filtered_subcelltype_markers$gene, feature_split = TOP_dataset_filtered_subcelltype_markers$cluster, cell_split_by = "SubCellType",
  species = "Mus_musculus", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  row_title_size = 0, height = 8, width = 10)
ht


pdf(file = paste0('plots/heatmap_Subcelltype_pathways_sc','.pdf'), width = 20, height = 15)
print(ht$plot)
dev.off()

ggsave(paste0('plots/heatmap_pathways_SubCellType' , '.jpeg'),
  plot = ht$plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 500,
  height = 300,
  units = c("mm"),
  dpi = 300)

######## VOLCANO ######## 


selected_label <- c('Anxa5','Cct5','Cd36','Cct7','C5ar1','C3ar1','Il10ra','Tanc2','Lpl','Blnk',
'Arf4','Tgfb1','Cd68','Psen2','Flt1','Msr1','Ms4a7','Clec7a','Lgals3','Ccl3',
'Sorl1','CD33', 'App','Inpp5d','Apoe','Ch25h','Cx3cr1','Pfdn1','Ccl9','Fcer1g',
'Ccl4', 'Cd33', 'Il13ra1','Mef2a', 'Tpcn1', 'Cd69', 'Cd300lb','Lpl','Fcgr4',
'Axl','Lilrb4a','Ctsb')

pdf(file = paste0('plots/', '.pdf'), pointsize = 10, width = 5, height = 6)
EnhancedVolcano(All_microglia_xenon_vs_control,
                lab = All_microglia_xenon_vs_control$gene,
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c(),
                title = paste0('Xenon vs Ctrl'),
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 0.5,
                labSize = 4,
                # xlim = c(-3,3),
                # ylim = c(0,6),
                legendLabSize = 8,
                legendPosition = 'top',
                drawConnectors = T,
                widthConnectors = 0.2,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                max.overlaps = Inf,
                col=c('grey', 'grey', 'grey', 'red3'))
dev.off()

######## SCATTER PLOT ####################
FeaturePlot(dataset_filtered_no_macrophage, features = 'Stat2')
volcano <- All_microglia_xenon_vs_control
volcano$log_pval <- -log10(volcano$p_val)

sig_data <- volcano %>% filter(p_val < 0.05)
UP_data <- sig_data %>% filter(avg_log2FC > 0) 
DOWN_data <- sig_data %>% filter(avg_log2FC < 0)

##Theme
theme_scatter <- theme(aspect.ratio = 1, 
                       panel.background = element_blank(),
                       panel.border=element_rect(fill=NA),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.background=element_blank(),
                       axis.title  = element_text(size = 20, colour = "black"),
                       axis.text.x=element_text(colour="black",size =20),
                       axis.text.y=element_text(colour="black",size =20),
                       axis.ticks=element_line(colour="black"),
                       plot.margin=unit(c(1,1,1,1),"line"),
                       plot.title = element_text(colour="black",size =20))

#Select data
select_DOWN <- c('Ccl3', 'Nlrp3', 'Atf3', 'Junb', 'Egr1', 'Itgb5', 'Ccl4', 'Ccl12', 'Nfkbiz', 'Lyz2', 'Ccrl2', 'Tyrobp', 'Ier2', 'Sgk1','H3f3b')

select_UP <- c('Ifi27l2a', 'Mid1', 'Hspaa1', 'Hspa1a', 'Hspa1b', 'Rps27', 'Rps28', 'H2-Q7', 'Ifitm3', 'Klf2','Dnaja1','Klf4')

select_UP_data <- sig_data %>% filter(gene %in% select_UP)
select_DOWN_data <- sig_data %>% filter(gene %in% select_DOWN)


scatter_volcano <- ggplot(volcano, aes(x= avg_log2FC, y= log_pval, label = gene)) +
  geom_point(aes(fill = avg_log2FC),size = 3, shape = 21, color = 'black') +
  scale_fill_gradient2(limits=c(-0.5, 0.53), low="blue", mid="whitesmoke", high = "red", na.value = 'blue') +
  theme_scatter +
  xlim(-1,1)+
  # ylim(0,9)+
  geom_vline(xintercept = 0.1,colour="grey", linetype = "longdash") +
  geom_vline(xintercept = -0.1,colour="grey", linetype = "longdash") +
  geom_hline(yintercept = 1.3, colour="grey", linetype = "longdash") +
  labs(title = str_wrap("All microglia Xenon vs Air",60),
       x = "Log2FC",
       y = "-log(pvalue)") +
  geom_label_repel(data = select_UP_data,
                   color = 'black', size = 8,fontface = 'italic',
                   min.segment.length = 0,
                   nudge_y = 30,
                   nudge_x = 0.2,
                   max.overlaps = Inf,
                   point.padding = unit(0.5, 'mm')) +
  geom_label_repel(data = select_DOWN_data,
                   color = 'black',size = 8, fontface = 'italic',
                   min.segment.length =0,
                   nudge_y = 30,
                   nudge_x = -0.2,
                   max.overlaps = Inf,
                   point.padding = unit(0.5, 'mm'))
scatter_volcano 


ggsave(paste0('plots/volcano','All_microglia_xenon_vs_air', '.pdf'),
  plot = scatter_volcano,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 20,
  units = c("mm"),
  dpi = 600)

######## DOTPLOT CELLTYPES ##########
TOP_dataset_filtered_intermediate_no_macrophage_markers
dataset_filtered_no_macrophage@active.ident <- factor(dataset_filtered_no_macrophage@active.ident,
                                                        levels=c( "M0",
                                                                  "Ribosome-Intermediate",
                                                                  "HSP-Intermediate",
                                                                  "Cytokine-Intermediate",
                                                                  "Interferon",
                                                                  "Cycling-S",
                                                                  "Cycling-G2M",
                                                                  'MGnD',
                                                                  'MGnD-Antigen'))
dataset_filtered_no_macrophage@active.ident <- factor(dataset_filtered_no_macrophage@active.ident,
                                                      levels=c( 'MGnD-Antigen',
                                                                'MGnD',
                                                                "Cycling-S",
                                                                "Cycling-G2M",
                                                                "Ribosome-Intermediate",
                                                                "Interferon",
                                                                "HSP-Intermediate",
                                                                "Cytokine-Intermediate",
                                                                "M0"))
dataset_filtered_no_macrophage@active.ident <- factor(dataset_filtered_no_macrophage@active.ident,
                                                      levels=c( 'MGnD',
                                                                "Cycling",
                                                                "Intermediate",
                                                                "Cytokine-Intermediate"))

Dotplot <- DotPlot(dataset_filtered_intermediate, features = c(
  #'Tmem119', 'P2ry12', 'Malat1', 'Maf', 'Jun', 'Gpr34',                  #M0
   'Ifit1' ,'Ifit3', 'Stat1', 'Irf7', 'Isg15', 'Clec2d', 'Ccl12',             #Interferon
   'Ccl2', 'Ccl3','Ccl4','Sgk1','Camk2n1','Ier2','Nfkbia', 'Tnf',        #Cytokine
  'Hsp90aa1','Hspa8', 'Hspa5','Hspa1a','Hspa1b',                     #HSP-Intermediate
  'Rps29','Rps27', 'Rpl30','Rplp1','Fau'))  +                              #Ribosome-Intermediate
  #'Top2a', 'Mki67', 'Ube2c','Nusap1','Birc5',                            #Cycling-G2M
  #'Mcm5','Mcm6','Lig1','Mcm4','Mcm7',                                    #Cycling-S
  #'Lpl', 'Clec7a','Cst7', 'Lyz2', 'Tyrobp','Apoe',                       #MGnD
  #'H2-Aa','H2-Eb1', 'H2-Ab1')) +                                         #MGnD-Antigen
  theme(axis.text.x = element_text(angle = 60, face = 'italic', vjust= 0.01, size = 15),
        axis.text.y  = element_text(size = 15 , colour = "black"))+
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Dotplot

ggsave(
  paste0('plots/Dotplot', 'dataset_filtered_intermediate_SUBCELLTYES', '.pdf'),
  plot = Dotplot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 400,
  height = 300,
  units = c( "mm"),
  dpi = 600
)

FeaturePlot(dataset_filtered_no_macrophage, features = 'H2-Q7')
#################### FUNCTIONAL ANALYSIS ##################
######## STATS CELL SUBSETS #########
devtools::install_github("KlugerLab/DAseq")
library(DAseq)

python2use <- "/Users/kiliankleemann/opt/anaconda3/bin/python3"


labels_atms <- unique(barplot_data[barplot_data$treatment == "atmo", "grouping"])
labels_xenon <- unique(barplot_data[barplot_data$treatment == "xen", "grouping"])

labels_res <- X.label.info[X.label.info$condition == "R", "label"]
labels_nonres <- X.label.info[X.label.info$condition == "NR", "label"]

X.melanoma
X.label.melanoma
cell_labels
cell_labels <- as.character.factor(dataset_filtered_no_macrophage@meta.data$Treatment)
dataset_filtered_no_macrophage@reductions$umap
X_DAseq <- dataset_filtered_no_macrophage@reductions[["umap"]]@cell.embeddings
X_DAseq 

da_cells <- getDAcells(
  X = X_DAseq,
  cell.labels = cell_labels,
  labels.1 = labels_atms,
  labels.2 = labels_xenon,
  k.vector = seq(50, 500, 50),
  plot.embedding = dataset_filtered_no_macrophage_DA)


#Regression analysis
barplot_data_edit <- read.xlsx(paste0('results/barplot_counts_subbcelltype_per_origident.xlsx'))
barplot_data_edit <- separate(barplot_data_edit, col=grouping, into= c('grouping', 'treatment'), sep='_')
print(barplot_data_edit)

barplot_data_M0 <- barplot_data_edit %>% filter(label == 'Cytokine-Intermediate')
M0_glm <- glm(count ~ log2(total) + treatment, data = barplot_data_M0)


anova(M0_glm)
summary(M0_glm)

counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))

#Example
## an example with offsets from Venables & Ripley (2002, p.189)
utils::data(anorexia, package = "MASS")

anorex.1 <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
                family = gaussian, data = anorexia)
summary(anorex.1)
# }
# NOT RUN {
# A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
  u = c(5,10,15,20,30,40,60,80,100),
  lot1 = c(118,58,42,35,27,25,21,19,18),
  lot2 = c(69,35,26,21,18,16,13,12,12))
summary(glm(lot1 ~ log(u), data = clotting, family = Gamma))
summary(glm(lot2 ~ log(u), data = clotting, family = Gamma))
## Aliased ("S"ingular) -> 1 NA coefficient
(fS <- glm(lot2 ~ log(u) + log(u^2), data = clotting, family = Gamma))
tools::assertError(update(fS, singular.ok=FALSE), verbose=interactive())


######## DETERMINE INTERMEDIATE / TRANSITIONARY ####################################
Idents(dataset_filtered) = 'seurat_clusters'

#Intermediate vs M0
Clust5_vs_Clust0 <- FindMarkers(dataset_filtered, ident.1 = '5', ident.2 = '0', 
                               min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')

Clust5_vs_Clust0 <- FindMarkers(dataset_filtered, ident.1 = '5', ident.2 = '0', 
                                min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')
write.xlsx(Clust5_vs_Clust0, 'results/Clust5_vs_Clust0_dataset_filtered.xlsx')

Clust5_vs_Clust2 <- FindMarkers(dataset_filtered, ident.1 = '5', ident.2 = '2', 
                                min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')
write.xlsx(Clust5_vs_Clust2, 'results/Clust5_vs_Clust2_dataset_filtered.xlsx')


Clust4_vs_Clust0 <- FindMarkers(dataset_filtered, ident.1 = '4', ident.2 = '0', 
                                min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')
write.xlsx(Clust4_vs_Clust0, 'results/Clust4_vs_Clust0_dataset_filtered.xlsx')

Clust4_vs_Clust2 <- FindMarkers(dataset_filtered, ident.1 = '4', ident.2 = '2', 
                                min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')
write.xlsx(Clust4_vs_Clust2, 'results/Clust4_vs_Clust2_dataset_filtered.xlsx')


Clust4_vs_Clust5 <- FindMarkers(dataset_filtered, ident.1 = '4', ident.2 = '5', 
                                min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')
write.xlsx(Clust4_vs_Clust5, 'results/Clust4_vs_Clust5_dataset_filtered.xlsx')


#MGnD
Clust11_vs_Clust1 <- FindMarkers(dataset_filtered, ident.1 = '11', ident.2 = '1', 
                                min.pct = 0.05, logfc.threshold = 0.5) %>% rownames_to_column('gene')
write.xlsx(Clust4_vs_Clust5, 'results/Clust4_vs_Clust5_dataset_filtered.xlsx')

#Compare Intermediate 
Cytokine_vs_HSP <- FindMarkers(dataset_filtered_intermediate, ident.1 = 'Cytokine-Intermediate', ident.2 = 'HSP-Intermediate', group.by = 'SubCellType', min.pct = 0.05, logfc.threshold = 0.25, min.diff.pct = 0.0, ) %>% rownames_to_column('gene')
Cytokine_vs_Interferon <- FindMarkers(dataset_filtered_intermediate, ident.1 = 'Cytokine-Intermediate', ident.2 = 'HSP-Intermediate', group.by = 'SubCellType', min.pct = 0.05, logfc.threshold = 0.25, min.diff.pct = 0.0, ) %>% rownames_to_column('gene')

Intermediate_INF_HSP_Ribosome_vs_M0 <- FindMarkers(dataset_filtered, ident.1 = 'Intermediate', ident.2 = 'M0', min.pct = 0.05, logfc.threshold = 0, min.diff.pct = 0.0, ) %>% rownames_to_column('gene')
Intermediate_INF_HSP_Ribosome_vs_MGnD <- FindMarkers(dataset_filtered, ident.1 = 'Intermediate', ident.2 = 'MGnD', min.pct = 0.05, logfc.threshold = 0, min.diff.pct = 0.0, ) %>% rownames_to_column('gene')


#Macrophage
Macrophage_ <- FindAllMarkers(dataset_filtered, ident.1 = 'Macrophage', ident.2 = '5', 
                                min.pct = 0.05, logfc.threshold = 0.25) %>% rownames_to_column('gene')
write.xlsx(Clust4_vs_Clust5, 'results/Clust4_vs_Clust5_dataset_filtered.xlsx')



######## TRAJECTORY ASSESSMENT ####################################
Idents(dataset_filtered) <- 'SubCellType'
dataset_filtered_subset <- subset(dataset_filtered_no_macrophage, idents = c('M0','MGnD','Cytokine-Intermediate',
                                                               'Intermediate'))
Idents(dataset_filtered_no_macrophage_no_cycling) <- 'Treatment'
dataset_filtered_subset_xenon <- subset(dataset_filtered_no_macrophage_no_cycling, idents = 'xenon')
dataset_filtered_subset_air <- subset(dataset_filtered_no_macrophage_no_cycling, idents = 'atmosphere')

dataset_filtered_subset_air.slingshot <- RunSlingshot(srt = dataset_filtered_subset_air, group.by = "SubCellType_PreMGnD", reduction = "umap", start = 'M0')

# dataset_filtered_subset_slingshot@active.ident <- factor(dataset_filtered_subset_slingshot@active.ident,
#                                                         levels=c('M0','Ribosome-Intermediate','HSP-Intermediate',
#                                                                  'Interferon','MGnD-Antigen',
#                                                                  'MGnD','Cytokine-Intermediate'))
# levels(dataset_filtered_subset_slingshot)

pdf(file = paste0('plots/Trajectory_plot_','dataset_filtered_subset_air.slingshot','_SubCellType_PreMGnD' , projectname, '.pdf'), width = 10, height = 10)
Trajectory_plot <- ClassDimPlot(dataset_filtered_subset_xenon.slingshot, group.by = "SubCellType_PreMGnD", reduction = "umap", lineages = paste0("Lineage", 1:2), lineages_span = 0.3,lineages_line_bg = "white", palette = "Paired")
dev.off()

Trajectory_plot

ggsave(
  paste0('plots/Trajectory_plot_','dataset_filtered_subset_air.slingshot','_SubCellType_PreMGnD' , '.pdf'),
  plot = Trajectory_plot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c("mm"),
  dpi = 300)


#Dynamic Plot
dynamicplot <- DynamicPlot(
  srt = dataset_filtered_subset_xenon.slingshot, lineages = c("Lineage1", "Lineage2"), group.by = "SubCellType_PreMGnD",
  exp_method = 'log1p',
  features = c('Tmem119','P2ry12', 'Stat1','Clec2d'))

dynamicplot

ggsave(
  paste0('plots/dataset_filtered_subset_XENON_slingshot_','Cst7_Clec7a','PreMGnD.pdf'),
  plot = dynamicplot,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c( "mm"),
  dpi = 300)

#Expression Plots
ExpDimPlot_trajectory <- ExpDimPlot(dataset_filtered_subset_xenon.slingshot, features = paste0("Lineage", 1:2), reduction = "UMAP", theme_use = "theme_blank")
ExpDimPlot_trajectory

ggsave(
  paste0('plots/ExpDimPlot_trajectory_','dataset_filtered_subset_slingshot_xenon' ,'.pdf'),
  plot = ExpDimPlot_trajectory,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 300,
  height = 300,
  units = c( "mm"),
  dpi = 300)


######## Monocle 3 analysis####
gene_annotation <- as.data.frame(rownames(dataset_filtered_no_macrophage@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(dataset_filtered_no_macrophage@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(dataset_filtered_no_macrophage@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = dataset_filtered_no_macrophage@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


New_matrix <- dataset_filtered_no_macrophage@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(dataset_filtered_no_macrophage@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

library(monocle3)

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition



list_cluster <- dataset_filtered_no_macrophage@active.ident
names(list_cluster) <- dataset_filtered_no_macrophage@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-dataset_filtered_no_macrophage@reductions[["umap"]]@cell.embeddings


cds_from_seurat@preprocess_aux$gene_loadings <- dataset_filtered_no_macrophage@reductions[["pca"]]@feature.loadings
dataset_filtered_no_macrophage@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=4)


######## PATHWAY ANALYSIS ########
#Setting significant genes
DEgenes <- volcano %>% filter(p_val < 0.05) %>% pull('gene')
DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
DEG_Log2FC <-  volcano %>% filter(p_val < 0.05) %>%  pull('avg_log2FC')
names(DEG_Log2FC) = DEgenes
DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
DEG_Log2FC
###

pathway_data <- Sub_celltype_markers_filtered
rm(df_total_pathway)
df_total_pathway <- data.frame()
for (i in SubCellTypes_no_macrophage){
  pathway_analysis <- pathway_data %>% filter(cluster == paste0(i))
  DEgenes <- pathway_analysis %>% filter(p_val_adj < 0.05) %>% pull('gene')
  DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
  go <- enrichGO(gene = DEgenes,
                 'org.Mm.eg.db',
                 keyType = "ENTREZID",
                 ont = "BP",
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
  p1 <- as.data.frame(go@result) %>% filter(pvalue < 0.05)
  p1$log_p <- -log10(p1$pvalue)
  p1$group <- paste0(i)
  df <- data.frame(p1)
  df_total_pathway <- rbind(df_total_pathway, df)
}
write.xlsx(df_total_pathway, paste0('results/total_pathways_GO_SubCellType.xlsx'))
selected_pathways <- read.xlsx(paste0('results/total_pathways_GO_SubCellType.xlsx'), sheet = 2) %>% pull(selected)

df_total_pathway_selected <- df_total_pathway %>% filter(ID%in% selected_pathways)

df_total_pathway_selected_zscore

df_total_pathway_selected$group <- factor(df_total_pathway_selected$group,
                                                      levels=c( "M0",
                                                                "HSP-Intermediate",
                                                                "Cytokine-Intermediate",
                                                                "Interferon",
                                                                "Ribosome-Intermediate",
                                                                "Cycling-S",
                                                                "Cycling-G2M",
                                                                "MGnD",
                                                                "MGnD-Antigen"))
pointSize=10
lineWidth=1
dotplot_GO <- ggplot(df_total_pathway_selected,  # you can replace the numbers to the row number of pathway of your interest
                     aes(x = group, y = Description)) + 
  geom_point(aes(size = Count, color = log_p)) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x =  element_text(size = pointSize , colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90, vjust=0.1),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  scale_colour_gradient2(limits=c(0, 65), low="navyblue", mid="antiquewhite", high = "red") +
  labs(title = str_wrap("GO pathway enrichment", 20))
dotplot_GO


#GSEA GO
pathway_data <- Sub_celltype_markers_filtered
rm(df_total_pathway)
df_total_enrichment_pathway <- data.frame()
for (i in SubCellTypes_no_macrophage){
  pathway_analysis <- pathway_data %>% filter(cluster == paste0(i))
  DEgenes <- pathway_analysis %>% filter(p_val_adj < 0.05) %>% pull('gene')
  DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
  DEG_Log2FC <-  pathway_analysis %>% filter(p_val < 0.05) %>%  pull('avg_log2FC')
  names(DEG_Log2FC) = DEgenes
  DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
  go_result_BP <- gseGO(geneList = DEG_Log2FC,
                        OrgDb    = org.Mm.eg.db,
                        ont      = "BP",
                        # minGSSize = 15, 
                        # maxGSSize = 500,
                        nPermSimple = 1000,
                        pvalueCutoff = 0.1,
                        scoreType = 'pos')
  p1 <- as.data.frame(go_result_BP@result) %>% filter(pvalue < 0.1)
  p1$log_p <- -log10(p1$pvalue)
  p1$group <- paste0(i)
  df <- data.frame(p1)
  df_total_enrichment_pathway <- rbind(df_total_enrichment_pathway, df)
}


#GSEA KEGG
kegg_result <- gseKEGG(geneList = DEG_Log2FC,
                       organism     = 'mmu',
                       pAdjustMethod = "BH",
                       nPermSimple = 10000,
                       pvalueCutoff = 0.05)

head(kegg_result)
kegg_result <- setReadable(kegg_result, 'org.Mm.eg.db', 'ENTREZID')
kegg_result_dataframe <- as.data.frame(kegg_result)
write.xlsx(kegg_result_dataframe, paste0('results/KEGG_enrichment_All_microglia_xenon_vs_air', '.xlsx'))

cols_SubCellType <- c("#A6CEE3", #M0
                      "#1F78B4", #MGnD
                      "#B2DF8A", #Ribosome
                      "#33A02C", #HSP
                      "#FDBF6F", #Cytokine
                      "#FF7F00", #Interferon
                      "#FB9A99", #Cycling-S
                      "#E31A1C", #Cycling-G2M
                      "#6A3D9A") #MGnD-Antigen
#UPREGULATED KEGG
gseasplot2_KEGG <- gseaplot2(kegg_result,geneSetID = c(1,5,9),color = c("#6A3D9A","#33A02C","#B2DF8A"), pvalue_table = TRUE,
                             base_size = 20, subplots = 1:3)
gseasplot2_KEGG2 <- gseaplot2(kegg_result,geneSetID = 1,color = "#B2DF8A", pvalue_table = TRUE,
                             base_size = 20)

#DOWNREGULATED KEGG
gseasplot2_KEGG3 <- gseaplot2(kegg_result,geneSetID = c(2,4,7,17),color = c("#FDBF6F","#E31A1C","#FB9A99",'#FF7F00'), pvalue_table = TRUE,
                             base_size = 20)

gseasplot2_KEGG3
ggsave(paste0('plots/GseaPlot_','All_microglia_xenon_vs_air_','DOWN','.pdf'),
       gseasplot2_KEGG3,
       units = c("mm"),
       width = 300,
       height = 200,
       dpi = 300)

#Hallmark Pathways
pathways.hallmark <- gmtPathways("gsea_pathways/mh.all.v2022.1.Mm.entrez.gmt") 
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=DEG_Log2FC)

fgseaResTidy_hallmark <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Hallmark_pathways.xlsx'))


#Curated gene sets
pathways.curated <- gmtPathways("gsea_pathways/m2.all.v2022.1.Mm.entrez.gmt") 
fgseaRes <- fgsea(pathways=pathways.curated, stats=DEG_Log2FC)

fgseaResTidy_curated <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Hallmark_pathways.xlsx'))

#GO gene sets
pathways.go <- gmtPathways("gsea_pathways/m5.go.v2022.1.Mm.entrez.gmt") 
fgseaRes <- fgsea(pathways=pathways.go, stats=DEG_Log2FC)

fgseaResTidy_go <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)

write.xlsx(fgseaResTidy, paste0('results/Pathway_analysis/FGSEA_enrichment/',file_export, 'Hallmark_pathways.xlsx'))





#Canonical gene sets
pathways.canoical <- gmtPathways("gsea_pathways/m2.cp.v2022.1.Mm.entrez.gmt") 
rm(canonical_data_total)
canonical_data_total <- data.frame()
for (i in SubCellTypes_no_macrophage){
  pathway_analysis <- pathway_data %>% filter(cluster == paste0(i))
  DEgenes <- pathway_analysis %>% filter(p_val < 0.05) %>% pull('gene')
  DEgenes <- mapIds(org.Mm.eg.db, DEgenes, 'ENTREZID', 'SYMBOL')
  DEG_Log2FC <-  pathway_analysis %>% filter(p_val < 0.05) %>%  pull('avg_log2FC')
  names(DEG_Log2FC) = DEgenes
  DEG_Log2FC <- sort(DEG_Log2FC, decreasing = TRUE)
  fgseaRes <- fgsea(pathways=pathways.canoical, stats=DEG_Log2FC, scoreType ='pos')
  fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  filter(pval<0.05)
  write.xlsx(fgseaResTidy, paste0('results/Hallmark_pathways_',i, 'Hallmark_pathways.xlsx'))
  fgseaResTidy$log_p <- -log10(fgseaResTidy$pval)
  fgseaResTidy$group <- paste0(i)
  canonical_data_total <- rbind(canonical_data_total, fgseaResTidy)
}
write.xlsx(canonical_data_total, paste0('results/canonical_data_total.xlsx'))
canonical_selected_pathways <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 2)%>% pull(pathway)
canonical_data_selected <- canonical_data_total %>% filter(pathway %in% canonical_selected_pathways)
write.xlsx(canonical_data_selected_pval_mt, paste0('results/canonical_data_selected_pval_mt','.xlsx'))

canonical_data_selected$group <- factor(canonical_data_selected$group,
                                          levels=c( "M0",
                                                    "HSP-Intermediate",
                                                    "Cytokine-Intermediate",
                                                    "Interferon",
                                                    "Ribosome-Intermediate",
                                                    "Cycling-S",
                                                    "Cycling-G2M",
                                                    "MGnD",
                                                    "MGnD-Antigen"))

dotplot_GO <- ggplot(canonical_data_selected,  #[c(2,3,4,9,32,37,61,76,78,97),] you can replace the numbers to the row number of pathway of your interest
                     aes(x = group, y = reorder(pathway,group))) + 
  geom_point(aes(size = log_p, color = NES)) +
  theme(text = element_text(size = pointSize, colour = "black"),
        rect = element_blank(),
        line = element_line(size = lineWidth, colour = "black"),
        plot.title  = element_text(color="black", size=20),
        axis.title.x =  element_text(size = pointSize , colour = "black"),
        axis.text.x  = element_text(size = pointSize , colour = "black", angle = 90),
        axis.text.y  = element_text(size = pointSize , colour = "black"),
        axis.line = element_line(size = lineWidth, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = pointSize , colour = "black"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"))+
  scale_colour_gradient2(limits=c(0,3), low="navyblue", mid="antiquewhite", high = "red") +
  labs(title = str_wrap("GO pathway enrichment", 20))
dotplot_GO

ggsave(paste0('plots/doplot_','canonical_selected', '_' ,'SubCellType','.pdf'),
       dotplot_GO,
       units = c("mm"),
       width = 400,
       height = 300,
       dpi = 300)

## Heatmap
canonical_selected_pathways <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 2)%>% pull(pathway)
canonical_data_selected <- canonical_data_total %>% filter(pathway %in% canonical_selected_pathways)

canonical_data_selected_pval <- canonical_data_selected %>% select(c('pathway','group', 'log_p'))
canonical_data_selected_NES <- canonical_data_selected %>% select(c('pathway','group', 'NES'))


canonical_data_selected_pval_mt <- canonical_data_selected_pval %>% 
  as_tibble() %>% pivot_wider(id_cols = pathway,
                              names_from = group,
                              values_from = log_p,
                              values_fill = 0)%>%
  column_to_rownames('pathway')

canonical_data_selected_NES_mt <- canonical_data_selected_NES %>% 
  as_tibble() %>% pivot_wider(id_cols = pathway,
                              names_from = group,
                              values_from = NES,
                              values_fill = 0)%>%
  column_to_rownames('pathway')

#Change Column order
mt_colum_order <- c( "M0","Ribosome-Intermediate","HSP-Intermediate","Cytokine-Intermediate","Interferon","Cycling-G2M","Cycling-S","MGnD","MGnD-Antigen")
canonical_data_selected_NES_mt <- canonical_data_selected_NES_mt[,mt_colum_order]
canonical_data_selected_pval_mt <- canonical_data_selected_pval_mt[,mt_colum_order]

mt_row_order <- read.xlsx(paste0('results/canonical_data_total.xlsx'), sheet = 3)%>% pull(ordered_pathways)
canonical_data_selected_NES_mt <- canonical_data_selected_NES_mt[mt_row_order,]
canonical_data_selected_pval_mt <- canonical_data_selected_pval_mt[mt_row_order,]

#Edit rownames
rownames(canonical_data_selected_pval_mt) <- sub('REACTOME_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- sub('BIOCARTA_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- sub('WP_','',rownames(canonical_data_selected_pval_mt))
rownames(canonical_data_selected_pval_mt) <- gsub('_',' ',rownames(canonical_data_selected_pval_mt))

rownames(canonical_data_selected_NES_mt) <- sub('REACTOME_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- sub('BIOCARTA_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- sub('WP_','',rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_NES_mt) <- gsub('_',' ',rownames(canonical_data_selected_NES_mt))

rownames(canonical_data_selected_NES_mt) <- tolower(rownames(canonical_data_selected_NES_mt))
rownames(canonical_data_selected_pval_mt) <- tolower(rownames(canonical_data_selected_pval_mt))

canonical_data_selected_pval_mt <- as.matrix(canonical_data_selected_pval_mt)
canonical_data_selected_NES_mt <- as.matrix(canonical_data_selected_NES_mt)




#Set up legend
col_fun = colorRamp2(c(1.3, 2, 2.6), c("white", "orange", "red"))
col_fun(seq(1.3, 2.6))
lgd = Legend(col_fun = col_fun, title = "NES")

htmp_pval <- Heatmap(canonical_data_selected_NES_mt,
        col = col_fun,
        rect_gp = gpar(col = "black", lwd = 0.2),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        width = ncol(canonical_data_selected_pval_mt)*unit(5, "mm"), 
        height = nrow(canonical_data_selected_pval_mt)*unit(5, "mm"),
        show_heatmap_legend = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(canonical_data_selected_pval_mt[i, j] > 2) {
            grid.text("***", x, y)
          } else if(canonical_data_selected_pval_mt[i, j] > 1.6) {
            grid.text("**", x, y)
          } else if(canonical_data_selected_pval_mt[i, j] > 1.3 ) {
            grid.text("*", x, y)
          }else if(canonical_data_selected_pval_mt[i, j] < 1.3) {
            grid.text("", x, y)
          }
        }) 



pdf(file = paste0('plots/heatmap_canonical_pathways_', 'SubCellType', '.pdf'), width = 10, height = 12)
htmp_pval
draw(lgd, x = unit(0.5, "cm"), just = c("left"))
dev.off()





