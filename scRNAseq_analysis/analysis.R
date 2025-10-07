library(Seurat)
library(dplyr)
library(patchwork)
library(readxl)
library(stringr)
library(ggplot2)
library(gridExtra)
library(cowplot)

# https://satijalab.org/seurat/articles/integration_introduction.html

setwd("/Users/michaelcheng/Desktop/Yang_Lab/Ulmert_collaboration/scRNAseq_analysis")
# modify R memory allocation
# library(usethis) 
# usethis::edit_r_environ()
# set memory limit higher
options(future.globals.maxSize = 10000 * 1024^2)

data_dir <- '../ULMERT LAB_LRRC15 SINGLE CELL RNASEQ/'
datasets <- list.files(data_dir, pattern="bc_matrix")
metadata <- read_excel("../ULMERT LAB_LRRC15 SINGLE CELL RNASEQ/Cell line characteristics.xlsx")
seurat_obj_list <- c()


# merge objects and process each layer individually
# each layer in the merged object is one seurat object
read_and_combine_obj <- function(){
  for(f in datasets){
    sample_pref <- gsub('_filtered_feature_bc_matrix','',f)
    counts_data <- Read10X(paste0(data_dir,'/',f))
    seurat_obj <- CreateSeuratObject(counts = counts_data, project=gsub('\\d+_(.*)','\\1',sample_pref), min.cells=100, min.feature=200)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    # nFeature_RNA_max <- min(max(seurat_obj$nFeature_RNA), as.integer(readline('nFeature_RNA max: ')))
    nFeature_RNA_max <- 5000
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < nFeature_RNA_max & percent.mt < 5)
    
    seurat_obj[['cell_line']] <- metadata %>% filter(`Sample Name`==sample_pref) %>% pull(`Cell Line Name`)
    seurat_obj$cell_line <- gsub('_|-','',seurat_obj$cell_line)
    seurat_obj[['LRRC15_surfexp']] <- metadata %>% filter(`Sample Name`==sample_pref) %>% pull(`LRRC15 surface expression`)
    seurat_obj[['TGFB_inducible']] <- metadata %>% filter(`Sample Name`==sample_pref) %>% pull(`LRRC15 can be induced by TGFB`)
    seurat_obj_list <- c(seurat_obj_list, seurat_obj)
  }
  merged_obj <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1], add.cell.ids = sapply(seurat_obj_list,function(x){unique(x$orig.ident)}), project = "all_cell_lines")
  return(merged_obj)
}
# function for splitting data and performing single cell analysis
analyze_by_group <- function(cl_obj_list){
  # unintegrated analysis
  for (i in 1:length(cl_obj_list)){
    seurat_obj = cl_obj_list[[i]]
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj)
    # KNN and UMAP
    ElbowPlot(seurat_obj, ndims=40)
    n_pcs <- 20
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, reduction='pca', graph.name="knn")
    seurat_obj <- FindClusters(seurat_obj, resolution = 1, algorithm=1, graph.name = "knn", cluster.name = "clusters")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, reduction="pca", reduction.name='umap')
    cl_obj_list[[i]] <- seurat_obj
  }
  return(cl_obj_list)
}
# integrate each group separately
# - integration method: CCAIntegration, HarmonyIntegration, JointPCAIntegration, RPCAIntegration
integrate_by_group <- function(cl_obj_list, int_method, method_name, already_split = T){
  for (i in 1:length(cl_obj_list)){
    seurat_obj <- cl_obj_list[[i]]
    # split counts into batches if not already done
    if (!already_split){
      seurat_obj[['RNA']] <- split(seurat_obj[['RNA']], f=seurat_obj$orig.ident)
    }
    # integrate using cca
    seurat_obj <- IntegrateLayers(object = seurat_obj, method = int_method, orig.reduction = "pca", new.reduction = paste0(method_name,".pca"),
                                  verbose = FALSE)
    # re-join layers after integration
    seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
    n_pcs <- 20
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, reduction=paste0(method_name,".pca"), graph.name = paste0(method_name,"_nn"))
    seurat_obj <- FindClusters(seurat_obj, resolution = 1, algorithm=1, graph.name=paste0(method_name,"_nn"), cluster.name = paste0(method_name,"_clusters"))
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, reduction = paste0(method_name,".pca"), reduction.name = paste0(method_name,".umap"))
    cl_obj_list[[i]] <- seurat_obj
    rm(seurat_obj)
  }
  return(cl_obj_list)
}


##########################################################
# cell-line specific analysis
merged_obj <- read_and_combine_obj()
cl_obj_list <- SplitObject(merged_obj, split.by="cell_line")
rm(merged_obj)

cl_obj_list <- analyze_by_group(cl_obj_list)
plot_list <- list()
for (i in 1:length(cl_obj_list)){
  # print(DimPlot(cl_obj_list[[i]], group.by=c('LRRC15_surfexp')) + xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp')))
  plot_list[[i]] <- DimPlot(cl_obj_list[[i]], group.by=c('LRRC15_surfexp')) + xlab('UMAP1') + ylab('UMAP2') + 
                   ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp'))
}
# plot in grid
plot_grid(plotlist=plot_list, nrow=2, ncol=3)

# integrate each cell line
cl_obj_list <- integrate_by_group(cl_obj_list, HarmonyIntegration, "harmony")
# save list object
saveRDS(cl_obj_list, "./object_list_by_cell_line.RDS")

dir.create("umaps")
plot_list <- list()
for (i in 1:length(cl_obj_list)){
  plot_list[[i]] <- (DimPlot(cl_obj_list[[i]], group.by=c('LRRC15_surfexp')) + xlab('UMAP1') + ylab('UMAP2') + ggtitle('uncorrected')) +
    (DimPlot(cl_obj_list[[i]], reduction =paste0("harmony",".umap"), group.by=c('LRRC15_surfexp')) + xlab('UMAP1') + ylab('UMAP2') + ggtitle("harmony"))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
}
p <- plot_grid(plotlist=plot_list, ncol=3, labels = paste0(names(cl_obj_list), ' LRRC15_surfexp'), label_x = -0.1, label_y=1.01)
save_plot("umaps/cell_line_LRRC_surfexp_umap.pdf", p, base_width = 20, base_height=10)

plot_list <- list()
for (i in 1:length(cl_obj_list)){
  plot_list[[i]] <- FeaturePlot(cl_obj_list[[i]], features=c('LRRC15','TGFB1'),combine = T, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
}
p <- plot_grid(plotlist=plot_list, ncol=3, labels = names(cl_obj_list), label_y=1.01)
save_plot("umaps/cell_line_LRRC15_TGFB1_featureplot.pdf", p, base_width = 20, base_height=10)

plot_list <- list()
for (i in 1:length(cl_obj_list)){
  plot_list[[i]] <- FeaturePlot(cl_obj_list[[i]], reduction=paste0("harmony",".umap"), features=c('LRRC15','TGFB1'),combine = T, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
}
p <- plot_grid(plotlist=plot_list, ncol=3, labels = names(cl_obj_list), label_y=1.01)
save_plot("umaps/cell_line_LRRC15_TGFB1_featureplot_harmony.pdf", p, base_width = 20, base_height=10)


##########################################################
# LRRC15 surface expression analysis
merged_obj <- read_and_combine_obj()
cl_obj_list <- SplitObject(merged_obj, split.by="LRRC15_surfexp")
rm(merged_obj)

cl_obj_list <- analyze_by_group(cl_obj_list)
# integrate each cell line
cl_obj_list <- integrate_by_group(cl_obj_list, HarmonyIntegration, "harmony")
# save list object
saveRDS(cl_obj_list, "./object_list_by_LRRC15_surfexp.RDS")

plot_list <- list()
plot_list2 <- list()
for (i in 1:length(cl_obj_list)){
  plot_list[[i]] <- DimPlot(cl_obj_list[[i]], group.by=c('cell_line', 'TGFB_inducible')) + xlab('UMAP1') + ylab('UMAP2')
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$LRRC15_surfexp), ' LRRC15_surfexp_harmony_umap.png'), plot=p, width=15, height=10,units='in')
  plot_list2[[i]] <- DimPlot(cl_obj_list[[i]], reduction=paste0("harmony",".umap"), group.by=c('cell_line', 'TGFB_inducible')) + xlab('UMAP1') + ylab('UMAP2')
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$LRRC15_surfexp), ' LRRC15_surfexp_harmony_umap.png'), plot=p, width=15, height=10,units='in')
}
p <- plot_grid(plotlist=c(plot_list,plot_list2), ncol=2, labels = paste0(names(cl_obj_list), ' LRRC15 Surfexp'), label_y=1.02, label_x=-0.01)
save_plot("umaps/LRRC15_surfexp_umap.pdf", p, base_width = 20, base_height=10)


plot_list3 <- list()
plot_list4 <- list()
for (i in 1:length(cl_obj_list)){
  plot_list3[[i]] <- FeaturePlot(cl_obj_list[[i]], features=c('LRRC15','TGFB1'),combine = T, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
  plot_list4[[i]] <- FeaturePlot(cl_obj_list[[i]], reduction=paste0("harmony",".umap"), features=c('LRRC15','TGFB1'),combine = T, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
}
p <- plot_grid(plotlist=c(plot_list3,plot_list4), ncol=2, labels = paste0(names(cl_obj_list), ' LRRC15 Surfexp'), label_y=1.02, label_x=-0.01)
save_plot("umaps/LRRC15_surfexp_LRRC15_TGFB1_featureplot.pdf", p, base_width = 20, base_height=10)

p <- plot_grid(plotlist=c(plot_list,plot_list3,plot_list2,plot_list4), ncol=2, labels = paste0(names(cl_obj_list), ' LRRC15 Surfexp'), label_y=1.02, label_x=-0.01)
save_plot("umaps/LRRC15_surfexp_umap_featureplots.pdf", p, base_width = 20, base_height=20)
##########################################################

##########################################################
# TGFB inducibility analysis
merged_obj <- read_and_combine_obj()
cl_obj_list <- SplitObject(merged_obj, split.by="TGFB_inducible")
rm(merged_obj)

cl_obj_list <- analyze_by_group(cl_obj_list)
# integrate each cell line
cl_obj_list <- integrate_by_group(cl_obj_list, HarmonyIntegration, "harmony")
# save list object
saveRDS(cl_obj_list, "./object_list_by_TGFB_inducible.RDS")

dir.create("umaps")

plot_list <- list()
plot_list2 <- list()
for (i in 1:length(cl_obj_list)){
  plot_list[[i]] <- DimPlot(cl_obj_list[[i]], group.by=c('cell_line', 'LRRC15_surfexp')) + xlab('UMAP1') + ylab('UMAP2')
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$LRRC15_surfexp), ' LRRC15_surfexp_harmony_umap.png'), plot=p, width=15, height=10,units='in')
  plot_list2[[i]] <- DimPlot(cl_obj_list[[i]], reduction=paste0("harmony",".umap"), group.by=c('cell_line', 'LRRC15_surfexp')) + xlab('UMAP1') + ylab('UMAP2')
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$LRRC15_surfexp), ' LRRC15_surfexp_harmony_umap.png'), plot=p, width=15, height=10,units='in')
}
p <- plot_grid(plotlist=c(plot_list,plot_list2), ncol=2, labels = paste0(names(cl_obj_list), ' TGFB Inducible'), label_y=1.02, label_x=-0.01)
save_plot("umaps/TGFB_inducible_umap.pdf", p, base_width = 20, base_height=10)

plot_list3 <- list()
plot_list4 <- list()
for (i in 1:length(cl_obj_list)){
  plot_list3[[i]] <- FeaturePlot(cl_obj_list[[i]], features=c('LRRC15','TGFB1'),combine = T, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
  plot_list4[[i]] <- FeaturePlot(cl_obj_list[[i]], reduction=paste0("harmony",".umap"), features=c('LRRC15','TGFB1'),combine = T, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
  # ggsave(paste0("umaps/",unique(cl_obj_list[[i]]$cell_line), ' LRRC15_surfexp_harmony_umap.png'),width=15, height=10,units='in')
}
p <- plot_grid(plotlist=c(plot_list3,plot_list4), ncol=2, labels = paste0(names(cl_obj_list), ' TGFB Inducible'), label_y=1.02, label_x=-0.01)
save_plot("umaps/TGFB_inducible_LRRC15_TGFB1_featureplot.pdf", p, base_width = 20, base_height=10)

p <- plot_grid(plotlist=c(plot_list,plot_list3,plot_list2,plot_list4), ncol=2, labels = paste0(names(cl_obj_list), ' TGFB Inducible'), label_y=1.02, label_x=-0.01)
save_plot("umaps/TGFB_inducible_umap_featureplots.pdf", p, base_width = 20, base_height=20)
##########################################################
rm(cl_obj_list)

##########################################################
## Object used for downstream analysis
merged_obj <- read_and_combine_obj()
cl_obj_list <- list(merged_obj)
rm(merged_obj)
# unintegrated analysis
cl_obj_list <- analyze_by_group(cl_obj_list)
# integrate each cell line
cl_obj_list <- integrate_by_group(cl_obj_list, HarmonyIntegration, "harmony")
saveRDS(cl_obj_list, "object_list_all.RDS")

cl_obj_list[[1]][['RNA']] <- as(cl_obj_list[[1]][['RNA']], Class = "Assay")
saveRDS(cl_obj_list[[1]], 'merged_object_v4.RDS')

plot_list <- DimPlot(cl_obj_list[[1]], group.by=c('cell_line', 'LRRC15_surfexp', 'TGFB_inducible'),combine = F)
plot_list2 <- DimPlot(cl_obj_list[[1]], reduction=paste0("harmony",".umap"), group.by=c('cell_line', 'LRRC15_surfexp', 'TGFB_inducible'), combine=F)

p <- plot_grid(plotlist=c(plot_list,plot_list2), nrow=3, byrow=F, labels = c('UMAP','Harmony'), label_y=1.02, label_x=-0.01)
save_plot("umaps/all_umap.pdf", p, base_width = 20, base_height=20)

plot_list3 <- FeaturePlot(cl_obj_list[[1]], features=c('LRRC15','TGFB1'),combine = F, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
plot_list4 <- FeaturePlot(cl_obj_list[[1]], reduction=paste0("harmony",".umap"), features=c('LRRC15','TGFB1'),combine = F, keep.scale = "all") #+ xlab('UMAP1') + ylab('UMAP2') + ggtitle(paste0(unique(cl_obj_list[[i]]$cell_line), ' LRRC15'))
p <- plot_grid(plotlist=c(plot_list3,plot_list4), ncol=2, byrow=F, labels = c('UMAP','Harmony'), label_y=1.02, label_x=-0.01)
save_plot("umaps/all_LRRC15_TGFB1_featureplot.pdf", p, base_width = 20, base_height=10)

p <- plot_grid(plotlist=c(plot_list,plot_list3,plot_list2,plot_list4), nrow=5, byrow=F, labels = c('UMAP','Harmony'), label_y=1.02, label_x=-0.01)
save_plot("umaps/all_umap_featureplots.pdf", p, base_width = 20, base_height=20)




##########################################################
# Old script
merged_obj <- read_and_combine_obj()
# unintegrated analysis

merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
# KNN and UMAP
ElbowPlot(merged_obj, ndims=40)
n_pcs <- 20
merged_obj <- FindNeighbors(merged_obj, dims = 1:n_pcs, reduction='pca', graph.name="unintegrated_nn")
merged_obj <- FindClusters(merged_obj, resolution = 1, algorithm=1, graph.name = "unintegrated_nn", cluster.name = "unintegrated_clusters")
merged_obj <- RunUMAP(merged_obj, dims = 1:n_pcs, reduction="pca", reduction.name='umap.unintegrated')
DimPlot(merged_obj, reduction = "umap.unintegrated")
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = "orig.ident")
DimPlot(merged_obj, reduction = "umap.unintegrated", group.by = c("cell_line", "LRRC15_surfexp", "TGFB_inducible"))


# integrate using rpca
merged_obj <- IntegrateLayers(object = merged_obj, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                              verbose = FALSE)
# now that integration is complete, rejoin layers
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
n_pcs <- 20
merged_obj <- FindNeighbors(merged_obj, dims = 1:n_pcs, reduction='integrated.rpca', graph.name = "rpca_nn")
merged_obj <- FindClusters(merged_obj, resolution = 1, algorithm=1, graph.name="rpca_nn", cluster.name = "rpca_clusters")
merged_obj <- RunUMAP(merged_obj, dims = 1:n_pcs, reduction = 'integrated.rpca', reduction.name = "umap.rpca")
DimPlot(merged_obj, reduction = "umap.rpca")
DimPlot(merged_obj, reduction = "umap.rpca", group.by = "orig.ident")
DimPlot(merged_obj, reduction = "umap.rpca", group.by = c("cell_line", "LRRC15_surfexp", "TGFB_inducible"))



# integrate using cca
merged_obj[['RNA']] <- split(merged_obj[['RNA']], f=merged_obj$orig.ident)
merged_obj <- IntegrateLayers(object = merged_obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE)
# re-join layers after integration
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
n_pcs <- 20
merged_obj <- FindNeighbors(merged_obj, dims = 1:n_pcs, reduction='integrated.cca', graph.name = "cca_nn")
merged_obj <- FindClusters(merged_obj, resolution = 1, algorithm=1, graph.name="cca_nn", cluster.name = "cca_clusters")
merged_obj <- RunUMAP(merged_obj, dims = 1:n_pcs, reduction = 'integrated.cca', reduction.name = "umap.cca")
DimPlot(merged_obj, reduction = "umap.cca")
DimPlot(merged_obj, reduction = "umap.cca", group.by = "orig.ident")
DimPlot(merged_obj, reduction = "umap.cca", group.by = c("cell_line", "LRRC15_surfexp", "TGFB_inducible"))

# saveRDS(merged_obj, "merged_object.RDS")

merged_obj <- readRDS("merged_object.RDS")
# find DEGs between conditions
Idents(merged_obj) <- "TGFB_inducible"
degs <- FindMarkers(merged_obj, ident.1 = "yes", ident.2 = "no", verbose = T)

degs %>% filter(pct.1>0 | pct.2>0)

# aggregate cell lines
agg_merged_obj <- AggregateExpression(merged_obj, group.by = c("orig.ident", "LRRC15_surfexp", "TGFB_inducible"), return.seurat = TRUE)
Idents(agg_merged_obj) <- "TGFB_inducible"
degs <- FindMarkers(agg_merged_obj, ident.1 = "yes", ident.2 = "no", verbose = T)


# markers for each cell lines
Idents(merged_obj) <- "cell_line"
cl.markers <- FindConservedMarkers(merged_obj, ident.1 = "CALU1", grouping.var = "LRRC15_surfexp")
head(nk.markers)

saveRDS(merged_obj, "merged_object.RDS")
# save as a seurat version 3 or 4 object
merged_obj[['RNA']] <- as(merged_obj[['RNA']], Class = "Assay")
saveRDS(merged_obj, 'merged_object_v4.RDS')
