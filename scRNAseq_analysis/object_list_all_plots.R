library(Seurat)
library(dplyr)
library(patchwork)
library(readxl)
library(stringr)
library(ggplot2)



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
    # nFeature_RNA_max <- 5000
    # seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < nFeature_RNA_max & percent.mt < 5)
    
    seurat_obj[['cell_line']] <- metadata %>% filter(`Sample Name`==sample_pref) %>% pull(`Cell Line Name`)
    seurat_obj$cell_line <- gsub('_|-','',seurat_obj$cell_line)
    seurat_obj[['LRRC15_surfexp']] <- metadata %>% filter(`Sample Name`==sample_pref) %>% pull(`LRRC15 surface expression`)
    seurat_obj[['TGFB_inducible']] <- metadata %>% filter(`Sample Name`==sample_pref) %>% pull(`LRRC15 can be induced by TGFB`)
    seurat_obj_list <- c(seurat_obj_list, seurat_obj)
  }
  merged_obj <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1], add.cell.ids = sapply(seurat_obj_list,function(x){unique(x$orig.ident)}), project = "all_cell_lines")
  return(merged_obj)
}

seurat_obj <- read_and_combine_obj()

# qc plots by sample
seurat_obj$orig.ident <- as.factor(seurat_obj$orig.ident)
Idents(seurat_obj) <- 'orig.ident'
levels(Idents(seurat_obj)) 
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=F, alpha=0.5)
g <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), raster=F, combine = F)
g <- lapply(g, function(x){x$layers[[2]] <- NULL; x <- x+geom_boxplot(); return(x)})

g_combined <- wrap_plots(g) + plot_layout(axes = "collect", guides = "collect", nrow=3)
dir.create('./qcplots', showWarnings = F)
ggsave('./qcplots/qc_boxplots.pdf', g_combined, width=10, height = 15, units = "in")

f1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = F, combine = F) + facet_wrap(.~colors) + 
  ggtitle('') + theme(axis.text.x = element_text(angle=45, vjust = 0.5))
f2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = F, combine = F) + facet_wrap(.~colors) + 
  ggtitle('') + theme(axis.text.x = element_text(angle=45, vjust = 0.5))

p <- f1/f2 + plot_layout(nrow=2, axes="collect", guides="collect")
ggsave('./qcplots/feature_scatters.pdf', p, width=10, height = 10, units = "in")



# qc plots filtered
seurat_obj <- readRDS("object_list_all.RDS")[[1]]
seurat_obj@meta.data$orig.ident <- factor(seurat_obj@meta.data$orig.ident)
levels(seurat_obj@meta.data$orig.ident)
Idents(seurat_obj) <- 'orig.ident'
levels(Idents(seurat_obj))
g <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), raster=F, combine = F)
g <- lapply(g, function(x){x$layers[[2]] <- NULL; x <- x+geom_boxplot(); return(x)})

g_combined <- wrap_plots(g) + plot_layout(axes = "collect", guides = "collect", nrow=3)
ggsave('./qcplots/qc_boxplots_filtered.pdf', g_combined, width=10, height = 15, units = "in")

