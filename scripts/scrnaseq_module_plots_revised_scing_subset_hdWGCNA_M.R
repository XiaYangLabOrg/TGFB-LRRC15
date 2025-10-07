library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
# library(readxl)
library(stringr)
library(ggplot2)
# library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(WGCNA)
library(hdWGCNA)
library(scCustomize)
library(parallel)
library(PCAtools)

setwd("~/Desktop/Yang_Lab/Ulmert_collaboration/scripts")
### hdwgcna metacell functions
# create seurat object for module expression
create_seurat_modexp <- function(mod_exp, seurat_obj, pca_embed=NULL){
  # make sure cells in mod_exp and seurat object overlap
  if(sum(rownames(mod_exp)%in%Cells(seurat_obj))==0){stop("Cell names have 0 match")}
  # create new seurat object
  mod_obj <- CreateSeuratObject(counts=t(mod_exp),data=t(mod_exp))
  mod_obj[['RNA']]$scale.data <- mod_obj[['RNA']]$counts
  mod_obj <- FindVariableFeatures(mod_obj)
  if(!is.null(pca_embed)){
    mod_obj[['pca']] <- CreateDimReducObject(embeddings = as.matrix(pca_embed),assay=DefaultAssay(seurat_obj)) 
  }
  # add umap 
  mod_obj[['umap']] <- subset(seurat_obj, cells = Cells(mod_obj))[['umap']]
  
  mod_obj@meta.data <- merge(mod_obj@meta.data, subset(seurat_obj, cells = Cells(mod_obj))@meta.data, by=0)
  rownames(mod_obj@meta.data) <- mod_obj@meta.data$Row.names
  mod_obj@meta.data <- mod_obj@meta.data[,-1]
  return(mod_obj)
}

setup_wgcna <- function(seurat_obj, wgcna_name='wgcna'){
  # Set Up Seurat Object for WGCNA
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "custom", # the gene selection approach
    gene_list = Features(seurat_obj), # list of genes to include
    wgcna_name = wgcna_name # the name of the hdWGCNA experiment
  )
  return(seurat_obj)
}
get_metacells <- function(seurat_obj, celltype_col, sample_col, metacell_size=25, max_metacell_overlap=10, n_metacell=20){
  # construct metacells in each group
  seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c(celltype_col, sample_col), # specify the columns in seurat_obj@meta.data to group by
    reduction = 'pca', # select the dimensionality reduction to perform KNN on
    target_metacells = n_metacell,
    k = metacell_size, # nearest-neighbors parameter
    max_shared = max_metacell_overlap, # maximum number of shared cells between two metacells
    ident.group = celltype_col # set the Idents of the metacell seurat object
  )
  # seurat_obj[['RNA']]$data <- seurat_obj[['RNA']]$counts
  # normalize metacell expression matrix:
  # seurat_obj <- NormalizeMetacells(seurat_obj)
  # saveRDS(seurat_obj, paste0(outdir,'/',outpref,'_metacells.rds'))
  return(seurat_obj)
}
###

### Sample-wise DEG
# function for statistical test between group
test_groups <- function(x, group, test_function){
  res <- test_function(x[group],x[!group])
  avg_diff <- mean(x[group]) - mean(x[!group])
  med_diff <- median(x[group]) - median(x[!group])
  return(c(statistic = unname(res$statistic), p_val=unname(res$p.value),avg_diff=avg_diff, med_diff=med_diff))
}
# function to run t-test based on group
run_test <- function(df, groups, test_function,column_prefix=''){
  system.time({
    cl <- makeCluster(detectCores())
    res <- mclapply(1:ncol(df), FUN=function(i){ test_groups(df[,i],group=groups, test_function=test_function)})
    res <- as.data.frame(do.call(rbind,res))
    stopCluster(cl)
  })
  rownames(res) <- colnames(df)
  res$p_val_adj <- p.adjust(res$p, method='BH')
  colnames(res) <- paste0(column_prefix, colnames(res))
  return(res)
}

sample_wise_degs <- function(seurat_obj, ident, ident.1, ident.2,slot='scale.data', fc.slot='scale.data', method=NULL){
  Idents(seurat_obj) <- ident
  seurat_obj$tmp_sample_group <- NA
  seurat_obj$tmp_sample_group[seurat_obj@meta.data[,ident]%in%ident.1] <- 0
  seurat_obj$tmp_sample_group[seurat_obj@meta.data[,ident]%in%ident.2] <- 1
  # print(table(seurat_obj$tmp_sample_group, exclude=F))
  Idents(seurat_obj) <- "tmp_sample_group"
  if (class(method)=="function"){
    test_obj <- subset(seurat_obj, cells = rownames(seurat_obj@meta.data)[!is.na(seurat_obj$tmp_sample_group)])
    sample_degs <- run_test(t(as.data.frame(test_obj[['RNA']]$scale.data)), groups=test_obj$tmp_sample_group==0, test_function=method)
  } else {
    sample_degs <- FindMarkers(seurat_obj, ident.1=0, ident.2=1, slot=slot, fc.slot=fc.slot)
  }
  return(sample_degs)
}
run_samplewise_degs <- function(seurat_obj, ident, condition, method=t.test){
  # markers for each cell lines
  Idents(seurat_obj) <- ident
  s_degs <- data.frame()
  # identify conserved DEGs across cell lines
  for(cl in unique(seurat_obj@meta.data[[ident]])){
    # get cells belonging to opposite condition
    ident2 <- seurat_obj@meta.data %>% 
      filter(!!as.name(condition) != unique(seurat_obj@meta.data %>% filter(!!as.name(ident)==cl) %>% pull(!!as.name(condition)))) %>%
      pull(!!as.name(ident)) %>% unique()
    # run samplewise_degs
    cl.markers <- sample_wise_degs(seurat_obj, ident, ident.1 = cl, ident.2=ident2, method=method)
    cl.markers$cell_line <- cl
    cl.markers <- cl.markers %>% rownames_to_column('gene')
    s_degs <- rbind(s_degs, cl.markers)
  }
  return(s_degs)
}
get_common_degs <- function(sig_degs, cell_line_cutoff, same_direction=F){
  common_degs <- table(sig_degs$gene)
  common_degs <- names(common_degs[common_degs>=cell_line_cutoff])
  if(same_direction){
    same_dir_sig_genes <- sapply(common_degs, function(x){length(unique(sign(sig_degs[sig_degs$gene==x,"avg_diff"])))})
    same_dir_sig_genes <- names(same_dir_sig_genes[same_dir_sig_genes==1])
    return(same_dir_sig_genes)
  }
  return(common_degs)
}
######################################################
### create module expression objects
# original scrnaseq object
seurat_obj <- readRDS("../scRNAseq_analysis/object_list_all.RDS")
seurat_obj <- seurat_obj[[1]]
seurat_obj <- RenameCells(seurat_obj, new.names = gsub('_','',Cells(seurat_obj)))

# create new seurat object
hdwgcna_mod <- read.delim('../scRNAseq_analysis/hdwgcna/final_all_hdWGCNA_object_MEs.txt') # module expression
hdwgcna_obj <- create_seurat_modexp(mod_exp = hdwgcna_mod, seurat_obj = seurat_obj, pca_embed = hdwgcna_mod)
Idents(hdwgcna_obj) <- hdwgcna_obj$cell_line
# add metacell objects
hdwgcna_obj <- setup_wgcna(hdwgcna_obj, wgcna_name = 'modules')
hdwgcna_obj <- get_metacells(hdwgcna_obj, celltype_col = 'cell_line', sample_col='orig.ident.y')
hdwgcna_metacell <- GetMetacellObject(hdwgcna_obj)
hdwgcna_metacell[['RNA']]$data <- hdwgcna_metacell[['RNA']]$counts
hdwgcna_metacell[['RNA']]$scale.data <- hdwgcna_metacell[['RNA']]$counts
metadata <- hdwgcna_obj@meta.data %>% select(orig.ident.y, TGFB_inducible,LRRC15_surfexp) %>% unique() %>% remove_rownames() %>% column_to_rownames('orig.ident.y')
# add TGFB inducible metadata column
hdwgcna_metacell$TGFB_inducible <- metadata[hdwgcna_metacell$orig.ident.y,'TGFB_inducible']
hdwgcna_metacell$LRRC15_surfexp <- metadata[hdwgcna_metacell$orig.ident.y,"LRRC15_surfexp"]

# get samplewise_degs
hdwgcna_deg <- run_samplewise_degs(hdwgcna_metacell, ident='cell_line', condition='TGFB_inducible', method = wilcox.test)
# flip direction of fold change for TGFB-noninducible cell lines
tgfb_no <- hdwgcna_metacell@meta.data %>% filter(TGFB_inducible=='no') %>% pull(cell_line) %>% unique()
tgfb_yes <- hdwgcna_metacell@meta.data %>% filter(TGFB_inducible=='yes') %>% pull(cell_line) %>% unique()
hdwgcna_deg$flipped_avg_diff <- ifelse(hdwgcna_deg$cell_line %in% tgfb_no, -1*hdwgcna_deg$avg_diff, hdwgcna_deg$avg_diff)
sig_degs <- hdwgcna_deg %>% filter(p_val_adj <= 0.05)

## SCING
# read module expression
scing_mod <- read.delim('../pathway_enrichment/pathway_enrichments/scing_ME.50.txt')
rownames(scing_mod) <- gsub('_','',rownames(scing_mod))
colnames(scing_mod) <- gsub('all_','SCING_',colnames(scing_mod))
scing_mod <- scing_mod[,order(as.numeric(gsub('SCING_M','',colnames(scing_mod))))]
# subsitute scing module
scing0_subset <- read.delim("../pathway_enrichment/pathway_enrichments/scing_ME.50.subsetted.txt")
rownames(scing0_subset) <- gsub('_','',rownames(scing0_subset))
method <- 'hdWGCNA_M' # hdWGCNA_M, hdWGCNA_all, deg_0$, degree
outfile_suffix=paste0('with_SCING_',method,'_subsets')
scing0_subset <- scing0_subset[,grep(method,colnames(scing0_subset))]
scing_mod <- merge(scing_mod[,colnames(scing_mod)!='SCING_M0'], scing0_subset, by=0) %>% column_to_rownames("Row.names")
pca_scing <-  scing_mod
colnames(pca_scing) <- paste0('SCING_M',c(1:ncol(scing_mod)))
# create new seurat object
scing_obj <- create_seurat_modexp(mod_exp = scing_mod, seurat_obj = seurat_obj, pca_embed = pca_scing)
scing_obj[['RNA']]$new_scale.data <- t(scale(t(as.matrix(scing_obj[['RNA']]$scale.data))))
Idents(scing_obj) <- scing_obj$cell_line
# check degs
Idents(scing_obj) <- scing_obj$TGFB_inducible
scing_deg <- FindMarkers(scing_obj, ident.1 = 'yes',test.use='t', slot='scale.data',fc.slot='scale.data', min.pct=0)
scing_deg %>% arrange(avg_diff)

# add metacell objects
scing_obj <- setup_wgcna(scing_obj, wgcna_name = 'modules')
scing_obj <- get_metacells(scing_obj, celltype_col = 'cell_line', sample_col='orig.ident.y')
scing_metacell <- GetMetacellObject(scing_obj)
scing_metacell[['RNA']]$data <- scing_metacell[['RNA']]$counts
scing_metacell[['RNA']]$scale.data <- scale(scing_metacell[['RNA']]$counts)
scing_metacell[['RNA']]$new_scale.data <- t(scale(t(as.matrix(scing_metacell[['RNA']]$scale.data))))


# ggplot(scing_metacell_pca$x)
metadata <- scing_obj@meta.data %>% select(orig.ident.y, TGFB_inducible,LRRC15_surfexp) %>% unique() %>% remove_rownames() %>% column_to_rownames('orig.ident.y')
# add TGFB inducible metadata column
scing_metacell$TGFB_inducible <- metadata[scing_metacell$orig.ident.y,'TGFB_inducible']
scing_metacell$LRRC15_surfexp <- metadata[scing_metacell$orig.ident.y,"LRRC15_surfexp"]
# save
saveRDS(scing_metacell, paste0("../pathway_enrichment/pathway_enrichments/ME_metacells_",outfile_suffix,".rds"))
# pca
scing_metacell_pca <- pca(as.matrix(scing_metacell[['RNA']]$scale.data), metadata = scing_metacell@meta.data[,c("cell_line","TGFB_inducible","LRRC15_surfexp")], scale = T)
screeplot(scing_metacell_pca, axisLabSize = 18, titleLabSize = 22)
biplot(scing_metacell_pca, x='PC1', y='PC4',lab=NULL, colby="TGFB_inducible", legendPosition = 'right', showLoadings = T)
pairsplot(scing_metacell_pca, colby = 'TGFB_inducible')
scing_metacell[['pca']] <- CreateDimReducObject(embeddings = as.matrix(scing_metacell_pca$rotated), loadings = as.matrix(scing_metacell_pca$loadings), key='PCA')
scing_metacell <- FindNeighbors(scing_metacell, dims = 1:20, reduction='pca', graph.name="knn")
scing_metacell <- FindClusters(scing_metacell, resolution = 1, algorithm=1, graph.name = "knn", cluster.name = "clusters")
scing_metacell <- RunUMAP(scing_metacell, dims = 1:20, reduction="pca", reduction.name='umap')

# group-wise degs
Idents(scing_metacell) <- 'TGFB_inducible'
# Idents(scing_metacell) <- 'cell_line'
# levels(Idents(scing_metacell)) <- c('CALU1','KASUMI2','SAOS9', 'HuO9','RPMI7951','U118')
scing_deg <- run_test(df=t(as.data.frame(scing_metacell[['RNA']]$scale.data)), groups = scing_metacell$TGFB_inducible=='yes',test_function = t.test)
sig_degs <- scing_deg %>% filter(p_val_adj <= 0.05)
# scing_deg <- FindMarkers(scing_metacell, ident.1 = 'yes',test.use='t', slot='scale.data',fc.slot='scale.data', min.pct=0)

Idents(scing_metacell) <- 'cell_line'
Clustered_DotPlot(scing_metacell, features=rownames(sig_degs),colors_use_exp = colorRampPalette(c('#377EB8','white','#E41A1C'))(20), ggplot_default_colors = T,row_label_size = 15,)
DotPlot(scing_metacell, features = rownames(sig_degs), cols = 'RdBu', scale.by = 'size')


# get samplewise_degs
scing_deg <- run_samplewise_degs(scing_metacell, ident='cell_line', condition='TGFB_inducible', method = wilcox.test)
# flip direction of fold change for TGFB-noninducible cell lines
tgfb_no <- scing_metacell@meta.data %>% filter(TGFB_inducible=='no') %>% pull(cell_line) %>% unique()
tgfb_yes <- scing_metacell@meta.data %>% filter(TGFB_inducible=='yes') %>% pull(cell_line) %>% unique()
scing_deg$flipped_avg_diff <- ifelse(scing_deg$cell_line %in% tgfb_no, -1*scing_deg$avg_diff, scing_deg$avg_diff)
scing_deg$flipped_med_diff <- ifelse(scing_deg$cell_line %in% tgfb_no, -1*scing_deg$med_diff, scing_deg$med_diff)
sig_degs <- scing_deg %>% filter(p_val_adj <= 0.05)

common_degs <- get_common_degs(sig_degs, length(unique(sig_degs$cell_line)))
frequent_degs <- get_common_degs(sig_degs, length(unique(sig_degs$cell_line))-1)
tgfb_yes_degs <- get_common_degs(sig_degs[sig_degs$cell_line%in%tgfb_yes,], length(tgfb_yes)-1,same_direction = F)
tgfb_yes_degs_dir <- get_common_degs(sig_degs[sig_degs$cell_line%in%tgfb_yes,], length(tgfb_yes),same_direction = T)
tgfb_no_degs <- get_common_degs(sig_degs[sig_degs$cell_line%in%tgfb_no,], length(tgfb_no)-1, same_direction = F)
tgfb_no_degs_dir <- get_common_degs(sig_degs[sig_degs$cell_line%in%tgfb_no,], length(tgfb_no), same_direction = T)
common_degs_dir <- intersect(tgfb_yes_degs_dir, tgfb_no_degs_dir)

sig_degs[sig_degs$gene %in% common_degs_dir,] %>% arrange(gene)

# levels(Idents(scing_metacell)) <- c('CALU1','KASUMI2','SAOS9', 'HuO9','RPMI7951','U118')
# scing_metacell$cell_line <- factor(scing_metacell$cell_line, levels=c('CALU1','KASUMI2','SAOS9', 'HuO9','RPMI7951','U118'))
Idents(scing_metacell) <- 'cell_line'
DotPlot(scing_metacell, features = tgfb_yes_degs, cols = 'RdBu')
Idents(scing_metacell) <- 'cell_line'
Clustered_DotPlot(scing_metacell, features=common_degs_dir %>% unique(),colors_use_exp = colorRampPalette(c('#377EB8','white','#E41A1C'))(20), ggplot_default_colors = T,row_label_size = 15,)

p <- FeaturePlot(scing_metacell, features=c(common_degs_dir,'SCING-M0-hdWGCNA-M2','SCING-M0-hdWGCNA-M18'), combine = T,label = F, reduction = 'pca', split.by = 'TGFB_inducible', slot = 'new_scale.data') &
  theme(legend.position="right")&
  # scale_color_gradient2(low="blue",mid="white", high = "red")
  scale_colour_gradientn(limits=c(-4,4),colors = rev(brewer.pal(n = 11, name = "RdBu")))
save_plot(paste0("../scRNAseq_analysis/module_umaps/SCING_module_umaps_metacells_",outfile_suffix,".pdf"), p, base_width = 20, base_height=20)

Idents(scing_obj) <- 'cell_line'
p <- FeaturePlot(scing_obj, features=c(common_degs_dir), combine = T, split.by = 'TGFB_inducible',label = F, slot = 'scale.data', pt.size = 1) +
  plot_layout(axis= "collect", guides="collect") & xlab('UMAP1') & ylab('UMAP2')& ggtitle(NULL)&
  theme(legend.position="right", axis.text = element_text(size = 20), axis.title = element_text(size=20)) &
  scale_color_gradient2(low="#377EB8",mid="white", high = "#E41A1C",limits=c(-4,4), oob = scales::squish)
save_plot(paste0("../scRNAseq_analysis/module_umaps/SCING_module_umaps_singlecells_",outfile_suffix,".pdf"), p, base_width = 10, base_height=10)

##### Heatmap of Samplewise DEGs
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

hdwgcna_plot_df <- hdwgcna_deg %>% select(gene, med_diff, p_val_adj, cell_line) %>% pivot_wider(names_from = cell_line, values_from = c(p_val_adj,med_diff)) %>% mutate(method='hdWGCNA')
scing_plot_df <- scing_deg %>% select(gene, med_diff, p_val_adj, cell_line) %>% pivot_wider(names_from = cell_line, values_from = c(p_val_adj,med_diff)) %>% mutate(method='SCING')

comb_plot_df <- rbind(hdwgcna_plot_df, scing_plot_df) %>% column_to_rownames('gene')
colnames(comb_plot_df) <- gsub('HuO9','HUO9',colnames(comb_plot_df))

# convert encodings all to ASCII (for some reason HUO9 and U118 are in UTF-8)
library(stringi)
stri_enc_mark(colnames(comb_plot_df))
all(stri_enc_isutf8(colnames(comb_plot_df)))
colnames(comb_plot_df) <- stri_enc_toascii(colnames(comb_plot_df))
colnames(comb_plot_df) <- gsub('\032','',colnames(comb_plot_df), fixed = T)

# median difference
med_diff_df <- comb_plot_df[,grepl('med_diff', colnames(comb_plot_df))]
colnames(med_diff_df) <- gsub('med_diff_','',colnames(med_diff_df))

pval_marker_cell_fun <- function(df){
    return(function(j, i, x, y, w, h, fill) {
      if(df[i, j] < 0.001) {
        grid.text("***", x, y, vjust = 0.7)
      } else if(df[i, j] < 0.01) {
        grid.text("**", x, y, vjust = 0.7)
      } else if(df[i, j] < 0.05){
        grid.text("*", x, y, vjust = 0.7)
      }
    })
}
text_marker_cell_fun <- function(j, i, x, y, w, h, fill) {
  sub_df <- comb_plot_df[,grepl('med_diff', colnames(comb_plot_df))]
  grid.text(format(comb_plot_df[,grepl('med_diff', colnames(comb_plot_df))][i,j], digits=3), x, y)
}
color_pal <- brewer.pal(9,"Set1")[c(2,9,1)]
color_pal[2] <- 'white'
marker_col_fun <- function(min, mid, max, color_pal){
  return(colorRamp2(c(min,mid,max),color_pal))
}

column_split_vals <- metadata %>% filter(LRRC15_surfexp=='no') %>% pull(TGFB_inducible)
column_split_vals <- c(yes='TGFB_ind  ',no='  TGFB_nind')[column_split_vals]
names(column_split_vals) <- gsub('_NEG','',metadata %>% filter(LRRC15_surfexp=='no') %>% row.names())
column_split_vals <- column_split_vals[colnames(med_diff_df)]
# png('../scRNAseq_analysis/module_umaps/SCING_hdWGCNA_combined_samplewise_ME_deg_heatmap.png', width=10, height=15, units='in',res=600)
# heatmap_width <- 8
# h <- Heatmap(med_diff_df, name='Cross-Group\nMedian Diff', width = 2*unit(heatmap_width/9, "in"),
#              cell_fun = pval_marker_cell_fun(comb_plot_df[,grepl('p_val_adj',colnames(comb_plot_df))]), col = marker_col_fun(-7,0,7, color_pal),
#              row_split = comb_plot_df$method, cluster_row_slices = T,row_names_gp=gpar(fontsize=15), row_title_gp=gpar(fontsize=20),
#              column_split = c(column_split_vals), cluster_column_slices = T, border = T,gap = unit(2,'mm'), column_gap = unit(2,'mm'), column_title_gp=gpar(fontsize=20),column_names_rot = 45,)
# ht <- draw(h,
#            heatmap_legend_side="left", annotation_legend_side="left")
# dev.off()

### heatmap of scaled difference
# scaled expression
hdwgcna_metacell[['RNA']]$new_scale.data <- t(scale(t(as.matrix(hdwgcna_metacell[['RNA']]$scale.data))))
hdwgcna_value_df <- as.matrix(AverageExpression(hdwgcna_metacell,assays = "RNA", layer = 'new_scale.data',group.by = c('cell_line'), return.seurat = F)$RNA)
scing_metacell[['RNA']]$new_scale.data <- t(scale(t(as.matrix(scing_metacell[['RNA']]$scale.data))))
scing_value_df <- as.matrix(AverageExpression(scing_metacell,assays = "RNA", layer = 'new_scale.data',group.by = c('cell_line'), return.seurat = F)$RNA)
value_df <- rbind(hdwgcna_value_df, scing_value_df)
value_df <- value_df[rownames(comb_plot_df),]
colnames(value_df) <- gsub('HuO9','HUO9',colnames(value_df))
stri_enc_mark(colnames(value_df))
all(stri_enc_isutf8(colnames(value_df)))
colnames(value_df) <- stri_enc_toascii(colnames(value_df))
colnames(value_df) <- gsub('\032','',colnames(value_df), fixed = T)

keep_rows <- c(rownames(comb_plot_df)[rowSums(comb_plot_df[,grepl('p_val_adj', colnames(comb_plot_df))] <=0.05) >=4],'SCING-M0-hdWGCNA-M18')
column_split_vals <- column_split_vals[colnames(value_df)]
pdf(paste0('../scRNAseq_analysis/module_umaps/SCING_hdWGCNA_combined_samplewise_ME_deg_scaled_exp_heatmap',outfile_suffix,'.pdf'), width=8, height=10)
heatmap_width <- 8
h <- Heatmap(value_df[keep_rows,], name='Expression', width = 2*unit(heatmap_width/9, "in"),
             cell_fun = pval_marker_cell_fun(comb_plot_df[keep_rows,grepl('p_val_adj',colnames(comb_plot_df))]), col = marker_col_fun(-2,0,2, color_pal),
             row_split = comb_plot_df[keep_rows,]$method, cluster_row_slices = T,row_names_gp=gpar(fontsize=15), row_title_gp=gpar(fontsize=20),
             column_split = c(column_split_vals), cluster_column_slices = T, border = T,gap = unit(2,'mm'), column_gap = unit(2,'mm'), column_title_gp=gpar(fontsize=15),column_names_rot = 45,)
ht <- draw(h,
           heatmap_legend_side="left", annotation_legend_side="left")
dev.off()


# save data for meta heatmap
saveRDS(list(comb_plot_df=comb_plot_df, value_df=value_df, column_split_vals=column_split_vals), "./scrnaseq_heatmap_hdWGCNA_M.rds")



