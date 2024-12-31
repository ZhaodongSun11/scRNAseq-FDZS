###############Figure1###############################
library(Seurat)
library(ggsci)
celltype_color = c(
  "Mast_cells" = "#b3d1e3","Epithelial_cells" = "#7ec6da",
  "Plasma_cells"="#398caa",
  "Myeloid_cells"="#63b246",
  "Endothelial_cells"="#e13d28","Fibroblasts"="#91656d",
  "Smooth_muscle_cells"="#c0a36f",
  "T_cells_NK_cells"="#e7692f",
  "B_cells"="#cfb3cb"
)

col_epi <- pal_npg()(6)
col_caf <- c(pal_aaas()(8),"#c0a36f")
col_mye <- pal_d3("category20")(12)
col_T <- pal_d3("category20b")(15)
col_B <- pal_nejm()(7)
col_endo <- pal_bmj()(5)
col_pla <- pal_jco()(8)
col_ma <- c("#b3d1e3")
clusterCols <- c(col_epi,
                 col_caf,
                 col_mye,
                 col_T, 
                 col_B, 
                 col_endo, 
                 col_pla,
                 col_ma
)

DimPlot(sce, group.by='celltype1', cols=celltype_color, pt.size=0.01, raster=F,
        label = F)+NoLegend()
DimPlot(sce, group.by='celltype1', cols=celltype_color, pt.size=0.01, raster=F, 
        label = F,split.by = "location",ncol = 2)+NoLegend()

cols <- clusterCols
names(cols) <- c(
  "Epi_c1_CAPN8","Epi_c6_TM4SF1","Epi_c2_TFF1","Epi_c4_TFF3","Epi_c3_PGA3","Epi_c5_UBE2C",   
  "Mesothelial",
  "Fib_c1_CD74", "Fib_c2_TAGLN", "Fib_c3_NEAT1", 
  "Fib_c4_HSPA1B", "Fib_c5_CXCL14", "Fib_c6_CST1",
  "Fib_c7_CLEC3B", 
  "SMC", 
  "Mono_c1_S100A8", "Mono_c2_CD300E", "Mono_c3_IL1B", 
  "Macro_c1_CCL4L2","Macro_c2_FCN1",   "Macro_c3_SPP1", "Macro_c4_APOC1",
  "Macro_c5_LYVE1", "Macro_c6_CXCL10",
  "cDC1_c1_IRF8", "cDC2_c2_CD1C","cDC3_c3_LAMP3", 
  "NKT_c1_FGFBP2", "NK_c1_FCGR3A", "NK_c2_AREG","γδT", "T_STMN1",
  "CD8_c1_GZMK",   "CD8_c2_ZNF683",  "CD8_c3_ISG15", "CD8_c4_TNFSF9",
  "CD8_c5_CEBPD", "CD4_c1_LTB","CD4_c4_FOXP3", "Tn_c1_SELL",  
  "CD4_c2_FOS", "CD4_c3_CXCL13",
  "naiveB_c1_TCL1A", "SwBm_c2_ITGB1", "B_c3_IFIT3", "AtM_c4_DUSP4", 
  "ACB2_c7_NR4A2", "B_c5_HSP", "ACB1_c6_EGR1", 
  "Vein_c1_ACKR1",  "Tip_c2_SPARC",    "Cap_c3_CA4",  "Atr_c4_FBLN5",  "Lym_c5_CCL21", 
  "Pla_c1_CD83,", "Pla_c2_IGHJ5", "Pla_c3_DUSP5", "Pla_c4_IGLC2", 
  "Pla_c5_IGKC", "Pla_c6_IFI6", "Pla_c7_NEAT1", "Pla_c8_PDK1", 
  "Mast"
) 

library(ggplot2)
library(dplyr)
library(tidyverse)

nature_col = clusterCols
df <- data.frame(sce@meta.data[,c('seurat_clusters','celltype','celltype1')], 
                 sce@reductions$umap@cell.embeddings[,1:2])

colnames(df)
cluster_order <- c(
  "Epi_c1_CAPN8","Epi_c6_TM4SF1","Epi_c2_TFF1","Epi_c4_TFF3","Epi_c3_PGA3","Epi_c5_UBE2C",   
  "Mesothelial",
  "Fib_c1_CD74", "Fib_c2_TAGLN", "Fib_c3_NEAT1", 
  "Fib_c4_HSPA1B", "Fib_c5_CXCL14", "Fib_c6_CST1",
  "Fib_c7_CLEC3B", 
  "SMC", 
  "Mono_c1_S100A8", "Mono_c2_CD300E", "Mono_c3_IL1B", 
  "Macro_c1_CCL4L2","Macro_c2_FCN1",   "Macro_c3_SPP1", "Macro_c4_APOC1",
  "Macro_c5_LYVE1", "Macro_c6_CXCL10",
  "cDC1_c1_IRF8", "cDC2_c2_CD1C","cDC3_c3_LAMP3", 
  "NKT_c1_FGFBP2", "NK_c1_FCGR3A", "NK_c2_AREG","γδT", "T_STMN1",
  "CD8_c1_GZMK",   "CD8_c2_ZNF683",  "CD8_c3_ISG15", "CD8_c4_TNFSF9",
  "CD8_c5_CEBPD", "CD4_c1_LTB","CD4_c4_FOXP3", "Tn_c1_SELL",  
  "CD4_c2_FOS", "CD4_c3_CXCL13",
  "naiveB_c1_TCL1A", "SwBm_c2_ITGB1", "B_c3_IFIT3", "AtM_c4_DUSP4", 
  "ACB2_c7_NR4A2", "B_c5_HSP", "ACB1_c6_EGR1", 
  "Vein_c1_ACKR1",  "Tip_c2_SPARC",    "Cap_c3_CA4",  "Atr_c4_FBLN5",  "Lym_c5_CCL21", 
  "Pla_c1_CD83,", "Pla_c2_IGHJ5", "Pla_c3_DUSP5", "Pla_c4_IGLC2", 
  "Pla_c5_IGKC", "Pla_c6_IFI6", "Pla_c7_NEAT1", "Pla_c8_PDK1", 
  "Mast"
)                
scales::show_col(nature_col)
col_cluster <- setNames(clusterCols,
                        c(
                          "Epi_c1_CAPN8","Epi_c6_TM4SF1","Epi_c2_TFF1","Epi_c4_TFF3","Epi_c3_PGA3","Epi_c5_UBE2C",   
                          "Mesothelial",
                          "Fib_c1_CD74", "Fib_c2_TAGLN", "Fib_c3_NEAT1", 
                          "Fib_c4_HSPA1B", "Fib_c5_CXCL14", "Fib_c6_CST1",
                          "Fib_c7_CLEC3B", 
                          "SMC", 
                          "Mono_c1_S100A8", "Mono_c2_CD300E", "Mono_c3_IL1B", 
                          "Macro_c1_CCL4L2","Macro_c2_FCN1",   "Macro_c3_SPP1", "Macro_c4_APOC1",
                          "Macro_c5_LYVE1", "Macro_c6_CXCL10",
                          "cDC1_c1_IRF8", "cDC2_c2_CD1C","cDC3_c3_LAMP3", 
                          "NKT_c1_FGFBP2", "NK_c1_FCGR3A", "NK_c2_AREG","γδT", "T_STMN1",
                          "CD8_c1_GZMK",   "CD8_c2_ZNF683",  "CD8_c3_ISG15", "CD8_c4_TNFSF9",
                          "CD8_c5_CEBPD", "CD4_c1_LTB","CD4_c4_FOXP3", "Tn_c1_SELL",  
                          "CD4_c2_FOS", "CD4_c3_CXCL13",
                          "naiveB_c1_TCL1A", "SwBm_c2_ITGB1", "B_c3_IFIT3", "AtM_c4_DUSP4", 
                          "ACB2_c7_NR4A2", "B_c5_HSP", "ACB1_c6_EGR1", 
                          "Vein_c1_ACKR1",  "Tip_c2_SPARC",    "Cap_c3_CA4",  "Atr_c4_FBLN5",  "Lym_c5_CCL21", 
                          "Pla_c1_CD83,", "Pla_c2_IGHJ5", "Pla_c3_DUSP5", "Pla_c4_IGLC2", 
                          "Pla_c5_IGKC", "Pla_c6_IFI6", "Pla_c7_NEAT1", "Pla_c8_PDK1", 
                          "Mast"
                        )                )

ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(celltype, levels = cluster_order))) + 
  geom_point(size = 0.1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.8,'cm'),
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))

ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(celltype, levels = cluster_order))) + 
  geom_point(size = 0.1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.key=element_rect(fill='white'),
        legend.text = element_text(size=15), 
        legend.key.size=unit(0.6,'cm'),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color = guide_legend(override.aes = list(size=6)))+
  geom_segment(aes(x = min(UMAP_1) , y = min(UMAP_2),xend = min(UMAP_1)+3, yend = min(UMAP_2)),
               colour = "black", size=0.5,
               arrow = arrow(length = unit(0.2,"cm"), 
                             type = "closed"))+
  geom_segment(aes(x = min(UMAP_1), y = min(UMAP_2),xend = min(UMAP_1),yend = min(UMAP_2)+1.5),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.2,"cm"), 
                                                        type = "closed")) +
  annotate("text", x = min(df$UMAP_1) +1.4, y = min(df$UMAP_2) -0.3, label = "UMAP1",
           color="black",size = 5) + 
  annotate("text", x = min(df$UMAP_1) -0.8, y = min(df$UMAP_2) + 0.8, label = "UMAP2",
           color="black",size = 5,angle=90) 


cell <- df %>%group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))

rownames(cell) <- cell$celltype
A <- cell[cluster_order,]
a <- c(1:59)
A$ID <- a
A$ID <- as.factor(A$ID)

library(ggrepel)
p <- ggplot(df,aes(x= UMAP_1 , y = UMAP_2 ,col=factor(celltype, levels = cluster_order))) + 
  geom_point(size = 0.1, shape=16)+
  scale_color_manual("",values = col_cluster)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = 'none',
        axis.title = element_text(colour = 'black', size = 15, hjust = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  geom_text_repel(data=A, aes(label=ID),color="black", size=5, point.padding = 0.3)


ggplot(df, aes(x = UMAP_1, y = UMAP_2, col = factor(celltype, levels = cluster_order))) + 
  geom_point(size = 0.1, shape = 16) +
  scale_color_manual("", values = col_cluster) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'none',
    axis.title = element_text(colour = 'black', size = 15, hjust = 0),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  geom_text_repel(
    data = A, 
    aes(label = ID), 
    color = "black", 
    size = 5, 
    point.padding = 0.3,
    max.overlaps = Inf  
  )

p

B <- A
B$x <- 1
B$lab <- c(59:1)
leg <- B %>%
  mutate(Response_round = round(5 * lab) / 5) %>%
  group_by(x, Response_round) %>% 
  mutate(x = 0.1 * (seq_along(Response_round) - (0.5 * (n() + 1)))) %>%
  ungroup() %>%
  mutate(x = x + as.numeric(as.factor(x))) %>%
  ggplot(aes(x = x, y = lab)) +
  geom_point(shape = 21, size = 8, aes(x = x, y = Response_round, fill=celltype)) +
  geom_text(aes(label = ID, x = x, y = Response_round), size = 6)+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values = col_cluster)+
  annotate("text", x = 1, y = 16, label = "Cluster",
           color="black",size = 6)+
  geom_text(aes(label = celltype, x = x+0.001, 
                y = Response_round), size = 5, hjust=0)+
  scale_x_continuous(expand=c(-0.01,0.01))

library(cowplot)
plotlist <- list(p, leg)
plot_grid(plotlist = plotlist, ncol = 2, align="hv", rel_widths = c(5,2))


library(ggsci)
cl=pal_nejm("default",alpha = 0.5)(10)

markers <- FindAllMarkers(sce, logfc.threshold = 1, min.pct = 0.25, only.pos = T)
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
library(ggplot2)
p <- DotPlot(sce, features = top5$gene,
             cols = cl, group.by = "celltype1")
exp <- p$data
library(forcats)
exp$features.plot <- as.factor(exp$features.plot)
exp$features.plot <- fct_inorder(exp$features.plot)

exp$id <- as.factor(exp$id)
exp$id <- fct_inorder(exp$id)


ggplot(exp,aes(x=features.plot,y= id))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  geom_point(aes(size=`pct.exp`,color=`avg.exp.scaled`),
             shape=21,color="black",stroke=1)+
  theme(panel.background =element_blank(),
        axis.line=element_line(colour="black", size = 1),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=11,color="black",angle=90), 
        axis.text.y=element_text(size=11,color="black"))+
  scale_color_gradientn(colors = colorRampPalette(c("white", "#00C1D4", "#FFED99","#FF7600"))(10))+
  labs(x=NULL,y=NULL)

library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
in.dat <- sce@meta.data
R_oe <- calTissueDist(in.dat,
                      byPatient = F,
                      colname.cluster = "celltype",
                      colname.patient = "patientID",
                      colname.tissue = "location",
                      method = "chisq", 
                      min.rowSum = 0) 
R_oe
R_oe_adjusted <- pmin(R_oe, 2)
col_fun <- colorRamp2(c(0, 1,2), c("#f2f4c9","#349eb7", "#1a2852"))
Heatmap(as.matrix(R_oe_adjusted),
        show_heatmap_legend = TRUE, 
        cluster_rows = TRUE, 
        cluster_columns = TRUE,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = col_fun, 
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0, 2, by = 0.5),
          labels =  seq(0, 2, by = 0.5),
          legend_gp = gpar(fill = col_fun(seq(0, 2, by = 0.5)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", R_oe_adjusted[i, j]), x, y, gp = gpar(fontsize = 8))
        }
)

library(Seurat)
library(dplyr)
library(reshape2)
library(plyr)
library(ggplot2)
prop.table(table(Idents(sce)))
table(sce$celltype, sce$orig.ident)
Cellratio <- prop.table(table(sce$celltype1, sce$orig.ident), margin = 2)
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
sub <- factor(Cellratio$Var2)
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =sub, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  scale_fill_manual(values = pal)+
  labs(x='Tissue',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

###################Figure2################################
library(Seurat)
library(Nebulosa)
library(ggnetwork)
library(dplyr)

plot_density(sce, features = gene,
             pal = 'magma', raster = T, size = 0.8,reduction = "umap") &
  theme_blank()&
  theme(legend.frame = element_rect(colour = "black"),
        legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.8, "cm"),
        legend.title = element_text(color = 'black', face = "bold", size=8))

library(CytoTRACE2) 
results <- cytotrace2(sce,
                      species = "human",
                      is_seurat = TRUE,
                      full_model = FALSE,
                      slot_type = "counts"
)

annotation1 <- data.frame(phenotype = sce@meta.data$location) %>% 
  set_rownames(., colnames(sce))

plots <- plotData(cytotrace2_result = results, 
                  annotation = annotation1, 
                  is_seurat = TRUE)

emb_1 <- all@reductions$umap@cell.embeddings[,1]
emb_2 <- all@reductions$umap@cell.embeddings[,2] 

plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2

plots[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_Potency_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2

plots[["CytoTRACE2_Relative_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["CytoTRACE2_Relative_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2

plots[["Phenotype_UMAP"]][[1]][["data"]]["UMAP_1"] <- emb_1
plots[["Phenotype_UMAP"]][[1]][["data"]]["UMAP_2"]<- emb_2

plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$Phenotype_UMAP
plots$CytoTRACE2_Boxplot_byPheno

library(monocle3)
library(Seurat)
library(ggplot2)
library(dplyr)

expression_matrix <- as(as.matrix(data@assays$RNA@counts), 'sparseMatrix')
cell_metadata <- data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,reduction_method='UMAP',
                        preprocess_method = 'PCA')
plot_cells(cds, color_cells_by="celltype",cell_size=0.5,group_label_size=5)
plot_cells(cds, color_cells_by="orig.ident",cell_size=0.5,group_label_size=5)
cds <- cluster_cells(cds)
cds <- cluster_cells(cds, resolution =0.0025) 
plot_cells(cds,cell_size=0.5,group_label_size=5)
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed 
plot_cells(cds, color_cells_by="seurat_clusters",
           cell_size=0.5,group_label_size=4) 
mycds <- cds
mycds <- learn_graph(mycds,
                     verbose=T,use_partition = T,
                     learn_graph_control=list(minimal_branch_len=30,
                                              euclidean_distance_ratio=20
                     ))
plot_cells(mycds, 
           color_cells_by = 'seurat_clusters',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)

plot_cells(mycds, 
           color_cells_by = "celltype", 
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE, 
           label_branch_points=TRUE,
           graph_label_size=4)
mycds1 <- mycds
mycds1 <- order_cells(mycds1)
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 4, 
           cell_size=0.5, 
           trajectory_graph_segment_size = 2)

pd <- pseudotime(mycds1, reduction_method = 'UMAP')
data <- AddMetaData(data,metadata = pd,col.name = 'pseudotime')

mycds2_res <- graph_test(mycds1, 
                         neighbor_graph="principal_graph", cores=4)

genes_sig <- mycds2_res %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()
genes_sig <- c("RACK1","SELENOF","SELENOM","SELENOH","SELENOP",
               "SELENOW","SELENOK","SELENOS","SELENOT","GPX4"
)
plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="seurat_clusters", 
                         min_expr=0.5, ncol = 2)


DotPlot(sce,features = trans,
        group.by = "celltype1",
)+coord_flip()+theme_bw()+  
  theme(panel.grid = element_blank(),  
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+  
  scale_color_gradientn(values = seq(0,1,0.2),colors = c('#330066','#336699','#66CC66','#FFCC33'))+  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 

genes <- c("FLT1","KDR","NRP1","FLT4","NRP2","TIE1","TEK","PDGFB",#Tumor vascularization
           "JAG1","EFNB2","HEY1","EPHB4","HES1",#Morphogenesis
           "NT5E","STAB1","CCL2",#Inflammation
           "HLA-A","HLA-B","HLA-C",#MHC-I
           "HLA-DMA","HLA-DMB","HLA-DOA","HLA-DPA1",
           "HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA",
           "HLA-DRB1","HLA-DRB5",#MHC-II
           "ROBO1","ROBO4","PLXND1","ROBO3","UNC5B"#Guidance receptors
)

anergy <- c("ICAM1","VCAM1","SELE","SELP",#Endothelial Adhesion Molecules
            "IDO1","HAVCR2","CD274","PDCD1LG2",
            "FASLG","LGALS1"#Immunosuppressive Molecules
)

library(Seurat)
library(ggplot2)
library(tidyverse)
markers <- genes

markers <- as.data.frame(markers)
markers$x <- 1 
markers$y <- seq(1:nrow(markers))

markers$group <- c(rep('Tumor vascularization',8),rep('Morphogenesis',5),rep('Inflammation',3),rep('MHC-I',3),
                   rep("MHC-II",11),rep("Guidance receptors",5))

markers$group <- c(rep('Endothelial Adhesion Molecules',5),rep('Immunosuppressive Molecules',6))

marker_ave_exp <- AverageExpression(sce,assays = "RNA",
                                    features = markers$markers,
                                    group.by = "celltype",
                                    slot="data")
markers$markers=factor(markers$markers,levels = markers$markers)

data=base::apply(marker_ave_exp$RNA,1,
                 function(x) (x-mean(x))/sd(x))%>%t()%>%as.data.frame()%>%rownames_to_column('Gene')%>%reshape2::melt()

data$Gene=factor(data$Gene,levels = rev(unique(data$Gene)))

library(RColorBrewer)
color_palette <- colorRampPalette(brewer.pal(11, "RdYlBu"))(100)  

p1=ggplot(data,aes(x = variable,y = Gene,fill=value))+
  geom_tile()+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradientn(colors=rev(colorRampPalette(color_palette)(500)),
                       limits=c(-2,2),name="Z Score",
                       oob = scales::squish
  )+
  geom_vline(xintercept=as.numeric(cumsum(table(unique(data$variable)))+0.5),linetype=2)+
  geom_hline(yintercept=as.numeric(cumsum(c(5.5,11,3,3,5,8))),linetype=2)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text = element_text(face="bold"),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 10))
p1
palette1 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99")
p2=ggplot(markers,aes(x,reorder(y,-y),fill=group))+
  geom_tile()+
  geom_text(aes(label=markers),size=3)+
  scale_fill_manual(values = palette1)+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')+
  scale_x_continuous(expand = c(0,0))
color_labels <- c("Tumor vascularization", "Morphogenesis",'Inflammation', "MHC-I", "MHC-II", "Guidance receptors")
label_data <- data.frame(y = color_labels)
label_data$y=factor(label_data$y,levels = c("Tumor vascularization", "Morphogenesis",'Inflammation', "MHC-I", "MHC-II", "Guidance receptors"))
label_data$y=factor(label_data$y,levels = c("Endothelial Adhesion Molecules", "Immunosuppressive Molecules"))
p2_color_annotations <- ggplot(label_data, aes(x = 0, y = y)) +
  geom_tile(aes(fill = y), color = "white") +
  scale_fill_manual(values = palette1) +
  geom_text(aes(label = y), hjust = 0, size = 4, color = "black") +
  theme_void()
p2_color_annotations
heatmap_with_color_annotations = p2_color_annotations +p2 +p1 +
  plot_layout(ncol = 2, widths  = c(1, 3))
heatmap_with_color_annotations & theme(plot.margin = margin(0,0,0,0))

library(patchwork)
heatmap1 =p2+p1+plot_layout(ncol = 2, widths  = c(1, 3))
heatmap1 & theme(plot.margin = margin(0,0,0,0))

###################################Figure3######################################
library(scRepertoire)
library(immunarch)
library(Seurat)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggplot2)

for(i in 1:length(f)){
  tmp = read.csv(paste0("TCR/",f[[i]],"/",s[[i]],"/filtered_contig_annotations.csv"))
  t[[i]] = paste0(f[[i]],"_",s[[i]])
  BCR_list[[i]] <- tmp
}

names(BCR_list)  = t

f = dir("~/LN/TCR")
s = str_split(f, "-",simplify = T)[,1]
t = str_split(f, "-",simplify = T)[,2]
a = c('Pt1')
for(i in 1:length(f)){
  tmp = read.csv(paste0(f[[i]]))
  a[[i]] = paste0(s[[i]],"-",t[[i]])
  BCR_list1[[i]] <- tmp
}
names(BCR_list1)  = a

data_bcr1 <- combineTCR(BCR_list,
                        ID = t,
                        samples=s
)

data_bcr <- addVariable(data_bcr, name = "location", 
                        variables = c(c("Ascites", "M", "M", "Ascites", "M", "M", "Ascites", "M", "M", 
                                        "M", "PT", "Ascites", "M", "M", "PT", "Ascites", "M", "PT", 
                                        "Ascites", "PT", "Ascites", "M", "PBMC", "M", "PT", "Ascites", "N","M", "PBMC", "M", "PT" )))
library(ggsci)
library(RColorBrewer)
display.brewer.all()  
getPalette = colorRampPalette(brewer.pal(9,"Set1"))
p1 <- quantContig(data_bcr, 
                  cloneCall="strict", 
                  scale = T,
                  chain = "both")+
  scale_fill_manual(values = getPalette(39))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

p2 <- quantContig(data_bcr, 
                  cloneCall="strict", 
                  scale = T,
                  chain = "TRA")+
  scale_fill_manual(values = getPalette(39))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

quantContig(data_bcr, 
            cloneCall="strict", 
            scale = F,
            chain = "both",
            exportTable=T)

data <- quantContig(data_bcr, 
                    cloneCall="strict", 
                    scale = T,
                    chain = "both",
                    exportTable=T)
data$group <- str_split(data$values,"_",simplify = T)[,1]
library(ggplot2)
library(ggpubr)
ggplot(data, aes(fill=group, y=scaled, x=group))+
  geom_bar(position=position_dodge(),
           stat="summary",
           width=0.7,
           colour = "black",
           size=1)+
  stat_summary(fun.data = 'mean_se', 
               geom = "errorbar", 
               colour = "black",
               width = 0.2,
               position=position_dodge(0.7))+
  theme_classic()

quantContig(data_bcr,
            cloneCall="strict",
            scale = T,
            chain = "both",
            group.by = 'celltype1')

p3 <- abundanceContig(data_bcr, cloneCall = "strict", scale = F)

p4 <- abundanceContig(data_bcr, cloneCall = "strict", scale = F, group.by = 'sample')

p5 <- compareClonotypes(data_bcr, 
                        numbers = 10, 
                        cloneCall="gene+nt", 
                        samples = c("PT","Ascites","M"),
                        graph = "alluvial") + scale_colour_brewer(palette = "Greens")

p6 <- compareClonotypes(data_bcr, 
                        cloneCall="gene+nt", 
                        numbers = 5,
                        graph = "alluvial",
                        samples = c("OM","PM"))

p7 <- clonalHomeostasis(data_bcr, cloneCall = "aa", 
                        cloneTypes = c(Rare = 1e-04, 
                                       Small = 0.001, 
                                       Medium = 0.01, 
                                       Large = 0.1, 
                                       Hyperexpanded = 1))

p8 <- clonalHomeostasis(data_bcr, cloneCall = "aa", 
                        cloneTypes = c(Rare = 1e-04, 
                                       Small = 0.001, 
                                       Medium = 0.01, 
                                       Large = 0.1, 
                                       Hyperexpanded = 1),
                        group.by = 'sample')

p9 <- clonalProportion(data_bcr, cloneCall = "gene",
                       split = c(10, 100, 1000, 10000, 30000, 1e+05)) 

p10 <- clonalDiversity(data_bcr, 
                       cloneCall = "aa", 
                       n.boots = 1000, 
                       x.axis='sample', 
                       group='ID')


for (i in seq_along(data_bcr)) {
  data_bcr[[i]] <- stripBarcode(data_bcr[[i]], 
                                column = 1, connector = "_", num_connects = 4)
}

SeuratObj_combine <- sce
scBCR_RNA <- combineExpression(data_bcr, 
                               SeuratObj_combine, 
                               cloneCall="aa", 
                               cloneTypes=c(Single=1, Small=3, Medium=10, Large=30, Hyperexpanded=100),
                               proportion = FALSE)

scBCR_RNA@meta.data$TCR.known <- "yes"
scBCR_RNA@meta.data$TCR.known[is.na(scBCR_RNA@meta.data$cloneType)] <-  "no"
scBCR_RNA@meta.data$both_chains <- "yes"
scBCR_RNA@meta.data$both_chains[grep("NA",scBCR_RNA@meta.data$CTgene)] <- "no"
scBCR_RNA@meta.data$both_chains[is.na(scBCR_RNA@meta.data$CTaa)]<- "no"
scBCR_RNA@meta.data$TYPEs = scBCR_RNA@meta.data$CTgene
scBCR_RNA@meta.data$TYPEs[is.na(scBCR_RNA@meta.data$TYPEs)] <-  "NA_NA"
heavychains <- unlist(strsplit(scBCR_RNA@meta.data$TYPEs, "[_]"))[seq(1, length(unlist(strsplit(scBCR_RNA@meta.data$TYPEs, "[_]"))), 2)]
heavychains <- as.data.frame(heavychains)
heavychains$Isotype <- "unknown"
heavychains$Isotype[grep("TRAV",heavychains$heavychains)] <- "TRAV"
scBCR_RNA@meta.data$Isotype <- heavychains$Isotype
remove(heavychains)

DimPlot(scBCR_RNA, group.by = "cloneType",label = T)
Idents(scBCR_RNA) <- "celltype1"
p11 <- clonalOverlay(scBCR_RNA, reduction = "umap", freq.cutpoint = 10, bins = 10, facet = "location") +
  guides(color = "none")+
  scale_color_manual("",values = cols)      

library(ggsci)
scBCR_RNA$cloneType <- factor(scBCR_RNA$cloneType, levels = c("Hyperexpanded (30 < X <= 100)",
                                                              "Large (10 < X <= 30)",
                                                              "Medium (3 < X <= 10)",
                                                              "Small (1 < X <= 3)",
                                                              "Single (0 < X <= 1)"))
Idents(scBCR_RNA) <- "celltype1"
pal = c("#8B0000", "#FF4500", "#FFA07A", "#FFD700", "#FFFACD")
 occupiedscRepertoire(scBCR_RNA, x.axis = "ident",label = F,proportion=T)+
  scale_fill_manual(values = pal)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

StartracDiversity(scBCR_RNA, 
                         type = "celltype1", 
                         sample = "location", 
                         by = "overall")
scBCR_RNA <- scBCR_RNA[,scBCR_RNA$celltype1 %in% c("NKT_c1_FGFBP2", "NK_c1_FCGR3A", "NK_c2_AREG","γδT", "T_STMN1",
                                                   "CD8_c1_GZMK",   "CD8_c2_ZNF683",  "CD8_c3_ISG15", "CD8_c4_TNFSF9",
                                                   "CD8_c5_CEBPD", "CD4_c1_LTB","CD4_c4_FOXP3", "Tn_c1_SELL",  
                                                   "CD4_c2_FOS", "CD4_c3_CXCL13")]
p17 <- p17 + 
  scale_fill_manual(values = cols) +
  theme_minimal() 


Singlecellratio <- function (seu, by = "cell.type",meta.include = NULL, 
                                      group_by = NULL, shape_by = NULL,
                                      custom_fill_colors = NULL, group_by.point = NULL, color_by = NULL, 
                                      pb = FALSE, comparisons = my_comparisons, 
                                      ncol = NULL, label = c("p.format","p.signif"), 
                                      label.x = NA, pt.size = NA) 
{
  
  by <- match.arg(by)  
  if (is.null(group_by)){ 
    group_by <- "null.group" 
  } 
  shapes <- NULL 
  if (!is.null(shape_by)) {
    shapes <- c(16,15,3,7,8,18,5,6,2,4,1,17)
  }
  
  fq <- prop.table(table(seu@meta.data$celltype, seu@meta.data[,"orig.ident"]), margin=2) *100
  df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", "orig.ident")) 
  
  uniques <- apply(seu@meta.data, 2, function(x) length(unique(x))) 
  ei <- unique(seu@meta.data[, names(uniques[uniques<=100])])
  ei <- unique(ei[,colnames(ei) %in% meta.include])
  df <- merge(df, ei, by = "orig.ident")
  df <- cbind(df, null.group = paste("1"))
  df$orig.ident <- as.factor(df$orig.ident)
  
  
  if (is.null(x = ncol)) {
    ncol <- 3
    
    if (length(unique(df$celltype)) > 9) {
      ncol <- 4
    }
    if (length(unique(df$celltype)) > 20) {
      ncol <- 5
    }
  }
  
  custom_fill_colors = c(ggsci::pal_d3()(4),
                         RColorBrewer::brewer.pal(9, "Oranges")[2], 
                         RColorBrewer::brewer.pal(9, "Reds")[6], 
                         RColorBrewer::brewer.pal(9, "Oranges")[5], 
                         RColorBrewer::brewer.pal(9, "Blues")[4:9])
  
  
  
  p <- ggplot(df, aes_string(y = "freq", x = group_by)) + labs(x = NULL,y = "Proportion (%)") + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                                                                                                   panel.grid.major = element_blank(), strip.background = element_rect(fill = NA, 
                                                                                                                                                                                       color = NA), 
                                                                                                                   strip.text = element_text(face = "bold", size = 14), 
                                                                                                                   axis.ticks.x = element_blank(), axis.text = element_text(color = "black"), 
                                                                                                                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,color = 'black', size = 12),
                                                                                                                   axis.text.y = element_text(color = 'black', hjust = 1, vjust = 0.5, size = 12),
                                                                                                                   axis.title.y = element_text(color = 'black', size = 14))
  
  
  if(by=="cell.type" && color_by=="cell.type") {
    p + facet_wrap(group_by, scales = "free_x") + 
      geom_bar(aes_string(x = "orig.ident", fill = "factor(cell.type)"), position = "fill", stat = "identity") + 
      scale_fill_manual("cell.type", values = c("#FB8072", "#1965B0", "#7BAFDE", "#882E72","#B17BA6", 
                                                "#FF7F00", "#FDB462", "#E7298A", "#E78AC3","#33A02C", 
                                                "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D","#E6AB02")) + 
      scale_y_continuous(expand = c(0, 0), labels = seq(0, 100, 25)) + theme(panel.border = element_blank())
  }
  else {
    switch(by, cell.type = p + facet_wrap("cell.type", scales = "free_y", 
                                          ncol = ncol) + guides(fill = FALSE) + geom_boxplot(aes_string(x = group_by), 
                                                                                             alpha = 0.25, outlier.color = NA) + geom_point(size = 4, position = position_jitter(width = 0.25), 
                                                                                                                                            aes_string(x = group_by, y = "freq", color = color_by, 
                                                                                                                                                       shape = shape_by)) + scale_shape_manual(values = shapes) + 
             theme(panel.grid.major = element_line(color = "grey", 
                                                   size = 0.25)) + scale_color_manual(values = custom_fill_colors) + scale_fill_manual(values = custom_fill_colors)) + 
      scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.2)))+
      ggpubr::stat_compare_means(mapping = aes_string(group_by), comparisons = comparisons, label = label,method = "wilcox.test")
  }
  
}
Singlecellratio(sce, group_by = "location",
                meta.include = c("location","orig.ident"),
                comparisons = my_comparisons, color_by = 'location',
                group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                label = 'p.format', ncol =3)


library(CellChat)
library(Seurat)

HD_input <- sce[["RNA"]]@data
HD_meta <- sce@meta.data[,c("location","celltype","celltype1")]
colnames(HD_meta) <-  c("location","celltype","celltype1")
identical(colnames(HD_input),rownames(HD_meta)) 
HD.cellchat <- createCellChat(object = HD_input, meta = HD_meta, group.by = "celltype")
levels(HD.cellchat@idents) 
groupSize <- as.numeric(table(HD.cellchat@idents)) 

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

HD.cellchat@DB <- CellChatDB
HD.cellchat <- subsetData(HD.cellchat) 
future::plan("multisession", workers = 2) 
HD.cellchat <- identifyOverExpressedGenes(HD.cellchat)
HD.cellchat <- identifyOverExpressedInteractions(HD.cellchat)
HD.cellchat@idents <- droplevels(HD.cellchat@idents)
HD.cellchat <- computeCommunProb(HD.cellchat, type = "triMean")
HD.cellchat <- filterCommunication(HD.cellchat, min.cells = 10)
HD.cellchat <- computeCommunProbPathway(HD.cellchat)
HD.cellchat <- aggregateNet(HD.cellchat)
df.net <- subsetCommunication(HD.cellchat)

groupSize <- as.numeric(table(HD.cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(HD.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(HD.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pheatmap::pheatmap(HD.cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")

pathways.show <- c("CCL") 
netVisual_aggregate(HD.cellchat, signaling = pathways.show,targets.use = "LYVE1_Macro")

netVisual_aggregate(HD.cellchat, signaling = pathways.show, layout = "chord")

p1 = netVisual_heatmap(HD.cellchat, signaling = pathways.show,color.heatmap = "Reds")
p2 = netVisual_heatmap(HD.cellchat,color.heatmap = "Reds")
p1+p2


##########################################Figure4###################################
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(ggbeeswarm)
library(scater)
library(scales)
library(forcats)
library(data.table)
library(stringr)
library(dplyr)


sce <- as.SingleCellExperiment(sce)
sce_milo <- miloR::Milo(sce)
sce_milo <- miloR::buildGraph(sce_milo, k = 30, d = 50)
sce_milo <- makeNhoods(sce_milo, 
                           prop = 0.2, 
                           k = 30, 
                           d=50, 
                           refined = TRUE)

plotNhoodSizeHist(sce_milo) 
sce_milo <- countCells(sce_milo, meta.data = data.frame(colData(sce_milo)), 
                           sample="orig.ident")

traj_design <- data.frame(colData(sce_milo))[,c("orig.ident", "group")]
traj_design$orig.ident <- as.factor(traj_design$orig.ident)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident

sce_milo <- calcNhoodDistance(sce_milo, d=50)
da_results <- testNhoods(sce_milo, 
                         design = ~ group, 
                         design.df = traj_design)
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
sce_milo <- buildNhoodGraph(sce_milo)

plotUMAP(PMN_sce_milo, colour_by = "group")

pdf('milo.pdf', width=8, height=8)
plotNhoodGraphDA(sce_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) 


library(pheatmap)
library(RColorBrewer)
library(scales)

all.score.df=data.frame()
all.score.topn.df=data.frame()

for (i in dirs) {
  score.file=paste(dir_output,"/",i,"_program.Zscore.txt",sep = "")
  score.df=read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
  colnames(score.df) <- gsub("^X", "", colnames(score.df))
  if (i==dirs[1]) {all.score.df=score.df}
  if (i!=dirs[1]) {
    all.score.df=all.score.df%>%inner_join(score.df,by="gene")
  }
  
  score.topn.file=paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = "")
  score.topn.df=read.table(score.topn.file,header = T,sep = "\t",stringsAsFactors = F)
  if (i==dirs[1]) {all.score.topn.df=score.topn.df}
  if (i!=dirs[1]) {
    all.score.topn.df=cbind(all.score.topn.df,score.topn.df)
  }
}

rownames(all.score.df)=all.score.df$gene
all.score.df$gene=NULL
all.score.df=all.score.df[rowSums(is.na(all.score.df)) == 0,] 
all.score.rm.df=all.score.df[,setdiff(colnames(all.score.df),maybe.bg)] 
all.score.rm.df.cor=cor(all.score.rm.df,method = "pearson")

all.score.rm.df.cor[all.score.rm.df.cor < cor_min]=cor_min
all.score.rm.df.cor[all.score.rm.df.cor > cor_max]=cor_max

colanno=as.data.frame(colnames(all.score.rm.df.cor))
colnames(colanno)="colnames"
colanno$sample <- str_replace(colanno$colnames, "_\\d+$", "")
colanno$sample <- str_replace_all(colanno$sample, "\\.", "-")
rownames(colanno)=colanno$colnames
colanno$colnames=NULL

rowanno=as.data.frame(rownames(all.score.rm.df.cor))
colnames(rowanno)="rownames"
rowanno$sample <- str_replace(rowanno$rownames, "_\\d+$", "")
rowanno$sample <- str_replace_all(rowanno$sample, "\\.", "-")
rownames(rowanno)=rowanno$rownames
rowanno$rownames=NULL

if (is.null(color)){
  color_v=colorRampPalette(RColorBrewer::brewer.pal(39, "Set1"))(length(dirs))
  names(color_v)=cleaned_dirs
}else{
  color_v=color
}
ann_colors = list(sample = color_v)

tmpp=pheatmap(all.score.rm.df.cor,cluster_rows = T,cluster_cols = T,
              clustering_method = cluster_method, 
              show_colnames = F,
              treeheight_row=30,treeheight_col=0,
              border_color=NA,
              annotation_row = rowanno,annotation_col = colanno,
              annotation_names_row = F,annotation_names_col = F,
              annotation_colors = ann_colors,
              color = colorRampPalette(c("white","yellow", "red","#67001F"))(50),
              fontsize_row=12,
              width = 30,height = 25,
              filename = paste(dir_output,"/","program_pearson_cor.",cluster_method,".heatmap.pdf",sep = "")
)
colanno$sample <- factor(colanno$sample)
print(levels(colanno$sample))
print(names(ann_colors$sample))

library(monocle)
library(Seurat)
library(tidyverse)
library(RColorBrewer)
data <- GetAssayData(sce, assay = "RNA", slot = 'counts')
data <- data[rowSums(as.matrix(data)) != 0,]
pd <- new("AnnotatedDataFrame", data = sce@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)
monocds <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily=negbinomial.size())
monocds <- estimateSizeFactors(monocds)
monocds <- estimateDispersions(monocds)
monocds <- detectGenes(monocds, min_expr = 0.1)
test <- monocds
marker_gene = markers %>% group_by(cluster) %>% top_n(n =5, wt = avg_log2FC)
test_ordering_genes=unique(marker_gene$gene)
test=setOrderingFilter(test,ordering_genes = test_ordering_genes)
test=reduceDimension(test,reduction_method = "DDRTree",max_components = 2, norm_method = 'log',residualModelFormulaStr = "~num_genes_expressed")
test=orderCells(test)
meta <- sce@meta.data
phenoData(test) <- new("AnnotatedDataFrame", data = meta)

plot_cell_trajectory(test,color_by = "location") +
  scale_color_manual(breaks = c("PT", "Peri_Met", "Ascites"), values= c("#d62728", "#ff7f0e", "#1f77b4")) + theme(legend.position = "right")
ggsave("tissue.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))

plot_cell_trajectory(test,color_by = "orig.ident")
ggsave("State.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
plot_cell_trajectory(test,color_by = "Pseudotime")
ggsave("Pseudotime.pdf",device = "pdf",width = 8,height = 9,units = c("cm"))
test=orderCells(test,root_state = 6) 
plot_cell_trajectory(test,color_by = "Pseudotime")
expressed_genes=row.names(subset(fData(test),num_cells_expressed>=10))
pseudotime_de <- differentialGeneTest(test[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
pseudotime_de <- pseudotime_de[order(pseudotime_de$qval), ]


library(Seurat)
library(RColorBrewer)
library(GSVA)
library(GSEABase)
library(msigdbr)

msigdbr_species()
human <- msigdbr(species = "Homo sapiens")
human[1:5,1:5]
table(human$gs_cat)
table(human$gs_subcat)
genesets = msigdbr(species = "Homo sapiens",
                   category = "H") %>% 
  dplyr::select(gs_name,gene_symbol)
genesets = genesets %>% split(x = .$gene_symbol, f = .$gs_name)
names(genesets) <- sub("^HALLMARK_", "", names(genesets))

gene_matrix <- as.data.frame(human_data@assays$RNA@data)
GSVA <- gsva(as.matrix(gene_matrix), genesets, min.sz=10, max.sz=1000,verbose=TRUE)
head(rownames(GSVA))
rownames(GSVA) <- tolower(rownames(GSVA))
rownames(GSVA) <- gsub('_',' ',rownames(GSVA))
library(Hmisc)
rownames(GSVA) <- capitalize(rownames(GSVA))

meta=human_data@meta.data[,c("orig.ident","celltype1","group")]
row.names(meta)=colnames(GSVA)
GSVA_Seurat <- CreateSeuratObject(counts = GSVA, meta.data = meta, project = "GSVA_singleCell")
Idents(GSVA_Seurat) <- "celltype1"
GSVA_celltype=FindAllMarkers(GSVA_Seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
GSVA_exp <- AverageExpression(GSVA_Seurat,features = rownames(GSVA),group.by = 'celltype1', assays = 'RNA', slot = "counts") 
GSVA_exp <- as.data.frame(GSVA_exp$RNA)

pathways <- rownames(GSVA)
gsva_plot <- GSVA_exp[pathways,]

library(ComplexHeatmap)
pheatmap(gsva_plot,cluster_rows = T,cluster_cols = F,scale = 'none',
         colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
         border=FALSE,cellwidth = 9, cellheight = 9,heatmap_legend_param = list(title="GSVA score"))

##################################Figure5######################################
library(Seurat)
library(dplyr)

adj_flux <- read.csv("~/scFEA/scFEA/output/adj_flux.csv",row.names = 1)
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))

adj_bala <- read.csv("~/scFEA/scFEA/output/adj_balance.csv",row.names = 1)
adj_scRNA$group_cells <- paste0(adj_scRNA$orig.ident)
cell_anno <- data.frame(cellid=rownames(adj_scRNA@meta.data),
                        group_cells=adj_scRNA$group_cells)
cell_anno <- cell_anno[order(cell_anno$group_cells),]
adj_flux <- adj_flux[cell_anno$cellid,]
df_averages <- adj_flux %>%
  group_by(group = cell_anno$group_cells) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::select(-group)

rownames(df_averages) <- unique(cell_anno$group_cells)
df_averages <- t(df_averages)%>% as.data.frame()


df_flux <- df_averages[,c("Epi_c1_CAPN8","Epi_c6_TM4SF1","Epi_c2_TFF1",
                          "Epi_c4_TFF3","Epi_c3_PGA3","Epi_c5_UBE2C" )]


colnames(df_flux) <- c("Epi_c1_CAPN8","Epi_c6_TM4SF1","Epi_c2_TFF1",
                       "Epi_c4_TFF3","Epi_c3_PGA3","Epi_c5_UBE2C" )

df_flux = df_flux[apply(df_flux, 1, function(x) sd(x)!=0),]

library(ComplexHeatmap)
df_flux[is.na(df_flux)] <- 0
pheatmap(df_flux, cluster_cols = F, cluster_rows = T,
         show_rownames = F, scale = "row",
         colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
         border=F,heatmap_legend_param = list(title="Flux"),
         fontsize_col=8,
         treeheight_row=20)

annotation_col = data.frame(
  celltype = c(rep("Epi",6)))
row.names(annotation_col) <- colnames(df_flux)
human_moduleInfo <- read.csv("scFEA.human.moduleinfo.csv", header = T, row.names = 1)

annotation_row = human_moduleInfo[rownames(df_flux),]
annotation_row  = as.data.frame(annotation_row[,c("SM_anno")])
rownames(annotation_row) = rownames(df_flux)
colnames(annotation_row) = c("SM_anno")
cellcolor <- c('#66C5CC')
names(cellcolor) <- c("Epi")
modulecolor <- c("#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00",
                 "#808000","#FF00FF","#FA8072","#800080","#87CEEB","#40E0D0","#5F9EA0",
                 "#008B8B","#FFE4B5","#228B22","#4682B4","#32CD32","#F0E68C","#FFFFE0",
                 "#FF6347")
names(modulecolor) <- unique(annotation_row$SM_anno)
ann_colors <- list(celltype1=cellcolor,  SM_anno=modulecolor) 

pheatmap(df_flux, scale = "row",show_rownames = T,
         show_colnames = T,
         cluster_cols = F, cluster_rows = T,
         col = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_names_row = F,
         annotation_names_col = F ,
         column_title = NULL,
         row_title = NULL,
         fontsize_col=8,
         treeheight_row=20,
         heatmap_legend_param = list(title="Flux"),
         annotation_colors = ann_colors)

colnames(human_moduleInfo)
human_moduleInfo$module_name <- paste0(human_moduleInfo$Module_id,": ",
                                       human_moduleInfo$Compound_IN_name,
                                       "_", human_moduleInfo$Compound_OUT_name)

select_moduleInfo = human_moduleInfo[rownames(df_flux),]

df_flux_new <- df_flux
rownames(df_flux_new) <- select_moduleInfo$module_name
row.names(annotation_col) <- colnames(df_flux_new)
rownames(annotation_row) = rownames(df_flux_new)
modules <- c("M_61: Putrescine_GABA",
             "M_62: Ornithine_Putrescine",
             "M_63: Arginine_Ornithine",
             "M_64: Arginine_Putrescine",
             "M_65: Argininosuccinate_Arginine + Fumarate",
             "M_66: Citruline + Aspartate_Argininosuccinate",
             "M_67: Ornithine_Proline",
             "M_150: PRPP_UMP",
             "M_151: Orotidylic acid_UMP",
             "M_152: UMP_B-Alanine",
             "M_153: UMP_CDP",
             "M_154: Uracil_B-Alanine",
             "M_155: UTP_CDP",
             "M_156: CDP_Cytidine",
             "M_157: CDP_dCDP",
             "M_158: dCDP_Deoxycytidine",
             "M_159: dCMP_Deoxycytidine",
             "M_160: dCDP_dUMP",
             "M_161: dCDP_dCTP",
             "M_162: dUMP_Deoxyuridine",
             "M_163: dUMP_dTMP",
             "M_164: dTMP_Succinyl-CoA",
             "M_165: dTMP_dTTP",
             "M_166: Thymine_Succinyl-CoA",
             "M_171: dCDP_dCMP"
             
)
modules <- as.data.frame(modules)
new_at <- which(rownames(df_flux_new) %in% modules$modules)
pdf("./3-anno-celltypeFlux.pdf", width = 9, height = 7)
pheatmap(df_flux_new, scale = "row",show_rownames = T,
         show_colnames = T,
         cluster_cols = F, cluster_rows = T,
         col = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100),
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_names_row = F,
         annotation_names_col = F ,
         column_title = NULL,
         row_title = NULL,
         fontsize_col=8,
         treeheight_row=20,
         heatmap_legend_param = list(title="Flux"),
         annotation_colors = ann_colors)+
  rowAnnotation(link = anno_mark(at = new_at, 
                                 labels = modules$modules, labels_gp = gpar(fontsize = 9)))


dev.off()
adj_flux <- read.csv("~/scFEA/scFEA/output/adj_flux.csv",row.names = 1)
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))
predFlux <- t(data.matrix(adj_flux))
adj_scRNA[["FLUX"]] <- CreateAssayObject(counts = predFlux)
DefaultAssay(adj_scRNA) <- 'FLUX'

module_deg = FindMarkers(adj_scRNA, ident.1="nonPM", 
                         ident.2="PM",group.by = "group",
                         test.use = "LR",
                         min.cells.group = 1,
                         min.cells.feature = 1,
                         min.pct = 0,
                         logfc.threshold = 0)
DEG_FLUX <- list()

Idents(adj_scRNA) = 'celltype11'
cells <- c('Epi_c6_TM4SF1', 'Others')

for (i in 1:length(cells)){
  
  obj = subset(adj_scRNA, idents = cells[i])
  module_deg = FindMarkers(adj_scRNA, ident.1="Epi_c6_TM4SF1", 
                           ident.2="Others",
                           test.use = "LR",
                           min.cells.group = 1,
                           min.cells.feature = 1,
                           min.pct = 0,
                           logfc.threshold = 0)
  rownames(module_deg) <-  gsub("-", "_", rownames(module_deg))
  module_deg$MN <- rownames(module_deg)
  
  module_deg <- module_deg[human_moduleInfo$Module_id, ]
  module_deg <- cbind(module_deg, human_moduleInfo)
  
  DEG_FLUX[[1]] <- module_deg
  names(DEG_FLUX)[1] <- cells[1]
  
}
Epi_c6_TM4SF1 <- DEG_FLUX[[1]]
Epi_c6_TM4SF1$M_name <- paste(Epi_c6_TM4SF1$Compound_IN_name,"\u2192",Epi_c6_TM4SF1$Compound_OUT_name)

DEG_FLUX <- list(Epi_c6_TM4SF1)
names(DEG_FLUX) <- c("Epi_c6_TM4SF1")
adj_flux$celltype <- "Epi_c6_TM4SF1"
adj_flux$group <- adj_scRNA$celltype11

cohens_d <- function(x, y) {
  pooled_std <- sqrt(((length(x)-1) * var(x) + (length(y)-1) * var(y)) / (length(x) + length(y) - 2))
  return ((mean(x) - mean(y)) / pooled_std)
}

for (i in 1:length(DEG_FLUX)) {
  cells_flux <- adj_flux[adj_flux$celltype == names(DEG_FLUX)[1], ]
  cells_SD <- rownames(cells_flux)[cells_flux$group == "TM4SF1_Epi"]
  cells_HC <- rownames(cells_flux)[cells_flux$group == "Others"]
  
  df <- as.data.frame(t(cells_flux))
  df <- df[-c(169,170),]
  
  df_SD <- df[,cells_SD]
  df_HC <- df[,cells_HC]
  
  
  for (flux_id in rownames(DEG_FLUX[[1]])) {
    A <- as.numeric(df_SD[flux_id, ])
    B <- as.numeric(df_HC[flux_id, ])
    c_d <- cohens_d(A, B)
    DEG_FLUX[[1]][flux_id, 'cohens_d'] <- c_d
  }
  
}
library(ggplot2)
library(ggrepel)
library(cowplot)
zero_count <- sum(-log10(DEG_FLUX[[1]]$p_val_adj) > 285)
random_pvalues <- runif(zero_count, min=290, max=320)
DEG_FLUX[[1]]$p_val_adj[-log10(DEG_FLUX[[1]]$p_val_adj) > 285] <- 10^(-random_pvalues)
print(DEG_FLUX[[1]]$p_val_adj)
p = ggplot(DEG_FLUX[[1]], aes(x=cohens_d, y=-log10(p_val_adj))) +
  geom_hline(aes(yintercept=1.3), color = "#999999", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0), color = "#999999", linetype="dashed", size=1) + 
  geom_point(size=2,color="grey80")+
  geom_point(data = DEG_FLUX[[1]][DEG_FLUX[[1]]$SM_anno=="Pyrimidine_synthesis",], stroke = 0.5, size=2, shape=16, color="#B51F29") + 
  labs(x = "Cohen's D",y = "-Log10(P-value)", title = "") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.title =element_text(size = 14),axis.text =element_text(size = 8, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 12)) +
  theme(plot.margin=unit(c(0, 1, 2, 1),'cm'))+
  geom_text_repel(data=DEG_FLUX[[1]][DEG_FLUX[[1]]$SM_anno=="Pyrimidine_synthesis" & DEG_FLUX[[1]]$p_val_adj<=0.05 & abs(DEG_FLUX[[1]]$cohens_d)>0.2, ], 
                  aes(label=M_name), color="black", size=4, 
                  arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3,
                  max.overlaps = Inf)+
  ggtitle("TM4SF1_Epi \nPyrimidine_synthesis")

ggdraw(xlim = c(0, 1), ylim = c(0,1.1))+ 
  draw_plot(p,x = 0, y =0) +  
  draw_line(x = c(0.55,0.75), 
            y = c(0.1,0.1),
            lineend = "round",
            size =1, col = "#B51F29",  
            arrow=arrow(angle = 15, 
                        length = unit(0.1,"inches"),
                        type = "closed"))+
  draw_line(x = c(0.25,0.45), 
            y = c(0.1,0.1),
            lineend = "round",
            size =1, col = "#006699",  
            arrow=arrow(angle = 15, 
                        length = unit(0.1,"inches"),
                        type = "closed",
                        ends="first"))+
  draw_text(text = "Activate in TM4SF1_Epi", size = 12,
            x = 0.88, y = 0.1,
            color="black",fontface = "italic")+
  draw_text(text = "Activate in Others", size = 12,
            x = 0.15, y = 0.1,
            color="black",fontface = "italic")



library(ggrepel)
DEG_FLUX[[1]] <- DEG_FLUX[[1]] %>%
  mutate(Difference = pct.1 - pct.2) %>% 
  rownames_to_column("gene")

ggplot(DEG_FLUX[[1]], aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=1.2) + 
  geom_label_repel(data=subset(DEG_FLUX[[1]], 
                               avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05), 
                   aes(label=gene),  
                   color="black", 
                   segment.colour = "black",
                   label.padding = 0.2, 
                   segment.size = 0.3,  
                   size=4,
                   max.overlaps = 20) +  
  geom_label_repel(data=subset(DEG_FLUX[[1]], 
                               avg_log2FC <= -1 & Difference <= -0.1 & p_val_adj <= 0.05), 
                   aes(label=gene), 
                   label.padding = 0.2, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4,
                   max.overlaps = 20) + 
  geom_vline(xintercept = 0,linetype = 2) +
  geom_hline(yintercept = 0,linetype = 2) +
  labs(x="△ Percentage difference",y="Log-Fold change") + 
  theme_bw()

p = ggplot(DEG_FLUX[[1]], aes(x=cohens_d, y=-log10(p_val_adj))) +
  geom_hline(aes(yintercept=0), color = "#999999", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0), color = "#999999", linetype="dashed", size=1) + 
  geom_point(size=2,color="grey80")+
  geom_point(data = DEG_FLUX[[1]][DEG_FLUX[[1]]$SM_anno=="Glycolysis_TCA_cycle",], stroke = 0.5, size=2, shape=16, color="#B51F29") + 
  labs(x = "Cohen's D",y = "△ Percentage difference", title = "") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size=1, colour = "black"),
        axis.title =element_text(size = 14),axis.text =element_text(size = 8, color = "black"),
        plot.title =element_text(hjust = 0.5, size = 12)) +
  theme(plot.margin=unit(c(0, 1, 2, 1),'cm'))+
  geom_text_repel(data=DEG_FLUX[[1]][DEG_FLUX[[1]]$SM_anno=="Glycolysis_TCA_cycle" & DEG_FLUX[[1]]$p_val_adj<=0.05 & abs(DEG_FLUX[[1]]$cohens_d)>0.2, ], 
                  aes(label=M_name), color="black", size=4, 
                  arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                  point.padding = 0.3, segment.color = 'black', segment.size = 0.3, force = 1, max.iter = 3e3,
                  max.overlaps = Inf)+
  ggtitle("TM4SF1_Epi \nGlycolysis_TCA_cycle")

