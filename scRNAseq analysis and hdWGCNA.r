library(harmony)
library(Seurat)
library(ggplot2)
library(ggsci)
library(clustree)
library(cowplot)
library(dplyr)
library(tidydr)
# Load 10x scRNAseq files
# https://doi.org/10.1016/j.devcel.2020.12.015

filedir = c("ear1", "ear2", "ear3")
ear = list()
for (i in 1:length(filedir)){
    counts <- Read10X(paste0("../", filedir[i]), gene.column = 1)
    seuratobj = CreateSeuratObject(counts = counts, project = filedir[i], min.cells = 3)
    ear[[i]] = subset(seuratobj, subset= nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA < 15000 & nCount_RNA > 1000)
}
obj <- merge(x = ear[[1]], y = ear[-1], project = "integrated")
obj[["percent.PRG"]] <- PercentageFeatureSet(object = obj, features = protoplast_genes)
obj <- obj %>% NormalizeData(verbose = FALSE)
VariableFeatures(obj) <- split(row.names(obj@meta.data), obj$orig.ident) %>% lapply(function(cells_use) {
    obj[,cells_use] %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 1000) %>% 
        VariableFeatures()
}) %>% unlist %>% unique
# Scale protoplast responses genes  
obj <- obj %>% 
    ScaleData(vars.to.regress='percent.PRG') %>%      
    RunPCA(features = VariableFeatures(obj), npcs = 50)
# Harmony integrated multiple samples
obj <- obj %>% 
    RunHarmony("orig.ident", plot_convergence = FALSE, nclust = 50, max_iter = 10, early_stop = T)

p1 <- DimPlot(object = obj, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = obj, features = "harmony_1", group.by = "orig.ident",  pt.size = .1)
plot_grid(p1,p2) -> p
ggsave(p, filename='ear_harmony.png', width=12, height= 5)
# Resolutions check
obj = FindNeighbors(obj, reduction = "harmony", dims = 1:50)
obj = FindClusters(object = obj, resolution = c(seq(0.3, 0.8, by=0.1)))
clustree(obj@meta.data, prefix = 'RNA_snn_res.') -> p
ggsave(p, filename = './ear_clustree.png', width=12, height =15)
obj <- obj %>%
    RunUMAP(reduction = "harmony",  dims = 1:50)

p1 <- DimPlot(obj, reduction = "umap", group.by = "orig.ident", pt.size = .1)
p2 <- DimPlot(obj, reduction = "umap", label = TRUE,  pt.size = .1)
plot_grid(p1, p2) -> p
ggsave(p, filename='ear_umap.png', width=12, height= 5)

col_tab20 <- c(
  "#1f77b4",
  "#17becf",
  "#d62728",
  "#9467bd",
  "#8c564b",
  "#e377c2",
  "#ff7f0e",
  "#aec7e8",
  "#ffbb78",
  "#98df8a",
  "#ff9896",
  "#c5b0d5",
  "#c49c94",
  "#9edae5",
  "#003399"
)

DimPlot(obj, reduction='umap', label=T, pt.size = 0.1, cols = col_tab20) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
           labs(x = "UMAP1", y = "UMAP2") + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03)) -> p

ggsave(p, filename='ear_umap-2.png', width=8, height= 6)
# Find marker genes in all clusters
obj = JoinLayers(obj)
markers = FindAllMarkers(obj, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
# Dotplot visualization
top3 = markers %>% group_by(cluster) %>% top_n(3, wt = avg_log2FC)
p <- DotPlot(obj, features = unique(top3$gene), assay='RNA' ) +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend("Percent Expression"))+
  scale_color_gradientn(colours = c('#330066', '#336699', '#66CC66', '#FFCC33'))
ggsave(p, filename = 'Dotplot.pdf', width = 10, height =8)
write.table(markers, file='./ear_markers.txt', quote=F, append=F, sep='\t')
saveRDS(obj, file = 'ear.rds')

# hdWGCNA analysis
library(hdWGCNA)
library(harmony)
library(cowplot)
library(patchwork)
enableWGCNAThreads(nThreads = 10)
seurat_obj = readRDS('ear.rds')
seurat_obj <- RenameIdents(seurat_obj,
                "0" = "cellcycle",
                "1" = "epidermis",
                "2" = "lateral_organ",
                "3" = "rib_meristem",
                "4" = "meristem_internal",
                "5" = "pith",
                "6" = "meristem_periphery",
                "7" = "meristem_tip",
                "8" = "cellcycle",
                "9" = "unknown1",
                "10" = "xylem",
                "11" = "unknown2",
                "12" = "xylem",
                "13" = "epidermis",
                "14" = "cellcycle"
)
seurat_obj$cell_type <- Idents(seurat_obj)
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select = "fraction", fraction = 0.01, wgcna_name = "ear")
seurat_obj <- MetacellsByGroups(seurat_obj = seurat_obj, group.by = "cell_type", k = 25, max_shared = 15, ident.group = 'cell_type', reduction = 'harmony')
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(seurat_obj, 
                group_name = c("cellcycle", "epidermis", "lateral_organ", "rib_meristem", "meristem_internal", "pith","meristem_periphery", "meristem_tip", "unknown1","xylem", "unknown2"), 
                group.by='cell_type', assay = 'RNA', layer = 'data')
seurat_obj <- TestSoftPowers(seurat_obj, networkType = 'unsigned')
seurat_obj <- ConstructNetwork(seurat_obj, tom_name = 'ear', networkType = "unsigned", TOMType = "unsigned", overwrite_tom=TRUE, soft_power=2)
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars="cell_type")
hMEs <- GetMEs(seurat_obj)
MEs <- GetMEs(seurat_obj, harmonized=FALSE)
seurat_obj <- ModuleConnectivity(seurat_obj, 
                group_name = c("cellcycle", "epidermis", "lateral_organ", "rib_meristem", "meristem_internal", "pith","meristem_periphery", "meristem_tip", "unknown1","xylem", "unknown2"), 
                group.by='cell_type')
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
write.table(modules, file="modules.txt",quote=F, sep="\t", row.names=T, col.names=T)
hub_df <- GetHubGenes(seurat_obj, n_hubs = 20)
write.table(hub_df, file="hub_genes.txt",quote=F, sep="\t", row.names=T, col.names=T)
saveRDS(seurat_obj, file = 'ear_hdWGCNA.rds')
# module visualization
ModuleNetworkPlot(
    seurat_obj, 
    outdir='ModuleNetworks',
    n_inner = 20,
    n_outer = 30,
    n_conns = Inf,
    plot_size=c(10,10),
    vertex.label.cex=1
)