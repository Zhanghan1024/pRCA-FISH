# pciseq analysis
import os
import ast
import numpy as np
import pandas as pd
import skimage.color
from aicsimageio import AICSImage
from cellpose import models, io, plot, utils
from cellpose.io import imread
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import load_npz, coo_matrix
from skimage.segmentation import mark_boundaries
import pciSeq
import seaborn as sns
# Load segmentation out
coo = load_npz('segmentation.npz')
print('The image has %d cells' % len(set(coo.data)))
# Load fluorescence signal, calculate by U-FISH workflow
iss_spots = pd.read_csv('exp-full-genename.csv')
iss = iss_spots.rename(columns={
    "gene": "Gene",
    "dim_2": "x",
    "dim_1": "y"
})[["Gene", "x", "y"]]
# Load scRNAseq counts matrix
# R: remove cell_cycle cells
# ear <- readRDS(ear.rds)
# Idents(ear) <- ear$seurat_clusters
# ear_rm <- subset(ear, idents = setdiff(unique(Idents(ear)), c(0, 8, 14)))
# ear_rm$seurat_clusters <- paste0("celltype_", ear_rm$seurat_clusters)
# Idents(ear_rm) <- ear_rm$seurat_clusters
# expr_mat <- GetAssayData(ear_rm, layer = "counts")
# celltypes <- Idents(ear_rm)
# colnames(expr_mat) <- as.character(celltypes)
# pciSeq_scRNA <- as.data.frame(expr_mat)
# write.csv(pciSeq_scRNA, 'pciSeq_scRNA.csv')
scRNAseq = pd.read_csv('pciSeq_scRNA.csv', header = None, index_col = 0, dtype = object)
scRNAseq = scRNAseq.rename(columns=scRNAseq.iloc[0], copy = False).iloc[1:]
scRNAseq = scRNAseq.astype(float).astype(np.uint32)
#cell typing
pciSeq.setup_logger()
cellData, geneData = pciSeq.fit(spots = iss, coo = coo, scRNAseq = scRNAseq)
cellData['ClassName'] = cellData['ClassName'].apply(ast.literal_eval)
cellData = cellData[cellData['Genenames'].apply(len) > 0]
cellData['celltype'] = cellData['ClassName'].apply(lambda x: x[0])
final_df = cellData[cellData['celltype'].apply(lambda x: x[0] != 'Zero')]
celltype = cellData['ClassName'].apply(lambda x: x[0])
unique_classes = pd.unique(celltype)
num_classes = len(unique_classes)
palette = sns.color_palette("tab20", n_colors=num_classes)
palette[0] = "lightgray"
class_to_color = {cls: palette[i] for i, cls in enumerate(unique_classes)}
colors = [class_to_color[cls] for cls in celltype]
plt.figure(figsize=(2.8, 7.3))
plt.scatter(cellData['X'], cellData['Y'], alpha=0.5, s=3, c=colors)
legend_elements = []
for i, cls in enumerate(unique_classes):
    legend_elements.append(plt.Line2D(
        [0], [0],
        marker='o',
        color='w',
        markerfacecolor=palette[i],
        markersize=3,
        label=cls
    ))
plt.legend(handles=legend_elements, title='Cell Types', loc='upper left')
plt.title('Cell Distribution by Type')
plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.savefig('pciseq.pdf', dpi = 200)

color_dict = {cell: palette[i] for i, cell in enumerate(cellData['celltype'].unique())}
all_types = cellData['celltype'].unique()

for i, celltype in enumerate(all_types):
    plt.figure(figsize=(2.8, 7.3))
    bg = cellData[cellData['celltype'] != celltype]
    plt.scatter(
        bg['X'], bg['Y'],
        c="lightgray", alpha=0.5, s=3, label="other")
    target = cellData[cellData['celltype'] == celltype]
    plt.scatter(
        target['X'], target['Y'],
        c=[color_dict[celltype]], alpha=0.6, s=4, edgecolor="black", linewidth=0.2,label=celltype)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.legend(loc="best", fontsize=9)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f'{celltype}_highlight.pdf', dpi=200, bbox_inches='tight')
    plt.close()

cellData['Genenames'] = cellData['Genenames'].apply(ast.literal_eval)
cellData['Gene'] = cellData['Genenames'].apply(lambda x: x[0])
palette = sns.color_palette("husl", n_colors=len(cellData['Gene'].unique()))
color_dict = {gene: palette[i] for i, gene in enumerate(cellData['Gene'].unique())}
all_genes = cellData['Gene'].unique()

for i, geneID in enumerate(all_genes):
    plt.figure(figsize=(2.8, 7.3))
    bg = cellData[cellData['Gene'] != geneID]
    plt.scatter(
        bg['X'], bg['Y'],
        c="lightgray", alpha=0.5, s=3)
    target = cellData[cellData['Gene'] == geneID]
    plt.scatter(
        target['X'], target['Y'],
        c=[color_dict[geneID]], alpha=0.5, s=1.5, edgecolor=None, linewidth=0.2,label=geneID)
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.legend(loc="best", fontsize=9)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(f'{geneID}.png', dpi=200, bbox_inches='tight')
    plt.close()