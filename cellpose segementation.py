# cellpose v2 segementation
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
# Loading the DAPI channel image
img = imread('3D-20Xtile_405.tif')
img.shape #(11564, 3972, 3)
model = models.Cellpose(gpu=True, model_type='nuclei')
masks, flows, styles, _ = model.eval(img, diameter = 15, flow_threshold = 0.4, cellprob_threshold = -2.0, min_size = 15)
img_plot = np.transpose(img, (2, 0, 1))
img_gray = img_plot[0]
overlay = mark_boundaries(img_gray, masks, color=(1, 0, 0))
plt.imsave('segmentation_overlay.pdf', overlay)
# Writing segmentation to disk
segmentation = sparse.coo_matrix(masks)
sparse.save_npz("segmentation.npz", segmentation)
# FM (a small image) segmentation
img = imread('ear-FM-3D_AF405-DAPI.tif')
masks, flows, styles, _ = model.eval(img, diameter = 15, flow_threshold = 0.4, cellprob_threshold = -2.0, min_size = 15)
fig = plt.figure(figsize=(12,5))
plot.show_segmentation(fig, img, masks, flows[0])
plt.tight_layout()
plt.savefig("FM_segmentation.pdf", bbox_inches='tight')
plt.close(fig)