This repository is for development of multi-resolution heatmapping in Galaxy.

The multi-resolution heatmap (MRH) file format is as follows:
Header: 32 bytes
Lowest resolution in MRH - int32
Highest resolution in MRH - int32
Zoom factor - int32
Minimum observation cutoff - int32
First bin start coordinate - int32
Lowest resolution nuber of bins - int32
Total number of data bins - int32
Total number of bins (data + indices) - int32

Data bins: 1D float32 array
complete upper triangle 1D row-major array (including diagonal) of lowest resolution data
sets of zoom ** 2 or zoom * (zoom + 1) / 2 bins for each off-diagonal or on-diagonal lower-resolution bin, respectively if at least on higher resolution bin has sufficient numbers of reads.

Indices: 1D int32 array
One bin for every data bin in the same order, excepting the highest resolution data bins. These contain the index in the data array corresponding to the first bin in the subdivided next highest resolution bins corresponding to the region covered by the current bin.