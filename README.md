## Script to reproduce the specificity score analysis in [Javierre, Burren, Wilder, Kreuzhuber, Hill et al. Cell 167, 1369-1384 (2016)](http://www.cell.com/cell/fulltext/S0092-8674(16)31322-8).

### Overview
The specificity score analysis described in the paper quantifies the cell type-specificity of each gene’s interactions with active enhancers, or each gene's expression, through calculation of gene specificity scores (see the paper for full details).

The script `specificityScoreAnalysis.R` contains:
- a function `specificityScore` that calculates cell type-specificity scores given a vector of values (one value for each cell type);
- code to calculate gene specificity scores (based on PCHi-C data and expression data);
- code to reproduce Figures 4B,C,D,E and S4B,C.

### Required R packages
- `gplots`
- `ggplot2`
- `RColorBrewer`
- `gridExtra`

### Required data files
The script reads in the following data files:
**PCHiC_peak_matrix_cutoff5.txt**
Available to download via the Open Science Framework
