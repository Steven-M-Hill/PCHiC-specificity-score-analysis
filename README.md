## Script to reproduce the specificity score analysis in [Javierre, Burren, Wilder, Kreuzhuber, Hill et al. Cell 167, 1369-1384 (2016)](http://www.cell.com/cell/fulltext/S0092-8674(16)31322-8).

### Overview
The specificity score analysis described in the paper quantifies the cell type-specificity of each gene’s interactions with active enhancers, or each gene's expression, through calculation of gene specificity scores (see the paper for full details).

The script `specificityScoreAnalysis.R` contains:
- a function `specificityScore` that calculates cell type-specificity scores given a vector of values (one value for each cell type);
- code to calculate gene specificity scores (based on Promoter Capture Hi-C data and expression data);
- code to reproduce Figures 4B,C,D,E and S4B,C.

### Required R packages
- `gplots`
- `ggplot2`
- `RColorBrewer`
- `gridExtra`

### Required data files
The script reads in the following data files:

**PCHiC_peak_matrix_cutoff5.txt**<br />
Available to download from the Open Science Framework repository associated with the paper.<br />
Direct link to file: https://osf.io/63hh4/.<br />
This file contains the Promoter Capture Hi-C peak matrix consisting of CHiCAGO scores for all interactions that pass a cutoff of >=5 in at least one cell type. For further details of the contents of formatting of this file, see https://osf.io/cn4k8/.

**GeneExpressionMatrix.txt**<br />
Available to download from the Open Science Framework repository associated with the paper.<br />
Direct link to file: https://osf.io/wpjy8/.<br />
This file contains a matrix of gene expression data used in the study, containing expression quantifications generated with `MMSEQ` (Turro et al., Genome Biol 2011).

**PIRactivity.Rds**<br />
Available to download from this Github repository.
This R data file contains a data frame providing the activity statuses for the promoter interacting region (PIR) of each interaction, in each cell type. These activities are defined on the basis of chromHMM segmentations of BLUEPRINT histone modification ChIP data.

**baitAnnotations.Rds**<br />
Available to download from this Github repository.
This R data file contains a data frame providing Ensembl Regulatory Build features mapping to each baited promoter fragment.
