# Erica's Functions #

A comprehensive set of functions written by or frequently used by Erica DePasquale in her bioinformatics work. Use at your own risk!

# Installation #

Run the following code to install the package using devtools:

```
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github('EDePasquale/Ericas_Functions')
```

# Dependencies #
 
DoubletDecon requires the following R packages:
 
* Matrix
* Matrix.utils
* data.table
* gdata
* Rtsne
* umap
* dplyr
* DropletUtils
 
# Usage #
 
## qc_and_filter() ##
 
```javascript
qc_and_filter(directories, genome = "GRCh38", suffix = "default", min.umis=500, min.genes=250, max.mito=0.3, sig=T, dense=F)
```

#### Arguments ####

* directories - List of directories containing CellRanger output for QC and filtering. Example single directory path is "/Volumes/broad_vangalenlab/depasquale/CellRanger/210628_GeneCounts"
* genome - Desired genome: "GRCh38" and "GRCm38" are accepted, default is "GRCh38"
* suffix - Suffix for directory where the output of this function will be stored and where the input sample sheet is located. Default is the "QC" suffix.
* min.umis - QC filter for minimum number of Unique Molecular Identifiers (UMIs) detected per cell. Default is 500
* min.genes - QC filter for minimum number of genes detected per cell. Default is 250
* max.mito - QC filter to remove cells that exceed this fraction of mitochondrial alignments. Default is 0.3
* sig - Boolean for whether or not to plot cell type signatures. Default is T
* dense - Boolean for whether or not to save dense matrices alongside sparse matrices. Default is F

#### Value ####

No return value


## cell_signatures() ##

```javascript
cell_signatures(coords, expression, signatures, file.path.prefix, cexsize = 1, type = "UMAP", subset = 30:45)
```

#### Arguments ####

* coords - Dataframe containing 2 columns, each for X and Y coordinates with rows as cells
* expression - Dataframe of gene expression values, rows as genes and cells as columns
* signatures - The signatures file, with each column as a seperate list of genes (signature) with headers
* file.path.prefix - File path, including sample name, for plot to be saved as with the appropriate suffix. Example: "/Users/Documents/BPDCN123" 
* cexsize - The value that determines "cex" in the plot. Generated as part of the QC and Filter workflow but can be supplied independently. Default is 1
* type - The type of plots to generate/the type of coordinates given. "UMAP" (default) and "tSNE" accepted
* subset - Range of numeric values indicating which 16 plots to generate as the seperate "selected" plots. Default is 30:45 for the Griffin signatures

#### Value ####

No return value
