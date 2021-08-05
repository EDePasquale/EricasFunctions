#' Cell_Signatures
#'
#' This function allows you to provide dimensionality coordinates and a gene by cell matrix to get cell signature overlays.
#' @param coords Dataframe containing 2 columns, each for X and Y coordinates with rows as cells
#' @param expression Dataframe of gene expression values, rows as genes and cells as columns
#' @param signatures The signatures file, with each column as a seperate list of genes (signature) with headers
#' @param file.path.prefix File path, including sample name, for plot to be saved as with the appropriate suffix. Example: "/Users/Documents/BPDCN123" 
#' @param cexsize The value that determines "cex" in the plot. Generated as part of the QC and Filter workflow but can be supplied independently. Default is 1
#' @param type The type of plots to generate/the type of coordinates given. "UMAP" (default) and "tSNE" accepted
#' @param subset Range of numeric values indicating which 16 plots to generate as the seperate "selected" plots. Default is 30:45 for the Griffin signatures
#' @keywords cell signature UMAP tSNE plot
#' @export
#' @examples
#' cell_signatures()

cell_signatures <- function(coords, expression, signatures, file.path.prefix, cexsize = 1, type = "UMAP", subset = 30:45){
  
  # Full signatures
  signatures <- lapply(as.list(signatures), function(y) y[!is.na(y)])
  signatures <- lapply(signatures, intersect, rownames(expression))
  exp.mean <- rowMeans(expression)
  signScore <- lapply(names(signatures), function(g) {
    message(g)
    scoreSignature(CM = as.matrix(expression), signatures = signatures[[g]], CM.mean = exp.mean, verbose = T)
  })
  names(signScore) <- names(signatures)
  
  # Plot
  pdf(file = paste0(file.path.prefix, "_", type, ".pdf"), width = 6, height = 6)
  par(mar=c(4, 4, 4, 4))
  for (n in names(signScore)) {
    mycol <- colItay(signScore[[n]])
    if(type == "UMAP"){
      plotUMAP(coords, pch = 16, cex = cexsize, col = mycol, main = n)
    }else if (type == "tSNE"){
      plotTSNE(coords, pch = 16, cex = cexsize, col = mycol, main = n)
    }else{
      message("This function only supports UMAP and tSNE at this time")
    }
  }
  dev.off()
  
  # Bonus plot of selected signatures (Griffin)
  pdf(file = paste0(file.path.prefix, "_", type, "_selected.pdf"), width = 11, height = 8.5)
  par(mar=c(2, 2, 2, 2))
  
  layout(matrix(1:16, ncol=4, nrow=4, byrow=TRUE))
  
  for(s in 30:45){
    mycol <- colItay(signScore[[s]])
    plot(coords, xlab = NA, ylab = NA, tck = F, yaxt = "n", xaxt = "n", pch = 16, cex = cexsize, col = mycol, main = names(signatures)[s])
  }
  dev.off()
  
}

