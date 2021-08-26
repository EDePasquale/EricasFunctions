#' Random Forest
#'
#' This function allows you to classify cells in the list of datasets using the reference with the randomForest package
#' @param reference Path to a Seurat object that is an integrated reference
#' @param datasets List of directories containing a Seurat object Example single directory path is "~/Documents/Projects/Data_files_temp/Donor/SBM1014/Seurat.rds"
#' @param results_name The name you would prefer to call the output directory. Default is "results"
#' @param ordering_vector A numeric vector containing the order you would like your clusters to appear in the plot, starting with 1. Default is "none", which retains the original order.
#' @param cluster_names A character vector containing the names of the clusters, in the original order. Default is "none" or cluster numbers to be used.
#' @param cluster_colors A character vector containing hex codes for the colors of cluster, in the original order. Default is "none" for random colors.
#' @keywords random forest classification
#' @export
#' @examples
#' random_forest()

random_forest <- function(reference, datasets, results_name="results", ordering_vector="none", cluster_names="none", cluster_colors="none"){
  
  # Load required packages
  library(Seurat)
  library(dplyr)
  library(randomcoloR)
  library(gplots)
  library(gdata)
  library(randomForest)
  library(tidyverse)
  library(data.table)
  
  # Read in reference
  ref=readRDS(reference)
  
  # Define cluster names
  if(cluster_names!="none"){
    clus_names=as.data.frame(cbind(0:(length(unique(ref@active.ident))-1), cluster_names))
  }else{
    clus_names=as.data.frame(cbind(0:(length(unique(ref@active.ident))-1), 0:(length(unique(ref@active.ident))-1)))
  }
  
  # Make Seurat's count matrix easier to access (log, transcript per 10K)
  bm.cm=as.matrix(ref@assays[["integrated"]]@data)
  
  # Create tibble with information on BPDCN samples that were analyzed with 10X
  TenX_samples.tib=NULL
  for(i in datasets){
    name=basename(sub("/[^/]*$", "", i)) #"SBM1064"
    TenX_samples.tib=rbind(TenX_samples.tib, c(name, i))
  }
  colnames(TenX_samples.tib)=c("Patient_ID", "Rds")
  TenX_samples.tib=as_tibble(TenX_samples.tib)
  
  # Create list of processed Seurat objects
  seu_10X.ls <- vector(mode = "list", length = nrow(TenX_samples.tib))
  for (x in 1:nrow(TenX_samples.tib)) {
    seu_10X.ls[[x]] <- readRDS(as.character(TenX_samples.tib[x,2]))
  }
  names(seu_10X.ls) <- TenX_samples.tib$Patient_ID
  cm.ls <- lapply(seu_10X.ls, function(x) as.matrix(x@assays[["RNA"]]@data))
  
  # Merge all gene expression data & check if it all makes sense.
  all(colnames(bm.cm) == names(ref@active.ident))
  
  # Limit to common genes
  select_genes_common=Reduce(intersect, c(list(row.names(bm.cm)), lapply(cm.ls, row.names)))
  
  # Create results directory
  parent.dir=sub("/[^/]*$", "", reference) #"~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem"
  dir.create(paste0(parent.dir, "/", results_name))
  setwd(paste0(parent.dir, "/", results_name))
  
  #==================================
  # Build classifier
  #==================================
  
  # Plant classification trees
  set.seed(123)
  
  rf <- randomForest(x = t(bm.cm[select_genes_common,]),    # matrix of predictors
                     y = ref@active.ident, # response vector
                     sampsize = rep(50, length(levels(ref@active.ident))),
                     ntree = 1000,
                     do.trace = 100)
  
  # Plot confusion matrix (based on out-of-bag data), with colors normalized for the total cell number in each population
  Conf.mat <<- rf$confusion[, rownames(rf$confusion)]
  NormConf.mat <<- Conf.mat / rowSums(Conf.mat)
  
  pdf(paste0("1_ConfusionMatrix.pdf"), width = 8, height = 8)
  heatmap.2(NormConf.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
            col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
            colsep = c(0, ncol(NormConf.mat)), rowsep = c(0, nrow(NormConf.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
            main = paste0("Confusion matrix, ", round(sum(diag(Conf.mat)) / sum(Conf.mat)*100, 2), "% accurate"),
            add.expr = text(rep(1:ncol(NormConf.mat), each=nrow(NormConf.mat)),
                            rep(ncol(NormConf.mat):1, nrow(NormConf.mat)), Conf.mat))
  dev.off()
  
  
  #==================================
  # Five-fold cross-validation
  #==================================
  
  # Split dataset in five parts
  cv <- split(colnames(bm.cm), rep(1:5, 1E6)[1:ncol(bm.cm)])
  
  # Build five forests, each with 4/5 of the data
  rf.cv <- lapply(cv, function(n) {
    set.seed(123)
    randomForest(x = t(bm.cm[select_genes_common,setdiff(colnames(bm.cm), n)]),
                 y = ref@active.ident[! colnames(bm.cm) %in% n],
                 sampsize = sapply(table(ref@active.ident[! colnames(bm.cm) %in% n]), min, 50), # check that it's not too low
                 ntree = 1000,
                 do.trace = 100)
  })
  
  # Predict the sets that were not used for training
  rf.cv.predict.prob <- lapply(rf.cv, function(rf) {
    predict(rf, t(bm.cm[select_genes_common, setdiff(colnames(bm.cm), names(rf$y))]), type = "prob")
  })
  
  # Maximum probability
  rf.cv.predict <- lapply(rf.cv.predict.prob, function(x) {
    y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))  # maximum sore
    names(y) <- rownames(x)
    y
  })
  
  # Confusion matrix
  Conf.cv.mat <<- table(ref@active.ident[unlist(lapply(rf.cv.predict, names))], unlist(rf.cv.predict))
  NormConf.cv.mat <<- Conf.cv.mat / rowSums(Conf.cv.mat)
  
  pdf(paste0("2_CrossValidation.pdf"), width = 8, height = 8)
  heatmap.2(NormConf.cv.mat, Rowv = F, Colv = F, dendrogram = "none", scale = "none", zlim = c(0, 1),
            col = colCustom(seq(0, 1, 0.01), color = c("white", "red")), trace = "none", density.info = "none",
            colsep = c(0, ncol(NormConf.cv.mat)), rowsep = c(0, nrow(NormConf.cv.mat)), sepcolor = "black",sepwidth = rep(0.01, 4),
            main = paste0("Confusion matrix, ", round(sum(diag(Conf.cv.mat)) / sum(Conf.cv.mat)*100, 2), "% accurate"),
            add.expr = text(rep(1:ncol(NormConf.cv.mat), each=nrow(NormConf.cv.mat)),
                            rep(ncol(NormConf.cv.mat):1, nrow(NormConf.cv.mat)), Conf.cv.mat))
  dev.off()
  
  #==================================
  # Classify cells from patients
  #==================================
  
  # Reduce to common genes
  bm.cm2 <- bm.cm[select_genes_common,]
  cm.ls2=lapply(cm.ls, function(x) x[select_genes_common, ])
  
  # Predict cell types in each of the BPDCN samples (ties are broken at random)
  predictions.mat.ls <- lapply(c(BM = list(bm.cm2), cm.ls), function(x) predict(rf, t(x[select_genes_common,]), type = "prob"))
  
  # Maximum probability
  CellTypes.ls <- lapply(predictions.mat.ls, function(x) {
    y <- factor(colnames(x)[apply(x, 1, which.max)], colnames(x))
    names(y) <- rownames(x)
    y
  })
  
  for(i in 1:length(CellTypes.ls)){
    write.table(CellTypes.ls[[i]], names(CellTypes.ls[i]), sep="\t", col.names = F, quote=F)
  }
  
  # Plot the clustered cell type frequencies in BM (per donor) and the predicted cell type frequencies in BPDCN
  BMfreq.mat <- do.call(cbind, lapply(split(ref@meta.data, f = cutf(ref@meta.data$orig.ident, d = "\\.")), function(x) table(x$CellType)))
  PredictFreq.mat <- do.call(cbind, lapply(CellTypes.ls, table))
  PlotFreq.mat<-PredictFreq.mat
  
  # Normalize to 100
  PlotFreqNorm.mat <- sweep(PlotFreq.mat, 2, colSums(PlotFreq.mat), "/")*100
  
  
  #==================================
  # Plot results
  #==================================
  
  # Write tabular results
  write.table(PlotFreqNorm.mat, "PlotFreqNorm.txt", sep="\t", quote=F)
  PlotFreqNorm.mat[,1]=as.numeric(row.names(PlotFreqNorm.mat))
  PlotFreqNorm.mat=as.data.frame(PlotFreqNorm.mat)
  PlotFreqNorm.mat=PlotFreqNorm.mat[order(PlotFreqNorm.mat[,1]),]
  PlotFreqNorm.mat=PlotFreqNorm.mat[,-1]
  
  # Reorder results
  if(ordering_vector!="none"){
    ordering_vec=ordering_vector
  }else{
    ordering_vec=c(1:length(unique(ref@active.ident)))
  }
  renamed=PlotFreqNorm.mat
  row.names(renamed)=1:length(unique(ref@active.ident))
  renamed=renamed[ordering_vec,]
  
  # Define cluster names
  if(cluster_names!="none"){
    clus_names=as.data.frame(cbind(0:(length(unique(ref@active.ident))-1), cluster_names))
  }else{
    clus_names=as.data.frame(cbind(0:(length(unique(ref@active.ident))-1), 0:(length(unique(ref@active.ident))-1)))
  }
  
  # Define cluster colors
  if(cluster_colors!="none"){
    CellTypeCol.ch=cluster_colors
  }else{
    CellTypeCol.ch=distinctColorPalette(length(unique(ref@active.ident)))
  }
  clus_names_col=cbind(clus_names, CellTypeCol.ch)
  renamed2=clus_names_col$CellTypeCol.ch
  renamed2=renamed2[ordering_vec]
  
  # Plot
  pdf(paste0("3_CellTypeFrequencies.pdf"), width = 8, height = 6)
  par(mar = c(8,4,8,12), xpd = T)
  
  barplot(as.matrix(renamed[nrow(renamed):1,]), col = rev(renamed2), xaxt = "n", ylab = "Population frequency (%)", border = NA)
  
  axis(side = 1, at = seq(1,ncol(renamed))*1.2-0.5, labels = colnames(renamed), las = 2)
  legend(x = ncol(renamed)*1.2+0.5, y = 100, legend = clus_names[,2][ordering_vec], fill = renamed2, bty = "n", border = NA)
  dev.off()
  
}