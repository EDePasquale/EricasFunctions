#' scPred
#'
#' This function allows you to classify cells in the list of datasets using the reference with the scPred package
#' @param reference Path to a Seurat object that is an integrated reference
#' @param datasets List of directories containing a Seurat object Example single directory path is "~/Documents/Projects/Data_files_temp/Donor/SBM1014/Seurat.rds"
#' @param results_name The name you would prefer to call the output directory. Default is "results"
#' @param ordering_vector A numeric vector containing the order you would like your clusters to appear in the plot, starting with 1. Default is "none", which retains the original order.
#' @param cluster_names A character vector containing the names of the clusters, in the original order. Default is "none" or cluster numbers to be used.
#' @param cluster_colors A character vector containing hex codes for the colors of cluster, in the original order. Default is "none" for random colors.
#' @keywords scpred classification
#' @export
#' @examples
#' scpred()

scpred <- function(reference, datasets, results_name="results", ordering_vector="none", cluster_names="none", cluster_colors="none"){
  
  # Load required packages
  library(Seurat)
  library(scPred)
  library(dplyr)
  library(randomcoloR)
  library(gplots)
  
  # Read in reference and process
  ref=readRDS(reference)
  ref@meta.data[["new_clusters"]] <- ref@active.ident
  ref <- getFeatureSpace(ref, "new_clusters")
  ref <- trainModel(ref)
  get_probabilities(ref) %>% head()
  get_scpred(ref)
  
  # Classify datasets
  PlotFreq=data.frame(x=0:(length(unique(ref@active.ident))-1), y=0:(length(unique(ref@active.ident))-1))
  
  for(i in datasets){
    
    message(paste0("\nClassifying dataset: ", i))
    
    # read in query and save name
    query=readRDS(i) #"~/Documents/Projects/Data_files_temp/Donor/SBM1064/Seurat.rds"
    name=basename(sub("/[^/]*$", "", i)) #"SBM1064"
    parent.dir=sub("/[^/]*$", "", reference) #"~/Documents/Projects/Data_files_temp/Donor/Seurat_Integration_0.5_ccrem"
    
    # make predictions
    query <- scPredict(query, ref, threshold=0)
    query <- RunUMAP(query, reduction = "scpred", dims = 1:30)
    
    #assign(name, query@meta.data[["predicted.id"]])
    assign(name,  table(query@meta.data[["scpred_prediction"]]))
    
    # add to PlotFreq table
    temp1=as.data.frame(table(query@meta.data[["scpred_prediction"]]))
    temp1[,1]=as.numeric(as.character(temp1[,1]))
    colnames(temp1)[2]=name
    PlotFreq=left_join(PlotFreq, temp1, by=c("x"="Var1"))
    
  }
  
  message("\nCreating plots")
  
  # Clean up PlotFreq table
  PlotFreq=PlotFreq[,-(1:2)]
  PlotFreq[is.na(PlotFreq)] = 0
  PlotFreqNorm.mat <- sweep(PlotFreq, 2, colSums(PlotFreq), "/")*100
  
  # Create results directory
  dir.create(paste0(parent.dir, "/", results_name))
  setwd(paste0(parent.dir, "/", results_name))
  
  # Write tabular results
  write.table(PlotFreqNorm.mat, "PlotFreqNorm.txt", sep="\t", quote=F)
  
  # Save training plot
  pdf(paste0("1_plot_probabilities.pdf"), width = 8, height = 6)
  plot1=plot_probabilities(ref) 
  par(mar = c(8,4,8,12), xpd = T)
    print(plot1)
  dev.off()
  
  
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