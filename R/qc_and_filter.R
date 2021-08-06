#' QC_and_Filter
#'
#' This function allows you to provide a list of paths leading to the QC sample sheet (Erica's pipeline Step 4) to perform quality control and filtering.
#' @param directories List of directories containing CellRanger output for QC and filtering. Example single directory path is "/Volumes/broad_vangalenlab/depasquale/CellRanger/210628_GeneCounts"
#' @param genome Desired genome: "GRCh38" and "GRCm38" are accepted, default is "GRCh38"
#' @param suffix Suffix for directory where the output of this function will be stored and where the input sample sheet is located. Default is the "QC" suffix.
#' @param min.umis QC filter for minimum number of Unique Molecular Identifiers (UMIs) detected per cell. Default is 500
#' @param min.genes QC filter for minimum number of genes detected per cell. Default is 250
#' @param max.mito QC filter to remove cells that exceed this fraction of mitochondrial alignments. Default is 0.3
#' @param sig Boolean for whether or not to plot cell type signatures. Default is T
#' @param dense Boolean for whether or not to save dense matrices alongside sparse matrices. Default is F
#' @keywords qc filter
#' @export
#' @examples
#' qc_and_filter()

qc_and_filter <- function(directories, genome = "GRCh38", suffix = "default", min.umis=500, min.genes=250, max.mito=0.3, sig=T, dense=F){
  
  # Options
  options(max.print = 500)
  options(stringsAsFactors = FALSE)
  options(scipen = 999)
  
  # Load libraries
  library(Matrix)
  library(Matrix.utils)
  library(data.table)
  library(gdata)
  library(Rtsne)
  library(umap)
  library(dplyr)
  library(DropletUtils)
  
  # Source files 
  data(popcol)
  if (sig == T) {data(signatures)}
  if (genome == "GRCh38") {
    data(geneInfoH)
  } else {
    data(geneInfoM)
  }
  
  for(i in 1:length(directories)){
    
    #=================================================
    # Set Up
    #=================================================
    
    # Progress update
    parent.dir=sub("_[^_]+$", "", directories[i])
    name=gsub(".*/", "", parent.dir)
    print(paste0("Running experiment ", name))
    print(paste0("...loading files ", name))
          
    # Other parameters
    if(suffix == "default"){
      out.dir=paste0(parent.dir, "_QC")
    }else{
      out.dir=paste0(parent.dir, "_", suffix)  
    }
    setwd(out.dir) #folder that contains outs
    
    # Set up
    sampleSheet <- read.csv(file = paste0(out.dir, "/", list.files(path = out.dir)[grepl(".csv$", list.files(path = out.dir))]))
    sampleSheet$WellCol <- popcol$mycol[1:nrow(sampleSheet)]
    sampleSheet$DonorCol <- popcol$mycol[as.numeric(as.factor(sampleSheet$Donor))+6]
   
    #=================================================
    # Load count matrix from CellRanger
    #=================================================
    
    samples=NULL
    sample_matrices=list(rep(NA, nrow(sampleSheet)))
    for(x in 1:nrow(sampleSheet)){
      print(paste0("Loading matrix from ", sampleSheet[x,1], "/strainedCounts"))
      barcode.path = paste0(directories[i], "/", sampleSheet[x,1], "/strainedCounts/", "barcodes.tsv.gz")
      features.path = paste0(directories[i], "/", sampleSheet[x,1], "/strainedCounts/", "features.tsv.gz")
      matrix.path = paste0(directories[i], "/", sampleSheet[x,1], "/strainedCounts/", "matrix.mtx.gz")
      CM2=readMM(file = matrix.path)
      feature.names = read.delim(features.path, header = F)
      barcode.names = read.delim(barcode.path, header = F)
      colnames(CM2) = barcode.names$V1
      rownames(CM2) = cutf(feature.names$V2, d = "\\.", f = 1)
      CM1 <- aggregate.Matrix(CM2, groupings = rownames(CM2), fun = "sum")
      sample_matrices[x]=CM1
      rm(CM1, CM2, barcode.path, features.path, matrix.path, feature.names, barcode.names)
      samples=c(samples, sampleSheet[x,1])
    }
    names(sample_matrices)=paste0(samples, "_CM1")
    
    #=================================================
    # Make statistics file & QC plots
    #=================================================
    
    print(paste0("...QC plots ", name))
    
    # Make vectors with genes of interest
    if (genome == "GRCh38") {
      Ychr.ch <- geneInfoH$gene_name[geneInfoH$chr == "Y"]
      for(x in 2:length(sample_matrices)){       #check to make sure all gene names are the same
        if(all(sample_matrices[[x]]@Dimnames[[1]]!=sample_matrices[[x]]@Dimnames[[1]])){
          stop("gene names don't match") 
        }
      }
      Ychr.ch <- intersect(Ychr.ch, rownames(sample_matrices[[x]]))
      chrMgenes.ch <- geneInfoH$gene_name[geneInfoH$chr == "MT"]
    } else if (genome == "GRCm38") {
      Ychr.ch <- geneInfoM$gene_name[geneInfoM$chr == "chrY"]
      for(x in 2:length(sample_matrices)){       #check to make sure all gene names are the same
        if(all(sample_matrices[[x]]@Dimnames[[1]]!=sample_matrices[[x]]@Dimnames[[1]])){
          stop("gene names don't match") 
        }
      }
      Ychr.ch <- intersect(Ychr.ch, rownames(sample_matrices[[x]]))
      chrMgenes.ch <- geneInfoM$gene_name[geneInfoM$chr == "chrM"]
    }
    
    # Make file for plotting
    stats.dt=list(rep(NA, nrow(sampleSheet)))
    for(x in 1:length(sample_matrices)){
      print(dim(sample_matrices[[x]]))
      print(class(sample_matrices[[x]]))
      test<<-sample_matrices[[x]]
      print("umis")
      a<<-Matrix::colSums(sample_matrices[[x]])
      print("genes")
      b<<-colSums(sample_matrices[[x]] != 0)
      print("chrYfraction")
      c<<-colSums(sample_matrices[[x]][Ychr.ch,]) / colSums(sample_matrices[[x]])
      print("chrMfraction")
      d<<-colSums(sample_matrices[[x]][chrMgenes.ch,]) / colSums(sample_matrices[[x]])
      print("chrMcount")
      e<<-colSums(sample_matrices[[x]][chrMgenes.ch,])
      
      stats1.dt <- data.table(cell = colnames(sample_matrices[[x]]),
                              umis = colSums(sample_matrices[[x]]),
                              genes = colSums(sample_matrices[[x]] != 0),
                              chrYfraction = colSums(sample_matrices[[x]][Ychr.ch,]) / colSums(sample_matrices[[x]]),
                              chrMfraction = colSums(sample_matrices[[x]][chrMgenes.ch,]) / colSums(sample_matrices[[x]]),
                              chrMcount = colSums(sample_matrices[[x]][chrMgenes.ch,]),
                              wellcol = sampleSheet$WellCol[[x]],
                              donorcol = sampleSheet$DonorCol[[x]])
      setorder(stats1.dt, -umis)
      stats.dt[[x]]=stats1.dt
      rm(stats1.dt)
    }
    names(stats.dt)=paste0(samples, "_stats")
    
    # Define cells that pass umi and gene thresholds
    values=list(rep(NA, nrow(sampleSheet)))
    for(x in 1:nrow(sampleSheet)){
      umiGenePass.log <- stats.dt[[x]][,"genes"] >= min.genes & stats.dt[[x]][,"umis"] >= min.umis
      mitoPass.log <- stats.dt[[x]][,"chrMfraction"] <= max.mito
      mycol <- ifelse(umiGenePass.log, yes = ifelse(mitoPass.log, yes = "black", no = "red"), no = "grey")
      cexsize <- round(max(c(0.3, 0.5-nrow(stats.dt[[x]])/40000)),2)
      values[[x]]=list(umiGenePass.log=umiGenePass.log, mitoPass.log=mitoPass.log, mycol=mycol, cexsize=cexsize)
    }
    names(values)=paste0(samples, "_values")
    
    # Generate QC plots
    for(x in 1:nrow(sampleSheet)){
      
      print(paste0("...QC ", x, "/", nrow(sampleSheet), " ", name))
      pdf(file = paste0(out.dir, "/", samples[x], "_QC.pdf"), width = 6, height = 6)
      par(oma=c(0,0,2,0))
      par(mfrow=c(1,1),mar=c(4,4,4,4))
      
      
      # cells vs umis
      plot(stats.dt[[x]]$umis, pch = 16, cex = values[[x]]$cexsize, log = "xy", col = values[[x]]$mycol, ylim = c(1, max(stats.dt[[x]]$umis)),
           xlab = "Cells", ylab = "umis", main=samples[x])
      abline(h = min.umis, col = "black", lty = 2)
      title("Cell vs UMI", outer=TRUE)
      
      # cells vs genes
      plot(stats.dt[[x]]$genes, pch = 16, cex = values[[x]]$cexsize, log = "xy", col = values[[x]]$mycol, ylim = c(1, max(stats.dt[[x]]$genes)),
           xlab = "Cells", ylab = "genes", main=samples[x])
      abline(h = min.genes, col = "black", lty = 2)
      title("Cell vs Gene", outer=TRUE)
      
      # Fraction of reads that align to chrY (color by CM well)
      plot(stats.dt[[x]]$chrYfraction, pch = 16, cex = values[[x]]$cexsize, col = stats.dt[[x]]$donorcol, 
           ylab = "Fraction chrY umis", xlab = "Cells", log = "x", main=samples[x])
      text(x=1, y=max(stats.dt[[x]]$chrYfraction)-max(stats.dt[[x]]$chrYfraction)*0.05, labels = sampleSheet$library_id[[x]], col = sampleSheet$DonorCol[[x]], pos = 4)
      title("Fraction of reads that align to chrY", outer=TRUE)
      
      # Fraction of reads that align to chrM
      plot(stats.dt[[x]]$chrMfraction,  pch = 16, cex = values[[x]]$cexsize, col = values[[x]]$mycol, ylab = "Fraction chrM umis", xlab = "Cells", log = "x", ylim = c(0, 1), main=samples[x])
      abline(h = max.mito, col = "black", lty = 2)
      title("Fraction of reads that align to chrM", outer=TRUE)
      
      # Total reads that align to chrM
      plot(stats.dt[[x]]$chrMcount, pch = 16, cex = values[[x]]$cexsize, col = values[[x]]$mycol, ylab = "chrM umi count", xlab = "Cells", log = "x", main=samples[x])
      title("Total reads that align to chrM", outer=TRUE)
      
      # Plot: genes vs UMIs
      plot(stats.dt[[x]][,c("umis", "genes")], pch = 16, cex = values[[x]]$cexsize, col = values[[x]]$mycol, main=samples[x])
      abline(v = min.umis, col = "black", lty = 2)
      abline(h = min.genes, col = "black", lty = 2)
      text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.6, labels = sum(values[[x]]$mycol == "black"), col = "black")
      text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.55, labels = sum(values[[x]]$mycol == "green"), col = "green")
      text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.5, labels = sum(values[[x]]$mycol == "red"), col = "red")
      text(x = max(stats.dt[[x]]$umis)*0.8, y = max(stats.dt[[x]]$genes)*0.45, labels = sum(values[[x]]$mycol == "grey"), col = "grey")
      title("Gene vs UMI", outer=TRUE)
      
      dev.off()
      
    }
    
    #=================================================
    # Filters
    #=================================================
    
    print(paste0("...filters ", name))
    sample_matrices_filtered=list(rep(NA, nrow(sampleSheet)))
    for(x in 1:nrow(sampleSheet)){
      
      # Filter cells by UMI and number of genes
      CM1=sample_matrices[[x]]
      message(paste0("\n", ncol(CM1), " cells and ", nrow(CM1), " genes in count matrix."))
      CM_XYM.dgm <- CM1[,stats.dt[[x]]$cell[values[[x]][["umiGenePass.log"]] == T]]
      message(paste0(ncol(CM_XYM.dgm), " cells and ", nrow(CM_XYM.dgm), " genes remaining after filtering for ", min.umis, " UMIs and ", min.genes, " genes."))
      
      # Filter cells based on chrM alignments
      CM.dgm <- CM_XYM.dgm[,stats.dt[[x]]$cell[values[[x]][["umiGenePass.log"]] == T & values[[x]][["mitoPass.log"]] == T]]
      CM.dgm <- CM.dgm[! rownames(CM.dgm) %in% chrMgenes.ch,]
      message(paste0(ncol(CM.dgm), " cells and ", nrow(CM.dgm), " genes remaining after filtering for cells with <", max.mito, " chrM alignments and chrM removal."))
      
      # Filter genes that map to X or Y
      if (genome == "GRCh38") {
        XYchr.ch <- geneInfoH$gene_name[geneInfoH$chr == "Y" | geneInfoH$chr == "X"]
      } else if (genome == "GRCm38") {
        XYchr.ch <- geneInfoM$gene_name[geneInfoM$chr == "chrY" | geneInfoM$chr == "chrX"]
      }
      CM.dgm <- CM.dgm[! rownames(CM.dgm) %in% XYchr.ch,]
      message(paste0(ncol(CM.dgm), " cells and ", nrow(CM.dgm), " genes remaining after chrXY removal."))
      
      # Subset stats table
      stats.dt[[x]] <- stats.dt[[x]][stats.dt[[x]]$cell %in% colnames(CM.dgm),]
      
      # Save CM.dgm for future use
      sample_matrices_filtered[x] <- CM.dgm
      
    }
    names(sample_matrices_filtered)=paste0(samples, "_filtered")
    
    #=================================================
    # Save sparse matrices (MM format)
    #=================================================
    
    #Save files - 1 per donor
    dir.create(paste0(out.dir, "/Sparse"))
    setwd(paste0(out.dir, "/Sparse"))
    
    for(x in 1:length(unique(sampleSheet$Donor))){
      
      # Combine expression matrices for (genes and order match between datasets)
      temp_i=which(sampleSheet$Donor == unique(sampleSheet$Donor)[x])
      temp_m=NULL
      for(j in temp_i){
        temp_m2=sample_matrices_filtered[[j]]
        colnames(temp_m2)=paste0(colnames(sample_matrices_filtered[[j]]), "_", j)
        temp_m=cbind(temp_m, temp_m2)
      }
      
      # Write out expression matrix (10x format)
      DropletUtils:::write10xCounts(paste0("./", unique(sampleSheet$Donor)[x], "_FLT3"), temp_m, version="3")
    }
    
    #=================================================
    # Save dense matrices (optional)
    #=================================================
    
    if(dense == T){
      
      #Save files - 1 per donor
      dir.create(paste0(out.dir, "/Dense"))
      setwd(paste0(out.dir, "/Dense"))
      
      for(x in 1:length(unique(sampleSheet$Donor))){
        
        # Combine expression matrices for (genes and order match between datasets)
        temp_i=which(sampleSheet$Donor == unique(sampleSheet$Donor)[x])
        temp_m=NULL
        for(j in temp_i){
          temp_m2=sample_matrices_filtered[[j]]
          colnames(temp_m2)=paste0(colnames(sample_matrices_filtered[[j]]), "_", j)
          temp_m=cbind(temp_m, temp_m2)
        }
        
        # Write out expression matrix
        write.table(as.matrix(temp_m), file=paste0(unique(sampleSheet$Donor)[x], "_dense.txt"), sep="\t")
      }
    }
    
    #Reset wd
    setwd(out.dir) #folder that contains outs
    
    #=================================================
    # Generate Dimensionality Reduction Plots
    #=================================================
    
    for(x in 1:nrow(sampleSheet)){
      
      ##################
      # Variable Genes #
      ##################
      
      print(paste0("...variable genes ", x, "/", nrow(sampleSheet), " ", name))
      
      # Normalize count matrix
      CM.norm <- t(t(as.matrix(sample_matrices_filtered[[x]]))/colSums(as.matrix(sample_matrices_filtered[[x]]))) * 10000   # Normalize count matrix to 10,000 transcripts per cell
      
      # The following is basically identical to the function variableGenes that is used for Seq-Well
      CM.norm.min <- 0.01
      sd.cutoff <- 1
      
      # Calculate mean and coefficient of variation per gene, remove genes that are not expressed at all
      CM.norm.mean <- rowMeans(CM.norm)
      CM.norm.cv <- apply(CM.norm, 1, sd) / CM.norm.mean
      CM.norm.mean <- sort(CM.norm.mean[!is.na(CM.norm.cv)])
      CM.norm.cv <- CM.norm.cv[names(CM.norm.mean)]
      
      # Plot
      message(paste("\nGenerating", paste0(out.dir, "/", samples[x], "_VarGenes.pdf")))
      pdf(file = paste0(out.dir, "/", samples[x], "_VarGenes.pdf"), width = 6, height = 6)
      par(mar=c(4,4,4,4))
      
      plot(CM.norm.mean, CM.norm.cv, pch = ".", log = "xy", xlab = "Mean expression", ylab = "Coefficient of variation")
      abline(v = CM.norm.min, col = "blue")
      
      # Linear fit
      CM.norm.lm <- lm(y~x, subset = CM.norm.mean>=CM.norm.min, data.frame(x = log10(CM.norm.mean), y = log10(CM.norm.cv)))
      CM.norm.cv.predict <- predict(CM.norm.lm, data.frame(x = log10(CM.norm.mean)))
      lines(CM.norm.mean, 10^CM.norm.cv.predict, col = "red")
      
      # This is from Volker's function, I do not have a complete understanding
      CM.norm.cv.residue <- log10(CM.norm.cv) - CM.norm.cv.predict
      plot(CM.norm.mean, CM.norm.cv.residue, log="x", pch=".", xlab="Mean", ylab="Residue")
      abline(v=CM.norm.min, col="green")
      abline(h=0, col="red")
      abline(h=sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])*c(1, 2, 3), col="red", lty=2)
      axis(4, at = round(sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])*c(1, 2, 3), 3), las = 2)
      
      plot(seq(-0.25+0.005, 1-0.005, 0.01), as.vector(table(cut(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min], seq(-0.25, 1, 0.01)))), xlab="Residue", ylab="Frequency", type="b")
      abline(v = 0, col = "red")
      abline(v = sd(CM.norm.cv.residue[CM.norm.mean >= CM.norm.min])*c(1, 2, 3), col="red", lty = 2)
      axis(3, at = round(sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])*c(1, 2, 3), 3), las = 2)
      
      plot(ecdf(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min]), main="ECDF")
      qqnorm(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])
      
      # This is the key object
      vargenes.num <- sort((CM.norm.cv.residue)[CM.norm.mean>=CM.norm.min], decreasing = T)
      attr(vargenes.num, which = "sd") <- sd(CM.norm.cv.residue[CM.norm.mean>=CM.norm.min])
      VarGenes1.ch <- names(which(vargenes.num > sd.cutoff*attr(vargenes.num, "sd")))
      # Remove cell cycle genes to generate VarGenes.ch
      VarGenes.ch <- setdiff(VarGenes1.ch, c("ASPM", "CENPE", "CENPF", "DLGAP5", "MKI67", "NUSAP1", "PCLAF", "STMN1", "TOP2A", "TUBB"))
      
      plot(CM.norm.mean, CM.norm.cv, log="xy", pch=".", col=ifelse(names(CM.norm.cv) %in% names(which(vargenes.num > sd.cutoff*attr(vargenes.num, which = "sd"))), "red", "grey"), xlab="Mean", ylab="CV")
      abline(v=CM.norm.min, col="blue")
      lines(CM.norm.mean, 10^CM.norm.cv.predict, col="red")
      text(x = 50, y = 20, labels = paste(length(VarGenes.ch), "genes"), col = "red")
      
      dev.off()
      
      ########
      # tSNE #
      ########
      
      print(paste0("...tSNE ", x, "/", nrow(sampleSheet), " ", name))
      
      # Adjustments to count matrix
      CM.norm.log <- log(CM.norm+1)                         # Normalize & log
      CM.norm.log.var <- CM.norm.log[VarGenes.ch,]          # Normalize & log and only keep variable genes 
      
      # Based on 2_tSNE.PvG190810.R and 1d_runTSNE.vh20190325.R
      message("\nCalculating tSNE coordinates")
      
      # Make tSNE with default values & add values to stats.dt
      set.seed(42)
      tryCatch(
        expr = {
          CM.tsne <- Rtsne(as.matrix(t(CM.norm.log.var)), max_iter = 1000, theta = 0.5, perplexity = 30, verbose = TRUE)
        },
        error = function(e){ 
          tryCatch(
            expr={
              CM.tsne <- Rtsne(as.matrix(t(CM.norm.log.var)), max_iter = 1000, theta = 0.5, perplexity = 5, verbose = TRUE)
            },
            error = function(e){
              CM.tsne <- Rtsne(as.matrix(t(CM.norm.log.var)), max_iter = 1000, theta = 0.5, perplexity = 1, verbose = TRUE)
            }
          )
        }
      )
      stopifnot( all(colnames(CM.norm.log) == stats.dt[[x]]$cell)) # just checking...
      stats.dt$tSNEx <- CM.tsne$Y[,1]
      stats.dt$tSNEy <- CM.tsne$Y[,2]
      
      message(paste("\nGenerating", paste0(out.dir, "/", samples[x], "_tSNE.pdf")))
      pdf(file = paste0(out.dir, "/", samples[x], "_tSNE.pdf"), width = 6, height = 6)
      par(mar=c(4, 4, 4, 4))
      
      # Color by UMI
      plotTSNE(cbind(stats.dt$tSNEx, stats.dt$tSNEy), pch = 16, cex = cexsize, col = colCustom(stats.dt[[x]]$umis, c(2000, 20000)),
               main = paste0("UMIs (", min(stats.dt[[x]]$umis), "-", max(stats.dt[[x]]$umis), ", mean ", round(mean(stats.dt[[x]]$umis), 0), ")"))
      axis(side = 1, at = 0, label = paste(nrow(stats.dt[[x]]), "cells"), tick = F)
      
      # Color by Donor
      plotTSNE(cbind(stats.dt$tSNEx, stats.dt$tSNEy), cex = cexsize, col = stats.dt[[x]]$donorcol, main = "Donor")
      sapply(1:length(unique(sampleSheet$Donor)), function(z) axis(side = 1, at = ifelse(length(unique(sampleSheet$Donor)) == 1,
             yes = 0, no = min(stats.dt$tSNEx) + (max(stats.dt$tSNEx)-min(stats.dt$tSNEx))/(length(unique(sampleSheet$Donor))-1)*(z-1)), 
             tick = F, col.axis = unique(sampleSheet$DonorCol)[z], labels = unique(sampleSheet$Donor)[z]))
      
      # Color by CM well
      stats.dt[[x]]$wellcol=gsub(".*-","",stats.dt[[x]][["cell"]])
      stats.dt[[x]]$wellcol=gsub("1", sampleSheet$WellCol[1], stats.dt[[x]]$wellcol)
      stats.dt[[x]]$wellcol=gsub("2", sampleSheet$WellCol[2], stats.dt[[x]]$wellcol)
      stats.dt[[x]]$wellcol=gsub("3", sampleSheet$WellCol[3], stats.dt[[x]]$wellcol)
      plotTSNE(cbind(stats.dt$tSNEx, stats.dt$tSNEy), cex = cexsize, col = stats.dt[[x]]$wellcol, main = "Library")
      sapply(1:nrow(sampleSheet), function(z) axis(side = 1, at = min(stats.dt$tSNEx) + (max(stats.dt$tSNEx)-min(stats.dt$tSNEx))/(nrow(sampleSheet)-1)*(z-1),
                                                   tick = F, col.axis = sampleSheet$WellCol[z], labels = z))
      plot(x = rep(1, nrow(sampleSheet)), y = nrow(sampleSheet):1, pch = 21, bg = sampleSheet$WellCol, cex = 3, xlim = c(0,15), axes = F, ylab = "", xlab = "")
      text(x = rep(1.5, nrow(sampleSheet)), y = nrow(sampleSheet):1, labels = sampleSheet$library_id, pos = 4)
      
      dev.off()
      
      # Color by signature score
      if (sig == T) {
        cell_signatures(cbind(stats.dt$tSNEx, stats.dt$tSNEy), CM.norm.log, signatures = signatures, file.path.prefix = paste0(out.dir, "/", samples[x]),  cexsize, type = "tSNE", subset = 30:45)
      }
      
      
      ########
      # UMAP #
      ########
      
      message("\nCalculating UMAP coordinates")
      
      # Make UMAP with default values & add values to stats.dt
      set.seed(42)
      tryCatch(
        expr={
          CM.umap <- umap(as.matrix(t(CM.norm.log.var)))
        },
        error = function(e){
          print("NO")
        }
      )
      
      stopifnot( all(colnames(CM.norm.log) == stats.dt$cell) )# just checking...
      stats.dt$umapx <- CM.umap$layout[,1]
      stats.dt$umapy <- CM.umap$layout[,2]
      
      message(paste("\nGenerating", paste0(out.dir, "/", samples[x], "_UMAP.pdf")))
      pdf(file = paste0(out.dir, "/", samples[x], "_UMAP.pdf"), width = 6, height = 6)
      par(mar=c(4, 4, 4, 4))
      
      # Color by UMI
      plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), pch = 16, cex = cexsize, col = colCustom(stats.dt[[x]]$umis, c(2000, 20000)),
               main = paste0("UMIs (", min(stats.dt[[x]]$umis), "-", max(stats.dt[[x]]$umis), ", mean ", round(mean(stats.dt[[x]]$umis), 0), ")"))
      axis(side = 1, at = 0, label = paste(nrow(stats.dt[[x]]), "cells"), tick = F)
      
      # Color by Donor
      plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), cex = cexsize, col = stats.dt[[x]]$donorcol, main = "Donor")
      sapply(1:length(unique(sampleSheet$Donor)), function(z) axis(side = 1, at = ifelse(length(unique(sampleSheet$Donor)) == 1,
                                                                                         yes = 0, no = min(stats.dt$umapx) + (max(stats.dt$umapx)-min(stats.dt$umapx))/(length(unique(sampleSheet$Donor))-1)*(z-1)), 
                                                                   tick = F, col.axis = unique(sampleSheet$DonorCol)[z], labels = unique(sampleSheet$Donor)[z]))
      
      # Color by CM well
      stats.dt[[x]]$wellcol=gsub(".*-","",stats.dt[[x]][["cell"]])
      stats.dt[[x]]$wellcol=gsub("1", sampleSheet$WellCol[1], stats.dt[[x]]$wellcol)
      stats.dt[[x]]$wellcol=gsub("2", sampleSheet$WellCol[2], stats.dt[[x]]$wellcol)
      stats.dt[[x]]$wellcol=gsub("3", sampleSheet$WellCol[3], stats.dt[[x]]$wellcol)
      plotUMAP(cbind(stats.dt$umapx, stats.dt$umapy), cex = cexsize, col = stats.dt[[x]]$wellcol, main = "Library")
      sapply(1:nrow(sampleSheet), function(z) axis(side = 1, at = min(stats.dt$umapx) + (max(stats.dt$umapx)-min(stats.dt$umapx))/(nrow(sampleSheet)-1)*(z-1),
                                                   tick = F, col.axis = sampleSheet$WellCol[z], labels = z))
      plot(x = rep(1, nrow(sampleSheet)), y = nrow(sampleSheet):1, pch = 21, bg = sampleSheet$WellCol, cex = 3, xlim = c(0,15), axes = F, ylab = "", xlab = "")
      text(x = rep(1.5, nrow(sampleSheet)), y = nrow(sampleSheet):1, labels = sampleSheet$library_id, pos = 4)
      
      dev.off()
      
      # Color by signature score
      if (sig == T) {
        cell_signatures(cbind(stats.dt$umapx, stats.dt$umapy), CM.norm.log, signatures = signatures, file.path.prefix = paste0(out.dir, "/", samples[x]), cexsize, type = "UMAP", subset = 30:45)
      }
     
  } # plotting close
    
  } # main for loop close
  
  rm(list=setdiff(ls(), c("i", "directories", "genome", "suffix", "min.umis", "min.genes", "max.mito", "sig", "dense")))
  
} # function close
