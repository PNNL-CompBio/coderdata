#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# If BiocManager not installed, install it
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

library(GenomicRanges)
library(Homo.sapiens)


# function to return gene mappings for each coordinate range (row) in the segfile
splitColumnByOverlap <-
  function(query, subject, column="ENTREZID", ...)
  {
    olaps <- findOverlaps(query, subject, ...)
    f1 <- factor(subjectHits(olaps),
                 levels=seq_len(subjectLength(olaps)))
    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
  }

segfile = read.csv(args[1])

# create genomic ranges object from segfile for use with findOverlaps()
gr <- GRanges(seqnames = Rle(paste0('chr', segfile$chrom)), 
              ranges = IRanges(segfile$loc.start, end = segfile$loc.end), 
              strand = Rle(c("-", "0", "+")[sign(segfile$loc.start) +2]), 
              score = segfile$seg.mean, 
              ID = segfile$ID)

# create genes GRanges obj from database
genes<- genes(Homo.sapiens, columns =c("ENTREZID"))

geneHitsByRow <- splitColumnByOverlap(genes, gr, "ENTREZID")

# create matrices of annotations and scores/patient IDs in segfile
# with a catch if there are no gene hits found
matrixresults <- list()
noresults <- c()
for (i in 1:length(geneHitsByRow)){
  if(length(geneHitsByRow[[i]])>0){
    
    matrixresults[[i]] <- data.frame(ENTREZID = as.matrix(geneHitsByRow[[i]]), rep(mcols(gr[i,]), length(geneHitsByRow[[i]])))
  }else{
    noresults <- c(noresults, i)
  matrixresults[[i]] <- data.frame(ENTREZID = NA, mcols(gr)[i,])
}
}
# concatenate annotated matrices
allCNV <- do.call(rbind, matrixresults)

# drop NAs
completeAllCNV <- allCNV[complete.cases(allCNV),]
# aggregate scores for the same genes (that came from different regions) for the same patient ID 
aggregatedAllCNV <- aggregate(completeAllCNV$score, by = list(ENTREZID =completeAllCNV$ENTREZID, ID = completeAllCNV$ID), FUN=mean)
names(aggregatedAllCNV)[names(aggregatedAllCNV)=='x'] <- "score"
# write results to csv for further processing in Python
write.csv(aggregatedAllCNV, args[2], row.names=F, quote =F)