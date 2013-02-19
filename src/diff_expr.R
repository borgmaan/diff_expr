#!/usr/bin/env Rscript
# Andrew Borgman
# 12/1/2012
# EdgeR, DEseq, and SAMseq analysis 
# Sources in file formatted like comparisons.R to 
# specify groups of individuals to be compared. Runs all
# differential expression analyses found in diff_expr_functions.R
################################################################################################                                           
library(DESeq) 
library(edgeR) 
library(samr)
library(VennDiagram)
library(RColorBrewer)
library(gplots)

# Donde Esta?
setwd("/home/andrew/Dropbox/bernie_project/hi_seq/")

# Read in data and comparisons from data script
source("comparisons.R")

# Read in differential expression functions
source("diff_expr_functions.R")

# Move to analysis folder
setwd("/home/andrew/Dropbox/dinesh/analysis/")

# Running all comparisons
for (i in 1:length(comparisons)){
  comp_dir <- paste("./",i,"_filt",sep="")
  print(i)
  print(comp_dir)
  system(paste("rm -r",comp_dir,sep=" "))
  system(paste("mkdir",comp_dir,sep=" "))
  setwd(comp_dir)
  print("made ")
  # Writing out trimmed table with counts
  #write.csv(comparisons[[i]][[1]],"comparison_counts.csv")
  
  # Filtering out genes with no data mapping to them
  genewise_sums <- apply(comparisons[[i]][[1]],1,sum)
  print('f')
  sel <- which(genewise_sums > 6)
  trimmed <- comparisons[[i]][[1]][sel,]
  
  # Running different analysis
  
  # Trys and fails on GLM analysis if there is no
  # extra experimental factor to consider. [think about cleaning this up]
  try(DeseqAnalysis(trimmed, comparisons[[i]][[2]],i),silent=T)
  try(EdgerAnalysis(trimmed, comparisons[[i]][[2]],i),silent=T)
  SamseqAnalysis(trimmed, comparisons[[i]][[2]],i)
  DeseqSimple(trimmed, comparisons[[i]][[2]],i)
  EdgerSimple(trimmed, comparisons[[i]][[2]],i)
  
  # Memory management
  rm(trimmed)
  comparisons[[i]][[1]] <- c()
  
  setwd("../")  
}
