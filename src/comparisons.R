# Andrew Borgman
# 12/21/2012
# Comparisons to consider:
#
# The description should be clear. We have 8 groups, 3 each.
#  
# Wildtype:       Series_1_1_WT, 2_1_WT,  3_1_WT
# Mutant 1393:    Series_1_2_1393, 2_2_1393, 3_2_1393
# Mutant 1394:    Series_1_5_1394, 2_5_1394, 3_5_1394
# Mutant 1436:    Series_1_4_1436, 2_4_1436, 3_4_1436
# Mutant 1513:    Series_1_3_1513, 2_3_1513, 3_3_1513
# Lam control:    Series_1_Lam, 2_Lam, 3_Lam
# T7 control:     Series_1_T7_5, 2_T7, 3_T7
# Pool control:   Pool_1, Pool_2, Pool_3
# 
# What I did with the microarray and what needs to done with the HiSeq data are the following comparisons:
#  
# Wildtype over T7, Lam, and Pool control                      3 comparisons
# Wildtype over 1393, 1394, 1436, and 1513                     4 comparisons
# Mutant 1393 over T7, Lam, and Pool controls                 3 comparisons
# Mutant 1394 over T7, Lam, and Pool controls                 3 comparisons
# Mutant 1436 over T7, Lam, and Pool controls                 3 comparisons
# Mutant 1513 over T7, Lam, and Pool controls                 3 comparisons
#  
# So this is quite a number of comparisons. That will then be the core of the paper.
# 
# I attached also the stats. Sample 3-5_1394 had only few reads. The results look ok though. 
# They could not tell me what actually the reason is for this.
################################################################################################                                           
library(DESeq) 
library(edgeR) 
library(samr)
library(VennDiagram)
library(RColorBrewer)
library(gplots)

# Read in the data making the row names the first column
dat <- read.csv("/home/andrew/Dropbox/bernie_project/hi_seq/clean_counts.csv",row.names=1,check.names=F)

png("library_sizes.png", width=1200, height=800)
barplot(colSums(dat)[order(colSums(dat))], las=2, main="Library Size by Sample", ylab="# of Mapped Reads")
dev.off()

# Gene,Sample_Series_1_5_1394,Sample_3_Lam,Sample_2_4_1436,Sample_Series_1_Lam,Sample_3_T7,Sample_2_5_1394,Sample_3_2_1393,Sample_Pool_2,Sample_lane3-Undetermined_indices,Sample_2_Lam,Sample_3_3_1513,Sample_3_4_1436,Sample_Series_1_T7_5,Sample_lane2-Undetermined_indices,Sample_Series_1_1_WT,Sample_2_3_1513,Sample_2_T7,Sample_2_1_WT,Sample_3_5_1394,Sample_3_1_WT,Sample_Series_1_2_1393,Sample_Series_1_4_1436,Sample_Pool_1,Sample_Series_1_3_1513,Sample_2_2_1393,Sample_Pool_3

# Sample Groups
wildtype <- c("Sample_Series_1_1_WT", "Sample_2_1_WT",  "Sample_3_1_WT")
mutant_1393 <- c("Sample_Series_1_2_1393", "Sample_2_2_1393", "Sample_3_2_1393")
mutant_1394 <- c("Sample_Series_1_5_1394", "Sample_2_5_1394", "Sample_3_5_1394")
mutant_1436 <- c("Sample_Series_1_4_1436", "Sample_2_4_1436", "Sample_3_4_1436")
mutant_1513 <- c("Sample_Series_1_3_1513", "Sample_2_3_1513", "Sample_3_3_1513")
lam <- c("Sample_Series_1_Lam", "Sample_2_Lam", "Sample_3_Lam")
t7 <- c("Sample_Series_1_T7_5", "Sample_2_T7", "Sample_3_T7")
pool <- c("Sample_Pool_1", "Sample_Pool_2", "Sample_Pool_3")

# Store all comparisons here
comparisons <- list()

# Comparison 1
# Wildtype over T7
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% t7)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","T7","T7","T7"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="T7")

comparisons[[1]] <- list(counttable, meta)

# Comparison 2
# Wildtype over Lam
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% lam)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Lam","Lam","Lam"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Lam")

comparisons[[2]] <- list(counttable, meta)

# Comparison 3
# Wildtype over Pool
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% pool)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Pool","Pool","Pool"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Pool")

comparisons[[3]] <- list(counttable, meta)

# Comparison 4
# Wildtype over Mutant 1393
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% mutant_1393)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Mutant_1393","Mutant_1393","Mutant_1393"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Mutant_1393")

comparisons[[4]] <- list(counttable, meta)

# Comparison 5
# Wildtype over Mutant 1394
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% mutant_1394)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Mutant_1394","Mutant_1394","Mutant_1394"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Mutant_1394")

comparisons[[5]] <- list(counttable, meta)

# Comparison 6
# Wildtype over Mutant 1436
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% mutant_1436)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Mutant_1436","Mutant_1436","Mutant_1436"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Mutant_1436")

comparisons[[6]] <- list(counttable, meta)

# Comparison 7
# Wildtype over Mutant 1513
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% wildtype)], dat[,which(names(dat) %in% mutant_1513)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Mutant_1513","Mutant_1513","Mutant_1513"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Mutant_1513")

comparisons[[7]] <- list(counttable, meta)

# Comparison 8
# Mutant 1393 over T7
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1393)], dat[,which(names(dat) %in% t7)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","T7","T7","T7"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="T7")

comparisons[[8]] <- list(counttable, meta)

# Comparison 9
# Mutant 1393 over Lam
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1393)], dat[,which(names(dat) %in% lam)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Lam","Lam","Lam"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Lam")

comparisons[[9]] <- list(counttable, meta)

# Comparison 10
# Mutant 1393 over Pool
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1393)], dat[,which(names(dat) %in% pool)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Pool","Pool","Pool"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Pool")

comparisons[[10]] <- list(counttable, meta)

# Comparison 11
# Mutant 1394 over T7
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1394)], dat[,which(names(dat) %in% t7)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","T7","T7","T7"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="T7")

comparisons[[11]] <- list(counttable, meta)

# Comparison 12
# Mutant 1394 over Lam
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1394)], dat[,which(names(dat) %in% lam)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Lam","Lam","Lam"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Lam")

comparisons[[12]] <- list(counttable, meta)

# Comparison 13
# Mutant 1394 over Pool
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1394)], dat[,which(names(dat) %in% pool)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Pool","Pool","Pool"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Pool")

comparisons[[13]] <- list(counttable, meta)

# Comparison 14
# Mutant 1436 over T7
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1436)], dat[,which(names(dat) %in% t7)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","T7","T7","T7"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="T7")

comparisons[[14]] <- list(counttable, meta)

# Comparison 15
# Mutant 1436 over Lam
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1436)], dat[,which(names(dat) %in% lam)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Lam","Lam","Lam"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Lam")

comparisons[[15]] <- list(counttable, meta)

# Comparison 16
# Mutant 1436 over Pool
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1436)], dat[,which(names(dat) %in% pool)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Pool","Pool","Pool"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Pool")

comparisons[[16]] <- list(counttable, meta)

# Comparison 17
# Mutant 1513 over T7
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1513)], dat[,which(names(dat) %in% t7)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","T7","T7","T7"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="T7")

comparisons[[17]] <- list(counttable, meta)

# Comparison 18
# Mutant 1513 over Lam
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1513)], dat[,which(names(dat) %in% lam)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Lam","Lam","Lam"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Lam")

comparisons[[18]] <- list(counttable, meta)

# Comparison 19
# Mutant 1513 over Pool
# Parse out cases and controls
counttable <- data.frame(cbind(dat[,which(names(dat) %in% mutant_1513)], dat[,which(names(dat) %in% pool)]))

## Make metadata data.frame -- identify condition label and  
meta <- data.frame(
  row.names=colnames(counttable),
  condition=c("treated", "treated", "treated","Pool","Pool","Pool"),
  libType=c("single", "single", "single", "single", "single", "single"))
meta$condition <- relevel(meta$condition, ref="Pool")

comparisons[[19]] <- list(counttable, meta)
