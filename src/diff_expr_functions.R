#!/usr/bin/env Rscript
# Andrew Borgman
# 12/7/2012
# EdgeR, DEseq, and SAMseq Functions for differential expression
# All functions are fed identical input parameters

# DEseq analysis function  
DeseqAnalysis <- function(count_data_frame, meta_data_frame,comparison_number){
  # Creating DEseq base object
  d <- newCountDataSet(count_data_frame, meta_data_frame)
  
  # Quality Control Steps
  # The idea of independent filtering is to filter out those tests from the procedure that 
  # have no, or little chance of showing significant evidence, without even looking at their
  # test statistic. Typically, this results in increased detection power at the same experiment-wide type I
  # error. Here, we measure experiment-wide type I error in terms of the false discovery rate.
  gene_sums <- rowSums(counts(d))
    
  # Normalization of counts
  d <- estimateSizeFactors(d)
  
  # Estimate dispersions
  d <- estimateDispersions(d ,fitType="local")
  
  dispersions <- data.frame(cbind(
    fData(estimateDispersions(d ,fitType="local"))[,1],
    fData(estimateDispersions(d ,fitType="local", sharingMode="fit-only"))[,1],
    fData(estimateDispersions(d ,fitType="local", sharingMode="maximum"))[,1],
    fData(estimateDispersions(d ,fitType="local", sharingMode="gene-est-only"))[,1],
    fData(estimateDispersions(d ,fitType="local",method="blind"))[,1],
    fData(estimateDispersions(d ,fitType="local",method="blind", sharingMode="fit-only"))[,1],
    fData(estimateDispersions(d ,fitType="local",method="blind", sharingMode="maximum"))[,1],
    fData(estimateDispersions(d ,fitType="local",method="blind", sharingMode="gene-est-only"))[,1]
  ))
  
  png(paste("DEseq_dispersion_plot_",comparison_number,".png",sep=""),width=1200,height=600)
  boxplot(dispersions,main="Dispersions for Assorted Estimates", col="grey", names=c("local", "local;fit-only", "local;maximum",
                                                                         "local;gene-est-only", "blind", "blind;fit-only", 
                                                                         "blind;maximum","blind;gene-est-only"), las=3)
  dev.off()
    
  # Plotting dispersion estimates for each
  png(paste("DEseq_dispersion_plot_",comparison_number,".png",sep=""),width=600,height=600)
  plotDispEsts(d, main="DESeq: Per-gene dispersion\nestimates", col="#00000040")
  dev.off()
    
  ## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
  png(paste("DEseq_pca_cluster_",comparison_number,".png",sep=""),width=600,height=600)
  try(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition", "libType")),silent=T)
  dev.off()
    
  # Multi-Factor Analysis
  # For inference, we now specify two models by formulas. The full model regresses the genes’ expression on both the
  # library type and the treatment condition, the reduced model regresses them only on the library type. For each gene,
  # we fit generalized linear models (GLMs) according to the two models, and then compare them in order to infer whether
  # the additional specification of the treatment improves the fit and hence, whether the treatment has significant effect.
  
  # Fitting model with libtype and condition
  # The first three columns show the fitted coefficients, converted to a logarithm base 2 scale. The log2 fold change
  # due to the condition is shown in the third column.
  if(length(unique(meta_data_frame$libType)) == 1){ 
    d_fit1 <- fitNbinomGLMs(d, count ~  condition)
    d_fit0 <- fitNbinomGLMs(d, count ~ 1)
  }else {  
    d_fit1 <- fitNbinomGLMs(d, count ~ libType + condition)
    d_fit0 <- fitNbinomGLMs(d, count ~ libType)
  }
  
  # Check if inclusion of condition significantly improve fit
  dpval <- nbinomGLMTest(d_fit1, d_fit0)
  
  # P-value adjustment 
  dpadj <- p.adjust(dpval, method="BH")
  
  ## Writing out top hits  
  # Adding in p-values 
  d_fit1$p_value <- dpval
  d_fit1$p_adjust <- dpadj
  
  # Adding counts to output
  count_data_frame$n_reads <- apply(count_data_frame, 1, sum)
  d_fit1 <- merge(d_fit1,count_data_frame,by='row.names')
  d_fit1 <- d_fit1[order(d_fit1$p_adjust),]
  d_fit1 <- d_fit1[!is.na(d_fit1[,3]),]
  
  write.csv(d_fit1,paste("DEseq_sig_genes_",comparison_number,".csv",sep=""))
  
  # Writing out normalized counts for good measure
  write.csv(counts( d, normalized=TRUE ), paste("DEseq_normalized_counts_",comparison_number,".csv",sep=""))
  
  # Heatmap of read counts
  d_full_blind <- estimateDispersions(d, fitType="local", method = "blind")
  vsdFull <- varianceStabilizingTransformation(d_full_blind)
  
  # Plotting heatmap for top genes -- work on getting these for genes with top P-Values
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  png(paste("DEseq_all_genes_",comparison_number,".png",sep=""),width=600,height=2600)
  heatmap.2(exprs(vsdFull), Colv=F, col = hmcol, trace="none", margin=c(10, 6))
  dev.off()
  
  # Subsetting top genes
  if(length(meta_data_frame[,1]) == 6) {
    avg_case <- apply(exprs(vsdFull)[ ,c(1,2,3)], 1, mean)
    avg_con <- apply(exprs(vsdFull)[ ,c(4,5,6)], 1, mean)
    avg_diff <- avg_case - avg_con
    
    top_genes <- names(sort(avg_diff, decreasing=T))[1:100]
          
    png(paste("DEseq_enriched_genes_",comparison_number,".png",sep=""),width=600,height=1500)
    heatmap.2(exprs(vsdFull)[which(row.names(exprs(vsdFull)) %in% top_genes), ], Colv=F, col = hmcol, trace="none", margin=c(10, 6))
    dev.off()  
  }
}

DeseqSimple <- function(count_data_frame, meta_data_frame, comparison_number){
  # Creating DEseq base object
  d <- newCountDataSet(count_data_frame, meta_data_frame$condition)
  
  # Quality Control Steps
  # The idea of independent filtering is to filter out those tests from the procedure that 
  # have no, or little chance of showing significant evidence, without even looking at their
  # test statistic. Typically, this results in increased detection power at the same experiment-wide type I
  # error. Here, we measure experiment-wide type I error in terms of the false discovery rate.
  gene_sums <- rowSums(counts(d))
  
  # Normalization of counts
  d <- estimateSizeFactors(d)
  
  # Estimate dispersions
  d <- estimateDispersions(d, fitType="local")
  
  # Plotting dispersion estimates for each
  png(paste("DEseq_dispersion_plot_",comparison_number,".png",sep=""),width=600,height=600)
  plotDispEsts(d, main="DESeq: Per-gene dispersion\nestimates",col="#00000040")
  dev.off()
    
  ## Principal components biplot on variance stabilized data, color-coded by condition-librarytype
  png(paste("DEseq_simple_pca_cluster_",comparison_number,".png",sep=""),width=600,height=600)
  try(plotPCA(varianceStabilizingTransformation(d), intgroup=c("condition")),silent=T)
  dev.off()
    
  # Significance testing
  res <- nbinomTest(d, levels(meta_data_frame$condition)[1], levels(meta_data_frame$condition)[2])
  
  ## Writing out top hits  
  write.csv(res,paste("DEseq_simple_sig_genes_",comparison_number,".csv",sep=""))
  
  png(paste("DEseq_simple_pvalue",comparison_number,".png",sep=""),width=600,height=600)  
  hist(res$pval, breaks=100, col="skyblue", border="slateblue",main="Distribution of P-Values")
  dev.off()
  
  # Heatmap of read counts
  d_full_blind <- estimateDispersions(d, fitType="local", method = "blind")
  vsdFull <- varianceStabilizingTransformation(d_full_blind)
  
  # Plotting heatmap for top genes -- work on getting these for genes with top P-Values
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  png(paste("DEseq_simple_all_genes_",comparison_number,".png",sep=""),width=600,height=2600)
  heatmap.2(exprs(vsdFull), Colv=F, col = hmcol, trace="none", margin=c(10, 6))
  dev.off()
  
  # Subsetting top genes
  if(length(meta_data_frame[,1]) == 6) {
    avg_case <- apply(exprs(vsdFull)[ ,c(1,2,3)], 1, mean)
    avg_con <- apply(exprs(vsdFull)[ ,c(4,5,6)], 1, mean)
    avg_diff <- avg_case - avg_con
    
    top_genes <- names(sort(avg_diff, decreasing=T))[1:100]
    
    png(paste("DEseq_simple_enriched_genes_",comparison_number,".png",sep=""),width=600,height=1500)
    heatmap.2(exprs(vsdFull)[which(row.names(exprs(vsdFull)) %in% top_genes), ], Colv=F, col = hmcol, trace="none", margin=c(10, 6))
    dev.off()  
  }
  
}


# edgeR -------------------------------------------------------------------
EdgerAnalysis <- function(count_data_frame, meta_data_frame,comparison_number){
  
  ## Make design matrix  
  # Grabbing appropriate reference factor level
  if ("T7" %in% meta_data_frame$condition){
    condition <- relevel(factor(meta_data_frame$condition), ref="T7")
  } 
  if ("Lam" %in% meta_data_frame$condition) {
    condition <- relevel(factor(meta_data_frame$condition), ref="Lam")
  }  
  if ("pool" %in% meta_data_frame$condition) {
    condition <- relevel(factor(meta_data_frame$condition), ref="pool")
  }  
  
  if(length(unique(meta_data_frame$libType)) == 1) {
    edesign <- model.matrix(~condition)
  } else {    
    libType <- factor(meta_data_frame$libType)
    edesign <- model.matrix(~libType+condition)
  }
    
  ## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
  e <- DGEList(counts=count_data_frame)
  e <- calcNormFactors(e)
  e <- estimateGLMCommonDisp(e, edesign)
  e <- estimateGLMTrendedDisp(e, edesign) 
  e <- estimateGLMTagwiseDisp(e, edesign)
  
  ## MDS Plot
  png(paste("edgeR_mds_cluster_",comparison_number,".png",sep=""),width=600,height=600)
  plotMDS(e, main="edgeR MDS Plot")
  dev.off()
  
  ## Biological coefficient of variation plot
  png(paste("edgeR_bcv",comparison_number,".png",sep=""),width=600,height=600)
  plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance",col="#00000040")
  dev.off() 
  
  ## Fit the model, testing the coefficient for the treated vs untreated comparison
  efit <- glmFit(e, edesign)
  efit <- glmLRT(efit, coef="conditiontreated")
  
  ## Make a table of results
  etable <- topTags(efit, n=nrow(e))$table
  etable <- etable[order(etable$FDR), ]
  count_data_frame$n_reads <- apply(count_data_frame, 1, sum)
  mm <- merge(etable,count_data_frame,by='row.names')
  mm <- mm[order(mm$PValue),]
  write.csv(etable, paste("edgeR_sig_genes",comparison_number,".csv",sep=""))
  
  ## MA Plot
  png(paste("edgeR_top_genes_",comparison_number,".png",sep=""),width=800,height=600)
  with(etable, plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
  with(subset(etable, FDR<0.05), points(logCPM, logFC, pch=20, col="red"))
  abline(h=c(-1,1), col="blue")
  dev.off()   
  
  ## Heatmap
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  y <- predFC(e, prior.count.total=2*ncol(e))
  png(paste("edgeR_GLM_heatmap_all_genes_",comparison_number,".png",sep=""),width=600,height=2600)
  heatmap.2(y, Colv=F, col = hmcol, trace="none", margin=c(10, 6))
  dev.off()
  
  # Subsetting top genes
  if(length(meta_data_frame[,1]) == 6) {
    avg_case <- apply(y[ ,c(1,2,3)], 1, mean)
    avg_con <- apply(y[ ,c(4,5,6)], 1, mean)
    avg_diff <- avg_case - avg_con
    
    top_genes <- names(sort(avg_diff, decreasing=T))[1:100]
    
    png(paste("edgeR_GLM_heatmap_enriched_genes_",comparison_number,".png",sep=""),width=600,height=1500)
    heatmap.2(y[which(row.names(y) %in% top_genes), ], Colv=F, col = hmcol, trace="none", margin=c(10, 6))
    dev.off()  
  }     
}

EdgerSimple <- function(count_data_frame, meta_data_frame,comparison_number){
  
  ## Make design matrix  
  # Grabbing appropriate reference factor level
  if ("T7" %in% meta_data_frame$condition){
    condition <- relevel(factor(meta_data_frame$condition), ref="T7")
  } 
  if ("Lam" %in% meta_data_frame$condition) {
    condition <- relevel(factor(meta_data_frame$condition), ref="Lam")
  }  
  if ("pool" %in% meta_data_frame$condition) {
    condition <- relevel(factor(meta_data_frame$condition), ref="pool")
  }  
  
  ## Make new DGEList, normalize by library size, and estimate dispersion allowing possible trend with average count size
  e <- DGEList(counts=count_data_frame, group=meta_data_frame$condition)
  e <- calcNormFactors(e)
  e <- estimateCommonDisp(e)
  e <- estimateTagwiseDisp(e) 
  
  ## Biological coefficient of variation plot
  png(paste("edgeR_bcv_simple_",comparison_number,".png",sep=""),width=600,height=600)
  plotBCV(e, cex=0.4, main="edgeR: Biological coefficient of variation (BCV) vs abundance",col="#00000040")
  dev.off() 
  
  ## Fit the model, testing the coefficient for the treated vs untreated comparison
  et <- exactTest(e)
  etable <- topTags(et,length(et$table[,1]))
  
  ## Make a table of results
  count_data_frame$n_reads <- apply(count_data_frame, 1, sum)
  mm <- merge(etable,count_data_frame,by='row.names')
  mm <- mm[order(mm$PValue),]
  write.csv(etable, paste("edgeR_simple_sig_genes",comparison_number,".csv",sep=""))
  
  ## MA Plot
  png(paste("edgeR_simple_top_genes_",comparison_number,".png",sep=""),width=800,height=600)
  with(mm, plot(logCPM, logFC, pch=20, main="edgeR: Fold change vs abundance"))
  with(subset(mm, FDR<0.05), points(logCPM, logFC, pch=20, col="red"))
  abline(h=c(-1,1), col="blue")
  dev.off()   
  
  ## Heatmap
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  y <- predFC(e, prior.count.total=2*ncol(e))
  png(paste("edgeR_simple_heatmap_all_genes_",comparison_number,".png",sep=""),width=600,height=2600)
  heatmap.2(y, Colv=F, col = hmcol, trace="none", margin=c(10, 6))
  dev.off()
  
  # Subsetting top genes
  if(length(meta_data_frame[,1]) == 6) {
    avg_case <- apply(y[ ,c(1,2,3)], 1, mean)
    avg_con <- apply(y[ ,c(4,5,6)], 1, mean)
    avg_diff <- avg_case - avg_con
    
    top_genes <- names(sort(avg_diff, decreasing=T))[1:100]
    
    png(paste("edgeR_simple_heatmap_enriched_genes_",comparison_number,".png",sep=""),width=600,height=1500)
    heatmap.2(y[which(row.names(y) %in% top_genes), ], Colv=F, col = hmcol, trace="none", margin=c(10, 6))
    dev.off()  
  }
  
  
}



#########################################################################################
## Checking results with non-parametric SAMseq 
# This method might be preferred here given Bernie's weird data
# generated by his pulldown method.

SamseqAnalysis <- function(count_data_frame, meta_data_frame,comparison_number){
  
  # Assessing significance with non-parametric model
  y <- as.numeric(meta_data_frame$condition)
  samfit <- SAMseq(count_data_frame, y, resp.type = "Two class unpaired",  geneid = row.names(count_data_frame), 
                   genenames = row.names(count_data_frame), nperms = 10000, fdr.output = 1)
  
  # Sort, merge, write out top hits ## TODO: See if thie ever spits out down regulated genes...
  count_data_frame$n_reads <- apply(count_data_frame, 1, sum)
  up_sig_genes <- data.frame(samfit$siggenes.table$genes)
  row.names(up_sig_genes) <- up_sig_genes$Gene.ID
  sig_table <- merge(count_data_frame, up_sig_genes, by='row.names')
  #sig_table <- sig_table[order(sig_table$q.value...),]
  tmp <- data.frame(cbind(as.character(sig_table$Gene.Name),as.numeric(as.character(sig_table$q.value...))))
  names(tmp) <- c("gene","fdr_q_value")
  
  # Building new, better samseq table
  d <- samfit$samr.obj$x
  d$wilcoxon_stat <- samfit$samr.obj$tt
  d$fold_change <- samfit$samr.obj$foldchange
  d$gene <- row.names(d)
  
  mm <- merge(d,tmp,by.x = length(d),by.y=1,all.x=T)
  mm <- mm[order(mm$fold_change,decreasing=T),]
  write.csv(mm, paste("SAMseq_sig_genes",comparison_number,".csv",sep=""))
  
  # For heatmap
  #mapper <- sig_table
  #row.names(mapper) <- mapper[,1]
  #mapper <- mapper[,c(2:7)]
  #mapper <- newCountDataSet(mapper,c(1,1,1,2,2,2))
  #mapper <- estimateSizeFactors(mapper)
  #mapper <- estimateDispersions(mapper, fitType="local", method = "blind")
  #vsdFull <- varianceStabilizingTransformation(mapper)
  #hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  #png(paste("SAMseq_top_genes_",comparison_number,".png",sep=""),width=600,height=2600)
  #heatmap.2(exprs(vsdFull), Colv=F, col = hmcol, trace="none", margin=c(10, 6))
  #dev.off()
}
