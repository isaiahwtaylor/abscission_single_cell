edgeR_2_sample <- function(reads_df, s1_label, s2_label, s1_indices, s2_indices, annotations) {
  
  
  #make DGElist
  x3 <- DGEList(counts = reads_df, genes = rownames(reads_df))
  
  
  #reads per library uniquely mapped to a gene:
  
  #make cpm and lcpm
  cpm <- cpm(x3)
  lcpm <- cpm(x3, log=TRUE)
  
  #keep only genes that are expressed. "Expressed" here means counts observed in at least 3 samples
  dim(x3)
  keep.exprs <- rowSums(cpm>0)>=3
  x3 <- x3[keep.exprs,, keep.lib.sizes=FALSE]
  dim(x3) #compare to dim(x) above
  
  #normalize data after removing low expressed genes
  x3 <- calcNormFactors(x3)
  
  #cpm, lcpm of normalized values
  cpm <- cpm(x3)
  lcpm <- cpm(x3, log=TRUE)
  
  
  #account for factors in experiment
  treatment=as.factor(c(rep("condition_1",length(s1_indices)), rep("condition_2",length(s2_indices))))
  
  #assign factors to DGElist
  x3$samples$treatment=treatment
  
  #make design matrix of experimental factors
  design <- model.matrix(~0+treatment)
  colnames(design) <- levels(treatment)
  design #check matrix
  
  #fit normalized dge list to model matrix
  x3 <- estimateDisp(x3, design)
  fit <- glmQLFit(x3, design)

  #tr <- glmTreat(fit, coef=2, lfc=1)
  #topTags(tr)
  
  #making contrast matrix for tests of interest
  my.contrasts <- makeContrasts(s1_v_s2=condition_1-condition_2, levels=design)
  
  #QLF tests 
  s1_v_s2 <- glmQLFTest(fit, contrast=my.contrasts[,"s1_v_s2"]) 
 
   
  #extract results. toptags adds FDR column.
  sec2 <- data.frame(topTags(s1_v_s2, n = "Inf")$table)
  #sec2=sec2[c(1,7,2,3,4,5,6)]
  sec2[[s1_label]] <- apply(cpm[sec2[["genes"]],s1_indices], 1,mean)
  sec2[[s2_label]] <- apply(cpm[sec2[["genes"]],s2_indices], 1,mean)
  
  #return(sec2)
  sec2 = cbind(sec2, cpm[sec2$genes,c(s1_indices, s2_indices)])
  #return(sec2)
  sec2=merge(sec2, annotations, by.x=c("genes"), by.y=c("gene_id"), all.x = T, all.y = F)
  sec2= sec2[!duplicated(sec2$genes), ]
  #rownames(sec2)=sec2$genes
  return(sec2)

}


edgeR_2_sample_other_factors <- function(reads_df, s1_label, s2_label, s1_indices, s2_indices, annotations) {
  
  
  #make DGElist
  x3 <- DGEList(counts = reads_df, genes = rownames(reads_df))
  
  
  #reads per library uniquely mapped to a gene:
  
  #make cpm and lcpm
  cpm <- cpm(x3)
  lcpm <- cpm(x3, log=TRUE)
  
  #keep only genes that are expressed. "Expressed" here means counts observed in at least 3 samples
  dim(x3)
  keep.exprs <- rowSums(cpm>0)>=3
  x3 <- x3[keep.exprs,, keep.lib.sizes=FALSE]
  dim(x3) #compare to dim(x) above
  
  #normalize data after removing low expressed genes
  x3 <- calcNormFactors(x3)
  
  #cpm, lcpm of normalized values
  cpm <- cpm(x3)
  lcpm <- cpm(x3, log=TRUE)
  
  
  #account for factors in experiment
  phenotype=as.factor(c("wt", "wt", "wt", "wt", "mut", "mut", "mut", "mut"))
  method=as.factor(c("total", "total", "sort", "sort", "total", "total", "sort", "sort"))
  #insertion=as.factor(c("0","0","1","1","0","0","2","2"))
  
  #assign factors to DGElist
  x3$samples$phenotype=phenotype
  x3$samples$method=method
  #x3$samples$insertion=insertion
  
  #make design matrix of experimental factors
  design <- model.matrix(~0+phenotype+method)#+insertion)
  #colnames(design) <- c(levels(phenotype), levels(method))#, levels(insertion))
  #design #check matrix
  
  #fit normalized dge list to model matrix
  x3 <- estimateDisp(x3, design)
  fit <- glmQLFit(x3, design)

  #tr <- glmTreat(fit, coef=2, lfc=1)
  #topTags(tr)
  
  #making contrast matrix for tests of interest
  my.contrasts <- makeContrasts(s1_v_s2=phenotypewt-phenotypemut, levels=design)
  
  #QLF tests 
  s1_v_s2 <- glmQLFTest(fit, contrast=my.contrasts[,"s1_v_s2"]) 
 
   
  #extract results. toptags adds FDR column.
  sec2 <- data.frame(topTags(s1_v_s2, n = "Inf")$table)
  #sec2=sec2[c(1,7,2,3,4,5,6)]
  sec2[[s1_label]] <- apply(cpm[sec2[["genes"]],s1_indices], 1,mean)
  sec2[[s2_label]] <- apply(cpm[sec2[["genes"]],s2_indices], 1,mean)
  
  #return(sec2)
  sec2 = cbind(sec2, cpm[sec2$genes,c(s1_indices, s2_indices)])
  #return(sec2)
  sec2=merge(sec2, annotations, by.x=c("genes"), by.y=c("gene_id"), all.x = T, all.y = F)
  sec2= sec2[!duplicated(sec2$genes), ]
  #rownames(sec2)=sec2$genes
  return(sec2)

}