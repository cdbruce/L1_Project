#!/bin/env Rscript
require(IWTomics)
require(dendextend)
#install.packages("dendextend")

setwd("~/Desktop/cleanControl_pvalues/")


#################### test ####################
#### selecting largest pvalues from 5 random results ####
##############################################
#load('L1_autosomes_results_smoothed_mean_1.RData')
#random1_pval<-as.data.frame(adjusted_pval(result_mean))

for (r in 1:5){
  file=paste('L1_autosomes_results_smoothed_mean_', r, '.RData', sep="")
  load(file)
  out <- as.data.frame(adjusted_pval(result_mean))
  assign(paste('pval_random_',r,sep=""),out)
}






pdf('IWT_autosomes_smoothed_mean.pdf',width=7,height=10)
plotTest(result_mean,col=c('red','blue','green','black'),
         scale_threshold=unlist(lapply(result_mean@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
dev.off()
plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',scale_threshold=10,
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),"_scale10.pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),".pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',scale_threshold=10,
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),"_scale10.pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
