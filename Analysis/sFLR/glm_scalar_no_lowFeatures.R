## Logistic regression on L1 data: scalar predictors 
setwd('~/Desktop/regression/10_randoms/')
# Load the transformed data
# source("https://bioconductor.org/biocLite.R")
# biocLite("IWTomics")
require(IWTomics)
library(fda)



# Select functional predictors by the localization table created based on the bosplots & IWT results
localization_all<-read.table("localization_table.txt",header=TRUE)
localization_all$test1
func_1<-as.vector(localization_all$test1)
func_2<-as.vector(localization_all$test2)
func_3<-as.vector(localization_all$test3)
func_4<-as.vector(localization_all$test4)
func_5<-as.vector(localization_all$test5)
func_6<-as.vector(localization_all$test6)

##Previous version
# y=1
# n=0
# func_1<-as.vector(c(y,y,n,y,n,y,y,n,y,n,y,n,y,y,y,n,n,n,y,n,n,y,n,y,y,y,n,n,n,n,y,n,y,y,n,y,y,y,n,y,n,n,n,n,y))
# func_2<-as.vector(c(y,n,n,y,y,n,n,n,n,y,n,n,n,n,y,n,n,y,y,y,n,n,n,n,n,n,n,n,n,n,n,y,n,n,n,y,n,n,n,n,n,n,n,n,n))
# func_3<-as.vector(c(y,y,y,y,y,n,n,n,y,y,n,n,y,y,n,n,y,n,y,n,y,y,n,n,n,y,n,n,n,n,y,y,y,n,n,y,y,y,n,n,n,n,y,n,n))
# func_4<-as.vector(c(y,y,n,n,y,y,y,n,y,y,y,y,y,y,y,y,n,n,n,y,y,y,n,y,y,y,n,n,y,n,y,n,n,y,y,y,y,n,n,n,n,n,n,n,y))
# func_5<-as.vector(c(y,y,y,y,y,y,y,y,y,y,y,y,y,y,y,y,y,n,y,n,y,y,n,y,y,y,n,n,y,y,y,y,y,y,y,y,y,y,n,n,n,n,n,n,y))
# func_6<-as.vector(c(y,y,y,y,y,n,n,n,y,y,n,n,y,y,n,n,y,n,y,n,y,n,n,n,n,y,n,n,y,y,n,n,n,n,n,n,y,y,n,n,y,n,n,n,n))
# 

####Here for comparisons 4,5 and 6, we also changed the comparisons' directions to L1pol vs denovo, L1hs vs denovo, L1pol vs L1hs
for (r in 1:10){
  #r=1
  setwd(paste("~/Desktop/regression/10_randoms/random",r,sep=''))
  load(file = paste("L1_transformed_random_",r,'.RData',sep=''))
  # test<-result_mean@features$H2AFZ_signal[[1]]
  # dim(test)
  
  #Create a data table with feature names
  # feature_names<-as.vector(c('H2AFZ_signal','H3K27ac_signal','H4K20me1_signal','H3K36me3_signal','H3K4me1_signal','H3K4me2_signal',
  #                            'H3K4me3_signal','H3K79me2_signal','H3K9ac_signal','H3K9me3_signal','H3K27me3_signal','CTCF_signal',
  #                            'DNase_DHS_signal','RNA_PolII','Quadruplex','A_Phased','Direct_Repeats','Inverted_Repeats','Mirror_Repeats',
  #                            'Z_DNA','Most_Conserved','Exons','Introns','GC_Content','AT_Content','Mononucleotide','Morethan1_nuc',
  #                            'DNA_Transposons','SINE_Alu','SINE_MIR','LTR_Elements','L1_Targets','LINE_L2&L3','CpG_Islands','5hMc',
  #                            'Sperm_hypometh','Rep_origin','Recombination_Hot','Exon_Expression','Gene_Expression','Tracript_Expression',
  #                            'Testis_Expression','chh_meth','chg_meth','cpg_meth'))
  # 
  feature_names<-as.vector(result_mean@metadata$feature_datasets$name)
  feature_num<- as.vector(c(1:45))
  feature<-as.data.frame(matrix(ncol=1))
  for (i in 1:45){
    feature<-rbind(feature,paste("feature_", i, sep = ""))
  }
  feature<-feature[-1,]
  feature_table<-data.frame(feature_num,feature,feature_names)
  
  
  #######Logsitic regression for scalar predictors
  # Comparison 1 (denovo vs control): Logistic regression for scalar predictors
  comp1_scalar<-which(func_1=="n")
  comp1_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp1_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp1_scalar){
    #i=comp1_scalar[1]
    l1<-result_mean@features[[i]][[1]]
    ct<-result_mean@features[[i]][[4]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("feature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp1_scalar_logit_sum<-rbind(comp1_scalar_logit_sum,fit_sum)
  }
  comp1_scalar_logit_sum<-comp1_scalar_logit_sum[-1,]
  rownames(comp1_scalar_logit_sum)<-c(1:length(comp1_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp1_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp1_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp1_scalar_logit_sum, file='comp1_scalar_logit_sum.csv')
  
  
  ##########
  # Comparison 2 (l1pol vs control): Logistic regression for scalar predictors
  comp2_scalar<-which(func_2=="n")
  comp2_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp2_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp2_scalar){
    #i=comp1_scalar[1]
    l1<-result_mean@features[[i]][[2]]
    ct<-result_mean@features[[i]][[4]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("feature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp2_scalar_logit_sum<-rbind(comp2_scalar_logit_sum,fit_sum)
  }
  comp2_scalar_logit_sum<-comp2_scalar_logit_sum[-1,]
  rownames(comp2_scalar_logit_sum)<-c(1:length(comp2_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp2_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp2_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp2_scalar_logit_sum, file='comp2_scalar_logit_sum.csv')
  
  
  ##########
  # Comparison 3 (l1HS vs control): Logistic regression for scalar predictors
  comp3_scalar<-which(func_3=="n")
  comp3_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp3_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp3_scalar){
    #i=comp1_scalar[1]
    l1<-result_mean@features[[i]][[3]]
    ct<-result_mean@features[[i]][[4]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("feature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp3_scalar_logit_sum<-rbind(comp3_scalar_logit_sum,fit_sum)
  }
  comp3_scalar_logit_sum<-comp3_scalar_logit_sum[-1,]
  rownames(comp3_scalar_logit_sum)<-c(1:length(comp3_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp3_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp3_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp3_scalar_logit_sum, file='comp3_scalar_logit_sum.csv')
  
  
  ##########
  # Comparison 4 (l1pol vs denovo): Logistic regression for scalar predictors
  comp4_scalar<-which(func_4=="n")
  comp4_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp4_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp4_scalar){
    #i=comp1_scalar[1]
    l1<-result_mean@features[[i]][[2]]
    ct<-result_mean@features[[i]][[1]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("feature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp4_scalar_logit_sum<-rbind(comp4_scalar_logit_sum,fit_sum)
  }
  comp4_scalar_logit_sum<-comp4_scalar_logit_sum[-1,]
  rownames(comp4_scalar_logit_sum)<-c(1:length(comp4_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp4_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp4_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp4_scalar_logit_sum, file='comp4_scalar_logit_sum.csv')
  
  ##########
  # Comparison 5 (l1hs vs denovo): Logistic regression for scalar predictors
  comp5_scalar<-which(func_5=="n")
  comp5_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp5_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp5_scalar){
    #i=comp1_scalar[1]
    l1<-result_mean@features[[i]][[3]]
    ct<-result_mean@features[[i]][[1]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("feature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp5_scalar_logit_sum<-rbind(comp5_scalar_logit_sum,fit_sum)
  }
  comp5_scalar_logit_sum<-comp5_scalar_logit_sum[-1,]
  rownames(comp5_scalar_logit_sum)<-c(1:length(comp5_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp5_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp5_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp5_scalar_logit_sum, file='comp5_scalar_logit_sum.csv')
  
  ##########
  # Comparison 6 (l1hs vs l1pol): Logistic regression for scalar predictors
  comp6_scalar<-which(func_6=="n")
  comp6_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp6_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp6_scalar){
    #i=comp1_scalar[1]
    l1<-result_mean@features[[i]][[3]]
    ct<-result_mean@features[[i]][[2]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("feature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp6_scalar_logit_sum<-rbind(comp6_scalar_logit_sum,fit_sum)
  }
  comp6_scalar_logit_sum<-comp6_scalar_logit_sum[-1,]
  rownames(comp6_scalar_logit_sum)<-c(1:length(comp6_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp6_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp6_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp6_scalar_logit_sum, file='comp6_scalar_logit_sum.csv')
  
  
  # Plot the pseudo-R2 for individual logistic regressions of all comparisons 
  #par(mfrow=c(2,3))
  # Comparison1
  # plot(comp1_scalar_logit_sum$pseudoR2, xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='de novo L1s vs Control')
  # text(comp1_scalar_logit_sum$pseudoR2,labels=comp1_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=0.2,col=3)
  # # Comparison2
  # plot(comp2_scalar_logit_sum$pseudoR2,xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='Polmorphic L1s vs Control')
  # text(comp2_scalar_logit_sum$pseudoR2,labels=comp2_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=0.2,col=3)
  # # Comparison3
  # plot(comp3_scalar_logit_sum$pseudoR2,xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='Human-specific L1s vs Control')
  # text(comp3_scalar_logit_sum$pseudoR2,labels=comp3_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=0.2,col=3)
  # 
  # # Comparison4
  # plot(comp4_scalar_logit_sum$pseudoR2,xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='de novo L1s vs Human-specific L1s')
  # text(comp4_scalar_logit_sum$pseudoR2,labels=comp4_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=0.2,col=3)
  # 
  # # Comparison5
  # plot(comp5_scalar_logit_sum$pseudoR2,xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='de novo L1s vs Human-specific L1s')
  # text(comp5_scalar_logit_sum$pseudoR2,labels=comp5_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=0.2,col=3)
  # 
  # # Comparison6
  # plot(comp6_scalar_logit_sum$pseudoR2,xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='Polymorphic L1s vs Human-specific L1s')
  # text(comp6_scalar_logit_sum$pseudoR2,labels=comp6_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=0.2,col=3)
  # 
  # # Plot the (log) significance for individual logistic regressions of all comparisons 
  # #par(mfrow=c(3,2))
  # # Comparison1
  # plot(log(as.numeric(comp1_scalar_logit_sum$significance)), xlab='Feature Number', ylab='P-Value (Log)', ylim=c(-250,50),cex=0.6,main='de novo L1s vs Control')
  # text(log(as.numeric(comp1_scalar_logit_sum$significance)),labels=comp1_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=log(0.05),col=2)
  # # Comparison2
  # plot(log(as.numeric(comp2_scalar_logit_sum$significance)), xlab='Feature Number', ylab='P-Value (Log)',ylim=c(-250,50),cex=0.6,main='Polymorphic L1s vs Control')
  # text(log(as.numeric(comp2_scalar_logit_sum$significance)),labels=comp1_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=log(0.05),col=2)
  # 
  # # Comparison3
  # plot(log(as.numeric(comp3_scalar_logit_sum$significance)),xlab='Feature Number', ylab='P-Value (Log)',ylim=c(-250,50), cex=0.6,main='Human-specific L1s vs Control')
  # text(log(as.numeric(comp3_scalar_logit_sum$significance)),labels=comp3_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=log(0.05),col=2)
  # 
  # # Comparison4
  # plot(log(as.numeric(comp4_scalar_logit_sum$significance)),xlab='Feature Number', ylab='P-Value (Log)',ylim=c(-250,50),cex=0.6, main='de novo L1s vs Polymorphic L1s')
  # text(log(as.numeric(comp4_scalar_logit_sum$significance)),labels=comp4_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=log(0.05),col=2)
  # 
  # # Comparison5
  # plot(log(as.numeric(comp5_scalar_logit_sum$significance)),xlab='Feature Number', ylab='P-Value (Log)',ylim=c(-250,50),cex=0.6, main='de novo L1s vs Human-specific L1s')
  # text(log(as.numeric(comp5_scalar_logit_sum$significance)),labels=comp5_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=log(0.05),col=2)
  # 
  # # Comparison6
  # plot(log(as.numeric(comp6_scalar_logit_sum$significance)),xlab='Feature Number', ylab='P-Value (Log)',ylim=c(-250,50),cex=0.6, main='Polymorphic L1s vs Human-specific L1s')
  # text(log(as.numeric(comp6_scalar_logit_sum$significance)),labels=comp6_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # abline(h=log(0.05),col=2)
  
  ##Save the image
  #I=paste('scalar_noLowres_individual_regression_random',i,'.RData',sep='')
  save.image(paste('scalar_noLowres_individual_regression_random',r,'.RData',sep=''))
}

