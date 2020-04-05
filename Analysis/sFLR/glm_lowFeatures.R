## Logistic regression on L1 data: scalar predictors 
setwd('~/Google_Drive/L1/L1_Project/Analysis/sFLR/')
# Load the transformed data
# source("https://bioconductor.org/biocLite.R")
# biocLite("IWTomics")
require(IWTomics)
#library(fda)

####Here for comparisons 4,5 and 6, we also changed the comparisons' directions to L1pol vs denovo, L1hs vs denovo, L1pol vs L1hs

for (r in 1:10){
  r=1#Since random 1 was selected in our analysis, here we only load random 1 for illustration purpose
  setwd(paste("~/Google_Drive/L1/L1_Project/Analysis/sFLR/random",r,sep=''))
  load(file = paste("L1_lowFeatures_transformed_random",r,'.RData',sep=''))
  # Use the similar pipepline for selecting scalar predictors among high-res features, 
  # only here all the 7 low-res features are "scalar" 
  y=1
  n=0
  func_1<-as.vector(c(n,n,n,n,n,n,n))
  func_2<-as.vector(c(n,n,n,n,n,n,n))
  func_3<-as.vector(c(n,n,n,n,n,n,n))
  func_4<-as.vector(c(n,n,n,n,n,n,n))
  func_5<-as.vector(c(n,n,n,n,n,n,n))
  func_6<-as.vector(c(n,n,n,n,n,n,n))
  # test<-regionsFeatures@features$H2AFZ_signal[[1]]
  # dim(test)
  
  #Create a data table with feature names
  feature_names<-as.vector(c('Dis_Telo', 'Dis_Cent', 'Rep_Timing', 'Female_Recom_Rate', 
                             'Male_Recom_Rate', 'Average_Recom_Rate', 'Telomere_Hexamer'))
  
  feature_num<- as.vector(c(1:7))
  feature<-as.data.frame(matrix(ncol=1))
  for (i in 1:7){
    feature<-rbind(feature,paste("lowfeature_", i, sep = ""))
  }
  feature<-feature[-1,]
  feature_table<-data.frame(feature_num,feature,feature_names)
  
  #######Logsitic regression for scalar predictors
  # Comparison 1 (denovo vs control): Logistic regression for scalar predictors
  comp1_scalar<-which(func_1==0)
  comp1_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp1_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp1_scalar){
    #i=comp1_scalar[1]
    l1<-regionsFeatures@features[[i]][[1]]
    ct<-regionsFeatures@features[[i]][[4]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    # #Remove rows with NA or Infinite values (e.g. feature_52 telomere_hexamer has log(0) = -infinite after transformation )
    # DAT <- DAT[is.finite(rowSums(DAT)),]
    
    colnames(DAT)<-c("l1","pred")
    ########Something wrong with glm 
    fit<-glm(l1~pred, family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("lowfeature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp1_scalar_logit_sum<-rbind(comp1_scalar_logit_sum,fit_sum)
  }
  comp1_scalar_logit_sum<-comp1_scalar_logit_sum[-1,]
  rownames(comp1_scalar_logit_sum)<-c(1:length(comp1_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp1_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp1_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp1_scalar_logit_sum, file='comp1_lowfeature_logit_sum.csv')
  
  
  ##########
  # Comparison 2 (l1pol vs control): Logistic regression for scalar predictors
  comp2_scalar<-which(func_2==0)
  comp2_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp2_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp2_scalar){
    #i=comp1_scalar[1]
    l1<-regionsFeatures@features[[i]][[2]]
    ct<-regionsFeatures@features[[i]][[4]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    # #Remove rows with NA or Infinite values (e.g. feature_52 telomere_hexamer has log(0) = -infinite after transformation )
    # DAT <- DAT[is.finite(rowSums(DAT)),]
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("lowfeature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp2_scalar_logit_sum<-rbind(comp2_scalar_logit_sum,fit_sum)
  }
  comp2_scalar_logit_sum<-comp2_scalar_logit_sum[-1,]
  rownames(comp2_scalar_logit_sum)<-c(1:length(comp2_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp2_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp2_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp2_scalar_logit_sum, file='comp2_lowfeature_logit_sum.csv')
  
  
  ##########
  # Comparison 3 (l1HS vs control): Logistic regression for scalar predictors
  comp3_scalar<-which(func_3==0)
  comp3_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp3_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp3_scalar){
    #i=comp1_scalar[1]
    l1<-regionsFeatures@features[[i]][[3]]
    ct<-regionsFeatures@features[[i]][[4]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    # #Remove rows with NA or Infinite values (e.g. feature_52 telomere_hexamer has log(0) = -infinite after transformation )
    # DAT <- DAT[is.finite(rowSums(DAT)),]
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("lowfeature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp3_scalar_logit_sum<-rbind(comp3_scalar_logit_sum,fit_sum)
  }
  comp3_scalar_logit_sum<-comp3_scalar_logit_sum[-1,]
  rownames(comp3_scalar_logit_sum)<-c(1:length(comp3_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp3_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp3_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp3_scalar_logit_sum, file='comp3_lowfeature_logit_sum.csv')
  
  
  ##########
  # Comparison 4 (de novo vs l1pol): Logistic regression for scalar predictors
  comp4_scalar<-which(func_4==0)
  comp4_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp4_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp4_scalar){
    #i=comp1_scalar[1]
    l1<-regionsFeatures@features[[i]][[2]]
    ct<-regionsFeatures@features[[i]][[1]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    # #Remove rows with NA or Infinite values (e.g. feature_52 telomere_hexamer has log(0) = -infinite after transformation )
    # DAT <- DAT[is.finite(rowSums(DAT)),]
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("lowfeature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp4_scalar_logit_sum<-rbind(comp4_scalar_logit_sum,fit_sum)
  }
  comp4_scalar_logit_sum<-comp4_scalar_logit_sum[-1,]
  rownames(comp4_scalar_logit_sum)<-c(1:length(comp4_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp4_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp4_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp4_scalar_logit_sum, file='comp4_lowfeature_logit_sum.csv')
  
  ##########
  # Comparison 5 (de novo vs l1hs): Logistic regression for scalar predictors
  comp5_scalar<-which(func_5==0)
  comp5_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp5_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp5_scalar){
    #i=comp1_scalar[1]
    l1<-regionsFeatures@features[[i]][[3]]
    ct<-regionsFeatures@features[[i]][[1]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    # #Remove rows with NA or Infinite values (e.g. feature_52 telomere_hexamer has log(0) = -infinite after transformation )
    # DAT <- DAT[is.finite(rowSums(DAT)),]
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("lowfeature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp5_scalar_logit_sum<-rbind(comp5_scalar_logit_sum,fit_sum)
  }
  comp5_scalar_logit_sum<-comp5_scalar_logit_sum[-1,]
  rownames(comp5_scalar_logit_sum)<-c(1:length(comp5_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp5_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp5_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp5_scalar_logit_sum, file='comp5_lowfeature_logit_sum.csv')
  
  ##########
  # Comparison 6 (l1pol vs l1hs): Logistic regression for scalar predictors
  comp6_scalar<-which(func_6==0)
  comp6_scalar_logit_sum <- as.data.frame(matrix(ncol = 4))
  colnames(comp6_scalar_logit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
  for (i in comp6_scalar){
    #i=comp1_scalar[1]
    l1<-regionsFeatures@features[[i]][[3]]
    ct<-regionsFeatures@features[[i]][[2]]
    # dim(l1)
    # dim(ct)
    l1_mean<-colMeans(l1)
    ct_mean<-colMeans(ct)
    y1<-rep(1,length(l1_mean))
    y2<-rep(0,length(ct_mean))
    DAT<-data.frame(rbind(cbind(y1,l1_mean),cbind(y2,ct_mean)))
    # #Remove rows with NA or Infinite values (e.g. feature_52 telomere_hexamer has log(0) = -infinite after transformation )
    # DAT <- DAT[is.finite(rowSums(DAT)),]
    colnames(DAT)<-c("l1","pred")
    fit<-glm(l1~., family = binomial(link = logit),data=DAT)
    # fit$deviance
    # fit$null.deviance
    pseudoR2<-1-(fit$deviance/fit$null.deviance)
    coefficient<-fit$coefficients[2]
    significance<-coef(summary(fit))[,4][2]
    fit_sum<-cbind(paste("lowfeature_", i, sep = ""), coefficient, significance, pseudoR2)
    colnames(fit_sum)<-c("feature","coefficient", "significance", "pseudoR2")
    comp6_scalar_logit_sum<-rbind(comp6_scalar_logit_sum,fit_sum)
  }
  comp6_scalar_logit_sum<-comp6_scalar_logit_sum[-1,]
  rownames(comp6_scalar_logit_sum)<-c(1:length(comp6_scalar))
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp6_scalar_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp6_scalar_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp6_scalar_logit_sum, file='comp6_lowfeature_logit_sum.csv')
  
  
  # Plot the pseudo-R2 for individual logistic regressions of all comparisons
  pdf(paste('pseudoR2_glm_LowFeature_random_', r, '.pdf', sep=""))
  par(mfrow=c(3,2))

  # Comparison1: plot only 5 out of the 7 features (excluding female and male recom_rate)
  #plot(comp1_scalar_logit_sum$pseudoR2, xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='de novo L1s vs Control')
  # plot(comp1_scalar_logit_sum$pseudoR2[c(1:3,6,7)], xlab='Feature Number', ylab='PseudoR^2',ylim=c(0,0.5), main='de novo L1s vs Control')
  # text(comp1_scalar_logit_sum$pseudoR2,labels=comp1_scalar_logit_sum$feature_num, pos = 3,cex = 0.5)
  # 
  a=comp1_scalar_logit_sum$feature_names[c(1:3,6,7)]
  a<-as.vector(a)
  b<-as.vector(comp1_scalar_logit_sum[,6][c(1:3,6,7)])
  plot(b,type='n',main="Variances explained by individual regressions\nde novo L1s vs Control", 
       xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  axis(side=1,at=1:length(a),label=a,las = 3,font.axis=1.5,cex.axis=0.5)
  points(b,cex=0.4)
  
  # Plot with feature names
  # plot(comp1_scalar_logit_sum$pseudoR2,type='n',main="Variances explained by individual regressions\nde novo L1s vs Control",
  #      xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  # axis(side=1,at=1:7,label=a,las = 3,font.axis=1,cex.axis=0.4)
  # points(comp1_scalar_logit_sum$pseudoR2,cex=0.4)
  #abline(h=0.2,col=3)

  # Comparison2
  a=comp2_scalar_logit_sum$feature_names[c(1:3,6,7)]
  a<-as.vector(a)
  b<-as.vector(comp2_scalar_logit_sum[,6][c(1:3,6,7)])
  plot(b,type='n',main="Variances explained by individual regressions\nPolymorphic L1s vs Control", 
       xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  axis(side=1,at=1:length(a),label=a,las = 3,font.axis=1.5,cex.axis=0.5)
  points(b,cex=0.4)
  
  #abline(h=0.2,col=3)
  # Comparison3
  a=comp3_scalar_logit_sum$feature_names[c(1:3,6,7)]
  a<-as.vector(a)
  b<-as.vector(comp3_scalar_logit_sum[,6][c(1:3,6,7)])
  plot(b,type='n',main="Variances explained by individual regressions\nHuman-specific L1s vs Control", 
       xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  axis(side=1,at=1:length(a),label=a,las = 3,font.axis=1.5,cex.axis=0.5)
  points(b,cex=0.4)
  #abline(h=0.2,col=3)

  # Comparison4
  a=comp4_scalar_logit_sum$feature_names[c(1:3,6,7)]
  a<-as.vector(a)
  b<-as.vector(comp4_scalar_logit_sum[,6][c(1:3,6,7)])
  plot(b,type='n',main="Variances explained by individual regressions\nPolmorphic L1s vs de novo L1s", 
       xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  axis(side=1,at=1:length(a),label=a,las = 3,font.axis=1.5,cex.axis=0.5)
  points(b,cex=0.4)
  #abline(h=0.2,col=3)

  # Comparison5
  a=comp5_scalar_logit_sum$feature_names[c(1:3,6,7)]
  a<-as.vector(a)
  b<-as.vector(comp5_scalar_logit_sum[,6][c(1:3,6,7)])
  plot(b,type='n',main="Variances explained by individual regressions\nHuman-specific L1s vs de novo L1s",
       xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  axis(side=1,at=1:length(a),label=a,las = 3,font.axis=1.5,cex.axis=0.5)
  points(b,cex=0.4)
  #abline(h=0.2,col=3)

  # Comparison6
  a=comp6_scalar_logit_sum$feature_names[c(1:3,6,7)]
  a<-as.vector(a)
  b<-as.vector(comp6_scalar_logit_sum[,6][c(1:3,6,7)])
  plot(b,type='n',main="Variances explained by individual regressions\nPolmorphic L1s vs Human-specific L1s",
       xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
  axis(side=1,at=1:length(a),label=a,las = 3,font.axis=1.5,cex.axis=0.5)
  points(b,cex=0.4)
  #abline(h=0.2,col=3)
  dev.off()
  # 
  # # Plot the (log) significance for individual logistic regressions of all comparisons 
  # pdf(paste('Significance_glm_LowFeature_random_', r, '.pdf', sep=""))
  # par(mfrow=c(3,2))
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
  # 
  # dev.off()
  
  ##Save the image
  #I=paste('scalar_noLowres_individual_regression_random',i,'.RData',sep='')
  save.image(paste('low_features_individual_regression_random',r,'.RData',sep=''))
} 

# r=1
# setwd(paste("~/Desktop/regression/10_randoms/random",r,sep=''))
# load(paste('low_features_individual_regression_random',r,'.RData',sep=''))
# a=comp3_scalar_logit_sum$feature_names
# a<-as.vector(a)
