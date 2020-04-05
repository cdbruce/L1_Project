setwd('~/Google_Drive/l1/L1_Project/Analysis/sFLR/')
library(IWTomics)
require(dendextend)

########## Transform all 45 high-resolution features

load('L1_autosomes_zero_fixed_with_strand_smoothed.RData')
#load('L1_autosomes_results_smoothed_mean_10.RData')
#load('L1_autosomes_results_smoothed_mean_1.RData')

# Histogram of the data before transformation
pdf('Histogram_pre_transformation.pdf')
par(mfrow=c(2,2))
for (i in 1:45) {
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[1]]), main = paste('L1denovo_feature_',i), xlab = paste('Feature',i))
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[2]]), main = paste('L1Pol_feature_',i), xlab = paste('Feature',i))
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[3]]), main = paste('L1HS_feature_',i), xlab = paste('Feature',i))
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[4]]), main = paste('Control_feature_',i), xlab = paste('Feature',i))
}
dev.off()


#length(result_mean@features$H2AFZ_signal$L1denovo)
#class(result_mean@features$H2AFZ_signal$L1denovo)
#dim(result_mean@features$H2AFZ_signal$L1denovo)
#test<-as.vector(result_mean@features$H2AFZ_signal$L1denovo)
#class(test)
#hist(test)
#test2<-log(test+0.001)
#test_new<-matrix(data=test2,nrow=100)
#dim(test_new)


#Shapiro-Wilk Normality Test (limited no more than 5000 samples)
#do one shift, test across four datasets, choose lowest pval
#try all possible shifts, choose the highest pval 
# length(test)
# a<-shapiro.test(head(test))
# a$p.value

#Try Kolmogorov-Smirnov Tests (ks.test) -> does not allow repeated values
# b<-ks.test(test,'pnorm')
# b$p.value

#Try Anderson-Darling test
#install.packages("nortest")
library(nortest)
#c<-ad.test(test)
#c$p.value
#Probably gonna work!

#Loop through each of the 45 features
trans_summary_ad <- as.data.frame(matrix('shift',nrow = 12))
trans_summary_sw <- as.data.frame(matrix('shift',nrow = 12))

for (i in 1:45){
  #i=1
  #Create two empty dataframes, df_sum for each feature, trans_summary to append all features
  df_sum_ad<- data.frame()
  df_sum_sw<- data.frame()
  
  ##No transformation, subsample 5000 values from each vector
  dv<-lapply(1:10,function(k) sample(as.vector(regionsFeatures_smoothed@features[[i]][[1]]), replace=FALSE, 5000))
  pl<-lapply(1:10,function(k) sample(as.vector(regionsFeatures_smoothed@features[[i]][[2]]), replace=FALSE, 5000))
  hs<-lapply(1:10,function(k) sample(as.vector(regionsFeatures_smoothed@features[[i]][[3]]), replace=FALSE, 5000))
  ct<-lapply(1:10,function(k) sample(as.vector(regionsFeatures_smoothed@features[[i]][[4]]), replace=FALSE, 5000))
  
 # Test the pipeline on the first 100 values of each set  
 # dv<-head(as.vector(regionsFeatures_smoothed@features[[i]][[1]]),100)
 # pl<-head(as.vector(regionsFeatures_smoothed@features[[i]][[2]]),100)
 # hs<-head(as.vector(regionsFeatures_smoothed@features[[i]][[3]]),100)
 # ct<-head(as.vector(regionsFeatures_smoothed@features[[i]][[4]]),100)
  
  # Anderson-Darling test
  dv.ad.ts<-lapply(dv,ad.test)
  pl.ad.ts<-lapply(pl,ad.test)
  hs.ad.ts<-lapply(hs,ad.test)
  ct.ad.ts<-lapply(ct,ad.test)
  
  # Shapiro-Wilk test
   dv.sw.ts<-lapply(dv,shapiro.test)
   pl.sw.ts<-lapply(pl,shapiro.test)
   hs.sw.ts<-lapply(hs,shapiro.test)
   ct.sw.ts<-lapply(ct,shapiro.test)
  
  # Anderson-Darling
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.ad.ts),ncol=4,byrow=TRUE)[,2]))
  D1<-data.frame(data,pval)
  D1[,1]
  D1[,2]
  min_data<-as.character(with(D1, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-0
  df_sum_ad <- rbind(df_sum_ad, s)
  
  # Shapiro-Wilk
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.sw.ts),ncol=4,byrow=TRUE)[,2]))
  D2<-data.frame(data,pval)
  D2[,1]
  D2[,2]
  min_data<-as.character(with(D2, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-0
  df_sum_sw <- rbind(df_sum_sw, s)

##Log transformation with different shifts from 1 to 10^-10
  for (j in (c(1/(10^(0:10))))){
    #j=1
    dv_log<-lapply(dv,function(dvk) log(dvk+j))
    pl_log<-lapply(pl,function(plk) log(plk+j))
    hs_log<-lapply(hs,function(hsk) log(hsk+j))
    ct_log<-lapply(ct,function(ctk) log(ctk+j))
    
    # Anderson-Darling 
    dv_log.ad.ts<-lapply(dv_log,ad.test)
    pl_log.ad.ts<-lapply(pl_log,ad.test)
    hs_log.ad.ts<-lapply(hs_log,ad.test)
    ct_log.ad.ts<-lapply(ct_log,ad.test)
    
    # Shapiro-Wilk test
    dv_log.sw.ts<-lapply(dv_log,shapiro.test)
    pl_log.sw.ts<-lapply(pl_log,shapiro.test)
    hs_log.sw.ts<-lapply(hs_log,shapiro.test)
    ct_log.sw.ts<-lapply(ct_log,shapiro.test)
    
    # Anderson-Darling
    data<-c('dv','pl','hs','ct')
    pval<-c(min(matrix(unlist(dv_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.ad.ts),ncol=4,byrow=TRUE)[,2]))
    D1<-data.frame(data,pval)
    D1[,1]
    D1[,2]
    min_data<-as.character(with(D1, data[which.min(pval)]))
    min_pval<-min(pval)
    s<-data.frame(min_data, min_pval)
    colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
    row.names(s)<-j
    df_sum_ad <- rbind(df_sum_ad, s)
    
    # Shapiro-Wilk
    data<-c('dv','pl','hs','ct')
    pval<-c(min(matrix(unlist(dv_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.sw.ts),ncol=4,byrow=TRUE)[,2]))
    D2<-data.frame(data,pval)
    D2[,1]
    D2[,2]
    min_data<-as.character(with(D2, data[which.min(pval)]))
    min_pval<-min(pval)
    s<-data.frame(min_data, min_pval)
    colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
    row.names(s)<-j
    df_sum_sw <- rbind(df_sum_sw, s)
  }
  trans_summary_ad<-cbind(trans_summary_ad,df_sum_ad)
  trans_summary_sw<-cbind(trans_summary_sw,df_sum_sw)
  assign(paste("df_sum_ad_", i, sep = ""), df_sum_ad )
  assign(paste("df_sum_sw_", i, sep = ""), df_sum_sw )
}  
save.image("transformation_final.RData")

#write.table(trans_summary_ad, file='adtest_transformation_summary_random10.txt', sep='\t')
write.table(trans_summary_ad, file='adtest_transformation_summary_random5000.txt', sep='\t')
#write.table(trans_summary_sw, file='swtest_transformation_summary_random10.txt', sep='\t')
write.table(trans_summary_sw, file='swtest_transformation_summary_random5000.txt', sep='\t')
#row.names(trans_summary_sw)<-c(0,1/(10^(0:10)))


##Select the transformation with largest p-value for each feature
#load('transformation_final.RData')
head(trans_summary_sw)
which.max(trans_summary_sw[,24])
s=seq(from = 3,to = 91, by=2)
pvals<-trans_summary_sw[,s]
row.names(pvals)
trans_shift<-rownames(pvals)[apply(pvals, 2, which.max)]
# Make a matrix out of the trans_shift
trans <- as.vector(trans_shift, mode = "numeric")
shift_table<-cbind(c(1:45),trans)

test<-log(regionsFeatures_smoothed@features[[1]][[1]]+1)
head(as.vector(test))
head(as.vector(regionsFeatures_smoothed@features[[1]][[1]]))

# Variable transformation for the complete dataset (with all de novo L1s) 
# trans[1:12]
for (i in c(1:24, 26:43)){
  print(i)
  print(trans[i])
}

for (i in c(1:24, 26:43)){
  regionsFeatures_smoothed@features[[i]][[1]]<-log(regionsFeatures_smoothed@features[[i]][[1]]+trans[i])
  regionsFeatures_smoothed@features[[i]][[2]]<-log(regionsFeatures_smoothed@features[[i]][[2]]+trans[i])
  regionsFeatures_smoothed@features[[i]][[3]]<-log(regionsFeatures_smoothed@features[[i]][[3]]+trans[i])
  regionsFeatures_smoothed@features[[i]][[4]]<-log(regionsFeatures_smoothed@features[[i]][[4]]+trans[i])
}

save(regionsFeatures_smoothed,file='L1_allFeatures_transformed.RData')

# Histogram of the data after transformation
# load('L1_allFeatures_transformed.RData')
pdf('Histogram_post_transformation.pdf')
par(mfrow=c(2,2))
for (i in 1:45) {
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[1]]), main = paste('L1denovo_feature_',i), xlab = paste('Feature',i))
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[2]]), main = paste('L1Pol_feature_',i), xlab = paste('Feature',i))
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[3]]), main = paste('L1HS_feature_',i), xlab = paste('Feature',i))
  hist(as.vector(regionsFeatures_smoothed@features[[i]][[4]]), main = paste('Control_feature_',i), xlab = paste('Feature',i))
}
dev.off()


# Variable transformation for the 10 random samples (of de novo L1s)

for (i in c(1:24, 26:43)){
  regionsFeatures_smoothed@features[[i]][[1]]<-log(regionsFeatures_smoothed@features[[i]][[1]]+trans[i])
  regionsFeatures_smoothed@features[[i]][[2]]<-log(regionsFeatures_smoothed@features[[i]][[2]]+trans[i])
  regionsFeatures_smoothed@features[[i]][[3]]<-log(regionsFeatures_smoothed@features[[i]][[3]]+trans[i])
  regionsFeatures_smoothed@features[[i]][[4]]<-log(regionsFeatures_smoothed@features[[i]][[4]]+trans[i])
}

save(regionsFeatures_smoothed,file='L1_allFeatures_transformed.RData')

setwd("10_randoms/")
for (r in 1:10){
  file=paste('L1_autosomes_results_smoothed_mean_', r, '.RData', sep="")
  load(file)
  # Plot the variables before transformation 
  pdf(paste('Histogram_pre_transformation_random_', r, '.pdf', sep=""))
  par(mfrow=c(2,2))
  for (i in 1:45) {
    hist(as.vector(result_mean@features[[i]][[1]]), main = paste('L1denovo_feature_',i), xlab = paste('Feature',i))
    hist(as.vector(result_mean@features[[i]][[2]]), main = paste('L1Pol_feature_',i), xlab = paste('Feature',i))
    hist(as.vector(result_mean@features[[i]][[3]]), main = paste('L1HS_feature_',i), xlab = paste('Feature',i))
    hist(as.vector(result_mean@features[[i]][[4]]), main = paste('Control_feature_',i), xlab = paste('Feature',i))
  }
  dev.off()
    
  # Transform the variables based on the vector 'trans'
  for (i in c(1:24, 26:43)){
    result_mean@features[[i]][[1]]<-log(result_mean@features[[i]][[1]]+trans[i])
    result_mean@features[[i]][[2]]<-log(result_mean@features[[i]][[2]]+trans[i])
    result_mean@features[[i]][[3]]<-log(result_mean@features[[i]][[3]]+trans[i])
    result_mean@features[[i]][[4]]<-log(result_mean@features[[i]][[4]]+trans[i])
  }
  pdf(paste('Histogram_post_transformation_random_', r, '.pdf', sep=""))
  par(mfrow=c(2,2))
  for (i in 1:45) {
    hist(as.vector(result_mean@features[[i]][[1]]), main = paste('L1denovo_feature_',i), xlab = paste('Feature',i))
    hist(as.vector(result_mean@features[[i]][[2]]), main = paste('L1Pol_feature_',i), xlab = paste('Feature',i))
    hist(as.vector(result_mean@features[[i]][[3]]), main = paste('L1HS_feature_',i), xlab = paste('Feature',i))
    hist(as.vector(result_mean@features[[i]][[4]]), main = paste('Control_feature_',i), xlab = paste('Feature',i))
  }
  dev.off()
  save(result_mean,file=paste('L1_transformed_random_', r, '.RData', sep=""))
}

# Write a table of "shifts" selected for the 45 "high-resolution" features
shift_table
regionsFeatures_smoothed@metadata$feature_datasets 
f<-as.vector(regionsFeatures_smoothed@metadata$feature_datasets[,1])
shift_sum<-cbind(f,shift_table)
shift_sum<-data.frame(shift_sum)
write.table(shift_sum,file = "~/Google_Drive/l1/L1_Project/Analysis/sFLR/shift_table.txt",sep="\t",quote=FALSE)








# #pval_random_median[[i]][[j]]$adjusted_pval_matrix <- matrix(mapply(function(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10){median(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))},pval_random_1[[i]][[j]]$adjusted_pval_matrix,pval_random_2[[i]][[j]]$adjusted_pval_matrix,
#                                                                    pval_random_3[[i]][[j]]$adjusted_pval_matrix,pval_random_4[[i]][[j]]$adjusted_pval_matrix,
#                                                                    pval_random_5[[i]][[j]]$adjusted_pval_matrix,pval_random_6[[i]][[j]]$adjusted_pval_matrix,
#                                                                    pval_random_7[[i]][[j]]$adjusted_pval_matrix,pval_random_8[[i]][[j]]$adjusted_pval_matrix,
#                                                                    pval_random_9[[i]][[j]]$adjusted_pval_matrix,pval_random_10[[i]][[j]]$adjusted_pval_matrix),nrow=100)




###################################################################################################
################Transform the 7 low resolution features: 10 ramdom samples with the same de novo regions ##########
###Note that feature 3 has negative values, need to be transformed seperately!!!!! Here only tested for features 1,2,4:7

## Transform the 7 low resolution features base on the original dataset without subsampling denovos
setwd('~/Google_Drive/l1/L1_Project/Analysis/sFLR/')
load('L1_low_resolution_autosomes.RData')
library(nortest)


# Histogram of the low resolution features before transformation
pdf('LowRes_Histogram_pre_transformation.pdf')
par(mfrow=c(2,2))
for (i in 1:7) {
  hist(as.vector(regionsFeatures@features[[i]][[1]]), main = paste('L1denovo_feature_',i+45), xlab = paste('Feature',i+45))
  hist(as.vector(regionsFeatures@features[[i]][[2]]), main = paste('L1Pol_feature_',i+45), xlab = paste('Feature',i+45))
  hist(as.vector(regionsFeatures@features[[i]][[3]]), main = paste('L1HS_feature_',i+45), xlab = paste('Feature',i+45))
  hist(as.vector(regionsFeatures@features[[i]][[4]]), main = paste('Control_feature_',i+45), xlab = paste('Feature',i+45))
}
dev.off()


#Loop through each of the 7 features: use shift 1...e^-10 for features 1:2, 4:7 
#and shift 2:20 for feature 3!!!!

trans_summary_ad <- as.data.frame(matrix('shift',nrow = 12))
trans_summary_sw <- as.data.frame(matrix('shift',nrow = 12))

#First 2 features: feature 1 and 2
for (i in c(1,2)){
  #i=1
  #Create two empty dataframes, df_sum for each feature, trans_summary to append all features
  df_sum_ad<- data.frame()
  df_sum_sw<- data.frame()
  
  ##No transformation, subsample 800 values from each vector
  dv<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[1]]), replace=FALSE, 800))
  pl<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[2]]), replace=FALSE, 800))
  hs<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[3]]), replace=FALSE, 800))
  ct<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[4]]), replace=FALSE, 800))
  
  # Test the pipeline on the first 100 values of each set  
  # dv<-head(as.vector(regionsFeatures@features[[i]][[1]]),100)
  # pl<-head(as.vector(regionsFeatures@features[[i]][[2]]),100)
  # hs<-head(as.vector(regionsFeatures@features[[i]][[3]]),100)
  # ct<-head(as.vector(regionsFeatures@features[[i]][[4]]),100)
  
  # Anderson-Darling test
  dv.ad.ts<-lapply(dv,ad.test)
  pl.ad.ts<-lapply(pl,ad.test)
  hs.ad.ts<-lapply(hs,ad.test)
  ct.ad.ts<-lapply(ct,ad.test)
  
  # Shapiro-Wilk test
  dv.sw.ts<-lapply(dv,shapiro.test)
  pl.sw.ts<-lapply(pl,shapiro.test)
  hs.sw.ts<-lapply(hs,shapiro.test)
  ct.sw.ts<-lapply(ct,shapiro.test)
  
  # Anderson-Darling
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.ad.ts),ncol=4,byrow=TRUE)[,2]))
  D1<-data.frame(data,pval)
  D1[,1]
  D1[,2]
  min_data<-as.character(with(D1, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-0
  df_sum_ad <- rbind(df_sum_ad, s)
  
  # Shapiro-Wilk
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.sw.ts),ncol=4,byrow=TRUE)[,2]))
  D2<-data.frame(data,pval)
  D2[,1]
  D2[,2]
  min_data<-as.character(with(D2, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-0
  df_sum_sw <- rbind(df_sum_sw, s)
  
  ##Log transformation with different shifts from 1 to 10^-10
  for (j in (c(1/(10^(0:10))))){
    #j=1
    dv_log<-lapply(dv,function(dvk) log(dvk+j))
    pl_log<-lapply(pl,function(plk) log(plk+j))
    hs_log<-lapply(hs,function(hsk) log(hsk+j))
    ct_log<-lapply(ct,function(ctk) log(ctk+j))
    
    # Anderson-Darling 
    dv_log.ad.ts<-lapply(dv_log,ad.test)
    pl_log.ad.ts<-lapply(pl_log,ad.test)
    hs_log.ad.ts<-lapply(hs_log,ad.test)
    ct_log.ad.ts<-lapply(ct_log,ad.test)
    
    # Shapiro-Wilk test
    dv_log.sw.ts<-lapply(dv_log,shapiro.test)
    pl_log.sw.ts<-lapply(pl_log,shapiro.test)
    hs_log.sw.ts<-lapply(hs_log,shapiro.test)
    ct_log.sw.ts<-lapply(ct_log,shapiro.test)
    
    # Anderson-Darling
    data<-c('dv','pl','hs','ct')
    pval<-c(min(matrix(unlist(dv_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.ad.ts),ncol=4,byrow=TRUE)[,2]))
    D1<-data.frame(data,pval)
    D1[,1]
    D1[,2]
    min_data<-as.character(with(D1, data[which.min(pval)]))
    min_pval<-min(pval)
    s<-data.frame(min_data, min_pval)
    colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
    row.names(s)<-j
    df_sum_ad <- rbind(df_sum_ad, s)
    
    # Shapiro-Wilk
    data<-c('dv','pl','hs','ct')
    pval<-c(min(matrix(unlist(dv_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.sw.ts),ncol=4,byrow=TRUE)[,2]))
    D2<-data.frame(data,pval)
    D2[,1]
    D2[,2]
    min_data<-as.character(with(D2, data[which.min(pval)]))
    min_pval<-min(pval)
    s<-data.frame(min_data, min_pval)
    colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
    row.names(s)<-j
    df_sum_sw <- rbind(df_sum_sw, s)
  }
  trans_summary_ad<-cbind(trans_summary_ad,df_sum_ad)
  trans_summary_sw<-cbind(trans_summary_sw,df_sum_sw)
  assign(paste("df_sum_ad_", i, sep = ""), df_sum_ad )
  assign(paste("df_sum_sw_", i, sep = ""), df_sum_sw )
}  


## Feature 3
#######Replication Timing: Feature 3
i=3
## for replication timing, since there are negative values no smaller than '-2',  try shifts from 2:20 (step2)
# min(regionsFeatures@features$Rep_Timing[[1]])
# min(regionsFeatures@features$Rep_Timing[[2]])
# min(regionsFeatures@features$Rep_Timing[[3]])
# min(regionsFeatures@features$Rep_Timing[[4]])
#Use the previous dataframes: trans_summary to append all features
trans_summary_ad
trans_summary_sw
##No transformation, subsample 800 values from each vector
df_sum_ad<- data.frame()
df_sum_sw<- data.frame()

dv<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[1]]), replace=FALSE, 800))
pl<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[2]]), replace=FALSE, 800))
hs<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[3]]), replace=FALSE, 800))
ct<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[4]]), replace=FALSE, 800))

# Test the pipeline on the first 100 values of each set  
# dv<-head(as.vector(regionsFeatures@features[[i]][[1]]),100)
# pl<-head(as.vector(regionsFeatures@features[[i]][[2]]),100)
# hs<-head(as.vector(regionsFeatures@features[[i]][[3]]),100)
# ct<-head(as.vector(regionsFeatures@features[[i]][[4]]),100)

# Anderson-Darling test
dv.ad.ts<-lapply(dv,ad.test)
pl.ad.ts<-lapply(pl,ad.test)
hs.ad.ts<-lapply(hs,ad.test)
ct.ad.ts<-lapply(ct,ad.test)

# Shapiro-Wilk test
dv.sw.ts<-lapply(dv,shapiro.test)
pl.sw.ts<-lapply(pl,shapiro.test)
hs.sw.ts<-lapply(hs,shapiro.test)
ct.sw.ts<-lapply(ct,shapiro.test)

# Anderson-Darling
data<-c('dv','pl','hs','ct')
pval<-c(min(matrix(unlist(dv.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.ad.ts),ncol=4,byrow=TRUE)[,2]))
D1<-data.frame(data,pval)
D1[,1]
D1[,2]
min_data<-as.character(with(D1, data[which.min(pval)]))
min_pval<-min(pval)
s<-data.frame(min_data, min_pval)
colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
row.names(s)<-0
df_sum_ad <- rbind(df_sum_ad, s)

# Shapiro-Wilk
data<-c('dv','pl','hs','ct')
pval<-c(min(matrix(unlist(dv.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.sw.ts),ncol=4,byrow=TRUE)[,2]))
D2<-data.frame(data,pval)
D2[,1]
D2[,2]
min_data<-as.character(with(D2, data[which.min(pval)]))
min_pval<-min(pval)
s<-data.frame(min_data, min_pval)
colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
row.names(s)<-0
df_sum_sw <- rbind(df_sum_sw, s)

##Log transformation with different shifts from 2 to 20
for (j in (seq(2,22,2))){
  #j=1
  dv_log<-lapply(dv,function(dvk) log(dvk+j))
  pl_log<-lapply(pl,function(plk) log(plk+j))
  hs_log<-lapply(hs,function(hsk) log(hsk+j))
  ct_log<-lapply(ct,function(ctk) log(ctk+j))
  
  # Anderson-Darling 
  dv_log.ad.ts<-lapply(dv_log,ad.test)
  pl_log.ad.ts<-lapply(pl_log,ad.test)
  hs_log.ad.ts<-lapply(hs_log,ad.test)
  ct_log.ad.ts<-lapply(ct_log,ad.test)
  
  # Shapiro-Wilk test
  dv_log.sw.ts<-lapply(dv_log,shapiro.test)
  pl_log.sw.ts<-lapply(pl_log,shapiro.test)
  hs_log.sw.ts<-lapply(hs_log,shapiro.test)
  ct_log.sw.ts<-lapply(ct_log,shapiro.test)
  
  # Anderson-Darling
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.ad.ts),ncol=4,byrow=TRUE)[,2]))
  D1<-data.frame(data,pval)
  D1[,1]
  D1[,2]
  min_data<-as.character(with(D1, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-j
  df_sum_ad <- rbind(df_sum_ad, s)
  
  # Shapiro-Wilk
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.sw.ts),ncol=4,byrow=TRUE)[,2]))
  D2<-data.frame(data,pval)
  D2[,1]
  D2[,2]
  min_data<-as.character(with(D2, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-j
  df_sum_sw <- rbind(df_sum_sw, s)
}
trans_summary_ad<-cbind(trans_summary_ad,df_sum_ad)
trans_summary_sw<-cbind(trans_summary_sw,df_sum_sw)
assign(paste("df_sum_ad_", i, sep = ""), df_sum_ad )
assign(paste("df_sum_sw_", i, sep = ""), df_sum_sw )


## Feature 4-7, use the same strategy as feature 1 and feature 2
for (i in c(4:7)){
  #i=1
  #Create two empty dataframes, df_sum for each feature, trans_summary to append all features
  df_sum_ad<- data.frame()
  df_sum_sw<- data.frame()
  
  ##No transformation, subsample 800 values from each vector
  dv<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[1]]), replace=FALSE, 800))
  pl<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[2]]), replace=FALSE, 800))
  hs<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[3]]), replace=FALSE, 800))
  ct<-lapply(1:10,function(k) sample(as.vector(regionsFeatures@features[[i]][[4]]), replace=FALSE, 800))
  
  # Test the pipeline on the first 100 values of each set  
  # dv<-head(as.vector(regionsFeatures@features[[i]][[1]]),100)
  # pl<-head(as.vector(regionsFeatures@features[[i]][[2]]),100)
  # hs<-head(as.vector(regionsFeatures@features[[i]][[3]]),100)
  # ct<-head(as.vector(regionsFeatures@features[[i]][[4]]),100)
  
  # Anderson-Darling test
  dv.ad.ts<-lapply(dv,ad.test)
  pl.ad.ts<-lapply(pl,ad.test)
  hs.ad.ts<-lapply(hs,ad.test)
  ct.ad.ts<-lapply(ct,ad.test)
  
  # Shapiro-Wilk test
  dv.sw.ts<-lapply(dv,shapiro.test)
  pl.sw.ts<-lapply(pl,shapiro.test)
  hs.sw.ts<-lapply(hs,shapiro.test)
  ct.sw.ts<-lapply(ct,shapiro.test)
  
  # Anderson-Darling
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.ad.ts),ncol=4,byrow=TRUE)[,2]))
  D1<-data.frame(data,pval)
  D1[,1]
  D1[,2]
  min_data<-as.character(with(D1, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-0
  df_sum_ad <- rbind(df_sum_ad, s)
  
  # Shapiro-Wilk
  data<-c('dv','pl','hs','ct')
  pval<-c(min(matrix(unlist(dv.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct.sw.ts),ncol=4,byrow=TRUE)[,2]))
  D2<-data.frame(data,pval)
  D2[,1]
  D2[,2]
  min_data<-as.character(with(D2, data[which.min(pval)]))
  min_pval<-min(pval)
  s<-data.frame(min_data, min_pval)
  colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
  row.names(s)<-0
  df_sum_sw <- rbind(df_sum_sw, s)
  
  ##Log transformation with different shifts from 1 to 10^-10
  for (j in (c(1/(10^(0:10))))){
    #j=1
    dv_log<-lapply(dv,function(dvk) log(dvk+j))
    pl_log<-lapply(pl,function(plk) log(plk+j))
    hs_log<-lapply(hs,function(hsk) log(hsk+j))
    ct_log<-lapply(ct,function(ctk) log(ctk+j))
    
    # Anderson-Darling 
    dv_log.ad.ts<-lapply(dv_log,ad.test)
    pl_log.ad.ts<-lapply(pl_log,ad.test)
    hs_log.ad.ts<-lapply(hs_log,ad.test)
    ct_log.ad.ts<-lapply(ct_log,ad.test)
    
    # Shapiro-Wilk test
    dv_log.sw.ts<-lapply(dv_log,shapiro.test)
    pl_log.sw.ts<-lapply(pl_log,shapiro.test)
    hs_log.sw.ts<-lapply(hs_log,shapiro.test)
    ct_log.sw.ts<-lapply(ct_log,shapiro.test)
    
    # Anderson-Darling
    data<-c('dv','pl','hs','ct')
    pval<-c(min(matrix(unlist(dv_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.ad.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.ad.ts),ncol=4,byrow=TRUE)[,2]))
    D1<-data.frame(data,pval)
    D1[,1]
    D1[,2]
    min_data<-as.character(with(D1, data[which.min(pval)]))
    min_pval<-min(pval)
    s<-data.frame(min_data, min_pval)
    colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
    row.names(s)<-j
    df_sum_ad <- rbind(df_sum_ad, s)
    
    # Shapiro-Wilk
    data<-c('dv','pl','hs','ct')
    pval<-c(min(matrix(unlist(dv_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(pl_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(hs_log.sw.ts),ncol=4,byrow=TRUE)[,2]), min(matrix(unlist(ct_log.sw.ts),ncol=4,byrow=TRUE)[,2]))
    D2<-data.frame(data,pval)
    D2[,1]
    D2[,2]
    min_data<-as.character(with(D2, data[which.min(pval)]))
    min_pval<-min(pval)
    s<-data.frame(min_data, min_pval)
    colnames(s)<- c(paste("dataset_feature_",i, sep='' ), paste("pval_feature_",i, sep=''))
    row.names(s)<-j
    df_sum_sw <- rbind(df_sum_sw, s)
  }
  trans_summary_ad<-cbind(trans_summary_ad,df_sum_ad)
  trans_summary_sw<-cbind(trans_summary_sw,df_sum_sw)
  assign(paste("df_sum_ad_", i, sep = ""), df_sum_ad )
  assign(paste("df_sum_sw_", i, sep = ""), df_sum_sw )
}  

#save.image("low_resolution_transformation_final.RData")

save.image("low_resolution_transformation_prep.RData")




######## Transform the 7 low-resolution variables
## Normality test
##Select the transformation with largest p-value for each feature
load("low_resolution_transformation_prep.RData")
head(trans_summary_sw)
which.max(trans_summary_sw[,3])
s=seq(from = 3,to = 15, by=2)
pvals<-trans_summary_sw[,s]
row.names(pvals)
trans_shift<-rownames(pvals)[apply(pvals, 2, which.max)]
# Make a matrix out of the trans_shift
trans <- as.vector(trans_shift, mode = "numeric")
shift_table<-cbind(c(1:7),trans)

## Here to change the shift value for feature3 (not ideal though) 
# to the real shift value from the seq (2:22, step=2)
# The shift "1.0" of feature 3 is acutally shift "4" seq(2,22,2) 
shift_table[3,]<-4
trans[3] <- 4

# test<-log(regionsFeatures@features[[1]][[1]]+1)
# head(as.vector(test))
# head(as.vector(regionsFeatures@features[[1]][[1]]))

# Variable transformation for the complete dataset (with all de novo L1s)
# Only choose the non-zero values from the vector "trans", 
# since the "0" means no transformation
# > trans[1:7]
# [1] 0.0 0.0 4.0 0.1 0.1 0.1 0.0

for (i in c(3,4,5,6)){
  #print(i)
  print(trans[i])
}

for (i in c(3,4,5,6)){
  regionsFeatures@features[[i]][[1]]<-log(regionsFeatures@features[[i]][[1]]+trans[i])
  regionsFeatures@features[[i]][[2]]<-log(regionsFeatures@features[[i]][[2]]+trans[i])
  regionsFeatures@features[[i]][[3]]<-log(regionsFeatures@features[[i]][[3]]+trans[i])
  regionsFeatures@features[[i]][[4]]<-log(regionsFeatures@features[[i]][[4]]+trans[i])
}

#save(regionsFeatures,file='L1_Lowr_Features_transformed.RData')

# Histogram of the data after transformation
# load('L1_allFeatures_transformed.RData')
pdf('LowRes_Histogram_post_transformation.pdf')
par(mfrow=c(2,2))
for (i in 1:7) {
  hist(as.vector(regionsFeatures@features[[i]][[1]]), main = paste('L1denovo_feature_',i+45), xlab = paste('Feature',i+45))
  hist(as.vector(regionsFeatures@features[[i]][[2]]), main = paste('L1Pol_feature_',i+45), xlab = paste('Feature',i+45))
  hist(as.vector(regionsFeatures@features[[i]][[3]]), main = paste('L1HS_feature_',i+45), xlab = paste('Feature',i+45))
  hist(as.vector(regionsFeatures@features[[i]][[4]]), main = paste('Control_feature_',i+45), xlab = paste('Feature',i+45))
}
dev.off()


#write.table(trans_summary_ad, file='adtest_transformation_summary_random10.txt', sep='\t')
write.table(trans_summary_ad, file='low_resolution_adtest_transformation_summary_random800.txt', sep='\t')
#write.table(trans_summary_sw, file='swtest_transformation_summary_random10.txt', sep='\t')
write.table(trans_summary_sw, file='low_resolution_swtest_transformation_summary_random800.txt', sep='\t')
#row.names(trans_summary_sw)<-c(0,1/(10^(0:10)))



#############################################
################# 10 Random Samples ######
## Create (denovo) subsamples for the 7 low resolution features 
setwd('~/Google_Drive/l1/L1_Project/Analysis/sFLR/')
#load('L1_low_resolution_autosomes.RData')

# Check the transformation shift table, should be the same as above
trans

#Transformation for each random sample
for (r in 1:10){
  r=1#Since random 1 was selected in our analysis, here we only load random 1 for illustration purpose
  #load('L1_low_resolution_autosomes.RData')
  file=paste('~/Google_Drive/l1/L1_Project/Analysis/IWTomics/high_resolution/L1_autosomes_results_smoothed_mean_', r, '.RData', sep="")
  load(file)
  # result_mean@regions$L1denovo
  # regionsFeatures@regions$L1denovo
  ## number of windows in each dataset
  # lengthRegions(regionsFeatures)
  # lengthRegions(result_mean)
  
  # Create an index from the de novo regions from the 'r-th' random sample, and extract the regions for the low resolution features
  x=findOverlaps(result_mean@regions$L1denovo,regionsFeatures@regions$L1denovo,type='equal')
  #head(x)
  index=as.matrix(x)[,2]
  
  ##Use the index to match the 1000 denovo regions 
  
    ##Rearrange region index
    regionsFeatures@metadata$region_datasets['L1denovo','size']<-length(index)
    regionsFeatures@regions$L1denovo<-regionsFeatures@regions$L1denovo[index]
    ##Rearrange features
    regionsFeatures@features=lapply(regionsFeatures@features,
                                             function(feature){
                                               feature$L1denovo=matrix(feature$L1denovo[,index],nrow=1)
                                               return(feature)
                                             })
    regionsFeatures@length_features=lapply(regionsFeatures@length_features,
                                                    function(feature){
                                                      feature$L1denovo=feature$L1denovo[index]
                                                      return(feature)
                                                    })
  ## check again number of windows in each dataset
  #validObject(regionsFeatures)
  #lengthRegions(regionsFeatures)
  save(regionsFeatures,file=paste('L1_lowFeatures_raw_random',r,'.RData',sep=''))
  
  # Plot the variables before transformation
  pdf(paste('Histogram_lowFeatures_pre_transformation_random_', r, '.pdf', sep=""))
  par(mfrow=c(2,2))
  for (i in 1:7) {
    hist(as.vector(regionsFeatures@features[[i]][[1]]), main = paste('L1denovo_feature_',i+45), xlab = paste('Feature',i+45))
    hist(as.vector(regionsFeatures@features[[i]][[2]]), main = paste('L1Pol_feature_',i+45), xlab = paste('Feature',i+45))
    hist(as.vector(regionsFeatures@features[[i]][[3]]), main = paste('L1HS_feature_',i+45), xlab = paste('Feature',i+45))
    hist(as.vector(regionsFeatures@features[[i]][[4]]), main = paste('Control_feature_',i+45), xlab = paste('Feature',i+45))
  }
  dev.off()
  
  
  ## Transform the 7 low resolution features (one random sample per loop)
  for (i in 3:6){
    #print(i)
    print(trans[i])
  }
  
  for (i in 3:6){
    regionsFeatures@features[[i]][[1]]<-log(regionsFeatures@features[[i]][[1]]+trans[i])
    regionsFeatures@features[[i]][[2]]<-log(regionsFeatures@features[[i]][[2]]+trans[i])
    regionsFeatures@features[[i]][[3]]<-log(regionsFeatures@features[[i]][[3]]+trans[i])
    regionsFeatures@features[[i]][[4]]<-log(regionsFeatures@features[[i]][[4]]+trans[i])
  }
  
  save(regionsFeatures,file=paste('L1_lowFeatures_transformed_random',r,'.RData',sep=''))
  # Copy the saved Rdata to the "regression" directory, sorted by random samples 
  
  # Histogram of the data after transformation
  pdf(paste('Histogram_lowFeatures_post_transformation_random_', r, '.pdf', sep=""))
  par(mfrow=c(2,2))
  for (i in 1:7) {
    hist(as.vector(regionsFeatures@features[[i]][[1]]), main = paste('L1denovo_feature_',i+45), xlab = paste('Feature',i+45))
    hist(as.vector(regionsFeatures@features[[i]][[2]]), main = paste('L1Pol_feature_',i+45), xlab = paste('Feature',i+45))
    hist(as.vector(regionsFeatures@features[[i]][[3]]), main = paste('L1HS_feature_',i+45), xlab = paste('Feature',i+45))
    hist(as.vector(regionsFeatures@features[[i]][[4]]), main = paste('Control_feature_',i+45), xlab = paste('Feature',i+45))
  }
  dev.off()
  
}

##Double check
# regionsFeatures@features[[1]][[1]]
# regionsFeatures@features[[7]][[1]]

