test_scalar <- function(data1,data2,mu=0,statistics='mean',probs=0.5,paired=FALSE,B=1000){
  # data1 and data2 matrices with 1 row
  
  if(statistics=='median'){
    statistics='quantile'
    probs=0.5
  }
  n1=ncol(data1)
  n2=ncol(data2)
  n=n1+n2
  data1=data1-mu
  if(paired){
    exact=(B>=(2^n1))
  }else{
    exact=(B>=choose(n,n1))
  }

  p=nrow(data1)
  result=list(test='2pop',mu=mu)
  data=cbind(data1,data2)
  
  # Univariate permutations
  message('    Point-wise tests...')
  if(statistics=='mean'){
    T0_plot=rowMeans(data1,na.rm=TRUE)-rowMeans(data2,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(statistics=='quantile'){
    T0_plot=apply(data1,1,quantile,probs=probs,na.rm=TRUE)-apply(data2,1,quantile,probs=probs,na.rm=TRUE)
    T0=(T0_plot)^2
    if(is.matrix(T0_plot)){
      T0_plot=colSums(T0_plot)
      T0=colSums(T0)
    }
  }
  if(statistics=='variance'){
    T0_plot=apply(data1,1,var,na.rm=TRUE)/apply(data2,1,var,na.rm=TRUE)
    T0=(T0_plot)^2
  }
  if(exact){
    if(paired){
      T_perm=do.call(cbind,lapply(seq_len(n1+1)-1,
                                  function(m){
                                    group_change=combn(n1,m)
                                    T_perm=apply(group_change,2,
                                                 function(change){
                                                   data_perm=data
                                                   data_perm[,c(change,n1+change)]=matrix(data_perm[,c(n1+change,change)],nrow=1)
                                                   if(statistics=='mean')
                                                     return((rowMeans(matrix(data_perm[,seq.int(n1)],nrow=1),na.rm=TRUE)-rowMeans(matrix(data_perm[,n1+seq.int(n2)],nrow=1),na.rm=TRUE))^2)
                                                   if(statistics=='quantile'){
                                                     T_perm=(apply(matrix(data_perm[,seq.int(n1)],nrow=1),1,quantile,probs=probs,na.rm=TRUE)-apply(matrix(data_perm[,n1+seq.int(n2)],nrow=1),1,quantile,probs=probs,na.rm=TRUE))^2
                                                     if(is.matrix(T_perm))
                                                       return(colSums(T_perm))
                                                     return(T_perm)
                                                   }
                                                   if(statistics=='variance')
                                                     return((apply(matrix(data_perm[,seq.int(n1)],nrow=1),1,var,na.rm=TRUE)/apply(matrix(data_perm[,n1+seq.int(n2)],nrow=1),1,var,na.rm=TRUE))^2)
                                                 })
                                    return(T_perm)
                                  }))
    }else{
      first_group=combn(n,n1)
      T_perm=apply(first_group,2,
                   function(group){
                     data_perm=matrix(data[,c(group,setdiff(seq_len(n1+n2),group))],nrow=1)
                     if(statistics=='mean')
                       return((rowMeans(matrix(data_perm[,seq.int(n1)],nrow=1),na.rm=TRUE)-rowMeans(matrix(data_perm[,n1+seq.int(n2)],nrow=1),na.rm=TRUE))^2)
                     if(statistics=='quantile'){
                       T_perm=(apply(matrix(data_perm[,seq.int(n1)],nrow=1),1,quantile,probs=probs,na.rm=TRUE)-apply(matrix(data_perm[,n1+seq.int(n2)],nrow=1),1,quantile,probs=probs,na.rm=TRUE))^2
                       if(is.matrix(T_perm))
                         return(colSums(T_perm))
                       return(T_perm)
                     }
                     if(statistics=='variance')
                       return((apply(matrix(data_perm[,seq.int(n1)],nrow=1),1,var,na.rm=TRUE)/apply(matrix(data_perm[,n1+seq.int(n2)],nrow=1),1,var,na.rm=TRUE))^2)
                   })
    }
  }else{
    T_perm=do.call(cbind,lapply(seq.int(B-1),
                                function(perm){
                                  if(paired){
                                    couple.perm=rbinom(n1,1,0.5)
                                    data_perm=matrix(data[,c(n1*couple.perm,-n1*couple.perm)+seq.int(2*n1)],nrow=1)
                                  }else{
                                    permutation=sample(n,n1)
                                    data_perm=matrix(data[,c(permutation,setdiff(seq_len(n),permutation))],nrow=1)
                                  }
                                  if(statistics=='mean')
                                    return((rowMeans(matrix(data_perm[,seq.int(n1)],nrow=1),na.rm=TRUE)-rowMeans(matrix(data_perm[,n1+seq.int(n2)],nrow=1),na.rm=TRUE))^2)
                                  if(statistics=='quantile'){
                                    T_perm=(apply(matrix(data_perm[,seq.int(n1)],nrow=1),1,quantile,probs=probs,na.rm=TRUE)-apply(matrix(data_perm[,n1+seq.int(n2)],nrow=1),1,quantile,probs=probs,na.rm=TRUE))^2
                                    if(is.matrix(T_perm))
                                      return(colSums(T_perm))
                                    return(T_perm)
                                  }
                                  if(statistics=='variance')
                                    return((apply(matrix(data_perm[,seq.int(n1)],nrow=1),1,var,na.rm=TRUE)/apply(matrix(data_perm[,n1+seq.int(n2)],nrow=1),1,var,na.rm=TRUE))^2)
                                }))
    T_perm=cbind(T_perm,T0)
  }
  
  # Not fully computable p-values (some permutations do not produce any test statistics because of the NAs)
  if(statistics=='variance'){
    no.pval=rowSums(is.na(data))>=(min(n1,n2)-1)
  }else{
    no.pval=rowSums(is.na(data))>=min(n1,n2)
  }
  #T_perm[no.pval,]=NaN # do not compute any p-value when it is not fully computable
  # do not compute any p-value if NaN is in T_perm
  #if(statistics!='variance'){
  #  pval=rowSums(T_perm>=T0)/B
  #}else{
  #  pval=pmin(2*rowSums(T_perm>=T0)/B,2*rowSums(T_perm<=T0)/B)
  #}
  # compute p-value omitting NaN
  if(statistics!='variance'){
    pval=rowSums(T_perm>=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm))
  }else{
    pval=pmin(2*rowSums(T_perm>=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm)),2*rowSums(T_perm<=T0,na.rm=TRUE)/rowSums(!is.nan(T_perm)))
  }
  
  result$T0_plot=T0_plot
  result$unadjusted_pval=pval
  result$exact=exact
  class(result)='ITWomics.2pop'
  return(list(result=result,no.pval=no.pval))
}

#setwd("C:/Users/MarziaAngela/Google Drive/L1/IWTomics analysis")
setwd("~/Downloads/IWT_lowres/")
require(IWTomics)


source('IWTomicsData_low_resolution.r') # function to create IWTomicsData object with low resolution features

# files with datasets and features
datasets=read.table("~/Downloads/IWT_lowres/datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasets=read.table("~/Downloads/IWT_lowres/features_datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)

#datasets=read.table("C:/Users/MarziaAngela/Google Drive/L1/datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
#features_datasets=read.table("C:/Users/MarziaAngela/Google Drive/L1/features_datasets.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)

# select low-resoution features
features_datasets=features_datasets[56:62,]

regionsFeatures_new=IWTomicsData(datasets$file,features_datasets[,datasets$id],'center',
                                 datasets$id,datasets$name,features_datasets$id,features_datasets$name,
                                 path='~/Downloads/IWT_lowres/files')
rm(list=setdiff(ls(),c('regionsFeatures_new','test_scalar')))
require(IWTomics)
regionsFeatures=IWTomicsData(regions(regionsFeatures_new),features(regionsFeatures_new),alignment='center',
                             id_regions=idRegions(regionsFeatures_new),name_regions=nameRegions(regionsFeatures_new),
                             id_features=idFeatures(regionsFeatures_new),name_features=nameFeatures(regionsFeatures_new),
                             length_features=lengthFeatures(regionsFeatures_new))
regionsFeatures@metadata=regionsFeatures_new@metadata
save(regionsFeatures,file='L1_low_resolution_complete.RData')

# select only autosomes
load('L1_low_resolution_complete.RData')
index=lapply(regionsFeatures@regions,function(region) seqnames(region)!='chrX')
regionsFeatures@metadata$region_datasets$size=unlist(lapply(index,sum))
regionsFeatures@regions=GRangesList(mapply(function(region,ind) region[ind,],regionsFeatures@regions,index,SIMPLIFY=FALSE))
regionsFeatures@features=lapply(regionsFeatures@features,function(feature) mapply(function(feat,ind) t(as.matrix(feat[,which(ind)])),feature,index,SIMPLIFY=FALSE))
regionsFeatures@length_features=lapply(regionsFeatures@length_features,function(feature) mapply(function(feat,ind) feat[which(ind)],feature,index,SIMPLIFY=FALSE))
validObject(regionsFeatures)
# number of windows in each dataset
lengthRegions(regionsFeatures)
save(regionsFeatures,file='L1_low_resolution_autosomes.RData')

# # select only chrX
# load('L1_low_resolution_complete.RData')
# index=lapply(regionsFeatures@regions,function(region) seqnames(region)=='chrX')
# regionsFeatures@metadata$region_datasets$size=unlist(lapply(index,sum))
# regionsFeatures@regions=GRangesList(mapply(function(region,ind) region[ind,],regionsFeatures@regions,index,SIMPLIFY=FALSE))
# regionsFeatures@features=lapply(regionsFeatures@features,function(feature) mapply(function(feat,ind) t(as.matrix(feat[,which(ind)])),feature,index,SIMPLIFY=FALSE))
# regionsFeatures@length_features=lapply(regionsFeatures@length_features,function(feature) mapply(function(feat,ind) feat[which(ind)],feature,index,SIMPLIFY=FALSE))
# validObject(regionsFeatures)
# # number of windows in each dataset
# lengthRegions(regionsFeatures)
# save(regionsFeatures,file='L1_low_resolution_chrX.RData')













######################
### ONLY AUTOSOMES ###
######################

load('L1_low_resolution_autosomes.RData')
idFeatures_select=c("Dis_Telo","Dis_Cent","Rep_Timing","Average_Recom_Rate","Telomere_Hexamer")
regionsFeatures=regionsFeatures[,idFeatures_select]
regionsFeatures@features=lapply(regionsFeatures@features,function(feat) lapply(feat,function(fe) t(fe)))

pdf('boxplot_low_resolution.pdf',width=10,height=7)
for(idFeature in idFeatures(regionsFeatures)){
  means=unlist(lapply(features(regionsFeatures)[[idFeature]],mean,na.rm=TRUE))
  boxplot(features(regionsFeatures)[[idFeature]],outline=FALSE,range=0.0000001,
          col=c('red','blue','green','gray'),names=nameRegions(regionsFeatures),main=nameFeatures(regionsFeatures)[idFeature])
  points(1:4,means)
}
dev.off()





# clustering using spearman correlation
clustering_spearman <- function(data,file_name,h,labels){
  correlation <- cor(data,method="spearman",use="pairwise.complete.obs");
  abs.corr=as.dist(1-abs(correlation))
  fit <- hclust(abs.corr)
  pdf(paste0("clustering_spearman",file_name,".pdf"),width=12,height=12)
  par(mar=c(4,1,1,14)+0.1)
  library("dendextend")
  dend=as.dendrogram(fit,hang=0.05)
  dend_cut=cutree(dend,h=h)[order.dendrogram(dend)]
  dend_col=rep(1,length(dend_cut))
  cluster=which(table(dend_cut)>1)
  for(k in seq_along(cluster)){
    dend_col[dend_cut==cluster[k]]=k+2
  }
  labels_colors(dend)=dend_col
  labels(dend)=labels[labels(dend)]
  plot(dend,xlab= "1-|Spearman's correlation|", ylab="",main=NULL,cex=0.7,cex.lab=1.5,horiz=TRUE,xlim=c(1,0))
  lines(c(h,h),c(-2,nrow(correlation)+1),col='red',lty='dashed',lwd=2)
  dev.off()
}

# all
features_plot=matrix(unlist(lapply(regionsFeatures@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_spearman(features_plot,"_low_resolution_all",h=0.2,labels=nameFeatures(regionsFeatures))
# de novo
features_plot=matrix(unlist(lapply(regionsFeatures["L1denovo",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_spearman(features_plot,"_low_resolution_denovo",h=0.2,labels=nameFeatures(regionsFeatures))
# polymorphic
features_plot=matrix(unlist(lapply(regionsFeatures["L1Pol",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_spearman(features_plot,"_low_resolution_pol",h=0.2,labels=nameFeatures(regionsFeatures))
# human specific
features_plot=matrix(unlist(lapply(regionsFeatures["L1HS",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_spearman(features_plot,"_low_resolution_hs",h=0.2,labels=nameFeatures(regionsFeatures))
# control
features_plot=matrix(unlist(lapply(regionsFeatures["Control",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_spearman(features_plot,"_low_resolution_controls",h=0.2,labels=nameFeatures(regionsFeatures))


# clustering using pearson correlation
clustering_pearson <- function(data,file_name,h,labels){
  correlation <- cor(data,method="pearson",use="pairwise.complete.obs");
  abs.corr=as.dist(1-abs(correlation))
  fit <- hclust(abs.corr)
  pdf(paste0("clustering_pearson",file_name,".pdf"),width=12,height=12)
  par(mar=c(4,1,1,14)+0.1)
  library("dendextend")
  dend=as.dendrogram(fit,hang=0.05)
  dend_cut=cutree(dend,h=h)[order.dendrogram(dend)]
  dend_col=rep(1,length(dend_cut))
  cluster=which(table(dend_cut)>1)
  for(k in seq_along(cluster)){
    dend_col[dend_cut==cluster[k]]=k+2
  }
  labels_colors(dend)=dend_col
  labels(dend)=labels[labels(dend)]
  plot(dend,xlab= "1-|Pearson's correlation|", ylab="",main=NULL,cex=0.7,cex.lab=1.5,horiz=TRUE,xlim=c(1,0))
  lines(c(h,h),c(-2,nrow(correlation)+1),col='red',lty='dashed',lwd=2)
  dev.off()
}

# all
features_plot=matrix(unlist(lapply(regionsFeatures@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_pearson(features_plot,"_low_resolution_all",h=0.2,labels=nameFeatures(regionsFeatures))
# de novo
features_plot=matrix(unlist(lapply(regionsFeatures["L1denovo",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_pearson(features_plot,"_low_resolution_denovo",h=0.2,labels=nameFeatures(regionsFeatures))
# polymorphic
features_plot=matrix(unlist(lapply(regionsFeatures["L1Pol",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_pearson(features_plot,"_low_resolution_pol",h=0.2,labels=nameFeatures(regionsFeatures))
# human specific
features_plot=matrix(unlist(lapply(regionsFeatures["L1HS",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_pearson(features_plot,"_low_resolution_hs",h=0.2,labels=nameFeatures(regionsFeatures))
# control
features_plot=matrix(unlist(lapply(regionsFeatures["Control",]@features,unlist)),ncol=nFeatures(regionsFeatures))
colnames(features_plot)=idFeatures(regionsFeatures)
clustering_pearson(features_plot,"_low_resolution_controls",h=0.2,labels=nameFeatures(regionsFeatures))










#################### test ####################
M=matrix(NA,nrow=nRegions(regionsFeatures),ncol=nRegions(regionsFeatures))
rownames(M)=idRegions(regionsFeatures)
colnames(M)=idRegions(regionsFeatures)
M=list(M=M)

result_mean=rep(M,nFeatures(regionsFeatures))
names(result_mean)=idFeatures(regionsFeatures)
for(idFeature in idFeatures(regionsFeatures)){
  message(nameFeatures(regionsFeatures)[idFeature])
  message('De novo vs Control')
  result_mean[[idFeature]][1,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$Control,B=10000)$result$unadjusted_pval
  message(result_mean[[idFeature]][1,4])
  message('Polymorphic vs Control')
  result_mean[[idFeature]][2,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$Control,B=10000)$result$unadjusted_pval
  message(result_mean[[idFeature]][2,4])
  message('Human specific vs Control')
  result_mean[[idFeature]][3,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1HS,features(regionsFeatures)[[idFeature]]$Control,B=10000)$result$unadjusted_pval
  message(result_mean[[idFeature]][3,4])
  message('De novo vs Polymorphic')
  result_mean[[idFeature]][1,2]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1Pol,B=10000)$result$unadjusted_pval
  message(result_mean[[idFeature]][1,2])
  message('De novo vs Human specific')
  result_mean[[idFeature]][1,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1HS,B=10000)$result$unadjusted_pval
  message(result_mean[[idFeature]][1,3])
  message('Polymorphic vs Human specific')
  result_mean[[idFeature]][2,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$L1HS,B=10000)$result$unadjusted_pval
  message(result_mean[[idFeature]][2,3])
}
save(result_mean,file='L1_low_resolution_results_mean.RData')

result_median=rep(M,nFeatures(regionsFeatures))
names(result_median)=idFeatures(regionsFeatures)
for(idFeature in idFeatures(regionsFeatures)){
  message(nameFeatures(regionsFeatures)[idFeature])
  message('De novo vs Control')
  result_median[[idFeature]][1,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$Control,statistics="median",B=10000)$result$unadjusted_pval
  message(result_median[[idFeature]][1,4])
  message('Polymorphic vs Control')
  result_median[[idFeature]][2,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$Control,statistics="median",B=10000)$result$unadjusted_pval
  message(result_median[[idFeature]][2,4])
  message('Human specific vs Control')
  result_median[[idFeature]][3,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1HS,features(regionsFeatures)[[idFeature]]$Control,statistics="median",B=10000)$result$unadjusted_pval
  message(result_median[[idFeature]][3,4])
  message('De novo vs Polymorphic')
  result_median[[idFeature]][1,2]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1Pol,statistics="median",B=10000)$result$unadjusted_pval
  message(result_median[[idFeature]][1,2])
  message('De novo vs Human specific')
  result_median[[idFeature]][1,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1HS,statistics="median",B=10000)$result$unadjusted_pval
  message(result_median[[idFeature]][1,3])
  message('Polymorphic vs Human specific')
  result_median[[idFeature]][2,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$L1HS,statistics="median",B=10000)$result$unadjusted_pval
  message(result_median[[idFeature]][2,3])
}
save(result_median,file='L1_low_resolution_results_median.RData')

result_variance=rep(M,nFeatures(regionsFeatures))
names(result_variance)=idFeatures(regionsFeatures)
for(idFeature in idFeatures(regionsFeatures)){
  message(nameFeatures(regionsFeatures)[idFeature])
  message('De novo vs Control')
  result_variance[[idFeature]][1,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$Control,statistics="variance",B=10000)$result$unadjusted_pval
  message(result_variance[[idFeature]][1,4])
  message('Polymorphic vs Control')
  result_variance[[idFeature]][2,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$Control,statistics="variance",B=10000)$result$unadjusted_pval
  message(result_variance[[idFeature]][2,4])
  message('Human specific vs Control')
  result_variance[[idFeature]][3,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1HS,features(regionsFeatures)[[idFeature]]$Control,statistics="variance",B=10000)$result$unadjusted_pval
  message(result_variance[[idFeature]][3,4])
  message('De novo vs Polymorphic')
  result_variance[[idFeature]][1,2]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1Pol,statistics="variance",B=10000)$result$unadjusted_pval
  message(result_variance[[idFeature]][1,2])
  message('De novo vs Human specific')
  result_variance[[idFeature]][1,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1HS,statistics="variance",B=10000)$result$unadjusted_pval
  message(result_variance[[idFeature]][1,3])
  message('Polymorphic vs Human specific')
  result_variance[[idFeature]][2,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$L1HS,statistics="variance",B=10000)$result$unadjusted_pval
  message(result_variance[[idFeature]][2,3])
}
save(result_variance,file='L1_low_resolution_results_variance.RData')

result_quantile=rep(M,nFeatures(regionsFeatures))
names(result_quantile)=idFeatures(regionsFeatures)
for(idFeature in idFeatures(regionsFeatures)){
  message(nameFeatures(regionsFeatures)[idFeature])
  message('De novo vs Control')
  result_quantile[[idFeature]][1,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$Control,statistics="quantile",probs=0.9,B=10000)$result$unadjusted_pval
  message(result_quantile[[idFeature]][1,4])
  message('Polymorphic vs Control')
  result_quantile[[idFeature]][2,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$Control,statistics="quantile",probs=0.9,B=10000)$result$unadjusted_pval
  message(result_quantile[[idFeature]][2,4])
  message('Human specific vs Control')
  result_quantile[[idFeature]][3,4]=test_scalar(features(regionsFeatures)[[idFeature]]$L1HS,features(regionsFeatures)[[idFeature]]$Control,statistics="quantile",probs=0.9,B=10000)$result$unadjusted_pval
  message(result_quantile[[idFeature]][3,4])
  message('De novo vs Polymorphic')
  result_quantile[[idFeature]][1,2]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1Pol,statistics="quantile",probs=0.9,B=10000)$result$unadjusted_pval
  message(result_quantile[[idFeature]][1,2])
  message('De novo vs Human specific')
  result_quantile[[idFeature]][1,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1denovo,features(regionsFeatures)[[idFeature]]$L1HS,statistics="quantile",probs=0.9,B=10000)$result$unadjusted_pval
  message(result_quantile[[idFeature]][1,3])
  message('Polymorphic vs Human specific')
  result_quantile[[idFeature]][2,3]=test_scalar(features(regionsFeatures)[[idFeature]]$L1Pol,features(regionsFeatures)[[idFeature]]$L1HS,statistics="quantile",probs=0.9,B=10000)$result$unadjusted_pval
  message(result_quantile[[idFeature]][2,3])
}
save(result_quantile,file='L1_low_resolution_results_quantile_90.RData')





