## Test fda.usc
#install.packages('fda.usc')
require(fda.usc)
# Load the transformed data
# source("https://bioconductor.org/biocLite.R")
# biocLite("IWTomics")
require(IWTomics)

setwd('~/Desktop/regression/10_randoms/')

# Select functional predictors by the localization table created based on the bosplots & IWT results
localization_all<-read.table("localization_table.txt",header=TRUE)
localization_all$test1

func_1<-as.vector(localization_all$test1)
func_2<-as.vector(localization_all$test2)
func_3<-as.vector(localization_all$test3)
func_4<-as.vector(localization_all$test4)
func_5<-as.vector(localization_all$test5)
func_6<-as.vector(localization_all$test6)


# y=1
# n=0
# func_1<-as.vector(c(y,y,n,y,n,y,y,n,y,n,y,n,y,y,y,n,n,n,y,n,n,y,n,y,y,y,n,n,n,n,y,n,y,y,n,y,y,y,n,y,n,n,n,n,y))
# func_2<-as.vector(c(y,n,n,y,y,n,n,n,n,y,n,n,n,n,y,n,n,y,y,y,n,n,n,n,n,n,n,n,n,n,n,y,n,n,n,y,n,n,n,n,n,n,n,n,n))
# func_3<-as.vector(c(y,y,y,y,y,n,n,n,y,y,n,n,y,y,n,n,y,n,y,n,y,y,n,n,n,y,n,n,n,n,y,y,y,n,n,y,y,y,n,n,n,n,y,n,n))
# func_4<-as.vector(c(y,y,n,n,y,y,y,n,y,y,y,y,y,y,y,y,n,n,n,y,y,y,n,y,y,y,n,n,y,n,y,n,n,y,y,y,y,n,n,n,n,n,n,n,y))
# func_5<-as.vector(c(y,y,y,y,y,y,y,y,y,y,y,y,y,y,y,y,y,n,y,n,y,y,n,y,y,y,n,n,y,y,y,y,y,y,y,y,y,y,n,n,n,n,n,n,y))
# func_6<-as.vector(c(y,y,y,y,y,n,n,n,y,y,n,n,y,y,n,n,y,n,y,n,y,n,n,n,n,y,n,n,y,y,n,n,n,n,n,n,y,y,n,n,y,n,n,n,n))

####Here for comparisons 4,5 and 6, we also changed the comparisons' directions to L1pol vs denovo, L1hs vs denovo, L1pol vs L1hs
for (r in 1:10){
  #r=1
  setwd(paste("~/Desktop/regression/10_randoms/random",r,sep=''))
  load(file = paste("L1_transformed_random_",r,'.RData',sep=''))
  #result_mean
  # test<-result_mean@features$H2AFZ_signal[[1]]
  # dim(test)
  
  #Create a data table with feature names
  # feature_names<-as.vector(c('H2AFZ signal','H3K27ac_signal','H4K20me1_signal','H3K36me3_signal','H3K4me1_signal','H3K4me2_signal',
  #                            'H3K4me3_signal','H3K79me2_signal','H3K9ac_signal','H3K9me3_signal','H3K27me3_signal','CTCF_signal',
  #                            'DNase_DHS_signal','RNA_PolII','Quadruplex','A_Phased','Direct_Repeats','Inverted_Repeats','Mirror_Repeats',
  #                            'Z_DNA','Most_Conserved','Exons','Introns','GC_Content','AT_Content','Mononucleotide','Morethan1_nuc',
  #                            'DNA_Transposons','SINE_Alu','SINE_MIR','LTR_Elements','L1_Targets','LINE_L2&L3','CpG_Islands','5hMc',
  #                            'Sperm_hypometh','Rep_origin','Recombination_Hot','Exon_Expression','Gene_Expression','Tracript_Expression',
  #                            'Testis_Expression','chh_meth','chg_meth','cpg_meth'))
  
  ## Change the feature names to the same as IWT plot
  feature_names<-as.vector(result_mean@metadata$feature_datasets$name)
  
  feature_num<- as.vector(c(1:45))
  feature<-as.data.frame(matrix(ncol=1))
  for (i in 1:45){
    feature<-rbind(feature,paste("feature_", i, sep = ""))
  }
  feature<-feature[-1,]
  feature_table<-data.frame(feature_num,feature,feature_names)
  
  
  #######Logsitic regression for func predictors
  # Comparison 1 (denovo vs control): Logistic regression for func predictors
  comp1_func<-which(func_1=="y")
  comp1_func_logit_sum <- as.data.frame(matrix(ncol = 3))
  colnames(comp1_func_logit_sum)<-c("feature", "significance", "pseudoR2")
  
  comp1_func_logit_7coeff <- as.data.frame(matrix(nrow = 8))
  comp1_func_logit_7coeff<-unname(comp1_func_logit_7coeff)
  for (i in comp1_func){
    #i=comp1_func[1]
    l1<-result_mean@features[[i]][[1]]
    ct<-result_mean@features[[i]][[4]]
    dim(l1)
    dim(ct)
    # Convert  from class fda to fdata
    x<-rbind(t(l1),t(ct))
    # bsp1 <- create.bspline.basis(c(1:100),10)
    # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
    # fdataobj=fdata(fd1)
    x_f<-fdata(x)
    #feature1<-x_f$data
    y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
    dataf=as.data.frame(y)
    x=paste('feature',i,sep='')
    tdata=list("df"=dataf,'x'=x_f)
    #x_f[["argvals"]]=x_f[["argvals"]]-50.5
    tt=x_f[["argvals"]]
    # define the same basis for both x and b
    # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
    basis_common=create.bspline.basis(rangeval = range(tt),norder=3,breaks=seq(from =0.5,to=100.5, length.out=6))
    basis.x=list("x"=basis_common)
    basis.b=list("x"=basis_common)
    f=y~x
    res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
    coef(res)
    # res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
    #                basis.b=basis.b)
    #summary(res)
    pseudoR2<-1-(res$deviance/res$null.deviance)
    #coefficient<-coef(summary(res))
    #summary(res)[min(coef(res)[,4])]
    #coef(summary(res))[which.min(coef(summary(res))[,4])]
    #min(coef(summary(res))[,4])
    significance<-min(coef(summary(res))[,4])
    res_sum<-cbind(paste("feature_", i, sep = ""), significance, pseudoR2)
    colnames(res_sum)<-c("feature", "significance", "pseudoR2")
    comp1_func_logit_sum<-rbind(comp1_func_logit_sum,res_sum)
    # save the 20+1 coefficient and their pvals with feature number
    coef<- as.matrix(coef(summary(res))[,1])
    colnames(coef)<-paste('coef_',feature_table$feature_names[i],sep='')
    # colnames(coef)<-paste('coef_feature',i,sep='')
    pval<- as.matrix(coef(summary(res))[,4])
    #colnames(pval)<-paste('pval_feature',i,sep='')
    colnames(pval)<-paste('pval_',feature_table$feature_names[i],sep='')
    comp1_func_logit_7coeff<-cbind(comp1_func_logit_7coeff,coef)
    comp1_func_logit_7coeff<-cbind(comp1_func_logit_7coeff,pval)
    
    #Plot the beta function of each feature
    betafd=fd(res$coefficients[-1],basisobj=basis_common)
    pdf(paste('comp1_beta_function_',feature_table$feature_names[i],'.pdf',sep=''))
    plot(1:100,eval.fd(1:100,betafd),main=paste('Comparison 1: beta function\n',feature_table$feature_names[i],sep=''))
    dev.off()
  }
  comp1_func_logit_sum<-comp1_func_logit_sum[-1,]
  #rownames(comp1_func_logit_sum)<-c(1:length(comp1_func))
  comp1_func_logit_7coeff<-comp1_func_logit_7coeff[,-1]
  
  
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp1_func_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp1_func_logit_sum<-sum_table[order(sum_table$feature_num),]
  # write.table(comp1_func_logit_sum, file='comp1_func_logit_sum.txt',sep = '\t')
  # write.table(comp1_func_logit_7coeff, file='comp1_func_logit_coefficients.txt',sep = '\t')
  write.csv(comp1_func_logit_sum, file='comp1_func_logit_sum.csv')
  write.csv(comp1_func_logit_7coeff, file='comp1_func_logit_coefficients.csv')
  
  ###########
  # Comparison 2 (L1pol VS control): Logistic regression for func predictors
  comp2_func<-which(func_2=="y")
  comp2_func_logit_sum <- as.data.frame(matrix(ncol = 3))
  colnames(comp2_func_logit_sum)<-c("feature", "significance", "pseudoR2")
  
  comp2_func_logit_7coeff <- as.data.frame(matrix(nrow = 8))
  comp2_func_logit_7coeff<-unname(comp2_func_logit_7coeff)
  for (i in comp2_func){
    #i=comp2_func[1]
    l1<-result_mean@features[[i]][[2]]
    ct<-result_mean@features[[i]][[4]]
    dim(l1)
    dim(ct)
    # Convert  from class fda to fdata
    x<-rbind(t(l1),t(ct))
    # bsp1 <- create.bspline.basis(c(1:100),10)
    # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
    # fdataobj=fdata(fd1)
    x_f<-fdata(x)
    #feature1<-x_f$data
    y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
    dataf=as.data.frame(y)
    x=paste('feature',i,sep='')
    tdata=list("df"=dataf,'x'=x_f)
    #x_f[["argvals"]]=x_f[["argvals"]]-50.5
    tt=x_f[["argvals"]]
    # define the same basis for both x and b
    # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
    basis_common=create.bspline.basis(rangeval = range(tt),norder=3,breaks=seq(from =0.5,to=100.5, length.out=6))
    basis.x=list("x"=basis_common)
    basis.b=list("x"=basis_common)
    f=y~x
    res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
    coef(res)
    # res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
    #                basis.b=basis.b)
    #summary(res)
    pseudoR2<-1-(res$deviance/res$null.deviance)
    #coefficient<-coef(summary(res))
    #summary(res)[min(coef(res)[,4])]
    #coef(summary(res))[which.min(coef(summary(res))[,4])]
    #min(coef(summary(res))[,4])
    significance<-min(coef(summary(res))[,4])
    res_sum<-cbind(paste("feature_", i, sep = ""), significance, pseudoR2)
    colnames(res_sum)<-c("feature", "significance", "pseudoR2")
    comp2_func_logit_sum<-rbind(comp2_func_logit_sum,res_sum)
    # save the 20+1 coefficient and their pvals with feature number
    coef<- as.matrix(coef(summary(res))[,1])
    colnames(coef)<-paste('coef_',feature_table$feature_names[i],sep='')
    pval<- as.matrix(coef(summary(res))[,4])
    colnames(pval)<-paste('pval_',feature_table$feature_names[i],sep='')
    
    comp2_func_logit_7coeff<-cbind(comp2_func_logit_7coeff,coef)
    comp2_func_logit_7coeff<-cbind(comp2_func_logit_7coeff,pval)
    
    #Plot the beta function of each feature
    betafd=fd(res$coefficients[-1],basisobj=basis_common)
    pdf(paste('comp2_beta_function_',feature_table$feature_names[i],'.pdf',sep=''))
    plot(1:100,eval.fd(1:100,betafd),main=paste('Comparison 2: beta function\n',feature_table$feature_names[i],sep=''))
    dev.off()
  }
  comp2_func_logit_sum<-comp2_func_logit_sum[-1,]
  #rownames(comp2_func_logit_sum)<-c(1:length(comp2_func))
  comp2_func_logit_7coeff<-comp2_func_logit_7coeff[,-1]
  
  
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp2_func_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp2_func_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp2_func_logit_sum, file='comp2_func_logit_sum.csv')
  write.csv(comp2_func_logit_7coeff, file='comp2_func_logit_coefficients.csv')
  
  
  ##########
  # Comparison 3 ( L1HS vs Control): Logistic regression for func predictors
  comp3_func<-which(func_3=="y")
  comp3_func_logit_sum <- as.data.frame(matrix(ncol = 3))
  colnames(comp3_func_logit_sum)<-c("feature", "significance", "pseudoR2")
  
  comp3_func_logit_7coeff <- as.data.frame(matrix(nrow = 8))
  comp3_func_logit_7coeff<-unname(comp3_func_logit_7coeff)
  for (i in comp3_func){
    #i=comp3_func[1]
    l1<-result_mean@features[[i]][[3]]
    ct<-result_mean@features[[i]][[4]]
    dim(l1)
    dim(ct)
    # Convert  from class fda to fdata
    x<-rbind(t(l1),t(ct))
    # bsp1 <- create.bspline.basis(c(1:100),10)
    # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
    # fdataobj=fdata(fd1)
    x_f<-fdata(x)
    #feature1<-x_f$data
    y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
    dataf=as.data.frame(y)
    x=paste('feature',i,sep='')
    tdata=list("df"=dataf,'x'=x_f)
    #x_f[["argvals"]]=x_f[["argvals"]]-50.5
    tt=x_f[["argvals"]]
    # define the same basis for both x and b
    # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
    basis_common=create.bspline.basis(rangeval = range(tt),norder=3,breaks=seq(from =0.5,to=100.5, length.out=6))
    basis.x=list("x"=basis_common)
    basis.b=list("x"=basis_common)
    f=y~x
    res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
    coef(res)
    # res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
    #                basis.b=basis.b)
    #summary(res)
    pseudoR2<-1-(res$deviance/res$null.deviance)
    #coefficient<-coef(summary(res))
    #summary(res)[min(coef(res)[,4])]
    #coef(summary(res))[which.min(coef(summary(res))[,4])]
    #min(coef(summary(res))[,4])
    significance<-min(coef(summary(res))[,4])
    res_sum<-cbind(paste("feature_", i, sep = ""), significance, pseudoR2)
    colnames(res_sum)<-c("feature", "significance", "pseudoR2")
    comp3_func_logit_sum<-rbind(comp3_func_logit_sum,res_sum)
    # save the 20+1 coefficient and their pvals with feature number
    coef<- as.matrix(coef(summary(res))[,1])
    colnames(coef)<-paste('coef_',feature_table$feature_names[i],sep='')
    pval<- as.matrix(coef(summary(res))[,4])
    colnames(pval)<-paste('pval_',feature_table$feature_names[i],sep='')
    
    comp3_func_logit_7coeff<-cbind(comp3_func_logit_7coeff,coef)
    comp3_func_logit_7coeff<-cbind(comp3_func_logit_7coeff,pval)
    
    #Plot the beta function of each feature
    betafd=fd(res$coefficients[-1],basisobj=basis_common)
    pdf(paste('comp3_beta_function_',feature_table$feature_names[i],'.pdf',sep=''))
    plot(1:100,eval.fd(1:100,betafd),main=paste('Comparison 3: beta function\n',feature_table$feature_names[i],sep=''))
    dev.off()
  }
  comp3_func_logit_sum<-comp3_func_logit_sum[-1,]
  #rownames(comp3_func_logit_sum)<-c(1:length(comp3_func))
  comp3_func_logit_7coeff<-comp3_func_logit_7coeff[,-1]
  
  
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp3_func_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp3_func_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp3_func_logit_sum, file='comp3_func_logit_sum.csv')
  write.csv(comp3_func_logit_7coeff, file='comp3_func_logit_coefficients.csv')
  
  ##########
  # Comparison 4 ( de novo vs L1Pol): Logistic regression for func predictors
  comp4_func<-which(func_4=="y")
  comp4_func_logit_sum <- as.data.frame(matrix(ncol = 3))
  colnames(comp4_func_logit_sum)<-c("feature", "significance", "pseudoR2")
  
  comp4_func_logit_7coeff <- as.data.frame(matrix(nrow = 8))
  comp4_func_logit_7coeff<-unname(comp4_func_logit_7coeff)
  for (i in comp4_func){
    #i=comp4_func[1]
    l1<-result_mean@features[[i]][[2]]
    ct<-result_mean@features[[i]][[1]]
    dim(l1)
    dim(ct)
    # Convert  from class fda to fdata
    x<-rbind(t(l1),t(ct))
    # bsp1 <- create.bspline.basis(c(1:100),10)
    # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
    # fdataobj=fdata(fd1)
    x_f<-fdata(x)
    #feature1<-x_f$data
    y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
    dataf=as.data.frame(y)
    x=paste('feature',i,sep='')
    tdata=list("df"=dataf,'x'=x_f)
    #x_f[["argvals"]]=x_f[["argvals"]]-50.5
    tt=x_f[["argvals"]]
    # define the same basis for both x and b
    # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
    basis_common=create.bspline.basis(rangeval = range(tt),norder=3,breaks=seq(from =0.5,to=100.5, length.out=6))
    basis.x=list("x"=basis_common)
    basis.b=list("x"=basis_common)
    f=y~x
    res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
    coef(res)
    # res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
    #                basis.b=basis.b)
    #summary(res)
    pseudoR2<-1-(res$deviance/res$null.deviance)
    #coefficient<-coef(summary(res))
    #summary(res)[min(coef(res)[,4])]
    #coef(summary(res))[which.min(coef(summary(res))[,4])]
    #min(coef(summary(res))[,4])
    significance<-min(coef(summary(res))[,4])
    res_sum<-cbind(paste("feature_", i, sep = ""), significance, pseudoR2)
    colnames(res_sum)<-c("feature", "significance", "pseudoR2")
    comp4_func_logit_sum<-rbind(comp4_func_logit_sum,res_sum)
    # save the 20+1 coefficient and their pvals with feature number
    coef<- as.matrix(coef(summary(res))[,1])
    colnames(coef)<-paste('coef_',feature_table$feature_names[i],sep='')
    pval<- as.matrix(coef(summary(res))[,4])
    colnames(pval)<-paste('pval_',feature_table$feature_names[i],sep='')
    
    comp4_func_logit_7coeff<-cbind(comp4_func_logit_7coeff,coef)
    comp4_func_logit_7coeff<-cbind(comp4_func_logit_7coeff,pval)
    
    #Plot the beta function of each feature
    betafd=fd(res$coefficients[-1],basisobj=basis_common)
    pdf(paste('comp4_beta_function_',feature_table$feature_names[i],'.pdf',sep=''))
    plot(1:100,eval.fd(1:100,betafd),main=paste('Comparison 4: beta function\n',feature_table$feature_names[i],sep=''))
    dev.off()
  }
  comp4_func_logit_sum<-comp4_func_logit_sum[-1,]
  #rownames(comp4_func_logit_sum)<-c(1:length(comp4_func))
  comp4_func_logit_7coeff<-comp4_func_logit_7coeff[,-1]
  
  
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp4_func_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp4_func_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp4_func_logit_sum, file='comp4_func_logit_sum.csv')
  write.csv(comp4_func_logit_7coeff, file='comp4_func_logit_coefficients.csv')
  
  ##########
  # Comparison 5 ( de novo vs L1HS): Logistic regression for func predictors
  comp5_func<-which(func_5=="y")
  comp5_func_logit_sum <- as.data.frame(matrix(ncol = 3))
  colnames(comp5_func_logit_sum)<-c("feature", "significance", "pseudoR2")
  
  comp5_func_logit_7coeff <- as.data.frame(matrix(nrow = 8))
  comp5_func_logit_7coeff<-unname(comp5_func_logit_7coeff)
  for (i in comp5_func){
    #i=comp5_func[1]
    l1<-result_mean@features[[i]][[3]]
    ct<-result_mean@features[[i]][[1]]
    dim(l1)
    dim(ct)
    # Convert  from class fda to fdata
    x<-rbind(t(l1),t(ct))
    # bsp1 <- create.bspline.basis(c(1:100),10)
    # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
    # fdataobj=fdata(fd1)
    x_f<-fdata(x)
    #feature1<-x_f$data
    y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
    dataf=as.data.frame(y)
    x=paste('feature',i,sep='')
    tdata=list("df"=dataf,'x'=x_f)
    #x_f[["argvals"]]=x_f[["argvals"]]-50.5
    tt=x_f[["argvals"]]
    # define the same basis for both x and b
    # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
    basis_common=create.bspline.basis(rangeval = range(tt),norder=3,breaks=seq(from =0.5,to=100.5, length.out=6))
    basis.x=list("x"=basis_common)
    basis.b=list("x"=basis_common)
    f=y~x
    res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
    coef(res)
    # res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
    #                basis.b=basis.b)
    #summary(res)
    pseudoR2<-1-(res$deviance/res$null.deviance)
    #coefficient<-coef(summary(res))
    #summary(res)[min(coef(res)[,4])]
    #coef(summary(res))[which.min(coef(summary(res))[,4])]
    #min(coef(summary(res))[,4])
    significance<-min(coef(summary(res))[,4])
    res_sum<-cbind(paste("feature_", i, sep = ""), significance, pseudoR2)
    colnames(res_sum)<-c("feature", "significance", "pseudoR2")
    comp5_func_logit_sum<-rbind(comp5_func_logit_sum,res_sum)
    # save the 20+1 coefficient and their pvals with feature number
    coef<- as.matrix(coef(summary(res))[,1])
    colnames(coef)<-paste('coef_',feature_table$feature_names[i],sep='')
    pval<- as.matrix(coef(summary(res))[,4])
    colnames(pval)<-paste('pval_',feature_table$feature_names[i],sep='')
    
    comp5_func_logit_7coeff<-cbind(comp5_func_logit_7coeff,coef)
    comp5_func_logit_7coeff<-cbind(comp5_func_logit_7coeff,pval)
    
    #Plot the beta function of each feature
    betafd=fd(res$coefficients[-1],basisobj=basis_common)
    pdf(paste('comp5_beta_function_',feature_table$feature_names[i],'.pdf',sep=''))
    plot(1:100,eval.fd(1:100,betafd),main=paste('Comparison 5: beta function\n',feature_table$feature_names[i],sep=''))
    dev.off()
  }
  comp5_func_logit_sum<-comp5_func_logit_sum[-1,]
  #rownames(comp5_func_logit_sum)<-c(1:length(comp5_func))
  comp5_func_logit_7coeff<-comp5_func_logit_7coeff[,-1]
  
  
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp5_func_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp5_func_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp5_func_logit_sum, file='comp5_func_logit_sum.csv')
  write.csv(comp5_func_logit_7coeff, file='comp5_func_logit_coefficients.csv')
  
  ##########
  # Comparison 6 ( L1HS vs L1Pol): Logistic regression for func predictors
  comp6_func<-which(func_6=="y")
  comp6_func_logit_sum <- as.data.frame(matrix(ncol = 3))
  colnames(comp6_func_logit_sum)<-c("feature", "significance", "pseudoR2")
  
  comp6_func_logit_7coeff <- as.data.frame(matrix(nrow = 8))
  comp6_func_logit_7coeff<-unname(comp6_func_logit_7coeff)
  for (i in comp6_func){
    #i=comp6_func[1]
    l1<-result_mean@features[[i]][[3]]
    ct<-result_mean@features[[i]][[2]]
    dim(l1)
    dim(ct)
    # Convert  from class fda to fdata
    x<-rbind(t(l1),t(ct))
    # bsp1 <- create.bspline.basis(c(1:100),10)
    # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
    # fdataobj=fdata(fd1)
    x_f<-fdata(x)
    #feature1<-x_f$data
    y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
    dataf=as.data.frame(y)
    x=paste('feature',i,sep='')
    tdata=list("df"=dataf,'x'=x_f)
    #x_f[["argvals"]]=x_f[["argvals"]]-50.5
    tt=x_f[["argvals"]]
    # define the same basis for both x and b
    # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
    basis_common=create.bspline.basis(rangeval = range(tt),norder=3,breaks=seq(from =0.5,to=100.5, length.out=6))
    basis.x=list("x"=basis_common)
    basis.b=list("x"=basis_common)
    f=y~x
    res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
    coef(res)
    # res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
    #                basis.b=basis.b)
    #summary(res)
    pseudoR2<-1-(res$deviance/res$null.deviance)
    #coefficient<-coef(summary(res))
    #summary(res)[min(coef(res)[,4])]
    #coef(summary(res))[which.min(coef(summary(res))[,4])]
    #min(coef(summary(res))[,4])
    significance<-min(coef(summary(res))[,4])
    res_sum<-cbind(paste("feature_", i, sep = ""), significance, pseudoR2)
    colnames(res_sum)<-c("feature", "significance", "pseudoR2")
    comp6_func_logit_sum<-rbind(comp6_func_logit_sum,res_sum)
    # save the 20+1 coefficient and their pvals with feature number
    coef<- as.matrix(coef(summary(res))[,1])
    colnames(coef)<-paste('coef_',feature_table$feature_names[i],sep='')
    pval<- as.matrix(coef(summary(res))[,4])
    colnames(pval)<-paste('pval_',feature_table$feature_names[i],sep='')
    
    comp6_func_logit_7coeff<-cbind(comp6_func_logit_7coeff,coef)
    comp6_func_logit_7coeff<-cbind(comp6_func_logit_7coeff,pval)
    
    #Plot the beta function of each feature
    betafd=fd(res$coefficients[-1],basisobj=basis_common)
    pdf(paste('comp6_beta_function_',feature_table$feature_names[i],'.pdf',sep=''))
    plot(1:100,eval.fd(1:100,betafd),main=paste('Comparison 6: beta function\n',feature_table$feature_names[i],sep=''))
    dev.off()
  }
  comp6_func_logit_sum<-comp6_func_logit_sum[-1,]
  #rownames(comp6_func_logit_sum)<-c(1:length(comp6_func))
  comp6_func_logit_7coeff<-comp6_func_logit_7coeff[,-1]
  
  
  
  #Add names of the features in the table 
  sum_table<-merge(feature_table, comp6_func_logit_sum,by.x=c("feature"),by.y=c("feature"))
  comp6_func_logit_sum<-sum_table[order(sum_table$feature_num),]
  write.csv(comp6_func_logit_sum, file='comp6_func_logit_sum.csv')
  write.csv(comp6_func_logit_7coeff, file='comp6_func_logit_coefficients.csv')
  
  save.image(paste('functional_individual_regression_random',r,'.RData',sep=''))
  
}




###DONE######





# 
# 
# 
# ###############
# ## Check regression results individually for the random samples (1 & 6)
# r=1
# setwd(paste("~/Desktop/regression/10_randoms/random",r,sep=''))
# load(file = paste("L1_transformed_random_",r,'.RData',sep=''))
# ## comparison 1
# comp1_func
# comp1_func_logit_sum
# #H3K4me3(activation)
# i=comp1_func[5]
# #CpG methylation
# #i=comp1_func[23]
# l1<-result_mean@features[[i]][[1]]
# ct<-result_mean@features[[i]][[4]]
# dim(l1)
# dim(ct)
# # Convert  from class fda to fdata
# x<-rbind(t(l1),t(ct))
# # bsp1 <- create.bspline.basis(c(1:100),10)
# # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
# # fdataobj=fdata(fd1)
# x_f<-fdata(x)
# #feature1<-x_f$data
# y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
# dataf=as.data.frame(y)
# x=paste('feature',i,sep='')
# tdata=list("df"=dataf,'x'=x_f)
# #x_f[["argvals"]]=x_f[["argvals"]]-50.5
# tt=x_f[["argvals"]]
# # define the same basis for both x and b
# # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
# basis_common=create.bspline.basis(rangeval = range(tt),norder=1,breaks=seq(from =0.5,to=100.5, length.out=21))
# basis.x=list("x"=basis_common)
# basis.b=list("x"=basis_common)
# f=y~x
# res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
# summary(res)
# 
# 
# ###Coparison 3
# #comp3_func
# #H3K4me3(activation)
# i=comp3_func[5]
# #CpG methylation
# #i=comp3_func[23]
# l1<-result_mean@features[[i]][[3]]
# ct<-result_mean@features[[i]][[4]]
# dim(l1)
# dim(ct)
# # Convert  from class fda to fdata
# x<-rbind(t(l1),t(ct))
# # bsp1 <- create.bspline.basis(c(1:100),10)
# # fd1 <- Data2fd(1:100,y=t(y),basisobj=bsp1)
# # fdataobj=fdata(fd1)
# x_f<-fdata(x)
# #feature1<-x_f$data
# y<-c(rep(1,length(l1[1,])),rep(0,length(ct[1,])))
# dataf=as.data.frame(y)
# x=paste('feature',i,sep='')
# tdata=list("df"=dataf,'x'=x_f)
# #x_f[["argvals"]]=x_f[["argvals"]]-50.5
# tt=x_f[["argvals"]]
# # define the same basis for both x and b
# # Instead of defining nbasis, here to define the number of nodes: norder=1 & breaks=21 
# basis_common=create.bspline.basis(rangeval = range(tt),norder=1,breaks=seq(from =0.5,to=100.5, length.out=21))
# basis.x=list("x"=basis_common)
# basis.b=list("x"=basis_common)
# f=y~x
# res=fregre.glm(f,family=binomial(link=logit),data=tdata,basis.x=basis.x,basis.b=basis.b)
# summary(res)
# # 
# # 
# # 
# #Number of functional predictors where significances were not captured by regression (by the cut-off: smallest pval=0.05)
# length(as.vector(which(as.vector(comp1_func_logit_sum$significance)>0.1)))
# length(as.vector(which(as.vector(comp2_func_logit_sum$significance)>0.1)))
# length(as.vector(which(as.vector(comp3_func_logit_sum$significance)>0.1)))
# length(as.vector(which(as.vector(comp4_func_logit_sum$significance)>0.1)))
# length(as.vector(which(as.vector(comp5_func_logit_sum$significance)>0.1)))
# length(as.vector(which(as.vector(comp6_func_logit_sum$significance)>0.5)))
# # 
# # length(as.vector(which(as.vector(comp5_func_logit_sum$significance)<0.05)))


# ## Create functional data object
# data(phoneme)
# mlearn<-phoneme$learn[1:4,1:150]
# # Center curves
# fdata.c=fdata.cen(mlearn)$Xcen
# par(mfrow=c(2,1))
# plot.fdata(mlearn,type="l")
# plot.fdata(fdata.c,type="l")
# 
# # Convert  from class fda to fdata
# bsp1 <- create.bspline.basis(c(1,150),21)
# fd1 <- Data2fd(1:150,y=t(mlearn$data),basisobj=bsp1)
# fdataobj=fdata(fd1)
# 
# # Convert  from class fds, fts or sfts to fdata
# #require(fds)
# #a=fds(x = 1:20, y = Simulationdata$y, xname = "x", 
# # yname = "Simulated value")
# #b=fts(x = 15:49, y = Australiasmoothfertility$y, xname = "Age",
# #    yname = "Fertility rate")
# #c=sfts(ts(as.numeric(ElNino$y), frequency = 12), xname = "Month",
# #yname = "Sea surface temperature")
# #class(a);class(b);class(c)
# #fdataobj=fdata(b)
# 
# data(tecator)
# x=tecator$absorp.fdata
# y=tecator$y$Fat
# tt=x[["argvals"]]
# dataf=as.data.frame(tecator$y)
# 
# nbasis.x=11 # define the same basis for both x and b
# nbasis.b=7
# basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)# Instead, define the number of nodes instead, order 1, 21 knots 
# basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)
# 
# f=Fat~Protein+x
# basis.x=list("x"=basis1)
# basis.b=list("x"=basis2)
# ldata=list("df"=dataf,"x"=x)
# res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,
#                basis.b=basis.b)
# summary(res)

# Load the transformed data