require(IWTomics)
require(dendextend)

setwd("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/")
#setwd("/Users/Bruce/Desktop/L1/IWTomics analysis_1000control_1000denovo")

# Read files with datasets and features
datasets=read.table("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/datasets_new_revision.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasets=read.table("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/features_datasets_new_revision.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# Exclude low-resoution features (analyzed in IWT_L1_low_resolution.r)
features_datasets=features_datasets[1:13,]
tail(features_datasets)

# Load data
regionsFeatures=IWTomicsData(datasets$file,features_datasets[,datasets$id],'center',
                             datasets$id,datasets$name,features_datasets$id,features_datasets$name,
                             path='/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/')
save(regionsFeatures,file='L1_complete_new_revision.RData')



# select only autosomes
load('L1_complete_new_revision.RData')
index=lapply(regionsFeatures@regions,function(region) seqnames(region)!='chrX')
regionsFeatures@metadata$region_datasets$size=unlist(lapply(index,sum))
regionsFeatures@regions=GRangesList(mapply(function(region,ind) region[ind,],regionsFeatures@regions,index,SIMPLIFY=FALSE))
regionsFeatures@features=lapply(regionsFeatures@features,function(feature) mapply(function(feat,ind) feat[,which(ind)],feature,index,SIMPLIFY=FALSE))
regionsFeatures@length_features=lapply(regionsFeatures@length_features,function(feature) mapply(function(feat,ind) feat[which(ind)],feature,index,SIMPLIFY=FALSE))
validObject(regionsFeatures)
# number of windows in each dataset
lengthRegions(regionsFeatures)
save(regionsFeatures,file='L1_autosomes_new_revision.RData')


######################
### ONLY AUTOSOMES ###
######################

load('L1_autosomes_new_revision.RData')
# select only feature we want to analyze
idFeatures_select=c("H2AFZ_signal","H3K27ac_signal","H4K20me1_signal","H3K36me3_signal","H3K4me1_signal",
                    "H3K4me2_signal","H3K4me3_signal","H3K79me2_signal","H3K9ac_signal","H3K9me3_signal",
                    "H3K27me3_signal","CTCF_signal","DNase_DHS_signal")#,"cpg_meth")
regionsFeatures=regionsFeatures[,idFeatures_select]

# number of 0 in the different features
# require that at least 10% of data are not zeros... 
# if this is not the case, smooth them with kernel smoothing and increasing bandwidth
zero_count=lapply(regionsFeatures@features,
                  function(feature){
                    counts=Reduce(rbind,lapply(feature,function(feat) c(sum(feat==0),sum(feat!=0),length(as.vector(feat)),length(unique(as.vector(feat))))))
                    colnames(counts)=c('0','>0','tot','distinct')
                    rownames(counts)=names(feature)
                    return(counts)
                  })
zero_count
zero_count_tot=Reduce(rbind,lapply(regionsFeatures@features,
                                   function(feature){
                                     feat=unlist(feature)
                                     count=c(sum(feat==0),sum(feat!=0),length(feat),length(unique(feat)))
                                     names(count)=c('0','>0','tot','distinct')
                                     return(count)
                                   }))
rownames(zero_count_tot)=names(zero_count)
zero_count_tot

many_zeros=names(which((zero_count_tot[,2]/zero_count_tot[,3]*100)<10))
many_zeros
##No actions needed to fix zereos
save(regionsFeatures,file='L1_autosomes_new_revision_zero_fixed.RData')




##################################################PLOT 
load('L1_autosomes_new_revision_zero_fixed.RData')
# plot
pdf('curves_autosomes.pdf',width=10,height=8)
plot(regionsFeatures,type='curves',id_regions_subset=c('L1denovoFlasch','L1denovoSultana','Control','L1denovo','L1Pol','L1HS'),
     col=c('orange','purple','black','red','blue','green'),xlab='kb',ask=FALSE)
dev.off()

pdf('boxplot_autosomes.pdf',width=10,height=7)
plot(regionsFeatures,type='boxplot',id_regions_subset=c('L1denovoFlasch','L1denovoSultana','Control','L1denovo','L1Pol','L1HS'),
     col=c('orange','purple','black','red','blue','green'),xlab='kb',ask=FALSE)
dev.off()

###########
# include strand information and annotations
# (reverse the 100-kb region measurements for minus strand)
load('L1_autosomes_new_revision_zero_fixed.RData')

# de novo
L1_flanking_100=read.table('/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/L1denovo_100kb_no_gaps_no_blacklisted_method2.interval',header=FALSE,sep='\t')
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1denovo,L1_flanking_100,ignore.strand=TRUE)]
L1denovo_new=L1_flanking_100

# polymorphic
L1_flanking_100=read.table('/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/L1Pol_100kb_no_gaps_no_blacklisted_method2.interval',header=FALSE,sep='\t')
L1_flanking_100=cbind(L1_flanking_100[,1:4],NA,NA,L1_flanking_100[,5:7])
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1Pol,L1_flanking_100,ignore.strand=TRUE)]
L1Pol_new=L1_flanking_100

# human specific
L1_flanking_100=read.table('/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/L1HS_100kb_no_gaps_no_blacklisted_method2_withL1seq.interval',header=FALSE,sep='\t')
L1_flanking_100=cbind(L1_flanking_100[,1:4],NA,NA,L1_flanking_100[,5:7])
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1HS,L1_flanking_100,ignore.strand=TRUE)]
L1HS_new=L1_flanking_100

Control_new=regionsFeatures@regions$Control
mcols(Control_new)=data.frame(barcode=NA,site=NA,annot=NA,overlap_same=0,overlap_other=0)

# de novo Flasch
L1_flanking_100=read.table('/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/Flasch_denovo_hESC_100k_no_gaps_no_blacklisted_method2.interval',header=FALSE,sep='\t')
L1_flanking_100=cbind(L1_flanking_100[,1:4],NA,NA,NA,L1_flanking_100[,5],"no")
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1denovoFlasch,L1_flanking_100,ignore.strand=TRUE)]
L1denovo_new_flasch=L1_flanking_100

# de novo Sultana
L1_flanking_100=read.table('/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/Sultana_denovo_Hela_100k_no_gaps_no_blacklisted_method2.interval',header=FALSE,sep='\t')
L1_flanking_100=cbind(L1_flanking_100[,1:4],NA,NA,NA,L1_flanking_100[,5],"no")
names(L1_flanking_100)=c('chr','start','end','strand','barcode','site','annot','overlap_same','overlap_other')
L1_flanking_100=makeGRangesFromDataFrame(L1_flanking_100,starts.in.df.are.0based=TRUE,keep.extra.columns=TRUE)
L1_flanking_100=L1_flanking_100[match(regionsFeatures@regions$L1denovoSultana,L1_flanking_100,ignore.strand=TRUE)]
L1denovo_new_sultana=L1_flanking_100


# reverse the 100-kb region measurements for minus strand
regionsFeatures@regions=GRangesList(L1denovoFlasch=L1denovo_new_flasch,L1denovoSultana=L1denovo_new_sultana, L1denovo=L1denovo_new,L1Pol=L1Pol_new,L1HS=L1HS_new,Control=Control_new)
regionsFeatures@features=lapply(regionsFeatures@features,
                                function(feature){
                                  for(region in idRegions(regionsFeatures)){
                                    minus=which(strand(regions(regionsFeatures)[[region]])=='-')
                                    feature[[region]][,minus]=feature[[region]][seq(nrow(feature[[region]]),1),minus]
                                  }
                                  return(feature)
                                })
validObject(regionsFeatures)
save(regionsFeatures,file='L1_autosomes_zero_fixed_with_strand_new_revision.RData')



load('L1_autosomes_zero_fixed_with_strand_new_revision.RData')
# plot
# pdf('curves_autosomes_with_strand.pdf',width=10,height=8)
# plot(regionsFeatures,type='curves',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#      col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()

##Visulazation with boxplot
pdf('boxplot_autosomes_with_strand.pdf',width=10,height=7)
plot(regionsFeatures,type='boxplot',id_regions_subset=c('L1denovoFlasch','L1denovoSultana','Control','L1denovo','L1Pol','L1HS'),
     col=c('orange','purple','black','red','blue','green'),xlab='kb',ask=FALSE)
dev.off()

# pdf('boxplot_autosomes_with_strand_denovo.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1denovo'),
#      col=c('black','red'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_pol.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1Pol'),
#      col=c('black','blue'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_hs.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1HS'),
#      col=c('black','green'),xlab='kb',ask=FALSE)
# dev.off()

# smooth all the features
load('L1_autosomes_zero_fixed_with_strand_new_revision.RData')
already_smoothed=c("")
# smooth measurements... and keep 100 windows
band=2
regionsFeatures_smoothed=smooth(regionsFeatures,id_features_subset=setdiff(idFeatures(regionsFeatures),already_smoothed),type='kernel',bandwidth=band)
save(regionsFeatures_smoothed,file='L1_autosomes_zero_fixed_with_strand_smoothed_new_revision.RData')


load('L1_autosomes_zero_fixed_with_strand_smoothed_new_revision.RData')
# plot
# pdf('curves_autosomes_with_strand_smoothed.pdf',width=10,height=8)
# plot(regionsFeatures_smoothed,type='curves',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#      col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()

pdf('boxplot_autosomes_with_strand_smoothed_new_revision.pdf',width=10,height=7)
plot(regionsFeatures,type='boxplot',id_regions_subset=c('L1denovoFlasch','L1denovoSultana','Control','L1denovo','L1Pol','L1HS'),
     col=c('orange','purple','black','red','blue','green'),xlab='kb',ask=FALSE)
dev.off()

# pdf('boxplot_autosomes_with_strand_smoothed_denovo.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1denovo'),
#      col=c('black','red'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_pol.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1Pol'),
#      col=c('black','blue'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_hs.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1HS'),
#      col=c('black','green'),xlab='kb',ask=FALSE)
# dev.off()


#################### IWT test ####################
#### all regions, zero fixed, with strand ####
##############################################
load('L1_autosomes_zero_fixed_with_strand_smoothed_new_revision.RData')

result_mean=IWTomicsTest(regionsFeatures_smoothed,
                         id_region1=c("L1denovoFlasch","L1denovoSultana","L1denovoFlasch","L1denovoSultana"),
                         id_region2=c("Control","Control","L1HS","L1HS"),
                         statistics='mean',B=10000)
save(result_mean,file='L1_autosomes_results_smoothed_mean_new_revision.RData')


##Mar22, 2020
load('L1_autosomes_results_smoothed_mean_new_revision.RData')

pdf('IWT_autosomes_smoothed_mean_new_revision.pdf',width=7,height=10)
plotTest(result_mean,col=c('orange','purple','black','green'),
         scale_threshold=unlist(lapply(result_mean@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
dev.off()

##Switch id_regions1 and 2 to present comparisons between L1HS vs L1denovoFlasch and L1HS vs L1denovoSultana
result_mean@test$input$id_region1 = c("L1denovoFlasch","L1denovoSultana","L1HS","L1HS")
result_mean@test$input$id_region2=c("Control","Control","L1denovoFlasch","L1denovoSultana")
##change signs in all the corresponding mean differences from the test results (test 3 and 4)
result_mean@test$result[[3]]$H2AFZ_signal$T0_plot <- -result_mean@test$result[[3]]$H2AFZ_signal$T0_plot
result_mean@test$result[[3]]$H3K27ac_signal$T0_plot <- -result_mean@test$result[[3]]$H3K27ac_signal$T0_plot
result_mean@test$result[[3]]$H4K20me1_signal$T0_plot <- -result_mean@test$result[[3]]$H4K20me1_signal$T0_plot
result_mean@test$result[[3]]$H3K36me3_signal$T0_plot <- -result_mean@test$result[[3]]$H3K36me3_signal$T0_plot
result_mean@test$result[[3]]$H3K4me1_signal$T0_plot <- -result_mean@test$result[[3]]$H3K4me1_signal$T0_plot
result_mean@test$result[[3]]$H3K4me2_signal$T0_plot <- -result_mean@test$result[[3]]$H3K4me2_signal$T0_plot
result_mean@test$result[[3]]$H3K4me3_signal$T0_plot <- -result_mean@test$result[[3]]$H3K4me3_signal$T0_plot
result_mean@test$result[[3]]$H3K79me2_signal$T0_plot <- -result_mean@test$result[[3]]$H3K79me2_signal$T0_plot
result_mean@test$result[[3]]$H3K9ac_signal$T0_plot <- -result_mean@test$result[[3]]$H3K9ac_signal$T0_plot
result_mean@test$result[[3]]$H3K9me3_signal$T0_plot <- -result_mean@test$result[[3]]$H3K9me3_signal$T0_plot
result_mean@test$result[[3]]$H3K27me3_signal$T0_plot <- -result_mean@test$result[[3]]$H3K27me3_signal$T0_plot
result_mean@test$result[[3]]$CTCF_signal$T0_plot <- -result_mean@test$result[[3]]$CTCF_signal$T0_plot
result_mean@test$result[[3]]$DNase_DHS_signal$T0_plot <- -result_mean@test$result[[3]]$DNase_DHS_signal$T0_plot

result_mean@test$result[[4]]$H2AFZ_signal$T0_plot <- -result_mean@test$result[[4]]$H2AFZ_signal$T0_plot
result_mean@test$result[[4]]$H3K27ac_signal$T0_plot <- -result_mean@test$result[[4]]$H3K27ac_signal$T0_plot
result_mean@test$result[[4]]$H4K20me1_signal$T0_plot <- -result_mean@test$result[[4]]$H4K20me1_signal$T0_plot
result_mean@test$result[[4]]$H3K36me3_signal$T0_plot <- -result_mean@test$result[[4]]$H3K36me3_signal$T0_plot
result_mean@test$result[[4]]$H3K4me1_signal$T0_plot <- -result_mean@test$result[[4]]$H3K4me1_signal$T0_plot
result_mean@test$result[[4]]$H3K4me2_signal$T0_plot <- -result_mean@test$result[[4]]$H3K4me2_signal$T0_plot
result_mean@test$result[[4]]$H3K4me3_signal$T0_plot <- -result_mean@test$result[[4]]$H3K4me3_signal$T0_plot
result_mean@test$result[[4]]$H3K79me2_signal$T0_plot <- -result_mean@test$result[[4]]$H3K79me2_signal$T0_plot
result_mean@test$result[[4]]$H3K9ac_signal$T0_plot <- -result_mean@test$result[[4]]$H3K9ac_signal$T0_plot
result_mean@test$result[[4]]$H3K9me3_signal$T0_plot <- -result_mean@test$result[[4]]$H3K9me3_signal$T0_plot
result_mean@test$result[[4]]$H3K27me3_signal$T0_plot <- -result_mean@test$result[[4]]$H3K27me3_signal$T0_plot
result_mean@test$result[[4]]$CTCF_signal$T0_plot <- -result_mean@test$result[[4]]$CTCF_signal$T0_plot
result_mean@test$result[[4]]$DNase_DHS_signal$T0_plot <- -result_mean@test$result[[4]]$DNase_DHS_signal$T0_plot

save(result_mean,file='L1_autosomes_results_smoothed_mean_new_revision_direction_changed.RData')


#Change order of features to the same in our current study
features_reordered_revision<-read.table("features_reordered_revision.txt",header=FALSE)$V1
features_reordered_revision_ready<-as.vector(features_reordered_revision)

#Summary by test, scale=100
plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("L1denovoFlasch_control","L1denovoSultana_control","hs_L1denovoFlasch","hs_L1denovoSultana"),".pdf"),
            id_features_subset=features_reordered_revision_ready,
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
#Summary by test, scale=10
plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',scale_threshold=10,
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("L1denovoFlasch_control","L1denovoSultana_control","hs_L1denovoFlasch","hs_L1denovoSultana"),"_scale10.pdf"),
            id_features_subset=features_reordered_revision_ready,
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)

#Summary by features, scale=100
plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),".pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
#Summary by features, scale=10
plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',scale_threshold=10,
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),"_scale10.pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)





