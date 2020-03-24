require(IWTomics)
require(dendextend)

setwd("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/")
#setwd("/Users/Bruce/Desktop/L1/IWTomics analysis_1000control_1000denovo")


#####Subsample de novo L1s by taking 1000 random regions including X
denovo_all<-read.table("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/new_L1_data/Flasch_denovo_hESC_100k_no_gaps_no_blacklisted_method2.bed",header=FALSE)
##Create first 1000 random denovo 
# denovo1000<-denovo_all[sample(nrow(denovo_all), 1000), ]
# write.table(denovo1000, file="~/Desktop/L1/files/denovo1000_A_unsorted.bed", sep = "\t",row.names = FALSE,
#             col.names = FALSE, quote = FALSE)
# #sort by bedtools
# ##sortBed -i denovo1000_A_unsorted.bed > denovo1000_A.bed
# 
# #####Subsample the original control by taking 1000 random regions including X
# control_all<-read.table("~/Desktop/L1/files/100k_Control_noHS_noPol_nodenovo_nodbRIP.bed",header=FALSE)
# ##Create first 1000 random control 
# control1000<-control_all[sample(nrow(control_all), 1000), ]
# write.table(control1000, file="~/Desktop/L1/files/control1000_A_unsorted.bed", sep = "\t",row.names = FALSE,
#             col.names = FALSE, quote = FALSE)
#sort by bedtools
##sortBed -i control1000_A_unsorted.bed > control1000_A.bed






# files with datasets and features
datasets=read.table("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/datasets_new_revision.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
features_datasets=read.table("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/IWT_Revision/features_datasets_new_revision.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
# exclude low-resoution features (analyzed in IWT_L1_low_resolution.r)
features_datasets=features_datasets[1:13,]
tail(features_datasets)

# load data
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

# # select only chrX
# load('L1_complete.RData')
# index=lapply(regionsFeatures@regions,function(region) seqnames(region)=='chrX')
# regionsFeatures@metadata$region_datasets$size=unlist(lapply(index,sum))
# regionsFeatures@regions=GRangesList(mapply(function(region,ind) region[ind,],regionsFeatures@regions,index,SIMPLIFY=FALSE))
# regionsFeatures@features=lapply(regionsFeatures@features,function(feature) mapply(function(feat,ind) feat[,which(ind)],feature,index,SIMPLIFY=FALSE))
# regionsFeatures@length_features=lapply(regionsFeatures@length_features,function(feature) mapply(function(feat,ind) feat[which(ind)],feature,index,SIMPLIFY=FALSE))
# validObject(regionsFeatures)
# # number of windows in each dataset
# lengthRegions(regionsFeatures)
# save(regionsFeatures,file='L1_chrX.RData')













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

# #many_zeros=c("CpG_Islands","5hMc","Sperm_hypometh")
# #######################
# ##### NOTE: before when we weren't selecting only controls that overlap <7% with non human-specific L1, 
# #####       we had more features with less than 10% non zeros:
# #####       many_zeros=c("Quadruplex","Exons","CpG_Islands","5hMc","Sperm_hypometh","Exon_Expression")
# # smooth measurements... and keep 100 windows
# band=rep(1,length(many_zeros))
# names(band)=many_zeros
# for(id in many_zeros){
#   fix=TRUE
#   while(fix){
#     band[id]=band[id]+1
#     regionsFeatures_new=smooth(regionsFeatures,id_features_subset=id,type='kernel',bandwidth=band[id])
#     zero_count=sum(unlist(regionsFeatures_new@features[[id]])>0)/length(unlist(regionsFeatures_new@features[[id]]))*100
#     fix=!(zero_count>10)
#   }
#   regionsFeatures=regionsFeatures_new
# }
# band
# #band=c(3,2,2)
# #######################
# ##### NOTE: before when we weren't selecting only controls that overlap <7% with non human-specific L1, 
# #####       we had to smooth more CpG_Islands in order to have more than 10% non zeros:
# #####       band=c(2,2,5,2,2,2)
# zero_count_tot=Reduce(rbind,lapply(regionsFeatures@features,
#                                    function(feature){
#                                      feat=unlist(feature)
#                                      count=c(sum(feat==0),sum(feat!=0),length(feat),length(unique(feat)))
#                                      names(count)=c('0','>0','tot','distinct')
#                                      return(count)
#                                    }))
# rownames(zero_count_tot)=idFeatures(regionsFeatures)
# zero_count_tot
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


# pdf('boxplot_autosomes_denovo.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1denovo'),
#      col=c('black','red'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_pol.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1Pol'),
#      col=c('black','blue'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_hs.pdf',width=10,height=7)
# plot(regionsFeatures,type='boxplot',id_regions_subset=c('Control','L1HS'),
#      col=c('black','green'),xlab='kb',ask=FALSE)
# dev.off()



###########

#index=sort(sample(regionsFeatures_smoothed@metadata$region_datasets['L1denovo','size'],1000))


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










# smooth everything but the features that were already smoothed because they had too many zeros
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






# 
# 
# 
# # select the "right" threshold
# load('L1_autosomes_results_smoothed_mean.RData')
# test=c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs")
# names(test)=test
# select_scale=lapply(test,
#                     function(test){
#                       x=rep(100,nFeatures(result_mean))
#                       names(x)=idFeatures(result_mean)
#                       return(x)})
# select_scale$denovo_control[c(19,28,38)]=c(4,4,10)
# select_scale$pol_control[c(1,4,6,7,8,9,11,12,13,21,23,24,25,29,30,32,34,38,41)]=c(8,20,10,8,10,5,4,2,8,20,4,5,5,5,4,4,4,4,2)
# select_scale$hs_control[c(11,14,15,16,19,31,32,33,34,38,43,44)]=c(20,10,4,4,8,2,2,4,10,4,10,2)
# pdf('IWT_autosomes_smoothed_mean_select_scale.pdf',width=7,height=10)
# plotTest(result_mean,col=c('red','blue','green','black'),
#          scale_threshold=select_scale,ask=FALSE)
# dev.off()
# plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',scale_threshold=select_scale,
#             filenames=paste0("IWT_autosomes_smoothed_mean_select_scale_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',scale_threshold=10,
#             filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),"_scale10.pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',scale_threshold=10,
#             filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),"_scale10.pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# result_median=IWTomicsTest(regionsFeatures_smoothed,
#                            id_region1=c("L1denovo","L1Pol","L1HS","L1denovo","L1denovo","L1Pol"),
#                            id_region2=c("Control","Control","Control","L1Pol","L1HS","L1HS"),
#                            statistics='median',B=10000)
# save(result_median,file='L1_autosomes_results_smoothed_median.RData')
# load('L1_autosomes_results_smoothed_median.RData')
# pdf('IWT_autosomes_smoothed_median.pdf',width=7,height=10)
# plotTest(result_median,col=c('red','blue','green','black'),
#          scale_threshold=unlist(lapply(result_median@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
# dev.off()
# plotSummary(result_median,groupby="test",only_significant=FALSE,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_median_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_median,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_median_summary_feature_",idFeatures(result_median),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# result_variance=IWTomicsTest(regionsFeatures_smoothed,
#                              id_region1=c("L1denovo","L1Pol","L1HS","L1denovo","L1denovo","L1Pol"),
#                              id_region2=c("Control","Control","Control","L1Pol","L1HS","L1HS"),
#                              statistics='variance',B=10000)
# save(result_variance,file='L1_autosomes_results_smoothed_variance.RData')
# load('L1_autosomes_results_smoothed_variance.RData')
# pdf('IWT_autosomes_smoothed_variance.pdf',width=7,height=10)
# plotTest(result_variance,col=c('red','blue','green','black'),
#          scale_threshold=unlist(lapply(result_variance@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
# dev.off()
# plotSummary(result_variance,groupby="test",only_significant=FALSE,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_variance_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_variance,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_variance_summary_feature_",idFeatures(result_variance),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# result_quantile=IWTomicsTest(regionsFeatures_smoothed,
#                              id_region1=c("L1denovo","L1Pol","L1HS","L1denovo","L1denovo","L1Pol"),
#                              id_region2=c("Control","Control","Control","L1Pol","L1HS","L1HS"),
#                              statistics='quantile',probs=0.9,B=10000)
# save(result_quantile,file='L1_autosomes_results_smoothed_quantile_90.RData')
# load('L1_autosomes_results_smoothed_quantile_90.RData')
# pdf('IWT_autosomes_smoothed_quantile_90.pdf',width=7,height=10)
# plotTest(result_quantile,col=c('red','blue','green','black'),
#          scale_threshold=unlist(lapply(result_quantile@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
# dev.off()
# plotSummary(result_quantile,groupby="test",only_significant=FALSE,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_quantile_90_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_quantile,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_quantile_90_summary_feature_",idFeatures(result_quantile),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# result_multi_quantile=IWTomicsTest(regionsFeatures_smoothed,
#                                    id_region1=c("L1denovo","L1Pol","L1HS","L1denovo","L1denovo","L1Pol"),
#                                    id_region2=c("Control","Control","Control","L1Pol","L1HS","L1HS"),
#                                    statistics='quantile',probs=c(0.05,0.25,0.5,0.75,0.95),B=10000)
# save(result_multi_quantile,file='L1_autosomes_results_smoothed_quantile_5_25_50_75_95.RData')
# load('L1_autosomes_results_smoothed_quantile_5_25_50_75_95.RData')
# pdf('IWT_autosomes_smoothed_quantile_5_25_50_75_95.pdf',width=7,height=10)
# plotTest(result_multi_quantile,col=c('red','blue','green','black'),
#          scale_threshold=unlist(lapply(result_multi_quantile@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
# dev.off()
# plotSummary(result_multi_quantile,groupby="test",only_significant=FALSE,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_quantile_5_25_50_75_95_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_multi_quantile,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_quantile_5_25_50_75_95_summary_feature_",idFeatures(result_multi_quantile),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# 
# 
# 
# #################### test ####################
# #### denovo different number of overlaps  ####
# ##############################################
# load('L1_autosomes_zero_fixed_with_strand_smoothed_denovo_overlaps.RData')
# result_mean=IWTomicsTest(regionsFeatures_smoothed,
#                          id_region1=c("L1denovo0","L1denovo1","L1denovo2","L1denovo>2"),
#                          id_region2=c("Control","Control","Control","Control"),
#                          statistics='mean',B=10000)
# save(result_mean,file='L1_autosomes_results_smoothed_denovo_overlaps_mean.RData')
# load('L1_autosomes_results_smoothed_denovo_overlaps_mean.RData')
# pdf('IWT_autosomes_smoothed_denovo_overlaps_mean.pdf',width=7,height=10)
# plotTest(result_mean,col=c('red','blue','green','black'),
#          scale_threshold=unlist(lapply(result_mean@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
# dev.off()
# result_mean_denovo_overlaps=result_mean
# load('L1_autosomes_results_smoothed_mean.RData')
# result_mean_others=result_mean
# result_mean=rbind(result_mean_denovo_overlaps,result_mean_others)
# result_mean@test=result_mean_denovo_overlaps@test
# result_mean@test$input$id_region1=c(result_mean@test$input$id_region1,result_mean_others@test$input$id_region1)
# result_mean@test$input$id_region2=c(result_mean@test$input$id_region2,result_mean_others@test$input$id_region2)
# result_mean@test$result=c(result_mean@test$result,result_mean_others@test$result)
# validObject(result_mean)
# plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',test=1:4,
#             filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo0_control","denovo1_control","denovo2_control","denovomore2_pol"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=c(4,7),xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_denovo_overlaps_mean_summary_feature_",idFeatures(result_mean),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #################### test ####################
# ####   regions with only one type of L1   ####
# ####        zero fixed, with strand       ####
# ##############################################
# load('L1_autosomes_zero_fixed_with_strand_smoothed.RData')
# no_overlaps=lapply(regionsFeatures_smoothed@regions[1:3],function(regions) regions$overlap_other=='no')
# regionsFeatures_smoothed@regions[1:3]=GRangesList(mapply(function(regions,no) regions[no],regionsFeatures_smoothed@regions[1:3],no_overlaps,SIMPLIFY=FALSE))
# regionsFeatures_smoothed@features=lapply(regionsFeatures_smoothed@features,
#                                          function(feature){
#                                            feature[1:3]=mapply(function(feat,no) feat[,no],feature[1:3],no_overlaps,SIMPLIFY=FALSE)
#                                            return(feature)
#                                          })
# regionsFeatures_smoothed@length_features=lapply(regionsFeatures_smoothed@length_features,
#                                                 function(feature){
#                                                   feature[1:3]=mapply(function(feat,no) feat[no],feature[1:3],no_overlaps,SIMPLIFY=FALSE)
#                                                   return(feature)
#                                                 })
# validObject(regionsFeatures_smoothed)
# save(regionsFeatures_smoothed,file='L1_autosomes_zero_fixed_with_strand_smoothed_only_one_L1.RData')
# 
# # plot
# pdf('curves_autosomes_with_strand_smoothed_only_one_L1.pdf',width=10,height=8)
# plot(regionsFeatures_smoothed,type='curves',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#      col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()
# 
# pdf('means_autosomes_with_strand_smoothed_only_one_L1.pdf',width=10,height=8)
# plot_only_means(regionsFeatures_smoothed,id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#                 col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('means_autosomes_with_strand_smoothed_only_one_L1_no_zero.pdf',width=10,height=8)
# plot_only_means(regionsFeatures_smoothed,id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),include_zero=FALSE,
#                 col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()
# 
# pdf('boxplot_autosomes_with_strand_smoothed_only_one_L1.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1denovo','L1Pol','L1HS'),
#      col=c('black','red','blue','green'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_only_one_L1_denovo.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1denovo'),
#      col=c('black','red'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_only_one_L1_pol.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1Pol'),
#      col=c('black','blue'),xlab='kb',ask=FALSE)
# dev.off()
# pdf('boxplot_autosomes_with_strand_smoothed_only_one_L1_hs.pdf',width=10,height=7)
# plot(regionsFeatures_smoothed,type='boxplot',id_regions_subset=c('Control','L1HS'),
#      col=c('black','green'),xlab='kb',ask=FALSE)
# dev.off()
# 
# 
# # test
# load('L1_autosomes_zero_fixed_with_strand_smoothed_only_one_L1.RData')
# result_mean=IWTomicsTest(regionsFeatures_smoothed,
#                          id_region1=c("L1denovo","L1Pol","L1HS","L1denovo","L1denovo","L1Pol"),
#                          id_region2=c("Control","Control","Control","L1Pol","L1HS","L1HS"),
#                          statistics='mean',B=10000)
# save(result_mean,file='L1_autosomes_results_smoothed_only_one_L1_mean.RData')
# load('L1_autosomes_results_smoothed_only_one_L1_mean.RData')
# pdf('IWT_autosomes_smoothed_only_one_L1_mean.pdf',width=7,height=10)
# plotTest(result_mean,col=c('red','blue','green','black'),
#          scale_threshold=unlist(lapply(result_mean@length_features,function(feat) unique(unlist(feat)))),ask=FALSE)
# dev.off()
# plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_only_one_L1_mean_summary_test_",c("denovo_control","pol_control","hs_control","denovo_pol","denovo_hs","pol_hs"),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',
#             filenames=paste0("IWT_autosomes_smoothed_only_one_L1_mean_summary_feature_",idFeatures(result_mean),".pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)
# 
# 
# 
# 





