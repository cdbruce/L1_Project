##Karyotype plots for L1s from three datasets
##Backup scripts for jupytor notebook (L1s_Karyotype_SizeCorrelation.ipynb)

setwd("~/Google_Drive/L1/L1_Project/Analysis/Karyotype_plots/")

library("rtracklayer")
library(karyoploteR)
library(GenomicRanges)

#Load the data for de novo L1s  
#L1_denovo<-read.table("L1denovo_BWA_merged4_undecidable_100kb_no_gaps_no_blacklisted.interval",sep="\t")[,1:3]
L1_denovo<-read.table("L1denovo_BWA_17037_for_karyotype.bed",sep="\t")[,1:3]
names(L1_denovo)=c('chr','start','end')
#class(G137_Br)
#head(G137_Br)
## Karyotype plot with "karyoploteR", 
gains_denovo <- makeGRangesFromDataFrame(L1_denovo) ## Import target positions
#head(gains)
length(gains_denovo)
## Plot de novo
pdf('L1_denovo.pdf',width=15, height=13)
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.1
pp$data1height <- 80
kp<-plotKaryotype(genome="hg19", main="de novo L1s",plot.type=1, plot.params = pp) ## Set genome assembly
#kpPlotRegions(kp, gains,col="red",avoid.overlapping=FALSE ) ## Choose color
#getCytobandColors(color.table=NULL, color.schema=c("only.centromeres"))
kpPlotRegions(kp, gains_denovo,col="red") ## Choose color 
kpAddCytobandLabels(kp,cex=0.5)
dev.off()

#Load the data for polymorpihc L1s 
#L1_pol<-read.table("L1_Pol_clean.bed",sep="\t")[,1:3]
L1_pol<-read.table("PolyL1_Ewing_LiftOver_18to19_1012_for_karyotype.bed",sep="\t")[,1:3]
names(L1_pol)=c('chr','start','end')
#class(G137_Br)
#head(G137_Br)
## Karyotype plot with "karyoploteR", 
gains_pol <- makeGRangesFromDataFrame(L1_pol) ## Import target positions
#head(gains)
length(gains_pol)
## Plot L1 pol
pdf('L1_pol.pdf',width=15, height=13)
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.1
pp$data1height <- 80
kp<-plotKaryotype(genome="hg19", main="Polymorphic L1s",plot.type=1, plot.params = pp) ## Set genome assembly
#kpPlotRegions(kp, gains,col="red",avoid.overlapping=FALSE ) ## Choose color
getCytobandColors(color.table=NULL, color.schema=c("only.centromeres"))
kpPlotRegions(kp, gains_pol,col="blue") ## Choose color 
kpAddCytobandLabels(kp,cex=0.5)
dev.off()

#Load the data for L1HS
#L1_hs<-read.table("L1HS_UCSC.bed",sep="\t")[,1:3]
L1_hs<-read.table("L1HS_clean_sorted_1205.bed",sep="\t")[,1:3]
names(L1_hs)=c('chr','start','end')
#class(G137_Br)
#head(G137_Br)
## Karyotype plot with "karyoploteR", 
gains_hs <- makeGRangesFromDataFrame(L1_hs) ## Import target positions
#head(gains)
length(gains_hs)
## Plot L1 hs
pdf('L1_hs.pdf',width=15, height=13)
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.1
pp$data1height <- 80
kp<-plotKaryotype(genome="hg19", main="Human-specific L1s",plot.type=1, plot.params = pp) ## Set genome assembly
#kpPlotRegions(kp, gains,col="red",avoid.overlapping=FALSE ) ## Choose color
getCytobandColors(color.table=NULL, color.schema=c("only.centromeres"))
kpPlotRegions(kp, gains_hs,col="green") ## Choose color 
kpAddCytobandLabels(kp,cex=0.5)
dev.off()

#Overlay all L1s in the same karyotype plot
pdf('L1_ALL.pdf',width=15, height=13)
pp <- getDefaultPlotParams(plot.type=1)
pp$leftmargin <- 0.1
pp$data1height <- 60
kp<-plotKaryotype(genome="hg19", main="Genome-wide distribution of all L1s in the study",plot.type=1, plot.params = pp, cex=1.6) ## Set genome assembly
#kpPlotRegions(kp, gains,col="red",avoid.overlapping=FALSE ) ## Choose color
getCytobandColors(color.table=NULL, color.schema=c("only.centromeres"))
kpPlotRegions(kp, r0 = 0, r1 = 2,gains_denovo,col="red") ## Choose color
kpPlotRegions(kp, r0 = 0, r1 = 0.8,gains_pol, col="blue") ## Choose color 
kpPlotRegions(kp, r0 = 0, r1 = 0.8,gains_hs, col="green") ## Choose color 
kpAddCytobandLabels(kp,cex=0.8)
dev.off()