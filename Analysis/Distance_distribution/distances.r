setwd("~/Google_Drive/L1/L1_Project/Analysis/Distance_distribution/")

require(GenomicRanges)



# DE NOVO
L1denovo  <- read.table('~/Google_Drive/L1/L1_Project/Datasets/L1_genomic_coordinates/L1denovo_BWA_17037_reads.bed', header = FALSE, sep = '\t')
names(L1denovo) <- c('chr', 'start', 'end', 'strand', 'barcode', 'site', 'annot')
L1denovo$start[L1denovo$strand == '+'] <- L1denovo$end[L1denovo$strand == '+']
L1denovo$end[L1denovo$strand == '-'] <- L1denovo$start[L1denovo$strand == '-']
L1denovo <- makeGRangesFromDataFrame(L1denovo, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
L1denovo <- sort(L1denovo, ignore.strand = TRUE)

# POLYMORPHIC
L1Pol <- read.table('~/Google_Drive/L1/L1_Project/Datasets/L1_genomic_coordinates/L1Pol_Ewing_LiftedFromHG18.interval', header = FALSE, sep = '\t')
names(L1Pol) <- c('chr', 'start', 'end', 'strand', 'annot')
L1Pol <- makeGRangesFromDataFrame(L1Pol, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
L1Pol <- sort(L1Pol, ignore.strand = TRUE)

# HUMAN SPECIFIC
L1HS <- read.table('~/Google_Drive/L1/L1_Project/Datasets/L1_genomic_coordinates/L1HS_clean_sorted_1205.bed', header = FALSE, sep = '\t')
names(L1HS) <- c('chr', 'start', 'end', 'strand', 'annot')
L1HS <- makeGRangesFromDataFrame(L1HS, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
L1HS <- sort(L1HS, ignore.strand = TRUE)


# Check that they don't overlap gaps
gaps <- read.table('gaps.bed', header = FALSE, sep = '\t')
names(gaps) <- c('chr', 'start', 'end')
gaps <- makeGRangesFromDataFrame(gaps, starts.in.df.are.0based = TRUE, keep.extra.columns = TRUE)
findOverlaps(L1denovo, gaps)
findOverlaps(L1Pol, gaps)
findOverlaps(L1HS, gaps)

# remove HS L1s that overlap gaps
L1HS <- L1HS[-as.matrix(findOverlaps(L1HS, gaps))[,1], ]

# removing entry 250 because for some reason it's too hard to shuffle it!
# shuffleBed error:
# Error, line 250: tried 1000 potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on.
L1HS <- L1HS[-250, ]





#### Distance of L1 elements of the same type: distribution ####
dist_denovo <- (start(L1denovo)[-1]-1)-end(L1denovo)[-length(L1denovo)]
dist_denovo <- dist_denovo[dist_denovo>=0]
dist_Pol <- (start(L1Pol)[-1]-1)-end(L1Pol)[-length(L1Pol)]
dist_Pol <- dist_Pol[dist_Pol>=0]
dist_HS <- (start(L1HS)[-1]-1)-end(L1HS)[-length(L1HS)]
dist_HS <- dist_HS[dist_HS>=0]


# histogram
pdf("distance_histogram_log.pdf", width = 6, height = 6)
hist(log(dist_Pol+0.01),col=rgb(0,0,1,0.3),probability=TRUE,breaks=40,xlim=c(-5,20),ylim=c(0,0.35),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_HS+0.01),col=rgb(0,1,0,0.3),probability=TRUE,breaks=40,add=TRUE)
hist(log(dist_denovo+0.01),col=rgb(1,0,0,0.3),probability=TRUE,breaks=40,add=TRUE)
legend('topleft',legend=c('De novo L1','Polymorphic L1','Human specific L1'),fill=c(rgb(1,0,0,0.3),rgb(0,0,1,0.3),rgb(0,1,0,0.3)),bg='white')
dev.off()


# cumulative distribution
ecdf_log_denovo=ecdf(log(dist_denovo+0.01))
ecdf_log_Pol=ecdf(log(dist_Pol+0.01))
ecdf_log_HS=ecdf(log(dist_HS+0.01))
pdf("distance_cumulative_log.pdf", width = 6, height = 6)
plot(ecdf_log_denovo,col='red',do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_Pol,col='blue',do.points=FALSE,lwd=2,add=TRUE)
plot(ecdf_log_HS,col='green',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('De novo L1','Polymorphic L1','Human specific L1'),col=c('red','blue','green'),lty=1,lwd=2,bg='white')
dev.off()
###-------


#### Distance of L1 elements of the same type: distribution vs random ####
shuffleBed <- function(x, genome_file = NULL, excl_file = NULL, chrom = FALSE, maxTries = 1000){
  require(rtracklayer)
  
  stopifnot(class(x) == "GRanges")
  if (nchar(Sys.which("shuffleBed")) == 0) stop("shuffleBed must be in the PATH")
  if(is.null(genome_file)){
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    genome_file <- tempfile()
    writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file)
  }
  bed_file <- tempfile()
  writeLines(paste(as.character(seqnames(x)), start(x)-1, end(x), sep="\t"), bed_file)
  
  shuffle_file <- tempfile()
  to_run="shuffleBed"
  if(!is.null(excl_file))
    to_run <- paste(to_run, "-excl", excl_file)
  if(chrom)
    to_run <- paste(to_run, "-chrom")
  to_run <- paste(to_run, "-maxTries", maxTries, "-i", bed_file, "-g", genome_file, ">", shuffle_file)
  system(to_run)
  shuffled <- import.bed(shuffle_file)
  seqlevels(shuffled) <- seqlevels(x)
  seqlengths(shuffled) <- seqlengths(x)
  
  #clean up
  unlink(c(bed_file, shuffle_file))
  shuffled
}
decdf <- function(x, a, b)  ecdf(a)(x) - ecdf(b)(x)

# Denovo not same chromosome
set.seed(2017)
shuffled <- shuffleBed(L1denovo, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
shuffled <- sort(shuffled, ignore.strand= TRUE)
dist_denovo_shuffled <- (start(shuffled)[-1]-1)-end(shuffled)[-length(shuffled)]
dist_denovo_shuffled <- dist_denovo_shuffled[dist_denovo_shuffled>=0]
require(Matching)
p.value=ks.boot(log(dist_denovo+0.01),log(dist_denovo_shuffled+0.01))$ks$p.value
if(p.value<1e-16)
  p.value='<1e-16'

# histogram
pdf("distance_denovo_histogram_log.pdf", width = 6, height = 6)
hist(log(dist_denovo+0.01),col=rgb(1,0,0,0.3),probability=TRUE,breaks=40,xlim=c(-5,20),ylim=c(0,0.35),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_denovo_shuffled+0.01),col=rgb(0,0,0,0.3),probability=TRUE,breaks=20,add=TRUE)
legend('topleft',legend=c('De novo L1','Random',paste0('p-value ',format(p.value,digits=2))),fill=c(rgb(1,0,0,0.3),rgb(0,0,0,0.3),'white'),border=c('black','black','white'),bg='white')
dev.off()


# cumulative distribution
ecdf_log_denovo=ecdf(log(dist_denovo+0.01))
ecdf_log_denovo_shuffled=ecdf(log(dist_denovo_shuffled+0.01))
pdf("distance_denovo_cumulative_log.pdf", width = 6, height = 6)
plot(ecdf_log_denovo,col='red',do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_denovo_shuffled,col='black',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('De novo L1','Random',paste0('p-value ',format(p.value,digits=2))),col=c('red','black','white'),lty=1,lwd=2,bg='white')
dev.off()


# q-q plot
pdf("distance_denovo_qqplot_log.pdf", width = 6, height = 6)
qqplot(log(dist_denovo_shuffled+0.01),log(dist_denovo+0.01),main='Log Distance between L1',xlab='Random',ylab='De novo L1')
abline(a=0,b=1,col='red')
legend('bottomright',legend=paste0('p-value ',format(p.value,digits=2)),col='white',lty=1,lwd=2,bg='white')
dev.off()




# Pol not same chromosome
set.seed(2018)
shuffled <- shuffleBed(L1Pol, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
shuffled <- sort(shuffled, ignore.strand= TRUE)
dist_Pol_shuffled <- (start(shuffled)[-1]-1)-end(shuffled)[-length(shuffled)]
dist_Pol_shuffled <- dist_Pol_shuffled[dist_Pol_shuffled>=0]
require(Matching)
p.value=ks.boot(log(dist_Pol+0.01),log(dist_Pol_shuffled+0.01))$ks$p.value
if(p.value<1e-16)
  p.value='<1e-16'

# histogram
pdf("distance_Pol_histogram_log.pdf", width = 6, height = 6)
hist(log(dist_Pol+0.01),col=rgb(0,0,1,0.3),probability=TRUE,breaks=40,xlim=c(-5,20),ylim=c(0,0.35),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_Pol_shuffled+0.01),col=rgb(0,0,0,0.3),probability=TRUE,breaks=20,add=TRUE)
legend('topleft',legend=c('Polymorphic L1','Random',paste0('p-value ',format(p.value,digits=2))),fill=c(rgb(0,0,1,0.3),rgb(0,0,0,0.3),'white'),border=c('black','black','white'),bg='white')
dev.off()


# cumulative distribution
ecdf_log_Pol=ecdf(log(dist_Pol+0.01))
ecdf_log_Pol_shuffled=ecdf(log(dist_Pol_shuffled+0.01))
pdf("distance_Pol_cumulative_log.pdf", width = 6, height = 6)
plot(ecdf_log_Pol,col='blue',do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_Pol_shuffled,col='black',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('Polymorphic L1','Random',paste0('p-value ',format(p.value,digits=2))),col=c('blue','black','white'),lty=1,lwd=2,bg='white')
dev.off()


# q-q plot
pdf("distance_Pol_qqplot_log.pdf", width = 6, height = 6)
qqplot(log(dist_Pol_shuffled+0.01),log(dist_Pol+0.01),main='Log Distance between L1',xlab='Random',ylab='Polymorphic L1')
abline(a=0,b=1,col='blue')
legend('bottomright',legend=paste0('p-value ',format(p.value,digits=2)),col='white',lty=1,lwd=2,bg='white')
dev.off()




# HS not same chromosome
set.seed(2018)
shuffled <- shuffleBed(L1HS, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
shuffled <- sort(shuffled, ignore.strand= TRUE)
dist_HS_shuffled <- (start(shuffled)[-1]-1)-end(shuffled)[-length(shuffled)]
dist_HS_shuffled <- dist_HS_shuffled[dist_HS_shuffled>=0]
require(Matching)
p.value=ks.boot(log(dist_HS+0.01),log(dist_HS_shuffled+0.01))$ks$p.value
if(p.value<1e-16)
  p.value='<1e-16'

# histogram
pdf("distance_HS_histogram_log.pdf", width = 6, height = 6)
hist(log(dist_HS+0.01),col=rgb(0,1,0,0.3),probability=TRUE,breaks=40,xlim=c(-5,20),ylim=c(0,0.35),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_HS_shuffled+0.01),col=rgb(0,0,0,0.3),probability=TRUE,breaks=20,add=TRUE)
legend('topleft',legend=c('Human specific L1','Random',paste0('p-value ',format(p.value,digits=2))),fill=c(rgb(0,1,0,0.3),rgb(0,0,0,0.3),'white'),border=c('black','black','white'),bg='white')
dev.off()


# cumulative distribution
ecdf_log_HS=ecdf(log(dist_HS+0.01))
ecdf_log_HS_shuffled=ecdf(log(dist_HS_shuffled+0.01))
pdf("distance_HS_cumulative_log.pdf", width = 6, height = 6)
plot(ecdf_log_HS,col='green',do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_HS_shuffled,col='black',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('Human specific L1','Random',paste0('p-value ',format(p.value,digits=2))),col=c('green','black','white'),lty=1,lwd=2,bg='white')
dev.off()


# q-q plot
pdf("distance_HS_qqplot_log.pdf", width = 6, height = 6)
qqplot(log(dist_HS_shuffled+0.01),log(dist_HS+0.01),main='Log Distance between L1',xlab='Random',ylab='Human specific L1')
abline(a=0,b=1,col='green')
legend('bottomright',legend=paste0('p-value ',format(p.value,digits=2)),col='white',lty=1,lwd=2,bg='white')
dev.off()


####--------


#### Distance of L1 elements of the same type: distribution NORMALIZED ####
decdf_log_denovo <- function(x) decdf(x,log(dist_denovo+0.01),log(dist_denovo_shuffled+0.01))
decdf_log_Pol <- function(x) decdf(x,log(dist_Pol+0.01),log(dist_Pol_shuffled+0.01))
decdf_log_HS <- function(x) decdf(x,log(dist_HS+0.01),log(dist_HS_shuffled+0.01))
pdf("distance_cumulative_log_diff.pdf", width = 6, height = 6)
curve(decdf_log_denovo,col='red',lwd=2,xlim=c(-5,20),ylim=c(-0.05,0.85),xlab='Log Distance',ylab='Fn(x) diff',main='Distance between L1')
curve(decdf_log_Pol,col='blue',lwd=2,add=TRUE)
curve(decdf_log_HS,col='green',lwd=2,add=TRUE)
legend('topleft',legend=c('De novo L1','Polymorphic L1','Human specific L1'),col=c('red','blue','green'),lty=1,lwd=2,bg='white')
dev.off()


####--------


#### SAME SAMPLE SIZE - MANY TIMES

# Maximum size
min( c(length(L1denovo), length(L1Pol), length(L1HS)) )

# SAMPLE SIZE 900
N <- 900

#### Distance of L1 elements of the same type: distribution ####
set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N)
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo <- (start(L1denovo_sample)[-1]-1)-end(L1denovo_sample)[-length(L1denovo_sample)]
dist_denovo <- dist_denovo[dist_denovo>=0]
dist_Pol <- (start(L1Pol_sample)[-1]-1)-end(L1Pol_sample)[-length(L1Pol_sample)]
dist_Pol <- dist_Pol[dist_Pol>=0]
dist_HS <- (start(L1HS_sample)[-1]-1)-end(L1HS_sample)[-length(L1HS_sample)]
dist_HS <- dist_HS[dist_HS>=0]

ecdf_log_denovo=ecdf(log(dist_denovo+0.01))
ecdf_log_Pol=ecdf(log(dist_Pol+0.01))
ecdf_log_HS=ecdf(log(dist_HS+0.01))
pdf(paste0(N, "_distance_cumulative_log_many_times.pdf"), width = 6, height = 6)
plot(ecdf_log_denovo,col='red',do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_Pol,col='blue',do.points=FALSE,add=TRUE)
plot(ecdf_log_HS,col='green',do.points=FALSE,add=TRUE)

for(i in 1:100){
  index_denovo <- sample(seq_along(L1denovo), N)
  index_Pol <- sample(seq_along(L1Pol), N)
  index_HS <- sample(seq_along(L1HS), N)
  L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
  L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
  L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)
  
  dist_denovo <- (start(L1denovo_sample)[-1]-1)-end(L1denovo_sample)[-length(L1denovo_sample)]
  dist_denovo <- dist_denovo[dist_denovo>=0]
  dist_Pol <- (start(L1Pol_sample)[-1]-1)-end(L1Pol_sample)[-length(L1Pol_sample)]
  dist_Pol <- dist_Pol[dist_Pol>=0]
  dist_HS <- (start(L1HS_sample)[-1]-1)-end(L1HS_sample)[-length(L1HS_sample)]
  dist_HS <- dist_HS[dist_HS>=0]
  
  ecdf_log_denovo=ecdf(log(dist_denovo+0.01))
  ecdf_log_Pol=ecdf(log(dist_Pol+0.01))
  ecdf_log_HS=ecdf(log(dist_HS+0.01))
  
  plot(ecdf_log_denovo,col='red',do.points=FALSE,add=TRUE)
  plot(ecdf_log_Pol,col='blue',do.points=FALSE,add=TRUE)
  plot(ecdf_log_HS,col='green',do.points=FALSE,add=TRUE)
}
legend('topleft',legend=c('De novo L1','Polymorphic L1','Human specific L1'),col=c('red','blue','green'),lty=1,lwd=2,bg='white')
dev.off()
####--------


#### Distance of L1 elements of the same type: distribution NORMALIZED ####
shuffleBed <- function(x, genome_file = NULL, excl_file = NULL, chrom = FALSE, maxTries = 1000){
  require(rtracklayer)
  
  stopifnot(class(x) == "GRanges")
  if (nchar(Sys.which("shuffleBed")) == 0) stop("shuffleBed must be in the PATH")
  if(is.null(genome_file)){
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    genome_file <- tempfile()
    writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file)
  }
  bed_file <- tempfile()
  writeLines(paste(as.character(seqnames(x)), start(x)-1, end(x), sep="\t"), bed_file)
  
  shuffle_file <- tempfile()
  to_run="shuffleBed"
  if(!is.null(excl_file))
    to_run <- paste(to_run, "-excl", excl_file)
  if(chrom)
    to_run <- paste(to_run, "-chrom")
  to_run <- paste(to_run, "-maxTries", maxTries, "-i", bed_file, "-g", genome_file, ">", shuffle_file)
  system(to_run)
  shuffled <- import.bed(shuffle_file)
  seqlevels(shuffled) <- seqlevels(x)
  seqlengths(shuffled) <- seqlengths(x)
  
  #clean up
  unlink(c(bed_file, shuffle_file))
  shuffled
}
decdf <- function(x, a, b)  ecdf(a)(x) - ecdf(b)(x)

set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N) 
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo <- (start(L1denovo_sample)[-1]-1)-end(L1denovo_sample)[-length(L1denovo_sample)]
dist_denovo <- dist_denovo[dist_denovo>=0]
dist_Pol <- (start(L1Pol_sample)[-1]-1)-end(L1Pol_sample)[-length(L1Pol_sample)]
dist_Pol <- dist_Pol[dist_Pol>=0]
dist_HS <- (start(L1HS_sample)[-1]-1)-end(L1HS_sample)[-length(L1HS_sample)]
dist_HS <- dist_HS[dist_HS>=0]

shuffled <- shuffleBed(L1denovo_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1denovo_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1Pol_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1Pol_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1HS_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1HS_shuffled <- sort(shuffled, ignore.strand = TRUE)

dist_denovo_shuffled <- (start(L1denovo_shuffled)[-1]-1)-end(L1denovo_shuffled)[-length(L1denovo_shuffled)]
dist_denovo_shuffled <- dist_denovo_shuffled[dist_denovo_shuffled>=0]
dist_Pol_shuffled <- (start(L1Pol_shuffled)[-1]-1)-end(L1Pol_shuffled)[-length(L1Pol_shuffled)]
dist_Pol_shuffled <- dist_Pol_shuffled[dist_Pol_shuffled>=0]
dist_HS_shuffled <- (start(L1HS_shuffled)[-1]-1)-end(L1HS_shuffled)[-length(L1HS_shuffled)]
dist_HS_shuffled <- dist_HS_shuffled[dist_HS_shuffled>=0]

decdf_log_denovo=function(x) decdf(x,log(dist_denovo+0.01),log(dist_denovo_shuffled+0.01))
decdf_log_Pol=function(x) decdf(x,log(dist_Pol+0.01),log(dist_Pol_shuffled+0.01))
decdf_log_HS=function(x) decdf(x,log(dist_HS+0.01),log(dist_HS_shuffled+0.01))

pdf(paste0(N,"_distance_cumulative_log_diff_many_times.pdf"), width = 6, height = 6)
curve(decdf_log_denovo,col='red',lwd=2,xlim=c(-5,20),ylim=c(-0.05,0.15),xlab='Log Distance',ylab='Fn(x) diff',main='Distance between L1')
curve(decdf_log_Pol,col='blue',add=TRUE)
curve(decdf_log_HS,col='green',add=TRUE)

for(i in 1:100){
  index_denovo <- sample(seq_along(L1denovo), N)
  index_Pol <- sample(seq_along(L1Pol), N)
  index_HS <- sample(seq_along(L1HS), N) 
  L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
  L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
  L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)
  
  dist_denovo <- (start(L1denovo_sample)[-1]-1)-end(L1denovo_sample)[-length(L1denovo_sample)]
  dist_denovo <- dist_denovo[dist_denovo>=0]
  dist_Pol <- (start(L1Pol_sample)[-1]-1)-end(L1Pol_sample)[-length(L1Pol_sample)]
  dist_Pol <- dist_Pol[dist_Pol>=0]
  dist_HS <- (start(L1HS_sample)[-1]-1)-end(L1HS_sample)[-length(L1HS_sample)]
  dist_HS <- dist_HS[dist_HS>=0]
  
  shuffled <- shuffleBed(L1denovo_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
  L1denovo_shuffled <- sort(shuffled, ignore.strand = TRUE)
  shuffled <- shuffleBed(L1Pol_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
  L1Pol_shuffled <- sort(shuffled, ignore.strand = TRUE)
  shuffled <- shuffleBed(L1HS_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
  L1HS_shuffled <- sort(shuffled, ignore.strand = TRUE)
  
  dist_denovo_shuffled <- (start(L1denovo_shuffled)[-1]-1)-end(L1denovo_shuffled)[-length(L1denovo_shuffled)]
  dist_denovo_shuffled <- dist_denovo_shuffled[dist_denovo_shuffled>=0]
  dist_Pol_shuffled <- (start(L1Pol_shuffled)[-1]-1)-end(L1Pol_shuffled)[-length(L1Pol_shuffled)]
  dist_Pol_shuffled <- dist_Pol_shuffled[dist_Pol_shuffled>=0]
  dist_HS_shuffled <- (start(L1HS_shuffled)[-1]-1)-end(L1HS_shuffled)[-length(L1HS_shuffled)]
  dist_HS_shuffled <- dist_HS_shuffled[dist_HS_shuffled>=0]
  
  decdf_log_denovo=function(x) decdf(x,log(dist_denovo+0.01),log(dist_denovo_shuffled+0.01))
  decdf_log_Pol=function(x) decdf(x,log(dist_Pol+0.01),log(dist_Pol_shuffled+0.01))
  decdf_log_HS=function(x) decdf(x,log(dist_HS+0.01),log(dist_HS_shuffled+0.01))
  
  curve(decdf_log_denovo,col='red',add=TRUE)
  curve(decdf_log_Pol,col='blue',add=TRUE)
  curve(decdf_log_HS,col='green',add=TRUE)
}
legend('topleft',legend=c('De novo L1','Polymorphic L1','Human specific L1'),col=c('red','blue','green'),lty=1,lwd=2,bg='white')
dev.off()


####--------








#### SAME SAMPLE SIZE

# Maximum size
min( c(length(L1denovo), length(L1Pol), length(L1HS)) )

# SAMPLE SIZE 900
N <- 900

#### Distance of L1 elements of different types: distribution ####
set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N)
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                     as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                    as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                 as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)


# histogram
pdf(paste0(N, "_distance_diff_type_histogram_log.pdf"), width = 6, height = 6)
hist(log(dist_Pol_HS+0.01),col=rgb(0,1,1,0.3),probability=TRUE,breaks=40,xlim=c(-5,20),ylim=c(0,0.4),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_denovo_Pol+0.01),col=rgb(1,0,1,0.3),probability=TRUE,breaks=20,add=TRUE)
hist(log(dist_denovo_HS+0.01),col=rgb(1,0.7,0,0.3),probability=TRUE,breaks=40,add=TRUE)
legend('topleft',legend=c('De novo - Polymorphic L1','De novo - Human specific L1','Polymorphic - Human specific L1'),fill=c(rgb(1,0,1,1),rgb(1,0.7,0,1),rgb(0,1,1,1)),bg='white')
dev.off()


# cumulative distribution
ecdf_log_denovo_Pol=ecdf(log(dist_denovo_Pol+0.01))
ecdf_log_denovo_HS=ecdf(log(dist_denovo_HS+0.01))
ecdf_log_Pol_HS=ecdf(log(dist_Pol_HS+0.01))
pdf(paste0(N, "_distance_diff_type_cumulative_log.pdf"), width = 6, height = 6)
plot(ecdf_log_Pol_HS,col=rgb(0,1,1,1),do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_denovo_Pol,col=rgb(1,0,1,1),do.points=FALSE,lwd=2,add=TRUE)
plot(ecdf_log_denovo_HS,col=rgb(1,0.7,0,1),do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('De novo - Polymorphic L1','De novo - Human specific L1','Polymorphic - Human specific L1'),col=c(rgb(1,0,1,1),rgb(1,0.7,0,1),rgb(0,1,1,1)),lty=1,lwd=2,bg='white')
dev.off()
####--------


#### Distance of L1 elements of the same type: distribution vs random ####
shuffleBed <- function(x, genome_file = NULL, excl_file = NULL, chrom = FALSE, maxTries = 1000){
  require(rtracklayer)
  
  stopifnot(class(x) == "GRanges")
  if (nchar(Sys.which("shuffleBed")) == 0) stop("shuffleBed must be in the PATH")
  if(is.null(genome_file)){
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    genome_file <- tempfile()
    writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file)
  }
  bed_file <- tempfile()
  writeLines(paste(as.character(seqnames(x)), start(x)-1, end(x), sep="\t"), bed_file)
  
  shuffle_file <- tempfile()
  to_run="shuffleBed"
  if(!is.null(excl_file))
    to_run <- paste(to_run, "-excl", excl_file)
  if(chrom)
    to_run <- paste(to_run, "-chrom")
  to_run <- paste(to_run, "-maxTries", maxTries, "-i", bed_file, "-g", genome_file, ">", shuffle_file)
  system(to_run)
  shuffled <- import.bed(shuffle_file)
  seqlevels(shuffled) <- seqlevels(x)
  seqlengths(shuffled) <- seqlengths(x)
  
  #clean up
  unlink(c(bed_file, shuffle_file))
  shuffled
}

set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N) 
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                     as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                    as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                 as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)

shuffled <- shuffleBed(L1denovo_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1denovo_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1Pol_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1Pol_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1HS_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1HS_shuffled <- sort(shuffled, ignore.strand = TRUE)

dist_denovo_Pol_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance,
                              as.data.frame(distanceToNearest(L1Pol_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
dist_denovo_HS_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                             as.data.frame(distanceToNearest(L1HS_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
dist_Pol_HS_shuffled <- c(as.data.frame(distanceToNearest(L1Pol_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                          as.data.frame(distanceToNearest(L1HS_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance)


require(Matching)
p.value=ks.boot(log(dist_denovo_Pol+0.01),log(dist_denovo_Pol_shuffled+0.01))$ks$p.value
if(p.value<1e-16)
  p.value='<1e-16'


# histogram
pdf(paste0(N, "_distance_denovo_Pol_histogram_log.pdf"), width = 6, height = 6)
hist(log(dist_denovo_Pol+0.01),col=rgb(1,0,1,0.3),probability=TRUE,breaks=20,xlim=c(-5,20),ylim=c(0,0.4),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_denovo_Pol_shuffled+0.01),col=rgb(0,0,0,0.3),probability=TRUE,breaks=20,add=TRUE)
legend('topleft',legend=c('De novo - Polymorphic L1','Random',paste0('p-value ',format(p.value,digits=2))),fill=c(rgb(1,0,1,1),rgb(0,0,0,1),'white'),border=c('black','black','white'),bg='white')
dev.off()


# cumulative distribution
ecdf_log_denovo_Pol=ecdf(log(dist_denovo_Pol+0.01))
ecdf_log_denovo_Pol_shuffled=ecdf(log(dist_denovo_Pol_shuffled+0.01))
pdf(paste0(N, "_distance_denovo_Pol_cumulative_log.pdf"), width = 6, height = 6)
plot(ecdf_log_denovo_Pol,col=rgb(1,0,1,1),do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_denovo_Pol_shuffled,col='black',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('De novo - Polymorphic L1','Random',paste0('p-value ',format(p.value,digits=2))),col=c(rgb(1,0,1,1),'black','white'),lty=1,lwd=2,bg='white')
dev.off()


# q-q plot
pdf(paste0(N, "_distance_denovo_Pol_qqplot_log.pdf"), width = 6, height = 6)
qqplot(log(dist_denovo_Pol_shuffled+0.01),log(dist_denovo_Pol+0.01),main='Log Distance between L1',xlab='Random',ylab='De novo - Polymorphic L1')
abline(a=0,b=1,col=rgb(1,0,1,1))
legend('bottomright',legend=paste0('p-value ',format(p.value,digits=2)),col='white',lty=1,lwd=2,bg='white')
dev.off()


require(Matching)
p.value=ks.boot(log(dist_denovo_HS+0.01),log(dist_denovo_HS_shuffled+0.01))$ks$p.value
if(p.value<1e-16)
  p.value='<1e-16'


# histogram
pdf(paste0(N, "_distance_denovo_HS_histogram_log.pdf"), width = 6, height = 6)
hist(log(dist_denovo_HS+0.01),col=rgb(1,0.7,0,0.3),probability=TRUE,breaks=20,xlim=c(-5,20),ylim=c(0,0.4),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_denovo_HS_shuffled+0.01),col=rgb(0,0,0,0.3),probability=TRUE,breaks=20,add=TRUE)
legend('topleft',legend=c('De novo - Human specific L1','Random',paste0('p-value ',format(p.value,digits=2))),fill=c(rgb(1,0.7,0,1),rgb(0,0,0,1),'white'),border=c('black','black','white'),bg='white')
dev.off()


# cumulative distribution
ecdf_log_denovo_HS=ecdf(log(dist_denovo_HS+0.01))
ecdf_log_denovo_HS_shuffled=ecdf(log(dist_denovo_HS_shuffled+0.01))
pdf(paste0(N, "_distance_denovo_HS_cumulative_log.pdf"), width = 6, height = 6)
plot(ecdf_log_denovo_HS,col=rgb(1,0.7,0,1),do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_denovo_HS_shuffled,col='black',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('De novo - Human specific L1','Random',paste0('p-value ',format(p.value,digits=2))),col=c(rgb(1,0.7,0,1),'black','white'),lty=1,lwd=2,bg='white')
dev.off()


# q-q plot
pdf(paste0(N, "_distance_denovo_HS_qqplot_log.pdf"), width = 6, height = 6)
qqplot(log(dist_denovo_HS_shuffled+0.01),log(dist_denovo_HS+0.01),main='Log Distance between L1',xlab='Random',ylab='De novo - Human specific L1')
abline(a=0,b=1,col=rgb(1,0.7,0,1))
legend('bottomright',legend=paste0('p-value ',format(p.value,digits=2)),col='white',lty=1,lwd=2,bg='white')
dev.off()


require(Matching)
p.value=ks.boot(log(dist_Pol_HS+0.01),log(dist_Pol_HS_shuffled+0.01))$ks$p.value
if(p.value<1e-16)
  p.value='<1e-16'


# histogram
pdf(paste0(N, "_distance_Pol_HS_histogram_log.pdf"), width = 6, height = 6)
hist(log(dist_Pol_HS+0.01),col=rgb(0,1,1,0.3),probability=TRUE,breaks=40,xlim=c(-5,20),ylim=c(0,0.4),xlab='Log Distance ',main='Distance between L1')
hist(log(dist_Pol_HS_shuffled+0.01),col=rgb(0,0,0,0.3),probability=TRUE,breaks=20,add=TRUE)
legend('topleft',legend=c('Polymorphic - Human specific L1','Random',paste0('p-value ',format(p.value,digits=2))),fill=c(rgb(0,1,1,1),rgb(0,0,0,1),'white'),border=c('black','black','white'),bg='white')
dev.off()


# cumulative distribution
ecdf_log_Pol_HS=ecdf(log(dist_Pol_HS+0.01))
ecdf_log_Pol_HS_shuffled=ecdf(log(dist_Pol_HS_shuffled+0.01))
pdf(paste0(N, "_distance_Pol_HS_cumulative_log.pdf"), width = 6, height = 6)
plot(ecdf_log_Pol_HS,col=rgb(0,1,1,1),do.points=FALSE,lwd=2,xlim=c(-5,20),xlab='Log Distance',main='Distance between L1')
plot(ecdf_log_Pol_HS_shuffled,col='black',do.points=FALSE,lwd=2,add=TRUE)
legend('topleft',legend=c('Polymorphic - Human specific L1','Random',paste0('p-value ',format(p.value,digits=2))),col=c(rgb(0,1,1,1),'black','white'),lty=1,lwd=2,bg='white')
dev.off()


# q-q plot
pdf(paste0(N, "_distance_Pol_HS_qqplot_log.pdf"), width = 6, height = 6)
qqplot(log(dist_Pol_HS_shuffled+0.01),log(dist_Pol_HS+0.01),main='Log Distance between L1',xlab='Random',ylab='Polymorphic - Human specific L1')
abline(a=0,b=1,col=rgb(0,1,1,1))
legend('bottomright',legend=paste0('p-value ',format(p.value,digits=2)),col='white',lty=1,lwd=2,bg='white')
dev.off()
####--------


#### Distance of L1 elements of the same type: distribution NORMALIZED ####
shuffleBed <- function(x, genome_file = NULL, excl_file = NULL, chrom = FALSE, maxTries = 1000){
  require(rtracklayer)
  
  stopifnot(class(x) == "GRanges")
  if (nchar(Sys.which("shuffleBed")) == 0) stop("shuffleBed must be in the PATH")
  if(is.null(genome_file)){
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    genome_file <- tempfile()
    writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file)
  }
  bed_file <- tempfile()
  writeLines(paste(as.character(seqnames(x)), start(x)-1, end(x), sep="\t"), bed_file)
  
  shuffle_file <- tempfile()
  to_run="shuffleBed"
  if(!is.null(excl_file))
    to_run <- paste(to_run, "-excl", excl_file)
  if(chrom)
    to_run <- paste(to_run, "-chrom")
  to_run <- paste(to_run, "-maxTries", maxTries, "-i", bed_file, "-g", genome_file, ">", shuffle_file)
  system(to_run)
  shuffled <- import.bed(shuffle_file)
  seqlevels(shuffled) <- seqlevels(x)
  seqlengths(shuffled) <- seqlengths(x)
  
  #clean up
  unlink(c(bed_file, shuffle_file))
  shuffled
}
decdf <- function(x, a, b)  ecdf(a)(x) - ecdf(b)(x)

set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N) 
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                     as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                    as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                 as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)

shuffled <- shuffleBed(L1denovo_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1denovo_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1Pol_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1Pol_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1HS_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1HS_shuffled <- sort(shuffled, ignore.strand = TRUE)

dist_denovo_Pol_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance,
                              as.data.frame(distanceToNearest(L1Pol_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
dist_denovo_HS_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                             as.data.frame(distanceToNearest(L1HS_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
dist_Pol_HS_shuffled <- c(as.data.frame(distanceToNearest(L1Pol_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                          as.data.frame(distanceToNearest(L1HS_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance)

decdf_log_denovo_Pol = function(x) decdf(x,log(dist_denovo_Pol+0.01),log(dist_denovo_Pol_shuffled+0.01))
decdf_log_denovo_HS = function(x) decdf(x,log(dist_denovo_HS+0.01),log(dist_denovo_HS_shuffled+0.01))
decdf_log_Pol_HS = function(x) decdf(x,log(dist_Pol_HS+0.01),log(dist_Pol_HS_shuffled+0.01))

pdf(paste0(N, "_distance_diff_type_cumulative_log_diff.pdf"), width = 6, height = 6)
curve(decdf_log_denovo_Pol,col=rgb(1,0,1,1),lwd=2,xlim=c(-5,20),ylim=c(-0.2,0.1),xlab='Log Distance',ylab='Fn(x) diff',main='Distance between L1')
curve(decdf_log_denovo_HS,col=rgb(1,0.7,0,1),lwd=2,add=TRUE)
curve(decdf_log_Pol_HS,col=rgb(0,1,1,1),lwd=2,add=TRUE)
legend('topright',legend=c('De novo - Polymorphic L1','De novo - Human specific L1','Polymorphic - Human specific L1'),col=c(rgb(1,0,1,1),rgb(1,0.7,0,1),rgb(0,1,1,1)),lty=1,lwd=2,bg='white')
dev.off()

####--------


#### SAME SAMPLE SIZE - MANY TIMES 

# Maximum size
min( c(length(L1denovo), length(L1Pol), length(L1HS)) )

# SAMPLE SIZE 900
N <- 900

#### Distance of L1 elements of different types: distribution ####
set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N)
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                     as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                    as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                 as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)

ecdf_log_denovo_Pol <- ecdf(log(dist_denovo_Pol+0.01))
ecdf_log_denovo_HS <- ecdf(log(dist_denovo_HS+0.01))
ecdf_log_Pol_HS <- ecdf(log(dist_Pol_HS+0.01))

pdf(paste0(N, "_distance_diff_type_cumulative_log_many_times.pdf"), width = 6, height = 6)
plot(ecdf_log_Pol_HS, col = rgb(0,1,1,1), do.points = FALSE, lwd = 2, xlim = c(-5,20), xlab = 'Log Distance', main = 'Distance between L1')
plot(ecdf_log_denovo_Pol, col = rgb(1,0,1,1), do.points = FALSE, lwd = 2, add = TRUE)
plot(ecdf_log_denovo_HS, col = rgb(1,0.7,0,1), do.points = FALSE, lwd = 2, add = TRUE)

for(i in 1:100){
  index_denovo <- sample(seq_along(L1denovo), N)
  index_Pol <- sample(seq_along(L1Pol), N)
  index_HS <- sample(seq_along(L1HS), N)
  L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
  L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
  L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)
  
  
  dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                       as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
  dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                      as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
  dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                   as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)
  
  ecdf_log_denovo_Pol <- ecdf(log(dist_denovo_Pol+0.01))
  ecdf_log_denovo_HS <- ecdf(log(dist_denovo_HS+0.01))
  ecdf_log_Pol_HS <- ecdf(log(dist_Pol_HS+0.01))
  
  plot(ecdf_log_Pol_HS, col = rgb(0,1,1,1), do.points = FALSE, lwd = 2, add = TRUE)
  plot(ecdf_log_denovo_Pol, col = rgb(1,0,1,1), do.points = FALSE, lwd = 2, add = TRUE)
  plot(ecdf_log_denovo_HS, col = rgb(1,0.7,0,1), do.points = FALSE, lwd = 2, add = TRUE)
}
legend('topleft',legend=c('De novo - Polymorphic L1', 'De novo - Human specific L1', 'Polymorphic - Human specific L1'),
       col = c(rgb(1,0,1,1), rgb(1,0.7,0,1), rgb(0,1,1,1)), lty = 1, lwd = 2, bg='white')
dev.off()
###-------


#### Distance of L1 elements of different types: distribution NORMALIZED ####
shuffleBed <- function(x, genome_file = NULL, excl_file = NULL, chrom = FALSE, maxTries = 1000){
  require(rtracklayer)
  
  stopifnot(class(x) == "GRanges")
  if (nchar(Sys.which("shuffleBed")) == 0) stop("shuffleBed must be in the PATH")
  if(is.null(genome_file)){
    if (any(is.na(seqlengths(x)))) stop("seqlengths of x must be defined")
    genome_file <- tempfile()
    writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file)
  }
  bed_file <- tempfile()
  writeLines(paste(as.character(seqnames(x)), start(x)-1, end(x), sep="\t"), bed_file)
  
  shuffle_file <- tempfile()
  to_run="shuffleBed"
  if(!is.null(excl_file))
    to_run <- paste(to_run, "-excl", excl_file)
  if(chrom)
    to_run <- paste(to_run, "-chrom")
  to_run <- paste(to_run, "-maxTries", maxTries, "-i", bed_file, "-g", genome_file, ">", shuffle_file)
  system(to_run)
  shuffled <- import.bed(shuffle_file)
  seqlevels(shuffled) <- seqlevels(x)
  seqlengths(shuffled) <- seqlengths(x)
  
  #clean up
  unlink(c(bed_file, shuffle_file))
  shuffled
}
decdf <- function(x, a, b)  ecdf(a)(x) - ecdf(b)(x)

set.seed(2307)
index_denovo <- sample(seq_along(L1denovo), N)
index_Pol <- sample(seq_along(L1Pol), N)
index_HS <- sample(seq_along(L1HS), N) 
L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)

dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                     as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                    as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                 as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)

shuffled <- shuffleBed(L1denovo_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1denovo_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1Pol_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1Pol_shuffled <- sort(shuffled, ignore.strand = TRUE)
shuffled <- shuffleBed(L1HS_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
L1HS_shuffled <- sort(shuffled, ignore.strand = TRUE)

dist_denovo_Pol_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance,
                              as.data.frame(distanceToNearest(L1Pol_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
dist_denovo_HS_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                             as.data.frame(distanceToNearest(L1HS_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
dist_Pol_HS_shuffled <- c(as.data.frame(distanceToNearest(L1Pol_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                          as.data.frame(distanceToNearest(L1HS_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance)

decdf_log_denovo_Pol = function(x) decdf(x,log(dist_denovo_Pol+0.01),log(dist_denovo_Pol_shuffled+0.01))
decdf_log_denovo_HS = function(x) decdf(x,log(dist_denovo_HS+0.01),log(dist_denovo_HS_shuffled+0.01))
decdf_log_Pol_HS = function(x) decdf(x,log(dist_Pol_HS+0.01),log(dist_Pol_HS_shuffled+0.01))

pdf(paste0(N, "_distance_diff_type_cumulative_log_diff_many_times.pdf"), width = 6, height = 6)
curve(decdf_log_denovo_Pol, col = rgb(1,0,1,1), lwd = 2, xlim = c(-5,20), ylim = c(-0.2,0.1), xlab = 'Log Distance', ylab = 'Fn(x) diff', main = 'Distance between L1')
curve(decdf_log_denovo_HS, col = rgb(1,0.7,0,1), lwd = 2, add = TRUE)
curve(decdf_log_Pol_HS, col = rgb(0,1,1,1), lwd = 2, add = TRUE)

for(i in 1:100){
  index_denovo <- sample(seq_along(L1denovo), N)
  index_Pol <- sample(seq_along(L1Pol), N)
  index_HS <- sample(seq_along(L1HS), N)
  L1denovo_sample <- sort(L1denovo[index_denovo], ignore.strand = TRUE)
  L1Pol_sample <- sort(L1Pol[index_Pol], ignore.strand = TRUE)
  L1HS_sample <- sort(L1HS[index_HS], ignore.strand = TRUE)
  
  dist_denovo_Pol <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1Pol_sample, ignore.strand = TRUE))$distance,
                       as.data.frame(distanceToNearest(L1Pol_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
  dist_denovo_HS <- c(as.data.frame(distanceToNearest(L1denovo_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                      as.data.frame(distanceToNearest(L1HS_sample, L1denovo_sample, ignore.strand = TRUE))$distance)
  dist_Pol_HS <- c(as.data.frame(distanceToNearest(L1Pol_sample, L1HS_sample, ignore.strand = TRUE))$distance,
                   as.data.frame(distanceToNearest(L1HS_sample, L1Pol_sample, ignore.strand = TRUE))$distance)
  
  shuffled <- shuffleBed(L1denovo_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
  L1denovo_shuffled <- sort(shuffled, ignore.strand = TRUE)
  shuffled <- shuffleBed(L1Pol_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
  L1Pol_shuffled<- sort(shuffled, ignore.strand = TRUE)
  shuffled <-shuffleBed(L1HS_sample, genome_file = 'chromsizes_hg19_noY.tab', excl_file = 'gaps.bed', chrom = FALSE)
  L1HS_shuffled <- sort(shuffled, ignore.strand = TRUE)
  
  dist_denovo_Pol_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance,
                                as.data.frame(distanceToNearest(L1Pol_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
  dist_denovo_HS_shuffled <- c(as.data.frame(distanceToNearest(L1denovo_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                               as.data.frame(distanceToNearest(L1HS_shuffled, L1denovo_shuffled, ignore.strand = TRUE))$distance)
  dist_Pol_HS_shuffled <- c(as.data.frame(distanceToNearest(L1Pol_shuffled, L1HS_shuffled, ignore.strand = TRUE))$distance,
                            as.data.frame(distanceToNearest(L1HS_shuffled, L1Pol_shuffled, ignore.strand = TRUE))$distance)
  
  decdf_log_denovo_Pol = function(x) decdf(x,log(dist_denovo_Pol+0.01),log(dist_denovo_Pol_shuffled+0.01))
  decdf_log_denovo_HS = function(x) decdf(x,log(dist_denovo_HS+0.01),log(dist_denovo_HS_shuffled+0.01))
  decdf_log_Pol_HS = function(x) decdf(x,log(dist_Pol_HS+0.01),log(dist_Pol_HS_shuffled+0.01))
  
  curve(decdf_log_denovo_Pol, col = rgb(1,0,1,1), lwd = 2, add = TRUE)
  curve(decdf_log_denovo_HS, col = rgb(1,0.7,0,1), lwd = 2, add = TRUE)
  curve(decdf_log_Pol_HS, col = rgb(0,1,1,1), lwd = 2, add = TRUE)
}
legend('bottomleft', legend=c('De novo - Polymorphic L1', 'De novo - Human specific L1', 'Polymorphic - Human specific L1'),
       col = c(rgb(1,0,1,1), rgb(1,0.7,0,1), rgb(0,1,1,1)), lty = 1, lwd = 2, bg = 'white')
dev.off()
####--------