setwd("/Users/Bruce/Google_Drive/L1/L1_Manuscript/L1_Draft/MBE_Revision/Analysis/Restriction_sites/")

#Read distance values (in bp) into R 
distance_MSPI_TaqI <- read.table("distance_MSPI_TaqI.bed", sep="\t", header=FALSE)
dim(distance_MSPI_TaqI)
distance <- distance_MSPI_TaqI[,13]
head(distance)

#Plot the distribution
#Violin plot of coverage distribution 
#install.packages("vioplot")
library(vioplot)
vioplot(log(distance+0.0001), names=c("Log Distance Distribution Between MSPI and TaqI Sites"), ylab="Log Distance")
abline(h=log(6000+0.0001), col="green")
#ylim = c(0,0.02))

#Histogram of coverage distribution 
#hist(distance, xlim=c(0,10000))
hist(log(distance+0.0001), main = "Log Distance Distribution Between MSPI and TaqI Sites", xlab="Log Distance")
abline(v=log(6000+0.0001), col="green")

#boxplot(distance)
#plot(density(distance, main="Distance between MSPI and TaqI sites", col="blue"))

#Calculate percentage of distance values under 6k (length of FL-L1)
length(distance)
distance_less_6k <- which(distance<=6000)
length(distance_less_6k)
length(distance_less_6k)/length(distance)
# =0.9998751ß
     