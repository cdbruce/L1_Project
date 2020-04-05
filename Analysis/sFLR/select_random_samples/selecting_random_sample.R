
setwd('~/Google_Drive/L1/L1_Project/Analysis/sFLR/select_random_samples/')
sum_comp1_randoms<-as.data.frame(matrix())
sum_comp4_randoms<-as.data.frame(matrix())
sum_comp5_randoms<-as.data.frame(matrix())

for (r in 1:10){
  r=1#Since random 1 was selected in our analysis, here we only load random 1 for illustration purpose
  setwd(paste("~/Google_Drive/L1/L1_Project/Analysis/sFLR/random",r,sep=''))
  r_comp1_scal<-read.csv(file='comp1_scalar_logit_sum.csv')
  r_comp1_func<-read.csv(file='comp1_func_logit_sum.csv')
  r_comp1<-rbind(r_comp1_scal[,c(3,4,7)],r_comp1_func[,c(3,4,6)])
  r_comp1_ordered<- r_comp1[order(r_comp1$feature_num),] 
  
  r_comp4_scal<-read.csv(file='comp4_scalar_logit_sum.csv')
  r_comp4_func<-read.csv(file='comp4_func_logit_sum.csv')
  r_comp4<-rbind(r_comp4_scal[,c(3,4,7)],r_comp4_func[,c(3,4,6)])
  r_comp4_ordered<- r_comp4[order(r_comp4$feature_num),] 
  
  r_comp5_scal<-read.csv(file='comp5_scalar_logit_sum.csv')
  r_comp5_func<-read.csv(file='comp5_func_logit_sum.csv')
  r_comp5<-rbind(r_comp5_scal[,c(3,4,7)],r_comp5_func[,c(3,4,6)])
  r_comp5_ordered<- r_comp5[order(r_comp5$feature_num),] 
  
  sum_comp1_randoms<-cbind(sum_comp1_randoms,r_comp1_ordered[,3])
  sum_comp4_randoms<-cbind(sum_comp4_randoms,r_comp4_ordered[,3])
  sum_comp5_randoms<-cbind(sum_comp5_randoms,r_comp5_ordered[,3])
}

sum_comp1_randoms<-sum_comp1_randoms[,-1]
sum_comp4_randoms<-sum_comp4_randoms[,-1]
sum_comp5_randoms<-sum_comp5_randoms[,-1]

sum_comp1_randoms<-cbind(r_comp1_ordered[,c(1,2)],sum_comp1_randoms,rowMeans(sum_comp1_randoms))
sum_comp4_randoms<-cbind(r_comp4_ordered[,c(1,2)],sum_comp4_randoms,rowMeans(sum_comp4_randoms))
sum_comp5_randoms<-cbind(r_comp5_ordered[,c(1,2)],sum_comp5_randoms,rowMeans(sum_comp5_randoms))
# Order (descending) by psuedo-R2
# sum_comp1_randoms<-sum_comp1_randoms[order(-sum_comp1_randoms[,13]),]
# sum_comp4_randoms<-sum_comp4_randoms[order(-sum_comp4_randoms[,13]),]
# sum_comp5_randoms<-sum_comp5_randoms[order(-sum_comp5_randoms[,13]),]

# # Take the features with the top 10 psuedo-R2
# sum_comp1_randoms<-sum_comp1_randoms[c(1:10),]
# sum_comp4_randoms<-sum_comp4_randoms[c(1:10),]
# sum_comp5_randoms<-sum_comp5_randoms[c(1:10),]

# 
# plot(sum_comp1_randoms[,3],type="b")
# plot(sum_comp1_randoms[,4],type="b")
# plot(sum_comp1_randoms[,12],type="b")


####Plot the pseudo-R2 for three comparisons involving de novo L1s 
setwd('~/Desktop/regression/10_randoms/select_random_samples/')

#Comp1 (L1denovo vs control)
pdf('comparisons_pseudo-R2_10randoms.pdf', width=15,height = 10)

# Unordered version
# pdf('comparisons_pseudo-R2_10randoms_unordered.pdf', width=15,height = 10)
# Top 10 pseudo-R2's
# pdf('comparisons_pseudo-R2_10randoms_top10.pdf', width=15,height = 10)

plot(sum_comp1_randoms[,3],type='n',main="Comparison of pseudo-R2 across 10 random samples",
     sub="de novo L1s vs control", xlab='',ylab='pseudo-R2',ylim=c(0,0.3))

points(sum_comp1_randoms[,3],cex=0.4, pch='1')
lines(sum_comp1_randoms[,3],lty=1,col = "red1")

points(sum_comp1_randoms[,4],cex=0.4,pch='2')
lines(sum_comp1_randoms[,4],lty=1,col = "blue1")

points(sum_comp1_randoms[,5],cex=0.4,pch='3')
lines(sum_comp1_randoms[,5],lty=1,col = "green1")

points(sum_comp1_randoms[,6],cex=0.4,pch='4')
lines(sum_comp1_randoms[,6],lty=1,col = "pink1")

points(sum_comp1_randoms[,7],cex=0.4,pch='5')
lines(sum_comp1_randoms[,7],lty=1,col = "purple1")

points(sum_comp1_randoms[,8],cex=0.4,pch='6')
lines(sum_comp1_randoms[,8],lty=1,col = "brown1")

points(sum_comp1_randoms[,9],cex=0.4,pch='7')
lines(sum_comp1_randoms[,9],lty=1,col = "orange1")

points(sum_comp1_randoms[,10],cex=0.4,pch='8')
lines(sum_comp1_randoms[,10],lty=1,col = "salmon1")

points(sum_comp1_randoms[,11],cex=0.4,pch='9')
lines(sum_comp1_randoms[,11],lty=1,col = "tan1")

points(sum_comp1_randoms[,12],cex=0.4,pch='0')
lines(sum_comp1_randoms[,12],lty=1,col = "khaki")

#dev.off()

#Comp4 (L1deonvo vs L1Pol)
#pdf('comp4_pseudo-R2_10randoms.pdf', width=15,height = 10)
plot(sum_comp4_randoms[,3],type='n',main="Comparison of pseudo-R2 across 10 random samples",
     sub="de novo L1s vs Polymorphic L1s", xlab='',ylab='pseudo-R2',ylim=c(0,0.3))

points(sum_comp4_randoms[,3],cex=0.4, pch='1')
lines(sum_comp4_randoms[,3],lty=1,col = "red1")

points(sum_comp4_randoms[,4],cex=0.4,pch='2')
lines(sum_comp4_randoms[,4],lty=1,col = "blue1")

points(sum_comp4_randoms[,5],cex=0.4,pch='3')
lines(sum_comp4_randoms[,5],lty=1,col = "green1")

points(sum_comp4_randoms[,6],cex=0.4,pch='4')
lines(sum_comp4_randoms[,6],lty=1,col = "pink1")

points(sum_comp4_randoms[,7],cex=0.4,pch='5')
lines(sum_comp4_randoms[,7],lty=1,col = "purple1")

points(sum_comp4_randoms[,8],cex=0.4,pch='6')
lines(sum_comp4_randoms[,8],lty=1,col = "brown1")

points(sum_comp4_randoms[,9],cex=0.4,pch='7')
lines(sum_comp4_randoms[,9],lty=1,col = "orange1")

points(sum_comp4_randoms[,10],cex=0.4,pch='8')
lines(sum_comp4_randoms[,10],lty=1,col = "salmon1")

points(sum_comp4_randoms[,11],cex=0.4,pch='9')
lines(sum_comp4_randoms[,11],lty=1,col = "tan1")

points(sum_comp4_randoms[,12],cex=0.4,pch='0')
lines(sum_comp4_randoms[,12],lty=1,col = "khaki")

#dev.off()

## Comp5 (L1denovo vs L1HS)
#pdf('comp5_pseudo-R2_10randoms.pdf', width=15,height = 10)
plot(sum_comp5_randoms[,3],type='n',main="Comparison of pseudo-R2 across 10 random samples",
     sub="de novo L1s vs Human-specific L1s", xlab='',ylab='pseudo-R2',ylim=c(0,0.3))

points(sum_comp5_randoms[,3],cex=0.4, pch='1')
lines(sum_comp5_randoms[,3],lty=1,col = "red1")

points(sum_comp5_randoms[,4],cex=0.4,pch='2')
lines(sum_comp5_randoms[,4],lty=1,col = "blue1")

points(sum_comp5_randoms[,5],cex=0.4,pch='3')
lines(sum_comp5_randoms[,5],lty=1,col = "green1")

points(sum_comp5_randoms[,6],cex=0.4,pch='4')
lines(sum_comp5_randoms[,6],lty=1,col = "pink1")

points(sum_comp5_randoms[,7],cex=0.4,pch='5')
lines(sum_comp5_randoms[,7],lty=1,col = "purple1")

points(sum_comp5_randoms[,8],cex=0.4,pch='6')
lines(sum_comp5_randoms[,8],lty=1,col = "brown1")

points(sum_comp5_randoms[,9],cex=0.4,pch='7')
lines(sum_comp5_randoms[,9],lty=1,col = "orange1")

points(sum_comp5_randoms[,10],cex=0.4,pch='8')
lines(sum_comp5_randoms[,10],lty=1,col = "salmon1")

points(sum_comp5_randoms[,11],cex=0.4,pch='9')
lines(sum_comp5_randoms[,11],lty=1,col = "tan1")

points(sum_comp5_randoms[,12],cex=0.4,pch='0')
lines(sum_comp5_randoms[,12],lty=1,col = "khaki")

dev.off()



############
## Color random1 and 6 differently, while fading other 8 random samples
## Color random1 differently, while fading other 9 random samples

#pdf('comparisons_pseudo-R2_r1_r6.pdf', width=15,height = 10)
#pdf('comparisons_pseudo-R2_r1_r10.pdf', width=15,height = 10)
pdf('comparisons_pseudo-R2_r1.pdf', width=15,height = 10)

# Unordered version
# pdf('comparisons_pseudo-R2_10randoms_unordered.pdf', width=15,height = 10)
# Top 10 pseudo-R2's
# pdf('comparisons_pseudo-R2_10randoms_top10.pdf', width=15,height = 10)
a=sum_comp1_randoms$feature_names
a<-as.vector(a)

plot(sum_comp1_randoms[,3],type='n',main="Comparison of pseudo-R2 across 10 random samples\nde novo L1s vs Control",
     xlab='',ylab='pseudo-R2',ylim=c(0,0.3),xaxt='n')
axis(side=1,at=1:45,label=a,las = 2,cex.axis=0.5)

points(sum_comp1_randoms[,3],cex=0.4, pch='1')
lines(sum_comp1_randoms[,3],lty=1,col = "red1")

points(sum_comp1_randoms[,4],cex=0.4,pch='2')
lines(sum_comp1_randoms[,4],lty=1,col = "lightgray")

points(sum_comp1_randoms[,5],cex=0.4,pch='3')
lines(sum_comp1_randoms[,5],lty=1,col = "lightgray")

points(sum_comp1_randoms[,6],cex=0.4,pch='4')
lines(sum_comp1_randoms[,6],lty=1,col = "lightgray")

points(sum_comp1_randoms[,7],cex=0.4,pch='5')
lines(sum_comp1_randoms[,7],lty=1,col = "lightgray")

points(sum_comp1_randoms[,8],cex=0.4,pch='6')
#lines(sum_comp1_randoms[,8],lty=1,col = "green1")
lines(sum_comp1_randoms[,8],lty=1,col = "lightgray")

points(sum_comp1_randoms[,9],cex=0.4,pch='7')
lines(sum_comp1_randoms[,9],lty=1,col = "lightgray")

points(sum_comp1_randoms[,10],cex=0.4,pch='8')
lines(sum_comp1_randoms[,10],lty=1,col = "lightgray")

points(sum_comp1_randoms[,11],cex=0.4,pch='9')
lines(sum_comp1_randoms[,11],lty=1,col = "lightgray")

points(sum_comp1_randoms[,12],cex=0.4,pch='0')
lines(sum_comp1_randoms[,12],lty=1,col = "lightgray")
#lines(sum_comp1_randoms[,12],lty=1,col = "blue")

#dev.off()

#Comp4 (L1deonvo vs L1Pol)
#pdf('comp4_pseudo-R2_10randoms.pdf', width=15,height = 10)
plot(sum_comp4_randoms[,3],type='n',main="Comparison of pseudo-R2 across 10 random samples\nde novo L1s vs Polymorphic L1s",
     xlab='',ylab='pseudo-R2',ylim=c(0,0.3),xaxt='n')
axis(side=1,at=1:45,label=a,las = 2,cex.axis=0.5)


points(sum_comp4_randoms[,3],cex=0.4, pch='1')
lines(sum_comp4_randoms[,3],lty=1,col = "red1")

points(sum_comp4_randoms[,4],cex=0.4,pch='2')
lines(sum_comp4_randoms[,4],lty=1,col = "lightgray")

points(sum_comp4_randoms[,5],cex=0.4,pch='3')
lines(sum_comp4_randoms[,5],lty=1,col = "lightgray")

points(sum_comp4_randoms[,6],cex=0.4,pch='4')
lines(sum_comp4_randoms[,6],lty=1,col = "lightgray")

points(sum_comp4_randoms[,7],cex=0.4,pch='5')
lines(sum_comp4_randoms[,7],lty=1,col = "lightgray")

points(sum_comp4_randoms[,8],cex=0.4,pch='6')
#lines(sum_comp4_randoms[,8],lty=1,col = "green1")
lines(sum_comp4_randoms[,8],lty=1,col = "lightgray")

points(sum_comp4_randoms[,9],cex=0.4,pch='7')
lines(sum_comp4_randoms[,9],lty=1,col = "lightgray")

points(sum_comp4_randoms[,10],cex=0.4,pch='8')
lines(sum_comp4_randoms[,10],lty=1,col = "lightgray")

points(sum_comp4_randoms[,11],cex=0.4,pch='9')
lines(sum_comp4_randoms[,11],lty=1,col = "lightgray")

points(sum_comp4_randoms[,12],cex=0.4,pch='0')
lines(sum_comp4_randoms[,12],lty=1,col = "lightgray")
#lines(sum_comp4_randoms[,12],lty=1,col = "blue")

#dev.off()

## Comp5 (L1denovo vs L1HS)
#pdf('comp5_pseudo-R2_10randoms.pdf', width=15,height = 10)
plot(sum_comp5_randoms[,3],type='n',main="Comparison of pseudo-R2 across 10 random samples\nde novo L1s vs Human-specific L1s",
     xlab='',ylab='pseudo-R2',ylim=c(0,0.3),xaxt='n')
axis(side=1,at=1:45,label=a,las = 2,cex.axis=0.5)

points(sum_comp5_randoms[,3],cex=0.4, pch='1')
lines(sum_comp5_randoms[,3],lty=1,col = "red1")

points(sum_comp5_randoms[,4],cex=0.4,pch='2')
lines(sum_comp5_randoms[,4],lty=1,col = "lightgray")

points(sum_comp5_randoms[,5],cex=0.4,pch='3')
lines(sum_comp5_randoms[,5],lty=1,col = "lightgray")

points(sum_comp5_randoms[,6],cex=0.4,pch='4')
lines(sum_comp5_randoms[,6],lty=1,col = "lightgray")

points(sum_comp5_randoms[,7],cex=0.4,pch='5')
lines(sum_comp5_randoms[,7],lty=1,col = "lightgray")

points(sum_comp5_randoms[,8],cex=0.4,pch='6')
#lines(sum_comp5_randoms[,8],lty=1,col = "green1")
lines(sum_comp5_randoms[,8],lty=1,col = "lightgray")

points(sum_comp5_randoms[,9],cex=0.4,pch='7')
lines(sum_comp5_randoms[,9],lty=1,col = "lightgray")

points(sum_comp5_randoms[,10],cex=0.4,pch='8')
lines(sum_comp5_randoms[,10],lty=1,col = "lightgray")

points(sum_comp5_randoms[,11],cex=0.4,pch='9')
lines(sum_comp5_randoms[,11],lty=1,col = "lightgray")

points(sum_comp5_randoms[,12],cex=0.4,pch='0')
lines(sum_comp5_randoms[,12],lty=1,col = "lightgray")
#lines(sum_comp5_randoms[,12],lty=1,col = "blue")


dev.off()

##Save the image
save.image("Select_random1.RData")

##################################
#######Selected Random1###########
##################################

#Plot all the peudoR2's for random1
# r=1
# setwd(paste("~/Desktop/regression/10_randoms/random",r,sep=''))
# sum_comp1_randoms<-as.data.frame(matrix())
# sum_comp3_randoms<-as.data.frame(matrix())
# 
# r_comp1_scal<-read.csv(file='comp1_scalar_logit_sum.csv')
# r_comp1_func<-read.csv(file='comp1_func_logit_sum.csv')
# r_comp1<-rbind(r_comp1_scal[,c(2,3,6)],r_comp1_func[,c(2,3,5)])
# r_comp1_ordered<- r_comp1[order(r_comp1$feature_num),] 
# 
# r_comp3_scal<-read.csv(file='comp3_scalar_logit_sum.csv')
# r_comp3_func<-read.csv(file='comp3_func_logit_sum.csv')
# r_comp3<-rbind(r_comp3_scal[,c(2,3,6)],r_comp3_func[,c(2,3,5)])
# r_comp3_ordered<- r_comp3[order(r_comp3$feature_num),] 
# 
# # Optional: Extract the same feature names as it was in IWT results, not applied
# # in the poster, since some feature names exceed the plot margin
# #r=1
# #setwd(paste("~/Desktop/regression/10_randoms/random",r,sep=''))
# #load(file = paste("L1_transformed_random_",r,'.RData',sep=''))
# #feature_names<-as.vector(result_mean@metadata$feature_datasets$name)
# # sum_comp1_randoms<-cbind(sum_comp1_randoms,r_comp1_ordered[,1],feature_names,r_comp1_ordered[,3])
# # sum_comp3_randoms<-cbind(sum_comp3_randoms,r_comp3_ordered[,1],feature_names,r_comp3_ordered[,3])
# 
# sum_comp1_randoms<-cbind(sum_comp1_randoms,r_comp1_ordered)
# sum_comp3_randoms<-cbind(sum_comp3_randoms,r_comp3_ordered)
# 
# sum_comp1_randoms<-sum_comp1_randoms[,-1]
# sum_comp3_randoms<-sum_comp3_randoms[,-1]
# 
# 
# # Order (descending) by psuedo-R2
# sum_comp1_randoms<-sum_comp1_randoms[order(-sum_comp1_randoms[,3]),]
# # sum_comp4_randoms<-sum_comp4_randoms[order(-sum_comp4_randoms[,13]),]
# sum_comp3_randoms<-sum_comp3_randoms[order(-sum_comp3_randoms[,3]),]
# 
# # Create the pdf file in the same directory
# #pdf('PsudoR2_Comp1_Comp3.pdf',height = 5, width = 10)
# # #Comp1
# pdf('PsudoR2_Comp1.pdf',height = 5, width = 10)
# a=sum_comp1_randoms$feature_names
# #a<-as.character(a)
# #a<-factor(a,levels = a)
# a<-as.vector(a)
# b<-as.vector(sum_comp1_randoms[,3])
# plot(b,type='n',main="Variances explained by individual regressions\nde novo L1 vs Control", 
#      xlab="",ylab='pseudo R2', cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
# axis(side=1,at=1:45,label=a,las = 2,font.axis=1.5,cex.axis=0.3)
# points(b,cex=0.4)
# dev.off()
# 
# #Comp3
# pdf('PsudoR2_Comp3.pdf',height = 5, width = 10)
# a=sum_comp3_randoms$feature_names
# a<-as.vector(a)
# #A<-factor(a,levels = a)
# b<-as.vector(sum_comp3_randoms[,3])
# plot(b,type='n',main="Variances explained by individual regressions\nHuman-specific L1 vs Control", 
#      xlab="",ylab='pseudo R2',cex.lab=0.9,ylim=c(0,0.3),xaxt='n')
# axis(side=1,at=1:45,label=a,las = 3,font.axis=1.5,cex.axis=0.5)
# points(b,cex=0.4)
# dev.off()




