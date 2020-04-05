setwd('/Users/Bruce/Desktop/L1/IWT/cleanControl_pvalues/median/')
load('L1_autosomes_results_smoothed_mean_alldennovo_median.RData')

setwd('Reorder_Final_L1Paper_Sept2019/')

## Set scales for each feature in each test




## Check list/order of features
idFeaturesTest(result_mean)
class(idFeaturesTest(result_mean))
write(idFeaturesTest(result_mean), file='feature_order_original.txt')

## Change order manually in feature_order_original.txt based on Table 1 and read back in
features_reordered<-read.table("feature_order_final.txt", stringsAsFactors=FALSE)
class(features_reordered)
features_reordered<-as.character(features_reordered$V1)
features_reordered


####### ORIGINAL SCALE (FEB 2018) ########################
# scale_threshold=list(test1=c(100,100,10,4,100,
#                              16,16,100,4,100,
#                              8,100,4,4,10,
#                              100,6,100,100,100,
#                              100,4,100,10,10,
#                              2,100,100,100,100,
#                              6,100,100,20,100,
#                              20,2,8,100,8,
#                              100,100,100,100,10),
#                      test2=c(10,100,100,100,100,
#                              100,100,100,100,100,
#                              8,100,100,100,100,
#                              100,8,12,4,4,
#                              100,100,100,100,100,
#                              100,2,100,100,100,
#                              100,4,4,6,100,
#                              2,100,4,6,100,
#                              100,100,2,100,100),
#                      test3=c(10,100,100,100,100,
#                              100,100,100,100,100,
#                              100,100,100,100,100,
#                              100,8,100,6,100,
#                              100,100,100,100,100,
#                              100,6,100,100,100,
#                              100,10,100,100,100,
#                              2,100,10,100,100,
#                              100,100,6,100,100),
#                      test4=c(100,100,100,100,100,
#                              100,100,100,100,4,
#                              8,100,100,100,100,
#                              100,100,100,100,2,
#                              100,100,100,100,100,
#                              100,100,100,100,6,
#                              100,100,100,100,100,
#                              100,100,6,6,100,
#                              8,100,100,100,6),
#                      test5=c(6,100,100,100,100,
#                              100,100,100,100,4,
#                              8,100,100,100,100,
#                              100,4,8,44,6,
#                              100,100,100,100,100,
#                              100,100,100,100,6,
#                              100,100,100,100,100,
#                              10,100,4,100,100,
#                              100,100,100,100,2),
#                      test6=c(24,8,30,4,100,
#                              8,8,10,8,2,
#                              10,2,10,2,100,
#                              100,2,6,6,100,
#                              2,2,100,100,100,
#                              100,100,100,6,6,
#                              100,100,100,6,8,
#                              4,20,4,10,16,
#                              100,100,100,100,100)
# )
# write.table(scale_threshold,file="scales_all_six_tests_original.txt", sep='\t')

## Change the scales for re-ordered features
scale_threshold=list(test1=c(4,4,100,16,4,	
                             16,100,100,10,100,	
                             4,100,8,100,8,	
                             20,10,100,100,100,	
                             10,100,6,100,100,	
                             100,2,100,10,100,	
                             100,100,100,100,6,	
                             2,8,100,20,4,	
                             100),
                     
                     test2=c(100,100,100,100,100,	
                             100,100,100,100,100,	
                             100,100,8,10,100,	
                             2,100,100,2,100,	
                             100,100,8,12,4,	
                             4,100,2,100,4,	
                             100,100,4,100,100,	
                             100,4,100,6,100,	
                             100),
                    
                     test3=c(100,100,100,100,100,	
                             100,100,100,100,100,	
                             100,100,100,10,100,	
                             2,100,100,6,100,	
                             100,100,8,100,6,	
                             100,100,6,100,10,	
                             100,100,100,100,100,	
                             100,10,100,100,100,	
                             100),
                     
                     test4=c(100,100,100,100,100,	
                             100,100,100,100,100,	
                             100,4,8,100,100,	
                             100,6,100,100,100,	
                             100,100,100,100,100,	
                             2,100,100,100,100,	
                             100,6,100,100,100,	
                             100,6,100,100,100,	
                             100),
                     
                     test5=c(100,100,100,100,100,	
                             100,100,100,100,100,	
                             100,4,8,6,100,	
                             10,2,100,100,100,	
                             100,100,4,8,44,	
                             6,100,100,100,100,	
                             100,6,100,100,100,	
                             100,4,100,100,100,	
                             100),
                     
                     test6=c(10,2,2,8,8,	
                             8,10,8,30,100,	
                             4,2,10,24,16,	
                             4,100,8,100,100,	
                             100,100,2,6,6,	
                             100,100,100,100,100,	
                             6,6,100,100,100,	
                             20,4,2,6,2,	
                             100)
)
write.table(scale_threshold,file="scales_all_six_tests_reordered.txt", sep='\t')


## Change direction of the tests
result_mean@test$input$id_region1=c("L1denovo","L1Pol","L1HS","L1Pol","L1HS","L1HS")
result_mean@test$input$id_region2=c("Control","Control","Control","L1denovo","L1denovo","L1Pol")
#Change sign of the test result when swapping the comparison direction
#result_mean@test$result[[i]][[j]]$T0_plot, i=1:6,j=1:45
for (i in 1:45){
  result_mean@test$result[[4]][[i]]$T0_plot=-result_mean@test$result[[4]][[i]]$T0_plot
  result_mean@test$result[[5]][[i]]$T0_plot=-result_mean@test$result[[5]][[i]]$T0_plot
  result_mean@test$result[[6]][[i]]$T0_plot=-result_mean@test$result[[6]][[i]]$T0_plot
  
}

# Plot the heatmap

plotSummary(result_mean,groupby="test",only_significant=FALSE,xlab='kb',
            scale_threshold=scale_threshold,
            # change order of features based on Table 1 ##
            id_features_subset=features_reordered,
            ##
            filenames=paste0("IWT_autosomes_smoothed_mean_summary_test_",c("denovo_control","pol_control","hs_control","pol_denovo","hs_denovo","hs_pol"),"_scaled.pdf"),
            align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)

# plotSummary(result_mean,groupby="feature",only_significant=FALSE,gaps_tests=3,xlab='kb',scale_threshold=scale_threshold,
#             filenames=paste0("IWT_autosomes_smoothed_mean_summary_feature_",idFeatures(result_mean),"_scaled.pdf"),
#             align_lab="Integration site",ask=FALSE,cellwidth=10,cellheight=15)



