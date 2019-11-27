In general, there are 10 folders (random 1-10) and 4 Rscript files:
transformation_all_dataset.R
func_log_regression_individual.R
glm_scalar_no_lowFeatures.R
glm_lowFeatures.R

transformation_all_dataset: R pipeline to transform all the predictors before performing any scalar/functional regression analysis. Transformation was performed for high-resolution (both functional&scalar) predictors, and low-resolution features separately as two files, generating two output files: “L1_transformed_random_r” (r=1:10) and “L1_lowFeatures_transformed_randomr” (r=1:10).

func_log_regression_individual: R pipeline of individual functional logistic regressions on functional predictors, using the transformed data “L1_transformed_random_r” (r=1:10) as input. The package ‘fda.usc’ is used here.

glm_scalar_no_lowFeatures: R pipeline of individual logistic regressions on scalar predictors for each of the 6 comparisons, using the transformed data “L1_transformed_random_r” (r=1:10) as input. At this stage, the 7 low-resolution features are not included (since they were transformed separately). Here ‘glm’ is used, since in the case of non-functional variables, fda.usc uses the same generalized linear models as glm functions.

glm_lowFeatures.R: R pipeline of individual logistic regressions on the 7 low-resolution features, using the transformed data “L1_lowFeatures_transformed_randomr” (r=1:10) as the input. Generally it is similar with the above pipeline for scalar predictors.


In each folder (taking random1 as an example), there are 4 RData and 18 tables, they are:
1.’L1_autosomes_results_smoothed_mean_1’: random sample generated based on de novo L1s, pre-transformation variables can be found in “result_mean” after loading the data.

2.’L1_transformed_random_1.RData’: transformed variables for random sample 1, post-transformation variables can be found in “result_mean” after loading the data.
(the raw files, plots and transformation pipeline can be find at folder :’L1/transformed_data’ on google drive).

3.’functional_individual_regression_random_1.RData’: data image of individual functional logistic regression for 6 comparisons. The summary of all individual regressions can be found in ‘comp1_func_logit_sum’ (summary of peudoR2 and smallest pvals of 20 coefficients for regressions on all functional variables, for comparison1 in this case) and ‘comp1_func_logit_20coeff’ (all 20 coefficients and their values for regressions on all functional variables,for comparison1 in this case), after loading the data.

4. ’L1_lowFeatures_transformed_random1.RData’: transformed variables for 7 low-resolution features random sample 1, post-transformation variables can be found in “regionsFeatures” after loading the data.

5.’scalar_noLowres_individual_regression_random1.RData’: data image of individual logistic regression on scalar variables (excluding 7 low-resolution features) for 6 comparisons. The summary of all individual regressions can be found in ‘comp1_scalar_logit_sum’ (summary of coefficient, p-values,peudoR2 for regressions on scalar variables, for comparison1 in this case)

6. ‘low_features_individual_regression_random1.RData’: data image of individual logistic regression on 7 low-resolution features for 6 comparisons. Similar pipeline and scheme as in “5” (the above one for scalar predictors). 

7. 24 tables are the written outputs from 3,4 and 6, with each comparison written separately : 

comp1_func_logit_sum.txt
comp2_func_logit_sum.txt
comp3_func_logit_sum.txt
comp4_func_logit_sum.txt
comp5_func_logit_sum.txt
comp6_func_logit_sum.txt

comp1_func_logit_coefficients.txt
comp2_func_logit_coefficients.txt
comp3_func_logit_coefficients.txt
comp4_func_logit_coefficients.txt
comp5_func_logit_coefficients.txt
comp6_func_logit_coefficients.txt

comp1_scalar_logit_sum.txt
comp2_scalar_logit_sum.txt
comp3_scalar_logit_sum.txt
comp4_scalar_logit_sum.txt
comp5_scalar_logit_sum.txt
comp6_scalar_logit_sum.txt

comp1_lowfeature_logit_sum.txt
comp2_lowfeature_logit_sum.txt
comp3_lowfeature_logit_sum.txt
comp4_lowfeature_logit_sum.txt
comp5_lowfeature_logit_sum.txt
comp6_lowfeature_logit_sum.txt

p.s. Further analysis should be focused on random sample 1 (files in the folder “random1”), since we selected the sample to be a representative one based on the individual regression results.
