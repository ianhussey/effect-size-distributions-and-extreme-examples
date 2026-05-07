################################################################
### COVIDSTRESS database
### 
### Blackburn et al. (2022) https://www.nature.com/articles/s41597-022-01383-6
### data available at https://osf.io/36tsd/
###
### Perceived Stress Scale
###
### Full dataset to determine 2-factor solution
###

Final_COVIDiSTRESS_Vol2_cleaned <- read.csv("Final_COVIDiSTRESS_Vol2_cleaned.csv")
colnames(Final_COVIDiSTRESS_Vol2_cleaned)
summary(Final_COVIDiSTRESS_Vol2_cleaned)
table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)

mydata <- Final_COVIDiSTRESS_Vol2_cleaned[,45:54]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
fit1 <- fa(mydata,2)
fit1
diagram(fit1)

EGAmodel= '
 factor1 =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6
 factor2 =~ perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'




################################################################
### COVIDSTRESS database
### 
### UserLanguage = AR (Arabic)
### 
###

mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="AR")
mydata <- mydata[,45:54]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .87, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 component
# Eigenvalue 1 = 4.11; eigenvalue 2 = .60

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 component
# Eigenvalue 1 = 4.61; eigenvalue 2 = .68

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.41, RMSEA=.12, RMSR=.09, TLI=.818

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.164, RMSR=.1, TLI=.742

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.08, RMSR=.05, TLI=.917
#      MR1   MR2
#MR1  1.00 -0.53
#MR2 -0.53  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.128, RMSR=.05, TLI=.842
#      MR1   MR2
#MR1  1.00 -0.53
#MR2 -0.53  1.00

# Single factor model lavaan
UNImodel= '
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.965       0.892
#Tucker-Lewis Index (TLI)                       0.954       0.856
#Robust Comparative Fit Index (CFI)                         0.799
#Robust Tucker-Lewis Index (TLI)                            0.732
#RMSEA                                          0.150       0.184
#Robust RMSEA                                               0.183
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .504

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.842       0.858
#Tucker-Lewis Index (TLI)                       0.789       0.811
#Robust Comparative Fit Index (CFI)                         0.857
#Robust Tucker-Lewis Index (TLI)                            0.810
#RMSEA                                          0.139       0.116
#Robust RMSEA                                               0.131
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .430

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.956
#Tucker-Lewis Index (TLI)                       0.987       0.939
#Robust Comparative Fit Index (CFI)                         0.886
#Robust Tucker-Lewis Index (TLI)                            0.842
#RMSEA                                          0.078       0.119
#Robust RMSEA                                               0.140
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .579

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.705    0.057  -12.420    0.000   -0.705   -0.705

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.919       0.927
#Tucker-Lewis Index (TLI)                       0.888       0.899
#Robust Comparative Fit Index (CFI)                         0.931
#Robust Tucker-Lewis Index (TLI)                            0.904
#RMSEA                                          0.102       0.085
#Robust RMSEA                                               0.093
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .516

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2         -0.713    0.099   -7.216    0.000   -0.713   -0.713

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.808       0.801
#Tucker-Lewis Index (TLI)                       0.744       0.735
#Robust Comparative Fit Index (CFI)                         0.816
#Robust Tucker-Lewis Index (TLI)                            0.754
#RMSEA                                          0.154       0.137
#Robust RMSEA                                               0.149
#SRMR                                           0.239       0.239

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .515


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6
 factor2 =~ perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
# warning
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       1.000
#Tucker-Lewis Index (TLI)                       1.012       1.030
#Robust Comparative Fit Index (CFI)                         1.000
#Robust Tucker-Lewis Index (TLI)                            1.027
#RMSEA                                          0.000       0.000
#Robust RMSEA                                               0.000
#SRMR                                           0.030       0.030

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .465

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    prcvd_strss__3   -0.044    0.096   -0.458    0.647   -0.044   -0.045
#    prcvd_strss__2   -0.736    0.472   -1.560    0.119   -0.736   -0.662
#    prcvd_strs__10    0.361    0.217    1.662    0.097    0.361    0.357
#    prcvd_strss__9    0.219    0.141    1.556    0.120    0.219    0.223
#    prcvd_strss__6   -0.014    0.126   -0.111    0.911   -0.014   -0.015
#  factor2 =~                                                            
#    prcvd_strss__4    0.407    0.089    4.552    0.000    0.407    0.444
#    prcvd_strss__8    0.531    0.110    4.812    0.000    0.531    0.526
#    prcvd_strss__5    0.234    0.108    2.171    0.030    0.234    0.262
#    prcvd_strss__7    0.755    0.137    5.513    0.000    0.755    0.764
#  general =~                                                            
#    prcvd_strss__3    0.671    0.083    8.088    0.000    0.671    0.693
#    prcvd_strss__2    0.935    0.102    9.128    0.000    0.935    0.842
#    prcvd_strs__10    0.817    0.096    8.548    0.000    0.817    0.806
#    prcvd_strss__9    0.611    0.092    6.623    0.000    0.611    0.622
#    prcvd_strss__6    0.656    0.078    8.377    0.000    0.656    0.697
#    prcvd_strss__4   -0.445    0.080   -5.597    0.000   -0.445   -0.486
#    prcvd_strss__8   -0.416    0.094   -4.447    0.000   -0.416   -0.412
#    prcvd_strss__5   -0.529    0.081   -6.534    0.000   -0.529   -0.591
#    prcvd_strss__7   -0.386    0.102   -3.788    0.000   -0.386   -0.391

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6748963      0.5555556      0.6634306      0.2929382 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1855028 0.1151313 0.8144972 0.8891858 0.00135989 0.4952088 1.0388113
#factor2 0.5534978 0.2099725 0.4465022 0.7927621 0.41997096 0.6776299 0.8654432
#general 0.6748963 0.6748963 0.6748963 0.6634306 0.29293817 0.8890623 0.9782090





################################################################
### COVIDSTRESS database
### 
### UserLanguage = BG (Bulgarian)
### 
###

mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="BG")
mydata <- mydata[,45:54]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .88, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 component
# Eigenvalue 1 = 4.39; eigenvalue 2 = .73

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 component
# Eigenvalue 1 = 4.89; eigenvalue 2 = .79

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.44, RMSEA=.133, RMSR=.09, TLI=.816

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.161, RMSR=.1, TLI=.787

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.062, RMSR=.03, TLI=.960
#      MR1   MR2
#MR1  1.00 -0.57
#MR2 -0.57  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.088, RMSR=.03, TLI=.936
#      MR1   MR2
#MR1  1.00 -0.59
#MR2 -0.59  1.00

# Single factor model lavaan
UNImodel= '
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.969       0.875
#Tucker-Lewis Index (TLI)                       0.958       0.833
#Robust Comparative Fit Index (CFI)                         0.827
#Robust Tucker-Lewis Index (TLI)                            0.770
#RMSEA                                          0.155       0.204
#Robust RMSEA                                               0.172
#SRMR                                           0.091       0.091

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .537

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.853       0.831
#Tucker-Lewis Index (TLI)                       0.804       0.775
#Robust Comparative Fit Index (CFI)                         0.856
#Robust Tucker-Lewis Index (TLI)                            0.807
#RMSEA                                          0.141       0.129
#Robust RMSEA                                               0.139
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.993       0.956
#Tucker-Lewis Index (TLI)                       0.991       0.939
#Robust Comparative Fit Index (CFI)                         0.942
#Robust Tucker-Lewis Index (TLI)                            0.920
#RMSEA                                          0.072       0.123
#Robust RMSEA                                               0.102
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .655

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.717    0.034  -21.105    0.000   -0.717   -0.717

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.963       0.962
#Tucker-Lewis Index (TLI)                       0.949       0.947
#Robust Comparative Fit Index (CFI)                         0.967
#Robust Tucker-Lewis Index (TLI)                            0.954
#RMSEA                                          0.072       0.062
#Robust RMSEA                                               0.068
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .578

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.718    0.059  -12.148    0.000   -0.718   -0.718

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.843       0.826
#Tucker-Lewis Index (TLI)                       0.790       0.769
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.795
#RMSEA                                          0.145       0.131
#Robust RMSEA                                               0.143
#SRMR                                           0.240       0.240

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .548


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6
 factor2 =~ perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.984       0.978
#Tucker-Lewis Index (TLI)                       0.969       0.955
#Robust Comparative Fit Index (CFI)                         0.984
#Robust Tucker-Lewis Index (TLI)                            0.968
#RMSEA                                          0.056       0.057
#Robust RMSEA                                               0.057
#SRMR                                           0.031       0.031

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .565

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5953601      0.5555556      0.7246279      0.1671787 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3636478 0.2347269 0.6363522 0.9214668 0.1541120 1.3462135 1.3218637
#factor2 0.4792743 0.1699130 0.5207257 0.8094608 0.3930995 0.5790785 0.7820612
#general 0.5953601 0.5953601 0.5953601 0.7246279 0.1671787 0.8804066 0.9285248







################################################################
### COVIDSTRESS database
### 
### UserLanguage = DE (German)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="DE")
mydata <- mydata[,45:54]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .9, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 4.82; eigenvalue 2 = .39

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 5.44; eigenvalue 2 = .41

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.48, RMSEA=.104, RMSR=.06, TLI=.896

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.129, RMSR=.06, TLI=.875

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 1 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.068, RMSR=.03, TLI=.955
#      MR1   MR2
#MR1  1.00 -0.72
#MR2 -0.72  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.60, RMSEA=.091, RMSR=.03, TLI=.937
#      MR1   MR2
#MR1  1.00 -0.74
#MR2 -0.74  1.00

# Single factor model lavaan
UNImodel= '
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.994       0.971
#Tucker-Lewis Index (TLI)                       0.993       0.961
#Robust Comparative Fit Index (CFI)                         0.925
#Robust Tucker-Lewis Index (TLI)                            0.899
#RMSEA                                          0.082       0.127
#Robust RMSEA                                               0.124
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .533

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.939       0.935
#Tucker-Lewis Index (TLI)                       0.919       0.913
#Robust Comparative Fit Index (CFI)                         0.943
#Robust Tucker-Lewis Index (TLI)                            0.924
#RMSEA                                          0.097       0.079
#Robust RMSEA                                               0.094
#SRMR                                           0.044       0.044

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .490

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.998       0.989
#Tucker-Lewis Index (TLI)                       0.998       0.984
#Robust Comparative Fit Index (CFI)                         0.966
#Robust Tucker-Lewis Index (TLI)                            0.953
#RMSEA                                          0.045       0.081
#Robust RMSEA                                               0.085
#SRMR                                           0.033       0.033

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .577

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.878    0.012  -71.038    0.000   -0.878   -0.878

# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.975       0.974
#Tucker-Lewis Index (TLI)                       0.965       0.964
#Robust Comparative Fit Index (CFI)                         0.978
#Robust Tucker-Lewis Index (TLI)                            0.969
#RMSEA                                          0.064       0.051
#Robust RMSEA                                               0.060
#SRMR                                           0.031       0.031

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .521

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.868    0.019  -45.327    0.000   -0.868   -0.868

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.789       0.757
#Tucker-Lewis Index (TLI)                       0.719       0.676
#Robust Comparative Fit Index (CFI)                         0.791
#Robust Tucker-Lewis Index (TLI)                            0.722
#RMSEA                                          0.182       0.152
#Robust RMSEA                                               0.180
#SRMR                                           0.314       0.314

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .528


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6
 factor2 =~ perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.974
#Tucker-Lewis Index (TLI)                       0.954       0.948
#Robust Comparative Fit Index (CFI)                         0.979
#Robust Tucker-Lewis Index (TLI)                            0.957
#RMSEA                                          0.074       0.061
#Robust RMSEA                                               0.070
#SRMR                                           0.028       0.028

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .546

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#  0.8577013076   0.5555555556   0.4319990488   0.0007272546 
#
#$FactorLevelIndices
#            ECV_SS     ECV_SG    ECV_GS     Omega       OmegaH          H        FD
#factor1 0.25383361 0.12708677 0.7461664 0.8268341 0.2075853191 0.42277405 0.6858517
#factor2 0.03046464 0.01521192 0.9695354 0.8612143 0.0050274762 0.07182248 0.3920095
#general 0.85770131 0.85770131 0.8577013 0.4319990 0.0007272546 0.90294090 0.9432194







################################################################
### COVIDSTRESS database
### 
### UserLanguage = JA (Japan)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="JA")
mydata <- mydata[,45:54]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .84, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 component
# Eigenvalue 1 = 3.65; eigenvalue 2 = .99

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 component
# Eigenvalue 1 = 4.05; eigenvalue 2 = 1.1

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.36, RMSEA=.158, RMSR=.12, TLI=.700

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.184, RMSR=.13, TLI=.675

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.075, RMSR=.03, TLI=.932
#      MR1   MR2
#MR1  1.00 -0.43
#MR2 -0.43  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.094, RMSR=.03, TLI=.915
#      MR1   MR2
#MR1  1.00 -0.43
#MR2 -0.43  1.00

# Single factor model lavaan
UNImodel= '
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.914       0.796
#Tucker-Lewis Index (TLI)                       0.885       0.728
#Robust Comparative Fit Index (CFI)                         0.726
#Robust Tucker-Lewis Index (TLI)                            0.634
#RMSEA                                          0.221       0.236
#Robust RMSEA                                               0.200
#SRMR                                           0.127       0.127

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .422

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.764       0.738
#Tucker-Lewis Index (TLI)                       0.685       0.650
#Robust Comparative Fit Index (CFI)                         0.765
#Robust Tucker-Lewis Index (TLI)                            0.686
#RMSEA                                          0.166       0.147
#Robust RMSEA                                               0.166
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .419

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.994       0.976
#Tucker-Lewis Index (TLI)                       0.992       0.966
#Robust Comparative Fit Index (CFI)                         0.958
#Robust Tucker-Lewis Index (TLI)                            0.942
#RMSEA                                          0.058       0.083
#Robust RMSEA                                               0.080
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .578

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.468    0.020  -23.425    0.000   -0.468   -0.468

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.970       0.967
#Tucker-Lewis Index (TLI)                       0.958       0.955
#Robust Comparative Fit Index (CFI)                         0.971
#Robust Tucker-Lewis Index (TLI)                            0.960
#RMSEA                                          0.061       0.053
#Robust RMSEA                                               0.060
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .520

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.473    0.033  -14.265    0.000   -0.473   -0.473

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.919       0.915
#Tucker-Lewis Index (TLI)                       0.893       0.886
#Robust Comparative Fit Index (CFI)                         0.921
#Robust Tucker-Lewis Index (TLI)                            0.894
#RMSEA                                          0.097       0.084
#Robust RMSEA                                               0.096
#SRMR                                           0.153       0.153

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .524


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6
 factor2 =~ perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
 general =~ perceived_stress_sca_3+perceived_stress_sca_2+perceived_stress_sca_2+
            perceived_stress_sca_10+perceived_stress_sca_9+perceived_stress_sca_6+
            perceived_stress_sca_4+perceived_stress_sca_8+perceived_stress_sca_5+
            perceived_stress_sca_7
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.990       0.986
#Tucker-Lewis Index (TLI)                       0.980       0.972
#Robust Comparative Fit Index (CFI)                         0.990
#Robust Tucker-Lewis Index (TLI)                            0.980
#RMSEA                                          0.042       0.042
#Robust RMSEA                                               0.042
#SRMR                                           0.018       0.018

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .519

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6575923      0.5555556      0.7360512      0.3873140 
#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.07578339 0.04730425 0.9242166 0.8713193 0.00493139 0.2030768 0.6776725
#factor2 0.78527427 0.29510341 0.2147257 0.7589005 0.60194257 0.6958192 0.8487897
#general 0.65759234 0.65759234 0.6575923 0.7360512 0.38731395 0.8824718 0.9470348






