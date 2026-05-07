################################################################
### COVIDSTRESS database
### 
### Blackburn et al. (2022) https://www.nature.com/articles/s41597-022-01383-6
### data available at https://osf.io/36tsd/
###
### Trust Scale (reappraisal and suppression)
###
### Full dataset to determine 2-factor solution
###

Final_COVIDiSTRESS_Vol2_cleaned <- read.csv("Final_COVIDiSTRESS_Vol2_cleaned.csv")
colnames(Final_COVIDiSTRESS_Vol2_cleaned)
table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)

mydata <- Final_COVIDiSTRESS_Vol2_cleaned[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
fit1 <- fa(mydata,2)
fit1
diagram(fit1)

EGAmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
'


################################################################
### COVIDSTRESS database
### 
### UserLanguage = DA (Dari)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="DA")
mydata <- mydata[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
omega(mydata) # alpha = .86, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.38; eigenvalue 2 = .59

rho <- polychoric(mydata)$rho
# No polychoric needed for 11 categories
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 1 components
# Eigenvalue 1 = 3.55; eigenvalue 2 = .49

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.48, RMSEA=.22, RMSR=.11, TLI=.679

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.23, RMSR=.1, TLI=.685

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.156, RMSR=.05, TLI=.84
#      MR1   MR2
#MR1  1.00  0.54
#MR2  0.54  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.207, RMSR=.06, TLI=.744
#      MR1   MR2
#MR1  1.00  0.62
#MR2  0.62  1.00


# Single factor model lavaan
UNImodel= '
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.978       0.932
#Tucker-Lewis Index (TLI)                       0.967       0.898
#Robust Comparative Fit Index (CFI)                         0.798
#Robust Tucker-Lewis Index (TLI)                            0.697
#RMSEA                                          0.173       0.230
#Robust RMSEA                                               0.236
#SRMR                                           0.091       0.091

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .607

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.799       0.738
#Tucker-Lewis Index (TLI)                       0.698       0.606
#Robust Comparative Fit Index (CFI)                         0.797
#Robust Tucker-Lewis Index (TLI)                            0.696
#RMSEA                                          0.220       0.209
#Robust RMSEA                                               0.219
#SRMR                                           0.106       0.106

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .488

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.935
#Tucker-Lewis Index (TLI)                       0.966       0.894
#Robust Comparative Fit Index (CFI)                         0.808
#Robust Tucker-Lewis Index (TLI)                            0.690
#RMSEA                                          0.176       0.234
#Robust RMSEA                                               0.239
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .610

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.843    0.055   15.361    0.000    0.843    0.843

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.828       0.785
#Tucker-Lewis Index (TLI)                       0.722       0.653
#Robust Comparative Fit Index (CFI)                         0.829
#Robust Tucker-Lewis Index (TLI)                            0.723
#RMSEA                                          0.211       0.197
#Robust RMSEA                                               0.209
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .595

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.648    0.186    3.487    0.000    0.648    0.648

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.745       0.742
#Tucker-Lewis Index (TLI)                       0.618       0.613
#Robust Comparative Fit Index (CFI)                         0.753
#Robust Tucker-Lewis Index (TLI)                            0.629
#RMSEA                                          0.248       0.208
#Robust RMSEA                                               0.242
#SRMR                                           0.244       0.244

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .616


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.945       0.928
#Tucker-Lewis Index (TLI)                       0.835       0.785
#Robust Comparative Fit Index (CFI)                         0.945
#Robust Tucker-Lewis Index (TLI)                            0.835
#RMSEA                                          0.163       0.155
#Robust RMSEA                                               0.161
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .625

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    trust_1           1.637    0.201    8.167    0.000    1.637    0.759
#    trust_2           1.698    0.220    7.727    0.000    1.698    0.823
#    trust_3           1.632    0.205    7.956    0.000    1.632    0.782
#    trust_6           0.946    0.375    2.521    0.012    0.946    0.446
#    trust_4           0.734    0.337    2.180    0.029    0.734    0.496
#  factor2 =~                                                            
#    trust_7          -1.182    3.202   -0.369    0.712   -1.182   -1.115
#    trust_5           1.425    1.395    1.022    0.307    1.425    0.622
#  general =~                                                            
#    trust_1           0.480    0.356    1.348    0.178    0.480    0.223
#    trust_2           0.575    0.359    1.603    0.109    0.575    0.279
#    trust_3           0.724    0.351    2.062    0.039    0.724    0.347
#    trust_4           0.737    0.502    1.468    0.142    0.737    0.499
#    trust_5           2.538    1.393    1.823    0.068    2.538    1.109
#    trust_6           0.932    0.595    1.566    0.117    0.932    0.440
#    trust_7           1.099    0.528    2.083    0.037    1.099    1.037

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.4316273      0.4761905      0.9975820      0.5789270 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H       FD
#factor1 0.7701676 0.3332695 0.2298324 0.8760743 0.67814079 0.8487979 1.000214
#factor2 0.4144418 0.2351032 0.5855582 1.6645211 0.08340889 1.2874946 1.701702
#general 0.4316273 0.4316273 0.4316273 0.9975820 0.57892700 1.0557966 1.347378





################################################################
### COVIDSTRESS database
### 
### UserLanguage = ES-HN (Spanish-Honduras)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="ES-HN")
mydata <- mydata[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
omega(mydata) # alpha = .86, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 3.6; eigenvalue 2 = .47

rho <- polychoric(mydata)$rho
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 4.02; eigenvalue 2 = .45

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.51, RMSEA=.185, RMSR=.09, TLI=.811

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.199, RMSR=.09, TLI=.83

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.62, RMSEA=.157, RMSR=.04, TLI=.864
#      MR1   MR2
#MR1  1.00  0.47
#MR2  0.47  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.177, RMSR=.03, TLI=.865
#      MR1   MR2
#MR1  1.00  0.49
#MR2  0.49  1.00


# Single factor model lavaan
UNImodel= '
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.985       0.946
#Tucker-Lewis Index (TLI)                       0.977       0.919
#Robust Comparative Fit Index (CFI)                         0.887
#Robust Tucker-Lewis Index (TLI)                            0.831
#RMSEA                                          0.171       0.235
#Robust RMSEA                                               0.200
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .727

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.876       0.865
#Tucker-Lewis Index (TLI)                       0.814       0.797
#Robust Comparative Fit Index (CFI)                         0.877
#Robust Tucker-Lewis Index (TLI)                            0.815
#RMSEA                                          0.185       0.170
#Robust RMSEA                                               0.184
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .601

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.982
#Tucker-Lewis Index (TLI)                       0.994       0.971
#Robust Comparative Fit Index (CFI)                         0.933
#Robust Tucker-Lewis Index (TLI)                            0.892
#RMSEA                                          0.090       0.141
#Robust RMSEA                                               0.160
#SRMR                                           0.045       0.045

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .736

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.502    0.048   10.534    0.000    0.502    0.502

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.923       0.912
#Tucker-Lewis Index (TLI)                       0.876       0.859
#Robust Comparative Fit Index (CFI)                         0.924
#Robust Tucker-Lewis Index (TLI)                            0.877
#RMSEA                                          0.151       0.142
#Robust RMSEA                                               0.150
#SRMR                                           0.051       0.051

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .627

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.478    0.056    8.543    0.000    0.478    0.478

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.872       0.856
#Tucker-Lewis Index (TLI)                       0.808       0.785
#Robust Comparative Fit Index (CFI)                         0.873
#Robust Tucker-Lewis Index (TLI)                            0.809
#RMSEA                                          0.188       0.176
#Robust RMSEA                                               0.187
#SRMR                                           0.185       0.185

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .615


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.992       0.987
#Tucker-Lewis Index (TLI)                       0.976       0.960
#Robust Comparative Fit Index (CFI)                         0.991
#Robust Tucker-Lewis Index (TLI)                            0.974
#RMSEA                                          0.067       0.075
#Robust RMSEA                                               0.070
#SRMR                                           0.019       0.019

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .720

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    trust_1           1.463    0.154    9.487    0.000    1.463    0.679
#    trust_2           0.937    0.308    3.042    0.002    0.937    0.410
#    trust_3           0.462    0.327    1.412    0.158    0.462    0.199
#    trust_6           1.179    0.182    6.478    0.000    1.179    0.554
#    trust_4           0.047    0.188    0.249    0.804    0.047    0.020
#  factor2 =~                                                            
#    trust_7           2.005    0.271    7.391    0.000    2.005    0.674
#    trust_5           1.376    0.038   35.951    0.000    1.376    0.483
#  general =~                                                            
#    trust_1           1.224    0.202    6.071    0.000    1.224    0.568
#    trust_2           1.735    0.198    8.773    0.000    1.735    0.760
#    trust_3           2.018    0.154   13.132    0.000    2.018    0.870
#    trust_4           1.995    0.132   15.129    0.000    1.995    0.848
#    trust_5           1.419    0.162    8.732    0.000    1.419    0.498
#    trust_6           1.323    0.184    7.192    0.000    1.323    0.622
#    trust_7           0.847    0.156    5.437    0.000    0.847    0.285

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6500861      0.4761905      0.9165285      0.7375500 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2611093 0.2052727 0.7388907 0.9307210 0.1907460 0.6067037 0.8150270
#factor2 0.6763864 0.1446412 0.3236136 0.6652364 0.4563582 0.5324382 0.7557956
#general 0.6500861 0.6500861 0.6500861 0.9165285 0.7375500 0.8954290 0.9313510





################################################################
### COVIDSTRESS database
### 
### UserLanguage = SV (Swedish)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="SV")
mydata <- mydata[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
omega(mydata) # alpha = .89, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 components
# Eigenvalue 1 = 3.8; eigenvalue 2 = .3

rho <- polychoric(mydata)$rho
# No polychoric needed for 11 categories
library(EGAnet)
rho <- polychoric.matrix((mydata))
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 factors and 1 components
# Eigenvalue 1 = 3.82; eigenvalue 2 = .27

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.54, RMSEA=.183, RMSR=.07, TLI=.825

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.147, RMSR=.06, TLI=.883

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 1 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.201, RMSR=.04, TLI=.788
#      MR1   MR2
#MR1  1.00  0.71
#MR2  0.71  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.65, RMSEA=.126, RMSR=.03, TLI=.913
#      MR1   MR2
#MR1  1.00  0.56
#MR2  0.56  1.00

# Single factor model lavaan
UNImodel= '
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.981
#Tucker-Lewis Index (TLI)                       0.994       0.971
#Robust Comparative Fit Index (CFI)                         0.931
#Robust Tucker-Lewis Index (TLI)                            0.897
#RMSEA                                          0.083       0.142
#Robust RMSEA                                               0.142
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .464

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.890       0.877
#Tucker-Lewis Index (TLI)                       0.836       0.816
#Robust Comparative Fit Index (CFI)                         0.894
#Robust Tucker-Lewis Index (TLI)                            0.841
#RMSEA                                          0.180       0.159
#Robust RMSEA                                               0.176
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .492

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.992
#Tucker-Lewis Index (TLI)                       1.000       0.988
#Robust Comparative Fit Index (CFI)                         0.966
#Robust Tucker-Lewis Index (TLI)                            0.944
#RMSEA                                          0.023       0.093
#Robust RMSEA                                               0.104
#SRMR                                           0.037       0.037

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .567

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.759    0.053   14.319    0.000    0.759    0.759

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.915       0.905
#Tucker-Lewis Index (TLI)                       0.862       0.847
#Robust Comparative Fit Index (CFI)                         0.918
#Robust Tucker-Lewis Index (TLI)                            0.868
#RMSEA                                          0.165       0.145
#Robust RMSEA                                               0.160
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .531

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.735    0.080    9.150    0.000    0.735    0.735

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.816       0.804
#Tucker-Lewis Index (TLI)                       0.724       0.706
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.730
#RMSEA                                          0.233       0.201
#Robust RMSEA                                               0.229
#SRMR                                           0.264       0.264

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .665


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.984       0.980
#Tucker-Lewis Index (TLI)                       0.952       0.939
#Robust Comparative Fit Index (CFI)                         0.985
#Robust Tucker-Lewis Index (TLI)                            0.954
#RMSEA                                          0.097       0.091
#Robust RMSEA                                               0.095
#SRMR                                           0.035       0.035

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .649

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    trust_1           2.264    2.519    0.899    0.369    2.264    0.818
#    trust_2          -0.145    0.226   -0.644    0.519   -0.145   -0.060
#    trust_3          -0.003    0.249   -0.014    0.989   -0.003   -0.001
#    trust_6           0.444    0.634    0.700    0.484    0.444    0.155
#    trust_4          -0.398    0.449   -0.885    0.376   -0.398   -0.199
#  factor2 =~                                                            
#    trust_7           0.602    0.125    4.813    0.000    0.602    0.377
#    trust_5           1.428    0.183    7.805    0.000    1.428    0.619
#  general =~                                                            
#    trust_1           2.472    0.286    8.637    0.000    2.472    0.893
#    trust_2           1.790    0.207    8.639    0.000    1.790    0.741
#    trust_3           2.144    0.207   10.333    0.000    2.144    0.915
#    trust_4           1.348    0.260    5.183    0.000    1.348    0.674
#    trust_5           1.215    0.230    5.286    0.000    1.215    0.526
#    trust_6           2.200    0.222    9.890    0.000    2.200    0.770
#    trust_7           0.848    0.215    3.947    0.000    0.848    0.531

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.7503346      0.4761905      0.9327209      0.8809640 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1856613 0.1458502 0.8143387 0.9410213 0.02906599 0.6768755 1.6440395
#factor2 0.4841476 0.1038152 0.5158524 0.6969148 0.32743674 0.4397630 0.7492231
#general 0.7503346 0.7503346 0.7503346 0.9327209 0.88096401 0.9303128 0.9718263





################################################################
### COVIDSTRESS database
### 
### UserLanguage = TR (Turkish)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="TR")
mydata <- mydata[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
omega(mydata) # alpha = .83, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 3.35; eigenvalue 2 = .5

rho <- polychoric(mydata)$rho
# No polychoric needed for 11 categories
library(EGAnet)
rho <- polychoric.matrix((mydata))
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 3.63; eigenvalue 2 = .57

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.48, RMSEA=.163, RMSR=.1, TLI=.833

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.199, RMSR=.11, TLI=.803

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.6, RMSEA=.117, RMSR=.03, TLI=.913
#      MR1   MR2
#MR1  1.00  0.26
#MR2  0.26  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.66, RMSEA=.148, RMSR=.03, TLI=.89
#      MR1   MR2
#MR1  1.00  0.20
#MR2  0.20  1.00

# Single factor model lavaan
UNImodel= '
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.982       0.932
#Tucker-Lewis Index (TLI)                       0.973       0.898
#Robust Comparative Fit Index (CFI)                         0.874
#Robust Tucker-Lewis Index (TLI)                            0.811
#RMSEA                                          0.186       0.269
#Robust RMSEA                                               0.199
#SRMR                                           0.101       0.101

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .714

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.894
#Tucker-Lewis Index (TLI)                       0.841       0.841
#Robust Comparative Fit Index (CFI)                         0.897
#Robust Tucker-Lewis Index (TLI)                            0.846
#RMSEA                                          0.162       0.147
#Robust RMSEA                                               0.159
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .67

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.979
#Tucker-Lewis Index (TLI)                       0.994       0.966
#Robust Comparative Fit Index (CFI)                         0.930
#Robust Tucker-Lewis Index (TLI)                            0.888
#RMSEA                                          0.085       0.155
#Robust RMSEA                                               0.154
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .755

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.275    0.081    3.418    0.001    0.275    0.275

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.942       0.940
#Tucker-Lewis Index (TLI)                       0.906       0.902
#Robust Comparative Fit Index (CFI)                         0.944
#Robust Tucker-Lewis Index (TLI)                            0.909
#RMSEA                                          0.124       0.115
#Robust RMSEA                                               0.122
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .669

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.292    0.100    2.912    0.004    0.292    0.292

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.930       0.933
#Tucker-Lewis Index (TLI)                       0.894       0.899
#Robust Comparative Fit Index (CFI)                         0.933
#Robust Tucker-Lewis Index (TLI)                            0.900
#RMSEA                                          0.132       0.117
#Robust RMSEA                                               0.128
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .664


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.988       0.989
#Tucker-Lewis Index (TLI)                       0.965       0.967
#Robust Comparative Fit Index (CFI)                         0.990
#Robust Tucker-Lewis Index (TLI)                            0.969
#RMSEA                                          0.076       0.067
#Robust RMSEA                                               0.071
#SRMR                                           0.028       0.028

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .666

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    trust_1           2.218    0.544    4.078    0.000    2.218    0.792
#    trust_2           2.183    0.515    4.239    0.000    2.183    0.707
#    trust_3           1.540    0.451    3.414    0.001    1.540    0.578
#    trust_6           1.830    0.615    2.974    0.003    1.830    0.580
#    trust_4           0.131    1.415    0.093    0.926    0.131    0.047
#  factor2 =~                                                            
#    trust_7           1.527    0.377    4.053    0.000    1.527    0.695
#    trust_5           1.288    0.114   11.294    0.000    1.288    0.503
#  general =~                                                            
#    trust_1           1.257    0.840    1.495    0.135    1.257    0.449
#    trust_2           1.474    0.786    1.876    0.061    1.474    0.478
#    trust_3           1.506    0.457    3.295    0.001    1.506    0.565
#    trust_4           2.797    0.849    3.295    0.001    2.797    1.006
#    trust_5           0.601    0.418    1.438    0.150    0.601    0.235
#    trust_6           1.847    0.585    3.158    0.002    1.847    0.586
#    trust_7           0.641    0.306    2.094    0.036    0.641    0.291

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.4695812      0.4761905      0.9074838      0.5429529 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4609867 0.3765616 0.5390133 0.9388596 0.4081773 0.7869179 0.9388882
#factor2 0.8401072 0.1538572 0.1598928 0.6033918 0.5058026 0.5596120 0.7732652
#general 0.4695812 0.4695812 0.4695812 0.9074838 0.5429529 1.0132690 1.0080187





################################################################
### COVIDSTRESS database
### 
### UserLanguage = UK (English)
### 
###
### not included because too many bad analyses

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="UK")
mydata <- mydata[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
omega(mydata) # alpha = .86, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.43; eigenvalue 2 = .41

rho <- polychoric(mydata)$rho
# No polychoric needed for 11 categories
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.55; eigenvalue 2 = .42

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.49, RMSEA=.148, RMSR=.09, TLI=.851

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.159, RMSR=.09, TLI=.845

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.086, RMSR=.02, TLI=.95
#      MR1   MR2
#MR1  1.00  0.49
#MR2  0.49  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.11, RMSR=.03, TLI=.926
#      MR1   MR2
#MR1  1.00  0.49
#MR2  0.49  1.00

# Single factor model lavaan
UNImodel= '
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
#warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.985       0.946
#Tucker-Lewis Index (TLI)                       0.977       0.918
#Robust Comparative Fit Index (CFI)                         0.898
#Robust Tucker-Lewis Index (TLI)                            0.847
#RMSEA                                          0.134       0.192
#Robust RMSEA                                               0.160
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .618

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.903       0.907
#Tucker-Lewis Index (TLI)                       0.855       0.861
#Robust Comparative Fit Index (CFI)                         0.906
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.148       0.136
#Robust RMSEA                                               0.146
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .602

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.991
#Tucker-Lewis Index (TLI)                       0.999       0.985
#Robust Comparative Fit Index (CFI)                         0.967
#Robust Tucker-Lewis Index (TLI)                            0.946
#RMSEA                                          0.027       0.081
#Robust RMSEA                                               0.095
#SRMR                                           0.035       0.035

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .647

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.556    0.057    9.732    0.000    0.556    0.556

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.971       0.973
#Tucker-Lewis Index (TLI)                       0.953       0.957
#Robust Comparative Fit Index (CFI)                         0.973
#Robust Tucker-Lewis Index (TLI)                            0.957
#RMSEA                                          0.085       0.076
#Robust RMSEA                                               0.081
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .612

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.558    0.060    9.321    0.000    0.558    0.558

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.900       0.906
#Tucker-Lewis Index (TLI)                       0.850       0.860
#Robust Comparative Fit Index (CFI)                         0.903
#Robust Tucker-Lewis Index (TLI)                            0.854
#RMSEA                                          0.151       0.136
#Robust RMSEA                                               0.148
#SRMR                                           0.199       0.199

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .617

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.900       0.906
#Tucker-Lewis Index (TLI)                       0.850       0.860
#Robust Comparative Fit Index (CFI)                         0.903
#Robust Tucker-Lewis Index (TLI)                            0.854
#RMSEA                                          0.151       0.136
#Robust RMSEA                                               0.148
#SRMR                                           0.199       0.199

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .617


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
# warning
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning with bad weight estimates; analysis stopped here





################################################################
### COVIDSTRESS database
### 
### UserLanguage = ZH-T (Taiwan)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="ZH-T")
mydata <- mydata[,111:117]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 11

library(psych)
omega(mydata) # alpha = .87, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.84; eigenvalue 2 = .34

rho <- polychoric(mydata)$rho
# No polychoric needed for 11 categories
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.61; eigenvalue 2 = .45

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.55, RMSEA=.189, RMSR=.07, TLI=.827

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.224, RMSR=.09, TLI=.738

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 1 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.089, RMSR=.03, TLI=.962
#      MR1   MR2
#MR1  1.00  0.68
#MR2  0.68  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.62, RMSEA=.089, RMSR=.03, TLI=.959
#      MR1   MR2
#MR1  1.00  0.64
#MR2  0.64  1.00

# Single factor model lavaan
UNImodel= '
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.983       0.943
#Tucker-Lewis Index (TLI)                       0.975       0.914
#Robust Comparative Fit Index (CFI)                         0.806
#Robust Tucker-Lewis Index (TLI)                            0.710
#RMSEA                                          0.176       0.241
#Robust RMSEA                                               0.240
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .673

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.888       0.871
#Tucker-Lewis Index (TLI)                       0.832       0.806
#Robust Comparative Fit Index (CFI)                         0.891
#Robust Tucker-Lewis Index (TLI)                            0.837
#RMSEA                                          0.188       0.161
#Robust RMSEA                                               0.185
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .658

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.984       0.940
#Tucker-Lewis Index (TLI)                       0.974       0.903
#Robust Comparative Fit Index (CFI)                         0.808
#Robust Tucker-Lewis Index (TLI)                            0.690
#RMSEA                                          0.181       0.256
#Robust RMSEA                                               0.247
#SRMR                                           0.079       0.079

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .673

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.716    0.114    6.277    0.000    0.716    0.716

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.891       0.876
#Tucker-Lewis Index (TLI)                       0.823       0.800
#Robust Comparative Fit Index (CFI)                         0.893
#Robust Tucker-Lewis Index (TLI)                            0.828
#RMSEA                                          0.193       0.163
#Robust RMSEA                                               0.189
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .684

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.640    0.175    3.656    0.000    0.640    0.640

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.822       0.813
#Tucker-Lewis Index (TLI)                       0.733       0.719
#Robust Comparative Fit Index (CFI)                         0.826
#Robust Tucker-Lewis Index (TLI)                            0.738
#RMSEA                                          0.237       0.193
#Robust RMSEA                                               0.233
#SRMR                                           0.205       0.205

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .650


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ trust_1+trust_2+trust_3+trust_6+trust_4
 factor2 =~ trust_7+trust_5
 general =~ trust_1+trust_2+trust_3+trust_4+trust_5+trust_6+trust_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.970       0.958
#Tucker-Lewis Index (TLI)                       0.909       0.875
#Robust Comparative Fit Index (CFI)                         0.970
#Robust Tucker-Lewis Index (TLI)                            0.910
#RMSEA                                          0.138       0.129
#Robust RMSEA                                               0.137
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .703

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    trust_1           0.128    0.533    0.240    0.810    0.128    0.054
#    trust_2           0.926    0.413    2.241    0.025    0.926    0.432
#    trust_3           0.954    0.395    2.413    0.016    0.954    0.503
#    trust_6          -0.463    0.739   -0.627    0.531   -0.463   -0.182
#    trust_4           0.347    0.393    0.884    0.377    0.347    0.169
#  factor2 =~                                                            
#    trust_7           0.934    1.352    0.691    0.490    0.934    0.523
#    trust_5           0.550    1.101    0.499    0.618    0.550    0.231
#  general =~                                                            
#    trust_1           2.084    0.166   12.522    0.000    2.084    0.872
#    trust_2           1.540    0.302    5.099    0.000    1.540    0.718
#    trust_3           1.498    0.312    4.807    0.000    1.498    0.790
#    trust_4           1.601    0.221    7.233    0.000    1.601    0.778
#    trust_5           0.221    0.277    0.798    0.425    0.221    0.093
#    trust_6           2.500    0.177   14.099    0.000    2.500    0.981
#    trust_7           0.947    0.155    6.127    0.000    0.947    0.530

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.8189552      0.4761905      0.9093879      0.8522511 
#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1269072 0.10988234 0.8730928 0.9462830 0.04981832 0.3883145 0.9091723
#factor2 0.5304610 0.07116246 0.4695390 0.4083864 0.24277617 0.3017621 0.6319865
#general 0.8189552 0.81895521 0.8189552 0.9093879 0.85225115 0.9706679 0.9955847





