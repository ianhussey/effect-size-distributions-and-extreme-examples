################################################################
### COVIDSTRESS database
### 
### Blackburn et al. (2022) https://www.nature.com/articles/s41597-022-01383-6
### data available at https://osf.io/36tsd/
###
### Emotion regulation Scale (reappraisal and suppression)
###
### Full dataset to determine 2-factor solution
###

Final_COVIDiSTRESS_Vol2_cleaned <- read.csv("Final_COVIDiSTRESS_Vol2_cleaned.csv")
colnames(Final_COVIDiSTRESS_Vol2_cleaned)
table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)

mydata <- Final_COVIDiSTRESS_Vol2_cleaned[,186:193]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
fit1 <- fa(mydata,2)
fit1
diagram(fit1)

EGAmodel= '
 factor1 =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral
 factor2 =~ emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'


################################################################
### COVIDSTRESS database
### 
### UserLanguage = EN (American English)
### 
###

mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="EN")
mydata <- mydata[,186:193]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .73, omega T = .82

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.1; eigenvalue 2 = .96

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.65; eigenvalue 2 = 1.04

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.26, RMSEA=.17, RMSR=.15, TLI=.561

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.197, RMSR=.16, TLI=.625

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.007, RMSR=.01, TLI=.999
#      MR1   MR2
#MR1  1.00  0.18
#MR2  0.18  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.025, RMSR=.01, TLI=.994
#      MR1   MR2
#MR1  1.00  0.23
#MR2  0.23  1.00


# Single factor model lavaan
UNImodel= '
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.907       0.809
#Tucker-Lewis Index (TLI)                       0.869       0.733
#Robust Comparative Fit Index (CFI)                         0.707
#Robust Tucker-Lewis Index (TLI)                            0.589
#RMSEA                                          0.252       0.282
#Robust RMSEA                                               0.207
#SRMR                                           0.149       0.149

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .229

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.707       0.680
#Tucker-Lewis Index (TLI)                       0.589       0.552
#Robust Comparative Fit Index (CFI)                         0.708
#Robust Tucker-Lewis Index (TLI)                            0.591
#RMSEA                                          0.165       0.141
#Robust RMSEA                                               0.164
#SRMR                                           0.133       0.133

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .124

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.998       0.995
#Tucker-Lewis Index (TLI)                       0.997       0.992
#Robust Comparative Fit Index (CFI)                         0.991
#Robust Tucker-Lewis Index (TLI)                            0.986
#RMSEA                                          0.038       0.048
#Robust RMSEA                                               0.038
#SRMR                                           0.032       0.032

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .483

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.261    0.030    8.691    0.000    0.261    0.261

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.993
#Tucker-Lewis Index (TLI)                       0.987       0.990
#Robust Comparative Fit Index (CFI)                         0.994
#Robust Tucker-Lewis Index (TLI)                            0.991
#RMSEA                                          0.029       0.021
#Robust RMSEA                                               0.024
#SRMR                                           0.031       0.031

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .417

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.196    0.042    4.671    0.000    0.196    0.196

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.980
#Tucker-Lewis Index (TLI)                       0.972       0.972
#Robust Comparative Fit Index (CFI)                         0.982
#Robust Tucker-Lewis Index (TLI)                            0.975
#RMSEA                                          0.043       0.035
#Robust RMSEA                                               0.040
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .416


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral
 factor2 =~ emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.998       0.998
#Tucker-Lewis Index (TLI)                       0.995       0.994
#Robust Comparative Fit Index (CFI)                         0.998
#Robust Tucker-Lewis Index (TLI)                            0.996
#RMSEA                                          0.018       0.016
#Robust RMSEA                                               0.016
#SRMR                                           0.013       0.013

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .417

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    emotns_8_0ntrl    1.434    0.315    4.555    0.000    1.434    0.728
#    emotns_6_0ntrl    1.394    0.409    3.407    0.001    1.394    0.691
#    emotns_7_0ntrl    1.493    0.538    2.776    0.006    1.493    0.747
#    emotns_3_0ntrl    0.630    1.005    0.627    0.531    0.630    0.324
#  factor2 =~                                                            
#    emotns_5_0ntrl    1.123    0.712    1.577    0.115    1.123    0.557
#    emotns_1_0ntrl    1.118    0.640    1.748    0.081    1.118    0.566
#    emotns_4_0ntrl    1.074    0.649    1.655    0.098    1.074    0.549
#    emotns_2_0ntrl    0.525    0.579    0.907    0.364    0.525    0.300
#  general =~                                                            
#    emotns_8_0ntrl    0.519    0.944    0.550    0.582    0.519    0.264
#    emotns_6_0ntrl    0.617    0.917    0.673    0.501    0.617    0.306
#    emotns_7_0ntrl    0.585    1.266    0.462    0.644    0.585    0.293
#    emotns_3_0ntrl    0.826    1.483    0.557    0.577    0.826    0.424
#    emotns_5_0ntrl    0.692    1.153    0.600    0.548    0.692    0.343
#    emotns_1_0ntrl    0.565    1.214    0.466    0.642    0.565    0.286
#    emotns_4_0ntrl    0.627    1.127    0.556    0.578    0.627    0.321
#    emotns_2_0ntrl    0.439    0.761    0.577    0.564    0.439    0.251

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.2278003      0.5714286      0.7828644      0.2976018 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.7956245 0.4788262 0.2043755 0.8049970 0.6352876 0.7735505 0.8519567
#factor2 0.7367943 0.2933735 0.2632057 0.6712591 0.4897979 0.5923640 0.7336124
#general 0.2278003 0.2278003 0.2278003 0.7828644 0.2976018 0.4724811 0.5758011





################################################################
### COVIDSTRESS database
### 
### UserLanguage = ES-GU (Spanish)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="ES-GU")
mydata <- mydata[,186:193]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .63, omega T = .78

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 1.89; eigenvalue 2 = .87

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.36; eigenvalue 2 = .82

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.24, RMSEA=.159, RMSR=.14, TLI=.549

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.29, RMSEA=.182, RMSR=.14, TLI=.617

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
# %variance explained=.39, RMSEA=.079, RMSR=.05, TLI=.889
#      MR1   MR2
#MR1  1.00  -0.10
#MR2  -0.10  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.133, RMSR=.06, TLI=.796
#      MR1   MR2
#MR1  1.00  -.15
#MR2  -.15  1.00


# Single factor model lavaan
UNImodel= '
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.908       0.850
#Tucker-Lewis Index (TLI)                       0.871       0.790
#Robust Comparative Fit Index (CFI)                         0.722
#Robust Tucker-Lewis Index (TLI)                            0.610
#RMSEA                                          0.216       0.231
#Robust RMSEA                                               0.188
#SRMR                                           0.127       0.127

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .136

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.689       0.684
#Tucker-Lewis Index (TLI)                       0.565       0.557
#Robust Comparative Fit Index (CFI)                         0.699
#Robust Tucker-Lewis Index (TLI)                            0.579
#RMSEA                                          0.158       0.132
#Robust RMSEA                                               0.153
#SRMR                                           0.126       0.126

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .033

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.978       0.966
#Tucker-Lewis Index (TLI)                       0.968       0.950
#Robust Comparative Fit Index (CFI)                         0.877
#Robust Tucker-Lewis Index (TLI)                            0.818
#RMSEA                                          0.108       0.112
#Robust RMSEA                                               0.129
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .433

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.221    0.068   -3.236    0.001   -0.221   -0.221

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.919       0.919
#Tucker-Lewis Index (TLI)                       0.881       0.881
#Robust Comparative Fit Index (CFI)                         0.930
#Robust Tucker-Lewis Index (TLI)                            0.896
#RMSEA                                          0.083       0.068
#Robust RMSEA                                               0.076
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .369

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.119    0.096   -1.241    0.215   -0.119   -0.119

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.917       0.916
#Tucker-Lewis Index (TLI)                       0.884       0.883
#Robust Comparative Fit Index (CFI)                         0.928
#Robust Tucker-Lewis Index (TLI)                            0.899
#RMSEA                                          0.082       0.068
#Robust RMSEA                                               0.075
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .369


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral
 factor2 =~ emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.940       0.772
#Tucker-Lewis Index (TLI)                       0.861       0.467
#Robust Comparative Fit Index (CFI)                         0.920
#Robust Tucker-Lewis Index (TLI)                            0.813
#RMSEA                                          0.089       0.145
#Robust RMSEA                                               0.102
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .264

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.20724997     0.57142857     0.73759587     0.02229881 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#factor1 0.7375474 0.5261382 0.26245264 0.8792733 0.75057955 1.9914822 1.2696359
#factor2 0.9301340 0.2666118 0.06986597 0.6017498 0.56869333 0.6726364 0.8322619
#general 0.2072500 0.2072500 0.20724997 0.7375959 0.02229881 0.6527743 1.0780581








################################################################
### COVIDSTRESS database
### 
### UserLanguage = IT (Italian)
### 
###

### Not included because problems with the BIF model

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="IT")
mydata <- mydata[,186:193]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .72, omega T = .81

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.04; eigenvalue 2 = 1.14

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.39; eigenvalue 2 = 1.22

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.209, RMSR=.17, TLI=.364

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.23, RMSR=.18, TLI=.44

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
# %variance explained=.45, RMSEA=.023, RMSR=.03, TLI=.992
#      MR1   MR2
#MR1  1.00   0.20
#MR2  0.20   1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.082, RMSR=.04, TLI=.929
#      MR1   MR2
#MR1  1.00  0.19
#MR2  0.19  1.00


# Single factor model lavaan
UNImodel= '
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.840       0.716
#Tucker-Lewis Index (TLI)                       0.777       0.603
#Robust Comparative Fit Index (CFI)                         0.584
#Robust Tucker-Lewis Index (TLI)                            0.417
#RMSEA                                          0.282       0.296
#Robust RMSEA                                               0.239
#SRMR                                           0.172       0.172

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .224

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.978       0.954
#Tucker-Lewis Index (TLI)                       0.968       0.932
#Robust Comparative Fit Index (CFI)                         0.922
#Robust Tucker-Lewis Index (TLI)                            0.885
#RMSEA                                          0.106       0.123
#Robust RMSEA                                               0.106
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .477

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.247    0.061    4.070    0.000    0.247    0.247

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.952       0.945
#Tucker-Lewis Index (TLI)                       0.933       0.923
#Robust Comparative Fit Index (CFI)                         0.913
#Robust Tucker-Lewis Index (TLI)                            0.878
#RMSEA                                          0.155       0.130
#Robust RMSEA                                               0.109
#SRMR                                           0.109       0.109

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .47


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral
 factor2 =~ emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
# warning
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.994       0.998
#Tucker-Lewis Index (TLI)                       0.985       0.994
#Robust Comparative Fit Index (CFI)                         0.998
#Robust Tucker-Lewis Index (TLI)                            0.995
#RMSEA                                          0.032       0.017
#Robust RMSEA                                               0.018
#SRMR                                           0.028       0.028

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .517

#Latent Variables:
#                   Estimate   Std.Err  z-value  P(>|z|)   Std.lv   Std.all
#  factor1 =~                                                              
#    emotns_8_0ntrl     1.604    0.121   13.295    0.000     1.604    0.779
#    emotns_6_0ntrl     1.565    0.138   11.349    0.000     1.565    0.715
#    emotns_7_0ntrl     1.593    0.111   14.377    0.000     1.593    0.779
#    emotns_3_0ntrl     0.492    0.133    3.687    0.000     0.492    0.255
#  factor2 =~                                                              
#    emotns_5_0ntrl     0.008    0.011    0.712    0.476     0.008    0.004
#    emotns_1_0ntrl     0.009    0.010    0.851    0.395     0.009    0.004
#    emotns_4_0ntrl     0.013    0.009    1.412    0.158     0.013    0.007
#    emotns_2_0ntrl    32.882    0.079  417.490    0.000    32.882   21.379
#  general =~                                                              
#    emotns_8_0ntrl     0.341    0.167    2.041    0.041     0.341    0.166
#    emotns_6_0ntrl     0.170    0.208    0.815    0.415     0.170    0.077
#    emotns_7_0ntrl     0.521    0.158    3.300    0.001     0.521    0.255
#    emotns_3_0ntrl     0.734    0.151    4.847    0.000     0.734    0.380
#    emotns_5_0ntrl     1.360    0.127   10.700    0.000     1.360    0.743
#    emotns_1_0ntrl     1.230    0.148    8.336    0.000     1.230    0.628
#    emotns_4_0ntrl     1.050    0.136    7.709    0.000     1.050    0.559
#    emotns_2_0ntrl     0.596    0.270    2.208    0.027     0.596    0.388

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#   0.003585112    0.571428571   21.741491283    0.467840687 
#$FactorLevelIndices
#         ECV_SS      ECV_SG      ECV_GS      Omega     OmegaH           H         FD
#factor1 0.880682363 0.003885639 0.119317637  0.7843248  0.6999483   0.8077004  0.9001723
#factor2 0.996927772 0.992529249 0.003072228 53.6991416 53.0764771 475.2719906 24.1157428
#general 0.003585112 0.003585112 0.003585112 21.7414913  0.4678407   0.7356792  0.8475070






################################################################
### COVIDSTRESS database
### 
### UserLanguage = RU (Russian)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="RU")
mydata <- mydata[,186:193]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .78, omega T = .83

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.54; eigenvalue 2 = .59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 3.07; eigenvalue 2 = .64

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.145, RMSR=.09, TLI=.697

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.178, RMSR=.1, TLI=.679

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
# %variance explained=.43, RMSEA=.04, RMSR=.02, TLI=.977
#      MR1   MR2
#MR1  1.00   0.48
#MR2  0.48   1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.5, RMSEA=.057, RMSR=.02, TLI=.967
#      MR1   MR2
#MR1  1.00  0.51
#MR2  0.51  1.00


# Single factor model lavaan
UNImodel= '
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.939       0.826
#Tucker-Lewis Index (TLI)                       0.915       0.756
#Robust Comparative Fit Index (CFI)                         0.756
#Robust Tucker-Lewis Index (TLI)                            0.659
#RMSEA                                          0.165       0.207
#Robust RMSEA                                               0.184
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .446

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.784       0.766
#Tucker-Lewis Index (TLI)                       0.697       0.672
#Robust Comparative Fit Index (CFI)                         0.785
#Robust Tucker-Lewis Index (TLI)                            0.698
#RMSEA                                          0.145       0.127
#Robust RMSEA                                               0.144
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .349

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.976       0.926
#Tucker-Lewis Index (TLI)                       0.964       0.891
#Robust Comparative Fit Index (CFI)                         0.906
#Robust Tucker-Lewis Index (TLI)                            0.861
#RMSEA                                          0.107       0.139
#Robust RMSEA                                               0.117
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .536

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.658    0.015   42.762    0.000    0.658    0.658

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.922       0.917
#Tucker-Lewis Index (TLI)                       0.885       0.877
#Robust Comparative Fit Index (CFI)                         0.923
#Robust Tucker-Lewis Index (TLI)                            0.887
#RMSEA                                          0.089       0.078
#Robust RMSEA                                               0.088
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .440

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.577    0.031   18.891    0.000    0.577    0.577

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.819       0.804
#Tucker-Lewis Index (TLI)                       0.747       0.726
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.748
#RMSEA                                          0.132       0.116
#Robust RMSEA                                               0.132
#SRMR                                           0.176       0.176

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .414


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral
 factor2 =~ emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
 general =~ emotions_8_0neutral+emotions_6_0neutral+emotions_7_0neutral+
            emotions_3_0neutral+
            emotions_5_0neutral+emotions_1_0neutral+emotions_4_0neutral+
            emotions_2_0neutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.992       0.993
#Tucker-Lewis Index (TLI)                       0.982       0.984
#Robust Comparative Fit Index (CFI)                         0.993
#Robust Tucker-Lewis Index (TLI)                            0.984
#RMSEA                                          0.036       0.028
#Robust RMSEA                                               0.033
#SRMR                                           0.015       0.015

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .391

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#   0.009864865    0.571428571    9.433617172    0.655961093 
#
#$FactorLevelIndices
#             ECV_SS      ECV_SG      ECV_GS      Omega     OmegaH           H         FD
#factor1 0.589044589 0.005153227 0.410955411  0.7716872  0.4289368   0.6701781  0.8203592
#factor2 0.993675024 0.984981907 0.006324976 27.2139312 26.5574978 220.9196636 18.5504843
#general 0.009864865 0.009864865 0.009864865  9.4336172  0.6559611   0.7694012  0.8560208







