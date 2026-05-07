################################################################
### COVIDSTRESS database
### 
### Blackburn et al. (2022) https://www.nature.com/articles/s41597-022-01383-6
### data available at https://osf.io/36tsd/
###
### Compliance Scale (reappraisal and suppression)
###
### Full dataset to determine 2-factor solution
###

Final_COVIDiSTRESS_Vol2_cleaned <- read.csv("Final_COVIDiSTRESS_Vol2_cleaned.csv")
colnames(Final_COVIDiSTRESS_Vol2_cleaned)
table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)

mydata <- Final_COVIDiSTRESS_Vol2_cleaned[,80:87]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
fit1 <- fa(mydata,2)
fit1
diagram(fit1)

EGAmodel= '
 factor1 =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3
 factor2 =~ compliance_5+compliance_7
'


################################################################
### COVIDSTRESS database
### 
### UserLanguage = ES-ES (Spanish)
### 
###

mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="ES-ES")
mydata <- mydata[,80:87]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .81, omega T = .86

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.93; eigenvalue 2 = .47

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.84; eigenvalue 2 = .47

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.37, RMSEA=.126, RMSR=.08, TLI=.806

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.185, RMSR=.08, TLI=.76

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
# %variance explained=.46, RMSEA=.034, RMSR=.02, TLI=.986
#      MR1   MR2
#MR1  1.00  0.45
#MR2  0.45  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.58, RMSEA=.088, RMSR=.02, TLI=.946
#      MR1   MR2
#MR1  1.00  0.49
#MR2  0.49  1.00


# Single factor model lavaan
UNImodel= '
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7

'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.930
#Tucker-Lewis Index (TLI)                       0.971       0.902
#Robust Comparative Fit Index (CFI)                         0.828
#Robust Tucker-Lewis Index (TLI)                            0.760
#RMSEA                                          0.094       0.132
#Robust RMSEA                                               0.186
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .555

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.862       0.847
#Tucker-Lewis Index (TLI)                       0.806       0.786
#Robust Comparative Fit Index (CFI)                         0.866
#Robust Tucker-Lewis Index (TLI)                            0.812
#RMSEA                                          0.126       0.106
#Robust RMSEA                                               0.124
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .366

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.990       0.961
#Tucker-Lewis Index (TLI)                       0.985       0.943
#Robust Comparative Fit Index (CFI)                         0.886
#Robust Tucker-Lewis Index (TLI)                            0.833
#RMSEA                                          0.069       0.101
#Robust RMSEA                                               0.155
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .511

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.730    0.032   22.907    0.000    0.730    0.730

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.916       0.913
#Tucker-Lewis Index (TLI)                       0.876       0.871
#Robust Comparative Fit Index (CFI)                         0.921
#Robust Tucker-Lewis Index (TLI)                            0.883
#RMSEA                                          0.101       0.082
#Robust RMSEA                                               0.098
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .382

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.640    0.070    9.199    0.000    0.640    0.640

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.733       0.733
#Tucker-Lewis Index (TLI)                       0.626       0.626
#Robust Comparative Fit Index (CFI)                         0.738
#Robust Tucker-Lewis Index (TLI)                            0.633
#RMSEA                                          0.176       0.139
#Robust RMSEA                                               0.173
#SRMR                                           0.189       0.189

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .419


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3
 factor2 =~ compliance_5+compliance_7
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.972
#Tucker-Lewis Index (TLI)                       0.953       0.936
#Robust Comparative Fit Index (CFI)                         0.981
#Robust Tucker-Lewis Index (TLI)                            0.955
#RMSEA                                          0.062       0.058
#Robust RMSEA                                               0.061
#SRMR                                           0.033       0.033

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .496

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    compliance_2      0.736    0.145    5.079    0.000    0.736    0.717
#    compliance_4      0.681    0.529    1.286    0.198    0.681    0.407
#    compliance_1      0.700    0.294    2.381    0.017    0.700    0.493
#    compliance_6      0.550    0.447    1.229    0.219    0.550    0.287
#    compliance_8      0.568    0.347    1.636    0.102    0.568    0.374
#    compliance_3      0.925    0.277    3.340    0.001    0.925    0.555
#  factor2 =~                                                            
#    compliance_5     -0.266    0.128   -2.077    0.038   -0.266   -0.131
#    compliance_7     -0.943    7.369   -0.128    0.898   -0.943   -0.444
#  general =~                                                            
#    compliance_2      0.264    0.124    2.128    0.033    0.264    0.257
#    compliance_4      0.981    0.462    2.123    0.034    0.981    0.586
#    compliance_1      0.458    0.294    1.561    0.119    0.458    0.322
#    compliance_6      0.643    0.399    1.612    0.107    0.643    0.336
#    compliance_8      0.739    0.309    2.394    0.017    0.739    0.487
#    compliance_3      0.792    0.300    2.637    0.008    0.792    0.475
#    compliance_5      1.970    0.898    2.194    0.028    1.970    0.971
#    compliance_7     -1.134    0.557   -2.037    0.042   -1.134   -0.535

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.5817046      0.4285714      0.8067648      0.4046897 
#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5715407 0.36442529 0.4284593 0.8028928 0.4571678 0.6944008 0.8551115
#factor2 0.1486560 0.05387008 0.8513440 0.4841921 0.3073415 0.2085941 0.5874555
#general 0.5817046 0.58170463 0.5817046 0.8067648 0.4046897 0.9487960 0.9771606





################################################################
### COVIDSTRESS database
### 
### UserLanguage = FI (Finnish)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="FI")
mydata <- mydata[,80:87]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .81, omega T = .87

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.99; eigenvalue 2 = .55

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 3.72; eigenvalue 2 = .54

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.37, RMSEA=.149, RMSR=.09, TLI=.754

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.202, RMSR=.09, TLI=.708

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
# %variance explained=.48, RMSEA=.056, RMSR=.02, TLI=.965
#      MR1   MR2
#MR1  1.00  0.45
#MR2  0.45  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.127, RMSR=.04, TLI=.885
#      MR1   MR2
#MR1  1.00  0.51
#MR2  0.51  1.00


# Single factor model lavaan
UNImodel= '
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7

'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.968       0.905
#Tucker-Lewis Index (TLI)                       0.955       0.867
#Robust Comparative Fit Index (CFI)                         0.792
#Robust Tucker-Lewis Index (TLI)                            0.709
#RMSEA                                          0.118       0.157
#Robust RMSEA                                               0.203
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .477

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.826       0.806
#Tucker-Lewis Index (TLI)                       0.756       0.728
#Robust Comparative Fit Index (CFI)                         0.827
#Robust Tucker-Lewis Index (TLI)                            0.757
#RMSEA                                          0.149       0.137
#Robust RMSEA                                               0.148
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .404

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.982       0.944
#Tucker-Lewis Index (TLI)                       0.973       0.918
#Robust Comparative Fit Index (CFI)                         0.854
#Robust Tucker-Lewis Index (TLI)                            0.784
#RMSEA                                          0.091       0.123
#Robust RMSEA                                               0.175
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .480

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.615    0.045   13.805    0.000    0.615    0.615

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.887       0.875
#Tucker-Lewis Index (TLI)                       0.833       0.816
#Robust Comparative Fit Index (CFI)                         0.888
#Robust Tucker-Lewis Index (TLI)                            0.835
#RMSEA                                          0.123       0.113
#Robust RMSEA                                               0.122
#SRMR                                           0.065       0.065

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .415

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.439    0.099    4.420    0.000    0.439    0.439

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.732       0.702
#Tucker-Lewis Index (TLI)                       0.624       0.583
#Robust Comparative Fit Index (CFI)                         0.732
#Robust Tucker-Lewis Index (TLI)                            0.625
#RMSEA                                          0.185       0.170
#Robust RMSEA                                               0.184
#SRMR                                           0.183       0.183

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .433


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3
 factor2 =~ compliance_5+compliance_7
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.975       0.972
#Tucker-Lewis Index (TLI)                       0.942       0.934
#Robust Comparative Fit Index (CFI)                         0.976
#Robust Tucker-Lewis Index (TLI)                            0.944
#RMSEA                                          0.073       0.068
#Robust RMSEA                                               0.071
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .519

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    compliance_2      1.120    0.101   11.039    0.000    1.120    0.696
#    compliance_4      0.452    0.106    4.242    0.000    0.452    0.262
#    compliance_1      0.862    0.099    8.680    0.000    0.862    0.701
#    compliance_6      0.536    0.130    4.143    0.000    0.536    0.252
#    compliance_8      0.755    0.101    7.437    0.000    0.755    0.442
#    compliance_3      0.340    0.091    3.739    0.000    0.340    0.162
#  factor2 =~                                                            
#    compliance_5      0.673    0.080    8.389    0.000    0.673    0.307
#    compliance_7     -1.136    0.199   -5.715    0.000   -1.136   -0.549
#  general =~                                                            
#    compliance_2      0.675    0.098    6.885    0.000    0.675    0.420
#    compliance_4      1.055    0.084   12.526    0.000    1.055    0.613
#    compliance_1      0.394    0.075    5.239    0.000    0.394    0.320
#    compliance_6      1.062    0.114    9.290    0.000    1.062    0.498
#    compliance_8      1.086    0.095   11.438    0.000    1.086    0.636
#    compliance_3      1.064    0.084   12.616    0.000    1.064    0.508
#    compliance_5      1.836    0.107   17.226    0.000    1.836    0.836
#    compliance_7     -0.545    0.140   -3.892    0.000   -0.545   -0.263

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.5749324      0.4285714      0.8291067      0.5520905 
#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.4594467 0.32766543 0.5405533 0.8312747 0.34391802 0.6986597 0.8329675
#factor2 0.3395863 0.09740217 0.6604137 0.3167243 0.04803833 0.3486153 0.6038527
#general 0.5749324 0.57493240 0.5749324 0.8291067 0.55209050 0.8243045 0.8919618





################################################################
### COVIDSTRESS database
### 
### UserLanguage = PT-BR (Portuguese)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="PT-BR")
mydata <- mydata[,80:87]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .74, omega T = .8

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.36; eigenvalue 2 = .5

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 3.29; eigenvalue 2 = .64

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.29, RMSEA=.117, RMSR=.08, TLI=.767

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.21, RMSR=.11, TLI=.629

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
# %variance explained=.4, RMSEA=.038, RMSR=.03, TLI=.975
#      MR1   MR2
#MR1  1.00  0.32
#MR2  0.32  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.083, RMSR=.03, TLI=.943
#      MR1   MR2
#MR1  1.00  0.38
#MR2  0.38  1.00


# Single factor model lavaan
UNImodel= '
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7

'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.957       0.894
#Tucker-Lewis Index (TLI)                       0.940       0.851
#Robust Comparative Fit Index (CFI)                         0.749
#Robust Tucker-Lewis Index (TLI)                            0.649
#RMSEA                                          0.111       0.142
#Robust RMSEA                                               0.206
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .417

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.837       0.801
#Tucker-Lewis Index (TLI)                       0.772       0.722
#Robust Comparative Fit Index (CFI)                         0.839
#Robust Tucker-Lewis Index (TLI)                            0.774
#RMSEA                                          0.117       0.108
#Robust RMSEA                                               0.115
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .245

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.974       0.933
#Tucker-Lewis Index (TLI)                       0.962       0.901
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.734
#RMSEA                                          0.088       0.116
#Robust RMSEA                                               0.179
#SRMR                                           0.084       0.084

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .459

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.702    0.046   15.403    0.000    0.702    0.702

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.885       0.865
#Tucker-Lewis Index (TLI)                       0.831       0.801
#Robust Comparative Fit Index (CFI)                         0.888
#Robust Tucker-Lewis Index (TLI)                            0.836
#RMSEA                                          0.101       0.091
#Robust RMSEA                                               0.098
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .204

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.714    0.100    7.162    0.000    0.714    0.714

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.604       0.576
#Tucker-Lewis Index (TLI)                       0.446       0.406
#Robust Comparative Fit Index (CFI)                         0.607
#Robust Tucker-Lewis Index (TLI)                            0.450
#RMSEA                                          0.182       0.157
#Robust RMSEA                                               0.179
#SRMR                                           0.178       0.178

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .344


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3
 factor2 =~ compliance_5+compliance_7
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.965       0.964
#Tucker-Lewis Index (TLI)                       0.919       0.916
#Robust Comparative Fit Index (CFI)                         0.969
#Robust Tucker-Lewis Index (TLI)                            0.928
#RMSEA                                          0.070       0.059
#Robust RMSEA                                               0.065
#SRMR                                           0.041       0.041

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .422

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    compliance_2      0.454    0.196    2.313    0.021    0.454    0.816
#    compliance_4      0.138    0.102    1.361    0.173    0.138    0.083
#    compliance_1      0.168    0.119    1.414    0.157    0.168    0.152
#    compliance_6      0.317    0.136    2.331    0.020    0.317    0.201
#    compliance_8      0.232    0.138    1.679    0.093    0.232    0.167
#    compliance_3      0.615    0.268    2.296    0.022    0.615    0.421
#  factor2 =~                                                            
#    compliance_5      0.438    0.169    2.590    0.010    0.438    0.241
#    compliance_7     -1.010    0.230   -4.400    0.000   -1.010   -0.491
#  general =~                                                            
#    compliance_2      0.137    0.072    1.901    0.057    0.137    0.247
#    compliance_4      1.141    0.116    9.871    0.000    1.141    0.684
#    compliance_1      0.403    0.093    4.339    0.000    0.403    0.365
#    compliance_6      0.322    0.131    2.464    0.014    0.322    0.204
#    compliance_8      0.873    0.132    6.630    0.000    0.873    0.626
#    compliance_3      0.684    0.140    4.869    0.000    0.684    0.468
#    compliance_5      1.480    0.124   11.909    0.000    1.480    0.815
#    compliance_7     -0.880    0.174   -5.061    0.000   -0.880   -0.428

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6354162      0.4285714      0.7283421      0.5248040 
#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.4171729 0.27657793 0.5828271 0.7297373 0.24414244 0.6975467 0.8432329
#factor2 0.2611307 0.08800585 0.7388693 0.1988577 0.05859831 0.2751627 0.5436695
#general 0.6354162 0.63541622 0.6354162 0.7283421 0.52480404 0.8100301 0.8897779





################################################################
### COVIDSTRESS database
### 
### UserLanguage = SK (Slovak)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="SK")
mydata <- mydata[,80:87]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .82, omega T = .86

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.91; eigenvalue 2 = .55

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 3.29; eigenvalue 2 = .64

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.36, RMSEA=.137, RMSR=.09, TLI=.777

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.204, RMSR=.1, TLI=.684

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
# %variance explained=.47, RMSEA=.066, RMSR=.03, TLI=.948
#      MR1   MR2
#MR1  1.00  0.49
#MR2  0.49  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.138, RMSR=.04, TLI=.855
#      MR1   MR2
#MR1  1.00  0.58
#MR2  0.58  1.00


# Single factor model lavaan
UNImodel= '
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.975       0.925
#Tucker-Lewis Index (TLI)                       0.964       0.895
#Robust Comparative Fit Index (CFI)                         0.776
#Robust Tucker-Lewis Index (TLI)                            0.687
#RMSEA                                          0.118       0.158
#Robust RMSEA                                               0.206
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .418

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.846       0.852
#Tucker-Lewis Index (TLI)                       0.784       0.793
#Robust Comparative Fit Index (CFI)                         0.858
#Robust Tucker-Lewis Index (TLI)                            0.801
#RMSEA                                          0.136       0.107
#Robust RMSEA                                               0.129
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .289

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.936
#Tucker-Lewis Index (TLI)                       0.968       0.906
#Robust Comparative Fit Index (CFI)                         0.814
#Robust Tucker-Lewis Index (TLI)                            0.725
#RMSEA                                          0.111       0.150
#Robust RMSEA                                               0.193
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .389

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.809    0.040   20.202    0.000    0.809    0.809

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.863       0.871
#Tucker-Lewis Index (TLI)                       0.797       0.810
#Robust Comparative Fit Index (CFI)                         0.875
#Robust Tucker-Lewis Index (TLI)                            0.816
#RMSEA                                          0.132       0.103
#Robust RMSEA                                               0.124
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .31

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.833    0.069   12.051    0.000    0.833    0.833

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.586       0.588
#Tucker-Lewis Index (TLI)                       0.421       0.423
#Robust Comparative Fit Index (CFI)                         0.594
#Robust Tucker-Lewis Index (TLI)                            0.431
#RMSEA                                          0.223       0.179
#Robust RMSEA                                               0.218
#SRMR                                           0.220       0.220

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .356


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3
 factor2 =~ compliance_5+compliance_7
 general =~ compliance_2+compliance_4+compliance_1+
            compliance_6+compliance_8+compliance_3+
            compliance_5+compliance_7
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.969       0.963
#Tucker-Lewis Index (TLI)                       0.928       0.913
#Robust Comparative Fit Index (CFI)                         0.972
#Robust Tucker-Lewis Index (TLI)                            0.934
#RMSEA                                          0.079       0.069
#Robust RMSEA                                               0.074
#SRMR                                           0.037       0.037

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .456

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    compliance_2      0.781    0.209    3.729    0.000    0.781    0.721
#    compliance_4      0.423    0.205    2.057    0.040    0.423    0.233
#    compliance_1      0.538    0.152    3.551    0.000    0.538    0.529
#    compliance_6      0.658    0.215    3.055    0.002    0.658    0.402
#    compliance_8      0.329    0.145    2.278    0.023    0.329    0.231
#    compliance_3      0.136    0.220    0.619    0.536    0.136    0.063
#  factor2 =~                                                            
#    compliance_5      0.069    0.357    0.193    0.847    0.069    0.033
#    compliance_7     -0.876    0.048  -18.284    0.000   -0.876   -0.452
#  general =~                                                            
#    compliance_2      0.416    0.098    4.229    0.000    0.416    0.384
#    compliance_4      1.296    0.124   10.452    0.000    1.296    0.713
#    compliance_1      0.255    0.078    3.269    0.001    0.255    0.251
#    compliance_6      0.592    0.135    4.397    0.000    0.592    0.362
#    compliance_8      0.713    0.105    6.772    0.000    0.713    0.500
#    compliance_3      1.440    0.139   10.361    0.000    1.440    0.671
#    compliance_5      1.901    0.141   13.503    0.000    1.901    0.913
#    compliance_7     -0.975    0.169   -5.776    0.000   -0.975   -0.503

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6734606      0.4285714      0.7941009      0.5459338 
#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4091241 0.27412482 0.5908759 0.7943841 0.2891324 0.6407278 0.8201870
#factor2 0.1588457 0.05241457 0.8411543 0.3267504 0.1665468 0.2047911 0.5038046
#general 0.6734606 0.67346061 0.6734606 0.7941009 0.5459338 0.8880971 0.9379404





