################################################################
### Emotional intelligence
### 
### One factor or several?
### 



################################################################
### 
### Brienza et al (2018) https://psycnet.apa.org/record/2017-42043-001
### Wong and Law’s (2002) 16-item Emotional Intelligence Scale measuring 4 dimensions of emotional intelligence
### data https://osf.io/7jhdn

library(haven)
EI_Brienza <- read_sav("EI_Brienza.sav")
colnames(EI_Brienza)
mydata <- EI_Brienza[,283:298]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
Names <- colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .91, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 6.52; eigenvalue 2 = 1.57

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 7.28; eigenvalue 2 = 1.64

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.41, RMSEA=.221, RMSR=.15, TLI=.472

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.248, RMSR=.15, TLI=.46

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.061, RMSR=.02, TLI=.96
#     MR4  MR1  MR2  MR3
#MR4 1.00 0.56 0.50 0.53
#MR1 0.56 1.00 0.31 0.50
#MR2 0.50 0.31 1.00 0.32
#MR3 0.53 0.50 0.32 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.74, RMSEA=.08, RMSR=.02, TLI=.944
#     MR4  MR1  MR2  MR3
#MR4 1.00 0.57 0.55 0.57
#MR1 0.57 1.00 0.34 0.54
#MR2 0.55 0.34 1.00 0.35
#MR3 0.57 0.54 0.35 1.00


# Single factor model lavaan
UNImodel= '
 general =~ EI_1+EI_2+EI_3+EI_4+
            EI_5+EI_6+EI_7+EI_8+
            EI_9+EI_10+EI_11+EI_12+
            EI_13+EI_14+EI_15+EI_16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.918       0.782
#Tucker-Lewis Index (TLI)                       0.905       0.748
#Robust Comparative Fit Index (CFI)                         0.444
#Robust Tucker-Lewis Index (TLI)                            0.359
#RMSEA                                          0.323       0.280
#Robust RMSEA                                               0.275
#SRMR                                           0.184       0.184

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .566

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.546       0.524
#Tucker-Lewis Index (TLI)                       0.476       0.450
#Robust Comparative Fit Index (CFI)                         0.548
#Robust Tucker-Lewis Index (TLI)                            0.479
#RMSEA                                          0.222       0.186
#Robust RMSEA                                               0.220
#SRMR                                           0.140       0.140

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .335

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors (theory based) 
EGAmodel= '
 self_appraisal =~ EI_1+EI_2+EI_3+EI_4
 other_appraisal =~ EI_5+EI_6+EI_7+EI_8
 use_emotion =~ EI_9+EI_10+EI_11+EI_12
 regulation =~ EI_13+EI_14+EI_15+EI_16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.989
#Tucker-Lewis Index (TLI)                       0.998       0.986
#Robust Comparative Fit Index (CFI)                         0.957
#Robust Tucker-Lewis Index (TLI)                            0.948
#RMSEA                                          0.041       0.066
#Robust RMSEA                                               0.079
#SRMR                                           0.033       0.033

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .776

#Covariances:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self_appraisal ~~                                                       
#    other_appraisl      0.567    0.025   22.838    0.000    0.567    0.567
#    use_emotion         0.608    0.024   25.204    0.000    0.608    0.608
#    regulation          0.581    0.024   24.433    0.000    0.581    0.581
#  other_appraisal ~~                                                      
#    use_emotion         0.397    0.032   12.531    0.000    0.397    0.397
#    regulation          0.355    0.032   11.048    0.000    0.355    0.355
#  use_emotion ~~                                                          
#    regulation          0.568    0.026   21.686    0.000    0.568    0.568

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.971       0.971
#Tucker-Lewis Index (TLI)                       0.964       0.965
#Robust Comparative Fit Index (CFI)                         0.974
#Robust Tucker-Lewis Index (TLI)                            0.969
#RMSEA                                          0.058       0.047
#Robust RMSEA                                               0.054
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .734

#Covariances:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self_appraisal ~~                                                       
#    other_appraisl      0.520    0.038   13.644    0.000    0.520    0.520
#    use_emotion         0.549    0.038   14.288    0.000    0.549    0.549
#    regulation          0.569    0.037   15.513    0.000    0.569    0.569
#  other_appraisal ~~                                                      
#    use_emotion         0.348    0.046    7.490    0.000    0.348    0.348
#    regulation          0.319    0.044    7.203    0.000    0.319    0.319
#  use_emotion ~~                                                          
#    regulation          0.528    0.036   14.853    0.000    0.528    0.528

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.875
#Tucker-Lewis Index (TLI)                       0.867       0.855
#Robust Comparative Fit Index (CFI)                         0.888
#Robust Tucker-Lewis Index (TLI)                            0.871
#RMSEA                                          0.112       0.095
#Robust RMSEA                                               0.110
#SRMR                                           0.288       0.288

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .740


# Bifactor model with 2 factors
BIFmodel= '
 self_appraisal =~ EI_1+EI_2+EI_3+EI_4
 other_appraisal =~ EI_5+EI_6+EI_7+EI_8
 use_emotion =~ EI_9+EI_10+EI_11+EI_12
 regulation =~ EI_13+EI_14+EI_15+EI_16
 general =~ EI_1+EI_2+EI_3+EI_4+
            EI_5+EI_6+EI_7+EI_8+
            EI_9+EI_10+EI_11+EI_12+
            EI_13+EI_14+EI_15+EI_16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.973       0.973
#Tucker-Lewis Index (TLI)                       0.963       0.963
#Robust Comparative Fit Index (CFI)                         0.976
#Robust Tucker-Lewis Index (TLI)                            0.968
#RMSEA                                          0.059       0.049
#Robust RMSEA                                               0.055
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .737

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5037498      0.8000000      0.9554907      0.7656417 
#
#$FactorLevelIndices
#                   ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#self_appraisal  0.2735603 0.07166935 0.7264397 0.9097962 0.2374149 0.5094227 0.7258496
#other_appraisal 0.6741966 0.16483257 0.3258034 0.8904076 0.5979407 0.7832788 0.8972075
#use_emotion     0.5082849 0.12073540 0.4917151 0.8816288 0.4444561 0.6795387 0.8383348
#regulation      0.5430397 0.13901292 0.4569603 0.9036319 0.4903253 0.7254706 0.8778725
#general         0.5037498 0.50374977 0.5037498 0.9554907 0.7656417 0.9054628 0.8967648







################################################################
### 
### Robinson et al (2023) https://www.sciencedirect.com/science/article/pii/S0191886923002246
### Wong and Law’s (2002) 16-item Emotional Intelligence Scale measuring 4 dimensions of emotional intelligence
### data sent by the author

library(readxl)
EI_Robinson <- read_excel("EI_Robinson.xlsx")
colnames(EI_Robinson)
mydata <- EI_Robinson[,2:17]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .92, omega T = .94

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 6.6; eigenvalue 2 = 1.48

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 7.0; eigenvalue 2 = 1.6

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.41, RMSEA=.194, RMSR=.14, TLI=.533

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.227, RMSR=.15, TLI=.473

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 0 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.67, RMSEA=.02, RMSR=.02, TLI=.994
#     MR1  MR4  MR2  MR3
#MR1 1.00 0.54 0.47 0.49
#MR4 0.54 1.00 0.32 0.41
#MR2 0.47 0.32 1.00 0.46
#MR3 0.49 0.41 0.46 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.72, RMSEA=.077, RMSR=.02, TLI=.938
#     MR4  MR2  MR1  MR3
#MR4 1.00 0.51 0.51 0.49
#MR2 0.51 1.00 0.32 0.37
#MR1 0.51 0.32 1.00 0.48
#MR3 0.49 0.37 0.48 1.00


# Single factor model lavaan
UNImodel= '
 general =~ EI_1+EI_2+EI_3+EI_4+
            EI_5+EI_6+EI_7+EI_8+
            EI_9+EI_10+EI_11+EI_12+
            EI_13+EI_14+EI_15+EI_16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.909       0.751
#Tucker-Lewis Index (TLI)                       0.895       0.713
#Robust Comparative Fit Index (CFI)                         0.526
#Robust Tucker-Lewis Index (TLI)                            0.453
#RMSEA                                          0.281       0.254
#Robust RMSEA                                               0.243
#SRMR                                           0.161       0.161

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .518

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.596       0.587
#Tucker-Lewis Index (TLI)                       0.533       0.523
#Robust Comparative Fit Index (CFI)                         0.600
#Robust Tucker-Lewis Index (TLI)                            0.539
#RMSEA                                          0.200       0.183
#Robust RMSEA                                               0.197
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .408

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors (theory based) 
EGAmodel= '
 self_appraisal =~ EI_1+EI_2+EI_3+EI_4
 other_appraisal =~ EI_5+EI_6+EI_7+EI_8
 use_emotion =~ EI_9+EI_10+EI_11+EI_12
 regulation =~ EI_13+EI_14+EI_15+EI_16
'
library(lavaan)
#CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.982       0.988
#Tucker-Lewis Index (TLI)                       0.978       0.985
#Robust Comparative Fit Index (CFI)                         0.989
#Robust Tucker-Lewis Index (TLI)                            0.986
#RMSEA                                          0.043       0.032
#Robust RMSEA                                               0.034
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .704

#Covariances:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self_appraisal ~~                                                       
#    other_appraisl      0.600    0.080    7.523    0.000    0.600    0.600
#    use_emotion         0.578    0.091    6.358    0.000    0.578    0.578
#    regulation          0.567    0.088    6.467    0.000    0.567    0.567
#  other_appraisal ~~                                                      
#    use_emotion         0.464    0.103    4.530    0.000    0.464    0.464
#    regulation          0.395    0.105    3.762    0.000    0.395    0.395
#  use_emotion ~~                                                          
#    regulation          0.559    0.086    6.471    0.000    0.559    0.559

EGAmodel= '
 self_appraisal =~ EI_1+EI_2+EI_3+EI_4
 other_appraisal =~ EI_5+EI_6+EI_7+EI_8
 use_emotion =~ EI_9+EI_10+EI_11+EI_12
 regulation =~ EI_13+EI_14+EI_15+EI_16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.984
#Tucker-Lewis Index (TLI)                       1.001       0.981
#Robust Comparative Fit Index (CFI)                         0.954
#Robust Tucker-Lewis Index (TLI)                            0.944
#RMSEA                                          0.000       0.065
#Robust RMSEA                                               0.078
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .730

#Covariances:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self_appraisal ~~                                                       
#    other_appraisl      0.588    0.051   11.421    0.000    0.588    0.588
#    use_emotion         0.584    0.055   10.531    0.000    0.584    0.584
#    regulation          0.600    0.050   11.900    0.000    0.600    0.600
#  other_appraisal ~~                                                      
#    use_emotion         0.420    0.058    7.229    0.000    0.420    0.420
#    regulation          0.374    0.059    6.311    0.000    0.374    0.374
#  use_emotion ~~                                                          
#    regulation          0.582    0.053   10.891    0.000    0.582    0.582

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.982       0.988
#Tucker-Lewis Index (TLI)                       0.978       0.985
#Robust Comparative Fit Index (CFI)                         0.989
#Robust Tucker-Lewis Index (TLI)                            0.986
#RMSEA                                          0.043       0.032
#Robust RMSEA                                               0.034
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .704

#Covariances:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self_appraisal ~~                                                       
#    other_appraisl      0.600    0.080    7.523    0.000    0.600    0.600
#    use_emotion         0.578    0.091    6.358    0.000    0.578    0.578
#    regulation          0.567    0.088    6.467    0.000    0.567    0.567
#  other_appraisal ~~                                                      
#    use_emotion         0.464    0.103    4.530    0.000    0.464    0.464
#    regulation          0.395    0.105    3.762    0.000    0.395    0.395
#  use_emotion ~~                                                          
#    regulation          0.559    0.086    6.471    0.000    0.559    0.559


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.877       0.875
#Tucker-Lewis Index (TLI)                       0.858       0.856
#Robust Comparative Fit Index (CFI)                         0.883
#Robust Tucker-Lewis Index (TLI)                            0.865
#RMSEA                                          0.110       0.101
#Robust RMSEA                                               0.107
#SRMR                                           0.300       0.300

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .707


# Bifactor model with 2 factors
BIFmodel= '
 self_appraisal =~ EI_1+EI_2+EI_3+EI_4
 other_appraisal =~ EI_5+EI_6+EI_7+EI_8
 use_emotion =~ EI_9+EI_10+EI_11+EI_12
 regulation =~ EI_13+EI_14+EI_15+EI_16
 general =~ EI_1+EI_2+EI_3+EI_4+
            EI_5+EI_6+EI_7+EI_8+
            EI_9+EI_10+EI_11+EI_12+
            EI_13+EI_14+EI_15+EI_16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.988       0.994
#Tucker-Lewis Index (TLI)                       0.984       0.991
#Robust Comparative Fit Index (CFI)                         0.994
#Robust Tucker-Lewis Index (TLI)                            0.992
#RMSEA                                          0.037       0.025
#Robust RMSEA                                               0.026
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .707

#Latent Variables:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self_appraisal =~                                                       
#    EI_1                0.799    0.169    4.725    0.000    0.799    0.533
#    EI_2                0.813    0.143    5.695    0.000    0.813    0.503
#    EI_3                0.681    0.185    3.691    0.000    0.681    0.441
#    EI_4                0.272    0.290    0.938    0.348    0.272    0.163
#  other_appraisal =~                                                      
#    EI_5                0.785    0.089    8.788    0.000    0.785    0.616
#    EI_6                0.972    0.102    9.538    0.000    0.972    0.677
#    EI_7                0.625    0.142    4.414    0.000    0.625    0.381
#    EI_8                0.912    0.089   10.302    0.000    0.912    0.657
#  use_emotion =~                                                          
#    EI_9                0.717    0.157    4.560    0.000    0.717    0.446
#    EI_10               0.691    0.163    4.236    0.000    0.691    0.450
#    EI_11               0.998    0.134    7.467    0.000    0.998    0.670
#    EI_12               0.844    0.138    6.095    0.000    0.844    0.598
#  regulation =~                                                           
#    EI_13               0.933    0.142    6.565    0.000    0.933    0.586
#    EI_14               0.958    0.128    7.509    0.000    0.958    0.635
#    EI_15               1.093    0.123    8.862    0.000    1.093    0.690
#    EI_16               0.950    0.133    7.134    0.000    0.950    0.649
#  general =~                                                              
#    EI_1                1.018    0.125    8.157    0.000    1.018    0.679
#    EI_2                1.214    0.116   10.451    0.000    1.214    0.751
#    EI_3                1.165    0.121    9.609    0.000    1.165    0.754
#    EI_4                1.177    0.181    6.500    0.000    1.177    0.704
#    EI_5                0.654    0.130    5.029    0.000    0.654    0.512
#    EI_6                0.873    0.161    5.418    0.000    0.873    0.608
#    EI_7                0.788    0.178    4.421    0.000    0.788    0.481
#    EI_8                0.865    0.153    5.659    0.000    0.865    0.623
#    EI_9                0.913    0.169    5.413    0.000    0.913    0.568
#    EI_10               0.778    0.151    5.153    0.000    0.778    0.506
#    EI_11               0.810    0.190    4.263    0.000    0.810    0.544
#    EI_12               0.836    0.197    4.249    0.000    0.836    0.593
#    EI_13               0.956    0.171    5.593    0.000    0.956    0.601
#    EI_14               0.968    0.151    6.402    0.000    0.968    0.642
#    EI_15               0.475    0.151    3.155    0.002    0.475    0.300
#    EI_16               0.764    0.144    5.312    0.000    0.764    0.522

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5316985      0.8000000      0.9532308      0.7815975 
#
#$FactorLevelIndices
#                   ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#self_appraisal  0.2661405 0.07063841 0.7338595 0.9053384 0.2206549 0.5008353 0.7308467
#other_appraisal 0.5303783 0.13180716 0.4696217 0.8860752 0.4636307 0.7045380 0.8689999
#use_emotion     0.4962705 0.11260877 0.5037295 0.8593699 0.4204936 0.6519219 0.8324338
#regulation      0.5913298 0.15324713 0.4086702 0.8986245 0.5443956 0.7391923 0.8725267
#general         0.5316985 0.53169853 0.5316985 0.9532308 0.7815975 0.9078463 0.9097970





################################################################
### 
### Anglim et al (2020) https://onlinelibrary.wiley.com/doi/full/10.1111/jopy.12493
### Assessing Emotions Scale (AES) developed by Schutte
### data https://osf.io/5wpj6

AES_Anglim <- read.csv("AES_Anglim.csv")
colnames(AES_Anglim)
mydata <- AES_Anglim[,254:286]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .87, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 5 components
# Eigenvalue 1 = 6.34; eigenvalue 2 = 1.87

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 5 components
# Eigenvalue 1 = 8.24; eigenvalue 2 = 2.16

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.082, RMSR=.08, TLI=.575

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.104, RMSR=.1, TLI=.552

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.045, RMSR=.03, TLI=.873
#      MR1   MR2   MR3   MR4
#MR1  1.00  0.32 -0.39 -0.04
#MR2  0.32  1.00 -0.12  0.00
#MR3 -0.39 -0.12  1.00 -0.09
#MR4 -0.04  0.00 -0.09  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.066, RMSR=.04, TLI=.821
#      MR1   MR2   MR3   MR4
#MR1  1.00  0.37 -0.40 -0.04
#MR2  0.37  1.00 -0.14 -0.02
#MR3 -0.40 -0.14  1.00 -0.10
#MR4 -0.04 -0.02 -0.10  1.00


# Single factor model lavaan
UNImodel= '
 general =~ aes11+aes13+aes16+aes26+aes1+aes4+aes24+aes30+aes10+aes28+aes3+aes21+
            aes31+aes2+aes12+aes14+aes23+aes5+aes18+aes19+aes32+aes33+aes9+aes15+
            aes22+aes25+aes29+aes6+aes7+aes8+aes17+aes20+aes27
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.887       0.741
#Tucker-Lewis Index (TLI)                       0.880       0.723
#Robust Comparative Fit Index (CFI)                         0.587
#Robust Tucker-Lewis Index (TLI)                            0.559
#RMSEA                                          0.108       0.101
#Robust RMSEA                                               0.104
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .254

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.603       0.597
#Tucker-Lewis Index (TLI)                       0.577       0.570
#Robust Comparative Fit Index (CFI)                         0.607
#Robust Tucker-Lewis Index (TLI)                            0.581
#RMSEA                                          0.082       0.075
#Robust RMSEA                                               0.082
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .186

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors (theory based) 
EGAmodel= '
 managing_others =~ aes11+aes13+aes16+aes26+aes1+aes4+aes24+aes30
 managing_own =~ aes10+aes28+aes3+aes21+aes31+aes2+aes12+aes14+aes23
 perception_emotion =~ aes5+aes18+aes19+aes32+aes33+aes9+aes15+aes22+aes25+aes29
 use_emotion =~ aes6+aes7+aes8+aes17+aes20+aes27
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.920       0.810
#Tucker-Lewis Index (TLI)                       0.914       0.795
#Robust Comparative Fit Index (CFI)                         0.697
#Robust Tucker-Lewis Index (TLI)                            0.672
#RMSEA                                          0.091       0.087
#Robust RMSEA                                               0.089
#SRMR                                           0.084       0.084

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .354

#Covariances:
#                        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  managing_others ~~                                                         
#    managing_own           0.760    0.020   37.106    0.000    0.760    0.760
#    perception_mtn         0.757    0.018   42.498    0.000    0.757    0.757
#    use_emotion           -0.611    0.028  -22.153    0.000   -0.611   -0.611
#  managing_own ~~                                                            
#    perception_mtn         0.667    0.018   36.217    0.000    0.667    0.667
#    use_emotion           -0.601    0.026  -23.376    0.000   -0.601   -0.601
#  perception_emotion ~~                                                      
#    use_emotion           -0.402    0.029  -13.870    0.000   -0.402   -0.402

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.723       0.718
#Tucker-Lewis Index (TLI)                       0.701       0.696
#Robust Comparative Fit Index (CFI)                         0.728
#Robust Tucker-Lewis Index (TLI)                            0.706
#RMSEA                                          0.069       0.063
#Robust RMSEA                                               0.068
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .257

#Covariances:
#                        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  managing_others ~~                                                         
#    managing_own           0.707    0.043   16.493    0.000    0.707    0.707
#    perception_mtn         0.700    0.028   24.946    0.000    0.700    0.700
#    use_emotion           -0.548    0.043  -12.685    0.000   -0.548   -0.548
#  managing_own ~~                                                            
#    perception_mtn         0.617    0.030   20.332    0.000    0.617    0.617
#    use_emotion           -0.466    0.071   -6.599    0.000   -0.466   -0.466
#  perception_emotion ~~                                                      
#    use_emotion           -0.315    0.046   -6.886    0.000   -0.315   -0.315

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (4 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.624       0.618
#Tucker-Lewis Index (TLI)                       0.599       0.593
#Robust Comparative Fit Index (CFI)                         0.628
#Robust Tucker-Lewis Index (TLI)                            0.603
#RMSEA                                          0.080       0.073
#Robust RMSEA                                               0.079
#SRMR                                           0.152       0.152

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .261


# Bifactor model with 4 factors
BIFmodel= '
 managing_others =~ aes11+aes13+aes16+aes26+aes1+aes4+aes24+aes30
 managing_own =~ aes10+aes28+aes3+aes21+aes31+aes2+aes12+aes14+aes23
 perception_emotion =~ aes5+aes18+aes19+aes32+aes33+aes9+aes15+aes22+aes25+aes29
 use_emotion =~ aes6+aes7+aes8+aes17+aes20+aes27
 general =~ aes11+aes13+aes16+aes26+aes1+aes4+aes24+aes30+aes10+aes28+aes3+aes21+
            aes31+aes2+aes12+aes14+aes23+aes5+aes18+aes19+aes32+aes33+aes9+aes15+
            aes22+aes25+aes29+aes6+aes7+aes8+aes17+aes20+aes27
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.807       0.803
#Tucker-Lewis Index (TLI)                       0.779       0.775
#Robust Comparative Fit Index (CFI)                         0.811
#Robust Tucker-Lewis Index (TLI)                            0.785
#RMSEA                                          0.060       0.054
#Robust RMSEA                                               0.059
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .314

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5216605      0.7651515      0.5676614      0.4183101 
#
#$FactorLevelIndices
#                      ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#managing_others    0.4003074 0.06660613 0.5996926 0.4988113 0.13161915 0.5001451 0.7300909
#managing_own       0.3447307 0.10285938 0.6552693 0.5023372 0.38820849 0.5843539 0.7833954
#perception_emotion 0.5436767 0.21110938 0.4563233 0.1126184 0.06747994 0.7673784 0.8691725
#use_emotion        0.6653503 0.09776462 0.3346497 0.2848187 0.17714407 0.5942615 0.7753700
#general            0.5216605 0.52166049 0.5216605 0.5676614 0.41831012 0.8734490 0.9159839







################################################################
### 
### AES dataset sent by Demetrovics based on study https://psycnet.apa.org/record/2022-01619-001
### for 3-factor model, see http://www.demetrovics.hu/dokumentumok/Kun_etal_2010_AES_BRM.pdf
### Assessing Emotions Scale (AES) developed by Schutte
### 

AES_Demetrovics <- read_sav("AES_Demetrovics.sav")
colnames(AES_Demetrovics)
mydata <- AES_Demetrovics[,38:70]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .92, omega T = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 5 components
# Eigenvalue 1 = 8.94; eigenvalue 2 = 1.44

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 5 components
# Eigenvalue 1 = 10.49; eigenvalue 2 = 1.62

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.27, RMSEA=.088, RMSR=.08, TLI=.653

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.108, RMSR=.08, TLI=.610

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 4 factors (theory based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.059, RMSR=.04, TLI=.844
#     MR1  MR2  MR4  MR3
#MR1 1.00 0.44 0.37 0.26
#MR2 0.44 1.00 0.33 0.31
#MR4 0.37 0.33 1.00 0.34
#MR3 0.26 0.31 0.34 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.077, RMSR=.04, TLI=.8
#     MR1  MR2  MR4  MR3
#MR1 1.00 0.46 0.35 0.27
#MR2 0.46 1.00 0.30 0.31
#MR4 0.35 0.30 1.00 0.33
#MR3 0.27 0.31 0.33 1.00


# Single factor model lavaan
UNImodel= '
 general =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30+AES10+AES28+AES3+AES21+
            AES31+AES2+AES12+AES14+AES23+AES5+AES18+AES19+AES32+AES33+AES9+AES15+
            AES22+AES25+AES29+AES6+AES7+AES8+AES17+AES20+AES27
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.933       0.765
#Tucker-Lewis Index (TLI)                       0.928       0.749
#Robust Comparative Fit Index (CFI)                         0.636
#Robust Tucker-Lewis Index (TLI)                            0.612
#RMSEA                                          0.105       0.106
#Robust RMSEA                                               0.109
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .319

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.675       0.675
#Tucker-Lewis Index (TLI)                       0.653       0.653
#Robust Comparative Fit Index (CFI)                         0.678
#Robust Tucker-Lewis Index (TLI)                            0.657
#RMSEA                                          0.089       0.079
#Robust RMSEA                                               0.088
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .264

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors (theory based) 
EGAmodel= '
 managing_others =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30
 managing_own =~ AES10+AES28+AES3+AES21+AES31+AES2+AES12+AES14+AES23
 perception_emotion =~ AES5+AES18+AES19+AES32+AES33+AES9+AES15+AES22+AES25+AES29
 use_emotion =~ AES6+AES7+AES8+AES17+AES20+AES27
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.959       0.847
#Tucker-Lewis Index (TLI)                       0.956       0.834
#Robust Comparative Fit Index (CFI)                         0.736
#Robust Tucker-Lewis Index (TLI)                            0.715
#RMSEA                                          0.082       0.086
#Robust RMSEA                                               0.093
#SRMR                                           0.070       0.070

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .398

#Covariances:
#                        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  managing_others ~~                                                         
#    managing_own           0.807    0.013   60.533    0.000    0.807    0.807
#    perception_mtn        -0.733    0.014  -54.125    0.000   -0.733   -0.733
#    use_emotion            0.758    0.017   44.910    0.000    0.758    0.758
#  managing_own ~~                                                            
#    perception_mtn        -0.678    0.015  -44.247    0.000   -0.678   -0.678
#    use_emotion            0.729    0.016   45.336    0.000    0.729    0.729
#  perception_emotion ~~                                                      
#    use_emotion           -0.618    0.017  -36.206    0.000   -0.618   -0.618

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.779       0.778
#Tucker-Lewis Index (TLI)                       0.761       0.761
#Robust Comparative Fit Index (CFI)                         0.783
#Robust Tucker-Lewis Index (TLI)                            0.765
#RMSEA                                          0.073       0.065
#Robust RMSEA                                               0.073
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .333

#Covariances:
#                        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  managing_others ~~                                                         
#    managing_own           0.781    0.020   38.153    0.000    0.781    0.781
#    perception_mtn        -0.724    0.020  -35.338    0.000   -0.724   -0.724
#    use_emotion            0.758    0.024   31.000    0.000    0.758    0.758
#  managing_own ~~                                                            
#    perception_mtn        -0.644    0.026  -25.124    0.000   -0.644   -0.644
#    use_emotion            0.739    0.024   31.050    0.000    0.739    0.739
#  perception_emotion ~~                                                      
#    use_emotion           -0.634    0.026  -24.822    0.000   -0.634   -0.634

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (4 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.648       0.644
#Tucker-Lewis Index (TLI)                       0.624       0.620
#Robust Comparative Fit Index (CFI)                         0.651
#Robust Tucker-Lewis Index (TLI)                            0.628
#RMSEA                                          0.092       0.082
#Robust RMSEA                                               0.091
#SRMR                                           0.215       0.215

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .338


# Bifactor model with 4 factors
BIFmodel= '
 managing_others =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30
 managing_own =~ AES10+AES28+AES3+AES21+AES31+AES2+AES12+AES14+AES23
 perception_emotion =~ AES5+AES18+AES19+AES32+AES33+AES9+AES15+AES22+AES25+AES29
 use_emotion =~ AES6+AES7+AES8+AES17+AES20+AES27
 general =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30+AES10+AES28+AES3+AES21+
            AES31+AES2+AES12+AES14+AES23+AES5+AES18+AES19+AES32+AES33+AES9+AES15+
            AES22+AES25+AES29+AES6+AES7+AES8+AES17+AES20+AES27
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.824       0.820
#Tucker-Lewis Index (TLI)                       0.798       0.794
#Robust Comparative Fit Index (CFI)                         0.827
#Robust Tucker-Lewis Index (TLI)                            0.802
#RMSEA                                          0.067       0.061
#Robust RMSEA                                               0.067
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .351

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6385550      0.7651515      0.9088545      0.8177412 
#
#$FactorLevelIndices
#                      ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#managing_others    0.3753338 0.08103646 0.6246662 0.7867649 0.2502284 0.5685095 0.7693126
#managing_own       0.2706118 0.06750248 0.7293882 0.7655041 0.2094718 0.4923340 0.7006781
#perception_emotion 0.3594860 0.13404522 0.6405140 0.7802080 0.1666733 0.6984852 0.8432680
#use_emotion        0.4874851 0.07886082 0.5125149 0.7218939 0.3473547 0.5864751 0.7696822
#general            0.6385550 0.63855502 0.6385550 0.9088545 0.8177412 0.9207538 0.9402321









################################################################
### 
### AES dataset sent by Demetrovics 2 based on study https://www.frontiersin.org/journals/psychiatry/articles/10.3389/fpsyt.2022.831992/full
### for 3-factor model, see http://www.demetrovics.hu/dokumentumok/Kun_etal_2010_AES_BRM.pdf
### Assessing Emotions Scale (AES) developed by Schutte

AES_Demetrovics2 <- read_sav("AES_Demetrovics2.sav")
colnames(AES_Demetrovics2)
mydata <- AES_Demetrovics2[,44:76]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .92, omega T = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 7.43; eigenvalue 2 = 1.52

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 6 components
# Eigenvalue 1 = 8.66; eigenvalue 2 = 1.71

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.23, RMSEA=.089, RMSR=.08, TLI=.58

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.113, RMSR=.09, TLI=.518

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 4 factors (theory based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.065, RMSR=.05, TLI=.778
#     MR1  MR3  MR2  MR4
#MR1 1.00 0.47 0.20 0.21
#MR3 0.47 1.00 0.29 0.19
#MR2 0.20 0.29 1.00 0.24
#MR4 0.21 0.19 0.24 1.0

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.089, RMSR=.05, TLI=.698
#     MR1  MR3  MR2  MR4
#MR1 1.00 0.43 0.24 0.21
#MR3 0.43 1.00 0.32 0.25
#MR2 0.24 0.32 1.00 0.30
#MR4 0.21 0.25 0.30 1.00


# Single factor model lavaan
UNImodel= '
 general =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30+AES10+AES28+AES3+AES21+
            AES31+AES2+AES12+AES14+AES23+AES5+AES18+AES19+AES32+AES33+AES9+AES15+
            AES22+AES25+AES29+AES6+AES7+AES8+AES17+AES20+AES27
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.898       0.739
#Tucker-Lewis Index (TLI)                       0.891       0.722
#Robust Comparative Fit Index (CFI)                         0.553
#Robust Tucker-Lewis Index (TLI)                            0.524
#RMSEA                                          0.111       0.102
#Robust RMSEA                                               0.114
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .252

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.606       0.608
#Tucker-Lewis Index (TLI)                       0.580       0.582
#Robust Comparative Fit Index (CFI)                         0.615
#Robust Tucker-Lewis Index (TLI)                            0.589
#RMSEA                                          0.091       0.081
#Robust RMSEA                                               0.089
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .201

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors (theory based) 
EGAmodel= '
 managing_others =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30
 managing_own =~ AES10+AES28+AES3+AES21+AES31+AES2+AES12+AES14+AES23
 perception_emotion =~ AES5+AES18+AES19+AES32+AES33+AES9+AES15+AES22+AES25+AES29
 use_emotion =~ AES6+AES7+AES8+AES17+AES20+AES27
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.930       0.813
#Tucker-Lewis Index (TLI)                       0.925       0.798
#Robust Comparative Fit Index (CFI)                         0.648
#Robust Tucker-Lewis Index (TLI)                            0.620
#RMSEA                                          0.093       0.087
#Robust RMSEA                                               0.101
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .335

#Covariances:
#                        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  managing_others ~~                                                         
#    managing_own           0.830    0.022   37.397    0.000    0.830    0.830
#    perception_mtn        -0.673    0.028  -24.076    0.000   -0.673   -0.673
#    use_emotion            0.712    0.033   21.464    0.000    0.712    0.712
#  managing_own ~~                                                            
#    perception_mtn        -0.701    0.023  -30.976    0.000   -0.701   -0.701
#    use_emotion            0.602    0.032   18.785    0.000    0.602    0.602
#  perception_emotion ~~                                                      
#    use_emotion           -0.544    0.032  -17.091    0.000   -0.544   -0.544

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.708       0.711
#Tucker-Lewis Index (TLI)                       0.685       0.688
#Robust Comparative Fit Index (CFI)                         0.718
#Robust Tucker-Lewis Index (TLI)                            0.695
#RMSEA                                          0.078       0.070
#Robust RMSEA                                               0.076
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .270

#Covariances:
#                        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  managing_others ~~                                                         
#    managing_own           0.788    0.043   18.289    0.000    0.788    0.788
#    perception_mtn        -0.642    0.049  -13.138    0.000   -0.642   -0.642
#    use_emotion            0.687    0.054   12.666    0.000    0.687    0.687
#  managing_own ~~                                                            
#    perception_mtn        -0.669    0.038  -17.521    0.000   -0.669   -0.669
#    use_emotion            0.588    0.055   10.767    0.000    0.588    0.588
#  perception_emotion ~~                                                      
#    use_emotion           -0.541    0.050  -10.815    0.000   -0.541   -0.541

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (4 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.595       0.595
#Tucker-Lewis Index (TLI)                       0.568       0.567
#Robust Comparative Fit Index (CFI)                         0.603
#Robust Tucker-Lewis Index (TLI)                            0.577
#RMSEA                                          0.092       0.082
#Robust RMSEA                                               0.090
#SRMR                                           0.181       0.181

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .284


# Bifactor model with 4 factors
BIFmodel= '
 managing_others =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30
 managing_own =~ AES10+AES28+AES3+AES21+AES31+AES2+AES12+AES14+AES23
 perception_emotion =~ AES5+AES18+AES19+AES32+AES33+AES9+AES15+AES22+AES25+AES29
 use_emotion =~ AES6+AES7+AES8+AES17+AES20+AES27
 general =~ AES11+AES13+AES16+AES26+AES1+AES4+AES24+AES30+AES10+AES28+AES3+AES21+
            AES31+AES2+AES12+AES14+AES23+AES5+AES18+AES19+AES32+AES33+AES9+AES15+
            AES22+AES25+AES29+AES6+AES7+AES8+AES17+AES20+AES27
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.775       0.774
#Tucker-Lewis Index (TLI)                       0.743       0.741
#Robust Comparative Fit Index (CFI)                         0.783
#Robust Tucker-Lewis Index (TLI)                            0.752
#RMSEA                                          0.071       0.063
#Robust RMSEA                                               0.069
#SRMR                                           0.063       0.063

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .318

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5940778      0.7651515      0.8847311      0.7899555 
#
#$FactorLevelIndices
#                      ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#managing_others    0.4114807 0.08605346 0.5885193 0.7415207 0.245056573 0.5866753 0.7820240
#managing_own       0.2018594 0.05113649 0.7981406 0.6959659 0.007647049 0.4073515 0.6694588
#perception_emotion 0.4336244 0.15903122 0.5663756 0.7515763 0.249236221 0.7197481 0.8515282
#use_emotion        0.6423050 0.10970099 0.3576950 0.6996469 0.477313459 0.6619797 0.8180958
#general            0.5940778 0.59407784 0.5940778 0.8847311 0.789955538 0.9012335 0.9344981













