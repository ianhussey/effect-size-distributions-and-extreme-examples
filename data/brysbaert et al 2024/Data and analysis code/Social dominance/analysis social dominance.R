################################################################
### Social dominance orientation scale
### 
### 16 items: either unidimensional or bidimensional (dominance vs. equality)
### Could be due to acquescence bias



################################################################
### 
### Roets (2019) https://www.sciencedirect.com/science/article/pii/S0191886919300777
### data study 1: https://osf.io/d4vym/?view_only=8da12c7350174552bce4b10c3ee3f246

library(haven)
SDO_Roets1 <- read_sav("SDO_Roets1.sav")
colnames(SDO_Roets1)
mydata <- SDO_Roets1[,15:30]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .96, omega T = .98

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 9.23; eigenvalue 2 = 2.57

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 11.03; eigenvalue 2 = 1.99

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.58, RMSEA=.25, RMSR=.17, TLI=.555

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.265, RMSR=.13, TLI=.618

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.76, RMSEA=.084, RMSR=.02, TLI=.949
#      MR1   MR2
#MR1  1.00   0.52
#MR2  0.52  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.83, RMSEA=.1, RMSR=.02, TLI=.945
#      MR1   MR2
#MR1  1.00   0.65
#MR2  0.65   1.00


# Single factor model lavaan
UNImodel= '
 general =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08+
            sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.940
#Tucker-Lewis Index (TLI)                       0.973       0.931
#Robust Comparative Fit Index (CFI)                         0.657
#Robust Tucker-Lewis Index (TLI)                            0.604
#RMSEA                                          0.296       0.231
#Robust RMSEA                                               0.277
#SRMR                                           0.145       0.145

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .775

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.613       0.526
#Tucker-Lewis Index (TLI)                       0.554       0.453
#Robust Comparative Fit Index (CFI)                         0.618
#Robust Tucker-Lewis Index (TLI)                            0.559
#RMSEA                                          0.254       0.193
#Robust RMSEA                                               0.249
#SRMR                                           0.178       0.178

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .546

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based) 
EGAmodel= '
 factor1 =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08
 factor2 =~ sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.989
#Tucker-Lewis Index (TLI)                       0.999       0.987
#Robust Comparative Fit Index (CFI)                         0.948
#Robust Tucker-Lewis Index (TLI)                            0.939
#RMSEA                                          0.065       0.100
#Robust RMSEA                                               0.108
#SRMR                                           0.042       0.042

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .84

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.691    0.026   26.241    0.000    0.691    0.691

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.951       0.971
#Tucker-Lewis Index (TLI)                       0.943       0.966
#Robust Comparative Fit Index (CFI)                         0.973
#Robust Tucker-Lewis Index (TLI)                            0.968
#RMSEA                                          0.091       0.048
#Robust RMSEA                                               0.067
#SRMR                                           0.049       0.049

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .769

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.553    0.061    9.098    0.000    0.553    0.553

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.932       0.951
#Tucker-Lewis Index (TLI)                       0.922       0.943
#Robust Comparative Fit Index (CFI)                         0.954
#Robust Tucker-Lewis Index (TLI)                            0.946
#RMSEA                                          0.106       0.062
#Robust RMSEA                                               0.087
#SRMR                                           0.291       0.291

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .766


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08
 factor2 =~ sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
 general =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08+
            sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.969       0.980
#Tucker-Lewis Index (TLI)                       0.958       0.973
#Robust Comparative Fit Index (CFI)                         0.983
#Robust Tucker-Lewis Index (TLI)                            0.977
#RMSEA                                          0.078       0.043
#Robust RMSEA                                               0.056
#SRMR                                           0.030       0.030

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .786

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
##$ModelLevelIndices
##ECV.general            PUC  Omega.general OmegaH.general 
##     0.6887961      0.5333333      0.9802025      0.7808144 
##$FactorLevelIndices
##        ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
##factor1 0.1782705 0.08947396 0.8217295 0.9678800 0.1356588 0.6089511 0.7939723
##factor2 0.4451515 0.22172997 0.5548485 0.9671586 0.3905000 0.8342927 0.9266390
##general 0.6887961 0.68879607 0.6887961 0.9802025 0.7808144 0.9657571 0.9563073





################################################################
### 
### Roets (2019) https://www.sciencedirect.com/science/article/pii/S0191886919300777
### data study 2: https://osf.io/d4vym/?view_only=8da12c7350174552bce4b10c3ee3f246

library(haven)
SDO_Roets2 <- read_sav("SDO_Roets2.sav")
colnames(SDO_Roets2)
mydata <- SDO_Roets2[,16:31]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .96, omega T = .97

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 9.56; eigenvalue 2 = 1.65

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 11.09; eigenvalue 2 = 1.37

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.60, RMSEA=.201, RMSR=.11, TLI=.687

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.218, RMSR=.1, TLI=.717

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.72, RMSEA=.091, RMSR=.03, TLI=.936
#      MR1   MR2
#MR1  1.00   0.65
#MR2  0.65   1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.79, RMSEA=.105, RMSR=.02, TLI=.934
#      MR1   MR2
#MR1  1.00   0.73
#MR2  0.73   1.00


# Single factor model lavaan
UNImodel= '
 general =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08+
            sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.987       0.941
#Tucker-Lewis Index (TLI)                       0.985       0.931
#Robust Comparative Fit Index (CFI)                         0.740
#Robust Tucker-Lewis Index (TLI)                            0.700
#RMSEA                                          0.205       0.209
#Robust RMSEA                                               0.228
#SRMR                                           0.101       0.101

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .747

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.729       0.722
#Tucker-Lewis Index (TLI)                       0.687       0.679
#Robust Comparative Fit Index (CFI)                         0.741
#Robust Tucker-Lewis Index (TLI)                            0.701
#RMSEA                                          0.203       0.135
#Robust RMSEA                                               0.196
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .627

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based) 
EGAmodel= '
 factor1 =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08
 factor2 =~ sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.986
#Tucker-Lewis Index (TLI)                       0.998       0.984
#Robust Comparative Fit Index (CFI)                         0.940
#Robust Tucker-Lewis Index (TLI)                            0.930
#RMSEA                                          0.067       0.102
#Robust RMSEA                                               0.111
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .806

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.777    0.017   44.641    0.000    0.777    0.777

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.938       0.950
#Tucker-Lewis Index (TLI)                       0.928       0.942
#Robust Comparative Fit Index (CFI)                         0.955
#Robust Tucker-Lewis Index (TLI)                            0.947
#RMSEA                                          0.097       0.057
#Robust RMSEA                                               0.082
#SRMR                                           0.051       0.051

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .729

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.688    0.036   18.976    0.000    0.688    0.688

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.903       0.910
#Tucker-Lewis Index (TLI)                       0.888       0.896
#Robust Comparative Fit Index (CFI)                         0.918
#Robust Tucker-Lewis Index (TLI)                            0.906
#RMSEA                                          0.122       0.077
#Robust RMSEA                                               0.110
#SRMR                                           0.344       0.344

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .723


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08
 factor2 =~ sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
 general =~ sdo01+sdo02+sdo03+sdo04+sdo05+sdo06+sdo07+sdo08+
            sdo09+sdo10+sdo11+sdo12+sdo13+sdo14+sdo15+sdo16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.961       0.972
#Tucker-Lewis Index (TLI)                       0.947       0.962
#Robust Comparative Fit Index (CFI)                         0.975
#Robust Tucker-Lewis Index (TLI)                            0.966
#RMSEA                                          0.083       0.047
#Robust RMSEA                                               0.066
#SRMR                                           0.040       0.040

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .742

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6878806      0.5333333      0.9746856      0.8106972 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1267686 0.06677108 0.8732314 0.9653293 0.05183868 0.6742776 1.0519039
#factor2 0.5183959 0.24534833 0.4816041 0.9511423 0.49208406 0.8321904 0.9445828
#general 0.6878806 0.68788059 0.6878806 0.9746856 0.81069718 0.9637267 0.9770829






################################################################
### 
### Perkins et al (2020) https://journals.sagepub.com/doi/full/10.1177/1368430220951621
### data https://osf.io/g3zyj

SDO_Perkins <- read.csv("SDO_Perkins.csv")
colnames(SDO_Perkins)
mydata <- SDO_Perkins[,8:23]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .95, omega T = .96

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 8.72; eigenvalue 2 = 1.12

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 10.57; eigenvalue 2 = .98

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.54, RMSEA=.149, RMSR=.08, TLI=.782

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.66, RMSEA=.173, RMSR=.07, TLI=.791

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.092, RMSR=.04, TLI=.917
#      MR1   MR2
#MR1  1.00   -.72
#MR2  -.72   1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.73, RMSEA=.113, RMSR=.03, TLI=.91
#      MR1   MR2
#MR1  1.00   -.80
#MR2  -.80   1.00


# Single factor model lavaan
UNImodel= '
 general =~ SDO_1+SDO_2+SDO_3+SDO_4+SDO_5+SDO_6+SDO_7+SDO_8+
            SDO_9+SDO_10+SDO_11+SDO_12+SDO_13+SDO_14+SDO_15+SDO_16
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.992       0.926
#Tucker-Lewis Index (TLI)                       0.991       0.914
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.792
#RMSEA                                          0.133       0.178
#Robust RMSEA                                               0.176
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .717

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.816       0.807
#Tucker-Lewis Index (TLI)                       0.788       0.777
#Robust Comparative Fit Index (CFI)                         0.827
#Robust Tucker-Lewis Index (TLI)                            0.800
#RMSEA                                          0.149       0.115
#Robust RMSEA                                               0.143
#SRMR                                           0.084       0.084

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .576

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based; notice different order) 
EGAmodel= '
 factor1 =~ SDO_1+SDO_2+SDO_3+SDO_4+SDO_9+SDO_10+SDO_11+SDO_12
 factor2 =~ SDO_5+SDO_6+SDO_7+SDO_8+SDO_13+SDO_14+SDO_15+SDO_16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.974
#Tucker-Lewis Index (TLI)                       0.998       0.970
#Robust Comparative Fit Index (CFI)                         0.922
#Robust Tucker-Lewis Index (TLI)                            0.910
#RMSEA                                          0.055       0.106
#Robust RMSEA                                               0.116
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .781

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.829    0.018  -45.886    0.000   -0.829   -0.829

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.925       0.927
#Tucker-Lewis Index (TLI)                       0.913       0.915
#Robust Comparative Fit Index (CFI)                         0.937
#Robust Tucker-Lewis Index (TLI)                            0.927
#RMSEA                                          0.096       0.071
#Robust RMSEA                                               0.087
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .670

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.762    0.039  -19.461    0.000   -0.762   -0.762

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.870       0.863
#Tucker-Lewis Index (TLI)                       0.849       0.842
#Robust Comparative Fit Index (CFI)                         0.881
#Robust Tucker-Lewis Index (TLI)                            0.862
#RMSEA                                          0.126       0.097
#Robust RMSEA                                               0.119
#SRMR                                           0.325       0.325

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .668


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ SDO_1+SDO_2+SDO_3+SDO_4+SDO_9+SDO_10+SDO_11+SDO_12
 factor2 =~ SDO_5+SDO_6+SDO_7+SDO_8+SDO_13+SDO_14+SDO_15+SDO_16
 general =~ SDO_1+SDO_2+SDO_3+SDO_4+SDO_5+SDO_6+SDO_7+SDO_8+
            SDO_9+SDO_10+SDO_11+SDO_12+SDO_13+SDO_14+SDO_15+SDO_16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.960       0.962
#Tucker-Lewis Index (TLI)                       0.946       0.949
#Robust Comparative Fit Index (CFI)                         0.969
#Robust Tucker-Lewis Index (TLI)                            0.958
#RMSEA                                          0.076       0.055
#Robust RMSEA                                               0.066
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .697

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7564323      0.5333333      0.7920060      0.1436422 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1409997 0.08256683 0.8590003 0.9596731 0.08955395 0.5279000 0.8100808
#factor2 0.3884984 0.16100085 0.6115016 0.8999084 0.35867986 0.6809995 0.8538202
#general 0.7564323 0.75643232 0.7564323 0.7920060 0.14364223 0.9546029 0.9638669





################################################################
### 
### Simon et al (2022) https://bpspsychub.onlinelibrary.wiley.com/doi/10.1111/bjop.12587
### data experiment 1 https://osf.io/z28aq/?view_only=

library(readxl)
SDO_Simon1 <- read_excel("SDO_Simon1.xlsx")
colnames(SDO_Simon1)
mydata <- SDO_Simon1[,585:600]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .93, omega T = .96

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 7.5; eigenvalue 2 = 2.07

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 9.83; eigenvalue 2 = 1.72

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.47, RMSEA=.189, RMSR=.14, TLI=.62

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.213, RMSR=.12, TLI=.678

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.62, RMSEA=.079, RMSR=.03, TLI=.932
#      MR1   MR2
#MR1  1.00   -.51
#MR2  -.51   1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.74, RMSEA=.099, RMSR=.03, TLI=.931
#      MR1   MR2
#MR1  1.00   -.64
#MR2  -.64   1.00


# Single factor model lavaan
UNImodel= '
 general =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8+
            SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.974       0.909
#Tucker-Lewis Index (TLI)                       0.970       0.895
#Robust Comparative Fit Index (CFI)                         0.712
#Robust Tucker-Lewis Index (TLI)                            0.668
#RMSEA                                          0.192       0.187
#Robust RMSEA                                               0.218
#SRMR                                           0.123       0.123

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .695

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.675       0.629
#Tucker-Lewis Index (TLI)                       0.624       0.572
#Robust Comparative Fit Index (CFI)                         0.680
#Robust Tucker-Lewis Index (TLI)                            0.631
#RMSEA                                          0.189       0.141
#Robust RMSEA                                               0.186
#SRMR                                           0.139       0.139

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .429

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based; notice different order) 
EGAmodel= '
 factor1 =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8
 factor2 =~ SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.974
#Tucker-Lewis Index (TLI)                       0.995       0.970
#Robust Comparative Fit Index (CFI)                         0.928
#Robust Tucker-Lewis Index (TLI)                            0.917
#RMSEA                                          0.075       0.100
#Robust RMSEA                                               0.109
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .759

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.713    0.019  -36.695    0.000   -0.713   -0.713

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.929       0.929
#Tucker-Lewis Index (TLI)                       0.917       0.917
#Robust Comparative Fit Index (CFI)                         0.938
#Robust Tucker-Lewis Index (TLI)                            0.928
#RMSEA                                          0.089       0.062
#Robust RMSEA                                               0.082
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .643

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.564    0.042  -13.303    0.000   -0.564   -0.564

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.901       0.896
#Tucker-Lewis Index (TLI)                       0.885       0.880
#Robust Comparative Fit Index (CFI)                         0.910
#Robust Tucker-Lewis Index (TLI)                            0.896
#RMSEA                                          0.104       0.074
#Robust RMSEA                                               0.098
#SRMR                                           0.242       0.242

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .644


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8
 factor2 =~ SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
 general =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8+
            SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.953       0.948
#Tucker-Lewis Index (TLI)                       0.936       0.929
#Robust Comparative Fit Index (CFI)                         0.959
#Robust Tucker-Lewis Index (TLI)                            0.945
#RMSEA                                          0.078       0.057
#Robust RMSEA                                               0.072
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .678

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.57855209     0.53333333     0.85687688     0.07361564 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.2125934 0.1047798 0.7874066 0.9259250 0.16774945 0.5686452 0.7278368
#factor2 0.6244253 0.3166682 0.3755747 0.9329389 0.59293659 0.8472047 0.9142966
#general 0.5785521 0.5785521 0.5785521 0.8568769 0.07361564 0.9201194 0.9260624




################################################################
### 
### Simon et al (2022) https://bpspsychub.onlinelibrary.wiley.com/doi/10.1111/bjop.12587
### data experiment 3 https://osf.io/z28aq/?view_only=

library(readxl)
SDO_Simon3 <- read_excel("SDO_Simon3.xlsx")
colnames(SDO_Simon3)
mydata <- SDO_Simon3[,579:594]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .92, omega T = .94

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 6.73; eigenvalue 2 = 1.58

rho <- polychoric(mydata)$rho
# error; so tried another polychor command
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
# warning
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 9.51; eigenvalue 2 = 1.34

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.42, RMSEA=.165, RMSR=.12, TLI=.633

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.23, RMSR=.1, TLI=.622

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.111, RMSR=.05, TLI=.832
#      MR1   MR2
#MR1  1.00   -.52
#MR2  -.52   1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.181, RMSR=.05, TLI=.763
#      MR1   MR2
#MR1  1.00   -.66
#MR2  -.66   1.00


# Single factor model lavaan
UNImodel= '
 general =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8+
            SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.920
#Tucker-Lewis Index (TLI)                       0.977       0.908
#Robust Comparative Fit Index (CFI)                         0.699
#Robust Tucker-Lewis Index (TLI)                            0.652
#RMSEA                                          0.137       0.153
#Robust RMSEA                                               0.227
#SRMR                                           0.103       0.103

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .665


CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.681       0.637
#Tucker-Lewis Index (TLI)                       0.632       0.581
#Robust Comparative Fit Index (CFI)                         0.699
#Robust Tucker-Lewis Index (TLI)                            0.653
#RMSEA                                          0.169       0.121
#Robust RMSEA                                               0.158
#SRMR                                           0.112       0.112

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .413

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based; notice different order) 
EGAmodel= '
 factor1 =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8
 factor2 =~ SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.972
#Tucker-Lewis Index (TLI)                       0.995       0.967
#Robust Comparative Fit Index (CFI)                         0.836
#Robust Tucker-Lewis Index (TLI)                            0.809
#RMSEA                                          0.063       0.091
#Robust RMSEA                                               0.168
#SRMR                                           0.063       0.063

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .722

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.764    0.028  -27.487    0.000   -0.764   -0.764

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.847       0.823
#Tucker-Lewis Index (TLI)                       0.822       0.794
#Robust Comparative Fit Index (CFI)                         0.868
#Robust Tucker-Lewis Index (TLI)                            0.846
#RMSEA                                          0.118       0.085
#Robust RMSEA                                               0.106
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .54

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.625    0.120   -5.200    0.000   -0.625   -0.625

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.807       0.777
#Tucker-Lewis Index (TLI)                       0.777       0.743
#Robust Comparative Fit Index (CFI)                         0.828
#Robust Tucker-Lewis Index (TLI)                            0.801
#RMSEA                                          0.132       0.095
#Robust RMSEA                                               0.120
#SRMR                                           0.236       0.236

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .549


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8
 factor2 =~ SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
 general =~ SDO1+SDO2+SDO3+SDO4+SDO5+SDO6+SDO7+SDO8+
            SDO9+SDO10+SDO11+SDO12+SDO13+SDO14+SDO15+SDO16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.917       0.909
#Tucker-Lewis Index (TLI)                       0.887       0.876
#Robust Comparative Fit Index (CFI)                         0.935
#Robust Tucker-Lewis Index (TLI)                            0.911
#RMSEA                                          0.093       0.066
#Robust RMSEA                                               0.080
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .591

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    SDO1              0.617    0.161    3.838    0.000    0.617    0.617
#    SDO2              0.710    0.123    5.763    0.000    0.710    0.620
#    SDO3              0.454    0.207    2.195    0.028    0.454    0.607
#    SDO4              0.614    0.156    3.943    0.000    0.614    0.586
#    SDO5              0.613    0.129    4.746    0.000    0.613    0.618
#    SDO6              0.298    0.112    2.655    0.008    0.298    0.310
#    SDO7              0.760    0.121    6.292    0.000    0.760    0.557
#    SDO8              0.562    0.167    3.367    0.001    0.562    0.435
#  factor2 =~                                                            
#    SDO9              0.095    0.211    0.449    0.653    0.095    0.068
#    SDO10             0.384    0.183    2.101    0.036    0.384    0.254
#    SDO11             0.258    0.209    1.235    0.217    0.258    0.200
#    SDO12             0.471    0.167    2.813    0.005    0.471    0.351
#    SDO13             0.961    0.133    7.240    0.000    0.961    0.644
#    SDO14             0.030    0.150    0.197    0.844    0.030    0.023
#    SDO15             0.873    0.191    4.568    0.000    0.873    0.646
#    SDO16             0.639    0.162    3.952    0.000    0.639    0.525
#  general =~                                                            
#    SDO1              0.505    0.095    5.318    0.000    0.505    0.506
#    SDO2              0.593    0.117    5.086    0.000    0.593    0.518
#    SDO3              0.312    0.090    3.479    0.001    0.312    0.417
#    SDO4              0.501    0.114    4.396    0.000    0.501    0.478
#    SDO5              0.432    0.134    3.223    0.001    0.432    0.435
#    SDO6              0.710    0.140    5.076    0.000    0.710    0.737
#    SDO7              0.648    0.137    4.728    0.000    0.648    0.475
#    SDO8              0.518    0.123    4.224    0.000    0.518    0.401
#    SDO9             -1.050    0.141   -7.453    0.000   -1.050   -0.752
#    SDO10            -0.560    0.146   -3.828    0.000   -0.560   -0.370
#    SDO11            -1.096    0.154   -7.137    0.000   -1.096   -0.850
#    SDO12            -0.683    0.158   -4.313    0.000   -0.683   -0.510
#    SDO13            -0.659    0.161   -4.103    0.000   -0.659   -0.441
#    SDO14            -1.083    0.133   -8.119    0.000   -1.083   -0.855
#    SDO15            -0.734    0.155   -4.731    0.000   -0.734   -0.544
#    SDO16            -0.771    0.153   -5.025    0.000   -0.771   -0.634

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#    0.58514995     0.53333333     0.79917244     0.02867906 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.5455087 0.2683176 0.4544913 0.9083126 0.49592057 0.7897063 0.8907481
#factor2 0.2883741 0.1465324 0.7116259 0.9049999 0.20847949 0.6734106 0.8393012
#general 0.5851500 0.5851500 0.5851500 0.7991724 0.02867906 0.9205949 0.9423457




################################################################
### 
### Fleeson et al. (2023) https://onlinelibrary.wiley.com/doi/full/10.1111/jopy.12867
### data study 1: https://osf.io/uay2f

library(haven)
SDO_Fleeson1 <- read_sav("SDO_Fleeson1.sav")
colnames(SDO_Fleeson1)
mydata <- SDO_Fleeson1[,37:52]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .96, omega T = .97

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 9.25; eigenvalue 2 = 1.69

rho <- polychoric(mydata)$rho
# error
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 10.69; eigenvalue 2 = 1.51

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.58, RMSEA=.187, RMSR=.12, TLI=.704

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.67, RMSEA=.212, RMSR=.11, TLI=.711

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.7, RMSEA=.082, RMSR=.03, TLI=.944
#      MR1   MR2
#MR1  1.00  -.59
#MR2  -.59  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.78, RMSEA=.104, RMSR=.03, TLI=.93
#      MR1   MR2
#MR1  1.00  -.65
#MR2  -.65  1.00


# Single factor model lavaan
UNImodel= '
 general =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8+
            sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.983       0.915
#Tucker-Lewis Index (TLI)                       0.980       0.902
#Robust Comparative Fit Index (CFI)                         0.740
#Robust Tucker-Lewis Index (TLI)                            0.700
#RMSEA                                          0.209       0.210
#Robust RMSEA                                               0.222
#SRMR                                           0.114       0.114

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .726

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.743       0.721
#Tucker-Lewis Index (TLI)                       0.703       0.678
#Robust Comparative Fit Index (CFI)                         0.755
#Robust Tucker-Lewis Index (TLI)                            0.717
#RMSEA                                          0.190       0.141
#Robust RMSEA                                               0.183
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .597

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based) 
EGAmodel= '
 factor1 =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8
 factor2 =~ sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.957
#Tucker-Lewis Index (TLI)                       0.995       0.950
#Robust Comparative Fit Index (CFI)                         0.916
#Robust Tucker-Lewis Index (TLI)                            0.902
#RMSEA                                          0.103       0.149
#Robust RMSEA                                               0.127
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .778

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.769    0.024  -31.761    0.000   -0.769   -0.769

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.926       0.932
#Tucker-Lewis Index (TLI)                       0.914       0.921
#Robust Comparative Fit Index (CFI)                         0.942
#Robust Tucker-Lewis Index (TLI)                            0.933
#RMSEA                                          0.103       0.070
#Robust RMSEA                                               0.089
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .696

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.698    0.052  -13.531    0.000   -0.698   -0.698

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.887       0.887
#Tucker-Lewis Index (TLI)                       0.870       0.870
#Robust Comparative Fit Index (CFI)                         0.903
#Robust Tucker-Lewis Index (TLI)                            0.888
#RMSEA                                          0.126       0.089
#Robust RMSEA                                               0.115
#SRMR                                           0.335       0.335

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .700


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8
 factor2 =~ sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
 general =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8+
            sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.961       0.971
#Tucker-Lewis Index (TLI)                       0.947       0.960
#Robust Comparative Fit Index (CFI)                         0.976
#Robust Tucker-Lewis Index (TLI)                            0.967
#RMSEA                                          0.080       0.049
#Robust RMSEA                                               0.063
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .719

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7148750      0.5333333      0.8451820      0.0521471 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4195414 0.20384712 0.5804586 0.9472877 0.3978341 0.7728337 0.8960419
#factor2 0.1580916 0.08127792 0.8419084 0.9544910 0.1058827 0.5555785 0.8087622
#general 0.7148750 0.71487496 0.7148750 0.8451820 0.0521471 0.9522165 0.9577990





################################################################
### 
### Fleeson et al. (2023) https://onlinelibrary.wiley.com/doi/full/10.1111/jopy.12867
### data study 2: https://osf.io/uay2f

library(haven)
SDO_Fleeson2 <- read_sav("SDO_Fleeson2.sav")
colnames(SDO_Fleeson2)
mydata <- SDO_Fleeson2[,85:100]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .95, omega T = .97

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 9.07; eigenvalue 2 = 1.93

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 10.89; eigenvalue 2 = 1.63

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.57, RMSEA=.215, RMSR=.13, TLI=.63

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.68, RMSEA=.249, RMSR=.11, TLI=.646

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.71, RMSEA=.09, RMSR=.03, TLI=.934
#      MR1   MR2
#MR1  1.00  -.56
#MR2  -.56  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.8, RMSEA=.125, RMSR=.02, TLI=.91
#      MR1   MR2
#MR1  1.00  -.65
#MR2  -.65  1.00


# Single factor model lavaan
UNImodel= '
 general =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8+
            sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.985       0.941
#Tucker-Lewis Index (TLI)                       0.982       0.932
#Robust Comparative Fit Index (CFI)                         0.664
#Robust Tucker-Lewis Index (TLI)                            0.612
#RMSEA                                          0.213       0.186
#Robust RMSEA                                               0.267
#SRMR                                           0.128       0.128

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .73

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.684       0.662
#Tucker-Lewis Index (TLI)                       0.636       0.610
#Robust Comparative Fit Index (CFI)                         0.696
#Robust Tucker-Lewis Index (TLI)                            0.649
#RMSEA                                          0.216       0.155
#Robust RMSEA                                               0.209
#SRMR                                           0.137       0.137

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .557

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (theory based) 
EGAmodel= '
 factor1 =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8
 factor2 =~ sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.964
#Tucker-Lewis Index (TLI)                       0.996       0.958
#Robust Comparative Fit Index (CFI)                         0.897
#Robust Tucker-Lewis Index (TLI)                            0.880
#RMSEA                                          0.104       0.147
#Robust RMSEA                                               0.148
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .806

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.765    0.027  -28.379    0.000   -0.765   -0.765

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.919       0.922
#Tucker-Lewis Index (TLI)                       0.905       0.909
#Robust Comparative Fit Index (CFI)                         0.935
#Robust Tucker-Lewis Index (TLI)                            0.924
#RMSEA                                          0.110       0.075
#Robust RMSEA                                               0.097
#SRMR                                           0.079       0.079

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .700

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.621    0.063   -9.851    0.000   -0.621   -0.621

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.891       0.889
#Tucker-Lewis Index (TLI)                       0.874       0.872
#Robust Comparative Fit Index (CFI)                         0.907
#Robust Tucker-Lewis Index (TLI)                            0.892
#RMSEA                                          0.127       0.089
#Robust RMSEA                                               0.116
#SRMR                                           0.320       0.320

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .702


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8
 factor2 =~ sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
 general =~ sdo_1+sdo_2+sdo_3+sdo_4+sdo_5+sdo_6+sdo_7+sdo_8+
            sdo_9+sdo_10+sdo_11+sdo_12+sdo_13+sdo_14+sdo_15+sdo_16
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.950       0.948
#Tucker-Lewis Index (TLI)                       0.932       0.929
#Robust Comparative Fit Index (CFI)                         0.961
#Robust Tucker-Lewis Index (TLI)                            0.947
#RMSEA                                          0.093       0.066
#Robust RMSEA                                               0.081
#SRMR                                           0.044       0.044

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .747

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6499655      0.5333333      0.8681453      0.0122098 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2754160 0.1391569 0.7245840 0.9523862 0.2331942 0.7236310 0.8435613
#factor2 0.4262398 0.2108776 0.5737602 0.9506894 0.3945632 0.7949130 0.9011570
#general 0.6499655 0.6499655 0.6499655 0.8681453 0.0122098 0.9407045 0.9377826






