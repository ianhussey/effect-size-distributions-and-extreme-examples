################################################################
### Analysis Edinburgh Handedness Inventory

### No clear indications for a two-factor scale
### split as suggested by Christman et al. (2015) https://www.sciencedirect.com/science/article/pii/S0278262615000676
### writing and drawing vs. the rest




################################################################

### Dragan et al. 2022 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0262803#sec013

Dragan_handedness <- read.delim("Dragan_handedness.dat")
colnames(Dragan_handedness)
mydata  <- as.data.frame(Dragan_handedness[,5:14])
colnames(mydata)
mydata <- na.omit(mydata)
colnames(mydata) <- c("EHI_1","EHI_2","EHI_3","EHI_4","EHI_5","EHI_6","EHI_7","EHI_8","EHI_9","EHI_10")
summary(mydata)
min(mydata)
max(mydata) # 3 response alternatives
names = colnames(mydata)


library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 4.39; eigenvalue 2 = 0.59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 5.92; eigenvalue 2 = 0.52

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.44, RMSEA=.133, RMSR=.07, TLI=.832

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.162, RMSR=.07, TLI=.849

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.34, RMSEA=.067, RMSR=.03, TLI=.958
#     MR1  MR2
#MR1 1.00 0.56
#MR2 0.56 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.67, RMSEA=.106, RMSR=.03, TLI=.935
#     MR1  MR2
#MR1 1.00 0.58
#MR2 0.58 1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~ EHI_1+EHI_2+EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.981
#Tucker-Lewis Index (TLI)                       0.994       0.976
#Robust Comparative Fit Index (CFI)                         0.897
#Robust Tucker-Lewis Index (TLI)                            0.868
#RMSEA                                          0.073       0.104
#Robust RMSEA                                               0.153
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .643

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.869
#Tucker-Lewis Index (TLI)                       0.854       0.832
#Robust Comparative Fit Index (CFI)                         0.895
#Robust Tucker-Lewis Index (TLI)                            0.865
#RMSEA                                          0.124       0.082
#Robust RMSEA                                               0.118
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .415

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 2 factors based on Christman et al.
EGAmodel= '
 factor1=~ EHI_1+EHI_2
 factor2=~ EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.982
#Tucker-Lewis Index (TLI)                       0.995       0.976
#Robust Comparative Fit Index (CFI)                         0.908
#Robust Tucker-Lewis Index (TLI)                            0.878
#RMSEA                                          0.072       0.103
#Robust RMSEA                                               0.147
#SRMR                                           0.065       0.065

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.965    0.007  129.494    0.000    0.965    0.965

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .650

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.892       0.869
#Tucker-Lewis Index (TLI)                       0.856       0.826
#Robust Comparative Fit Index (CFI)                         0.899
#Robust Tucker-Lewis Index (TLI)                            0.867
#RMSEA                                          0.123       0.083
#Robust RMSEA                                               0.117
#SRMR                                           0.071       0.071

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.958    0.028   34.053    0.000    0.958    0.958

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .424

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# does not converge either
summary(CFA_model3)
# give max values

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# min %variance given = .2


# Bifactor model
BIFmodel= '
 factor1=~ EHI_1+EHI_2
 factor2=~ EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
 general=~ EHI_1+EHI_2+EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.974       0.972
#Tucker-Lewis Index (TLI)                       0.952       0.949
#Robust Comparative Fit Index (CFI)                         0.980
#Robust Tucker-Lewis Index (TLI)                            0.963
#RMSEA                                          0.071       0.045
#Robust RMSEA                                               0.062
#SRMR                                           0.025       0.025
fitRsquares <- lavInspect(CFA_model4, what='rsquare')

summary(fitRsquares)[-4]
# median %variance explained = .424

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    EHI_1             0.087    0.074    1.165    0.244    0.087    0.181
#    EHI_2            -0.013    0.158   -0.082    0.935   -0.013   -0.027
#  factor2 =~                                                            
#    EHI_3             0.132    0.034    3.832    0.000    0.132    0.227
#    EHI_4             0.312    0.042    7.490    0.000    0.312    0.444
#    EHI_5             0.000    0.033    0.004    0.997    0.000    0.000
#    EHI_6             0.445    0.048    9.208    0.000    0.445    0.562
#    EHI_7             0.296    0.043    6.954    0.000    0.296    0.369
#    EHI_8             0.149    0.050    2.976    0.003    0.149    0.253
#    EHI_9             0.459    0.044   10.404    0.000    0.459    0.531
#    EHI_10            0.153    0.035    4.344    0.000    0.153    0.230
#  general =~                                                            
#    EHI_1             0.440    0.016   27.541    0.000    0.440    0.921
#    EHI_2             0.425    0.017   24.852    0.000    0.425    0.876
#    EHI_3             0.380    0.021   18.384    0.000    0.380    0.652
#    EHI_4             0.339    0.021   16.448    0.000    0.339    0.481
#    EHI_5             0.434    0.015   27.981    0.000    0.434    0.860
#    EHI_6             0.342    0.028   12.308    0.000    0.342    0.431
#    EHI_7             0.197    0.028    7.034    0.000    0.197    0.246
#    EHI_8             0.387    0.025   15.479    0.000    0.387    0.658
#    EHI_9             0.205    0.030    6.777    0.000    0.205    0.237
#    EHI_10            0.384    0.023   17.012    0.000    0.384    0.579

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.7827558      0.3555556      0.8981196      0.7518133 

#$FactorLevelIndices
#         ECV_SS      ECV_SG    ECV_GS     Omega      OmegaH          H        FD
#factor1 0.02038505 0.006444509 0.9796150 0.9025923 0.006652374 0.03356855 0.3843918
#factor2 0.30824937 0.210799721 0.6917506 0.8442285 0.240568536 0.58945913 0.7880430
#general 0.78275577 0.782755769 0.7827558 0.8981196 0.751813325 0.93492452 0.9654132




################################################################

### McManus sent by email

### select items in order of the EHI and replace two missing items
### holding spoon replaced by holding racket
### opening lid replaced by pull trigger of a rifle

library(haven)
McManus_handedness <- read_sav("McManus_handedness.sav")
colnames(McManus_handedness)
mydata  <- as.data.frame(McManus_handedness[,c(2,14,3,8,12,17,4,5,7,27)])
colnames(mydata) 
mydata <- na.omit(mydata)
colnames(mydata) <- c("EHI_1","EHI_2","EHI_3","EHI_4","EHI_5","EHI_6","EHI_7","EHI_8","EHI_9","EHI_10")
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 1 component
# Eigenvalue 1 = 7.24; eigenvalue 2 = 0.17

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 8.17; eigenvalue 2 = 0.11

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.72, RMSEA=.204, RMSR=.03, TLI=.844

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.82, RMSEA=.733, RMSR=.02, TLI=.302

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community without response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.75, RMSEA=.078, RMSR=.02, TLI=.977
#     MR1  MR2
#MR1 1.00 0.22
#MR2 0.22 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.87, RMSEA=.844, RMSR=.01, TLI=.075
#     MR1  MR2
#MR1 1.00 0.72
#MR2 0.72 1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~ EHI_1+EHI_2+EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
'

library(lavaan)
# ordered does not converge for all models; hence MLR
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.899       0.848
#Tucker-Lewis Index (TLI)                       0.870       0.805
#Robust Comparative Fit Index (CFI)                         0.902
#Robust Tucker-Lewis Index (TLI)                            0.874
#RMSEA                                          0.186       0.120
#Robust RMSEA                                               0.183
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .688

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on EGA analysis
EGAmodel= '
 factor1=~ EHI_1+EHI_2
 factor2=~ EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.981       0.980
#Tucker-Lewis Index (TLI)                       0.974       0.973
#Robust Comparative Fit Index (CFI)                         0.985
#Robust Tucker-Lewis Index (TLI)                            0.981
#RMSEA                                          0.083       0.045
#Robust RMSEA                                               0.071
#SRMR                                           0.018       0.018


#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.925    0.012   75.112    0.000    0.925    0.925


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .742

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.842       0.798
#Tucker-Lewis Index (TLI)                       0.797       0.741
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.802
#RMSEA                                          0.233       0.139
#Robust RMSEA                                               0.229
#SRMR                                           0.410       0.410

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .736

# Bifactor model
BIFmodel= '
 factor1=~ EHI_1+EHI_2
 factor2=~ EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
 general=~ EHI_1+EHI_2+EHI_3+EHI_4+EHI_5+EHI_6+EHI_7+EHI_8+EHI_9+EHI_10
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="WLSM",orthogonal=TRUE)
# no model with other variables converged either; so took max values
library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)

