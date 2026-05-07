################################################################
### analysis Sexual Desire Inventory (SDI-2) https://vk-kustannus.fi/wp-content/uploads/2021/09/The-sexual-desire-inventory-Development-factor-structure-and-evidence-of-reliability.pdf
### Jones & DeBruine https://osf.io/psmqg

SDI.2_Jones <- read.delim("SDI-2_Jones.txt")
colnames(SDI.2_Jones)
mydata <- as.data.frame(SDI.2_Jones[,5:18])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 10 response alternatives


library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 5.45; eigenvalue 2 = 1.62

rho <- polychoric(mydata)$rho
# does not run because 8 answer alternatives
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 5.69; eigenvalue 2 = 1.61


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.39, RMSEA=.212, RMSR=.14, TLI=.499

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.222, RMSR=.14, TLI=.492

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities 

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.124, RMSR=.06, TLI=.827
#     MR1  MR2
#MR1 1.00 0.41
#MR2 0.41 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.137, RMSR=.06, TLI=.807
#     MR1  MR2
#MR1 1.00 0.43
#MR2 0.43 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  SDI1+SDI2+SDI3+SDI4+SDI5+SDI6+SDI7+SDI8+SDI9+SDI10+SDI11+SDI12+SDI13+SDI14
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.907       0.750
#Tucker-Lewis Index (TLI)                       0.890       0.704
#Robust Comparative Fit Index (CFI)                         0.508
#Robust Tucker-Lewis Index (TLI)                            0.418
#RMSEA                                          0.254       0.248
#Robust RMSEA                                               0.237
#SRMR                                           0.163       0.163

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .447

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.579       0.557
#Tucker-Lewis Index (TLI)                       0.503       0.476
#Robust Comparative Fit Index (CFI)                         0.580
#Robust Tucker-Lewis Index (TLI)                            0.503
#RMSEA                                          0.211       0.190
#Robust RMSEA                                               0.211
#SRMR                                           0.137       0.137

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .337

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 dyadic =~  SDI1+SDI2+SDI3+SDI4+SDI5+SDI6+SDI7+SDI8+SDI9
 solitary=~ SDI10+SDI11+SDI12+SDI13+SDI14
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.965       0.879
#Tucker-Lewis Index (TLI)                       0.959       0.855
#Robust Comparative Fit Index (CFI)                         0.800
#Robust Tucker-Lewis Index (TLI)                            0.760
#RMSEA                                          0.156       0.173
#Robust RMSEA                                               0.152
#SRMR                                           0.095       0.095

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#dyadic ~~                                                             
#  solitary          0.531    0.008   66.683    0.000    0.531    0.531

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .530

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.847       0.841
#Tucker-Lewis Index (TLI)                       0.817       0.809
#Robust Comparative Fit Index (CFI)                         0.848
#Robust Tucker-Lewis Index (TLI)                            0.817
#RMSEA                                          0.128       0.115
#Robust RMSEA                                               0.128
#SRMR                                           0.096       0.096

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#dyadic ~~                                                             
#  solitary          0.442    0.011   38.769    0.000    0.442    0.442

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .529

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.826       0.818
#Tucker-Lewis Index (TLI)                       0.794       0.784
#Robust Comparative Fit Index (CFI)                         0.826
#Robust Tucker-Lewis Index (TLI)                            0.794
#RMSEA                                          0.136       0.122
#Robust RMSEA                                               0.136
#SRMR                                           0.201       0.201

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .529

# Bifactor model
BIFmodel= '
 dyadic =~  SDI1+SDI2+SDI3+SDI4+SDI5+SDI6+SDI7+SDI8+SDI9
 solitary=~ SDI10+SDI11+SDI12+SDI13+SDI14
 general=~  SDI1+SDI2+SDI3+SDI4+SDI5+SDI6+SDI7+SDI8+SDI9+SDI10+SDI11+SDI12+SDI13+SDI14
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.941       0.938
#Tucker-Lewis Index (TLI)                       0.915       0.910
#Robust Comparative Fit Index (CFI)                         0.941
#Robust Tucker-Lewis Index (TLI)                            0.915
#RMSEA                                          0.087       0.079
#Robust RMSEA                                               0.087
#SRMR                                           0.045       0.045

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .616

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6070171      0.4945055      0.9337695      0.8037939 
#
#$FactorLevelIndices
#            ECV_SS    ECV_SG    ECV_GS     Omega       OmegaH         H        FD
#dyadic   0.1789073 0.1093609 0.8210927 0.9072039 5.289019e-06 0.5659637 0.8839978
#solitary 0.7296141 0.2836220 0.2703859 0.8987600 6.422886e-01 0.8589365 0.9466365
#general  0.6070171 0.6070171 0.6070171 0.9337695 8.037939e-01 0.9090509 0.9610262



