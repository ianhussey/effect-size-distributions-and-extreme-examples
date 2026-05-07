################################################################
### Analysis Questionnaires DASS-21


################################################################
### DASS21 Hua Chen https://www.sciencedirect.com/science/article/pii/S0001691823002184?dgcid=rss_sd_all
### data https://osf.io/x3f78/
### Cross-section adults

library(haven)
DASS21_IHuaChen_osf <- read_sav("DASS21_IHuaChen_osf.sav")
colnames(DASS21_IHuaChen_osf)
mydata <- DASS21_IHuaChen_osf[DASS21_IHuaChen_osf$group==3,]
mydata <- mydata[,2:22]
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) #response alternatives = 4

colnames2 <- c("stress1","anx1","dep1","anx2","dep2","stress2","anx3","stress3","anx4","dep3",
               "stress4","stress5","dep4","stress6","anx5","dep5","dep6","stress7","anx6","anx7","dep7")
colnames(mydata) <- colnames2

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 8.5; eigenvalue 2 = 0.59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 4 component
# Eigenvalue 1 = 10.45; eigenvalue 2 = .68

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.40, RMSEA=.074, RMSR=.05, TLI=.874

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.114, RMSR=.06, TLI=.806

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities 

fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.051, RMSR=.03, TLI=.94
#     MR1  MR3  MR2
#MR1 1.00 0.74 0.31
#MR3 0.74 1.00 0.28
#MR2 0.31 0.28 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.091, RMSR=.04, TLI=.874
#     MR1  MR2  MR3
#MR1 1.00 0.74 0.26
#MR2 0.74 1.00 0.24
#MR3 0.26 0.24 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ anx1+anx2+anx3+anx4+anx5+anx6+anx7+
            dep1+dep2+dep3+dep4+dep5+dep6+dep7+
            stress1+stress2+stress3+stress4+stress5+stress6+stress7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.992       0.958
#Tucker-Lewis Index (TLI)                       0.991       0.953
#Robust Comparative Fit Index (CFI)                         0.847
#Robust Tucker-Lewis Index (TLI)                            0.830
#RMSEA                                          0.056       0.068
#Robust RMSEA                                               0.108
#SRMR                                           0.058       0.058

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .50

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.898
#Tucker-Lewis Index (TLI)                       0.873       0.887
#Robust Comparative Fit Index (CFI)                         0.900
#Robust Tucker-Lewis Index (TLI)                            0.888
#RMSEA                                          0.075       0.060
#Robust RMSEA                                               0.070
#SRMR                                           0.050       0.050

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .400

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on theoretical analysis
EGAmodel= '
 anxiety=~ anx1+anx2+anx3+anx4+anx5+anx6+anx7
 depression=~ dep1+dep2+dep3+dep4+dep5+dep6+dep7
 stress=~ stress1+stress2+stress3+stress4+stress5+stress6+stress7
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.993       0.961
#Tucker-Lewis Index (TLI)                       0.992       0.956
#Robust Comparative Fit Index (CFI)                         0.862
#Robust Tucker-Lewis Index (TLI)                            0.845
#RMSEA                                          0.053       0.066
#Robust RMSEA                                               0.103
#SRMR                                           0.055       0.055

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#anxiety ~~                                                            
#  depression        0.935    0.014   66.342    0.000    0.935    0.935
#  stress            0.961    0.013   75.388    0.000    0.961    0.961
#depression ~~                                                         
#  stress            0.955    0.012   79.001    0.000    0.955    0.955

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .53

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.893       0.904
#Tucker-Lewis Index (TLI)                       0.879       0.892
#Robust Comparative Fit Index (CFI)                         0.906
#Robust Tucker-Lewis Index (TLI)                            0.894
#RMSEA                                          0.074       0.059
#Robust RMSEA                                               0.069
#SRMR                                           0.049       0.049

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  anxiety ~~                                                            
#    depression        0.935    0.022   42.073    0.000    0.935    0.935
#    stress            0.950    0.020   47.696    0.000    0.950    0.950
#  depression ~~                                                         
#    stress            0.950    0.018   52.494    0.000    0.950    0.950


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.659       0.662
#Tucker-Lewis Index (TLI)                       0.621       0.625
#Robust Comparative Fit Index (CFI)                         0.668
#Robust Tucker-Lewis Index (TLI)                            0.632
#RMSEA                                          0.130       0.109
#Robust RMSEA                                               0.128
#SRMR                                           0.320       0.320

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .450

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 anxiety=~ anx1+anx2+anx3+anx4+anx5+anx6+anx7
 depression=~ dep1+dep2+dep3+dep4+dep5+dep6+dep7
 stress=~ stress1+stress2+stress3+stress4+stress5+stress6+stress7
 general=~ anx1+anx2+anx3+anx4+anx5+anx6+anx7+
            dep1+dep2+dep3+dep4+dep5+dep6+dep7+
            stress1+stress2+stress3+stress4+stress5+stress6+stress7
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.925       0.931
#Tucker-Lewis Index (TLI)                       0.907       0.914
#Robust Comparative Fit Index (CFI)                         0.936
#Robust Tucker-Lewis Index (TLI)                            0.919
#RMSEA                                          0.065       0.052
#Robust RMSEA                                               0.060
#SRMR                                           0.044       0.044

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .469

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.8492988      0.7000000      0.9400433      0.9141655 
#
#$FactorLevelIndices
#              ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#anxiety    0.1905409 0.06255145 0.8094591 0.8408250 0.09285392 0.4305565 0.7316500
#depression 0.1467947 0.05111880 0.8532053 0.8581284 0.06729732 0.3888428 0.7028256
#stress     0.1144759 0.03703099 0.8855241 0.8400153 0.03866522 0.3061205 0.6235535
#general    0.8492988 0.84929876 0.8492988 0.9400433 0.91416553 0.9386354 0.9654529









################################################################
### DASS21 Thiyagarajan  https://doi.org/10.7910/DVN/AJVLN6

library(readxl)
dass21_Thiyagarajan <- read_excel("dass21_Thiyagarajan.xlsx")
colnames(dass21_Thiyagarajan)
mydata <- dass21_Thiyagarajan[,9:29]
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 4 response alternatives

colnames2 <- c("stress1","anx1","dep1","anx2","dep2","stress2","anx3","stress3","anx4","dep3",
               "stress4","stress5","dep4","stress6","anx5","dep5","dep6","stress7","anx6","anx7","dep7")
colnames(mydata) <- colnames2
summary(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 11.17; eigenvalue 2 = 0.67

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 1 component
# Eigenvalue 1 = 12.74; eigenvalue 2 = .70

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.53, RMSEA=.088, RMSR=.04, TLI=.889

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.113, RMSR=.05, TLI=.858

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community

fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.057, RMSR=.02, TLI=.953
#     MR1  MR3  MR2
#MR1 1.00 0.74 0.52
#MR3 0.74 1.00 0.53
#MR2 0.52 0.53 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.67, RMSEA=.08, RMSR=.02, TLI=.929
#     MR1  MR3  MR2
#MR1 1.00 0.75 0.55
#MR3 0.75 1.00 0.57
#MR2 0.55 0.57 1.00

# Single factor model lavaan
UNImodel= '
 general=~ anx1+anx2+anx3+anx4+anx5+anx6+anx7+
            dep1+dep2+dep3+dep4+dep5+dep6+dep7+
            stress1+stress2+stress3+stress4+stress5+stress6+stress7
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.995       0.962
#Tucker-Lewis Index (TLI)                       0.995       0.958
#Robust Comparative Fit Index (CFI)                         0.877
#Robust Tucker-Lewis Index (TLI)                            0.864
#RMSEA                                          0.068       0.092
#Robust RMSEA                                               0.112
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .66

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.900       0.902
#Tucker-Lewis Index (TLI)                       0.889       0.891
#Robust Comparative Fit Index (CFI)                         0.903
#Robust Tucker-Lewis Index (TLI)                            0.892
#RMSEA                                          0.088       0.076
#Robust RMSEA                                               0.086
#SRMR                                           0.043       0.043

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .571

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on theoretical analysis
EGAmodel= '
 anxiety=~ anx1+anx2+anx3+anx4+anx5+anx6+anx7
 depression=~ dep1+dep2+dep3+dep4+dep5+dep6+dep7
 stress=~ stress1+stress2+stress3+stress4+stress5+stress6+stress7
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.997       0.973
#Tucker-Lewis Index (TLI)                       0.997       0.969
#Robust Comparative Fit Index (CFI)                         0.910
#Robust Tucker-Lewis Index (TLI)                            0.899
#RMSEA                                          0.054       0.078
#Robust RMSEA                                               0.097
#SRMR                                           0.038       0.038

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#anxiety ~~                                                            
#  depression        0.890    0.010   86.498    0.000    0.890    0.890
#  stress            0.990    0.006  168.611    0.000    0.990    0.990
#depression ~~                                                         
#  stress            0.935    0.008  123.211    0.000    0.935    0.935

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .687

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.930       0.932
#Tucker-Lewis Index (TLI)                       0.921       0.923
#Robust Comparative Fit Index (CFI)                         0.934
#Robust Tucker-Lewis Index (TLI)                            0.925
#RMSEA                                          0.074       0.064
#Robust RMSEA                                               0.072
#SRMR                                           0.037       0.037

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .588

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  anxiety ~~                                                            
#    depression        0.875    0.015   60.317    0.000    0.875    0.875
#    stress            0.985    0.008  119.692    0.000    0.985    0.985
#  depression ~~                                                         
#    stress            0.923    0.011   81.321    0.000    0.923    0.923


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.740       0.737
#Tucker-Lewis Index (TLI)                       0.711       0.708
#Robust Comparative Fit Index (CFI)                         0.743
#Robust Tucker-Lewis Index (TLI)                            0.714
#RMSEA                                          0.142       0.125
#Robust RMSEA                                               0.141
#SRMR                                           0.415       0.415

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .584

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 anxiety=~ NA*anx1+anx2+anx3+anx4+anx5+anx6+anx7
 depression=~ NA*dep1+dep2+dep3+dep4+dep5+dep6+dep7
 stress=~ NA*stress1+stress2+stress3+stress4+stress5+stress6+stress7
 general=~ NA*anx1+anx2+anx3+anx4+anx5+anx6+anx7+
            dep1+dep2+dep3+dep4+dep5+dep6+dep7+
            stress1+stress2+stress3+stress4+stress5+stress6+stress7
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
#does not converge
#leave stress3 out
BIFmodel= '
 anxiety=~ NA*anx1+anx2+anx3+anx4+anx5+anx6+anx7
 depression=~ NA*dep1+dep2+dep3+dep4+dep5+dep6+dep7
 stress=~ NA*stress1+stress2+stress4+stress5+stress6+stress7
 general=~ NA*anx1+anx2+anx3+anx4+anx5+anx6+anx7+
            dep1+dep2+dep3+dep4+dep5+dep6+dep7+
            stress1+stress2+stress4+stress5+stress6+stress7
'
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.956
#Tucker-Lewis Index (TLI)                       0.944       0.944
#Robust Comparative Fit Index (CFI)                         0.958
#Robust Tucker-Lewis Index (TLI)                            0.947
#RMSEA                                          0.063       0.055
#Robust RMSEA                                               0.061
#SRMR                                           0.029       0.029

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .588

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.8794299      0.7000000      0.9621109      0.9328501 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#anxiety    0.12852698 0.04052338 0.8714730 0.8815451 0.093985033 0.3540771 0.6690708
#depression 0.15289126 0.05877606 0.8471087 0.9240059 0.108230300 0.4688152 0.8037975
#stress     0.07083657 0.02127071 0.9291634 0.8882846 0.009947752 0.2387387 0.7761043
#general    0.87942985 0.87942985 0.8794299 0.9621109 0.932850116 0.9594325 0.9737070

