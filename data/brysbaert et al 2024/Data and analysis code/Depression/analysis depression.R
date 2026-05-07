################################################################
### Analysis Beck's Depression Inventory II https://naviauxlab.ucsd.edu/wp-content/uploads/2020/09/BDI21.pdf


### CFA based on the overall analysis Whisman et al. (2012)
### https://journals.sagepub.com/doi/10.1177/1073191112460273
### However because most BIF models did not converge, perfdif and somatic combined into one factor


################################################################
### https://osf.io/m4q8j/
### Perfectionism and depression Smith https://cruxpsychology.ca/wp-content/uploads/2017/07/Smithetal.2016.670-687.pdf 


library(haven)
BDI_II_Smith <- read_sav("BDI-II_Smith.sav")
colnames(BDI_II_Smith)
mydata <- as.data.frame(BDI_II_Smith[,49:68])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 component
# Eigenvalue 1 = 7.74; eigenvalue 2 = .62

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 component
# Eigenvalue 1 = 9.13; eigenvalue 2 = .77

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.39, RMSEA=.059, RMSR=.06, TLI=.914

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.122, RMSR=.07, TLI=.764

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communities with response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.046, RMSR=.04, TLI=.947
#     MR1  MR2
#MR1 1.00 0.57
#MR2 0.57 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.114, RMSR=.06, TLI=.791
#     MR1  MR2
#MR1 1.00 0.58
#MR2 0.58 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  BDI_II_1+BDI_II_2+BDI_II_3+BDI_iI_4+BDI_II_5+BDI_II_6+BDI_II_7+BDI_II_8+BDI_II_10+
            BDI_II_11+BDI_II_12+BDI_II_13+BDI_II_14+BDI_II_15+BDI_II_16+BDI_II_17+BDI_II_18+
            BDI_II_19+BDI_II_20+BDI_II_21
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.966
#Tucker-Lewis Index (TLI)                       0.999       0.962
#Robust Comparative Fit Index (CFI)                         0.857
#Robust Tucker-Lewis Index (TLI)                            0.840
#RMSEA                                          0.020       0.057
#Robust RMSEA                                               0.103
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .496

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.918       0.935
#Tucker-Lewis Index (TLI)                       0.909       0.927
#Robust Comparative Fit Index (CFI)                         0.936
#Robust Tucker-Lewis Index (TLI)                            0.928
#RMSEA                                          0.062       0.050
#Robust RMSEA                                               0.055
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .407

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 negatt =~  BDI_II_1+BDI_II_2+BDI_II_3+BDI_II_5+BDI_II_6+BDI_II_7+BDI_II_8+BDI_II_10+
            BDI_II_14
 perfsom=~  BDI_iI_4+BDI_II_11+BDI_II_12+BDI_II_13+BDI_II_17+BDI_II_19+BDI_II_21+
            BDI_II_15+BDI_II_16+BDI_II_18+BDI_II_20
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.971
#Tucker-Lewis Index (TLI)                       1.000       0.968
#Robust Comparative Fit Index (CFI)                         0.874
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.000       0.052
#Robust RMSEA                                               0.097
#SRMR                                           0.066       0.066

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negatt ~~                                                             
#       perfsom           0.921    0.022   42.810    0.000    0.921    0.921

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .511

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.931       0.948
#Tucker-Lewis Index (TLI)                       0.923       0.942
#Robust Comparative Fit Index (CFI)                         0.949
#Robust Tucker-Lewis Index (TLI)                            0.942
#RMSEA                                          0.057       0.045
#Robust RMSEA                                               0.049
#SRMR                                           0.052       0.052

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negatt ~~                                                             
#       perfsom           0.922    0.027   34.593    0.000    0.922    0.922

# median %variance explained = .420

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.799       0.810
#Tucker-Lewis Index (TLI)                       0.775       0.787
#Robust Comparative Fit Index (CFI)                         0.813
#Robust Tucker-Lewis Index (TLI)                            0.791
#RMSEA                                          0.098       0.086
#Robust RMSEA                                               0.093
#SRMR                                           0.263       0.263

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .405

# Bifactor model
BIFmodel= '
 negatt =~  BDI_II_1+BDI_II_2+BDI_II_3+BDI_II_5+BDI_II_6+BDI_II_7+BDI_II_8+BDI_II_10+
            BDI_II_14
 perfsom=~  BDI_iI_4+BDI_II_11+BDI_II_12+BDI_II_13+BDI_II_17+BDI_II_19+BDI_II_21+
            BDI_II_15+BDI_II_16+BDI_II_18+BDI_II_20
 general=~  BDI_II_1+BDI_II_2+BDI_II_3+BDI_iI_4+BDI_II_5+BDI_II_6+BDI_II_7+BDI_II_8+BDI_II_10+
            BDI_II_11+BDI_II_12+BDI_II_13+BDI_II_14+BDI_II_15+BDI_II_16+BDI_II_17+BDI_II_18+
            BDI_II_19+BDI_II_20+BDI_II_21
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.953       0.955
#Tucker-Lewis Index (TLI)                       0.941       0.943
#Robust Comparative Fit Index (CFI)                         0.959
#Robust Tucker-Lewis Index (TLI)                            0.948
#RMSEA                                          0.050       0.045
#Robust RMSEA                                               0.046
#SRMR                                           0.045       0.045

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .448

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.8482790      0.5210526      0.9280130      0.8870591 
#$FactorLevelIndices
#        ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#negatt  0.1625272 0.08213160 0.8374728 0.8831496 0.09856025 0.5082797 0.7817260
#perfsom 0.1406815 0.06958941 0.8593185 0.8468839 0.03739584 0.4858818 0.7852538
#general 0.8482790 0.84827899 0.8482790 0.9280130 0.88705911 0.9483822 0.9668151







################################################################
### Sandoval-Lentisco et al Spanish https://osf.io/x5gwn

library(readxl)
depression_Sandoval <- read_excel("depression_Sandoval.xlsx")
colnames(depression_Sandoval)
mydata <- as.data.frame(depression_Sandoval[,27:47])
colnames(mydata)
mydata[mydata<0] <- NA
mydata <- na.omit(mydata)
colnames(mydata) <- c("BDI_1","BDI_2","BDI_3","BDI_4","BDI_5","BDI_6","BDI_7","BDI_8","BDI_9","BDI_10",
                      "BDI_11","BDI_12","BDI_13","BDI_14","BDI_15","BDI_16","BDI_17","BDI_18","BDI_19",
                      "BDI_20","BDI_21")
summary(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 6.98; eigenvalue 2 = .82

rho <- polychoric(mydata)$rho
# Does not run because max is not the same in each column
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 component
# Eigenvalue 1 = 9.34; eigenvalue 2 = 1.04


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.33, RMSEA=.081, RMSR=.07, TLI=.805

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.171, RMSR=.08, TLI=.586

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities when no response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.067, RMSR=.05, TLI=.865
#     MR1  MR2
#MR1 1.00 0.64
#MR2 0.64 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.163, RMSR=.06, TLI=.621
#     MR1  MR2
#MR1 1.00 0.65
#MR2 0.65 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  BDI_1+BDI_2+BDI_3+BDI_4+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_11+BDI_12+
            BDI_13+BDI_14+BDI_15+BDI_16+BDI_17+BDI_18+BDI_19+BDI_20+BDI_21
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.982       0.933
#Tucker-Lewis Index (TLI)                       0.980       0.925
#Robust Comparative Fit Index (CFI)                         0.701
#Robust Tucker-Lewis Index (TLI)                            0.667
#RMSEA                                          0.061       0.072
#Robust RMSEA                                               0.155
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .465

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.822       0.840
#Tucker-Lewis Index (TLI)                       0.802       0.822
#Robust Comparative Fit Index (CFI)                         0.841
#Robust Tucker-Lewis Index (TLI)                            0.823
#RMSEA                                          0.083       0.070
#Robust RMSEA                                               0.078
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .337

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 negatt  =~  BDI_1+BDI_2+BDI_3+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_14
 perfsom =~  BDI_4+BDI_11+BDI_12+BDI_13+BDI_17+BDI_19+BDI_21+
             BDI_15+BDI_16+BDI_18+BDI_20
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.990       0.954
#Tucker-Lewis Index (TLI)                       0.989       0.949
#Robust Comparative Fit Index (CFI)                         0.743
#Robust Tucker-Lewis Index (TLI)                            0.713
#RMSEA                                          0.045       0.060
#Robust RMSEA                                               0.144
#SRMR                                           0.073       0.073

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#negatt ~~                                                             
#   perfsom           0.832    0.028   29.338    0.000    0.832    0.832

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .504

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.860       0.879
#Tucker-Lewis Index (TLI)                       0.844       0.864
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.866
#RMSEA                                          0.074       0.061
#Robust RMSEA                                               0.068
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .394

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#negatt ~~                                                             
#   perfsom           0.828    0.037   22.121    0.000    0.828    0.828

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.766       0.780
#Tucker-Lewis Index (TLI)                       0.740       0.756
#Robust Comparative Fit Index (CFI)                         0.783
#Robust Tucker-Lewis Index (TLI)                            0.759
#RMSEA                                          0.096       0.082
#Robust RMSEA                                               0.091
#SRMR                                           0.219       0.219

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .380

# Bifactor model
BIFmodel= '
 negatt  =~  BDI_1+BDI_2+BDI_3+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_14
 perfsom =~  BDI_4+BDI_11+BDI_12+BDI_13+BDI_17+BDI_19+BDI_21+
             BDI_15+BDI_16+BDI_18+BDI_20
 general=~  BDI_1+BDI_2+BDI_3+BDI_4+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_11+BDI_12+
            BDI_13+BDI_14+BDI_15+BDI_16+BDI_17+BDI_18+BDI_19+BDI_20+BDI_21
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.897       0.902
#Tucker-Lewis Index (TLI)                       0.872       0.878
#Robust Comparative Fit Index (CFI)                         0.909
#Robust Tucker-Lewis Index (TLI)                            0.887
#RMSEA                                          0.067       0.058
#Robust RMSEA                                               0.062
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .399

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7833577      0.5238095      0.9202935      0.8538541 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#negatt  0.1425221 0.07547654 0.8574779 0.8755581 0.03065486 0.4287210 0.7306149
#perfsom 0.3000832 0.14116579 0.6999168 0.8474289 0.19775047 0.6086122 0.7967014
#general 0.7833577 0.78335767 0.7833577 0.9202935 0.85385409 0.9106994 0.9475513







################################################################
### García-Batista et al Dominican Republic 
### #https://figshare.com/articles/dataset/Validity_and_reliability_of_the_Beck_Depression_Inventory_BDI-II_in_general_and_hospital_population_of_Dominican_Republic/6725552

library(haven)
BDI_II_Garcia_Batista <- read_sav("BDI-II_Garcia_Batista.sav")
colnames(BDI_II_Garcia_Batista)
mydata <- as.data.frame(BDI_II_Garcia_Batista[,3:23])
colnames(mydata)
mydata <- na.omit(mydata)
colnames(mydata) <- c("BDI_1","BDI_2","BDI_3","BDI_4","BDI_5","BDI_6","BDI_7","BDI_8","BDI_9","BDI_10",
                      "BDI_11","BDI_12","BDI_13","BDI_14","BDI_15","BDI_16","BDI_17","BDI_18","BDI_19",
                      "BDI_20","BDI_21")
summary(mydata)
min(mydata)
max(mydata) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 2 components
# Eigenvalue 1 = 6.03; eigenvalue 2 = .67

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Eigenvalue 1 = 8.37; eigenvalue 2 = .86


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.29, RMSEA=.063, RMSR=.05, TLI=.883

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.098, RMSR=.06, TLI=.797

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities when no response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.048, RMSR=.04, TLI=.913
#     MR1  MR2
#MR1 1.00 0.58
#MR2 0.58 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.082, RMSR=.05, TLI=.857
#     MR1  MR2
#MR1 1.0 0.6
#MR2 0.6 1.0

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  BDI_1+BDI_2+BDI_3+BDI_4+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_11+BDI_12+
            BDI_13+BDI_14+BDI_15+BDI_16+BDI_17+BDI_18+BDI_19+BDI_20+BDI_21
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.933
#Tucker-Lewis Index (TLI)                       0.978       0.925
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.829
#RMSEA                                          0.051       0.060
#Robust RMSEA                                               0.090
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .417

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.868       0.885
#Tucker-Lewis Index (TLI)                       0.853       0.873
#Robust Comparative Fit Index (CFI)                         0.886
#Robust Tucker-Lewis Index (TLI)                            0.873
#RMSEA                                          0.063       0.045
#Robust RMSEA                                               0.058
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .299

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 negatt  =~  BDI_1+BDI_2+BDI_3+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_14
 perfsom =~  BDI_4+BDI_11+BDI_12+BDI_13+BDI_17+BDI_19+BDI_21+
             BDI_15+BDI_16+BDI_18+BDI_20
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.986       0.949
#Tucker-Lewis Index (TLI)                       0.984       0.943
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.866
#RMSEA                                          0.043       0.052
#Robust RMSEA                                               0.080
#SRMR                                           0.055       0.055

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#negatt ~~                                                             
#  perfsom           0.857    0.016   54.957    0.000    0.857    0.857


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .45

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.895       0.913
#Tucker-Lewis Index (TLI)                       0.883       0.903
#Robust Comparative Fit Index (CFI)                         0.914
#Robust Tucker-Lewis Index (TLI)                            0.904
#RMSEA                                          0.056       0.039
#Robust RMSEA                                               0.050
#SRMR                                           0.044       0.044

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .322

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#negatt ~~                                                             
#  perfsom           0.850    0.021   41.088    0.000    0.850    0.850

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.777       0.791
#Tucker-Lewis Index (TLI)                       0.752       0.768
#Robust Comparative Fit Index (CFI)                         0.792
#Robust Tucker-Lewis Index (TLI)                            0.769
#RMSEA                                          0.082       0.061
#Robust RMSEA                                               0.078
#SRMR                                           0.188       0.188

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .302

# Bifactor model
BIFmodel= '
 negatt  =~  BDI_1+BDI_2+BDI_3+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_14
 perfsom =~  BDI_4+BDI_11+BDI_12+BDI_13+BDI_17+BDI_19+BDI_21+
             BDI_15+BDI_16+BDI_18+BDI_20
 general=~  BDI_1+BDI_2+BDI_3+BDI_4+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_11+BDI_12+
            BDI_13+BDI_14+BDI_15+BDI_16+BDI_17+BDI_18+BDI_19+BDI_20+BDI_21
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.939       0.954
#Tucker-Lewis Index (TLI)                       0.924       0.943
#Robust Comparative Fit Index (CFI)                         0.955
#Robust Tucker-Lewis Index (TLI)                            0.944
#RMSEA                                          0.045       0.030
#Robust RMSEA                                               0.038
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .331

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7936393      0.5238095      0.9011517      0.8363439 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#negatt  0.1121032 0.0583905 0.8878968 0.8422394 0.01029061 0.3174784 0.6580088
#perfsom 0.3088269 0.1479702 0.6911731 0.8195504 0.20644953 0.5977737 0.7814090
#general 0.7936393 0.7936393 0.7936393 0.9011517 0.83634389 0.8941481 0.9408224







################################################################
### Rabeya et al https://zenodo.org/record/5808314


library(readxl)
BDI_II_Rabeya <- read_excel("BDI-II_Rabeya.xlsx")
colnames(BDI_II_Rabeya)
mydata <- as.data.frame(BDI_II_Rabeya[,17:37])
colnames(mydata)
mydata <- na.omit(mydata)
colnames(mydata) <- c("BDI_1","BDI_2","BDI_3","BDI_4","BDI_5","BDI_6","BDI_7","BDI_8","BDI_9","BDI_10",
                      "BDI_11","BDI_12","BDI_13","BDI_14","BDI_15","BDI_16","BDI_17","BDI_18","BDI_19",
                      "BDI_20","BDI_21")
min(mydata)
max(mydata) #4 response alternatives


library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 5.73; eigenvalue 2 = 2.47

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Eigenvalue 1 = 6.56; eigenvalue 2 = 2.97


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.27, RMSEA=.12, RMSR=.13, TLI=.6

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.31, RMSEA=.151, RMSR=.16, TLI=.547

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communities when no response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.055, RMSR=.04, TLI=.917
#     MR1  MR2
#MR1 1.00 0.34
#MR2 0.34 1.00


fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.082, RMSR=.05, TLI=.867
#     MR1  MR2
#MR1 1.00 0.32
#MR2 0.32 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  BDI_1+BDI_2+BDI_3+BDI_4+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_11+BDI_12+
            BDI_13+BDI_14+BDI_15+BDI_16+BDI_17+BDI_18+BDI_19+BDI_20+BDI_21
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.903       0.812
#Tucker-Lewis Index (TLI)                       0.892       0.791
#Robust Comparative Fit Index (CFI)                         0.638
#Robust Tucker-Lewis Index (TLI)                            0.598
#RMSEA                                          0.160       0.145
#Robust RMSEA                                               0.145
#SRMR                                           0.155       0.155

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .214

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.667       0.666
#Tucker-Lewis Index (TLI)                       0.631       0.629
#Robust Comparative Fit Index (CFI)                         0.674
#Robust Tucker-Lewis Index (TLI)                            0.637
#RMSEA                                          0.117       0.107
#Robust RMSEA                                               0.115
#SRMR                                           0.135       0.135

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .096

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 negatt  =~  BDI_1+BDI_2+BDI_3+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_14
 perfsom =~  BDI_4+BDI_11+BDI_12+BDI_13+BDI_17+BDI_19+BDI_21+
             BDI_15+BDI_16+BDI_18+BDI_20
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.917       0.836
#Tucker-Lewis Index (TLI)                       0.908       0.817
#Robust Comparative Fit Index (CFI)                         0.555
#Robust Tucker-Lewis Index (TLI)                            0.503
#RMSEA                                          0.148       0.136
#Robust RMSEA                                               0.161
#SRMR                                           0.145       0.145

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .287

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negatt ~~                                                             
#   perfsom           0.718    0.031   23.403    0.000    0.718    0.718

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.761       0.763
#Tucker-Lewis Index (TLI)                       0.733       0.735
#Robust Comparative Fit Index (CFI)                         0.768
#Robust Tucker-Lewis Index (TLI)                            0.741
#RMSEA                                          0.100       0.091
#Robust RMSEA                                               0.098
#SRMR                                           0.128       0.128

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negatt ~~                                                             
#   perfsom           0.489    0.072    6.815    0.000    0.489    0.489

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .32

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.738       0.740
#Tucker-Lewis Index (TLI)                       0.708       0.711
#Robust Comparative Fit Index (CFI)                         0.745
#Robust Tucker-Lewis Index (TLI)                            0.717
#RMSEA                                          0.104       0.095
#Robust RMSEA                                               0.102
#SRMR                                           0.177       0.177

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .332

# Bifactor model
BIFmodel= '
 negatt  =~  BDI_1+BDI_2+BDI_3+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_14
 perfsom =~  BDI_4+BDI_11+BDI_12+BDI_13+BDI_17+BDI_19+BDI_21+
             BDI_15+BDI_16+BDI_18+BDI_20
 general=~  BDI_1+BDI_2+BDI_3+BDI_4+BDI_5+BDI_6+BDI_7+BDI_8+BDI_9+BDI_10+BDI_11+BDI_12+
            BDI_13+BDI_14+BDI_15+BDI_16+BDI_17+BDI_18+BDI_19+BDI_20+BDI_21
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.892       0.898
#Tucker-Lewis Index (TLI)                       0.865       0.872
#Robust Comparative Fit Index (CFI)                         0.900
#Robust Tucker-Lewis Index (TLI)                            0.875
#RMSEA                                          0.071       0.063
#Robust RMSEA                                               0.068
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .373

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negatt =~                                                             
#    BDI_1             0.235    0.075    3.139    0.002    0.235    0.240
#    BDI_2             0.242    0.093    2.606    0.009    0.242    0.243
#    BDI_3             0.355    0.086    4.146    0.000    0.355    0.337
#    BDI_5             0.552    0.088    6.250    0.000    0.552    0.497
#    BDI_6             0.219    0.080    2.754    0.006    0.219    0.197
#    BDI_7             0.392    0.075    5.209    0.000    0.392    0.359
#    BDI_8             0.493    0.074    6.663    0.000    0.493    0.471
#    BDI_9            -0.399    0.076   -5.277    0.000   -0.399   -0.376
#    BDI_10           -0.457    0.089   -5.156    0.000   -0.457   -0.391
#    BDI_14           -0.280    0.082   -3.406    0.001   -0.280   -0.274
#  perfsom =~                                                            
#    BDI_4            -0.251    0.071   -3.531    0.000   -0.251   -0.239
#    BDI_11            0.365    0.057    6.430    0.000    0.365    0.457
#    BDI_12            0.384    0.066    5.797    0.000    0.384    0.448
#    BDI_13            0.317    0.067    4.743    0.000    0.317    0.348
#    BDI_17            0.381    0.060    6.389    0.000    0.381    0.416
#    BDI_19            0.440    0.055    7.954    0.000    0.440    0.493
#    BDI_21            0.261    0.073    3.571    0.000    0.261    0.291
#    BDI_15            0.341    0.050    6.878    0.000    0.341    0.440
#    BDI_16            0.313    0.049    6.351    0.000    0.313    0.412
#    BDI_18            0.438    0.050    8.748    0.000    0.438    0.515
#    BDI_20            0.396    0.053    7.437    0.000    0.396    0.467
#  general =~                                                            
#    BDI_1             0.745    0.046   16.134    0.000    0.745    0.760
#    BDI_2             0.712    0.055   12.972    0.000    0.712    0.715
#    BDI_3             0.714    0.053   13.512    0.000    0.714    0.677
#    BDI_4             0.789    0.045   17.459    0.000    0.789    0.751
#    BDI_5             0.638    0.072    8.807    0.000    0.638    0.574
#    BDI_6             0.730    0.057   12.730    0.000    0.730    0.656
#    BDI_7             0.821    0.051   16.071    0.000    0.821    0.751
#    BDI_8             0.681    0.063   10.873    0.000    0.681    0.650
#    BDI_9             0.533    0.069    7.668    0.000    0.533    0.502
#    BDI_10            0.439    0.080    5.472    0.000    0.439    0.375
#    BDI_11            0.250    0.052    4.801    0.000    0.250    0.313
#    BDI_12            0.293    0.057    5.135    0.000    0.293    0.342
#    BDI_13            0.185    0.060    3.105    0.002    0.185    0.203
#    BDI_14            0.267    0.067    3.955    0.000    0.267    0.260
#    BDI_15            0.275    0.050    5.479    0.000    0.275    0.355
#    BDI_16            0.223    0.044    5.046    0.000    0.223    0.294
#    BDI_17            0.480    0.054    8.938    0.000    0.480    0.524
#    BDI_18            0.184    0.054    3.389    0.001    0.184    0.217
#    BDI_19            0.321    0.061    5.289    0.000    0.321    0.360
#    BDI_20            0.315    0.058    5.444    0.000    0.315    0.372
#    BDI_21            0.276    0.054    5.106    0.000    0.276    0.307

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6330714      0.5238095      0.9046200      0.7650418 
#$FactorLevelIndices
#        ECV_SS    ECV_SG    ECV_GS    Omega     OmegaH         H        FD
#negatt  0.2475438 0.1429384 0.7524562 0.879985 0.04068044 0.5964849 0.8045514
#perfsom 0.5300622 0.2239902 0.4699378 0.816441 0.40936385 0.7076554 0.8526086
#general 0.6330714 0.6330714 0.6330714 0.904620 0.76504184 0.9070386 0.9458346






