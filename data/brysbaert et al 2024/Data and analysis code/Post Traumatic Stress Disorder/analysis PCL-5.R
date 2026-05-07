################################################################
### Analysis PCL-5
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161645

#First determine overall 2 factor structure for comparison

library(haven)
Ashbaugh_PCL_5 <- read_sav("Ashbaugh_PCL-5.sav")
colnames(Ashbaugh_PCL_5)
mydata  <- as.matrix(Ashbaugh_PCL_5[,4:23])
mydata[mydata > 4] <- NA
mydata <- na.omit(mydata)
colnames(mydata) <- c("PCL1","PCL2","PCL3","PCL4","PCL5","PCL6","PCL7","PCL8","PCL9","PCL10",
                      "PCL11","PCL12","PCL13","PCL14","PCL15","PCL16","PCL17","PCL18","PCL19","PCL20")

Ashbaugh_PCL_5_FR <- read_sav("Ashbaugh_PCL-5.FR.sav")
mydata2  <- as.matrix(Ashbaugh_PCL_5_FR[,4:23])
mydata2 <- na.omit(mydata2)
colnames(mydata2) <- c("PCL1","PCL2","PCL3","PCL4","PCL5","PCL6","PCL7","PCL8","PCL9","PCL10",
                      "PCL11","PCL12","PCL13","PCL14","PCL15","PCL16","PCL17","PCL18","PCL19","PCL20")

Orovou_PCL_5 <- read_sav("Orovou_PCL-5.sav")
colnames(Orovou_PCL_5)
mydata3  <- as.matrix(Orovou_PCL_5[,26:45])
mydata3 <- na.omit(mydata3)
colnames(mydata3) <- c("PCL1","PCL2","PCL3","PCL4","PCL5","PCL6","PCL7","PCL8","PCL9","PCL10",
                      "PCL11","PCL12","PCL13","PCL14","PCL15","PCL16","PCL17","PCL18","PCL19","PCL20")

mydata_all <- rbind(mydata,mydata2,mydata3)

library(psych)
# Give solution with 2 factors
fit1 <- fa(mydata_all,2)
fit1
diagram(fit1)

#factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
#factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9


################################################################
### Ashbaugh et al. 2016 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161645
### English sample see also https://www.tandfonline.com/doi/epdf/10.1080/00223891.2018.1449116 for single factor model

min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 1 component
# Eigenvalue 1 = 10.05; eigenvalue 2 = .69

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 11.76; eigenvalue 2 = .67

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.50, RMSEA=.105, RMSR=.06, TLI=.838

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.129, RMSR=.06, TLI=.819

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities and response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.085, RMSR=.04, TLI=.895
#     MR1  MR2
#MR1 1.00 0.77
#MR2 0.77 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.109, RMSR=.04, TLI=.87
#    MR1 MR2
#MR1   1   0
#MR2   0   1

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ PCL1+PCL2+PCL3+PCL4+PCL5+PCL6+PCL7+PCL8+PCL9+PCL10+PCL11+PCL12+PCL13+
           PCL14+PCL15+PCL16+PCL17+PCL18+PCL19+PCL20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.948
#Tucker-Lewis Index (TLI)                       0.990       0.942
#Robust Comparative Fit Index (CFI)                         0.847
#Robust Tucker-Lewis Index (TLI)                            0.829
#RMSEA                                          0.078       0.100
#Robust RMSEA                                               0.127
#SRMR                                           0.056       0.056

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .621

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.856       0.860
#Tucker-Lewis Index (TLI)                       0.839       0.843
#Robust Comparative Fit Index (CFI)                         0.861
#Robust Tucker-Lewis Index (TLI)                            0.845
#RMSEA                                          0.105       0.086
#Robust RMSEA                                               0.103
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .515

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on EGA as the minimal multidimensional model
EGAmodel= '
 factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
 factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.994       0.959
#Tucker-Lewis Index (TLI)                       0.993       0.954
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.865
#RMSEA                                          0.067       0.089
#Robust RMSEA                                               0.113
#SRMR                                           0.050       0.050

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.900    0.010   90.558    0.000    0.900    0.900

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .641

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.898
#Tucker-Lewis Index (TLI)                       0.880       0.885
#Robust Comparative Fit Index (CFI)                         0.899
#Robust Tucker-Lewis Index (TLI)                            0.887
#RMSEA                                          0.091       0.074
#Robust RMSEA                                               0.088
#SRMR                                           0.049       0.049

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.870    0.015   59.076    0.000    0.870    0.870

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .548


CFA_model3 <- cfa(EGAmodel,mydata,estimator='MLR',std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.815       0.817
#Tucker-Lewis Index (TLI)                       0.793       0.795
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.799
#RMSEA                                          0.119       0.098
#Robust RMSEA                                               0.117
#SRMR                                           0.313       0.313

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .526

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
 factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9
 general=~ PCL1+PCL2+PCL3+PCL4+PCL5+PCL6+PCL7+PCL8+PCL9+PCL10+PCL11+PCL12+PCL13+
           PCL14+PCL15+PCL16+PCL17+PCL18+PCL19+PCL20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.932       0.933
#Tucker-Lewis Index (TLI)                       0.913       0.915
#Robust Comparative Fit Index (CFI)                         0.936
#Robust Tucker-Lewis Index (TLI)                            0.919
#RMSEA                                          0.077       0.063
#Robust RMSEA                                               0.074
#SRMR                                           0.039       0.039

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .575

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    PCL1              0.231    0.101    2.276    0.023    0.231    0.186
#    PCL4              0.095    0.120    0.790    0.429    0.095    0.075
#    PCL3              0.339    0.091    3.707    0.000    0.339    0.293
#    PCL5              0.203    0.084    2.401    0.016    0.203    0.161
#    PCL6             -0.043    0.123   -0.349    0.727   -0.043   -0.033
#    PCL2              0.495    0.057    8.694    0.000    0.495    0.459
#    PCL17             0.169    0.080    2.109    0.035    0.169    0.131
#    PCL7             -0.067    0.097   -0.694    0.488   -0.067   -0.051
#    PCL18             0.224    0.074    3.036    0.002    0.224    0.204
#    PCL20             0.313    0.059    5.307    0.000    0.313    0.252
#    PCL11            -0.369    0.068   -5.397    0.000   -0.369   -0.275
#    PCL10            -0.307    0.066   -4.665    0.000   -0.307   -0.227
#    PCL8              0.077    0.066    1.169    0.243    0.077    0.064
#    PCL19             0.270    0.053    5.124    0.000    0.270    0.217
#  factor2 =~                                                            
#    PCL13             0.629    0.040   15.649    0.000    0.629    0.525
#    PCL12             0.534    0.040   13.275    0.000    0.534    0.481
#    PCL14             0.519    0.046   11.190    0.000    0.519    0.411
#    PCL15             0.366    0.046    7.896    0.000    0.366    0.320
#    PCL16             0.293    0.045    6.582    0.000    0.293    0.291
#    PCL9              0.206    0.058    3.537    0.000    0.206    0.157
#  general =~                                                            
#    PCL1              0.958    0.040   23.916    0.000    0.958    0.774
#    PCL2              0.709    0.052   13.645    0.000    0.709    0.657
#    PCL3              0.796    0.046   17.446    0.000    0.796    0.689
#    PCL4              0.939    0.033   28.099    0.000    0.939    0.741
#    PCL5              0.908    0.040   22.885    0.000    0.908    0.719
#    PCL6              1.004    0.034   29.706    0.000    1.004    0.759
#    PCL7              0.953    0.040   23.676    0.000    0.953    0.722
#    PCL8              0.577    0.044   13.022    0.000    0.577    0.481
#    PCL9              0.974    0.041   23.543    0.000    0.974    0.746
#    PCL10             1.027    0.044   23.483    0.000    1.027    0.760
#    PCL11             1.115    0.036   30.711    0.000    1.115    0.830
#    PCL12             0.771    0.041   18.859    0.000    0.771    0.694
#    PCL13             0.844    0.041   20.747    0.000    0.844    0.704
#    PCL14             0.939    0.042   22.342    0.000    0.939    0.744
#    PCL15             0.792    0.040   19.739    0.000    0.792    0.692
#    PCL16             0.618    0.043   14.429    0.000    0.618    0.613
#    PCL17             0.693    0.044   15.733    0.000    0.693    0.538
#    PCL18             0.683    0.045   15.316    0.000    0.683    0.622
#    PCL19             0.904    0.042   21.678    0.000    0.904    0.725
#    PCL20             0.838    0.042   19.831    0.000    0.838    0.673

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
ECV.general            PUC  Omega.general OmegaH.general 
#0.8626520      0.4421053      0.9583141      0.9252021 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.08904057 0.05891162 0.9109594 0.9364774 0.02069523 0.4287446 0.7878490
#factor2 0.23180433 0.07843640 0.7681957 0.9114921 0.19471626 0.5275291 0.8234452
#general 0.86265199 0.86265199 0.8626520 0.9583141 0.92520213 0.9538537 0.9742858







################################################################
### Ashbaugh et al. 2016 https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0161645
### French sample

mydata <- mydata2
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 1 component
# Eigenvalue 1 = 9.30; eigenvalue 2 = .82

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 11.14; eigenvalue 2 = .82

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.47, RMSEA=.108, RMSR=.07, TLI=.808

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.143, RMSR=.07, TLI=.763

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communities and response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.09, RMSR=.05, TLI=.866
#     MR1  MR2
#MR1 1.00 0.68
#MR2 0.68 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.128, RMSR=.05, TLI=.81
#     MR1  MR2
#MR1 1.00 0.72
#MR2 0.72 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ PCL1+PCL2+PCL3+PCL4+PCL5+PCL6+PCL7+PCL8+PCL9+PCL10+PCL11+PCL12+PCL13+
           PCL14+PCL15+PCL16+PCL17+PCL18+PCL19+PCL20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.952
#Tucker-Lewis Index (TLI)                       0.990       0.946
#Robust Comparative Fit Index (CFI)                         0.825
#Robust Tucker-Lewis Index (TLI)                            0.805
#RMSEA                                          0.072       0.092
#Robust RMSEA                                               0.133
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .551

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.827       0.838
#Tucker-Lewis Index (TLI)                       0.806       0.819
#Robust Comparative Fit Index (CFI)                         0.842
#Robust Tucker-Lewis Index (TLI)                            0.823
#RMSEA                                          0.110       0.089
#Robust RMSEA                                               0.104
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .428

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on EGA as the minimal multidimensional model
EGAmodel= '
 factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
 factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.993       0.959
#Tucker-Lewis Index (TLI)                       0.992       0.953
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.831
#RMSEA                                          0.064       0.086
#Robust RMSEA                                               0.124
#SRMR                                           0.063       0.063

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.909    0.016   57.285    0.000    0.909    0.909

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .577

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.855       0.866
#Tucker-Lewis Index (TLI)                       0.837       0.850
#Robust Comparative Fit Index (CFI)                         0.871
#Robust Tucker-Lewis Index (TLI)                            0.855
#RMSEA                                          0.101       0.081
#Robust RMSEA                                               0.094
#SRMR                                           0.062       0.062


#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.869    0.029   30.344    0.000    0.869    0.869

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .538


CFA_model3 <- cfa(EGAmodel,mydata,estimator='MLR',std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.784       0.789
#Tucker-Lewis Index (TLI)                       0.755       0.762
#Robust Comparative Fit Index (CFI)                         0.797
#Robust Tucker-Lewis Index (TLI)                            0.770
#RMSEA                                          0.125       0.104
#Robust RMSEA                                               0.119
#SRMR                                           0.288       0.288

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .504

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
 factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9
 general=~ PCL1+PCL2+PCL3+PCL4+PCL5+PCL6+PCL7+PCL8+PCL9+PCL10+PCL11+PCL12+PCL13+
           PCL14+PCL15+PCL16+PCL17+PCL18+PCL19+PCL20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.919       0.927
#Tucker-Lewis Index (TLI)                       0.898       0.907
#Robust Comparative Fit Index (CFI)                         0.932
#Robust Tucker-Lewis Index (TLI)                            0.913
#RMSEA                                          0.080       0.064
#Robust RMSEA                                               0.073
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .515

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    PCL1              0.436    0.111    3.926    0.000    0.436    0.370
#    PCL4              0.566    0.125    4.514    0.000    0.566    0.445
#    PCL3              0.051    0.136    0.374    0.709    0.051    0.044
#    PCL5              0.229    0.145    1.584    0.113    0.229    0.183
#    PCL6              0.779    0.142    5.475    0.000    0.779    0.585
#    PCL2              0.042    0.080    0.527    0.598    0.042    0.044
#    PCL17            -0.132    0.172   -0.768    0.443   -0.132   -0.110
#    PCL7              0.600    0.172    3.482    0.000    0.600    0.441
#    PCL18            -0.140    0.188   -0.743    0.457   -0.140   -0.118
#    PCL20             0.123    0.088    1.392    0.164    0.123    0.105
#    PCL11             0.263    0.082    3.220    0.001    0.263    0.198
#    PCL10             0.121    0.092    1.315    0.188    0.121    0.098
#    PCL8              0.280    0.094    2.995    0.003    0.280    0.226
#    PCL19            -0.134    0.123   -1.088    0.277   -0.134   -0.109
#  factor2 =~                                                            
#    PCL13             0.894    0.123    7.293    0.000    0.894    0.711
#    PCL12             0.335    0.096    3.504    0.000    0.335    0.291
#    PCL14             0.477    0.122    3.901    0.000    0.477    0.408
#    PCL15             0.231    0.101    2.281    0.023    0.231    0.200
#    PCL16             0.134    0.078    1.716    0.086    0.134    0.163
#    PCL9              0.100    0.114    0.881    0.379    0.100    0.082
#  general =~                                                            
#    PCL1              0.799    0.072   11.048    0.000    0.799    0.678
#    PCL2              0.628    0.063   10.005    0.000    0.628    0.656
#    PCL3              0.753    0.071   10.527    0.000    0.753    0.651
#    PCL4              0.774    0.082    9.404    0.000    0.774    0.608
#    PCL5              0.870    0.072   12.082    0.000    0.870    0.694
#    PCL6              0.770    0.104    7.400    0.000    0.770    0.577
#    PCL7              0.827    0.102    8.102    0.000    0.827    0.607
#    PCL8              0.505    0.090    5.621    0.000    0.505    0.408
#    PCL9              0.838    0.077   10.940    0.000    0.838    0.687
#    PCL10             0.729    0.078    9.346    0.000    0.729    0.593
#    PCL11             1.001    0.056   17.876    0.000    1.001    0.754
#    PCL12             0.813    0.080   10.122    0.000    0.813    0.705
#    PCL13             0.918    0.084   10.938    0.000    0.918    0.730
#    PCL14             0.677    0.085    7.976    0.000    0.677    0.580
#    PCL15             0.869    0.074   11.707    0.000    0.869    0.752
#    PCL16             0.480    0.077    6.235    0.000    0.480    0.585
#    PCL17             0.841    0.074   11.430    0.000    0.841    0.700
#    PCL18             0.973    0.072   13.481    0.000    0.973    0.819
#    PCL19             1.013    0.064   15.793    0.000    1.013    0.822
#    PCL20             0.747    0.073   10.245    0.000    0.747    0.638

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general  Omega.general OmegaH.general 
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.8257625      0.4421053      0.9526550      0.9051307 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1457174 0.09758072 0.8542826 0.9305706 0.05931669 0.5757852 0.8243161
#factor2 0.2320523 0.07665680 0.7679477 0.8907470 0.15514108 0.5816774 1.0105095
#general 0.8257625 0.82576248 0.8257625 0.9526550 0.90513067 0.9476615 0.9681120






################################################################
### Orovou et al. 2021 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8363016/#sec017
### 

mydata <- mydata3
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 1 component
# Eigenvalue 1 = 12.0; eigenvalue 2 = .50

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 15.23; eigenvalue 2 = .36

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.60, RMSEA=.116, RMSR=.05, TLI=.851

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.76, RMSEA=.159, RMSR=.04, TLI=.833

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 communities and response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.106, RMSR=.04, TLI=.877
#     MR1  MR2
#MR1 1.00 0.77
#MR2 0.77 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.78, RMSEA=.149, RMSR=.03, TLI=.853
#     MR1  MR2
#MR1 1.00 0.86
#MR2 0.86 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ PCL1+PCL2+PCL3+PCL4+PCL5+PCL6+PCL7+PCL8+PCL9+PCL10+PCL11+PCL12+PCL13+
           PCL14+PCL15+PCL16+PCL17+PCL18+PCL19+PCL20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
# Comparative Fit Index (CFI)                    0.998       0.976
#Tucker-Lewis Index (TLI)                       0.998       0.973
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.866
#RMSEA                                          0.057       0.081
#Robust RMSEA                                               0.144
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .783

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.866       0.900
#Tucker-Lewis Index (TLI)                       0.851       0.888
#Robust Comparative Fit Index (CFI)                         0.901
#Robust Tucker-Lewis Index (TLI)                            0.890
#RMSEA                                          0.117       0.058
#Robust RMSEA                                               0.099
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .617

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on EGA as the minimal multidimensional model
EGAmodel= '
 factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
 factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.998       0.980
#Tucker-Lewis Index (TLI)                       0.998       0.977
#Robust Comparative Fit Index (CFI)                         0.896
#Robust Tucker-Lewis Index (TLI)                            0.883
#RMSEA                                          0.050       0.075
#Robust RMSEA                                               0.135
#SRMR                                           0.036       0.036

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.946    0.009  109.123    0.000    0.946    0.946

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .796

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.918
#Tucker-Lewis Index (TLI)                       0.869       0.908
#Robust Comparative Fit Index (CFI)                         0.919
#Robust Tucker-Lewis Index (TLI)                            0.909
#RMSEA                                          0.110       0.053
#Robust RMSEA                                               0.089
#SRMR                                           0.047       0.047


#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.911    0.020   45.159    0.000    0.911    0.911

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .644


CFA_model3 <- cfa(EGAmodel,mydata,estimator='MLR',std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.812       0.841
#Tucker-Lewis Index (TLI)                       0.790       0.822
#Robust Comparative Fit Index (CFI)                         0.844
#Robust Tucker-Lewis Index (TLI)                            0.826
#RMSEA                                          0.139       0.073
#Robust RMSEA                                               0.124
#SRMR                                           0.363       0.363

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .627

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ PCL1+PCL4+PCL3+PCL5+PCL6+PCL2+PCL17+PCL7+PCL18+PCL20+PCL11+PCL10+PCL8+PCL19
 factor2 =~ PCL13+PCL12+PCL14+PCL15+PCL16+PCL9
 general=~ PCL1+PCL2+PCL3+PCL4+PCL5+PCL6+PCL7+PCL8+PCL9+PCL10+PCL11+PCL12+PCL13+
           PCL14+PCL15+PCL16+PCL17+PCL18+PCL19+PCL20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.923       0.946
#Tucker-Lewis Index (TLI)                       0.903       0.932
#Robust Comparative Fit Index (CFI)                         0.952
#Robust Tucker-Lewis Index (TLI)                            0.939
#RMSEA                                          0.095       0.045
#Robust RMSEA                                               0.074
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .654

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    PCL1              0.378    0.126    2.993    0.003    0.378    0.328
#    PCL4              0.303    0.137    2.208    0.027    0.303    0.289
#    PCL3              0.323    0.154    2.104    0.035    0.323    0.313
#    PCL5              0.169    0.123    1.377    0.169    0.169    0.174
#    PCL6             -0.012    0.071   -0.166    0.869   -0.012   -0.012
#    PCL2              0.106    0.102    1.040    0.298    0.106    0.116
#    PCL17            -0.116    0.078   -1.482    0.138   -0.116   -0.124
#    PCL7              0.116    0.098    1.192    0.233    0.116    0.113
#    PCL18            -0.032    0.090   -0.349    0.727   -0.032   -0.037
#    PCL20            -0.188    0.160   -1.173    0.241   -0.188   -0.258
#    PCL11             0.009    0.108    0.085    0.932    0.009    0.011
#    PCL10            -0.050    0.121   -0.413    0.680   -0.050   -0.050
#    PCL8             -0.008    0.093   -0.091    0.928   -0.008   -0.009
#    PCL19            -0.174    0.181   -0.961    0.337   -0.174   -0.200
#  factor2 =~                                                            
#    PCL13             0.540    0.078    6.911    0.000    0.540    0.556
#    PCL12             0.285    0.070    4.065    0.000    0.285    0.310
#    PCL14             0.491    0.082    6.009    0.000    0.491    0.513
#    PCL15             0.237    0.069    3.433    0.001    0.237    0.219
#    PCL16             0.025    0.050    0.503    0.615    0.025    0.039
#    PCL9              0.122    0.078    1.571    0.116    0.122    0.118
#  general =~                                                            
#    PCL1              0.928    0.062   14.955    0.000    0.928    0.805
#    PCL2              0.679    0.065   10.453    0.000    0.679    0.744
#    PCL3              0.852    0.062   13.733    0.000    0.852    0.824
#    PCL4              0.892    0.062   14.381    0.000    0.892    0.852
#    PCL5              0.811    0.058   14.002    0.000    0.811    0.834
#    PCL6              0.784    0.060   13.093    0.000    0.784    0.809
#    PCL7              0.797    0.059   13.545    0.000    0.797    0.778
#    PCL8              0.641    0.058   11.109    0.000    0.641    0.714
#    PCL9              0.739    0.062   11.986    0.000    0.739    0.713
#    PCL10             0.756    0.063   11.976    0.000    0.756    0.756
#    PCL11             0.679    0.057   11.941    0.000    0.679    0.781
#    PCL12             0.707    0.056   12.716    0.000    0.707    0.769
#    PCL13             0.702    0.054   13.064    0.000    0.702    0.722
#    PCL14             0.661    0.060   11.061    0.000    0.661    0.691
#    PCL15             0.840    0.057   14.800    0.000    0.840    0.778
#    PCL16             0.389    0.066    5.882    0.000    0.389    0.606
#    PCL17             0.778    0.062   12.505    0.000    0.778    0.835
#    PCL18             0.691    0.062   11.117    0.000    0.691    0.806
#    PCL19             0.704    0.057   12.276    0.000    0.704    0.810
#    PCL20             0.538    0.056    9.621    0.000    0.538    0.739

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.9079641      0.4421053      0.9719222      0.9577021 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#factor1 0.05090099 0.03609700 0.9490990 0.9630979 0.003323899 0.3382801 0.7854154
#factor2 0.19233638 0.05593891 0.8076636 0.9068873 0.130572803 0.4939978 0.8488000
#general 0.90796409 0.90796409 0.9079641 0.9719222 0.957702080 0.9691045 0.9835181







