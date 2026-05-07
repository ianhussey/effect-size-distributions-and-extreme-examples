################################################################
### Brief Self control scale
### 
### Paap et al. (2023) https://link.springer.com/article/10.3758/s13428-023-02089-2
### data: https://osf.io/35y6k/?vie


####################################################################
### First determine factor solution on the basis of all datasets
###

library(haven)
BSC1_Paap <- read_sav("BSC1_Paap.sav")
colnames(BSC1_Paap)
mydata1 <- BSC1_Paap[,63:75]
mydata1 <- na.omit(mydata1)
colnames(mydata1)
min(mydata1)
max(mydata1) # response alternatives = 5

mydata2 <- BSC1_Paap[,79:91]
mydata2 <- na.omit(mydata2)
colnames(mydata2) <- colnames(mydata1)
colnames(mydata2)
min(mydata2)
max(mydata2) # response alternatives = 5

mydata3 <- BSC1_Paap[,94:106]
mydata3 <- na.omit(mydata3)
colnames(mydata3) <- colnames(mydata1)
colnames(mydata3)
min(mydata3)
max(mydata3) # response alternatives = 5

mydata4 <- BSC1_Paap[,109:121]
mydata4 <- na.omit(mydata4)
colnames(mydata4) <- colnames(mydata1)
colnames(mydata4)
min(mydata4)
max(mydata4) # response alternatives = 5

mydata5 <- BSC1_Paap[,124:136]
mydata5 <- na.omit(mydata5)
colnames(mydata5) <- colnames(mydata1)
colnames(mydata5)
min(mydata5)
max(mydata5) # response alternatives = 5

library(haven)
BSC2_Paap <- read_sav("BSC2_Paap.sav")
colnames(BSC2_Paap)
mydata6 <- BSC2_Paap[,78:90]
mydata6 <- na.omit(mydata6)
colnames(mydata6) <- colnames(mydata1)
colnames(mydata6)
min(mydata6)
max(mydata6) # response alternatives = 5

mydata <- rbind(mydata1,mydata2,mydata3,mydata4,mydata5,mydata6)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 2 components
# Eigenvalue 1 = 4.71; eigenvalue 2 = .6

# Give solution with 2 factors (general EFA-based)
fit2 <- fa(mydata,2)
fit2
diagram(fit2)

# Two components largely agree with reverse scoring; the latter distinction used
# also because EFA model otherwise does not always work

EGAmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
'

################################################################
### 
### Study 1 Group 1 Original
### 

mydata <- mydata1

library(psych)
omega(mydata) # alpha = .42, omega T = .72

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 4.74; eigenvalue 2 = 1.17

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 5.3; eigenvalue 2 = 1.28

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.36, RMSEA=.135, RMSR=.12, TLI=.736

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.166, RMSR=.13, TLI=.696

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
# %variance explained=.49, RMSEA=.088, RMSR=.05, TLI=.889
#      MR1   MR2
#MR1  1.00  -0.06
#MR2  -0.06  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.121, RMSR=.05, TLI=.838
#      MR1   MR2
#MR1  1.00  -0.08
#MR2  -0.08  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.916       0.741
#Tucker-Lewis Index (TLI)                       0.898       0.684
#Robust Comparative Fit Index (CFI)                         0.709
#Robust Tucker-Lewis Index (TLI)                            0.644
#RMSEA                                          0.194       0.229
#Robust RMSEA                                               0.182
#SRMR                                           0.128       0.128

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .496

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.742       0.734
#Tucker-Lewis Index (TLI)                       0.685       0.674
#Robust Comparative Fit Index (CFI)                         0.748
#Robust Tucker-Lewis Index (TLI)                            0.692
#RMSEA                                          0.148       0.133
#Robust RMSEA                                               0.146
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .433

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.928
#Tucker-Lewis Index (TLI)                       0.972       0.910
#Robust Comparative Fit Index (CFI)                         0.830
#Robust Tucker-Lewis Index (TLI)                            0.788
#RMSEA                                          0.102       0.122
#Robust RMSEA                                               0.141
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .526

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.239    0.061   -3.906    0.000   -0.239   -0.239

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.868       0.865
#Tucker-Lewis Index (TLI)                       0.835       0.832
#Robust Comparative Fit Index (CFI)                         0.875
#Robust Tucker-Lewis Index (TLI)                            0.844
#RMSEA                                          0.107       0.095
#Robust RMSEA                                               0.104
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .453

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.126    0.121   -1.045    0.296   -0.126   -0.126

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.867       0.867
#Tucker-Lewis Index (TLI)                       0.837       0.838
#Robust Comparative Fit Index (CFI)                         0.875
#Robust Tucker-Lewis Index (TLI)                            0.847
#RMSEA                                          0.107       0.094
#Robust RMSEA                                               0.103
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .452


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.913       0.901
#Tucker-Lewis Index (TLI)                       0.864       0.845
#Robust Comparative Fit Index (CFI)                         0.915
#Robust Tucker-Lewis Index (TLI)                            0.867
#RMSEA                                          0.098       0.092
#Robust RMSEA                                               0.096
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .479

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4659237      0.4848485      0.8383278      0.5015106 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4184610 0.3199520 0.58153897 0.9149856 0.3113558 1.4643683 1.2315391
#factor2 0.9095883 0.2141243 0.09041175 0.2007409 0.1174279 0.7569218 0.8734075
#general 0.4659237 0.4659237 0.46592373 0.8383278 0.5015106 0.8352627 0.9623136





################################################################
### 
### Study 1 Group 2 All pos
### 

### Not included because many analyses did not work well

mydata <- mydata2

library(psych)
omega(mydata) # alpha = .84, omega T = .87

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 3.98; eigenvalue 2 = .67

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 4.63; eigenvalue 2 = .76

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.31, RMSEA=.079, RMSR=.07, TLI=.854

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.104, RMSR=.08, TLI=.81

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
# %variance explained=.37, RMSEA=.049, RMSR=.05, TLI=.942
#      MR1   MR2
#MR1  1.00  0.62
#MR2  0.62  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.073, RMSR=.05, TLI=.905
#      MR1   MR2
#MR1  1.00  0.64
#MR2  0.64  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.973       0.906
#Tucker-Lewis Index (TLI)                       0.967       0.885
#Robust Comparative Fit Index (CFI)                         0.840
#Robust Tucker-Lewis Index (TLI)                            0.805
#RMSEA                                          0.094       0.123
#Robust RMSEA                                               0.116
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .423

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.877       0.871
#Tucker-Lewis Index (TLI)                       0.849       0.843
#Robust Comparative Fit Index (CFI)                         0.882
#Robust Tucker-Lewis Index (TLI)                            0.856
#RMSEA                                          0.087       0.081
#Robust RMSEA                                               0.085
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .343

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.974       0.910
#Tucker-Lewis Index (TLI)                       0.968       0.888
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.813
#RMSEA                                          0.093       0.121
#Robust RMSEA                                               0.113
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .401

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           1.076    0.023   46.710    0.000    1.076    1.076

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.883       0.878
#Tucker-Lewis Index (TLI)                       0.855       0.848
#Robust Comparative Fit Index (CFI)                         0.888
#Robust Tucker-Lewis Index (TLI)                            0.861
#RMSEA                                          0.086       0.080
#Robust RMSEA                                               0.083
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .338

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           1.076    0.034   31.843    0.000    1.076    1.076

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.569       0.546
#Tucker-Lewis Index (TLI)                       0.474       0.445
#Robust Comparative Fit Index (CFI)                         0.571
#Robust Tucker-Lewis Index (TLI)                            0.476
#RMSEA                                          0.163       0.152
#Robust RMSEA                                               0.162
#SRMR                                           0.234       0.234

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .312


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.928       0.879
#Tucker-Lewis Index (TLI)                       0.887       0.810
#Robust Comparative Fit Index (CFI)                         0.917
#Robust Tucker-Lewis Index (TLI)                            0.869
#RMSEA                                          0.076       0.089
#Robust RMSEA                                               0.081
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .407

#Latent Variables:
#                   Estimate  Std.Err  z-value   P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                             
#    BSC10R            0.166    0.294     0.566    0.572    0.166    0.176
#    BSC5R            -0.033    0.084    -0.386    0.699   -0.033   -0.032
#    BSC3R             0.359    0.560     0.642    0.521    0.359    0.287
#    BSC13R            0.309    0.542     0.571    0.568    0.309    0.293
#    BSC4R             0.613    0.822     0.746    0.456    0.613    0.503
#    BSC7R            -0.082    0.115    -0.710    0.478   -0.082   -0.088
#    BSC9R             0.450    0.793     0.568    0.570    0.450    0.376
#    BSC2R            -0.131    0.106    -1.236    0.216   -0.131   -0.117
#  factor2 =~                                                             
#    BSC1              0.003    0.002     1.768    0.077    0.003    0.003
#    BSC6              0.006    0.002     2.628    0.009    0.006    0.005
#    BSC11           -26.454    0.007 -3621.276    0.000  -26.454  -27.899
#    BSC8              0.004    0.002     2.109    0.035    0.004    0.004
#  general =~                                                             
#    BSC10R            0.476    0.061     7.778    0.000    0.476    0.505
#    BSC5R             0.657    0.069     9.468    0.000    0.657    0.643
#    BSC3R             0.598    0.084     7.152    0.000    0.598    0.478
#    BSC13R            0.257    0.092     2.800    0.005    0.257    0.244
#    BSC4R             0.496    0.096     5.186    0.000    0.496    0.407
#    BSC7R             0.681    0.059    11.512    0.000    0.681    0.730
#    BSC9R             0.305    0.098     3.112    0.002    0.305    0.255
#    BSC2R             0.717    0.071    10.144    0.000    0.717    0.643
#    BSC8              0.888    0.062    14.350    0.000    0.888    0.794
#    BSC1              0.603    0.061     9.832    0.000    0.603    0.638
#    BSC6              0.679    0.073     9.274    0.000    0.679    0.631
#    BSC11             0.588    0.072     8.139    0.000    0.588    0.621

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#    0.00504973     0.48484848    15.80354979     0.83332879 
#$FactorLevelIndices
#        ECV_SS       ECV_SG      ECV_GS     Omega      OmegaH           H         FD
#factor1 0.22410665 0.0007871512 0.775893355  0.766168  0.08708032   0.4259842  0.6856690
#factor2 0.99766732 0.9941631186 0.002332679 90.498642 89.66866227 813.2716431 33.5006177
#general 0.00504973 0.0050497302 0.005049730 15.803550  0.83332879   0.8780505  0.9482232





################################################################
### 
### Study 1 Group 3 All neg
### 

### Not included because many analyses did not work well

mydata <- mydata3

library(psych)
omega(mydata) # alpha = .44, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 7.12; eigenvalue 2 = .43

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 1 components
# Eigenvalue 1 = 7.84 eigenvalue 2 = .42

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.55, RMSEA=.089, RMSR=.05, TLI=.922

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.6, RMSEA=.111, RMSR=.05, TLI=.901

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 1 communities

# Give solution with 2 factors (EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.061, RMSR=.03, TLI=.963
#      MR1   MR2
#MR1  1.00  0.79
#MR2  0.79  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.64, RMSEA=.087, RMSR=.03, TLI=.938
#      MR1   MR2
#MR1  1.00  0.82
#MR2  0.82  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.968
#Tucker-Lewis Index (TLI)                       0.996       0.961
#Robust Comparative Fit Index (CFI)                         0.913
#Robust Tucker-Lewis Index (TLI)                            0.894
#RMSEA                                          0.064       0.109
#Robust RMSEA                                               0.119
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .59

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.930       0.940
#Tucker-Lewis Index (TLI)                       0.913       0.926
#Robust Comparative Fit Index (CFI)                         0.942
#Robust Tucker-Lewis Index (TLI)                            0.928
#RMSEA                                          0.096       0.073
#Robust RMSEA                                               0.086
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .518

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           1.024    0.019   53.469    0.000    1.024    1.024

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.472       0.107
#Tucker-Lewis Index (TLI)                       0.354      -0.092
#Robust Comparative Fit Index (CFI)                         0.717
#Robust Tucker-Lewis Index (TLI)                            0.654
#RMSEA                                          0.778       0.578
#Robust RMSEA                                               0.214
#SRMR                                           0.382       0.382

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .586


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
# warning
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.947       0.934
#Tucker-Lewis Index (TLI)                       0.917       0.896
#Robust Comparative Fit Index (CFI)                         0.949
#Robust Tucker-Lewis Index (TLI)                            0.920
#RMSEA                                          0.093       0.086
#Robust RMSEA                                               0.091
#SRMR                                           0.041       0.041

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .542

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    BSC10R           -0.265    0.297   -0.892    0.373   -0.265   -0.206
#    BSC5R             0.112    0.354    0.318    0.751    0.112    0.092
#    BSC3R            -0.023    0.378   -0.062    0.951   -0.023   -0.018
#    BSC13R           -0.264    0.309   -0.855    0.392   -0.264   -0.205
#    BSC4R             0.213    0.299    0.713    0.476    0.213    0.169
#    BSC7R             0.051    0.214    0.239    0.811    0.051    0.040
#    BSC9R            -0.298    0.268   -1.112    0.266   -0.298   -0.221
#    BSC2R             0.414    0.408    1.014    0.311    0.414    0.334
#  factor2 =~                                                            
#    BSC1             21.013    0.017 1232.481    0.000   21.013   17.129
#    BSC6              0.004    0.004    1.044    0.296    0.004    0.003
#    BSC11            -0.007    0.003   -2.660    0.008   -0.007   -0.005
#    BSC8             -0.003    0.003   -0.951    0.342   -0.003   -0.002
#  general =~                                                            
#    BSC10R            0.993    0.059   16.868    0.000    0.993    0.771
#    BSC5R             0.894    0.069   12.900    0.000    0.894    0.731
#    BSC3R             0.943    0.068   13.765    0.000    0.943    0.711
#    BSC13R            1.026    0.068   14.993    0.000    1.026    0.797
#    BSC4R             0.863    0.075   11.524    0.000    0.863    0.685
#    BSC7R             0.885    0.069   12.868    0.000    0.885    0.687
#    BSC9R             1.010    0.068   14.871    0.000    1.010    0.750
#    BSC2R             0.883    0.081   10.936    0.000    0.883    0.713
#    BSC8              0.917    0.070   13.049    0.000    0.917    0.702
#    BSC1              0.949    0.068   13.902    0.000    0.949    0.773
#    BSC6              0.774    0.073   10.632    0.000    0.774    0.649
#    BSC11             1.024    0.062   16.562    0.000    1.024    0.797

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#    0.02141916     0.48484848     4.51444715     0.93733448 
#$FactorLevelIndices
#        ECV_SS       ECV_SG      ECV_GS      Omega       OmegaH           H         FD
#factor1 0.06221022 0.0009463487 0.937789779  0.9086475 5.264709e-06   0.2338549  0.6520417
#factor2 0.99273610 0.9776344881 0.007263902 29.4847538 2.865109e+01 297.2783842 26.0156172
#general 0.02141916 0.0214191632 0.021419163  4.5144471 9.373345e-01   0.9350435  0.9706784






################################################################
### 
### Study 1 Group 4 Who
### 

mydata <- mydata4

library(psych)
omega(mydata) # alpha = .82, omega T = .88

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 4.19; eigenvalue 2 = 1.34

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 4.62 eigenvalue 2 = 1.46

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.135, RMSR=.13, TLI=.687

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.158, RMSR=.14, TLI=.653

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.07, RMSR=.04, TLI=.917
#      MR1   MR2
#MR1  1.00  -.03
#MR2  -.03  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.095, RMSR=.05, TLI=.876
#      MR1   MR2
#MR1  1.00  -.04
#MR2  -.04  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.872       0.649
#Tucker-Lewis Index (TLI)                       0.844       0.571
#Robust Comparative Fit Index (CFI)                         0.698
#Robust Tucker-Lewis Index (TLI)                            0.630
#RMSEA                                          0.195       0.223
#Robust RMSEA                                               0.169
#SRMR                                           0.140       0.140

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .419

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.720       0.718
#Tucker-Lewis Index (TLI)                       0.658       0.655
#Robust Comparative Fit Index (CFI)                         0.724
#Robust Tucker-Lewis Index (TLI)                            0.663
#RMSEA                                          0.145       0.134
#Robust RMSEA                                               0.143
#SRMR                                           0.128       0.128

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .369

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.936
#Tucker-Lewis Index (TLI)                       0.975       0.920
#Robust Comparative Fit Index (CFI)                         0.896
#Robust Tucker-Lewis Index (TLI)                            0.870
#RMSEA                                          0.079       0.096
#Robust RMSEA                                               0.100
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .554

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.128    0.070   -1.843    0.065   -0.128   -0.128

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.920       0.928
#Tucker-Lewis Index (TLI)                       0.901       0.910
#Robust Comparative Fit Index (CFI)                         0.928
#Robust Tucker-Lewis Index (TLI)                            0.911
#RMSEA                                          0.078       0.068
#Robust RMSEA                                               0.074
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .498

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.072    0.098   -0.730    0.465   -0.072   -0.072

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.921       0.929
#Tucker-Lewis Index (TLI)                       0.903       0.913
#Robust Comparative Fit Index (CFI)                         0.929
#Robust Tucker-Lewis Index (TLI)                            0.913
#RMSEA                                          0.077       0.067
#Robust RMSEA                                               0.073
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .498


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.957       0.956
#Tucker-Lewis Index (TLI)                       0.933       0.931
#Robust Comparative Fit Index (CFI)                         0.960
#Robust Tucker-Lewis Index (TLI)                            0.937
#RMSEA                                          0.064       0.060
#Robust RMSEA                                               0.062
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .524


library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.12872166     0.48484848     0.83693271     0.02151193 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#factor1 0.8498500 0.6046449 0.15015003 0.8885083 0.86803537 0.8749786 0.9505753
#factor2 0.9241180 0.2666334 0.07588197 0.3005975 0.29940577 0.7533131 0.8713123
#general 0.1287217 0.1287217 0.12872166 0.8369327 0.02151193 0.4808276 0.8411411






################################################################
### 
### Study 1 Group 5 Mirror
### 

mydata <- mydata5

library(psych)
omega(mydata) # alpha = .69, omega T = .76

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 2.33; eigenvalue 2 = 1.23

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 2.77 eigenvalue 2 = 1.34

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.118, RMSR=.12, TLI=.484

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.21, RMSEA=.14, RMSR=.13, TLI=.476

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.062, RMSR=.05, TLI=.858
#      MR1   MR2
#MR1  1.00  -.13
#MR2  -.13  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.088, RMSR=.06, TLI=.792
#      MR1   MR2
#MR1  1.00  -.14
#MR2  -.14  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.711       0.516
#Tucker-Lewis Index (TLI)                       0.646       0.409
#Robust Comparative Fit Index (CFI)                         0.561
#Robust Tucker-Lewis Index (TLI)                            0.464
#RMSEA                                          0.186       0.186
#Robust RMSEA                                               0.154
#SRMR                                           0.132       0.132

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .166

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.565       0.542
#Tucker-Lewis Index (TLI)                       0.469       0.440
#Robust Comparative Fit Index (CFI)                         0.569
#Robust Tucker-Lewis Index (TLI)                            0.474
#RMSEA                                          0.130       0.118
#Robust RMSEA                                               0.127
#SRMR                                           0.119       0.119

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .153

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Robust Comparative Fit Index (CFI)                         0.600
#Robust Tucker-Lewis Index (TLI)                            0.502
#RMSEA                                          0.164       0.164
#Robust RMSEA                                               0.148
#SRMR                                           0.127       0.127

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .36

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.411    0.063   -6.489    0.000   -0.411   -0.411

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.705       0.688
#Tucker-Lewis Index (TLI)                       0.632       0.611
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.640
#RMSEA                                          0.108       0.098
#Robust RMSEA                                               0.105
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .302

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2          -0.138    0.169   -0.814    0.415   -0.138   -0.138

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.703       0.699
#Tucker-Lewis Index (TLI)                       0.637       0.632
#Robust Comparative Fit Index (CFI)                         0.713
#Robust Tucker-Lewis Index (TLI)                            0.649
#RMSEA                                          0.107       0.096
#Robust RMSEA                                               0.104
#SRMR                                           0.118       0.118

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .296


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.891       0.874
#Tucker-Lewis Index (TLI)                       0.829       0.803
#Robust Comparative Fit Index (CFI)                         0.893
#Robust Tucker-Lewis Index (TLI)                            0.832
#RMSEA                                          0.073       0.070
#Robust RMSEA                                               0.072
#SRMR                                           0.070       0.070

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .385

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4238879      0.4848485      0.6764904      0.2885544 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5558501 0.3058753 0.4441499 0.7207763 0.4415196 0.6222740 0.7711973
#factor2 0.6009051 0.2702368 0.3990949 0.1684107 0.1432016 0.6489429 0.8089062
#general 0.4238879 0.4238879 0.4238879 0.6764904 0.2885544 0.7327607 0.8390490





################################################################
### 
### Study 2 Condition original
### 

mydata <- mydata6

library(psych)
omega(mydata) # alpha = .89, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 5.16; eigenvalue 2 = .84

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 5.83 eigenvalue 2 = .87

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.40, RMSEA=.095, RMSR=.08, TLI=.86

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.114, RMSR=.08, TLI=.841

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.045, RMSR=.03, TLI=.969
#      MR1   MR2
#MR1  1.00  0.46
#MR2  0.46  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.065, RMSR=.04, TLI=.948
#      MR1   MR2
#MR1  1.00  0.48
#MR2  0.48  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.928
#Tucker-Lewis Index (TLI)                       0.971       0.912
#Robust Comparative Fit Index (CFI)                         0.867
#Robust Tucker-Lewis Index (TLI)                            0.837
#RMSEA                                          0.112       0.135
#Robust RMSEA                                               0.122
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .497

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.879       0.889
#Tucker-Lewis Index (TLI)                       0.852       0.864
#Robust Comparative Fit Index (CFI)                         0.893
#Robust Tucker-Lewis Index (TLI)                            0.869
#RMSEA                                          0.103       0.082
#Robust RMSEA                                               0.096
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .438

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.961
#Tucker-Lewis Index (TLI)                       0.987       0.951
#Robust Comparative Fit Index (CFI)                         0.924
#Robust Tucker-Lewis Index (TLI)                            0.905
#RMSEA                                          0.076       0.101
#Robust RMSEA                                               0.093
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .528

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.700    0.038   18.233    0.000    0.700    0.700

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.936       0.950
#Tucker-Lewis Index (TLI)                       0.920       0.937
#Robust Comparative Fit Index (CFI)                         0.951
#Robust Tucker-Lewis Index (TLI)                            0.940
#RMSEA                                          0.076       0.056
#Robust RMSEA                                               0.065
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .478

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#   factor2           0.626    0.085    7.401    0.000    0.626    0.626

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.878       0.886
#Tucker-Lewis Index (TLI)                       0.851       0.860
#Robust Comparative Fit Index (CFI)                         0.891
#Robust Tucker-Lewis Index (TLI)                            0.867
#RMSEA                                          0.103       0.083
#Robust RMSEA                                               0.096
#SRMR                                           0.187       0.187

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .482


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.982       0.991
#Tucker-Lewis Index (TLI)                       0.971       0.986
#Robust Comparative Fit Index (CFI)                         0.992
#Robust Tucker-Lewis Index (TLI)                            0.987
#RMSEA                                          0.045       0.026
#Robust RMSEA                                               0.030
#SRMR                                           0.039       0.039

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .522

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6951871      0.4848485      0.9046465      0.7870211 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2457292 0.1851440 0.7542708 0.9038158 0.1293225 0.5805595 0.7918908
#factor2 0.4853686 0.1196688 0.5146314 0.6968133 0.3280787 0.4833346 0.7157781
#general 0.6951871 0.6951871 0.6951871 0.9046465 0.7870211 0.8799265 0.9268285






################################################################
### 
### Study 1 Sjastad (2024) https://osf.io/w24eb/download
### 
### Data https://osf.io/karuw/?view_only=81f20dfeacdc4e90aa83dee6e1738be5
###

library(readxl)
BSC_Sjastad1 <- read_excel("BSC_Sjastad1.xlsx")
colnames(BSC_Sjastad1)
mydata <- BSC_Sjastad1[,14:26]
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("BSC1","BSC2R","BSC3R","BSC4R","BSC5R","BSC6","BSC7R","BSC8","BSC9R","BSC10R","BSC11","BSC12R","BSC13R")
colnames(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
omega(mydata) # alpha = .91, omega T = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 5.66; eigenvalue 2 = 1.08

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 6.25 eigenvalue 2 = 1.12

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.44, RMSEA=.161, RMSR=.1, TLI=.705

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.18, RMSR=.11, TLI=.692

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.103, RMSR=.04, TLI=.88
#      MR1   MR2
#MR1  1.00  -.49
#MR2  -.49  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.118, RMSR=.04, TLI=.867
#      MR1   MR2
#MR1  1.00  -.51
#MR2  -.51  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.953       0.836
#Tucker-Lewis Index (TLI)                       0.942       0.799
#Robust Comparative Fit Index (CFI)                         0.722
#Robust Tucker-Lewis Index (TLI)                            0.660
#RMSEA                                          0.184       0.209
#Robust RMSEA                                               0.193
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .501

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.744       0.731
#Tucker-Lewis Index (TLI)                       0.688       0.671
#Robust Comparative Fit Index (CFI)                         0.746
#Robust Tucker-Lewis Index (TLI)                            0.690
#RMSEA                                          0.168       0.144
#Robust RMSEA                                               0.167
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .428

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
EGAmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.888
#Tucker-Lewis Index (TLI)                       0.971       0.860
#Robust Comparative Fit Index (CFI)                         0.842
#Robust Tucker-Lewis Index (TLI)                            0.803
#RMSEA                                          0.130       0.174
#Robust RMSEA                                               0.147
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .56

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.862       0.857
#Tucker-Lewis Index (TLI)                       0.829       0.822
#Robust Comparative Fit Index (CFI)                         0.865
#Robust Tucker-Lewis Index (TLI)                            0.832
#RMSEA                                          0.124       0.106
#Robust RMSEA                                               0.123
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .51

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2          -0.630    0.037  -16.953    0.000   -0.630   -0.630

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.804       0.800
#Tucker-Lewis Index (TLI)                       0.760       0.756
#Robust Comparative Fit Index (CFI)                         0.806
#Robust Tucker-Lewis Index (TLI)                            0.763
#RMSEA                                          0.147       0.124
#Robust RMSEA                                               0.146
#SRMR                                           0.228       0.228

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .522


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.914       0.906
#Tucker-Lewis Index (TLI)                       0.865       0.853
#Robust Comparative Fit Index (CFI)                         0.916
#Robust Tucker-Lewis Index (TLI)                            0.867
#RMSEA                                          0.110       0.096
#Robust RMSEA                                               0.109
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .554

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6138914      0.4848485      0.7811188      0.1710992 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3970659 0.2622920 0.6029341 0.8960795 0.3402858 0.7085735 0.8113708
#factor2 0.3647839 0.1238166 0.6352161 0.8315457 0.2996817 0.5137885 0.7386230
#general 0.6138914 0.6138914 0.6138914 0.7811188 0.1710992 0.8786563 0.9000815





################################################################
### 
### Study 2 Sjastad (2024) https://osf.io/w24eb/download
### 
### Data https://osf.io/karuw/?view_only=81f20dfeacdc4e90aa83dee6e1738be5
###

library(readxl)
BSC_Sjastad2 <- read_excel("BSC_Sjastad2.xlsx")
colnames(BSC_Sjastad2)
mydata <- BSC_Sjastad2[,14:26]
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- colnames(mydata) <- c("BSC1","BSC2R","BSC3R","BSC4R","BSC5R","BSC6","BSC7R","BSC8","BSC9R","BSC10R","BSC11","BSC12R","BSC13R")
colnames(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
omega(mydata) # alpha = .87, omega T = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 5.63; eigenvalue 2 = 1.7

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 6.01 eigenvalue 2 = 1.76

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.43, RMSEA=.181, RMSR=.16, TLI=.683

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.195, RMSR=.16, TLI=.678

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.60, RMSEA=.094, RMSR=.04, TLI=.914
#      MR1   MR2
#MR1  1.00  -.11
#MR2  -.11  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.64, RMSEA=.108, RMSR=.04, TLI=.901
#      MR1   MR2
#MR1  1.00  -.13
#MR2  -.13  1.00


# Single factor model lavaan
UNImodel= '
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.620
#Tucker-Lewis Index (TLI)                       0.859       0.535
#Robust Comparative Fit Index (CFI)                         0.684
#Robust Tucker-Lewis Index (TLI)                            0.614
#RMSEA                                          0.312       0.362
#Robust RMSEA                                               0.217
#SRMR                                           0.176       0.176

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .598

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.705       0.702
#Tucker-Lewis Index (TLI)                       0.639       0.635
#Robust Comparative Fit Index (CFI)                         0.706
#Robust Tucker-Lewis Index (TLI)                            0.641
#RMSEA                                          0.195       0.165
#Robust RMSEA                                               0.195
#SRMR                                           0.158       0.158

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .566

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
EGAmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.929
#Tucker-Lewis Index (TLI)                       0.973       0.912
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.850
#RMSEA                                          0.135       0.158
#Robust RMSEA                                               0.135
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .654

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2          -0.206    0.024   -8.428    0.000   -0.206   -0.206

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.882
#Tucker-Lewis Index (TLI)                       0.856       0.853
#Robust Comparative Fit Index (CFI)                         0.886
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.124       0.105
#Robust RMSEA                                               0.122
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .603

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2          -0.149    0.041   -3.634    0.000   -0.149   -0.149

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.882       0.880
#Tucker-Lewis Index (TLI)                       0.856       0.853
#Robust Comparative Fit Index (CFI)                         0.884
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.124       0.105
#Robust RMSEA                                               0.122
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .603


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R
 factor2 =~ BSC1+BSC6+BSC11+BSC8
 general =~ BSC10R+BSC5R+BSC3R+BSC13R+BSC4R+BSC7R+BSC9R+BSC2R+BSC8+BSC1+BSC6+BSC11
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# does not run
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.993       0.965
#Tucker-Lewis Index (TLI)                       0.989       0.945
#Robust Comparative Fit Index (CFI)                         0.909
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.085       0.124
#Robust RMSEA                                               0.132
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .682

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    BSC10R            0.740    0.023   31.621    0.000    0.740    0.740
#    BSC5R             0.737    0.022   33.827    0.000    0.737    0.737
#    BSC3R             0.706    0.025   28.510    0.000    0.706    0.706
#    BSC13R            0.809    0.019   43.384    0.000    0.809    0.809
#    BSC4R             0.784    0.019   40.954    0.000    0.784    0.784
#    BSC7R             0.437    0.041   10.753    0.000    0.437    0.437
#    BSC9R             0.697    0.025   27.621    0.000    0.697    0.697
#    BSC2R             0.527    0.036   14.621    0.000    0.527    0.527
#  factor2 =~                                                            
#    BSC1              0.737    0.016   47.105    0.000    0.737    0.737
#    BSC6              0.723    0.018   40.398    0.000    0.723    0.723
#    BSC11             0.477    0.022   22.000    0.000    0.477    0.477
#    BSC8              0.832    0.020   41.060    0.000    0.832    0.832
#  general =~                                                            
#    BSC10R            0.405    0.042    9.717    0.000    0.405    0.405
#    BSC5R             0.373    0.041    9.120    0.000    0.373    0.373
#    BSC3R             0.426    0.040   10.604    0.000    0.426    0.426
#    BSC13R            0.243    0.046    5.297    0.000    0.243    0.243
#    BSC4R             0.193    0.049    3.954    0.000    0.193    0.193
#    BSC7R             0.684    0.032   21.213    0.000    0.684    0.684
#    BSC9R             0.430    0.040   10.684    0.000    0.430    0.430
#    BSC2R             0.674    0.031   21.679    0.000    0.674    0.674
#    BSC8             -0.055    0.037   -1.481    0.138   -0.055   -0.055
#    BSC1             -0.437    0.029  -15.148    0.000   -0.437   -0.437
#    BSC6             -0.229    0.034   -6.820    0.000   -0.229   -0.229
#    BSC11            -0.487    0.031  -15.578    0.000   -0.487   -0.487

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.2725699      0.4848485      0.9082765      0.1323887 
#$FactorLevelIndices
#        ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.6930070 0.4784229 0.3069930 0.9321287 0.6655066 0.8941571 0.9166158
#factor2 0.8041777 0.2490072 0.1958223 0.8328274 0.6932209 0.8283268 0.9094618
#general 0.2725699 0.2725699 0.2725699 0.9082765 0.1323887 0.7634687 0.8379590




