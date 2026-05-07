################################################################
### Analysis trait anxiety
### Sundelin https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9332334/ https://osf.io/tcdne/ (Study 1)



library(readxl)
STAI_T_Sundelin <- read_excel("STAI-T_Sundelin.xlsx")
colnames(STAI_T_Sundelin)
mydata  <- as.data.frame(STAI_T_Sundelin[,3:22])
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 5.69; eigenvalue 2 = .82

rho <- polychoric(mydata)$rho
# not equal number of responses
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Eigenvalue 1 = 6.99; eigenvalue 2 = .99

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.28, RMSEA=.074, RMSR=.07, TLI=.8

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.127, RMSR=.09, TLI=.644

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.34, RMSEA=.055, RMSR=.06, TLI=.889
#      MR1   MR2
#MR1  1.00 -0.52
#MR2 -0.52  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.113, RMSR=.07, TLI=.72
#      MR1   MR2
#MR1  1.00 -0.53
#MR2 -0.53  1.00

# Single factor model lavaan
UNImodel= '
 general=~ STT_1+STT_2+STT_3+STT_4+STT_5+STT_6+STT_7+STT_8+STT_9+STT_10+STT_11+STT_12+STT_13+STT_14+STT_15+STT_16+STT_17+STT_18+STT_19+STT_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.972       0.884
#Tucker-Lewis Index (TLI)                       0.968       0.870
#Robust Comparative Fit Index (CFI)                         0.736
#Robust Tucker-Lewis Index (TLI)                            0.705
#RMSEA                                          0.067       0.080
#Robust RMSEA                                               0.119
#SRMR                                           0.087       0.087

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .393

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.814       0.819
#Tucker-Lewis Index (TLI)                       0.792       0.797
#Robust Comparative Fit Index (CFI)                         0.822
#Robust Tucker-Lewis Index (TLI)                            0.802
#RMSEA                                          0.078       0.073
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .313

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on theory
EGAmodel= '
 factor1=~ STT_9+STT_20+STT_17+STT_18+STT_2+STT_5+STT_8+STT_11+STT_12+STT_3+STT_14
 factor2=~ STT_10+STT_16+STT_13+STT_19+STT_7+STT_6+STT_1+STT_4+STT_15
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.984       0.917
#Tucker-Lewis Index (TLI)                       0.982       0.907
#Robust Comparative Fit Index (CFI)                         0.790
#Robust Tucker-Lewis Index (TLI)                            0.764
#RMSEA                                          0.051       0.068
#Robust RMSEA                                               0.107
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .437

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.796    0.037  -21.307    0.000   -0.796   -0.796

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.864       0.868
#Tucker-Lewis Index (TLI)                       0.847       0.852
#Robust Comparative Fit Index (CFI)                         0.872
#Robust Tucker-Lewis Index (TLI)                            0.856
#RMSEA                                          0.067       0.062
#Robust RMSEA                                               0.064
#SRMR                                           0.066       0.066

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.775    0.060  -12.920    0.000   -0.775   -0.775

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .34


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
# ordered=TRUE does not converge
#Comparative Fit Index (CFI)                    0.772       0.772
#Tucker-Lewis Index (TLI)                       0.745       0.745
#Robust Comparative Fit Index (CFI)                         0.778
#Robust Tucker-Lewis Index (TLI)                            0.752
#RMSEA                                          0.086       0.082
#Robust RMSEA                                               0.084
#SRMR                                           0.187       0.187

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .31

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1=~ STT_9+STT_20+STT_17+STT_18+STT_2+STT_5+STT_8+STT_11+STT_12+STT_3+STT_14
 factor2=~ STT_10+STT_16+STT_13+STT_19+STT_7+STT_6+STT_1+STT_4+STT_15
 general=~ STT_9+STT_20+STT_17+STT_18+STT_2+STT_5+STT_8+STT_11+STT_12+STT_3+STT_14+
           STT_10+STT_16+STT_13+STT_19+STT_7+STT_6+STT_1+STT_4+STT_15
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# ordered=TRUE and MLR give error message; item STT_6 is not well estimated
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .38

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


#Reverse code items to see whether the combination of negative and positive loadings complicates matters
mydata2 <- mydata
#Items that are worded negatively
reverse_cols = c("STT_10","STT_16","STT_19","STT_13","STT_7","STT_6","STT_2")
mydata2[ , reverse_cols] = 5 - mydata2[ , reverse_cols]
CFA_model4 <- cfa(BIFmodel, data = mydata2,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.7447797      0.5210526      0.8993651      0.8067584 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1614873 0.08509597 0.8385127 0.8339224 0.04734907 0.4733769 0.7331762
#factor2 0.3596339 0.17012434 0.6403661 0.8410204 0.26490494 0.7243994 0.8889807
#general 0.7447797 0.74477969 0.7447797 0.8993651 0.80675839 0.9139852 0.9456321




#############################################
### Analysis trait anxiety Werner et al. (2022)
### https://zenodo.org/record/6523803

load("STAI_T_Werner.RData")
colnames(B07_STAI_T_items)
mydata <- B07_STAI_T_items[,3:22]
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors 2 components
# Eigenvalue 1 = 7.41; eigenvalue 2 = 1.10

rho <- polychoric(mydata)$rho
# warning
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 2 components but warning about estimates
# Eigenvalue 1 = 9.06; eigenvalue 2 = 1.28

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.37, RMSEA=.109, RMSR=.08, TLI=.741

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.163, RMSR=.09, TLI=.635

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.082, RMSR=.05, TLI=.853
#      MR1   MR2
#MR2 1.00 0.64
#MR1 0.64 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.135, RMSR=.06, TLI=.75
#      MR2   MR1
#MR2 1.00 0.65
#MR1 0.65 1.00



# Single factor model lavaan
UNImodel= '
 general=~ STAI1+STAI2+STAI3+STAI4+STAI5+STAI6+STAI7+STAI8+STAI9+STAI10+STAI11+STAI12+STAI13+STAI14+STAI15+STAI16+STAI17+STAI18+STAI19+STAI20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.970       0.903
#Tucker-Lewis Index (TLI)                       0.966       0.892
#Robust Comparative Fit Index (CFI)                         0.703
#Robust Tucker-Lewis Index (TLI)                            0.668
#RMSEA                                          0.103       0.112
#Robust RMSEA                                               0.158
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .487

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.769       0.768
#Tucker-Lewis Index (TLI)                       0.742       0.741
#Robust Comparative Fit Index (CFI)                         0.772
#Robust Tucker-Lewis Index (TLI)                            0.745
#RMSEA                                          0.110       0.105
#Robust RMSEA                                               0.109
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .36

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on factor analysis
EGAmodel= '
 factor1=~ STAI9+STAI17+STAI14+STAI11+STAI18+STAI20+STAI8+STAI3+STAI12+STAI15+STAI2+STAI5
 factor2=~ STAI10+STAI16+STAI1+STAI19+STAI13+STAI6+STAI7+STAI14
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.987       0.952
#Tucker-Lewis Index (TLI)                       0.985       0.946
#Robust Comparative Fit Index (CFI)                         0.804
#Robust Tucker-Lewis Index (TLI)                            0.777
#RMSEA                                          0.069       0.082
#Robust RMSEA                                               0.133
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .537

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.748    0.027   27.958    0.000    0.748    0.748

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.865       0.865
#Tucker-Lewis Index (TLI)                       0.846       0.847
#Robust Comparative Fit Index (CFI)                         0.868
#Robust Tucker-Lewis Index (TLI)                            0.850
#RMSEA                                          0.088       0.083
#Robust RMSEA                                               0.087
#SRMR                                           0.055       0.055

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.738    0.031   23.660    0.000    0.738    0.738

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .41

semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.800       0.799
#Tucker-Lewis Index (TLI)                       0.774       0.772
#Robust Comparative Fit Index (CFI)                         0.803
#Robust Tucker-Lewis Index (TLI)                            0.777
#RMSEA                                          0.107       0.101
#Robust RMSEA                                               0.106
#SRMR                                           0.227       0.227

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .434

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1=~ STAI9+STAI17+STAI14+STAI11+STAI18+STAI20+STAI8+STAI3+STAI12+STAI15+STAI2+STAI5
 factor2=~ STAI10+STAI16+STAI1+STAI19+STAI13+STAI6+STAI7+STAI14
 general=~ STAI9+STAI17+STAI14+STAI11+STAI18+STAI20+STAI8+STAI3+STAI12+STAI15+STAI2+STAI5+
           STAI10+STAI16+STAI1+STAI19+STAI13+STAI6+STAI7+STAI14
'

CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.911       0.910
#Tucker-Lewis Index (TLI)                       0.885       0.884
#Robust Comparative Fit Index (CFI)                         0.914
#Robust Tucker-Lewis Index (TLI)                            0.888
#RMSEA                                          0.076       0.072
#Robust RMSEA                                               0.075
#SRMR                                           0.045       0.045

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .444

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    STAI9             0.441    0.045    9.695    0.000    0.441    0.563
#    STAI17            0.385    0.042    9.255    0.000    0.385    0.528
#    STAI14            0.395    0.041    9.574    0.000    0.395    0.531
#    STAI11            0.395    0.036   10.879    0.000    0.395    0.496
#    STAI18            0.356    0.044    8.125    0.000    0.356    0.430
#    STAI20            0.303    0.042    7.147    0.000    0.303    0.398
#    STAI8             0.252    0.034    7.392    0.000    0.252    0.377
#    STAI3             0.218    0.034    6.457    0.000    0.218    0.384
#    STAI12            0.319    0.049    6.525    0.000    0.319    0.403
#    STAI15            0.198    0.033    6.091    0.000    0.198    0.343
#    STAI2             0.214    0.046    4.690    0.000    0.214    0.281
#    STAI5             0.210    0.041    5.087    0.000    0.210    0.297
#  factor2 =~                                                            
#    STAI10            0.315    0.089    3.534    0.000    0.315    0.418
#    STAI16            0.198    0.066    3.000    0.003    0.198    0.263
#    STAI1             0.115    0.064    1.797    0.072    0.115    0.169
#    STAI19           -0.180    0.081   -2.222    0.026   -0.180   -0.218
#    STAI13            0.223    0.071    3.149    0.002    0.223    0.291
#    STAI6            -0.121    0.068   -1.794    0.073   -0.121   -0.160
#    STAI7            -0.261    0.089   -2.926    0.003   -0.261   -0.315
#    STAI14           -0.007    0.042   -0.156    0.876   -0.007   -0.009
#  general =~                                                            
#    STAI9             0.303    0.039    7.690    0.000    0.303    0.387
#    STAI17            0.351    0.036    9.782    0.000    0.351    0.482
#    STAI14            0.258    0.038    6.747    0.000    0.258    0.348
#    STAI11            0.413    0.039   10.653    0.000    0.413    0.519
#    STAI18            0.381    0.041    9.336    0.000    0.381    0.460
#    STAI20            0.407    0.039   10.444    0.000    0.407    0.534
#    STAI8             0.369    0.032   11.355    0.000    0.369    0.552
#    STAI3             0.219    0.030    7.324    0.000    0.219    0.385
#    STAI12            0.324    0.045    7.207    0.000    0.324    0.409
#    STAI15            0.298    0.029   10.114    0.000    0.298    0.516
#    STAI2             0.287    0.039    7.400    0.000    0.287    0.377
#    STAI5             0.237    0.037    6.485    0.000    0.237    0.334
#    STAI10            0.608    0.044   13.938    0.000    0.608    0.806
#    STAI16            0.614    0.036   17.218    0.000    0.614    0.816
#    STAI1             0.476    0.032   15.000    0.000    0.476    0.697
#    STAI19            0.690    0.038   17.942    0.000    0.690    0.835
#    STAI13            0.450    0.046    9.832    0.000    0.450    0.587
#    STAI6             0.473    0.038   12.347    0.000    0.473    0.625
#    STAI7             0.598    0.046   12.902    0.000    0.598    0.724

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general  Omega.general OmegaH.general 
#     0.6928149      0.9296604      0.7520416 
#$FactorLevelIndices
#        ECV_SS    ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#factor1 0.4777303 0.2477275 0.5222527 0.8785920 0.416164720 0.7392235 0.8643261
#factor2 0.1128711 0.0594577 0.8268867 0.9006397 0.005718377 0.3725500 0.7966198
#general 0.6928149 0.6928149 0.6928149 0.9296604 0.752041614 0.9257606 0.9625123






################################################################
### Analysis trait anxiety India
### Joshi et al. (2023) https://www.longdom.org/open-access/convergent-validation-of-the-statetrait-anxiety-inventory-with--measures-of-personality-affective-control-and-risk-prope.pdf
### Data sent via email by Bhoomika Kar


library(readxl)
STAI_Joshi <- read_excel("STAI_Joshi.xlsx")
colnames(STAI_Joshi)
mydata  <- as.data.frame(STAI_Joshi[,22:41])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 2 components
# Eigenvalue 1 = 3.94; eigenvalue 2 = .95

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components 
# Eigenvalue 1 = 4.69; eigenvalue 2 = 1.12

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.20, RMSEA=.067, RMSR=.07, TLI=.741

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.091, RMSR=.08, TLI=.668

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communities

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.048, RMSR=.05, TLI=.87
#     MR1  MR2
#MR1 1.00 0.53
#MR2 0.53 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.072, RMSR=.06, TLI=.793
#     MR2  MR1
#MR2 1.00 0.52
#MR1 0.52 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.928       0.844
#Tucker-Lewis Index (TLI)                       0.919       0.826
#Robust Comparative Fit Index (CFI)                         0.725
#Robust Tucker-Lewis Index (TLI)                            0.693
#RMSEA                                          0.071       0.075
#Robust RMSEA                                               0.089
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .24

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.767       0.781
#Tucker-Lewis Index (TLI)                       0.739       0.755
#Robust Comparative Fit Index (CFI)                         0.782
#Robust Tucker-Lewis Index (TLI)                            0.757
#RMSEA                                          0.068       0.060
#Robust RMSEA                                               0.065
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .190

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on EGA analysis
EGAmodel= '
 factor1=~ t2+t4+t5+t8+t9+t11+t12+t15+t17+t18+t20
 factor2=~ t1+t3+t6+t7+t10+t13+t14+t16+t19
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.976       0.935
#Tucker-Lewis Index (TLI)                       0.973       0.927
#Robust Comparative Fit Index (CFI)                         0.844
#Robust Tucker-Lewis Index (TLI)                            0.825
#RMSEA                                          0.041       0.049
#Robust RMSEA                                               0.067
#SRMR                                           0.059       0.059


#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.626    0.043   14.453    0.000    0.626    0.626

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .29

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.885       0.902
#Tucker-Lewis Index (TLI)                       0.871       0.890
#Robust Comparative Fit Index (CFI)                         0.903
#Robust Tucker-Lewis Index (TLI)                            0.891
#RMSEA                                          0.048       0.040
#Robust RMSEA                                               0.044
#SRMR                                           0.050       0.050

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .241

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.620    0.058   10.709    0.000    0.620    0.620


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.810       0.823
#Tucker-Lewis Index (TLI)                       0.787       0.803
#Robust Comparative Fit Index (CFI)                         0.826
#Robust Tucker-Lewis Index (TLI)                            0.805
#RMSEA                                          0.062       0.054
#Robust RMSEA                                               0.058
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .25

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1=~ t2+t4+t5+t8+t9+t11+t12+t15+t17+t18+t20
 factor2=~ t1+t3+t6+t7+t10+t13+t14+t16+t19
 general=~ t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.917       0.921
#Tucker-Lewis Index (TLI)                       0.895       0.900
#Robust Comparative Fit Index (CFI)                         0.926
#Robust Tucker-Lewis Index (TLI)                            0.907
#RMSEA                                          0.043       0.039
#Robust RMSEA                                               0.040
#SRMR                                           0.044       0.044

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .247

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    t2                0.044    0.078    0.565    0.572    0.044    0.070
#    t4                0.167    0.181    0.925    0.355    0.167    0.151
#    t5               -0.124    0.083   -1.491    0.136   -0.124   -0.153
#    t8                0.108    0.149    0.720    0.471    0.108    0.134
#    t9                0.481    0.118    4.070    0.000    0.481    0.512
#    t11               0.116    0.072    1.627    0.104    0.116    0.142
#    t12              -0.142    0.284   -0.500    0.617   -0.142   -0.153
#    t15              -0.010    0.212   -0.049    0.961   -0.010   -0.013
#    t17               0.352    0.179    1.964    0.050    0.352    0.409
#    t18               0.149    0.231    0.647    0.517    0.149    0.164
#    t20               0.092    0.259    0.357    0.721    0.092    0.103
#  factor2 =~                                                            
#    t1                0.369    0.055    6.711    0.000    0.369    0.474
#    t3                0.394    0.099    3.971    0.000    0.394    0.426
#    t6                0.305    0.055    5.539    0.000    0.305    0.377
#    t7                0.351    0.066    5.285    0.000    0.351    0.389
#    t10               0.466    0.058    8.080    0.000    0.466    0.587
#    t13               0.323    0.072    4.498    0.000    0.323    0.388
#    t14               0.234    0.080    2.906    0.004    0.234    0.260
#    t16               0.286    0.065    4.422    0.000    0.286    0.340
#    t19               0.296    0.065    4.554    0.000    0.296    0.335
#  general =~                                                            
#    t1                0.265    0.047    5.649    0.000    0.265    0.341
#    t2                0.288    0.044    6.533    0.000    0.288    0.455
#    t3                0.340    0.092    3.679    0.000    0.340    0.369
#    t4                0.218    0.091    2.397    0.017    0.218    0.197
#    t5                0.482    0.097    4.963    0.000    0.482    0.595
#    t6                0.061    0.050    1.237    0.216    0.061    0.076
#    t7                0.328    0.058    5.695    0.000    0.328    0.362
#    t8                0.297    0.062    4.791    0.000    0.297    0.371
#    t9                0.400    0.152    2.637    0.008    0.400    0.426
#    t10               0.308    0.051    6.065    0.000    0.308    0.388
#    t11               0.354    0.055    6.451    0.000    0.354    0.431
#    t12               0.580    0.054   10.727    0.000    0.580    0.628
#    t13               0.347    0.059    5.923    0.000    0.347    0.417
#    t14               0.301    0.062    4.864    0.000    0.301    0.335
#    t15               0.360    0.049    7.355    0.000    0.360    0.471
#    t16               0.272    0.049    5.530    0.000    0.272    0.325
#    t17               0.422    0.110    3.846    0.000    0.422    0.491
#    t18               0.406    0.082    4.926    0.000    0.406    0.447
#    t19               0.204    0.059    3.449    0.001    0.204    0.231
#    t20               0.454    0.075    6.058    0.000    0.454    0.508


library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.6216882      0.5210526      0.8403754      0.6792451 

#$FactorLevelIndices
#ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1933397 0.1058915 0.8066603 0.7717145 0.05321853 0.4150868 0.6749373
#factor2 0.6022953 0.2724203 0.3977047 0.7619623 0.46693279 0.6539613 0.7959797
#general 0.6216882 0.6216882 0.6216882 0.8403754 0.67924507 0.8162338 0.8894753





