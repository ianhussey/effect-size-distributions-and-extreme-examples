################################################################
### Analysis Moral Foundations Questionnaire
### Nottingham Trent University: Harper & Rhodes https://psyarxiv.com/p5fj8/

MoralF_Harper <- read.csv("MoralF_Harper.csv")
colnames(MoralF_Harper)
mydata  <- as.data.frame(MoralF_Harper[,25:65])
print(colnames(mydata))
mydata <- mydata[ , -which(names(mydata) %in% c("math","good","care"))]
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 6 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 7.07; eigenvalue 2 = 3.33

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components
# Eigenvalue 1 = 8.02; eigenvalue 2 = 3.73

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.102, RMSR=.12, TLI=.422

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.21, RMSEA=.124, RMSR=.14, TLI=.368

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities 

# Give solution with 5 factors
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.05, RMSR=.04, TLI=.846
#      MR1   MR4   MR2   MR3   MR5
#MR1  1.00  0.45 -0.03  0.00 -0.18
#MR4  0.45  1.00  0.01  0.05 -0.28
#MR2 -0.03  0.01  1.00  0.00  0.16
#MR3  0.00  0.05  0.00  1.00 -0.15
#MR5 -0.18 -0.28  0.16 -0.15  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.078, RMSR=.04, TLI=.752
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.03  0.42 -0.01  0.34
#MR2 -0.03  1.00 -0.03  0.03 -0.13
#MR4  0.42 -0.03  1.00 -0.05  0.21
#MR3 -0.01  0.03 -0.05  1.00  0.14
#MR5  0.34 -0.13  0.21  0.14  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ emotionally+treated+lovecountry+respect+decency+weak+unfairly+betray+traditions+
           disgusting+cruel+rights+loyalty+chaos+god+liberty_1+liberty_2+compassion+
           fairly+history+kidrespect+harmlessdg+liberty_3+animal+justice+family+sexroles+
           unnatural+liberty_4+liberty_5+kill+rich+team+liberty_6+soldier+chastity+liberty_7+liberty_8+liberty_9 
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.728       0.546
#Tucker-Lewis Index (TLI)                       0.713       0.521
#Robust Comparative Fit Index (CFI)                         0.413
#Robust Tucker-Lewis Index (TLI)                            0.381
#RMSEA                                          0.159       0.117
#Robust RMSEA                                               0.126
#SRMR                                           0.137       0.137

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .17

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.451       0.448
#Tucker-Lewis Index (TLI)                       0.421       0.417
#Robust Comparative Fit Index (CFI)                         0.455
#Robust Tucker-Lewis Index (TLI)                            0.425
#RMSEA                                          0.104       0.098
#Robust RMSEA                                               0.103
#SRMR                                           0.121       0.121

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .117

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors based on theoretical analysis
EGAmodel= '
 factor1 =~ loyalty+betray+traditions+lovecountry+god+family
 factor2 =~ cruel+rights+unfairly+treated+compassion+weak+emotionally+fairly+liberty_2+sexroles
 factor3 =~ harmlessdg+kidrespect+unnatural+kill+respect+decency+chastity+animal+disgusting+justice+soldier+team
 factor4 =~ liberty_4+liberty_5+liberty_8+liberty_6+liberty_9+liberty_7
 factor5 =~ rich+liberty_3+chaos+history
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.836       0.717
#Tucker-Lewis Index (TLI)                       0.824       0.696
#Robust Comparative Fit Index (CFI)                         0.625
#Robust Tucker-Lewis Index (TLI)                            0.597
#RMSEA                                          0.125       0.095
#Robust RMSEA                                               0.102
#SRMR                                           0.117       0.117

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.106    0.049   -2.143    0.032   -0.106   -0.106
#  factor3           0.774    0.023   33.066    0.000    0.774    0.774
#  factor4           0.185    0.048    3.885    0.000    0.185    0.185
#  factor5          -0.616    0.038  -16.165    0.000   -0.616   -0.616
#factor2 ~~                                                            
#  factor3          -0.071    0.050   -1.407    0.159   -0.071   -0.071
#  factor4          -0.076    0.057   -1.333    0.182   -0.076   -0.076
#  factor5           0.416    0.048    8.702    0.000    0.416    0.416
#factor3 ~~                                                            
#  factor4           0.050    0.049    1.025    0.305    0.050    0.050
#  factor5          -0.686    0.034  -20.340    0.000   -0.686   -0.686
#factor4 ~~                                                            
#  factor5          -0.494    0.043  -11.448    0.000   -0.494   -0.494

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .42


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.691       0.689
#Tucker-Lewis Index (TLI)                       0.668       0.667
#Robust Comparative Fit Index (CFI)                         0.696
#Robust Tucker-Lewis Index (TLI)                            0.674
#RMSEA                                          0.079       0.075
#Robust RMSEA                                               0.078
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .329

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.050    0.088   -0.567    0.571   -0.050   -0.050
#    factor3           0.758    0.034   22.004    0.000    0.758    0.758
#    factor4           0.178    0.067    2.654    0.008    0.178    0.178
#    factor5          -0.557    0.063   -8.820    0.000   -0.557   -0.557
#  factor2 ~~                                                            
#    factor3          -0.031    0.081   -0.379    0.705   -0.031   -0.031
#    factor4          -0.073    0.100   -0.732    0.464   -0.073   -0.073
#    factor5           0.423    0.085    4.977    0.000    0.423    0.423
#  factor3 ~~                                                            
#    factor4          -0.003    0.067   -0.042    0.966   -0.003   -0.003
#    factor5          -0.610    0.063   -9.751    0.000   -0.610   -0.610
#  factor4 ~~                                                            
#    factor5          -0.490    0.066   -7.412    0.000   -0.490   -0.490


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.605       0.603
#Tucker-Lewis Index (TLI)                       0.582       0.581
#Robust Comparative Fit Index (CFI)                         0.610
#Robust Tucker-Lewis Index (TLI)                            0.587
#RMSEA                                          0.089       0.084
#Robust RMSEA                                               0.088
#SRMR                                           0.158       0.158

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .327

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ loyalty+betray+traditions+lovecountry+god+family
 factor2 =~ cruel+rights+unfairly+treated+compassion+weak+emotionally+fairly+liberty_2+sexroles
 factor3 =~ harmlessdg+kidrespect+unnatural+kill+respect+decency+chastity+animal+disgusting+justice+soldier+team
 factor4 =~ liberty_4+liberty_5+liberty_8+liberty_6+liberty_9+liberty_7
 factor5 =~ rich+liberty_3+chaos+history
 general=~ emotionally+treated+lovecountry+respect+decency+weak+unfairly+betray+traditions+
           disgusting+cruel+rights+loyalty+chaos+god+liberty_1+liberty_2+compassion+
           fairly+history+kidrespect+harmlessdg+liberty_3+animal+justice+family+sexroles+
           unnatural+liberty_4+liberty_5+kill+rich+team+liberty_6+soldier+chastity+liberty_7+liberty_8+liberty_9 
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.738       0.736
#Tucker-Lewis Index (TLI)                       0.708       0.705
#Robust Comparative Fit Index (CFI)                         0.743
#Robust Tucker-Lewis Index (TLI)                            0.713
#RMSEA                                          0.074       0.070
#Robust RMSEA                                               0.073
#SRMR                                           0.087       0.087

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .363

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4528938      0.8016194      0.8531503      0.5783863 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#factor1 0.3352358 0.06723575 0.66476416 0.8263719 0.23426340 0.5936943 0.8050837
#factor2 0.8734398 0.19308234 0.12656022 0.7559033 0.75348432 0.8075299 0.9005254
#factor3 0.3364297 0.10268081 0.66357027 0.8408993 0.30803947 0.6600917 0.8171761
#factor4 0.9010449 0.12578050 0.09895514 0.6009910 0.59536088 0.7696759 0.8822281
#factor5 0.4944143 0.05832680 0.50558574 0.3617891 0.03557842 0.5712831 0.7974701
#general 0.4528938 0.45289381 0.45289381 0.8531503 0.57838631 0.9083153 0.9364476











################################################################
### Analysis Moral Foundations Questionnaire
### Zakharin & Bates, 2021, study 1 https://osf.io/bt9nh/

MoralF_Zakharin <- read.csv("MoralF_Zakharin.csv")
colnames(MoralF_Zakharin)
mydata  <- as.data.frame(MoralF_Zakharin[,c(3:7,9:23,25:34)])
print(colnames(mydata))
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 6 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 5 components
# Eigenvalue 1 = 5.41; eigenvalue 2 = 3.09

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 8 factprs a,d 5 components
# Eigenvalue 1 = 5.93; eigenvalue 2 = 3.58

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.103, RMSR=.12, TLI=.455

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.119, RMSR=.14, TLI=.417

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities 

# Give solution with 5 factors
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.053, RMSR=.03, TLI=.857
#      MR2   MR1   MR5   MR4   MR3
#MR2  1.00  0.10 -0.02  0.28  0.14
#MR1  0.10  1.00  0.47  0.41 -0.07
#MR5 -0.02  0.47  1.00  0.40 -0.05
#MR4  0.28  0.41  0.40  1.00 -0.02
#MR3  0.14 -0.07 -0.05 -0.02  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.064, RMSR=.03, TLI=.831
#      MR2   MR1   MR4  MR5   MR3
#MR2  1.00  0.09 -0.05 0.26  0.21
#MR1  0.09  1.00  0.44 0.41 -0.04
#MR4 -0.05  0.44  1.00 0.38 -0.04
#MR5  0.26  0.41  0.38 1.00  0.02
#MR3  0.21 -0.04 -0.04 0.02  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ H1+H2+H3+H4+H5+H6+
           F1+F2+F3+F4+F5+F6+
           I1+I2+I3+I4+I5+I6+
           A1+A2+A3+A4+A5+A6+
           P1+P2+P3+P4+P5+P6
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.704       0.381
#Tucker-Lewis Index (TLI)                       0.682       0.335
#Robust Comparative Fit Index (CFI)                         0.450
#Robust Tucker-Lewis Index (TLI)                            0.409
#RMSEA                                          0.160       0.144
#Robust RMSEA                                               0.120
#SRMR                                           0.136       0.136

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .179

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.494       0.486
#Tucker-Lewis Index (TLI)                       0.457       0.448
#Robust Comparative Fit Index (CFI)                         0.496
#Robust Tucker-Lewis Index (TLI)                            0.458
#RMSEA                                          0.103       0.095
#Robust RMSEA                                               0.103
#SRMR                                           0.116       0.116

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .179

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
 factor1 =~ H1+H2+H3+H4+H5+H6
 factor2 =~ F1+F2+F3+F4+F5+F6
 factor3 =~ I1+I2+I3+I4+I5+I6
 factor4 =~ A1+A2+A3+A4+A5+A6
 factor5 =~ P1+P2+P3+P4+P5+P6
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.885       0.743
#Tucker-Lewis Index (TLI)                       0.873       0.717
#Robust Comparative Fit Index (CFI)                         0.706
#Robust Tucker-Lewis Index (TLI)                            0.676
#RMSEA                                          0.101       0.094
#Robust RMSEA                                               0.089
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .381

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.984    0.012   79.817    0.000    0.984    0.984
#    factor3           0.317    0.024   13.239    0.000    0.317    0.317
#    factor4           0.176    0.026    6.803    0.000    0.176    0.176
#    factor5           0.231    0.024    9.450    0.000    0.231    0.231
#  factor2 ~~                                                            
#    factor3           0.279    0.025   11.304    0.000    0.279    0.279
#    factor4           0.129    0.027    4.879    0.000    0.129    0.129
#    factor5           0.100    0.025    4.053    0.000    0.100    0.100
#  factor3 ~~                                                            
#    factor4           0.883    0.012   74.502    0.000    0.883    0.883
#    factor5           0.775    0.015   50.757    0.000    0.775    0.775
#  factor4 ~~                                                            
#    factor5           0.866    0.012   72.075    0.000    0.866    0.866


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.736       0.733
#Tucker-Lewis Index (TLI)                       0.709       0.706
#Robust Comparative Fit Index (CFI)                         0.738
#Robust Tucker-Lewis Index (TLI)                            0.712
#RMSEA                                          0.075       0.070
#Robust RMSEA                                               0.075
#SRMR                                           0.078       0.078

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.950    0.018   52.785    0.000    0.950    0.950
#    factor3           0.332    0.037    8.912    0.000    0.332    0.332
#    factor4           0.203    0.037    5.538    0.000    0.203    0.203
#    factor5           0.209    0.036    5.782    0.000    0.209    0.209
#  factor2 ~~                                                            
#    factor3           0.267    0.041    6.439    0.000    0.267    0.267
#    factor4           0.134    0.040    3.372    0.001    0.134    0.134
#    factor5           0.084    0.038    2.226    0.026    0.084    0.084
#  factor3 ~~                                                            
#    factor4           0.885    0.021   41.627    0.000    0.885    0.885
#    factor5           0.799    0.022   36.353    0.000    0.799    0.799
#  factor4 ~~                                                            
#    factor5           0.855    0.018   46.355    0.000    0.855    0.855

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .31


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.521       0.514
#Tucker-Lewis Index (TLI)                       0.485       0.478
#Robust Comparative Fit Index (CFI)                         0.522
#Robust Tucker-Lewis Index (TLI)                            0.487
#RMSEA                                          0.100       0.093
#Robust RMSEA                                               0.100
#SRMR                                           0.172       0.172

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .282

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ H1+H2+H3+H4+H5+H6
 factor2 =~ F1+F2+F3+F4+F5+F6
 factor3 =~ I1+I2+I3+I4+I5+I6
 factor4 =~ A1+A2+A3+A4+A5+A6
 factor5 =~ P1+P2+P3+P4+P5+P6
 general=~ H1+H2+H3+H4+H5+H6+
           F1+F2+F3+F4+F5+F6+
           I1+I2+I3+I4+I5+I6+
           A1+A2+A3+A4+A5+A6+
           P1+P2+P3+P4+P5+P6
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.675       0.666
#Tucker-Lewis Index (TLI)                       0.624       0.613
#Robust Comparative Fit Index (CFI)                         0.677
#Robust Tucker-Lewis Index (TLI)                            0.625
#RMSEA                                          0.086       0.080
#Robust RMSEA                                               0.085
#SRMR                                           0.136       0.136

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .353

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    H1                0.608    0.018   34.508    0.000    0.608    0.608
#    H2                0.626    0.018   34.371    0.000    0.626    0.626
#    H3                0.587    0.019   31.003    0.000    0.587    0.587
#    H4                0.584    0.019   30.739    0.000    0.584    0.584
#    H5                0.297    0.023   12.927    0.000    0.297    0.297
#    H6                0.215    0.025    8.773    0.000    0.215    0.215
#  factor2 =~                                                            
#    F1                0.660    0.016   40.499    0.000    0.660    0.660
#    F2                0.620    0.017   36.118    0.000    0.620    0.620
#    F3                0.717    0.018   39.682    0.000    0.717    0.717
#    F4                0.488    0.021   23.270    0.000    0.488    0.488
#    F5                0.271    0.022   12.436    0.000    0.271    0.271
#    F6                0.248    0.023   10.585    0.000    0.248    0.248
#  factor3 =~                                                            
#    I1               -0.068    0.024   -2.811    0.005   -0.068   -0.068
#    I2                0.491    0.038   12.770    0.000    0.491    0.491
#    I3                0.510    0.040   12.630    0.000    0.510    0.510
#    I4               -0.200    0.028   -7.259    0.000   -0.200   -0.200
#    I5                0.011    0.027    0.397    0.692    0.011    0.011
#    I6               -0.111    0.029   -3.875    0.000   -0.111   -0.111
#  factor4 =~                                                            
#    A1                0.269    0.027    9.979    0.000    0.269    0.269
#    A2               -0.005    0.026   -0.197    0.844   -0.005   -0.005
#    A3               -0.023    0.027   -0.836    0.403   -0.023   -0.023
#    A4                0.636    0.048   13.141    0.000    0.636    0.636
#    A5                0.143    0.027    5.283    0.000    0.143    0.143
#    A6                0.294    0.030    9.970    0.000    0.294    0.294
#  factor5 =~                                                            
#    P1                0.061    0.025    2.482    0.013    0.061    0.061
#    P2                0.280    0.022   12.763    0.000    0.280    0.280
#    P3                0.475    0.025   19.230    0.000    0.475    0.475
#    P4                0.448    0.022   20.493    0.000    0.448    0.448
#    P5                0.375    0.023   16.649    0.000    0.375    0.375
#    P6                0.470    0.024   19.665    0.000    0.470    0.470
#  general =~                                                            
#    H1                0.203    0.021    9.838    0.000    0.203    0.203
#    H2                0.269    0.020   13.492    0.000    0.269    0.269
#    H3                0.298    0.020   14.842    0.000    0.298    0.298
#    H4                0.226    0.021   10.807    0.000    0.226    0.226
#    H5                0.204    0.022    9.308    0.000    0.204    0.204
#    H6                0.202    0.021    9.445    0.000    0.202    0.202
#    F1                0.194    0.021    9.149    0.000    0.194    0.194
#    F2                0.271    0.020   13.655    0.000    0.271    0.271
#    F3                0.171    0.021    8.050    0.000    0.171    0.171
#    F4                0.219    0.023    9.730    0.000    0.219    0.219
#    F5                0.361    0.020   18.168    0.000    0.361    0.361
#    F6               -0.031    0.022   -1.421    0.155   -0.031   -0.031
#    I1                0.651    0.013   49.450    0.000    0.651    0.651
#    I2                0.612    0.014   42.658    0.000    0.612    0.612
#    I3                0.622    0.014   44.708    0.000    0.622    0.622
#    I4                0.509    0.017   30.656    0.000    0.509    0.509
#    I5                0.464    0.017   27.144    0.000    0.464    0.464
#    I6                0.332    0.019   17.483    0.000    0.332    0.332
#    A1                0.716    0.011   62.485    0.000    0.716    0.716
#    A2                0.547    0.016   34.495    0.000    0.547    0.547
#    A3                0.535    0.016   33.720    0.000    0.535    0.535
#    A4                0.591    0.015   39.379    0.000    0.591    0.591
#    A5                0.400    0.019   21.520    0.000    0.400    0.400
#    A6                0.397    0.018   22.330    0.000    0.397    0.397
#    P1                0.650    0.013   50.094    0.000    0.650    0.650
#    P2                0.552    0.015   37.220    0.000    0.552    0.552
#    P3                0.468    0.019   24.858    0.000    0.468    0.468
#    P4                0.521    0.016   31.894    0.000    0.521    0.521
#    P5                0.504    0.016   31.200    0.000    0.504    0.504
#    P6                0.410    0.019   22.065    0.000    0.410    0.410

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.5237595      0.8275862      0.8801331      0.7722084 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.8250639 0.14201555 0.1749361 0.6383236 0.50866735 0.7073096 0.8457953
#factor2 0.8422638 0.15342867 0.1577362 0.6385893 0.54668726 0.7444607 0.8664530
#factor3 0.2398264 0.05013941 0.7601736 0.6649550 0.01952701 0.4211829 0.7402264
#factor4 0.2484952 0.05245109 0.7515048 0.6993219 0.10613366 0.4657121 0.7692544
#factor5 0.3467017 0.07820574 0.6532983 0.7586051 0.24599120 0.5190310 0.7339817
#general 0.5237595 0.52375953 0.5237595 0.8801331 0.77220837 0.8956256 0.9350770

