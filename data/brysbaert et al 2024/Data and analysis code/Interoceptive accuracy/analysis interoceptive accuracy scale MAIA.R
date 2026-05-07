################################################################
### Multidimensional assessment of interoceptive awareness (MAIA)
### 
### Theory-based 8-factor hierarchical model
### Various versions of the test with different numbers of items
###



################################################################
### 
### Rogowska et al (2023) https://www.nature.com/articles/s41598-023-48536-0#MOESM1
### data Study https://www.nature.com/articles/s41598-023-48536-0#MOESM1

library(readxl)
MAIA_Rogowska <- read_excel("MAIA_Rogowska.xlsx")
colnames(MAIA_Rogowska)
mydata <- MAIA_Rogowska[,6:42]
mydata <- na.omit(mydata)
colnames(mydata)
Names <- colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 6

library(psych)
omega(mydata) # alpha = .94, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 11.48; eigenvalue 2 = 3.08

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 5 components
# Eigenvalue 1 = 12.02; eigenvalue 2 = 3.33

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.31, RMSEA=.118, RMSR=.12, TLI=.535

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.137, RMSR=.13, TLI=.474

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 8 factors (theory-based)
fit4 <- fa(mydata,8)
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.046, RMSR=.02, TLI=.928
#      MR1   MR6   MR3   MR5   MR4   MR2   MR7   MR8
#MR1  1.00  0.54 -0.36  0.38  0.34  0.00  0.33  0.15
#MR6  0.54  1.00 -0.22  0.33  0.38  0.05  0.19  0.14
#MR3 -0.36 -0.22  1.00 -0.23 -0.12  0.10 -0.17 -0.05
#MR5  0.38  0.33 -0.23  1.00  0.26 -0.25  0.47  0.39
#MR4  0.34  0.38 -0.12  0.26  1.00 -0.03  0.20  0.11
#MR2  0.00  0.05  0.10 -0.25 -0.03  1.00 -0.35 -0.27
#MR7  0.33  0.19 -0.17  0.47  0.20 -0.35  1.00  0.32
#MR8  0.15  0.14 -0.05  0.39  0.11 -0.27  0.32  1.00

fit4 <- fa(rho,8,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.072, RMSR=.03, TLI=.854
#      MR1   MR6   MR3   MR5   MR4   MR2   MR7   MR8
#MR1  1.00  0.54 -0.36  0.39  0.33  0.00  0.27  0.12
#MR6  0.54  1.00 -0.20  0.36  0.37  0.04  0.11  0.08
#MR3 -0.36 -0.20  1.00 -0.20 -0.10  0.09 -0.14 -0.02
#MR5  0.39  0.36 -0.20  1.00  0.30 -0.27  0.42  0.36
#MR4  0.33  0.37 -0.10  0.30  1.00 -0.07  0.16  0.07
#MR2  0.00  0.04  0.09 -0.27 -0.07  1.00 -0.33 -0.26
#MR7  0.27  0.11 -0.14  0.42  0.16 -0.33  1.00  0.26
#MR8  0.12  0.08 -0.02  0.36  0.07 -0.26  0.26  1.00


# Single factor model lavaan
UNImodel= '
 general =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04+MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+
            MAIA2_09+MAIA2_10+MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15+MAIA2_16+
            MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22+MAIA2_23+MAIA2_24+
            MAIA2_25+MAIA2_26+MAIA2_27+MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31+MAIA2_32+
            MAIA2_33+MAIA2_34+MAIA2_35+MAIA2_36+MAIA2_37
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.881       0.679
#Tucker-Lewis Index (TLI)                       0.874       0.660
#Robust Comparative Fit Index (CFI)                         0.503
#Robust Tucker-Lewis Index (TLI)                            0.474
#RMSEA                                          0.194       0.147
#Robust RMSEA                                               0.143
#SRMR                                           0.132       0.132

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .404

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.561       0.573
#Tucker-Lewis Index (TLI)                       0.535       0.548
#Robust Comparative Fit Index (CFI)                         0.576
#Robust Tucker-Lewis Index (TLI)                            0.551
#RMSEA                                          0.121       0.103
#Robust RMSEA                                               0.117
#SRMR                                           0.116       0.116

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .314

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 8 factors (theory based) 
EGAmodel= '
 noticing    =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04
 distracting =~ MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+MAIA2_09+MAIA2_10
 worrying    =~ MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15
 attention   =~ MAIA2_16+MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22
 emotion     =~ MAIA2_23+MAIA2_24+MAIA2_25+MAIA2_26+MAIA2_27
 self_regul  =~ MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31
 body_listen =~ MAIA2_32+MAIA2_33+MAIA2_34
 trusting    =~ MAIA2_35+MAIA2_36+MAIA2_37
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.855
#Tucker-Lewis Index (TLI)                       0.951       0.839
#Robust Comparative Fit Index (CFI)                         0.752
#Robust Tucker-Lewis Index (TLI)                            0.726
#RMSEA                                          0.121       0.101
#Robust RMSEA                                               0.103
#SRMR                                           0.099       0.099

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .624

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting      -0.419    0.056   -7.547    0.000   -0.419   -0.419
#    worrying         -0.473    0.047  -10.095    0.000   -0.473   -0.473
#    attention         0.580    0.044   13.322    0.000    0.580    0.580
#    emotion           0.763    0.032   23.540    0.000    0.763    0.763
#    self_regul        0.497    0.045   11.084    0.000    0.497    0.497
#    body_listen       0.672    0.038   17.765    0.000    0.672    0.672
#    trusting          0.478    0.046   10.445    0.000    0.478    0.478
#  distracting ~~                                                        
#    worrying          0.363    0.045    8.079    0.000    0.363    0.363
#    attention        -0.512    0.041  -12.401    0.000   -0.512   -0.512
#    emotion          -0.373    0.049   -7.648    0.000   -0.373   -0.373
#    self_regul       -0.415    0.045   -9.298    0.000   -0.415   -0.415
#    body_listen      -0.347    0.048   -7.187    0.000   -0.347   -0.347
#    trusting         -0.357    0.044   -8.071    0.000   -0.357   -0.357
#  worrying ~~                                                           
#    attention        -0.267    0.045   -5.944    0.000   -0.267   -0.267
#    emotion          -0.456    0.044  -10.387    0.000   -0.456   -0.456
#    self_regul       -0.189    0.048   -3.929    0.000   -0.189   -0.189
#    body_listen      -0.394    0.049   -8.112    0.000   -0.394   -0.394
#    trusting         -0.218    0.045   -4.857    0.000   -0.218   -0.218
#  attention ~~                                                          
#    emotion           0.516    0.041   12.507    0.000    0.516    0.516
#    self_regul        0.743    0.026   29.003    0.000    0.743    0.743
#    body_listen       0.613    0.035   17.630    0.000    0.613    0.613
#    trusting          0.678    0.028   24.661    0.000    0.678    0.678
#  emotion ~~                                                            
#    self_regul        0.481    0.042   11.361    0.000    0.481    0.481
#    body_listen       0.753    0.028   26.565    0.000    0.753    0.753
#    trusting          0.429    0.045    9.470    0.000    0.429    0.429
#  self_regul ~~                                                         
#    body_listen       0.682    0.030   22.696    0.000    0.682    0.682
#    trusting          0.758    0.026   28.843    0.000    0.758    0.758
#  body_listen ~~                                                        
#    trusting          0.596    0.035   17.015    0.000    0.596    0.596


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.854       0.870
#Tucker-Lewis Index (TLI)                       0.838       0.855
#Robust Comparative Fit Index (CFI)                         0.873
#Robust Tucker-Lewis Index (TLI)                            0.860
#RMSEA                                          0.071       0.058
#Robust RMSEA                                               0.065
#SRMR                                           0.105       0.105

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .599

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting      -0.384    0.086   -4.482    0.000   -0.384   -0.384
#    worrying         -0.458    0.070   -6.566    0.000   -0.458   -0.458
#    attention         0.572    0.067    8.543    0.000    0.572    0.572
#    emotion           0.773    0.046   16.666    0.000    0.773    0.773
#    self_regul        0.453    0.073    6.226    0.000    0.453    0.453
#    body_listen       0.661    0.064   10.295    0.000    0.661    0.661
#    trusting          0.423    0.076    5.584    0.000    0.423    0.423
#  distracting ~~                                                        
#    worrying          0.166    0.085    1.954    0.051    0.166    0.166
#    attention        -0.469    0.074   -6.328    0.000   -0.469   -0.469
#    emotion          -0.352    0.077   -4.575    0.000   -0.352   -0.352
#    self_regul       -0.411    0.064   -6.442    0.000   -0.411   -0.411
#    body_listen      -0.302    0.081   -3.722    0.000   -0.302   -0.302
#    trusting         -0.318    0.071   -4.462    0.000   -0.318   -0.318
#  worrying ~~                                                           
#    attention        -0.028    0.080   -0.347    0.729   -0.028   -0.028
#    emotion          -0.384    0.069   -5.560    0.000   -0.384   -0.384
#    self_regul        0.020    0.077    0.257    0.797    0.020    0.020
#    body_listen      -0.348    0.073   -4.796    0.000   -0.348   -0.348
#    trusting          0.002    0.077    0.021    0.983    0.002    0.002
#  attention ~~                                                          
#    emotion           0.508    0.060    8.464    0.000    0.508    0.508
#    self_regul        0.743    0.037   19.874    0.000    0.743    0.743
#    body_listen       0.592    0.059    9.975    0.000    0.592    0.592
#    trusting          0.677    0.040   16.956    0.000    0.677    0.677
#  emotion ~~                                                            
#    self_regul        0.469    0.060    7.815    0.000    0.469    0.469
#    body_listen       0.718    0.056   12.795    0.000    0.718    0.718
#    trusting          0.407    0.064    6.325    0.000    0.407    0.407
#  self_regul ~~                                                         
#    body_listen       0.689    0.048   14.372    0.000    0.689    0.689
#    trusting          0.764    0.037   20.657    0.000    0.764    0.764
#  body_listen ~~                                                        
#    trusting          0.581    0.060    9.708    0.000    0.581    0.581

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.699       0.707
#Tucker-Lewis Index (TLI)                       0.681       0.690
#Robust Comparative Fit Index (CFI)                         0.715
#Robust Tucker-Lewis Index (TLI)                            0.698
#RMSEA                                          0.100       0.085
#Robust RMSEA                                               0.096
#SRMR                                           0.273       0.273

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .581


# Bifactor model with 8 factors
BIFmodel= '
 noticing    =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04
 distracting =~ MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+MAIA2_09+MAIA2_10
 worrying    =~ MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15
 attention   =~ MAIA2_16+MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22
 emotion     =~ MAIA2_23+MAIA2_24+MAIA2_25+MAIA2_26+MAIA2_27
 self_regul  =~ MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31
 body_listen =~ MAIA2_32+MAIA2_33+MAIA2_34
 trusting    =~ MAIA2_35+MAIA2_36+MAIA2_37
 general =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04+MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+
            MAIA2_09+MAIA2_10+MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15+MAIA2_16+
            MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22+MAIA2_23+MAIA2_24+
            MAIA2_25+MAIA2_26+MAIA2_27+MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31+MAIA2_32+
            MAIA2_33+MAIA2_34+MAIA2_35+MAIA2_36+MAIA2_37
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.862       0.877
#Tucker-Lewis Index (TLI)                       0.845       0.862
#Robust Comparative Fit Index (CFI)                         0.881
#Robust Tucker-Lewis Index (TLI)                            0.866
#RMSEA                                          0.070       0.057
#Robust RMSEA                                               0.064
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .590

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5183846      0.8888889      0.9293710      0.7425022 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#noticing    0.5396456 0.04891058 0.4603544 0.7567128 0.3996020 0.6293424 0.8126446
#distracting 0.7181999 0.09608557 0.2818001 0.8214203 0.5915097 0.7856007 0.8988870
#worrying    0.8368671 0.08647650 0.1631329 0.6757503 0.6707195 0.8238163 0.9113855
#attention   0.2953310 0.06173037 0.7046690 0.9159466 0.2658496 0.6154029 0.8027629
#emotion     0.5325080 0.07433419 0.4674920 0.8669192 0.4609017 0.7000545 0.8662368
#self_regul  0.2984533 0.03749971 0.7015467 0.8776916 0.2588111 0.4931059 0.7689142
#body_listen 0.4070543 0.03442270 0.5929457 0.7994980 0.3179896 0.4945344 0.7772532
#trusting    0.3716802 0.04215574 0.6283198 0.9089503 0.3198178 0.5989651 0.9111861
#general     0.5183846 0.51838465 0.5183846 0.9293710 0.7425022 0.9467646 0.9491958









################################################################
### 
### Desmedt et al (2022) https://www.sciencedirect.com/science/article/pii/S030105112200031X?via%3Dihub
### data Study https://osf.io/e2ax7/files/osfstorage

library(readxl)
MAIA_Desmedt <- read.csv("MAIA_Desmedt.csv", sep=";")
colnames(MAIA_Desmedt)
mydata <- MAIA_Desmedt[,27:63]
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 6

library(psych)
omega(mydata) # alpha = .94, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 3 components
# Eigenvalue 1 = 11.42; eigenvalue 2 = 2.51

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 4 components
# Eigenvalue 1 = 12.35; eigenvalue 2 = 2.72

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.31, RMSEA=.092, RMSR=.09, TLI=.662

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.103, RMSR=.1, TLI=.636

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 8 factors (theory-based)
fit4 <- fa(mydata,8)
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.033, RMSR=.02, TLI=.955
#      MR1   MR5   MR6   MR4   MR3   MR2   MR8   MR7
#MR1  1.00  0.52  0.44  0.36 -0.34 -0.37  0.02  0.01
#MR5  0.52  1.00  0.58  0.51 -0.36 -0.10 -0.15  0.13
#MR6  0.44  0.58  1.00  0.48 -0.33 -0.03 -0.08  0.07
#MR4  0.36  0.51  0.48  1.00 -0.31  0.00 -0.13  0.14
#MR3 -0.34 -0.36 -0.33 -0.31  1.00  0.13  0.46 -0.04
#MR2 -0.37 -0.10 -0.03  0.00  0.13  1.00 -0.17  0.18
#MR8  0.02 -0.15 -0.08 -0.13  0.46 -0.17  1.00 -0.05
#MR7  0.01  0.13  0.07  0.14 -0.04  0.18 -0.05  1.00

fit4 <- fa(rho,8,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.044, RMSR=.02, TLI=.933
#      MR1   MR5   MR3   MR6   MR4   MR2   MR8   MR7
#MR1  1.00  0.53 -0.31  0.46  0.37 -0.37  0.00  0.11
#MR5  0.53  1.00 -0.35  0.56  0.49 -0.15  0.21  0.07
#MR3 -0.31 -0.35  1.00 -0.29 -0.29  0.09 -0.36  0.22
#MR6  0.46  0.56 -0.29  1.00  0.48 -0.05  0.13  0.11
#MR4  0.37  0.49 -0.29  0.48  1.00 -0.03  0.20  0.08
#MR2 -0.37 -0.15  0.09 -0.05 -0.03  1.00  0.25 -0.06
#MR8  0.00  0.21 -0.36  0.13  0.20  0.25  1.00 -0.15
#MR7  0.11  0.07  0.22  0.11  0.08 -0.06 -0.15  1.00


# Single factor model lavaan
UNImodel= '
 general =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04+MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+
            MAIA2_09+MAIA2_10+MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15+MAIA2_16+
            MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22+MAIA2_23+MAIA2_24+
            MAIA2_25+MAIA2_26+MAIA2_27+MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31+MAIA2_32+
            MAIA2_33+MAIA2_34+MAIA2_35+MAIA2_36+MAIA2_37
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.920       0.680
#Tucker-Lewis Index (TLI)                       0.916       0.661
#Robust Comparative Fit Index (CFI)                         0.655
#Robust Tucker-Lewis Index (TLI)                            0.635
#RMSEA                                          0.133       0.123
#Robust RMSEA                                               0.105
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .397

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.680       0.682
#Tucker-Lewis Index (TLI)                       0.662       0.664
#Robust Comparative Fit Index (CFI)                         0.685
#Robust Tucker-Lewis Index (TLI)                            0.666
#RMSEA                                          0.093       0.085
#Robust RMSEA                                               0.092
#SRMR                                           0.091       0.091

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .340

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 8 factors (theory based) 
EGAmodel= '
 noticing    =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04
 distracting =~ MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+MAIA2_09+MAIA2_10
 worrying    =~ MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15
 attention   =~ MAIA2_16+MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22
 emotion     =~ MAIA2_23+MAIA2_24+MAIA2_25+MAIA2_26+MAIA2_27
 self_regul  =~ MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31
 body_listen =~ MAIA2_32+MAIA2_33+MAIA2_34
 trusting    =~ MAIA2_35+MAIA2_36+MAIA2_37
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.959       0.816
#Tucker-Lewis Index (TLI)                       0.954       0.796
#Robust Comparative Fit Index (CFI)                         0.800
#Robust Tucker-Lewis Index (TLI)                            0.778
#RMSEA                                          0.098       0.096
#Robust RMSEA                                               0.082
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .529

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting      -0.416    0.034  -12.136    0.000   -0.416   -0.416
#    worrying          0.733    0.032   23.075    0.000    0.733    0.733
#    attention         0.822    0.019   44.088    0.000    0.822    0.822
#    emotion           0.947    0.015   63.130    0.000    0.947    0.947
#    self_regul        0.734    0.023   31.350    0.000    0.734    0.734
#    body_listen       0.902    0.018   50.140    0.000    0.902    0.902
#    trusting          0.700    0.025   28.549    0.000    0.700    0.700
#  distracting ~~                                                        
#    worrying         -0.699    0.028  -24.940    0.000   -0.699   -0.699
#    attention        -0.478    0.030  -15.984    0.000   -0.478   -0.478
#    emotion          -0.373    0.033  -11.407    0.000   -0.373   -0.373
#    self_regul       -0.412    0.030  -13.563    0.000   -0.412   -0.412
#    body_listen      -0.342    0.034  -10.059    0.000   -0.342   -0.342
#    trusting         -0.354    0.032  -11.027    0.000   -0.354   -0.354
#  worrying ~~                                                           
#    attention         0.818    0.025   33.117    0.000    0.818    0.818
#    emotion           0.593    0.034   17.436    0.000    0.593    0.593
#    self_regul        0.716    0.029   24.851    0.000    0.716    0.716
#    body_listen       0.644    0.030   21.293    0.000    0.644    0.644
#    trusting          0.667    0.030   22.039    0.000    0.667    0.667
#  attention ~~                                                          
#    emotion           0.743    0.019   38.432    0.000    0.743    0.743
#    self_regul        0.865    0.013   67.411    0.000    0.865    0.865
#    body_listen       0.840    0.017   48.748    0.000    0.840    0.840
#    trusting          0.697    0.020   34.793    0.000    0.697    0.697
#  emotion ~~                                                            
#    self_regul        0.719    0.020   36.118    0.000    0.719    0.719
#    body_listen       0.910    0.014   64.537    0.000    0.910    0.910
#    trusting          0.589    0.026   22.946    0.000    0.589    0.589
#  self_regul ~~                                                         
#    body_listen       0.770    0.020   38.735    0.000    0.770    0.770
#    trusting          0.694    0.021   32.456    0.000    0.694    0.694
#  body_listen ~~                                                        
#    trusting          0.662    0.023   28.223    0.000    0.662    0.662

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.854       0.855
#Tucker-Lewis Index (TLI)                       0.838       0.840
#Robust Comparative Fit Index (CFI)                         0.858
#Robust Tucker-Lewis Index (TLI)                            0.843
#RMSEA                                          0.064       0.059
#Robust RMSEA                                               0.063
#SRMR                                           0.107       0.107

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .469

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting      -0.409    0.050   -8.105    0.000   -0.409   -0.409
#    worrying         -0.380    0.074   -5.166    0.000   -0.380   -0.380
#    attention         0.818    0.026   30.931    0.000    0.818    0.818
#    emotion           0.944    0.021   44.438    0.000    0.944    0.944
#    self_regul        0.729    0.032   22.555    0.000    0.729    0.729
#    body_listen       0.895    0.025   35.703    0.000    0.895    0.895
#    trusting          0.671    0.045   14.952    0.000    0.671    0.671
#  distracting ~~                                                        
#    worrying         -0.020    0.087   -0.230    0.818   -0.020   -0.020
#    attention        -0.463    0.043  -10.843    0.000   -0.463   -0.463
#    emotion          -0.354    0.049   -7.285    0.000   -0.354   -0.354
#    self_regul       -0.401    0.045   -8.856    0.000   -0.401   -0.401
#    body_listen      -0.323    0.049   -6.569    0.000   -0.323   -0.323
#    trusting         -0.343    0.049   -7.031    0.000   -0.343   -0.343
#  worrying ~~                                                           
#    attention        -0.087    0.084   -1.027    0.304   -0.087   -0.087
#    emotion          -0.356    0.062   -5.720    0.000   -0.356   -0.356
#    self_regul        0.009    0.079    0.119    0.905    0.009    0.009
#    body_listen      -0.325    0.065   -4.979    0.000   -0.325   -0.325
#    trusting          0.030    0.080    0.376    0.707    0.030    0.030
#  attention ~~                                                          
#    emotion           0.742    0.025   29.228    0.000    0.742    0.742
#    self_regul        0.863    0.019   44.810    0.000    0.863    0.863
#    body_listen       0.835    0.024   34.584    0.000    0.835    0.835
#    trusting          0.696    0.033   21.001    0.000    0.696    0.696
#  emotion ~~                                                            
#    self_regul        0.710    0.030   23.875    0.000    0.710    0.710
#    body_listen       0.903    0.019   47.429    0.000    0.903    0.903
#    trusting          0.568    0.042   13.509    0.000    0.568    0.568
#  self_regul ~~                                                         
#    body_listen       0.764    0.029   26.383    0.000    0.764    0.764
#    trusting          0.697    0.033   21.042    0.000    0.697    0.697
#  body_listen ~~                                                        
#    trusting          0.602    0.045   13.418    0.000    0.602    0.602

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.617       0.615
#Tucker-Lewis Index (TLI)                       0.594       0.593
#Robust Comparative Fit Index (CFI)                         0.620
#Robust Tucker-Lewis Index (TLI)                            0.598
#RMSEA                                          0.101       0.094
#Robust RMSEA                                               0.101
#SRMR                                           0.276       0.276

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .459


# Bifactor model with 8 factors
BIFmodel= '
 noticing    =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04
 distracting =~ MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+MAIA2_09+MAIA2_10
 worrying    =~ MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15
 attention   =~ MAIA2_16+MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22
 emotion     =~ MAIA2_23+MAIA2_24+MAIA2_25+MAIA2_26+MAIA2_27
 self_regul  =~ MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31
 body_listen =~ MAIA2_32+MAIA2_33+MAIA2_34
 trusting    =~ MAIA2_35+MAIA2_36+MAIA2_37
 general =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04+MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+
            MAIA2_09+MAIA2_10+MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15+MAIA2_16+
            MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22+MAIA2_23+MAIA2_24+
            MAIA2_25+MAIA2_26+MAIA2_27+MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31+MAIA2_32+
            MAIA2_33+MAIA2_34+MAIA2_35+MAIA2_36+MAIA2_37
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
#warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.896
#Tucker-Lewis Index (TLI)                       0.880       0.883
#Robust Comparative Fit Index (CFI)                         0.898
#Robust Tucker-Lewis Index (TLI)                            0.886
#RMSEA                                          0.055       0.050
#Robust RMSEA                                               0.054
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .473

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5656158      0.8888889      0.9308232      0.7914270 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#noticing    0.4074208 0.04704936 0.5925792 0.7953215 0.1829612 0.8809922 1.1039504
#distracting 0.7694612 0.10680341 0.2305388 0.8238150 0.6508627 0.7782083 0.8887511
#worrying    0.7688211 0.08765324 0.2311789 0.7526909 0.7521683 0.7255985 0.8700044
#attention   0.1729994 0.03034619 0.8270006 0.8673025 0.1497351 0.3939159 0.6637797
#emotion     0.2589744 0.03227377 0.7410256 0.8204337 0.2034427 0.4266975 0.7044714
#self_regul  0.2965075 0.03425776 0.7034925 0.8307199 0.2227750 0.4697195 0.7664642
#body_listen 0.4156340 0.04961797 0.5843660 0.8901114 0.2322955 0.9379802 1.2823611
#trusting    0.4793755 0.04638245 0.5206245 0.8262447 0.3784026 0.6050379 0.8395588
#general     0.5656158 0.56561585 0.5656158 0.9308232 0.7914270 0.9474579 0.9620201








################################################################
### 
### Ferentzi et al (2020) https://www.tandfonline.com/doi/full/10.1080/00223891.2020.1813147
### data https://osf.io/9mhry

library(readxl)
MAIA_Ferentzi <- read_excel("MAIA_Ferentzi.xlsx")
colnames(MAIA_Ferentzi)
mydata <- MAIA_Ferentzi[,9:40]
mydata <- na.omit(mydata)
colnames(mydata)
# Is MAIA 1 (32 questions)
Names2 <- colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 6

library(psych)
omega(mydata) # alpha = .93, omega T = .94

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 6 components
# Eigenvalue 1 = 10.19; eigenvalue 2 = 2.07

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 6 components
# Eigenvalue 1 = 11.13; eigenvalue 2 = 2.25

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.114, RMSR=.09, TLI=.6

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.129, RMSR=.1, TLI=.569

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 7 communities

# Give solution with 8 factors (theory-based)
fit4 <- fa(mydata,8)
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.045, RMSR=.02, TLI=.938
#     MR1   MR4  MR5  MR7  MR6   MR8   MR2   MR3
#MR1 1.00  0.41 0.50 0.56 0.59  0.44  0.22  0.08
#MR4 0.41  1.00 0.30 0.54 0.38  0.55 -0.03  0.08
#MR5 0.50  0.30 1.00 0.44 0.47  0.25  0.25  0.16
#MR7 0.56  0.54 0.44 1.00 0.52  0.51  0.00  0.17
#MR6 0.59  0.38 0.47 0.52 1.00  0.25  0.22  0.05
#MR8 0.44  0.55 0.25 0.51 0.25  1.00 -0.05  0.11
#MR2 0.22 -0.03 0.25 0.00 0.22 -0.05  1.00 -0.01
#MR3 0.08  0.08 0.16 0.17 0.05  0.11 -0.01  1.00

fit4 <- fa(rho,8,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.056, RMSR=.02, TLI=.92
#     MR1   MR4  MR5   MR7  MR6   MR8   MR2   MR3
#MR1 1.00  0.43 0.51  0.55 0.60  0.46  0.22  0.07
#MR4 0.43  1.00 0.32  0.55 0.39  0.58 -0.02  0.07
#MR5 0.51  0.32 1.00  0.43 0.47  0.26  0.26  0.16
#MR7 0.55  0.55 0.43  1.00 0.50  0.51 -0.01  0.17
#MR6 0.60  0.39 0.47  0.50 1.00  0.24  0.23  0.04
#MR8 0.46  0.58 0.26  0.51 0.24  1.00 -0.04  0.11
#MR2 0.22 -0.02 0.26 -0.01 0.23 -0.04  1.00 -0.01
#MR3 0.07  0.07 0.16  0.17 0.04  0.11 -0.01  1.00


# Single factor model lavaan
UNImodel= '
 general =~ maia_1+maia_2+maia_3+maia_4+maia_5+maia_6+maia_7+maia_8+maia_9+
            maia_10+maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17+
            maia_18+maia_19+maia_20+maia_21+maia_22+maia_23+maia_24+maia_25+
            maia_26+maia_27+maia_28+maia_29+maia_30+maia_31+maia_32
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.915       0.695
#Tucker-Lewis Index (TLI)                       0.910       0.674
#Robust Comparative Fit Index (CFI)                         0.583
#Robust Tucker-Lewis Index (TLI)                            0.554
#RMSEA                                          0.155       0.155
#Robust RMSEA                                               0.132
#SRMR                                           0.106       0.106

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .431

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.626       0.626
#Tucker-Lewis Index (TLI)                       0.600       0.601
#Robust Comparative Fit Index (CFI)                         0.628
#Robust Tucker-Lewis Index (TLI)                            0.603
#RMSEA                                          0.114       0.101
#Robust RMSEA                                               0.114
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .357

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 8 factors (theory based) 
EGAmodel= '
 noticing    =~ maia_1+maia_2+maia_3+maia_4
 distracting =~ maia_5+maia_6+maia_7
 worrying    =~ maia_8+maia_9+maia_10
 attention   =~ maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17
 emotion     =~ maia_18+maia_19+maia_20+maia_21+maia_22
 self_regul  =~ maia_23+maia_24+maia_25+maia_26
 body_listen =~ maia_27+maia_28+maia_29
 trusting    =~ maia_30+maia_31+maia_32
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.983       0.924
#Tucker-Lewis Index (TLI)                       0.980       0.913
#Robust Comparative Fit Index (CFI)                         0.863
#Robust Tucker-Lewis Index (TLI)                            0.844
#RMSEA                                          0.073       0.080
#Robust RMSEA                                               0.078
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .62

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting       0.420    0.030   14.161    0.000    0.420    0.420
#    worrying          0.073    0.027    2.760    0.006    0.073    0.073
#    attention         0.678    0.016   41.634    0.000    0.678    0.678
#    emotion           0.790    0.013   58.932    0.000    0.790    0.790
#    self_regul        0.584    0.019   30.331    0.000    0.584    0.584
#    body_listen       0.701    0.017   41.693    0.000    0.701    0.701
#    trusting          0.528    0.019   27.776    0.000    0.528    0.528
#  distracting ~~                                                        
#    worrying          0.026    0.024    1.096    0.273    0.026    0.026
#    attention         0.298    0.026   11.516    0.000    0.298    0.298
#    emotion           0.340    0.027   12.514    0.000    0.340    0.340
#    self_regul        0.260    0.026   10.052    0.000    0.260    0.260
#    body_listen       0.374    0.027   14.009    0.000    0.374    0.374
#    trusting          0.325    0.025   12.905    0.000    0.325    0.325
#  worrying ~~                                                           
#    attention         0.310    0.022   14.357    0.000    0.310    0.310
#    emotion          -0.021    0.024   -0.858    0.391   -0.021   -0.021
#    self_regul        0.296    0.022   13.246    0.000    0.296    0.296
#    body_listen       0.046    0.024    1.891    0.059    0.046    0.046
#    trusting          0.289    0.022   13.048    0.000    0.289    0.289
#  attention ~~                                                          
#    emotion           0.565    0.016   36.003    0.000    0.565    0.565
#    self_regul        0.747    0.011   68.758    0.000    0.747    0.747
#    body_listen       0.684    0.013   53.366    0.000    0.684    0.684
#    trusting          0.637    0.014   46.224    0.000    0.637    0.637
#  emotion ~~                                                            
#    self_regul        0.542    0.016   34.181    0.000    0.542    0.542
#    body_listen       0.707    0.013   56.372    0.000    0.707    0.707
#    trusting          0.423    0.019   22.717    0.000    0.423    0.423
#  self_regul ~~                                                         
#    body_listen       0.680    0.013   53.109    0.000    0.680    0.680
#    trusting          0.635    0.014   45.195    0.000    0.635    0.635
#  body_listen ~~                                                        
#    trusting          0.590    0.015   39.398    0.000    0.590    0.590

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.900       0.901
#Tucker-Lewis Index (TLI)                       0.886       0.888
#Robust Comparative Fit Index (CFI)                         0.903
#Robust Tucker-Lewis Index (TLI)                            0.890
#RMSEA                                          0.061       0.054
#Robust RMSEA                                               0.060
#SRMR                                           0.070       0.070

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .559

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting       0.140    0.045    3.090    0.002    0.140    0.140
#    worrying          0.009    0.032    0.271    0.786    0.009    0.009
#    attention         0.677    0.021   32.854    0.000    0.677    0.677
#    emotion           0.767    0.019   39.562    0.000    0.767    0.767
#    self_regul        0.582    0.024   24.537    0.000    0.582    0.582
#    body_listen       0.692    0.020   34.458    0.000    0.692    0.692
#    trusting          0.474    0.028   17.030    0.000    0.474    0.474
#  distracting ~~                                                        
#    worrying          0.020    0.029    0.707    0.480    0.020    0.020
#    attention         0.133    0.031    4.281    0.000    0.133    0.133
#    emotion           0.128    0.035    3.666    0.000    0.128    0.128
#    self_regul        0.107    0.031    3.424    0.001    0.107    0.107
#    body_listen       0.194    0.040    4.817    0.000    0.194    0.194
#    trusting          0.198    0.031    6.303    0.000    0.198    0.198
#  worrying ~~                                                           
#    attention         0.233    0.029    8.130    0.000    0.233    0.233
#    emotion          -0.051    0.028   -1.834    0.067   -0.051   -0.051
#    self_regul        0.241    0.029    8.448    0.000    0.241    0.241
#    body_listen       0.011    0.028    0.391    0.696    0.011    0.011
#    trusting          0.247    0.026    9.465    0.000    0.247    0.247
#  attention ~~                                                          
#    emotion           0.535    0.021   24.950    0.000    0.535    0.535
#    self_regul        0.750    0.015   49.592    0.000    0.750    0.750
#    body_listen       0.683    0.017   40.952    0.000    0.683    0.683
#    trusting          0.619    0.019   32.332    0.000    0.619    0.619
#  emotion ~~                                                            
#    self_regul        0.526    0.021   25.095    0.000    0.526    0.526
#    body_listen       0.670    0.018   37.021    0.000    0.670    0.670
#    trusting          0.375    0.025   14.788    0.000    0.375    0.375
#  self_regul ~~                                                         
#    body_listen       0.687    0.018   39.020    0.000    0.687    0.687
#    trusting          0.617    0.019   32.563    0.000    0.617    0.617
#  body_listen ~~                                                        
#    trusting          0.554    0.021   26.847    0.000    0.554    0.554

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.728       0.725
#Tucker-Lewis Index (TLI)                       0.709       0.706
#Robust Comparative Fit Index (CFI)                         0.730
#Robust Tucker-Lewis Index (TLI)                            0.712
#RMSEA                                          0.098       0.087
#Robust RMSEA                                               0.097
#SRMR                                           0.271       0.271

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .504


# Bifactor model with 8 factors
BIFmodel= '
 noticing    =~ maia_1+maia_2+maia_3+maia_4
 distracting =~ maia_5+maia_6+maia_7
 worrying    =~ maia_8+maia_9+maia_10
 attention   =~ maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17
 emotion     =~ maia_18+maia_19+maia_20+maia_21+maia_22
 self_regul  =~ maia_23+maia_24+maia_25+maia_26
 body_listen =~ maia_27+maia_28+maia_29
 trusting    =~ maia_30+maia_31+maia_32
 general =~ maia_1+maia_2+maia_3+maia_4+maia_5+maia_6+maia_7+maia_8+maia_9+
            maia_10+maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17+
            maia_18+maia_19+maia_20+maia_21+maia_22+maia_23+maia_24+maia_25+
            maia_26+maia_27+maia_28+maia_29+maia_30+maia_31+maia_32
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.916       0.917
#Tucker-Lewis Index (TLI)                       0.904       0.905
#Robust Comparative Fit Index (CFI)                         0.919
#Robust Tucker-Lewis Index (TLI)                            0.907
#RMSEA                                          0.056       0.049
#Robust RMSEA                                               0.055
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .583

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5328088      0.8891129      0.9542419      0.8557407 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#noticing    0.3887905 0.03613406 0.61120950 0.7231109 0.2476707 0.4975120 0.7466942
#distracting 0.8711895 0.06610576 0.12881049 0.6677837 0.6348859 0.8079098 0.9079705
#worrying    0.9460927 0.08572076 0.05390731 0.7679900 0.7318774 0.8697435 0.9368633
#attention   0.2871868 0.05941263 0.71281323 0.8835639 0.2398834 0.5792207 0.7929243
#emotion     0.4766034 0.07463760 0.52339657 0.8575932 0.3901558 0.6745700 0.8611676
#self_regul  0.3425494 0.04638940 0.65745060 0.8543862 0.2679929 0.5435752 0.8078155
#body_listen 0.3046916 0.03585489 0.69530837 0.8745747 0.2606740 0.4615788 0.7896553
#trusting    0.5075724 0.06293609 0.49242758 0.8920663 0.4430996 0.6814213 0.8962673
#general     0.5328088 0.53280882 0.53280882 0.9542419 0.8557407 0.9393940 0.9455462








################################################################
### 
### Reis (2017) https://econtent.hogrefe.com/doi/abs/10.1027/1015-5759/a000404?journalCode=jpa
### data https://osf.io/8hpb8/

MAIA_Reis <- read.csv("MAIA_Reis.csv")
colnames(MAIA_Reis)
mydata <- MAIA_Reis[,1:32]
mydata <- na.omit(mydata)
colnames(mydata)
# Is MAIA 1 (32 questions)
colnames(mydata) <- Names2
min(mydata)
max(mydata) # response alternatives = 6

library(psych)
omega(mydata) # alpha = .95, omega T = .96

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 2 components
# Eigenvalue 1 = 12.75; eigenvalue 2 = 1.95

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 13.79; eigenvalue 2 = 2.11

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.4, RMSEA=.102, RMSR=.09, TLI=.712

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.13, RMSR=.1, TLI=.632

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 8 factors (theory-based)
fit4 <- fa(mydata,8)
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.043, RMSR=.02, TLI=.947
#      MR2  MR7  MR6  MR3  MR1  MR8   MR5   MR4
#MR2  1.00 0.40 0.37 0.30 0.28 0.34 -0.13  0.17
#MR7  0.40 1.00 0.53 0.41 0.41 0.37  0.29  0.38
#MR6  0.37 0.53 1.00 0.32 0.35 0.36  0.15  0.19
#MR3  0.30 0.41 0.32 1.00 0.35 0.26  0.23  0.20
#MR1  0.28 0.41 0.35 0.35 1.00 0.36  0.09  0.26
#MR8  0.34 0.37 0.36 0.26 0.36 1.00  0.08  0.17
#MR5 -0.13 0.29 0.15 0.23 0.09 0.08  1.00 -0.02
#MR4  0.17 0.38 0.19 0.20 0.26 0.17 -0.02  1.00

fit4 <- fa(rho,8,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.68, RMSEA=.085, RMSR=.02, TLI=.839
#     MR1   MR4  MR5   MR7  MR6   MR8   MR2   MR3
#MR1 1.00  0.43 0.51  0.55 0.60  0.46  0.22  0.07
#MR4 0.43  1.00 0.32  0.55 0.39  0.58 -0.02  0.07
#MR5 0.51  0.32 1.00  0.43 0.47  0.26  0.26  0.16
#MR7 0.55  0.55 0.43  1.00 0.50  0.51 -0.01  0.17
#MR6 0.60  0.39 0.47  0.50 1.00  0.24  0.23  0.04
#MR8 0.46  0.58 0.26  0.51 0.24  1.00 -0.04  0.11
#MR2 0.22 -0.02 0.26 -0.01 0.23 -0.04  1.00 -0.01
#MR3 0.07  0.07 0.16  0.17 0.04  0.11 -0.01  1.00


# Single factor model lavaan
UNImodel= '
 general =~ maia_1+maia_2+maia_3+maia_4+maia_5+maia_6+maia_7+maia_8+maia_9+
            maia_10+maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17+
            maia_18+maia_19+maia_20+maia_21+maia_22+maia_23+maia_24+maia_25+
            maia_26+maia_27+maia_28+maia_29+maia_30+maia_31+maia_32
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.969       0.868
#Tucker-Lewis Index (TLI)                       0.966       0.859
#Robust Comparative Fit Index (CFI)                         0.667
#Robust Tucker-Lewis Index (TLI)                            0.644
#RMSEA                                          0.126       0.117
#Robust RMSEA                                               0.137
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .517

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.724       0.728
#Tucker-Lewis Index (TLI)                       0.705       0.709
#Robust Comparative Fit Index (CFI)                         0.733
#Robust Tucker-Lewis Index (TLI)                            0.714
#RMSEA                                          0.108       0.100
#Robust RMSEA                                               0.105
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .449

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 8 factors (theory based) 
EGAmodel= '
 noticing    =~ maia_1+maia_2+maia_3+maia_4
 distracting =~ maia_5+maia_6+maia_7
 worrying    =~ maia_8+maia_9+maia_10
 attention   =~ maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17
 emotion     =~ maia_18+maia_19+maia_20+maia_21+maia_22
 self_regul  =~ maia_23+maia_24+maia_25+maia_26
 body_listen =~ maia_27+maia_28+maia_29
 trusting    =~ maia_30+maia_31+maia_32
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.961
#Tucker-Lewis Index (TLI)                       0.996       0.955
#Robust Comparative Fit Index (CFI)                         0.835
#Robust Tucker-Lewis Index (TLI)                            0.812
#RMSEA                                          0.045       0.066
#Robust RMSEA                                               0.099
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .669

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting       0.422    0.069    6.097    0.000    0.422    0.422
#    worrying          0.241    0.067    3.602    0.000    0.241    0.241
#    attention         0.766    0.041   18.675    0.000    0.766    0.766
#    emotion           0.833    0.035   23.846    0.000    0.833    0.833
#    self_regul        0.711    0.045   15.977    0.000    0.711    0.711
#    body_listen       0.686    0.049   13.979    0.000    0.686    0.686
#    trusting          0.563    0.055   10.167    0.000    0.563    0.563
#  distracting ~~                                                        
#    worrying          0.570    0.055   10.386    0.000    0.570    0.570
#    attention         0.518    0.058    8.881    0.000    0.518    0.518
#    emotion           0.382    0.065    5.923    0.000    0.382    0.382
#    self_regul        0.506    0.056    8.966    0.000    0.506    0.506
#    body_listen       0.566    0.059    9.527    0.000    0.566    0.566
#    trusting          0.507    0.059    8.565    0.000    0.507    0.507
#  worrying ~~                                                           
#    attention         0.540    0.053   10.287    0.000    0.540    0.540
#    emotion           0.169    0.073    2.316    0.021    0.169    0.169
#    self_regul        0.522    0.060    8.764    0.000    0.522    0.522
#    body_listen       0.332    0.067    4.940    0.000    0.332    0.332
#    trusting          0.565    0.051   11.117    0.000    0.565    0.565
#  attention ~~                                                          
#    emotion           0.623    0.043   14.412    0.000    0.623    0.623
#    self_regul        0.872    0.021   41.351    0.000    0.872    0.872
#    body_listen       0.818    0.028   28.958    0.000    0.818    0.818
#    trusting          0.752    0.030   24.713    0.000    0.752    0.752
#  emotion ~~                                                            
#    self_regul        0.690    0.042   16.558    0.000    0.690    0.690
#    body_listen       0.677    0.045   15.169    0.000    0.677    0.677
#    trusting          0.615    0.049   12.512    0.000    0.615    0.615
#  self_regul ~~                                                         
#    body_listen       0.854    0.026   32.297    0.000    0.854    0.854
#    trusting          0.816    0.027   29.832    0.000    0.816    0.816
#  body_listen ~~                                                        
#    trusting          0.794    0.031   25.353    0.000    0.794    0.794

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.900       0.904
#Tucker-Lewis Index (TLI)                       0.887       0.891
#Robust Comparative Fit Index (CFI)                         0.908
#Robust Tucker-Lewis Index (TLI)                            0.895
#RMSEA                                          0.067       0.061
#Robust RMSEA                                               0.064
#SRMR                                           0.065       0.065

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .605

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting       0.460    0.094    4.908    0.000    0.460    0.460
#    worrying          0.240    0.107    2.247    0.025    0.240    0.240
#    attention         0.782    0.062   12.607    0.000    0.782    0.782
#    emotion           0.830    0.048   17.225    0.000    0.830    0.830
#    self_regul        0.727    0.063   11.482    0.000    0.727    0.727
#    body_listen       0.716    0.064   11.164    0.000    0.716    0.716
#    trusting          0.578    0.086    6.725    0.000    0.578    0.578
#  distracting ~~                                                        
#    worrying          0.615    0.132    4.653    0.000    0.615    0.615
#    attention         0.547    0.086    6.376    0.000    0.547    0.547
#    emotion           0.388    0.092    4.224    0.000    0.388    0.388
#    self_regul        0.538    0.077    6.995    0.000    0.538    0.538
#    body_listen       0.609    0.080    7.579    0.000    0.609    0.609
#    trusting          0.506    0.080    6.307    0.000    0.506    0.506
#  worrying ~~                                                           
#    attention         0.544    0.085    6.394    0.000    0.544    0.544
#    emotion           0.154    0.094    1.636    0.102    0.154    0.154
#    self_regul        0.510    0.084    6.099    0.000    0.510    0.510
#    body_listen       0.335    0.102    3.283    0.001    0.335    0.335
#    trusting          0.560    0.079    7.057    0.000    0.560    0.560
#  attention ~~                                                          
#    emotion           0.599    0.059   10.069    0.000    0.599    0.599
#    self_regul        0.865    0.031   27.819    0.000    0.865    0.865
#    body_listen       0.822    0.038   21.380    0.000    0.822    0.822
#    trusting          0.738    0.051   14.459    0.000    0.738    0.738
#  emotion ~~                                                            
#    self_regul        0.690    0.056   12.232    0.000    0.690    0.690
#    body_listen       0.669    0.056   11.887    0.000    0.669    0.669
#    trusting          0.578    0.075    7.729    0.000    0.578    0.578
#  self_regul ~~                                                         
#    body_listen       0.858    0.038   22.804    0.000    0.858    0.858
#    trusting          0.807    0.046   17.364    0.000    0.807    0.807
#  body_listen ~~                                                        
#    trusting          0.777    0.048   16.127    0.000    0.777    0.777


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.660       0.656
#Tucker-Lewis Index (TLI)                       0.637       0.632
#Robust Comparative Fit Index (CFI)                         0.666
#Robust Tucker-Lewis Index (TLI)                            0.643
#RMSEA                                          0.120       0.113
#Robust RMSEA                                               0.118
#SRMR                                           0.350       0.350

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .585


# Bifactor model with 8 factors
BIFmodel= '
 noticing    =~ maia_1+maia_2+maia_3+maia_4
 distracting =~ maia_5+maia_6+maia_7
 worrying    =~ maia_8+maia_9+maia_10
 attention   =~ maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17
 emotion     =~ maia_18+maia_19+maia_20+maia_21+maia_22
 self_regul  =~ maia_23+maia_24+maia_25+maia_26
 body_listen =~ maia_27+maia_28+maia_29
 trusting    =~ maia_30+maia_31+maia_32
 general =~ maia_1+maia_2+maia_3+maia_4+maia_5+maia_6+maia_7+maia_8+maia_9+
            maia_10+maia_11+maia_12+maia_13+maia_14+maia_15+maia_16+maia_17+
            maia_18+maia_19+maia_20+maia_21+maia_22+maia_23+maia_24+maia_25+
            maia_26+maia_27+maia_28+maia_29+maia_30+maia_31+maia_32
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.896
#Tucker-Lewis Index (TLI)                       0.878       0.881
#Robust Comparative Fit Index (CFI)                         0.901
#Robust Tucker-Lewis Index (TLI)                            0.886
#RMSEA                                          0.069       0.064
#Robust RMSEA                                               0.067
#SRMR                                           0.065       0.065

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .592

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing =~                                                           
#    maia_1            0.533    0.129    4.123    0.000    0.533    0.461
#    maia_2            0.508    0.102    4.990    0.000    0.508    0.536
#    maia_3            0.428    0.139    3.086    0.002    0.428    0.332
#    maia_4            0.489    0.133    3.680    0.000    0.489    0.401
#  distracting =~                                                        
#    maia_5            0.758    0.123    6.179    0.000    0.758    0.548
#    maia_6            0.579    0.121    4.782    0.000    0.579    0.502
#    maia_7            0.728    0.118    6.189    0.000    0.728    0.582
#  worrying =~                                                           
#    maia_8            0.669    0.135    4.956    0.000    0.669    0.435
#    maia_9            1.265    0.195    6.492    0.000    1.265    0.962
#    maia_10           0.554    0.118    4.712    0.000    0.554    0.454
#  attention =~                                                          
#    maia_11           0.499    0.150    3.336    0.001    0.499    0.381
#    maia_12           0.308    0.116    2.647    0.008    0.308    0.245
#    maia_13           0.488    0.142    3.434    0.001    0.488    0.404
#    maia_14           0.475    0.104    4.560    0.000    0.475    0.382
#    maia_15           0.375    0.127    2.964    0.003    0.375    0.315
#    maia_16           0.398    0.133    3.001    0.003    0.398    0.297
#    maia_17           0.335    0.117    2.871    0.004    0.335    0.243
#  emotion =~                                                            
#    maia_18           0.607    0.087    6.967    0.000    0.607    0.493
#    maia_19           0.544    0.095    5.733    0.000    0.544    0.453
#    maia_20           0.613    0.075    8.189    0.000    0.613    0.551
#    maia_21           0.627    0.091    6.855    0.000    0.627    0.551
#    maia_22           0.600    0.087    6.922    0.000    0.600    0.546
#  self_regul =~                                                         
#    maia_23           0.178    0.124    1.433    0.152    0.178    0.124
#    maia_24           0.239    0.134    1.788    0.074    0.239    0.188
#    maia_25           0.770    0.291    2.648    0.008    0.770    0.595
#    maia_26           0.358    0.163    2.202    0.028    0.358    0.260
#  body_listen =~                                                        
#    maia_27           0.332    0.089    3.737    0.000    0.332    0.252
#    maia_28           0.564    0.107    5.267    0.000    0.564    0.460
#    maia_29           0.523    0.119    4.417    0.000    0.523    0.392
#  trusting =~                                                           
#    maia_30           0.683    0.130    5.252    0.000    0.683    0.504
#    maia_31           0.755    0.144    5.256    0.000    0.755    0.554
#    maia_32           0.275    0.084    3.256    0.001    0.275    0.230
#  general =~                                                            
#    maia_1            0.614    0.089    6.899    0.000    0.614    0.531
#    maia_2            0.325    0.078    4.155    0.000    0.325    0.342
#    maia_3            0.875    0.084   10.441    0.000    0.875    0.678
#    maia_4            0.571    0.095    6.026    0.000    0.571    0.468
#    maia_5            0.735    0.095    7.734    0.000    0.735    0.532
#    maia_6            0.199    0.102    1.945    0.052    0.199    0.172
#    maia_7            0.496    0.099    4.998    0.000    0.496    0.397
#    maia_8            0.627    0.113    5.548    0.000    0.627    0.407
#    maia_9            0.195    0.110    1.767    0.077    0.195    0.148
#    maia_10           0.524    0.094    5.574    0.000    0.524    0.429
#    maia_11           0.904    0.084   10.757    0.000    0.904    0.690
#    maia_12           0.893    0.077   11.604    0.000    0.893    0.711
#    maia_13           0.632    0.079    8.011    0.000    0.632    0.524
#    maia_14           0.952    0.064   14.798    0.000    0.952    0.767
#    maia_15           0.864    0.072   11.955    0.000    0.864    0.725
#    maia_16           0.832    0.096    8.665    0.000    0.832    0.620
#    maia_17           1.060    0.088   12.104    0.000    1.060    0.767
#    maia_18           0.566    0.090    6.262    0.000    0.566    0.460
#    maia_19           0.523    0.101    5.196    0.000    0.523    0.436
#    maia_20           0.702    0.089    7.908    0.000    0.702    0.631
#    maia_21           0.646    0.096    6.751    0.000    0.646    0.568
#    maia_22           0.611    0.085    7.215    0.000    0.611    0.556
#    maia_23           1.079    0.085   12.681    0.000    1.079    0.753
#    maia_24           0.956    0.079   12.105    0.000    0.956    0.752
#    maia_25           0.913    0.082   11.132    0.000    0.913    0.705
#    maia_26           1.087    0.073   14.833    0.000    1.087    0.790
#    maia_27           1.094    0.067   16.384    0.000    1.094    0.831
#    maia_28           0.879    0.078   11.197    0.000    0.879    0.717
#    maia_29           0.984    0.081   12.218    0.000    0.984    0.737
#    maia_30           0.947    0.072   13.092    0.000    0.947    0.699
#    maia_31           0.998    0.076   13.127    0.000    0.998    0.732
#    maia_32           0.861    0.072   11.991    0.000    0.861    0.720

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6492131      0.8891129      0.9670360      0.9039289 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#noticing    0.4168015 0.04086569 0.5831985 0.7665206 0.3243545 0.4968098 0.7466981
#distracting 0.6548541 0.04731640 0.3451459 0.7029704 0.4832240 0.5613591 0.7886381
#worrying    0.7802272 0.07007365 0.2197728 0.7706876 0.6006977 0.9281860 0.9727149
#attention   0.1853495 0.04034622 0.8146505 0.9068369 0.1650956 0.4648680 0.7359467
#emotion     0.4859574 0.07182280 0.5140426 0.8613260 0.4213007 0.6522827 0.8544817
#self_regul  0.1732348 0.02506099 0.8267652 0.8905573 0.1171275 0.4020219 0.8092518
#body_listen 0.1970989 0.02274885 0.8029011 0.8863836 0.1678445 0.3411018 0.7180920
#trusting    0.2847198 0.03255230 0.7152802 0.8813474 0.2326746 0.4560735 0.8247706
#general     0.6492131 0.64921309 0.6492131 0.9670360 0.9039289 0.9609519 0.9675196








################################################################
### 
### Randelovic et al (2024) https://osf.io/preprints/psyarxiv/8fdma
### data https://osf.io/9fmc2?view_only=2c087da9c77a4db59d4b14e47037bc33

library(haven)
MAIA_Randelovic <- read_sav("MAIA_Randelovic.sav")
colnames(MAIA_Randelovic)
mydata <- MAIA_Randelovic[,4:40]
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 6

library(psych)
omega(mydata) # alpha = .88, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 6 components
# Eigenvalue 1 = 7.6; eigenvalue 2 = 2.86

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 6 components
# Eigenvalue 1 = 8.42; eigenvalue 2 = 3.26

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.21, RMSEA=.116, RMSR=.13, TLI=.382

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.141, RMSR=.14, TLI=.322

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities

# Give solution with 8 factors (theory-based)
fit4 <- fa(mydata,8)
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.062, RMSR=.04, TLI=.82
#      MR7   MR2   MR3   MR1   MR4   MR5   MR8   MR6
#MR7  1.00  0.14 -0.09  0.38  0.25  0.27  0.26  0.12
#MR2  0.14  1.00 -0.02 -0.08 -0.17 -0.06 -0.15 -0.07
#MR3 -0.09 -0.02  1.00 -0.03 -0.05  0.00  0.13 -0.09
#MR1  0.38 -0.08 -0.03  1.00  0.29  0.32  0.32  0.08
#MR4  0.25 -0.17 -0.05  0.29  1.00  0.38  0.25  0.20
#MR5  0.27 -0.06  0.00  0.32  0.38  1.00  0.38  0.12
#MR8  0.26 -0.15  0.13  0.32  0.25  0.38  1.00  0.09
#MR6  0.12 -0.07 -0.09  0.08  0.20  0.12  0.09  1.00

fit4 <- fa(rho,8,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.093, RMSR=.04, TLI=.697
#      MR2   MR3   MR1   MR7   MR4   MR5  MR6   MR8
#MR2  1.00 -0.05 -0.10  0.16 -0.17 -0.07 0.05 -0.13
#MR3 -0.05  1.00 -0.10 -0.09 -0.05  0.04 0.05  0.12
#MR1 -0.10 -0.10  1.00  0.34  0.37  0.31 0.20  0.25
#MR7  0.16 -0.09  0.34  1.00  0.24  0.21 0.28  0.19
#MR4 -0.17 -0.05  0.37  0.24  1.00  0.40 0.17  0.23
#MR5 -0.07  0.04  0.31  0.21  0.40  1.00 0.27  0.35
#MR6  0.05  0.05  0.20  0.28  0.17  0.27 1.00  0.15
#MR8 -0.13  0.12  0.25  0.19  0.23  0.35 0.15  1.00


# Single factor model lavaan
UNImodel= '
 general =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04+MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+
            MAIA2_09+MAIA2_10+MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15+MAIA2_16+
            MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22+MAIA2_23+MAIA2_24+
            MAIA2_25+MAIA2_26+MAIA2_27+MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31+MAIA2_32+
            MAIA2_33+MAIA2_34+MAIA2_35+MAIA2_36+MAIA2_37
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.702       0.423
#Tucker-Lewis Index (TLI)                       0.684       0.389
#Robust Comparative Fit Index (CFI)                         0.364
#Robust Tucker-Lewis Index (TLI)                            0.327
#RMSEA                                          0.210       0.152
#Robust RMSEA                                               0.148
#SRMR                                           0.146       0.146

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .302

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.412       0.408
#Tucker-Lewis Index (TLI)                       0.377       0.373
#Robust Comparative Fit Index (CFI)                         0.417
#Robust Tucker-Lewis Index (TLI)                            0.383
#RMSEA                                          0.122       0.114
#Robust RMSEA                                               0.120
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .246

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 8 factors (theory based) 
EGAmodel= '
 noticing    =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04
 distracting =~ MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+MAIA2_09+MAIA2_10
 worrying    =~ MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15
 attention   =~ MAIA2_16+MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22
 emotion     =~ MAIA2_23+MAIA2_24+MAIA2_25+MAIA2_26+MAIA2_27
 self_regul  =~ MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31
 body_listen =~ MAIA2_32+MAIA2_33+MAIA2_34
 trusting    =~ MAIA2_35+MAIA2_36+MAIA2_37
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
#warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.961       0.894
#Tucker-Lewis Index (TLI)                       0.957       0.883
#Robust Comparative Fit Index (CFI)                         0.751
#Robust Tucker-Lewis Index (TLI)                            0.724
#RMSEA                                          0.077       0.067
#Robust RMSEA                                               0.095
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .571

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting       0.028    0.075    0.378    0.705    0.028    0.028
#    worrying          0.101    0.073    1.392    0.164    0.101    0.101
#    attention         0.534    0.054    9.879    0.000    0.534    0.534
#    emotion           0.726    0.052   13.934    0.000    0.726    0.726
#    self_regul        0.416    0.060    6.883    0.000    0.416    0.416
#    body_listen       0.641    0.055   11.666    0.000    0.641    0.641
#    trusting          0.364    0.067    5.394    0.000    0.364    0.364
#  distracting ~~                                                        
#    worrying         -0.087    0.069   -1.265    0.206   -0.087   -0.087
#    attention         0.077    0.065    1.177    0.239    0.077    0.077
#    emotion          -0.063    0.067   -0.936    0.350   -0.063   -0.063
#    self_regul       -0.176    0.062   -2.859    0.004   -0.176   -0.176
#    body_listen      -0.043    0.072   -0.592    0.554   -0.043   -0.043
#    trusting         -0.074    0.064   -1.155    0.248   -0.074   -0.074
#  worrying ~~                                                           
#    attention        -0.192    0.064   -2.998    0.003   -0.192   -0.192
#    emotion           0.133    0.071    1.871    0.061    0.133    0.133
#    self_regul       -0.202    0.065   -3.118    0.002   -0.202   -0.202
#    body_listen      -0.122    0.070   -1.733    0.083   -0.122   -0.122
#    trusting         -0.216    0.064   -3.364    0.001   -0.216   -0.216
#  attention ~~                                                          
#    emotion           0.445    0.051    8.717    0.000    0.445    0.445
#    self_regul        0.590    0.040   14.873    0.000    0.590    0.590
#    body_listen       0.601    0.041   14.711    0.000    0.601    0.601
#    trusting          0.566    0.043   13.183    0.000    0.566    0.566
#  emotion ~~                                                            
#    self_regul        0.528    0.044   12.125    0.000    0.528    0.528
#    body_listen       0.688    0.045   15.328    0.000    0.688    0.688
#    trusting          0.421    0.056    7.491    0.000    0.421    0.421
#  self_regul ~~                                                         
#    body_listen       0.657    0.042   15.731    0.000    0.657    0.657
#    trusting          0.498    0.049   10.145    0.000    0.498    0.498
#  body_listen ~~                                                        
#    trusting          0.473    0.055    8.532    0.000    0.473    0.473

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.807       0.810
#Tucker-Lewis Index (TLI)                       0.786       0.790
#Robust Comparative Fit Index (CFI)                         0.817
#Robust Tucker-Lewis Index (TLI)                            0.797
#RMSEA                                          0.071       0.066
#Robust RMSEA                                               0.069
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .483

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing ~~                                                           
#    distracting       0.021    0.116    0.180    0.857    0.021    0.021
#    worrying          0.074    0.096    0.769    0.442    0.074    0.074
#    attention         0.536    0.082    6.527    0.000    0.536    0.536
#    emotion           0.673    0.095    7.102    0.000    0.673    0.673
#    self_regul        0.444    0.081    5.493    0.000    0.444    0.444
#    body_listen       0.629    0.079    7.959    0.000    0.629    0.629
#    trusting          0.285    0.084    3.391    0.001    0.285    0.285
#  distracting ~~                                                        
#    worrying         -0.056    0.125   -0.445    0.657   -0.056   -0.056
#    attention         0.090    0.120    0.748    0.455    0.090    0.090
#    emotion          -0.060    0.121   -0.491    0.623   -0.060   -0.060
#    self_regul       -0.185    0.109   -1.699    0.089   -0.185   -0.185
#    body_listen      -0.030    0.111   -0.266    0.790   -0.030   -0.030
#    trusting         -0.073    0.111   -0.656    0.512   -0.073   -0.073
#  worrying ~~                                                           
#    attention        -0.199    0.088   -2.260    0.024   -0.199   -0.199
#    emotion           0.121    0.084    1.438    0.150    0.121    0.121
#    self_regul       -0.203    0.084   -2.406    0.016   -0.203   -0.203
#    body_listen      -0.098    0.084   -1.173    0.241   -0.098   -0.098
#    trusting         -0.225    0.075   -3.017    0.003   -0.225   -0.225
#  attention ~~                                                          
#    emotion           0.443    0.079    5.600    0.000    0.443    0.443
#    self_regul        0.623    0.083    7.468    0.000    0.623    0.623
#    body_listen       0.615    0.064    9.603    0.000    0.615    0.615
#    trusting          0.513    0.070    7.364    0.000    0.513    0.513
#  emotion ~~                                                            
#    self_regul        0.568    0.066    8.548    0.000    0.568    0.568
#    body_listen       0.658    0.064   10.272    0.000    0.658    0.658
#    trusting          0.327    0.081    4.038    0.000    0.327    0.327
#  self_regul ~~                                                         
#    body_listen       0.672    0.066   10.215    0.000    0.672    0.672
#    trusting          0.498    0.085    5.879    0.000    0.498    0.498
#  body_listen ~~                                                        
#    trusting          0.407    0.071    5.713    0.000    0.407    0.407


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.691       0.691
#Tucker-Lewis Index (TLI)                       0.673       0.673
#Robust Comparative Fit Index (CFI)                         0.700
#Robust Tucker-Lewis Index (TLI)                            0.683
#RMSEA                                          0.088       0.082
#Robust RMSEA                                               0.086
#SRMR                                           0.183       0.183

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .49


# Bifactor model with 8 factors
BIFmodel= '
 noticing    =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04
 distracting =~ MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+MAIA2_09+MAIA2_10
 worrying    =~ MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15
 attention   =~ MAIA2_16+MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22
 emotion     =~ MAIA2_23+MAIA2_24+MAIA2_25+MAIA2_26+MAIA2_27
 self_regul  =~ MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31
 body_listen =~ MAIA2_32+MAIA2_33+MAIA2_34
 trusting    =~ MAIA2_35+MAIA2_36+MAIA2_37
 general =~ MAIA2_01+MAIA2_02+MAIA2_03+MAIA2_04+MAIA2_05+MAIA2_06+MAIA2_07+MAIA2_08+
            MAIA2_09+MAIA2_10+MAIA2_11+MAIA2_12+MAIA2_13+MAIA2_14+MAIA2_15+MAIA2_16+
            MAIA2_17+MAIA2_18+MAIA2_19+MAIA2_20+MAIA2_21+MAIA2_22+MAIA2_23+MAIA2_24+
            MAIA2_25+MAIA2_26+MAIA2_27+MAIA2_28+MAIA2_29+MAIA2_30+MAIA2_31+MAIA2_32+
            MAIA2_33+MAIA2_34+MAIA2_35+MAIA2_36+MAIA2_37
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.830       0.830
#Tucker-Lewis Index (TLI)                       0.809       0.808
#Robust Comparative Fit Index (CFI)                         0.838
#Robust Tucker-Lewis Index (TLI)                            0.817
#RMSEA                                          0.067       0.063
#Robust RMSEA                                               0.065
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .504

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  noticing =~                                                           
#    MAIA2_01          0.698    0.169    4.123    0.000    0.698    0.467
#    MAIA2_02          0.813    0.155    5.254    0.000    0.813    0.765
#    MAIA2_03          0.428    0.142    3.023    0.003    0.428    0.309
#    MAIA2_04          0.232    0.126    1.845    0.065    0.232    0.179
#  distracting =~                                                        
#    MAIA2_05          0.698    0.121    5.769    0.000    0.698    0.513
#    MAIA2_06          0.612    0.086    7.123    0.000    0.612    0.528
#    MAIA2_07          0.782    0.122    6.387    0.000    0.782    0.608
#    MAIA2_08          0.959    0.114    8.390    0.000    0.959    0.689
#    MAIA2_09          0.840    0.104    8.086    0.000    0.840    0.735
#    MAIA2_10          0.941    0.101    9.312    0.000    0.941    0.781
#  worrying =~                                                           
#    MAIA2_11          0.922    0.090   10.220    0.000    0.922    0.684
#    MAIA2_12          1.178    0.083   14.253    0.000    1.178    0.797
#    MAIA2_13         -0.867    0.090   -9.674    0.000   -0.867   -0.691
#    MAIA2_14         -0.949    0.083  -11.456    0.000   -0.949   -0.752
#    MAIA2_15          0.918    0.072   12.736    0.000    0.918    0.725
#  attention =~                                                          
#    MAIA2_16          0.448    0.202    2.215    0.027    0.448    0.335
#    MAIA2_17          0.387    0.169    2.286    0.022    0.387    0.311
#    MAIA2_18          0.462    0.109    4.219    0.000    0.462    0.369
#    MAIA2_19          0.753    0.122    6.160    0.000    0.753    0.710
#    MAIA2_20          0.629    0.096    6.527    0.000    0.629    0.559
#    MAIA2_21          0.299    0.134    2.230    0.026    0.299    0.248
#    MAIA2_22          0.275    0.135    2.031    0.042    0.275    0.225
#  emotion =~                                                            
#    MAIA2_23          0.377    0.143    2.636    0.008    0.377    0.267
#    MAIA2_24          0.374    0.140    2.668    0.008    0.374    0.285
#    MAIA2_25          0.550    0.105    5.233    0.000    0.550    0.479
#    MAIA2_26          0.728    0.114    6.397    0.000    0.728    0.508
#    MAIA2_27          0.861    0.099    8.687    0.000    0.861    0.686
#  self_regul =~                                                         
#    MAIA2_28          0.273    0.144    1.901    0.057    0.273    0.191
#    MAIA2_29          0.304    0.125    2.435    0.015    0.304    0.220
#    MAIA2_30          0.930    0.116    8.035    0.000    0.930    0.682
#    MAIA2_31          0.886    0.114    7.747    0.000    0.886    0.636
#  body_listen =~                                                        
#    MAIA2_32          0.500    0.132    3.779    0.000    0.500    0.369
#    MAIA2_33          0.547    0.126    4.348    0.000    0.547    0.407
#    MAIA2_34          0.637    0.147    4.328    0.000    0.637    0.461
#  trusting =~                                                           
#    MAIA2_35          0.873    0.105    8.299    0.000    0.873    0.665
#    MAIA2_36          1.045    0.117    8.921    0.000    1.045    0.814
#    MAIA2_37          0.402    0.098    4.122    0.000    0.402    0.368
#  general =~                                                            
#    MAIA2_01          0.775    0.116    6.705    0.000    0.775    0.519
#    MAIA2_02          0.304    0.087    3.479    0.001    0.304    0.286
#    MAIA2_03          0.602    0.111    5.445    0.000    0.602    0.434
#    MAIA2_04          0.383    0.113    3.378    0.001    0.383    0.295
#    MAIA2_05         -0.329    0.113   -2.908    0.004   -0.329   -0.242
#    MAIA2_06          0.037    0.108    0.346    0.729    0.037    0.032
#    MAIA2_07         -0.146    0.119   -1.232    0.218   -0.146   -0.114
#    MAIA2_08         -0.239    0.128   -1.872    0.061   -0.239   -0.172
#    MAIA2_09          0.172    0.098    1.760    0.078    0.172    0.150
#    MAIA2_10          0.063    0.109    0.582    0.561    0.063    0.052
#    MAIA2_11         -0.028    0.131   -0.215    0.830   -0.028   -0.021
#    MAIA2_12         -0.071    0.137   -0.520    0.603   -0.071   -0.048
#    MAIA2_13          0.146    0.108    1.353    0.176    0.146    0.117
#    MAIA2_14          0.205    0.103    1.991    0.047    0.205    0.162
#    MAIA2_15         -0.104    0.125   -0.830    0.407   -0.104   -0.082
#    MAIA2_16          0.627    0.136    4.591    0.000    0.627    0.468
#    MAIA2_17          0.658    0.124    5.307    0.000    0.658    0.529
#    MAIA2_18          0.354    0.100    3.531    0.000    0.354    0.282
#    MAIA2_19          0.399    0.101    3.959    0.000    0.399    0.377
#    MAIA2_20          0.604    0.086    6.993    0.000    0.604    0.537
#    MAIA2_21          0.620    0.111    5.585    0.000    0.620    0.514
#    MAIA2_22          0.764    0.106    7.221    0.000    0.764    0.624
#    MAIA2_23          0.521    0.117    4.463    0.000    0.521    0.368
#    MAIA2_24          0.630    0.119    5.314    0.000    0.630    0.479
#    MAIA2_25          0.641    0.098    6.546    0.000    0.641    0.559
#    MAIA2_26          0.699    0.118    5.918    0.000    0.699    0.489
#    MAIA2_27          0.636    0.096    6.594    0.000    0.636    0.507
#    MAIA2_28          0.707    0.135    5.235    0.000    0.707    0.494
#    MAIA2_29          1.003    0.096   10.475    0.000    1.003    0.726
#    MAIA2_30          0.697    0.103    6.751    0.000    0.697    0.511
#    MAIA2_31          0.833    0.093    8.947    0.000    0.833    0.598
#    MAIA2_32          0.876    0.101    8.688    0.000    0.876    0.646
#    MAIA2_33          0.878    0.089    9.913    0.000    0.878    0.653
#    MAIA2_34          0.811    0.088    9.198    0.000    0.811    0.587
#    MAIA2_35          0.641    0.093    6.858    0.000    0.641    0.488
#    MAIA2_36          0.645    0.093    6.928    0.000    0.645    0.502
#    MAIA2_37          0.544    0.096    5.675    0.000    0.544    0.497

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3824321      0.8888889      0.9150965      0.7397213 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#noticing    0.5976391 0.05119820 0.40236092 0.6850362 0.3814300 0.6465919 0.8239111
#distracting 0.9521542 0.13946621 0.04784583 0.8174264 0.8127371 0.8324327 0.9148568
#worrying    0.9818347 0.14685864 0.01816533 0.2081921 0.2025383 0.8558463 0.9258581
#attention   0.4336913 0.07002163 0.56630868 0.8214207 0.3337883 0.6644907 0.8261118
#emotion     0.4863525 0.06105186 0.51364746 0.7976610 0.3683747 0.6294766 0.8269586
#self_regul  0.4071619 0.05248710 0.59283807 0.8355143 0.2967951 0.6208712 0.8404075
#body_listen 0.3019361 0.02826389 0.69806391 0.7966937 0.2395178 0.3848094 0.6720200
#trusting    0.6271473 0.06822041 0.37285272 0.8462370 0.5133805 0.7444939 0.9254229
#general     0.3824321 0.38243207 0.38243207 0.9150965 0.7397213 0.9104722 0.9292819








