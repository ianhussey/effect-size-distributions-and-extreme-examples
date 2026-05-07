####### BFI2 Spanish validation study
####### Gallardo-Pujol et al. (2022) https://econtent.hogrefe.com/doi/10.1027/2698-1866/a000020
####### Data Study 1 https://osf.io/zmsr7

BFI.2_Gallardo.Pujol <- read.csv("BFI-2_Gallardo-Pujol.csv")
colnames(BFI.2_Gallardo.Pujol)
mydata <- BFI.2_Gallardo.Pujol[,5:64]
#mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
Names <- colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5


library(psych)
omega(mydata) # alpha = .89, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 10 factors and 8 components
# Eigenvalue1 = 7.64, eigenvalue2 = 4.43

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 11 factors and 8 components
# Eigenvalue 1 = 8.94; eigenvalue 2 = 5.32

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.13, RMSEA=.099, RMSR=.13, TLI=.238

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.15, RMSEA=.121, RMSR=.15, TLI=.201

# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.37, RMSEA=.056, RMSR=.04, TLI=.757
#     MR5   MR3   MR1   MR2   MR4
#MR5  1.00  0.11 -0.14  0.01 -0.13
#MR3  0.11  1.00 -0.16 -0.02 -0.18
#MR1 -0.14 -0.16  1.00  0.14  0.05
#MR2  0.01 -0.02  0.14  1.00  0.11
#MR4 -0.13 -0.18  0.05  0.11  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.44, RMSEA=.076, RMSR=.05, TLI=.68
#     MR2   MR3   MR5   MR1   MR4
#MR2  1.00  0.10  0.01 -0.13 -0.10
#MR3  0.10  1.00 -0.04 -0.16 -0.21
#MR5  0.01 -0.04  1.00  0.15  0.15
#MR1 -0.13 -0.16  0.15  1.00  0.07
#MR4 -0.10 -0.21  0.15  0.07  1.00


#unifactorial model with Lavaan
UNImodel= '
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
library(lavaan)
# ordered model
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.539       0.311
#Tucker-Lewis Index (TLI)                       0.523       0.287
#Robust Comparative Fit Index (CFI)                         0.220
#Robust Tucker-Lewis Index (TLI)                            0.193
#RMSEA                                          0.182       0.124
#Robust RMSEA                                               0.122
#SRMR                                           0.156       0.156

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .154

# MLR model
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.265       0.254
#Tucker-Lewis Index (TLI)                       0.239       0.228
#Robust Comparative Fit Index (CFI)                         0.266
#Robust Tucker-Lewis Index (TLI)                            0.240
#RMSEA                                          0.100       0.095
#Robust RMSEA                                               0.099
#SRMR                                           0.126       0.126

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .105



# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
 factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
 factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
 factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
 factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'

library(lavaan)
# ordered model
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#  Comparative Fit Index (CFI)                    0.846       0.738
#  Tucker-Lewis Index (TLI)                       0.840       0.727
#  Robust Comparative Fit Index (CFI)                         0.599
#  Robust Tucker-Lewis Index (TLI)                            0.583
#  RMSEA                                          0.105       0.077
#  Robust RMSEA                                               0.088
#  SRMR                                           0.102       0.102

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .40

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.125    0.014    8.672    0.000    0.296    0.296
#  factor3          -0.199    0.017  -11.521    0.000   -0.332   -0.332
#  factor4           0.161    0.014   11.858    0.000    0.343    0.343
#  factor5          -0.161    0.016  -10.371    0.000   -0.292   -0.292
#factor2 ~~                                                            
#  factor3          -0.142    0.014   -9.992    0.000   -0.322   -0.322
#  factor4           0.104    0.011    9.480    0.000    0.303    0.303
#  factor5          -0.085    0.013   -6.626    0.000   -0.209   -0.209
#  factor3 ~~                                                            
#  factor4          -0.078    0.014   -5.365    0.000   -0.159   -0.159
#  factor5           0.032    0.016    1.945    0.052    0.055    0.055
#factor4 ~~                                                            
#  factor5           0.009    0.013    0.651    0.515    0.019    0.019

# MLR model
CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.656       0.656
#Tucker-Lewis Index (TLI)                       0.642       0.642
#Robust Comparative Fit Index (CFI)                         0.660
#Robust Tucker-Lewis Index (TLI)                            0.646
#RMSEA                                          0.068       0.065
#Robust RMSEA                                               0.068
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .306

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.270    0.045    6.017    0.000    0.270    0.270
#    factor3          -0.227    0.040   -5.692    0.000   -0.227   -0.227
#    factor4           0.314    0.042    7.566    0.000    0.314    0.314
#    factor5          -0.173    0.039   -4.377    0.000   -0.173   -0.173
#  factor2 ~~                                                            
#    factor3          -0.292    0.035   -8.305    0.000   -0.292   -0.292
#    factor4           0.269    0.045    6.001    0.000    0.269    0.269
#    factor5          -0.200    0.035   -5.628    0.000   -0.200   -0.200
#  factor3 ~~                                                            
#    factor4          -0.149    0.037   -4.078    0.000   -0.149   -0.149
#    factor5          -0.009    0.035   -0.257    0.797   -0.009   -0.009
#  factor4 ~~                                                            
#    factor5           0.039    0.036    1.059    0.289    0.039    0.039

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# model with orthogonal factors
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.644       0.644
#Tucker-Lewis Index (TLI)                       0.631       0.631
#Robust Comparative Fit Index (CFI)                         0.648
#Robust Tucker-Lewis Index (TLI)                            0.635
#RMSEA                                          0.069       0.066
#Robust RMSEA                                               0.069
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .301


#bifactor model better MLR
BIFmodel= '
factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.729       0.728
#Tucker-Lewis Index (TLI)                       0.709       0.708
#Robust Comparative Fit Index (CFI)                         0.733
#Robust Tucker-Lewis Index (TLI)                            0.713
#RMSEA                                          0.062       0.058
#Robust RMSEA                                               0.061
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .329

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#   0.28176346    0.81355932    0.03049742    0.01516564 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS       Omega       OmegaH         H        FD
#factor1 0.3390446 0.07081309 0.6609554 0.022192133 5.826307e-06 0.7254342 0.8734864
#factor2 0.8420636 0.12150665 0.1579364 0.034276835 2.453118e-02 0.7724687 0.8763181
#factor3 0.8028471 0.17163270 0.1971529 0.053209804 4.490675e-02 0.8794420 0.9359599
#factor4 0.7756307 0.18573267 0.2243693 0.008835037 3.627250e-03 0.8667325 0.9282658
#factor5 0.8706044 0.16855143 0.1293956 0.002988062 8.739637e-04 0.8738259 0.9340896
#global  0.2817635 0.28176346 0.2817635 0.030497421 1.516564e-02 0.8895783 0.9246508


#Reverse code items to see whether the combination of negative and positive values drives the low
#omegaH estimates
mydata2 <- mydata
#Items that are worded negatively
reverse_cols = c("BFI11","BFI16","BFI26","BFI31","BFI36","BFI51","BFI12","BFI17","BFI22","BFI37","BFI42","BFI47",
                 "BFI13","BFI18","BFI33","BFI38","BFI43","BFI53","BFI14","BFI19","BFI34","BFI39","BFI54","BFI59",
                 "BFI10","BFI15","BFI20","BFI35","BFI40","BFI60")
mydata2[ , reverse_cols] = 6 - mydata2[ , reverse_cols]
#Factor that correlate negatively with the others (3 and 5)
reverse_cols = c("BFI3","BFI8","BFI13","BFI18","BFI23","BFI28","BFI33","BFI38","BFI43","BFI48","BFI53","BFI58",
                 "BFI5","BFI10","BFI15","BFI20","BFI25","BFI30","BFI35","BFI40","BFI45","BFI50","BFI55","BFI60")
mydata2[ , reverse_cols] = 6 - mydata2[ , reverse_cols]
CFA_model3b <- cfa(BIFmodel, data = mydata2,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3b, standardized = TRUE, fit.measures = TRUE)
semPaths(CFA_model3b,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")
bifactorIndices(CFA_model3b)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2817635     0.8135593     0.9146026     0.5633871 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3390446 0.07081309 0.6609554 0.8439553 0.1776706 0.7254342 0.8734864
#factor2 0.8420636 0.12150665 0.1579364 0.7964039 0.6814309 0.7724687 0.8763181
#factor3 0.8028471 0.17163270 0.1971529 0.8611444 0.6881740 0.8794420 0.9359599
#factor4 0.7756307 0.18573267 0.2243693 0.8874052 0.7429756 0.8667325 0.9282658
#factor5 0.8706044 0.16855143 0.1293956 0.8432228 0.7609520 0.8738259 0.9340896
#global  0.2817635 0.28176346 0.2817635 0.9146026 0.5633871 0.8895783 0.9246508

# So, no influence on ECV but big influence on omegaH







####### Big Five Personality Trait Short Questionnaire
####### Mezquita et al. (2019) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6917272/
####### Data https://osf.io/r65e7
####### Language = 1 (English)

library(haven)
Big_Five_Mezquita <- read_sav("Big_Five_Mezquita.sav")
colnames(Big_Five_Mezquita)
mydata <- subset(Big_Five_Mezquita,Language2==1)
mydata <- mydata[,5:54]
#mydata[mydata<=0] <- NA
colnames(mydata)
mydata <- na.omit(mydata)
colnames(mydata) <- paste0('BFI', 1:(ncol(mydata)))
min(mydata)
max(mydata) # response alternatives = 5


library(psych)
omega(mydata) # alpha = .84, omega T = .88

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 6 components
# Eigenvalue1 = 8.86, eigenvalue2 = 4.92

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 10.42; eigenvalue 2 = 5.56

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.113, RMSR=.14, TLI=.314

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.21, RMSEA=.134, RMSR=.16, TLI=.288

# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.44, RMSEA=.06, RMSR=.04, TLI=.802
#      MR1  MR5  MR3  MR2   MR4
#MR1  1.00 0.37 0.06 0.04 -0.13
#MR5  0.37 1.00 0.19 0.18  0.13
#MR3  0.06 0.19 1.00 0.27  0.02
#MR2  0.04 0.18 0.27 1.00  0.17
#MR4 -0.13 0.13 0.02 0.17  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.50, RMSEA=.081, RMSR=.05, TLI=.74
#      MR1  MR5  MR3   MR4  MR2
#MR1  1.00 0.38 0.08 -0.11 0.03
#MR5  0.38 1.00 0.21  0.15 0.17
#MR3  0.08 0.21 1.00  0.02 0.27
#MR4 -0.11 0.15 0.02  1.00 0.16
#MR2  0.03 0.17 0.27  0.16 1.00

colnames(mydata) <- Names[1:50]

#unifactorial model with Lavaan
UNImodel= '
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + 
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + 
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + 
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + 
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
'
library(lavaan)
# ordered model
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.674       0.443
#Tucker-Lewis Index (TLI)                       0.660       0.419
#Robust Comparative Fit Index (CFI)                         0.286
#Robust Tucker-Lewis Index (TLI)                            0.256
#RMSEA                                          0.202       0.126
#Robust RMSEA                                               0.139
#SRMR                                           0.167       0.167

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .238

# MLR model
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.343       0.329
#Tucker-Lewis Index (TLI)                       0.315       0.300
#Robust Comparative Fit Index (CFI)                         0.344
#Robust Tucker-Lewis Index (TLI)                            0.316
#RMSEA                                          0.114       0.105
#Robust RMSEA                                               0.113
#SRMR                                           0.138       0.138

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .149


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46
 factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47
 factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48
 factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49
 factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
'

library(lavaan)
# ordered model
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.789       0.659
#Tucker-Lewis Index (TLI)                       0.778       0.641
#Robust Comparative Fit Index (CFI)                         0.515
#Robust Tucker-Lewis Index (TLI)                            0.490
#RMSEA                                          0.163       0.099
#Robust RMSEA                                               0.115
#SRMR                                           0.140       0.140

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .444

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.468    0.022   21.637    0.000    0.468    0.468
#    factor3           0.551    0.023   24.156    0.000    0.551    0.551
#    factor4           0.460    0.024   19.079    0.000    0.460    0.460
#    factor5           0.136    0.028    4.918    0.000    0.136    0.136
#  factor2 ~~                                                            
#    factor3           0.435    0.023   19.285    0.000    0.435    0.435
#    factor4           0.455    0.023   19.557    0.000    0.455    0.455
#    factor5           0.429    0.024   17.988    0.000    0.429    0.429
#  factor3 ~~                                                            
#    factor4           0.618    0.021   29.380    0.000    0.618    0.618
#    factor5           0.310    0.026   11.713    0.000    0.310    0.310
#  factor4 ~~                                                            
#    factor5           0.409    0.024   16.965    0.000    0.409    0.409

# MLR model
CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.572       0.563
#Tucker-Lewis Index (TLI)                       0.550       0.541
#Robust Comparative Fit Index (CFI)                         0.576
#Robust Tucker-Lewis Index (TLI)                            0.554
#RMSEA                                          0.092       0.085
#Robust RMSEA                                               0.091
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .326

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.293    0.058    5.040    0.000    0.293    0.293
#    factor3           0.550    0.050   10.950    0.000    0.550    0.550
#    factor4           0.436    0.055    7.870    0.000    0.436    0.436
#    factor5           0.027    0.049    0.555    0.579    0.027    0.027
#  factor2 ~~                                                            
#    factor3           0.291    0.055    5.265    0.000    0.291    0.291
#    factor4           0.342    0.048    7.147    0.000    0.342    0.342
#    factor5           0.405    0.037   10.975    0.000    0.405    0.405
#  factor3 ~~                                                            
#    factor4           0.633    0.037   16.910    0.000    0.633    0.633
#    factor5           0.161    0.061    2.629    0.009    0.161    0.161
#  factor4 ~~                                                            
#    factor5           0.296    0.055    5.420    0.000    0.296    0.296


library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# MLR model with independent factors
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.535       0.527
#Tucker-Lewis Index (TLI)                       0.516       0.507
#Robust Comparative Fit Index (CFI)                         0.539
#Robust Tucker-Lewis Index (TLI)                            0.519
#RMSEA                                          0.095       0.088
#Robust RMSEA                                               0.095
#SRMR                                           0.171       0.171

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307


#bifactor MLR model
BIFmodel= '
factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46
factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47
factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48
factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49
factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.770       0.773
#Tucker-Lewis Index (TLI)                       0.749       0.753
#Robust Comparative Fit Index (CFI)                         0.776
#Robust Tucker-Lewis Index (TLI)                            0.756
#RMSEA                                          0.069       0.062
#Robust RMSEA                                               0.067
#SRMR                                           0.103       0.103

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .419

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.3377459     0.8163265     0.9036237     0.4702005 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4911571 0.08609684 0.5088429 0.8263469 0.4201172 0.7274141 0.8715803
#factor2 0.6625660 0.16532776 0.3374340 0.8982555 0.6250613 0.9208698 0.9605637
#factor3 0.6302846 0.11472706 0.3697154 0.8118904 0.6933319 0.7938121 0.8996473
#factor4 0.6444702 0.12470592 0.3555298 0.8369780 0.7284042 0.7903342 0.9001387
#factor5 0.8584676 0.17139648 0.1415324 0.8642260 0.8579452 0.8581695 0.9318518
#global  0.3377459 0.33774594 0.3377459 0.9036237 0.4702005 0.9092018 0.9542798






####### Big Five Personality Trait Short Questionnaire
####### Mezquita et al. (2019) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6917272/
####### Data https://osf.io/r65e7
####### Language = 2 (Spanish)

library(haven)
Big_Five_Mezquita <- read_sav("Big_Five_Mezquita.sav")
colnames(Big_Five_Mezquita)
mydata <- subset(Big_Five_Mezquita,Language2==2)
mydata <- mydata[,5:54]
#mydata[mydata<=0] <- NA
colnames(mydata)
colnames(mydata) <- paste0('BFI', 1:(ncol(mydata)))
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .89, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 6 components
# Eigenvalue1 = 7.9, eigenvalue2 = 3.86

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 7 components
# Eigenvalue 1 = 9.26, eigenvalue 2 = 4.46

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.108, RMSR=.13, TLI=.285

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.129, RMSR=.14, TLI=.257

# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.4, RMSEA=.063, RMSR=.04, TLI=.757
#     MR1  MR2  MR3  MR5  MR4
#MR1 1.00 0.05 0.19 0.28 0.08
#MR2 0.05 1.00 0.16 0.14 0.17
#MR3 0.19 0.16 1.00 0.10 0.11
#MR5 0.28 0.14 0.10 1.00 0.21
#MR4 0.08 0.17 0.11 0.21 1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.46, RMSEA=.083, RMSR=.05, TLI=.695
#     MR1  MR2  MR3  MR5  MR4
#MR1 1.00 0.04 0.21 0.29 0.09
#MR2 0.04 1.00 0.16 0.14 0.17
#MR3 0.21 0.16 1.00 0.11 0.11
#MR5 0.29 0.14 0.11 1.00 0.22
#MR4 0.09 0.17 0.11 0.22 1.00

colnames(mydata) <- Names[1:50]

#unifactorial model with Lavaan
UNImodel= '
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + 
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + 
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + 
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + 
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.678       0.489
#Tucker-Lewis Index (TLI)                       0.664       0.467
#Robust Comparative Fit Index (CFI)                         0.248
#Robust Tucker-Lewis Index (TLI)                            0.216
#RMSEA                                          0.179       0.122
#Robust RMSEA                                               0.134
#SRMR                                           0.155       0.155

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .224

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.315       0.302
#Tucker-Lewis Index (TLI)                       0.286       0.272
#Robust Comparative Fit Index (CFI)                         0.316
#Robust Tucker-Lewis Index (TLI)                            0.287
#RMSEA                                          0.109       0.103
#Robust RMSEA                                               0.109
#SRMR                                           0.124       0.124

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .156


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46
 factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47
 factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48
 factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49
 factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
'

library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.854       0.754
#Tucker-Lewis Index (TLI)                       0.847       0.742
#Robust Comparative Fit Index (CFI)                         0.586
#Robust Tucker-Lewis Index (TLI)                            0.565
#RMSEA                                          0.121       0.085
#Robust RMSEA                                               0.100
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .434

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.421    0.024   17.852    0.000    0.421    0.421
#    factor3           0.370    0.028   13.108    0.000    0.370    0.370
#    factor4           0.303    0.027   11.149    0.000    0.303    0.303
#    factor5           0.181    0.029    6.294    0.000    0.181    0.181
#  factor2 ~~                                                            
#    factor3           0.368    0.028   13.142    0.000    0.368    0.368
#    factor4           0.303    0.029   10.469    0.000    0.303    0.303
#    factor5           0.347    0.027   13.024    0.000    0.347    0.347
#  factor3 ~~                                                            
#    factor4           0.436    0.027   16.084    0.000    0.436    0.436
#    factor5           0.402    0.028   14.235    0.000    0.402    0.402
#  factor4 ~~                                                            
#    factor5           0.287    0.030    9.542    0.000    0.287    0.287

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.638       0.635
#Tucker-Lewis Index (TLI)                       0.619       0.616
#Robust Comparative Fit Index (CFI)                         0.641
#Robust Tucker-Lewis Index (TLI)                            0.623
#RMSEA                                          0.080       0.074
#Robust RMSEA                                               0.079
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.294    0.037    8.001    0.000    0.294    0.294
#    factor3           0.327    0.052    6.338    0.000    0.327    0.327
#    factor4           0.228    0.057    3.974    0.000    0.228    0.228
#    factor5           0.152    0.036    4.172    0.000    0.152    0.152
#  factor2 ~~                                                            
#    factor3           0.276    0.043    6.391    0.000    0.276    0.276
#    factor4           0.206    0.043    4.821    0.000    0.206    0.206
#    factor5           0.270    0.038    7.115    0.000    0.270    0.270
#  factor3 ~~                                                            
#    factor4           0.392    0.049    8.010    0.000    0.392    0.392
#    factor5           0.349    0.048    7.224    0.000    0.349    0.349
#  factor4 ~~                                                            
#    factor5           0.245    0.045    5.494    0.000    0.245    0.245

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.617       0.616
#Tucker-Lewis Index (TLI)                       0.601       0.599
#Robust Comparative Fit Index (CFI)                         0.621
#Robust Tucker-Lewis Index (TLI)                            0.605
#RMSEA                                          0.082       0.076
#Robust RMSEA                                               0.081
#SRMR                                           0.139       0.139

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .316

#bifactor model
BIFmodel= '
factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46
factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47
factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48
factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49
factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.759       0.759
#Tucker-Lewis Index (TLI)                       0.738       0.737
#Robust Comparative Fit Index (CFI)                         0.764
#Robust Tucker-Lewis Index (TLI)                            0.743
#RMSEA                                          0.066       0.062
#Robust RMSEA                                               0.065
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .373


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2900418     0.8163265     0.9020008     0.5105441 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5320557 0.1042470 0.4679443 0.8307188 0.3996471 0.8173991 0.9240313
#factor2 0.7260931 0.1705749 0.2739069 0.8686093 0.6353183 0.9059466 0.9511934
#factor3 0.6482503 0.1035983 0.3517497 0.7849134 0.6021918 0.7275107 0.8543859
#factor4 0.7181644 0.1347069 0.2818356 0.8184891 0.6946743 0.8074171 0.9015063
#factor5 0.8875755 0.1968311 0.1124245 0.8716958 0.8473248 0.8690459 0.9339460
#global  0.2900418 0.2900418 0.2900418 0.9020008 0.5105441 0.8789198 0.9269792



##################################################################################
### BFI2 Vermeiren study https://link.springer.com/article/10.3758/s13428-022-01856-x
### Summed observations of Study1 (N = 195), Study3 (N = 196)
### Data https://osf.io/ef3s4/

library(readxl)
BFI2_Vermeiren <- read_excel("BFI2_Vermeiren.xlsx")
BFI2_Vermeiren <- subset(BFI2_Vermeiren, Study=="S1")
colnames(BFI2_Vermeiren)
mydata <- as.data.frame(BFI2_Vermeiren[,2:61])
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
omega(mydata) # alpha = .91, omega = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 6 components
# Eigenvalue1 = 9.97, eigenvalue2 = 4.8

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 6 components 
# Eigenvalue 1 = 11.53; eigenvalue 2 = 5.71


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.17, RMSEA=.099, RMSR=.14, TLI=.283

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.206, RMSR=.16, TLI=.081


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.43, RMSEA=.057, RMSR=.05, TLI=.76
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.16 -0.22 -0.03  0.09
#MR2 -0.16  1.00  0.12  0.13 -0.20
#MR4 -0.22  0.12  1.00  0.20  0.02
#MR3 -0.03  0.13  0.20  1.00 -0.11
#MR5  0.09 -0.20  0.02 -0.11  1.00


fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.49, RMSEA=.193, RMSR=.06, TLI=.177
#     MR2   MR3   MR5   MR1   MR4
#MR1  1.00 -0.16  0.18  0.28 -0.01
#MR2 -0.16  1.00 -0.15  0.03  0.16
#MR5  0.18 -0.15  1.00  0.11 -0.17
#MR4  0.28  0.03  0.11  1.00 -0.07
#MR3 -0.01  0.16 -0.17 -0.07  1.00


# Single factor model
UNImodel= '
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
library(lavaan)
CFA_model4 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.672       0.485
#Tucker-Lewis Index (TLI)                       0.661       0.467
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.187       0.107
#Robust RMSEA                                                  NA
#SRMR                                           0.167       0.167

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .207

CFA_model4 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.300       0.294
#Tucker-Lewis Index (TLI)                       0.276       0.270
#Robust Comparative Fit Index (CFI)                         0.302
#Robust Tucker-Lewis Index (TLI)                            0.277
#RMSEA                                          0.109       0.105
#Robust RMSEA                                               0.108
#SRMR                                           0.141       0.141

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .123

semPaths(CFA_model3,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# model with theory-based structure
EGAmodel= '
 factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
 factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
 factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
 factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
 factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.891       0.807
#Tucker-Lewis Index (TLI)                       0.887       0.799
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.108       0.066
#Robust RMSEA                                                  NA
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .447

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.193    0.103    1.864    0.062    0.193    0.193
#    factor3          -0.214    0.094   -2.265    0.023   -0.214   -0.214
#    factor4           0.420    0.070    5.973    0.000    0.420    0.420
#    factor5          -0.334    0.083   -4.013    0.000   -0.334   -0.334
#  factor2 ~~                                                            
#    factor3          -0.518    0.088   -5.889    0.000   -0.518   -0.518
#    factor4           0.247    0.095    2.593    0.010    0.247    0.247
#    factor5          -0.228    0.085   -2.686    0.007   -0.228   -0.228
#  factor3 ~~                                                            
#    factor4          -0.289    0.080   -3.624    0.000   -0.289   -0.289
#    factor5           0.081    0.086    0.942    0.346    0.081    0.081
#  factor4 ~~                                                            
#    factor5          -0.044    0.092   -0.477    0.633   -0.044   -0.044

CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.647       0.649
#Tucker-Lewis Index (TLI)                       0.632       0.634
#Robust Comparative Fit Index (CFI)                         0.654
#Robust Tucker-Lewis Index (TLI)                            0.639
#RMSEA                                          0.077       0.074
#Robust RMSEA                                               0.076
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .363

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.193    0.103    1.864    0.062    0.193    0.193
#    factor3          -0.214    0.094   -2.265    0.023   -0.214   -0.214
#    factor4           0.420    0.070    5.973    0.000    0.420    0.420
#    factor5          -0.334    0.083   -4.013    0.000   -0.334   -0.334
#  factor2 ~~                                                            
#    factor3          -0.518    0.088   -5.889    0.000   -0.518   -0.518
#    factor4           0.247    0.095    2.593    0.010    0.247    0.247
#    factor5          -0.228    0.085   -2.686    0.007   -0.228   -0.228
#  factor3 ~~                                                            
#    factor4          -0.289    0.080   -3.624    0.000   -0.289   -0.289
#    factor5           0.081    0.086    0.942    0.346    0.081    0.081
#  factor4 ~~                                                            
#    factor5          -0.044    0.092   -0.477    0.633   -0.044   -0.044

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.628       0.631
#Tucker-Lewis Index (TLI)                       0.615       0.618
#Robust Comparative Fit Index (CFI)                         0.635
#Robust Tucker-Lewis Index (TLI)                            0.622
#RMSEA                                          0.079       0.076
#Robust RMSEA                                               0.078
#SRMR                                           0.142       0.142

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .362


#bifactor model better
BIFmodel= '
factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.705       0.695
#Tucker-Lewis Index (TLI)                       0.684       0.673
#Robust Comparative Fit Index (CFI)                         0.706
#Robust Tucker-Lewis Index (TLI)                            0.685
#RMSEA                                          0.072       0.070
#Robust RMSEA                                               0.071
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .394

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    BFI1              0.770    0.130    5.919    0.000    0.770    0.587
#    BFI6              0.496    0.116    4.281    0.000    0.496    0.417
#    BFI11            -0.207    0.167   -1.238    0.216   -0.207   -0.187
#    BFI16            -0.989    0.089  -11.128    0.000   -0.989   -0.813
#    BFI21             0.675    0.118    5.731    0.000    0.675    0.532
#    BFI26            -0.389    0.147   -2.645    0.008   -0.389   -0.306
#    BFI31            -0.849    0.077  -11.020    0.000   -0.849   -0.778
#    BFI36            -0.442    0.114   -3.865    0.000   -0.442   -0.400
#    BFI41             0.483    0.156    3.090    0.002    0.483    0.428
#    BFI46             0.755    0.114    6.605    0.000    0.755    0.616
#    BFI51            -0.570    0.114   -4.988    0.000   -0.570   -0.483
#    BFI56             0.314    0.150    2.091    0.037    0.314    0.321
#  factor2 =~                                                            
#    BFI2              0.394    0.092    4.280    0.000    0.394    0.459
#    BFI7              0.333    0.080    4.155    0.000    0.333    0.452
#    BFI12            -0.540    0.101   -5.339    0.000   -0.540   -0.474
#    BFI17            -0.457    0.116   -3.940    0.000   -0.457   -0.333
#    BFI22            -0.623    0.093   -6.728    0.000   -0.623   -0.547
#    BFI27             0.482    0.085    5.689    0.000    0.482    0.501
#    BFI32             0.215    0.103    2.093    0.036    0.215    0.223
#    BFI37            -0.816    0.090   -9.020    0.000   -0.816   -0.693
#    BFI42            -0.485    0.101   -4.780    0.000   -0.485   -0.432
#    BFI47            -0.749    0.087   -8.575    0.000   -0.749   -0.653
#    BFI52             0.325    0.078    4.148    0.000    0.325    0.471
#    BFI57             0.413    0.116    3.551    0.000    0.413    0.364
#  factor3 =~                                                            
#    BFI3              1.137    0.083   13.660    0.000    1.137    0.826
#    BFI8              0.515    0.138    3.743    0.000    0.515    0.442
#    BFI13            -0.379    0.087   -4.374    0.000   -0.379   -0.394
#    BFI18            -0.727    0.106   -6.829    0.000   -0.727   -0.605
#    BFI23             0.489    0.132    3.706    0.000    0.489    0.391
#    BFI28             0.765    0.114    6.702    0.000    0.765    0.636
#    BFI33            -0.802    0.078  -10.303    0.000   -0.802   -0.672
#    BFI38            -0.487    0.101   -4.839    0.000   -0.487   -0.489
#    BFI43            -0.367    0.101   -3.619    0.000   -0.367   -0.413
#    BFI48             0.651    0.083    7.846    0.000    0.651    0.580
#    BFI53            -0.507    0.100   -5.055    0.000   -0.507   -0.510
#    BFI58             0.638    0.118    5.391    0.000    0.638    0.534
#  factor4 =~                                                            
#    BFI4              0.800    0.100    7.973    0.000    0.800    0.632
#    BFI9              0.425    0.124    3.429    0.001    0.425    0.393
#    BFI14            -0.858    0.083  -10.340    0.000   -0.858   -0.671
#    BFI19            -0.734    0.073  -10.044    0.000   -0.734   -0.679
#    BFI24             0.577    0.116    4.966    0.000    0.577    0.462
#    BFI29             0.902    0.078   11.519    0.000    0.902    0.744
#    BFI34            -0.941    0.087  -10.828    0.000   -0.941   -0.724
#    BFI39            -1.007    0.081  -12.405    0.000   -1.007   -0.777
#    BFI44             0.668    0.074    8.983    0.000    0.668    0.615
#    BFI49             0.751    0.092    8.133    0.000    0.751    0.637
#    BFI54            -1.047    0.094  -11.135    0.000   -1.047   -0.779
#    BFI59            -0.936    0.063  -14.892    0.000   -0.936   -0.782
#  factor5 =~                                                            
#    BFI5              0.753    0.102    7.366    0.000    0.753    0.562
#    BFI10            -0.453    0.071   -6.341    0.000   -0.453   -0.508
#    BFI15            -0.326    0.108   -3.005    0.003   -0.326   -0.326
#    BFI20            -0.805    0.094   -8.546    0.000   -0.805   -0.657
#    BFI25             0.418    0.114    3.657    0.000    0.418    0.370
#    BFI30             0.735    0.100    7.310    0.000    0.735    0.647
#    BFI35            -0.748    0.104   -7.195    0.000   -0.748   -0.631
#    BFI40            -0.415    0.079   -5.248    0.000   -0.415   -0.438
#    BFI45             0.564    0.118    4.767    0.000    0.564    0.511
#    BFI50             0.830    0.090    9.254    0.000    0.830    0.614
#    BFI55             0.792    0.099    8.023    0.000    0.792    0.678
#    BFI60            -0.446    0.104   -4.297    0.000   -0.446   -0.438
#  global =~                                                             
#    BFI1              0.616    0.183    3.359    0.001    0.616    0.469
#    BFI6              0.524    0.150    3.482    0.000    0.524    0.440
#    BFI11            -0.459    0.234   -1.959    0.050   -0.459   -0.415
#    BFI16            -0.020    0.180   -0.109    0.913   -0.020   -0.016
#    BFI21             0.604    0.167    3.605    0.000    0.604    0.476
#    BFI26            -0.374    0.230   -1.624    0.104   -0.374   -0.294
#    BFI31            -0.040    0.202   -0.197    0.844   -0.040   -0.036
#    BFI36            -0.277    0.175   -1.586    0.113   -0.277   -0.251
#    BFI41             0.564    0.201    2.803    0.005    0.564    0.500
#    BFI46             0.370    0.170    2.186    0.029    0.370    0.302
#    BFI51            -0.422    0.137   -3.089    0.002   -0.422   -0.357
#    BFI56             0.639    0.143    4.451    0.000    0.639    0.654
#    BFI2              0.292    0.122    2.395    0.017    0.292    0.339
#    BFI7              0.273    0.103    2.642    0.008    0.273    0.372
#    BFI12            -0.072    0.187   -0.387    0.699   -0.072   -0.063
#    BFI17             0.047    0.194    0.242    0.809    0.047    0.034
#    BFI22             0.148    0.190    0.782    0.434    0.148    0.130
#    BFI27             0.229    0.097    2.371    0.018    0.229    0.238
#    BFI32             0.404    0.084    4.808    0.000    0.404    0.419
#    BFI37             0.042    0.282    0.151    0.880    0.042    0.036
#    BFI42            -0.136    0.169   -0.804    0.421   -0.136   -0.121
#    BFI47            -0.130    0.198   -0.659    0.510   -0.130   -0.114
#    BFI52             0.253    0.101    2.497    0.013    0.253    0.367
#    BFI57             0.380    0.141    2.698    0.007    0.380    0.336
#    BFI3             -0.272    0.181   -1.502    0.133   -0.272   -0.197
#    BFI8             -0.330    0.207   -1.592    0.111   -0.330   -0.283
#    BFI13             0.395    0.106    3.730    0.000    0.395    0.411
#    BFI18             0.371    0.136    2.727    0.006    0.371    0.309
#    BFI23            -0.264    0.210   -1.253    0.210   -0.264   -0.211
#    BFI28            -0.040    0.236   -0.171    0.864   -0.040   -0.033
#    BFI33             0.415    0.114    3.647    0.000    0.415    0.348
#    BFI38             0.582    0.112    5.193    0.000    0.582    0.584
#    BFI43             0.378    0.109    3.472    0.001    0.378    0.425
#    BFI48            -0.168    0.162   -1.035    0.301   -0.168   -0.150
#    BFI53             0.490    0.121    4.045    0.000    0.490    0.493
#    BFI58             0.012    0.220    0.057    0.955    0.012    0.010
#    BFI4              0.444    0.165    2.690    0.007    0.444    0.351
#    BFI9              0.495    0.134    3.692    0.000    0.495    0.458
#    BFI14            -0.060    0.230   -0.262    0.794   -0.060   -0.047
#    BFI19             0.064    0.214    0.300    0.764    0.064    0.059
#    BFI24             0.511    0.140    3.640    0.000    0.511    0.410
#    BFI29             0.329    0.188    1.752    0.080    0.329    0.272
#    BFI34            -0.193    0.198   -0.971    0.332   -0.193   -0.148
#    BFI39            -0.290    0.262   -1.106    0.269   -0.290   -0.224
#    BFI44             0.403    0.119    3.378    0.001    0.403    0.371
#    BFI49             0.341    0.181    1.888    0.059    0.341    0.290
#    BFI54            -0.375    0.271   -1.387    0.165   -0.375   -0.279
#    BFI59            -0.080    0.212   -0.378    0.706   -0.080   -0.067
#    BFI5              0.097    0.139    0.697    0.486    0.097    0.072
#    BFI10             0.416    0.086    4.857    0.000    0.416    0.466
#    BFI15             0.467    0.092    5.081    0.000    0.467    0.467
#    BFI20             0.193    0.127    1.527    0.127    0.193    0.158
#    BFI25            -0.331    0.136   -2.433    0.015   -0.331   -0.293
#    BFI30            -0.111    0.149   -0.750    0.453   -0.111   -0.098
#    BFI35             0.239    0.113    2.115    0.034    0.239    0.202
#    BFI40             0.294    0.075    3.908    0.000    0.294    0.310
#    BFI45            -0.210    0.159   -1.318    0.187   -0.210   -0.190
#    BFI50            -0.034    0.159   -0.216    0.829   -0.034   -0.025
#    BFI55            -0.175    0.189   -0.929    0.353   -0.175   -0.150
#    BFI60             0.371    0.118    3.156    0.002    0.371    0.364

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2379418     0.8135593     0.6417197     0.6260185 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega       OmegaH         H        FD
#factor1 0.6363225 0.1329343 0.3636775 0.2397243 0.0004755233 0.8610867 0.9227172
#factor2 0.7811013 0.1142663 0.2188987 0.3398013 0.0344728828 0.8034231 0.8967681
#factor3 0.7345918 0.1510643 0.2654082 0.3020706 0.0105807779 0.8685894 0.9281452
#factor4 0.8486938 0.2189751 0.1513062 0.3415251 0.0995525943 0.9171121 0.9557325
#factor5 0.7994695 0.1448182 0.2005305 0.1914350 0.0158029459 0.8464345 0.9179067
#global  0.2379418 0.2379418 0.2379418 0.6417197 0.6260184525 0.8791829 0.9212904


#Reverse code items to see whether the combination of negative and positive values drives the low
#omegaH estimates
mydata2 <- mydata
reverse_cols = c("BFI11","BFI16","BFI26","BFI31","BFI36","BFI51","BFI12","BFI17","BFI22","BFI37","BFI42","BFI47",
                 "BFI13","BFI18","BFI33","BFI38","BFI43","BFI53","BFI14","BFI19","BFI34","BFI39","BFI54","BFI59",
                 "BFI10","BFI15","BFI20","BFI35","BFI40","BFI60")
mydata2[ , reverse_cols] = 6 - mydata2[ , reverse_cols]
CFA_model3 <- cfa(BIFmodel, data = mydata2,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#   0.23794178    0.81355932    0.86151245    0.03876786 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.6363225 0.1329343 0.3636775 0.8834168 0.58305295 0.8610867 0.9227172
#factor2 0.7811013 0.1142663 0.2188987 0.8108545 0.70517617 0.8034231 0.8967681
#factor3 0.7345918 0.1510643 0.2654082 0.8856263 0.69212999 0.8685894 0.9281452
#factor4 0.8486938 0.2189751 0.1513062 0.9254058 0.81830117 0.9171121 0.9557325
#factor5 0.7994695 0.1448182 0.2005305 0.8631676 0.73604578 0.8464345 0.9179067
#global  0.2379418 0.2379418 0.2379418 0.8615125 0.03876786 0.8791829 0.9212904


semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")




##################################################################################
### BFI2 Vermeiren study https://link.springer.com/article/10.3758/s13428-022-01856-x
### Summed observations of Study3 (N = 196)
### Data https://osf.io/ef3s4/

library(readxl)
BFI2_Vermeiren <- read_excel("BFI2_Vermeiren.xlsx")
BFI2_Vermeiren <- subset(BFI2_Vermeiren, Study=="S3")
colnames(BFI2_Vermeiren)
mydata <- as.data.frame(BFI2_Vermeiren[,2:61])
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
omega(mydata) # alpha = .89, omega = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 6 components
# Eigenvalue1 = 8.45, eigenvalue2 = 4.44

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 6 components 
# Eigenvalue 1 = 9.65; eigenvalue 2 = 5.29


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.088, RMSR=.13, TLI=.287

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.16, RMSEA=.195, RMSR=.15, TLI=.067


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.37, RMSEA=.052, RMSR=.06, TLI=.738
#      MR1   MR2   MR3   MR4   MR5
#MR1  1.00 -0.09 -0.18 -0.17 -0.07
#MR2 -0.09  1.00  0.14  0.18  0.06
#MR3 -0.18  0.14  1.00  0.01  0.13
#MR4 -0.17  0.18  0.01  1.00  0.10
#MR5 -0.07  0.06  0.13  0.10  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.43, RMSEA=.186, RMSR=.07, TLI=.134
#      MR1   MR2   MR5   MR3   MR4
#MR1  1.00 -0.04 -0.18 -0.13 -0.13
#MR2 -0.04  1.00  0.16  0.10  0.09
#MR5 -0.18  0.16  1.00  0.11 -0.01
#MR3 -0.13  0.10  0.11  1.00  0.13
#MR4 -0.13  0.09 -0.01  0.13  1.00


# Single factor model
UNImodel= '
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
library(lavaan)
CFA_model4 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.642       0.505
#Tucker-Lewis Index (TLI)                       0.630       0.488
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.153       0.085
#Robust RMSEA                                                  NA
#SRMR                                           0.151       0.151

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .160

CFA_model4 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.298       0.289
#Tucker-Lewis Index (TLI)                       0.274       0.264
#Robust Comparative Fit Index (CFI)                         0.297
#Robust Tucker-Lewis Index (TLI)                            0.272
#RMSEA                                          0.097       0.095
#Robust RMSEA                                               0.096
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .128

semPaths(CFA_model3,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# model with theory-based structure
EGAmodel= '
 factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
 factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
 factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
 factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
 factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.865       0.787
#Tucker-Lewis Index (TLI)                       0.860       0.778
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.094       0.056
#Robust RMSEA                                                  NA
#SRMR                                           0.112       0.112

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .377

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.302    0.065    4.661    0.000    0.302    0.302
#    factor3           0.311    0.061    5.098    0.000    0.311    0.311
#    factor4          -0.391    0.059   -6.602    0.000   -0.391   -0.391
#    factor5           0.448    0.060    7.498    0.000    0.448    0.448
#  factor2 ~~                                                            
#    factor3           0.545    0.057    9.551    0.000    0.545    0.545
#    factor4          -0.327    0.069   -4.722    0.000   -0.327   -0.327
#    factor5           0.295    0.071    4.173    0.000    0.295    0.295
#  factor3 ~~                                                            
#    factor4          -0.275    0.065   -4.260    0.000   -0.275   -0.275
#    factor5           0.087    0.073    1.183    0.237    0.087    0.087
#  factor4 ~~                                                            
#    factor5          -0.046    0.076   -0.602    0.547   -0.046   -0.046

CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.616       0.615
#Tucker-Lewis Index (TLI)                       0.601       0.600
#Robust Comparative Fit Index (CFI)                         0.620
#Robust Tucker-Lewis Index (TLI)                            0.604
#RMSEA                                          0.072       0.070
#Robust RMSEA                                               0.071
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .304

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.237    0.121    1.968    0.049    0.237    0.237
#    factor3           0.272    0.095    2.869    0.004    0.272    0.272
#    factor4          -0.379    0.079   -4.802    0.000   -0.379   -0.379
#    factor5           0.426    0.082    5.221    0.000    0.426    0.426
#  factor2 ~~                                                            
#    factor3           0.512    0.090    5.654    0.000    0.512    0.512
#    factor4          -0.307    0.089   -3.463    0.001   -0.307   -0.307
#    factor5           0.280    0.100    2.810    0.005    0.280    0.280
#  factor3 ~~                                                            
#    factor4          -0.257    0.083   -3.101    0.002   -0.257   -0.257
#    factor5           0.063    0.109    0.581    0.561    0.063    0.063
#  factor4 ~~                                                            
#    factor5          -0.022    0.093   -0.237    0.813   -0.022   -0.022

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.592       0.591
#Tucker-Lewis Index (TLI)                       0.577       0.577
#Robust Comparative Fit Index (CFI)                         0.595
#Robust Tucker-Lewis Index (TLI)                            0.581
#RMSEA                                          0.074       0.072
#Robust RMSEA                                               0.073
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307


#bifactor model better
BIFmodel= '
factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.687       0.679
#Tucker-Lewis Index (TLI)                       0.664       0.656
#Robust Comparative Fit Index (CFI)                         0.687
#Robust Tucker-Lewis Index (TLI)                            0.664
#RMSEA                                          0.066       0.065
#Robust RMSEA                                               0.065
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .349

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    BFI1              0.865    0.092    9.370    0.000    0.865    0.676
#    BFI6              0.430    0.108    3.961    0.000    0.430    0.351
#    BFI11             0.258    0.114    2.260    0.024    0.258    0.239
#    BFI16             0.942    0.084   11.197    0.000    0.942    0.800
#    BFI21             0.714    0.100    7.146    0.000    0.714    0.574
#    BFI26             0.285    0.108    2.646    0.008    0.285    0.221
#    BFI31             0.805    0.074   10.815    0.000    0.805    0.730
#    BFI36             0.270    0.099    2.729    0.006    0.270    0.246
#    BFI41             0.369    0.108    3.424    0.001    0.369    0.321
#    BFI46             0.785    0.090    8.676    0.000    0.785    0.640
#    BFI51             0.607    0.109    5.565    0.000    0.607    0.518
#    BFI56             0.403    0.092    4.377    0.000    0.403    0.405
#  factor2 =~                                                            
#    BFI2              0.301    0.089    3.374    0.001    0.301    0.385
#    BFI7              0.141    0.064    2.202    0.028    0.141    0.243
#    BFI12             0.443    0.110    4.009    0.000    0.443    0.397
#    BFI17             0.334    0.151    2.210    0.027    0.334    0.221
#    BFI22             0.630    0.129    4.895    0.000    0.630    0.565
#    BFI27             0.447    0.124    3.597    0.000    0.447    0.396
#    BFI32             0.159    0.088    1.803    0.071    0.159    0.163
#    BFI37             0.677    0.098    6.881    0.000    0.677    0.610
#    BFI42             0.269    0.108    2.500    0.012    0.269    0.257
#    BFI47             0.657    0.127    5.182    0.000    0.657    0.549
#    BFI52             0.218    0.066    3.284    0.001    0.218    0.359
#    BFI57             0.480    0.143    3.364    0.001    0.480    0.423
#  factor3 =~                                                            
#    BFI3              0.789    0.093    8.518    0.000    0.789    0.645
#    BFI8             -0.696    0.101   -6.911    0.000   -0.696   -0.599
#    BFI13             0.122    0.105    1.160    0.246    0.122    0.138
#    BFI18             0.507    0.091    5.557    0.000    0.507    0.490
#    BFI23             0.576    0.098    5.870    0.000    0.576    0.509
#    BFI28             0.653    0.101    6.433    0.000    0.653    0.587
#    BFI33             0.626    0.100    6.253    0.000    0.626    0.561
#    BFI38             0.391    0.092    4.259    0.000    0.391    0.415
#    BFI43             0.134    0.097    1.391    0.164    0.134    0.161
#    BFI48             0.634    0.088    7.235    0.000    0.634    0.584
#    BFI53             0.374    0.115    3.247    0.001    0.374    0.332
#    BFI58             0.632    0.114    5.556    0.000    0.632    0.550
#  factor4 =~                                                            
#    BFI4              0.818    0.087    9.377    0.000    0.818    0.616
#    BFI9              0.386    0.094    4.106    0.000    0.386    0.320
#    BFI14             0.794    0.084    9.444    0.000    0.794    0.654
#    BFI19             0.606    0.068    8.881    0.000    0.606    0.586
#    BFI24             0.554    0.087    6.372    0.000    0.554    0.455
#    BFI29             0.904    0.075   12.006    0.000    0.904    0.713
#    BFI34             0.857    0.081   10.604    0.000    0.857    0.728
#    BFI39             0.883    0.087   10.113    0.000    0.883    0.742
#    BFI44             0.575    0.083    6.932    0.000    0.575    0.522
#    BFI49             0.644    0.079    8.171    0.000    0.644    0.623
#    BFI54             0.948    0.087   10.895    0.000    0.948    0.733
#    BFI59             0.799    0.088    9.049    0.000    0.799    0.641
#  factor5 =~                                                            
#    BFI5              0.822    0.099    8.281    0.000    0.822    0.617
#    BFI10             0.181    0.059    3.092    0.002    0.181    0.238
#    BFI15             0.307    0.079    3.909    0.000    0.307    0.328
#    BFI20             0.711    0.091    7.810    0.000    0.711    0.613
#    BFI25             0.505    0.092    5.503    0.000    0.505    0.428
#    BFI30             0.752    0.094    7.994    0.000    0.752    0.627
#    BFI35             0.522    0.075    6.943    0.000    0.522    0.548
#    BFI40             0.350    0.071    4.964    0.000    0.350    0.381
#    BFI45             0.563    0.087    6.470    0.000    0.563    0.519
#    BFI50             0.602    0.100    6.044    0.000    0.602    0.468
#    BFI55             0.614    0.081    7.587    0.000    0.614    0.566
#    BFI60             0.536    0.087    6.145    0.000    0.536    0.501
#  global =~                                                             
#    BFI1              0.412    0.148    2.778    0.005    0.412    0.323
#    BFI6              0.513    0.137    3.739    0.000    0.513    0.419
#    BFI11             0.307    0.150    2.046    0.041    0.307    0.286
#    BFI16             0.055    0.146    0.374    0.708    0.055    0.046
#    BFI21             0.439    0.128    3.435    0.001    0.439    0.352
#    BFI26             0.326    0.135    2.417    0.016    0.326    0.254
#    BFI31             0.139    0.112    1.244    0.213    0.139    0.126
#    BFI36             0.297    0.098    3.022    0.003    0.297    0.271
#    BFI41             0.609    0.113    5.392    0.000    0.609    0.531
#    BFI46             0.231    0.142    1.629    0.103    0.231    0.188
#    BFI51             0.275    0.105    2.615    0.009    0.275    0.235
#    BFI56             0.585    0.117    4.991    0.000    0.585    0.587
#    BFI2              0.306    0.089    3.423    0.001    0.306    0.391
#    BFI7              0.193    0.049    3.971    0.000    0.193    0.332
#    BFI12             0.345    0.128    2.687    0.007    0.345    0.310
#    BFI17             0.059    0.172    0.345    0.730    0.059    0.039
#    BFI22            -0.058    0.129   -0.449    0.654   -0.058   -0.052
#    BFI27             0.484    0.113    4.283    0.000    0.484    0.428
#    BFI32             0.323    0.076    4.224    0.000    0.323    0.330
#    BFI37             0.227    0.152    1.497    0.134    0.227    0.204
#    BFI42            -0.086    0.134   -0.639    0.523   -0.086   -0.082
#    BFI47             0.085    0.216    0.393    0.694    0.085    0.071
#    BFI52             0.193    0.070    2.750    0.006    0.193    0.319
#    BFI57             0.439    0.132    3.332    0.001    0.439    0.387
#    BFI3              0.017    0.127    0.133    0.894    0.017    0.014
#    BFI8             -0.313    0.149   -2.105    0.035   -0.313   -0.270
#    BFI13             0.493    0.082    6.011    0.000    0.493    0.559
#    BFI18             0.181    0.094    1.919    0.055    0.181    0.175
#    BFI23             0.231    0.154    1.497    0.134    0.231    0.204
#    BFI28             0.070    0.150    0.471    0.637    0.070    0.063
#    BFI33             0.277    0.107    2.583    0.010    0.277    0.249
#    BFI38             0.537    0.078    6.863    0.000    0.537    0.570
#    BFI43             0.512    0.078    6.567    0.000    0.512    0.615
#    BFI48            -0.031    0.126   -0.242    0.809   -0.031   -0.028
#    BFI53             0.612    0.097    6.304    0.000    0.612    0.543
#    BFI58             0.075    0.148    0.505    0.614    0.075    0.065
#    BFI4             -0.456    0.122   -3.727    0.000   -0.456   -0.343
#    BFI9             -0.514    0.124   -4.133    0.000   -0.514   -0.426
#    BFI14            -0.155    0.149   -1.038    0.299   -0.155   -0.127
#    BFI19             0.078    0.126    0.621    0.535    0.078    0.076
#    BFI24            -0.499    0.105   -4.748    0.000   -0.499   -0.410
#    BFI29            -0.423    0.146   -2.890    0.004   -0.423   -0.334
#    BFI34            -0.103    0.118   -0.876    0.381   -0.103   -0.088
#    BFI39            -0.202    0.137   -1.472    0.141   -0.202   -0.169
#    BFI44            -0.402    0.116   -3.454    0.001   -0.402   -0.365
#    BFI49            -0.158    0.099   -1.587    0.113   -0.158   -0.153
#    BFI54            -0.232    0.148   -1.570    0.116   -0.232   -0.180
#    BFI59            -0.088    0.145   -0.604    0.546   -0.088   -0.070
#    BFI5             -0.055    0.145   -0.381    0.703   -0.055   -0.041
#    BFI10             0.248    0.072    3.459    0.001    0.248    0.325
#    BFI15             0.348    0.084    4.162    0.000    0.348    0.372
#    BFI20             0.031    0.119    0.257    0.797    0.031    0.026
#    BFI25             0.082    0.113    0.723    0.470    0.082    0.070
#    BFI30             0.147    0.137    1.073    0.283    0.147    0.123
#    BFI35             0.244    0.093    2.622    0.009    0.244    0.257
#    BFI40             0.309    0.075    4.138    0.000    0.309    0.336
#    BFI45             0.135    0.097    1.387    0.165    0.135    0.124
#    BFI50             0.202    0.135    1.501    0.133    0.202    0.157
#    BFI55             0.092    0.121    0.760    0.447    0.092    0.085
#    BFI60             0.457    0.083    5.488    0.000    0.457    0.427

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2539226     0.8135593     0.8586922     0.2761630 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.7017326 0.15160665 0.2982674 0.8601463 0.6145646 0.8559916 0.9209719
#factor2 0.6715356 0.09314252 0.3284644 0.7552756 0.5621092 0.7233934 0.8485444
#factor3 0.6553453 0.13780184 0.3446547 0.7793331 0.5572763 0.8094962 0.8966246
#factor4 0.8485840 0.22107969 0.1514160 0.9028046 0.8026897 0.8948619 0.9441091
#factor5 0.8165399 0.14244673 0.1834601 0.8245611 0.7170017 0.8120682 0.9000530
#global  0.2539226 0.25392258 0.2539226 0.8586922 0.2761630 0.8711937 0.9182402


#Reverse code items to see whether the combination of negative and positive values drives the low
#omegaH estimates
mydata2 <- mydata
reverse_cols = c("BFI11","BFI16","BFI26","BFI31","BFI36","BFI51","BFI12","BFI17","BFI22","BFI37","BFI42","BFI47",
                 "BFI13","BFI18","BFI33","BFI38","BFI43","BFI53","BFI14","BFI19","BFI34","BFI39","BFI54","BFI59",
                 "BFI10","BFI15","BFI20","BFI35","BFI40","BFI60")
mydata2[ , reverse_cols] = 6 - mydata2[ , reverse_cols]
CFA_model3 <- cfa(BIFmodel, data = mydata2,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#   0.23794178    0.81355932    0.86151245    0.03876786 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.6363225 0.1329343 0.3636775 0.8834168 0.58305295 0.8610867 0.9227172
#factor2 0.7811013 0.1142663 0.2188987 0.8108545 0.70517617 0.8034231 0.8967681
#factor3 0.7345918 0.1510643 0.2654082 0.8856263 0.69212999 0.8685894 0.9281452
#factor4 0.8486938 0.2189751 0.1513062 0.9254058 0.81830117 0.9171121 0.9557325
#factor5 0.7994695 0.1448182 0.2005305 0.8631676 0.73604578 0.8464345 0.9179067
#global  0.2379418 0.2379418 0.2379418 0.8615125 0.03876786 0.8791829 0.9212904






#############################################################
### BFI2 Andrejevic et al. (2022) https://journals.sagepub.com/doi/full/10.1177/19485506211038295
### data https://osf.io/27HUQ/
###

BFI2.Andrejevic <- read.csv("BFI2-Andrejevic.csv")
colnames(BFI2.Andrejevic)
mydata <- as.data.frame(BFI2.Andrejevic[,5:64])
mydata <- na.omit(mydata)
library(readxl)
BFI2_data_Vermeiren <- read_excel("BFI2_Vermeiren.xlsx")
colnames(mydata) <- paste0('BFI', 1:(ncol(mydata)))
head(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
omega(mydata) # alpha = .9, omega = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 7 components
# Eigenvalue1 = 8.36, eigenvalue2 = 4.59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 8 components 
# Eigenvalue 1 = 9.69; eigenvalue 2 = 5.31


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.093, RMSR=.13, TLI=.271

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.16, RMSEA=.122, RMSR=.15, TLI=.202


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.37, RMSEA=.056, RMSR=.05, TLI=.733
#      MR1   MR2   MR3   MR4   MR5
#MR1  1.00 -0.17  0.17  0.28  0.00
#MR2 -0.17  1.00 -0.17  0.01  0.18
#MR3  0.17 -0.17  1.00  0.09 -0.15
#MR4  0.28  0.01  0.09  1.00 -0.05
#MR5  0.00  0.18 -0.15 -0.05  1.00


fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.43, RMSEA=.09, RMSR=.06, TLI=.561
#     MR2   MR3   MR5   MR1   MR4
#MR1  1.00 -0.16  0.18  0.28 -0.01
#MR2 -0.16  1.00 -0.15  0.03  0.16
#MR5  0.18 -0.15  1.00  0.11 -0.17
#MR4  0.28  0.03  0.11  1.00 -0.07
#MR3 -0.01  0.16 -0.17 -0.07  1.00


# Single factor model
UNImodel= '
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'

CFA_model4 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.602       0.417
#Tucker-Lewis Index (TLI)                       0.588       0.396
#Robust Comparative Fit Index (CFI)                         0.250
#Robust Tucker-Lewis Index (TLI)                            0.223
#RMSEA                                          0.166       0.100
#Robust RMSEA                                               0.126
#SRMR                                           0.147       0.147

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .15


CFA_model4 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.292       0.285
#Tucker-Lewis Index (TLI)                       0.267       0.260
#Robust Comparative Fit Index (CFI)                         0.293
#Robust Tucker-Lewis Index (TLI)                            0.268
#RMSEA                                          0.098       0.095
#Robust RMSEA                                               0.097
#SRMR                                           0.124       0.124

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .117

semPaths(CFA_model3,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with theory-based structure
EGAmodel= '
 factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
 factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
 factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
 factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
 factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'

library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.871       0.788
#Tucker-Lewis Index (TLI)                       0.866       0.779
#Robust Comparative Fit Index (CFI)                         0.574
#Robust Tucker-Lewis Index (TLI)                            0.557
#RMSEA                                          0.095       0.060
#Robust RMSEA                                               0.095
#SRMR                                           0.101       0.101

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .40

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.133    0.030    4.387    0.000    0.251    0.251
#  factor3          -0.214    0.033   -6.396    0.000   -0.333   -0.333
#  factor4           0.165    0.030    5.510    0.000    0.312    0.312
#  factor5          -0.165    0.026   -6.381    0.000   -0.391   -0.391
#factor2 ~~                                                            
#  factor3          -0.217    0.032   -6.793    0.000   -0.438   -0.438
#  factor4           0.114    0.023    5.007    0.000    0.281    0.281
#  factor5          -0.064    0.020   -3.198    0.001   -0.196   -0.196
#factor3 ~~                                                            
#  factor4          -0.147    0.028   -5.331    0.000   -0.298   -0.298
#  factor5           0.056    0.022    2.528    0.011    0.142    0.142
#factor4 ~~                                                            
#  factor5          -0.007    0.019   -0.363    0.717   -0.021   -0.021

CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.658       0.660
#Tucker-Lewis Index (TLI)                       0.644       0.646
#Robust Comparative Fit Index (CFI)                         0.664
#Robust Tucker-Lewis Index (TLI)                            0.651
#RMSEA                                          0.068       0.065
#Robust RMSEA                                               0.067
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .224

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.206    0.082    2.499    0.012    0.206    0.206
#    factor3          -0.253    0.076   -3.315    0.001   -0.253   -0.253
#    factor4           0.288    0.072    4.004    0.000    0.288    0.288
#    factor5          -0.350    0.063   -5.568    0.000   -0.350   -0.350
#  factor2 ~~                                                            
#    factor3          -0.391    0.072   -5.452    0.000   -0.391   -0.391
#    factor4           0.295    0.076    3.878    0.000    0.295    0.295
#    factor5          -0.173    0.078   -2.232    0.026   -0.173   -0.173
#  factor3 ~~                                                            
#    factor4          -0.275    0.068   -4.018    0.000   -0.275   -0.275
#    factor5           0.091    0.069    1.311    0.190    0.091    0.091
#  factor4 ~~                                                            
#    factor5          -0.027    0.073   -0.375    0.707   -0.027   -0.027

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.642       0.644
#Tucker-Lewis Index (TLI)                       0.629       0.632
#Robust Comparative Fit Index (CFI)                         0.648
#Robust Tucker-Lewis Index (TLI)                            0.636
#RMSEA                                          0.070       0.067
#Robust RMSEA                                               0.068
#SRMR                                           0.120       0.120

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .324


#bifactor model better?
BIFmodel= '
factor1=~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56
factor2=~ BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57
factor3=~ BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58
factor4=~ BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59
factor5=~ BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
global =~ BFI1 + BFI6 + BFI11 + BFI16 + BFI21 + BFI26 + BFI31 + BFI36 + BFI41 + BFI46 + BFI51 + BFI56 +
          BFI2 + BFI7 + BFI12 + BFI17 + BFI22 + BFI27 + BFI32 + BFI37 + BFI42 + BFI47 + BFI52 + BFI57 +
          BFI3 + BFI8 + BFI13 + BFI18 + BFI23 + BFI28 + BFI33 + BFI38 + BFI43 + BFI48 + BFI53 + BFI58 +
          BFI4 + BFI9 + BFI14 + BFI19 + BFI24 + BFI29 + BFI34 + BFI39 + BFI44 + BFI49 + BFI54 + BFI59 +
          BFI5 + BFI10 + BFI15 + BFI20 + BFI25 + BFI30 + BFI35 + BFI40 + BFI45 + BFI50 + BFI55 + BFI60
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.717       0.712
#Tucker-Lewis Index (TLI)                       0.696       0.691
#Robust Comparative Fit Index (CFI)                         0.720
#Robust Tucker-Lewis Index (TLI)                            0.699
#RMSEA                                          0.063       0.061
#Robust RMSEA                                               0.062
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .335

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2737761     0.8135593     0.3405045     0.3357389 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS      Omega       OmegaH         H        FD
#factor1 0.4309606 0.1028790 0.5690394 0.11951364 6.478784e-03 0.8031857 0.8864409
#factor2 0.8557522 0.1275298 0.1442478 0.20295916 4.029317e-03 0.7925168 0.8923038
#factor3 0.7938955 0.1704769 0.2061045 0.13099473 1.694929e-05 0.8762063 0.9341325
#factor4 0.8772696 0.1883551 0.1227304 0.05935753 1.289065e-02 0.8714213 0.9323587
#factor5 0.7493095 0.1369832 0.2506905 0.02451311 9.896632e-03 0.8240956 0.9028330
#global  0.2737761 0.2737761 0.2737761 0.34050451 3.357389e-01 0.8818048 0.9140725







####### HEXACO Henry
### https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0262465
### Time 1
library(readxl)
HEXACO_Henry_time1 <- read_excel("HEXACO_Henry_time1.xlsx")
mydata <- as.data.frame(HEXACO_Henry_time1[,1:96])
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 12 factors and 12 components
# Eigenvalue1 = 7.59, eigenvalue2 = 5.96

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 12 components 
# Eigenvalue 1 = 8.78; eigenvalue 2 = 7.00


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.08, RMSEA=.07, RMSR=.11, TLI=.178

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.09, RMSEA=.093, RMSR=.13, TLI=.125


# Give solution with 6 factors
fit5 <- fa(mydata,6)
fit5
diagram(fit5)
# %variance explained=.29, RMSEA=.045, RMSR=.05, TLI=.66
#      MR1   MR2  MR5   MR6   MR3   MR4
#MR1  1.00 -0.06 0.11  0.13 -0.05  0.06
#MR2 -0.06  1.00 0.16  0.20  0.06  0.06
#MR5  0.11  0.16 1.00  0.04  0.04  0.03
#MR6  0.13  0.20 0.04  1.00 -0.12 -0.01
#MR3 -0.05  0.06 0.04 -0.12  1.00  0.05
#MR4  0.06  0.06 0.03 -0.01  0.05  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.31, RMSEA=.074, RMSR=.06, TLI=.436
#     MR1   MR2  MR4   MR3  MR5
#MR1  1.00  0.03 0.07 -0.12 0.04
#MR2  0.03  1.00 0.14 -0.05 0.05
#MR4  0.07  0.14 1.00  0.06 0.09
#MR3 -0.12 -0.05 0.06  1.00 0.09
#MR5  0.04  0.05 0.09  0.09 1.00

UNImodel= '
 global =~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+
           C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16+
           E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16+
           H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16+
           O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16+
           X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
'
library(lavaan)
CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.392       0.242
#Tucker-Lewis Index (TLI)                       0.379       0.226
#Robust Comparative Fit Index (CFI)                         0.152
#Robust Tucker-Lewis Index (TLI)                            0.134
#RMSEA                                          0.135       0.068
#Robust RMSEA                                               0.098
#SRMR                                           0.126       0.126

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .084

CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.192       0.184
#Tucker-Lewis Index (TLI)                       0.175       0.167
#Robust Comparative Fit Index (CFI)                         0.192
#Robust Tucker-Lewis Index (TLI)                            0.175
#RMSEA                                          0.075       0.072
#Robust RMSEA                                               0.074
#SRMR                                           0.109       0.109

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .034


EGAmodel= '
 factor1=~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16
 factor2=~ C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16
 factor3=~ E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16
 factor4=~ H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16
 factor5=~ O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16
 factor6=~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
'

library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.758       0.682
#Tucker-Lewis Index (TLI)                       0.752       0.674
#Robust Comparative Fit Index (CFI)                         0.454
#Robust Tucker-Lewis Index (TLI)                            0.440
#RMSEA                                          0.086       0.044
#Robust RMSEA                                               0.079
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .30

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.022    0.012    1.936    0.053    0.095    0.095
#  factor3          -0.079    0.016   -4.819    0.000   -0.286   -0.286
#  factor4           0.095    0.017    5.692    0.000    0.273    0.273
#  factor5           0.036    0.025    1.445    0.149    0.081    0.081
#  factor6           0.118    0.023    5.209    0.000    0.280    0.280
#factor2 ~~                                                            
#  factor3           0.006    0.009    0.652    0.514    0.033    0.033
#  factor4           0.057    0.013    4.411    0.000    0.265    0.265
#  factor5           0.035    0.014    2.598    0.009    0.127    0.127
#  factor6           0.079    0.017    4.630    0.000    0.299    0.299
#factor3 ~~                                                            
#  factor4           0.055    0.013    4.268    0.000    0.218    0.218
#  factor5           0.021    0.016    1.271    0.204    0.064    0.064
#  factor6          -0.104    0.017   -6.128    0.000   -0.337   -0.337
#factor4 ~~                                                            
#  factor5           0.070    0.019    3.700    0.000    0.171    0.171
#  factor6          -0.019    0.019   -1.012    0.311   -0.050   -0.050
#factor5 ~~                                                            
#  factor6           0.084    0.025    3.337    0.001    0.169    0.169

CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.547       0.553
#Tucker-Lewis Index (TLI)                       0.535       0.542
#Robust Comparative Fit Index (CFI)                         0.556
#Robust Tucker-Lewis Index (TLI)                            0.545
#RMSEA                                          0.056       0.053
#Robust RMSEA                                               0.055
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .244

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.081    0.071    1.135    0.256    0.081    0.081
#    factor3          -0.292    0.069   -4.199    0.000   -0.292   -0.292
#    factor4           0.289    0.059    4.915    0.000    0.289    0.289
#    factor5           0.070    0.072    0.976    0.329    0.070    0.070
#    factor6           0.260    0.065    4.009    0.000    0.260    0.260
#  factor2 ~~                                                            
#    factor3           0.040    0.071    0.558    0.577    0.040    0.040
#    factor4           0.252    0.071    3.567    0.000    0.252    0.252
#    factor5           0.097    0.063    1.532    0.126    0.097    0.097
#    factor6           0.301    0.070    4.275    0.000    0.301    0.301
#  factor3 ~~                                                            
#    factor4           0.207    0.063    3.303    0.001    0.207    0.207
#    factor5           0.061    0.074    0.832    0.406    0.061    0.061
#    factor6          -0.287    0.086   -3.319    0.001   -0.287   -0.287
#  factor4 ~~                                                            
#    factor5           0.181    0.059    3.077    0.002    0.181    0.181
#    factor6          -0.079    0.081   -0.977    0.329   -0.079   -0.079
#  factor5 ~~                                                            
#    factor6           0.137    0.064    2.131    0.033    0.137    0.137

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.533       0.540
#Tucker-Lewis Index (TLI)                       0.523       0.531
#Robust Comparative Fit Index (CFI)                         0.543
#Robust Tucker-Lewis Index (TLI)                            0.533
#RMSEA                                          0.057       0.054
#Robust RMSEA                                               0.055
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .237


#bifactor model 
BIFmodel= '
 factor1=~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16
 factor2=~ C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16
 factor3=~ E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16
 factor4=~ H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16
 factor5=~ O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16
 factor6=~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
 global =~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+
           C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16+
           E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16+
           H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16+
           O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16+
           X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.610       0.617
#Tucker-Lewis Index (TLI)                       0.593       0.600
#Robust Comparative Fit Index (CFI)                         0.620
#Robust Tucker-Lewis Index (TLI)                            0.604
#RMSEA                                          0.052       0.050
#Robust RMSEA                                               0.051
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .262

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2638986     0.8421053     0.8545441     0.2793617 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#factor1 0.9029227 0.13753182 0.09707725 0.8329858 0.769958912 0.8310991 0.9117625
#factor2 0.8847731 0.14456778 0.11522694 0.8403770 0.792073742 0.8393866 0.9178304
#factor3 0.8016825 0.12810074 0.19831747 0.8201531 0.781070279 0.8180977 0.9079910
#factor4 0.9395558 0.15574475 0.06044423 0.8453916 0.842708703 0.8523857 0.9249031
#factor5 0.9550216 0.12437869 0.04497839 0.7937899 0.775128175 0.8235869 0.9078005
#factor6 0.2003435 0.04577761 0.79965649 0.8800351 0.004055442 0.5777983 0.8382801
#global  0.2638986 0.26389862 0.26389862 0.8545441 0.279361737 0.9072490 0.9579082





####### RIASEC openpsychometrics
### http://openpsychometrics.org/_rawdata/

data_RIASEC_openpsychometrics <- read.delim("data_RIASEC_openpsychometrics.csv")
mydata <- as.data.frame(data_RIASEC_openpsychometrics[,1:48])
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) # response alternatives = 6

# take random sample of 1000 participants
library(dplyr)
mydata2 <- mydata
set.seed(6465)
mydata <- sample_n(mydata2,1000)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 7 components
# Eigenvalue1 = 9.11, eigenvalue2 = 4.51

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 7 components 
# Eigenvalue 1 = 10.83; eigenvalue 2 = 4.89


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.121, RMSR=.14, TLI=.296

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.14, RMSR=.16, TLI=.284


# Give solution with 6 factors
fit5 <- fa(mydata,6)
fit5
diagram(fit5)
# %variance explained=.49, RMSEA=.059, RMSR=.03, TLI=.832
#      MR1  MR3  MR4   MR2  MR6  MR5
#MR1  1.00 0.11 0.36 -0.01 0.42 0.15
#MR3  0.11 1.00 0.27  0.27 0.01 0.13
#MR4  0.36 0.27 1.00  0.16 0.26 0.03
#MR2 -0.01 0.27 0.16  1.00 0.25 0.30
M#R6  0.42 0.01 0.26  0.25 1.00 0.29
M#R5  0.15 0.13 0.03  0.30 0.29 1.00

fit6 <- fa(rho,6,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.55, RMSEA=.078, RMSR=.03, TLI=.778
#     MR1  MR5  MR4  MR2  MR6  MR3
#MR1 1.00 0.39 0.12 0.00 0.42 0.16
#MR5 0.39 1.00 0.30 0.18 0.28 0.02
#MR4 0.12 0.30 1.00 0.28 0.01 0.13
#MR2 0.00 0.18 0.28 1.00 0.26 0.30
#MR6 0.42 0.28 0.01 0.26 1.00 0.29
#MR3 0.16 0.02 0.13 0.30 0.29 1.00

#Analysis lavaan
library(lavaan)

UNImodel= '
 global =~ R1+R2+R3+R4+R5+R6+R7+R8+
           I1+I2+I3+I4+I5+I6+I7+I8+
           A1+A2+A3+A4+A5+A6+A7+A8+
           S1+S2+S3+S4+S5+S6+S7+S8+
           E1+E2+E3+E4+E5+E6+E7+E8+
           C1+C2+C3+C4+C5+C6+C7+C8
'
CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.691       0.367
#Tucker-Lewis Index (TLI)                       0.678       0.339
#Robust Comparative Fit Index (CFI)                         0.297
#Robust Tucker-Lewis Index (TLI)                            0.266
#RMSEA                                          0.215       0.153
#Robust RMSEA                                               0.143
#SRMR                                           0.166       0.166

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .265

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.328       0.321
#Tucker-Lewis Index (TLI)                       0.298       0.291
#Robust Comparative Fit Index (CFI)                         0.330
#Robust Tucker-Lewis Index (TLI)                            0.300
#RMSEA                                          0.122       0.113
#Robust RMSEA                                               0.121
#SRMR                                           0.141       0.141

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .151

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


EGAmodel= '
 factor1=~ R1+R2+R3+R4+R5+R6+R7+R8
 factor2=~ I1+I2+I3+I4+I5+I6+I7+I8
 factor3=~ A1+A2+A3+A4+A5+A6+A7+A8
 factor4=~ S1+S2+S3+S4+S5+S6+S7+S8
 factor5=~ E1+E2+E3+E4+E5+E6+E7+E8
 factor6=~ C1+C2+C3+C4+C5+C6+C7+C8
'

CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.961       0.878
#Tucker-Lewis Index (TLI)                       0.959       0.870
#Robust Comparative Fit Index (CFI)                         0.760
#Robust Tucker-Lewis Index (TLI)                            0.745
#RMSEA                                          0.086       0.074
#Robust RMSEA                                               0.087
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .55

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.363    0.028   12.806    0.000    0.363    0.363
#    factor3           0.238    0.032    7.372    0.000    0.238    0.238
#    factor4           0.169    0.033    5.178    0.000    0.169    0.169
#    factor5           0.453    0.028   16.462    0.000    0.453    0.453
#    factor6           0.599    0.021   27.876    0.000    0.599    0.599
#  factor2 ~~                                                            
#    factor3           0.362    0.028   13.164    0.000    0.362    0.362
#    factor4           0.228    0.031    7.456    0.000    0.228    0.228
#    factor5           0.126    0.031    4.054    0.000    0.126    0.126
#    factor6           0.177    0.031    5.752    0.000    0.177    0.177
#  factor3 ~~                                                            
#    factor4           0.427    0.027   15.678    0.000    0.427    0.427
#    factor5           0.421    0.027   15.379    0.000    0.421    0.421
#    factor6           0.104    0.033    3.141    0.002    0.104    0.104
#  factor4 ~~                                                            
#    factor5           0.484    0.025   19.030    0.000    0.484    0.484
#    factor6           0.285    0.029    9.698    0.000    0.285    0.285
#  factor5 ~~                                                            
#    factor6           0.669    0.018   38.056    0.000    0.669    0.669

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.804       0.807
#Tucker-Lewis Index (TLI)                       0.793       0.796
#Robust Comparative Fit Index (CFI)                         0.810
#Robust Tucker-Lewis Index (TLI)                            0.798
#RMSEA                                          0.066       0.061
#Robust RMSEA                                               0.065
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .486

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.310    0.034    9.070    0.000    0.310    0.310
#    factor3           0.228    0.038    6.000    0.000    0.228    0.228
#    factor4           0.104    0.040    2.561    0.010    0.104    0.104
#    factor5           0.377    0.034   11.014    0.000    0.377    0.377
#    factor6           0.521    0.034   15.456    0.000    0.521    0.521
#  factor2 ~~                                                            
#    factor3           0.296    0.036    8.160    0.000    0.296    0.296
#    factor4           0.254    0.037    6.795    0.000    0.254    0.254
#    factor5           0.069    0.038    1.813    0.070    0.069    0.069
#    factor6           0.127    0.037    3.451    0.001    0.127    0.127
#  factor3 ~~                                                            
#    factor4           0.344    0.035    9.727    0.000    0.344    0.344
#    factor5           0.298    0.038    7.811    0.000    0.298    0.298
#    factor6           0.038    0.039    0.976    0.329    0.038    0.038
#  factor4 ~~                                                            
#    factor5           0.373    0.037   10.037    0.000    0.373    0.373
#    factor6           0.203    0.038    5.347    0.000    0.203    0.203
#  factor5 ~~                                                            
#    factor6           0.643    0.028   23.065    0.000    0.643    0.643

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.760       0.761
#Tucker-Lewis Index (TLI)                       0.749       0.750
#Robust Comparative Fit Index (CFI)                         0.765
#Robust Tucker-Lewis Index (TLI)                            0.754
#RMSEA                                          0.073       0.067
#Robust RMSEA                                               0.072
#SRMR                                           0.153       0.153

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .480


#bifactor model 
BIFmodel= '
 factor1=~ R1+R2+R3+R4+R5+R6+R7+R8
 factor2=~ I1+I2+I3+I4+I5+I6+I7+I8
 factor3=~ A1+A2+A3+A4+A5+A6+A7+A8
 factor4=~ S1+S2+S3+S4+S5+S6+S7+S8
 factor5=~ E1+E2+E3+E4+E5+E6+E7+E8
 factor6=~ C1+C2+C3+C4+C5+C6+C7+C8
 global =~ R1+R2+R3+R4+R5+R6+R7+R8+
           I1+I2+I3+I4+I5+I6+I7+I8+
           A1+A2+A3+A4+A5+A6+A7+A8+
           S1+S2+S3+S4+S5+S6+S7+S8+
           E1+E2+E3+E4+E5+E6+E7+E8+
           C1+C2+C3+C4+C5+C6+C7+C8
'

CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.816       0.816
#Tucker-Lewis Index (TLI)                       0.798       0.799
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.803
#RMSEA                                          0.065       0.060
#Robust RMSEA                                               0.064
#SRMR                                           0.087       0.087

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .490

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
CV.global           PUC  Omega.global OmegaH.global 
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2934886     0.8510638     0.9367839     0.6248751 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3188069 0.05834661 0.68119309 0.8890597 0.2361649 0.6764574 0.8146822
#factor2 0.8912496 0.16108992 0.10875042 0.8932867 0.8025371 0.8893148 0.9430886
#factor3 0.9095288 0.14319751 0.09047124 0.8613487 0.7971013 0.8646470 0.9296837
#factor4 0.9427448 0.15038905 0.05725522 0.8668566 0.8261169 0.8765661 0.9361255
#factor5 0.6923193 0.09650967 0.30768072 0.8387874 0.5837988 0.7657409 0.8711257
#factor6 0.5391479 0.09697862 0.46085212 0.8919071 0.4749921 0.7708059 0.8667217
#global  0.2934886 0.29348864 0.29348864 0.9367839 0.6248751 0.9068837 0.9159327





####### 16PF openpsychometrics
### http://openpsychometrics.org/_rawdata/

### Not included because took too long

data_16PF <- read.delim("~/Self regulation/data_16PF.csv")
mydata <- as.data.frame(data_16PF)
mydata <- na.omit(mydata[,1:163])
summary(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 21 components
# Eigenvalue1 = 18.1, eigenvalue2 = 14.94

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 21 components 
# Eigenvalue 1 = 20.16; eigenvalue 2 = 16.09


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.11, RMSEA=.068, RMSR=.13, TLI=.215

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.12, RMSEA=.076, RMSR=.14, TLI=.199


# Give solution with 16 factors
fit5 <- fa(mydata,16)
fit5
diagram(fit5)
#      MR16   MR1   MR9  MR12   MR2   MR7   MR5   MR4   MR3   MR6  MR11  MR15   MR8 MR10 MR14  MR13
#MR16  1.00 -0.32  0.17 -0.22 -0.19 -0.10  0.12  0.04  0.18 -0.34  0.09  0.24  0.20 0.09 0.02 -0.08
#MR1  -0.32  1.00 -0.40  0.22  0.23  0.38  0.19  0.03 -0.07  0.03 -0.08 -0.05 -0.01 0.14 0.11  0.18
#MR9   0.17 -0.40  1.00 -0.24  0.15 -0.16  0.06  0.08  0.15  0.02  0.39  0.25  0.13 0.12 0.22 -0.04
#MR12 -0.22  0.22 -0.24  1.00 -0.04  0.37  0.31 -0.30  0.02  0.33  0.04 -0.07 -0.04 0.06 0.00  0.15
#MR2  -0.19  0.23  0.15 -0.04  1.00  0.20  0.12  0.09 -0.07  0.09  0.26  0.15 -0.13 0.12 0.22  0.15
#MR7  -0.10  0.38 -0.16  0.37  0.20  1.00  0.31 -0.05 -0.02  0.06  0.11  0.07  0.01 0.10 0.05  0.13
#MR5   0.12  0.19  0.06  0.31  0.12  0.31  1.00 -0.17  0.03  0.07  0.23  0.14  0.06 0.27 0.19  0.07
#MR4   0.04  0.03  0.08 -0.30  0.09 -0.05 -0.17  1.00 -0.11 -0.09  0.02 -0.04  0.33 0.15 0.03 -0.16
#MR3   0.18 -0.07  0.15  0.02 -0.07 -0.02  0.03 -0.11  1.00 -0.04 -0.06  0.24  0.06 0.11 0.16  0.37
#MR6  -0.34  0.03  0.02  0.33  0.09  0.06  0.07 -0.09 -0.04  1.00  0.13 -0.16 -0.11 0.03 0.10  0.11
#MR11  0.09 -0.08  0.39  0.04  0.26  0.11  0.23  0.02 -0.06  0.13  1.00  0.13  0.07 0.14 0.21 -0.07
#MR15  0.24 -0.05  0.25 -0.07  0.15  0.07  0.14 -0.04  0.24 -0.16  0.13  1.00  0.18 0.21 0.30  0.15
#MR8   0.20 -0.01  0.13 -0.04 -0.13  0.01  0.06  0.33  0.06 -0.11  0.07  0.18  1.00 0.27 0.20 -0.09
#MR10  0.09  0.14  0.12  0.06  0.12  0.10  0.27  0.15  0.11  0.03  0.14  0.21  0.27 1.00 0.22  0.09
#MR14  0.02  0.11  0.22  0.00  0.22  0.05  0.19  0.03  0.16  0.10  0.21  0.30  0.20 0.22 1.00  0.20
#MR13 -0.08  0.18 -0.04  0.15  0.15  0.13  0.07 -0.16  0.37  0.11 -0.07  0.15 -0.09 0.09 0.20  1.00

fit6 <- fa(rho,16,n.obs=nrow(mydata))
fit6
diagram(fit6)
#MR1   MR9   MR2  MR12   MR7   MR3   MR4  MR14   MR6   MR5  MR16  MR13   MR8  MR10  MR11  MR15
#MR1   1.00  0.16 -0.27 -0.21 -0.12  0.05 -0.05 -0.25 -0.37  0.11  0.10  0.27  0.28  0.04  0.15 -0.03
#MR9   0.16  1.00 -0.34 -0.25 -0.22  0.09 -0.05  0.00 -0.02 -0.03  0.32  0.21  0.17  0.31  0.12 -0.04
#MR2  -0.27 -0.34  1.00  0.20  0.38  0.15 -0.01  0.30  0.11  0.19  0.01 -0.01 -0.04 -0.08  0.09  0.21
#MR12 -0.21 -0.25  0.20  1.00  0.36  0.07  0.30 -0.05  0.33  0.32 -0.02 -0.05 -0.05  0.00  0.04  0.02
#MR7  -0.12 -0.22  0.38  0.36  1.00  0.00  0.03  0.22  0.05  0.31  0.08  0.03  0.02  0.07  0.08  0.13
#MR3   0.05  0.09  0.15  0.07  0.00  1.00  0.14  0.00  0.08  0.01  0.12  0.36  0.07 -0.18  0.17 -0.05
#MR4  -0.05 -0.05 -0.01  0.30  0.03  0.14  1.00 -0.09  0.07  0.20  0.04  0.10 -0.31 -0.05 -0.11 -0.09
#MR14 -0.25  0.00  0.30 -0.05  0.22  0.00 -0.09  1.00  0.07  0.06  0.28  0.02 -0.09  0.16  0.06  0.20
#MR6  -0.37 -0.02  0.11  0.33  0.05  0.08  0.07  0.07  1.00  0.05  0.10 -0.16 -0.07  0.07  0.01  0.04
#MR5   0.11 -0.03  0.19  0.32  0.31  0.01  0.20  0.06  0.05  1.00  0.21  0.05  0.07  0.17  0.23  0.17
#MR16  0.10  0.32  0.01 -0.02  0.08  0.12  0.04  0.28  0.10  0.21  1.00  0.22  0.12  0.33  0.19  0.18
#MR13  0.27  0.21 -0.01 -0.05  0.03  0.36  0.10  0.02 -0.16  0.05  0.22  1.00  0.22 -0.06  0.23  0.07
#MR8   0.28  0.17 -0.04 -0.05  0.02  0.07 -0.31 -0.09 -0.07  0.07  0.12  0.22  1.00  0.14  0.29  0.27
#MR10  0.04  0.31 -0.08  0.00  0.07 -0.18 -0.05  0.16  0.07  0.17  0.33 -0.06  0.14  1.00  0.09  0.12
#MR11  0.15  0.12  0.09  0.04  0.08  0.17 -0.11  0.06  0.01  0.23  0.19  0.23  0.29  0.09  1.00  0.20
#MR15 -0.03 -0.04  0.21  0.02  0.13 -0.05 -0.09  0.20  0.04  0.17  0.18  0.07  0.27  0.12  0.20  1.00


library(EGAnet); library(psych); library(foreign); library(ggplot2)
riEGALV <- riEGA(mydata, algorithm="louvain")
riEGALV
summary(riEGALV)
# aborted after 6 hours


#Analysis lavaan
library(lavaan)

UNImodel= '
 global =~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+
            B1+B2+B3+B4+B5+B6+B7+B8+B9+B10+B11+B12+B13+
            C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+
            D1+D2+D3+D4+D5+D6+D7+D8+D9+D10+
            E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+
            F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+
            G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+
            H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+
            I1+I2+I3+I4+I5+I6+I7+I8+I9+I10+
            J1+J2+J3+J4+J5+J6+J7+J8+J9+J10+
            K1+K2+K3+K4+K5+K6+K7+K8+K9+K10+
            L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+
            M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+
            N1+N2+N3+N4+N5+N6+N7+N8+N9+N10+
            O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+
            P1+P2+P3+P4+P5+P6+P7+P8+P9+P10
'

CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
# aborted after 36+ hours


#Comparative Fit Index (CFI)                    0.702       0.234
#Tucker-Lewis Index (TLI)                       0.688       0.200
#Robust Comparative Fit Index (CFI)                         0.304
#Robust Tucker-Lewis Index (TLI)                            0.273
#RMSEA                                          0.217       0.166
#Robust RMSEA                                               0.141
#SRMR                                           0.168       0.168


fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .55

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


EGAmodel= '
 factor1=~ R1+R2+R3+R4+R5+R6+R7+R8
 factor2=~ I1+I2+I3+I4+I5+I6+I7+I8
 factor3=~ A1+A2+A3+A4+A5+A6+A7+A8
 factor4=~ S1+S2+S3+S4+S5+S6+S7+S8
 factor5=~ E1+E2+E3+E4+E5+E6+E7+E8
 factor6=~ C1+C2+C3+C4+C5+C6+C7+C8
'

CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.959       0.841
#Tucker-Lewis Index (TLI)                       0.957       0.832
#Robust Comparative Fit Index (CFI)                         0.785
#Robust Tucker-Lewis Index (TLI)                            0.772
#RMSEA                                          0.081       0.076
#Robust RMSEA                                               0.079
#SRMR                                           0.068       0.068

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#factor2           0.359    0.002  154.963    0.000    0.359    0.359
#factor3           0.262    0.003  103.172    0.000    0.262    0.262
#factor4           0.133    0.003   50.358    0.000    0.133    0.133
#factor5           0.431    0.002  184.504    0.000    0.431    0.431
#factor6           0.593    0.002  325.162    0.000    0.593    0.593
#factor2 ~~                                                            
#factor3           0.374    0.002  166.797    0.000    0.374    0.374
#factor4           0.256    0.002  106.520    0.000    0.256    0.256
#factor5           0.110    0.003   43.393    0.000    0.110    0.110
#factor6           0.157    0.002   64.081    0.000    0.157    0.157
#factor3 ~~                                                            
#factor4           0.399    0.002  174.187    0.000    0.399    0.399
#factor5           0.393    0.002  168.918    0.000    0.393    0.393
#factor6           0.089    0.003   34.250    0.000    0.089    0.089
#factor4 ~~                                                            
#factor5           0.488    0.002  229.175    0.000    0.488    0.488
#factor6           0.257    0.002  104.369    0.000    0.257    0.257
#factor5 ~~                                                            
#factor6           0.651    0.002  385.132    0.000    0.651    0.651

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .55

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.715       0.679
#Tucker-Lewis Index (TLI)                       0.702       0.664
#Robust Comparative Fit Index (CFI)                         0.743
#Robust Tucker-Lewis Index (TLI)                            0.732
#RMSEA                                          0.212       0.107
#Robust RMSEA                                               0.086
#SRMR                                           0.186       0.186

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .55


#bifactor model 
BIFmodel= '
 factor1=~ R1+R2+R3+R4+R5+R6+R7+R8
 factor2=~ I1+I2+I3+I4+I5+I6+I7+I8
 factor3=~ A1+A2+A3+A4+A5+A6+A7+A8
 factor4=~ S1+S2+S3+S4+S5+S6+S7+S8
 factor5=~ E1+E2+E3+E4+E5+E6+E7+E8
 factor6=~ C1+C2+C3+C4+C5+C6+C7+C8
 global =~ R1+R2+R3+R4+R5+R6+R7+R8+
           I1+I2+I3+I4+I5+I6+I7+I8+
           A1+A2+A3+A4+A5+A6+A7+A8+
           S1+S2+S3+S4+S5+S6+S7+S8+
           E1+E2+E3+E4+E5+E6+E7+E8+
           C1+C2+C3+C4+C5+C6+C7+C8
'

CFA_model4 <- cfa(BIFmodel, data = mydata,ordered=TRUE,orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.939       0.808
#Tucker-Lewis Index (TLI)                       0.934       0.791
#Robust Comparative Fit Index (CFI)                         0.793
#Robust Tucker-Lewis Index (TLI)                            0.773
#RMSEA                                          0.100       0.085
#Robust RMSEA                                               0.079
#SRMR                                           0.087       0.087

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    R1                0.273    0.003  107.287    0.000    0.273    0.273
#    R2                0.680    0.002  354.159    0.000    0.680    0.680
#    R3                0.506    0.002  207.850    0.000    0.506    0.506
#    R4                0.501    0.002  237.958    0.000    0.501    0.501
#    R5                0.576    0.002  252.128    0.000    0.576    0.576
#    R6                0.635    0.002  329.469    0.000    0.635    0.635
#    R7                0.488    0.002  221.441    0.000    0.488    0.488
#    R8                0.689    0.002  359.395    0.000    0.689    0.689
#  factor2 =~                                                            
#    I1                0.616    0.002  371.193    0.000    0.616    0.616
#    I2                0.742    0.001  582.886    0.000    0.742    0.742
#    I3                0.810    0.001  741.214    0.000    0.810    0.810
#    I4                0.620    0.002  397.411    0.000    0.620    0.620
#    I5                0.854    0.001  833.443    0.000    0.854    0.854
#    I6                0.751    0.001  599.107    0.000    0.751    0.751
#    I7                0.824    0.001  768.550    0.000    0.824    0.824
#    I8                0.515    0.002  285.334    0.000    0.515    0.515
#  factor3 =~                                                            
#    A1                0.651    0.002  404.243    0.000    0.651    0.651
#    A2                0.689    0.002  450.684    0.000    0.689    0.689
#    A3                0.572    0.002  314.815    0.000    0.572    0.572
#    A4                0.765    0.001  564.020    0.000    0.765    0.765
#    A5                0.734    0.001  510.276    0.000    0.734    0.734
#    A6                0.694    0.002  452.550    0.000    0.694    0.694
#    A7                0.315    0.002  134.474    0.000    0.315    0.315
#    A8                0.658    0.002  411.290    0.000    0.658    0.658
#  factor4 =~                                                            
#    S1                0.489    0.002  246.904    0.000    0.489    0.489
#    S2                0.535    0.002  269.096    0.000    0.535    0.535
#    S3                0.756    0.001  544.381    0.000    0.756    0.756
#    S4                0.331    0.002  149.160    0.000    0.331    0.331
#    S5                0.834    0.001  637.492    0.000    0.834    0.834
#    S6                0.609    0.002  349.892    0.000    0.609    0.609
#    S7                0.657    0.002  400.506    0.000    0.657    0.657
#    S8                0.607    0.002  360.474    0.000    0.607    0.607
#  factor5 =~                                                            
#    E1                0.430    0.003  164.681    0.000    0.430    0.430
#    E2                0.487    0.002  197.619    0.000    0.487    0.487
#    E3                0.367    0.002  147.353    0.000    0.367    0.367
#    E4                0.385    0.003  144.998    0.000    0.385    0.385
#    E5                0.289    0.003  106.375    0.000    0.289    0.289
#    E6                0.653    0.002  277.755    0.000    0.653    0.653
#    E7                0.427    0.003  162.641    0.000    0.427    0.427
#    E8                0.249    0.003   97.012    0.000    0.249    0.249
#  factor6 =~                                                            
#    C1                0.657    0.002  329.224    0.000    0.657    0.657
#    C2                0.416    0.002  187.290    0.000    0.416    0.416
#    C3                0.612    0.002  316.869    0.000    0.612    0.612
#    C4                0.600    0.002  308.526    0.000    0.600    0.600
#    C5                0.549    0.002  249.993    0.000    0.549    0.549
#    C6                0.430    0.002  187.192    0.000    0.430    0.430
#    C7                0.547    0.002  266.858    0.000    0.547    0.547
#    C8                0.567    0.002  285.581    0.000    0.567    0.567
#  global =~                                                             
#    R1                0.525    0.002  225.152    0.000    0.525    0.525
#    R2                0.422    0.003  156.699    0.000    0.422    0.422
#    R3                0.439    0.003  154.546    0.000    0.439    0.439
#    R4                0.528    0.002  227.354    0.000    0.528    0.528
#    R5                0.555    0.003  219.729    0.000    0.555    0.555
#    R6                0.493    0.002  201.472    0.000    0.493    0.493
#    R7                0.635    0.002  300.674    0.000    0.635    0.635
#    R8                0.494    0.003  194.253    0.000    0.494    0.494
#    I1                0.192    0.003   67.751    0.000    0.192    0.192
#    I2                0.240    0.003   88.804    0.000    0.240    0.240
#    I3                0.259    0.003   97.583    0.000    0.259    0.259
#    I4                0.313    0.003  121.250    0.000    0.313    0.313
#    I5                0.264    0.003  100.198    0.000    0.264    0.264
#    I6                0.310    0.003  122.340    0.000    0.310    0.310
#    I7                0.281    0.003  106.631    0.000    0.281    0.281
#    I8                0.433    0.002  183.975    0.000    0.433    0.433
#    A1                0.323    0.003  118.904    0.000    0.323    0.323
#    A2                0.342    0.003  132.195    0.000    0.342    0.342
#    A3                0.314    0.003  121.597    0.000    0.314    0.314
#    A4                0.258    0.003   95.954    0.000    0.258    0.258
#    A5                0.220    0.003   81.165    0.000    0.220    0.220
#    A6                0.313    0.003  119.123    0.000    0.313    0.313
#    A7                0.413    0.002  169.805    0.000    0.413    0.413
#    A8                0.416    0.002  174.630    0.000    0.416    0.416
#    S1                0.331    0.003  125.865    0.000    0.331    0.331
#    S2                0.281    0.003  104.090    0.000    0.281    0.281
#    S3                0.232    0.003   85.735    0.000    0.232    0.232
#    S4                0.441    0.002  187.667    0.000    0.441    0.441
#    S5                0.227    0.003   83.829    0.000    0.227    0.227
#    S6                0.397    0.002  164.478    0.000    0.397    0.397
#    S7                0.374    0.002  151.402    0.000    0.374    0.374
#    S8                0.386    0.002  156.455    0.000    0.386    0.386
#    E1                0.492    0.002  201.043    0.000    0.492    0.492
#    E2                0.479    0.002  199.317    0.000    0.479    0.479
#    E3                0.613    0.002  306.281    0.000    0.613    0.613
#    E4                0.454    0.003  181.258    0.000    0.454    0.454
#    E5                0.590    0.002  285.096    0.000    0.590    0.590
#    E6                0.522    0.002  228.928    0.000    0.522    0.522
#    E7                0.477    0.002  200.322    0.000    0.477    0.477
#    E8                0.627    0.002  310.423    0.000    0.627    0.627
#    C1                0.473    0.002  195.353    0.000    0.473    0.473
#    C2                0.589    0.002  286.408    0.000    0.589    0.589
#    C3                0.575    0.002  269.474    0.000    0.575    0.575
#    C4                0.532    0.002  240.263    0.000    0.532    0.532
#    C5                0.473    0.002  198.305    0.000    0.473    0.473
#    C6                0.537    0.002  242.210    0.000    0.537    0.537
#    C7                0.572    0.002  266.542    0.000    0.572    0.572
#    C8                0.630    0.002  316.372    0.000    0.630    0.630

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .55

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.global           PUC  Omega.global OmegaH.global 
#0.3541219     0.8510638     0.9511526     0.7066011 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5398665 0.09442682 0.4601335 0.8846658 0.4510117 0.8002384 0.8805007
#factor2 0.8587923 0.15951676 0.1412077 0.8967934 0.7630003 0.9154795 0.9552710
#factor3 0.7932316 0.12743494 0.2067684 0.8682402 0.6805386 0.8666988 0.9271886
#factor4 0.7666169 0.11645881 0.2333831 0.8451473 0.6285897 0.8646733 0.9270257
#factor5 0.3886993 0.05532853 0.6113007 0.8375343 0.3118422 0.6656586 0.8006244
#factor6 0.5026728 0.09271224 0.4973272 0.8994528 0.4326488 0.7867052 0.8653934
#global  0.3541219 0.35412191 0.3541219 0.9511526 0.7066011 0.9271746 0.9091705






####################################
####### HEXACO Anglim et al. https://osf.io/uwdgs/?view_only=    https://osf.io/5wpj6
### 
### not included in analysis because AES scale used

osf_5wpj6_HEXACO <- read.csv("~/Self regulation/osf_5wpj6_HEXACO.csv")
colnames(osf_5wpj6_HEXACO)
mydata <- as.data.frame(osf_5wpj6_HEXACO[,13:212])
mydata <- na.omit(mydata)
min(mydata)
max(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 20 components
# Eigenvalue1 = 19.33, eigenvalue2 = 9.42

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 20 components 
# Eigenvalue 1 = 23.74; eigenvalue 2 = 10.89


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.01, RMSEA=.054, RMSR=.09, TLI=.219

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.12, RMSEA=.07, RMSR=.11, TLI=.178


# Give solution with 6 factors
fit5 <- fa(mydata,6)
fit5
diagram(fit5)
# %variance explained=.28, RMSEA=.038, RMSR=.04, TLI=.621
#      MR2   MR6   MR1   MR5   MR4   MR3
#MR2  1.00  0.06  0.09 -0.12 -0.08 -0.04
#MR6  0.06  1.00  0.12 -0.13 -0.14 -0.07
#MR1  0.09  0.12  1.00 -0.19 -0.17 -0.20
#MR5 -0.12 -0.13 -0.19  1.00  0.15  0.12
#MR4 -0.08 -0.14 -0.17  0.15  1.00  0.07
#MR3 -0.04 -0.07 -0.20  0.12  0.07  1.00

fit6 <- fa(rho,6,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.34, RMSEA=.054, RMSR=.05, TLI=.5
#      MR2   MR3   MR1   MR4   MR6   MR5
#MR2  1.00  0.05  0.11 -0.10 -0.13 -0.05
#MR3  0.05  1.00  0.11 -0.12 -0.13 -0.06
#MR1  0.11  0.11  1.00 -0.19 -0.20 -0.22
#MR4 -0.10 -0.12 -0.19  1.00  0.17  0.10
#MR6 -0.13 -0.13 -0.20  0.17  1.00  0.13
#MR5 -0.05 -0.06 -0.22  0.10  0.13  1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~ hexaco1+hexaco2+hexaco3+hexaco4+hexaco5+hexaco6+hexaco7+hexaco8+hexaco9+hexaco10+
           hexaco11+hexaco12+hexaco13+hexaco14+hexaco15+hexaco16+hexaco17+hexaco18+hexaco19+hexaco20+
           hexaco21+hexaco22+hexaco23+hexaco24+hexaco25+hexaco26+hexaco27+hexaco28+hexaco29+hexaco30+
           hexaco31+hexaco32+hexaco33+hexaco34+hexaco35+hexaco36+hexaco37+hexaco38+hexaco39+hexaco40+
           hexaco41+hexaco42+hexaco43+hexaco44+hexaco45+hexaco46+hexaco47+hexaco48+hexaco49+hexaco50+
           hexaco51+hexaco52+hexaco53+hexaco54+hexaco55+hexaco56+hexaco57+hexaco58+hexaco59+hexaco60+
           hexaco61+hexaco62+hexaco63+hexaco64+hexaco65+hexaco66+hexaco67+hexaco68+hexaco69+hexaco70+
           hexaco71+hexaco72+hexaco73+hexaco74+hexaco75+hexaco76+hexaco77+hexaco78+hexaco79+hexaco80+
           hexaco81+hexaco82+hexaco83+hexaco84+hexaco85+hexaco86+hexaco87+hexaco88+hexaco89+hexaco90+
           hexaco91+hexaco92+hexaco93+hexaco94+hexaco95+hexaco96+hexaco97+hexaco98+hexaco99+hexaco100+
           hexaco101+hexaco102+hexaco103+hexaco104+hexaco105+hexaco106+hexaco107+hexaco108+hexaco109+hexaco110+
           hexaco111+hexaco112+hexaco113+hexaco114+hexaco115+hexaco116+hexaco117+hexaco118+hexaco119+hexaco120+
           hexaco121+hexaco122+hexaco123+hexaco124+hexaco125+hexaco126+hexaco127+hexaco128+hexaco129+hexaco130+
           hexaco131+hexaco132+hexaco133+hexaco134+hexaco135+hexaco136+hexaco137+hexaco138+hexaco139+hexaco140+
           hexaco141+hexaco142+hexaco143+hexaco144+hexaco145+hexaco146+hexaco147+hexaco148+hexaco149+hexaco150+
           hexaco151+hexaco152+hexaco153+hexaco154+hexaco155+hexaco156+hexaco157+hexaco158+hexaco159+hexaco160+
           hexaco161+hexaco162+hexaco163+hexaco164+hexaco165+hexaco166+hexaco167+hexaco168+hexaco169+hexaco170+
           hexaco171+hexaco172+hexaco173+hexaco174+hexaco175+hexaco176+hexaco177+hexaco178+hexaco179+hexaco180+
           hexaco181+hexaco182+hexaco183+hexaco184+hexaco185+hexaco186+hexaco187+hexaco188+hexaco189+hexaco190+
           hexaco191+hexaco192+hexaco193+hexaco194+hexaco195+hexaco196+hexaco197+hexaco198+hexaco199+hexaco200
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
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

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)



EGAmodel= '
 factor1=~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16
 factor2=~ C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16
 factor3=~ E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16
 factor4=~ H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16
 factor5=~ O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16
 factor6=~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
'

library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.758       0.682
#Tucker-Lewis Index (TLI)                       0.752       0.674
#Robust Comparative Fit Index (CFI)                         0.454
#Robust Tucker-Lewis Index (TLI)                            0.440
#RMSEA                                          0.086       0.044
#Robust RMSEA                                               0.079
#SRMR                                           0.092       0.092

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.022    0.012    1.936    0.053    0.095    0.095
#  factor3          -0.079    0.016   -4.819    0.000   -0.286   -0.286
#  factor4           0.095    0.017    5.692    0.000    0.273    0.273
#  factor5           0.036    0.025    1.445    0.149    0.081    0.081
#  factor6           0.118    0.023    5.209    0.000    0.280    0.280
#factor2 ~~                                                            
#  factor3           0.006    0.009    0.652    0.514    0.033    0.033
#  factor4           0.057    0.013    4.411    0.000    0.265    0.265
#  factor5           0.035    0.014    2.598    0.009    0.127    0.127
#  factor6           0.079    0.017    4.630    0.000    0.299    0.299
#factor3 ~~                                                            
#  factor4           0.055    0.013    4.268    0.000    0.218    0.218
#  factor5           0.021    0.016    1.271    0.204    0.064    0.064
#  factor6          -0.104    0.017   -6.128    0.000   -0.337   -0.337
#factor4 ~~                                                            
#  factor5           0.070    0.019    3.700    0.000    0.171    0.171
#  factor6          -0.019    0.019   -1.012    0.311   -0.050   -0.050
#factor5 ~~                                                            
#  factor6           0.084    0.025    3.337    0.001    0.169    0.169


fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .30

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.649       0.621
#Tucker-Lewis Index (TLI)                       0.641       0.612
#Robust Comparative Fit Index (CFI)                         0.444
#Robust Tucker-Lewis Index (TLI)                            0.432
#RMSEA                                          0.103       0.048
#Robust RMSEA                                               0.079
#SRMR                                           0.107       0.107


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .31


#bifactor model 
BIFmodel= '
 factor1=~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16
 factor2=~ C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16
 factor3=~ E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16
 factor4=~ H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16
 factor5=~ O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16
 factor6=~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
 global =~ A1+A2+A3+A4+A5+A6+A7+A8+A9+A10+A11+A12+A13+A14+A15+A16+
           C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16+
           E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+E16+
           H1+H2+H3+H4+H5+H6+H7+H8+H9+H10+H11+H12+H13+H14+H15+H16+
           O1+O2+O3+O4+O5+O6+O7+O8+O9+O10+O11+O12+O13+O14+O15+O16+
           X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16
'
CFA_model3 <- cfa(BIFmodel, data = mydata,ordered=TRUE,orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.790       0.713
#Tucker-Lewis Index (TLI)                       0.781       0.700
#Robust Comparative Fit Index (CFI)                         0.508
#Robust Tucker-Lewis Index (TLI)                            0.486
#RMSEA                                          0.080       0.042
#Robust RMSEA                                               0.075
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .33

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#ECV.global           PUC  Omega.global OmegaH.global 
#0.2443841     0.8421053     0.8766862     0.2969120 

#$FactorLevelIndices
#         ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8674604 0.13210199 0.13253957 0.8519069 0.7555523 0.8599401 0.9272945
#factor2 0.8765483 0.14506614 0.12345172 0.8548130 0.7990429 0.8757369 0.9372008
#factor3 0.7780358 0.12619331 0.22196420 0.8393431 0.7898870 0.8528415 0.9274576
#factor4 0.9447302 0.16608208 0.05526982 0.8672107 0.8650831 0.8946445 0.9476152
#factor5 0.9534080 0.12419576 0.04659202 0.8089543 0.7850749 0.8602615 0.9277829
#factor6 0.2896663 0.06197666 0.71033375 0.8829112 0.1229757 0.7334028 0.8823083
#global  0.2443841 0.24438407 0.24438407 0.8766862 0.2969120 0.9229271 0.9588863


#Reverse code items factor three to bring them in line with the rest
mydata2 <- mydata
reverse_cols = c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16")
mydata2[ , reverse_cols] = 6 - mydata2[ , reverse_cols]
CFA_model3 <- cfa(BIFmodel, data = mydata2,ordered=TRUE,orthogonal=TRUE,std.lv=TRUE)
options(max.print=999999)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.790       0.713
#Tucker-Lewis Index (TLI)                       0.781       0.700
#Robust Comparative Fit Index (CFI)                         0.508
#Robust Tucker-Lewis Index (TLI)                            0.486
#RMSEA                                          0.080       0.042
#Robust RMSEA                                               0.075
#SRMR                                           0.088       0.088

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    A1                0.547    0.037   14.667    0.000    0.547    0.547
#    A2                0.541    0.034   15.876    0.000    0.541    0.541
#    A3                0.401    0.040   10.011    0.000    0.401    0.401
#    A4                0.623    0.033   18.908    0.000    0.623    0.623
#    A5                0.653    0.030   21.880    0.000    0.653    0.653
#    A6                0.493    0.034   14.328    0.000    0.493    0.493
#    A7                0.364    0.042    8.570    0.000    0.364    0.364
#    A8                0.539    0.035   15.370    0.000    0.539    0.539
#    A9                0.310    0.047    6.535    0.000    0.310    0.310
#    A10               0.429    0.041   10.464    0.000    0.429    0.429
#    A11               0.443    0.039   11.414    0.000    0.443    0.443
#    A12               0.602    0.031   19.541    0.000    0.602    0.602
#    A13               0.549    0.036   15.125    0.000    0.549    0.549
#    A14               0.525    0.037   14.243    0.000    0.525    0.525
#    A15               0.430    0.040   10.883    0.000    0.430    0.430
#    A16               0.631    0.031   20.119    0.000    0.631    0.631
#  factor2 =~                                                            
#    C1                0.325    0.043    7.592    0.000    0.325    0.325
#    C2                0.526    0.039   13.502    0.000    0.526    0.526
#    C3                0.555    0.037   15.073    0.000    0.555    0.555
#    C4                0.544    0.035   15.669    0.000    0.544    0.544
#    C5                0.589    0.035   17.042    0.000    0.589    0.589
#    C6                0.631    0.032   19.542    0.000    0.631    0.631
#    C7                0.662    0.033   20.024    0.000    0.662    0.662
#    C8                0.614    0.033   18.749    0.000    0.614    0.614
#    C9                0.388    0.045    8.561    0.000    0.388    0.388
#    C10               0.546    0.034   16.004    0.000    0.546    0.546
#    C11               0.677    0.031   22.022    0.000    0.677    0.677
#    C12               0.329    0.040    8.239    0.000    0.329    0.329
#    C13               0.510    0.040   12.668    0.000    0.510    0.510
#    C14               0.609    0.033   18.369    0.000    0.609    0.609
#    C15               0.499    0.038   13.202    0.000    0.499    0.499
#    C16               0.455    0.040   11.304    0.000    0.455    0.455
#  factor3 =~                                                            
#    E1                0.442    0.040   10.974    0.000    0.442    0.442
#    E2                0.420    0.040   10.545    0.000    0.420    0.420
#    E3                0.610    0.032   18.907    0.000    0.610    0.610
#    E4                0.639    0.031   20.611    0.000    0.639    0.639
#    E5                0.355    0.043    8.222    0.000    0.355    0.355
#    E6                0.431    0.039   11.085    0.000    0.431    0.431
#    E7                0.489    0.039   12.567    0.000    0.489    0.489
#    E8                0.535    0.038   13.918    0.000    0.535    0.535
#    E9                0.440    0.040   11.099    0.000    0.440    0.440
#    E10               0.273    0.046    5.932    0.000    0.273    0.273
#    E11               0.504    0.037   13.494    0.000    0.504    0.504
#    E12               0.591    0.036   16.547    0.000    0.591    0.591
#    E13               0.571    0.033   17.267    0.000    0.571    0.571
#    E14               0.515    0.038   13.579    0.000    0.515    0.515
#    E15               0.419    0.042    9.910    0.000    0.419    0.419
#    E16               0.644    0.032   19.973    0.000    0.644    0.644
#  factor4 =~                                                            
#    H1                0.583    0.032   18.417    0.000    0.583    0.583
#    H2                0.651    0.031   21.087    0.000    0.651    0.651
#    H3                0.463    0.035   13.120    0.000    0.463    0.463
#    H4                0.453    0.037   12.189    0.000    0.453    0.453
#    H5                0.457    0.035   12.970    0.000    0.457    0.457
#    H6                0.554    0.035   15.641    0.000    0.554    0.554
#    H7                0.663    0.028   24.062    0.000    0.663    0.663
#    H8                0.458    0.038   12.015    0.000    0.458    0.458
#    H9                0.611    0.030   20.084    0.000    0.611    0.611
#    H10               0.603    0.030   19.864    0.000    0.603    0.603
#    H11               0.691    0.028   24.879    0.000    0.691    0.691
#    H12               0.504    0.034   14.979    0.000    0.504    0.504
#    H13               0.469    0.034   13.766    0.000    0.469    0.469
#    H14               0.665    0.030   22.384    0.000    0.665    0.665
#    H15               0.638    0.031   20.538    0.000    0.638    0.638
#    H16               0.659    0.035   18.655    0.000    0.659    0.659
#  factor5 =~                                                            
#    O1                0.731    0.029   25.144    0.000    0.731    0.731
#    O2                0.460    0.044   10.530    0.000    0.460    0.460
#    O3                0.429    0.042   10.104    0.000    0.429    0.429
#    O4                0.361    0.048    7.574    0.000    0.361    0.361
#    O5                0.521    0.041   12.812    0.000    0.521    0.521
#    O6                0.229    0.050    4.570    0.000    0.229    0.229
#    O7                0.689    0.031   21.995    0.000    0.689    0.689
#    O8                0.338    0.046    7.421    0.000    0.338    0.338
#    O9                0.602    0.035   17.398    0.000    0.602    0.602
#    O10               0.401    0.046    8.781    0.000    0.401    0.401
#    O11               0.466    0.041   11.487    0.000    0.466    0.466
#    O12               0.276    0.047    5.848    0.000    0.276    0.276
#    O13               0.395    0.045    8.856    0.000    0.395    0.395
#    O14               0.530    0.039   13.542    0.000    0.530    0.530
#    O15               0.605    0.036   16.773    0.000    0.605    0.605
#    O16               0.620    0.035   17.481    0.000    0.620    0.620
#  factor6 =~                                                            
#    X1               -0.159    0.057   -2.793    0.005   -0.159   -0.159
#    X2                0.360    0.045    7.946    0.000    0.360    0.360
#    X3                0.387    0.049    7.970    0.000    0.387    0.387
#    X4                0.105    0.058    1.819    0.069    0.105    0.105
#    X5                0.128    0.056    2.298    0.022    0.128    0.128
#    X6                0.624    0.045   13.970    0.000    0.624    0.624
#    X7                0.560    0.045   12.325    0.000    0.560    0.560
#    X8               -0.097    0.056   -1.747    0.081   -0.097   -0.097
#    X9                0.231    0.053    4.339    0.000    0.231    0.231
#    X10               0.564    0.048   11.798    0.000    0.564    0.564
#    X11               0.400    0.048    8.328    0.000    0.400    0.400
#    X12              -0.153    0.049   -3.110    0.002   -0.153   -0.153
#    X13              -0.152    0.059   -2.570    0.010   -0.152   -0.152
#    X14               0.184    0.048    3.858    0.000    0.184    0.184
#    X15               0.559    0.042   13.264    0.000    0.559    0.559
#    X16               0.133    0.056    2.383    0.017    0.133    0.133
#  global =~                                                             
#    A1                0.294    0.046    6.321    0.000    0.294    0.294
#    A2                0.084    0.049    1.722    0.085    0.084    0.084
#    A3                0.213    0.046    4.591    0.000    0.213    0.213
#    A4                0.058    0.050    1.166    0.244    0.058    0.058
#    A5                0.265    0.043    6.137    0.000    0.265    0.265
#    A6                0.220    0.045    4.910    0.000    0.220    0.220
#    A7                0.237    0.049    4.895    0.000    0.237    0.237
#    A8                0.243    0.045    5.451    0.000    0.243    0.243
#    A9                0.213    0.053    4.049    0.000    0.213    0.213
#    A10               0.099    0.048    2.051    0.040    0.099    0.099
#    A11               0.183    0.050    3.649    0.000    0.183    0.183
#    A12               0.178    0.044    4.024    0.000    0.178    0.178
#    A13               0.154    0.046    3.337    0.001    0.154    0.154
#    A14               0.048    0.050    0.950    0.342    0.048    0.048
#    A15               0.326    0.046    7.006    0.000    0.326    0.326
#    A16               0.136    0.047    2.876    0.004    0.136    0.136
#    C1                0.208    0.048    4.333    0.000    0.208    0.208
#    C2                0.227    0.050    4.573    0.000    0.227    0.227
#    C3               -0.096    0.050   -1.942    0.052   -0.096   -0.096
#    C4               -0.002    0.046   -0.039    0.969   -0.002   -0.002
#    C5                0.071    0.053    1.339    0.181    0.071    0.071
#    C6                0.272    0.050    5.465    0.000    0.272    0.272
#    C7                0.042    0.051    0.826    0.409    0.042    0.042
#    C8                0.180    0.049    3.720    0.000    0.180    0.180
#    C9                0.134    0.051    2.628    0.009    0.134    0.134
#    C10               0.444    0.042   10.530    0.000    0.444    0.444
#    C11               0.069    0.054    1.276    0.202    0.069    0.069
#    C12               0.242    0.048    5.090    0.000    0.242    0.242
#    C13               0.263    0.048    5.438    0.000    0.263    0.263
#    C14               0.294    0.050    5.842    0.000    0.294    0.294
#    C15               0.001    0.049    0.013    0.990    0.001    0.001
#    C16              -0.059    0.049   -1.212    0.226   -0.059   -0.059
#    E1                0.159    0.050    3.204    0.001    0.159    0.159
#    E2                0.480    0.041   11.666    0.000    0.480    0.480
#    E3               -0.079    0.050   -1.559    0.119   -0.079   -0.079
#    E4               -0.025    0.050   -0.494    0.622   -0.025   -0.025
#    E5                0.197    0.048    4.094    0.000    0.197    0.197
#    E6                0.445    0.040   11.215    0.000    0.445    0.445
#    E7                0.235    0.047    5.038    0.000    0.235    0.235
#    E8               -0.133    0.050   -2.637    0.008   -0.133   -0.133
#    E9                0.278    0.047    5.891    0.000    0.278    0.278
#    E10               0.312    0.046    6.801    0.000    0.312    0.312
#    E11              -0.203    0.048   -4.246    0.000   -0.203   -0.203
#    E12              -0.079    0.050   -1.592    0.111   -0.079   -0.079
#    E13               0.382    0.046    8.285    0.000    0.382    0.382
#    E14               0.389    0.045    8.640    0.000    0.389    0.389
#    E15              -0.248    0.049   -5.104    0.000   -0.248   -0.248
#    E16              -0.022    0.052   -0.432    0.666   -0.022   -0.022
#    H1                0.076    0.048    1.567    0.117    0.076    0.076
#    H2                0.228    0.051    4.500    0.000    0.228    0.228
#    H3                0.178    0.049    3.632    0.000    0.178    0.178
#    H4               -0.150    0.054   -2.809    0.005   -0.150   -0.150
#    H5                0.057    0.051    1.118    0.263    0.057    0.057
#    H6                0.142    0.056    2.549    0.011    0.142    0.142
#    H7                0.036    0.053    0.667    0.505    0.036    0.036
#    H8               -0.018    0.050   -0.360    0.719   -0.018   -0.018
#    H9                0.059    0.051    1.152    0.249    0.059    0.059
#    H10               0.168    0.053    3.189    0.001    0.168    0.168
#    H11              -0.037    0.050   -0.739    0.460   -0.037   -0.037
#    H12              -0.174    0.048   -3.628    0.000   -0.174   -0.174
#    H13               0.020    0.051    0.387    0.698    0.020    0.020
#    H14               0.198    0.053    3.697    0.000    0.198    0.198
#    H15              -0.011    0.050   -0.228    0.819   -0.011   -0.011
#    H16              -0.267    0.047   -5.651    0.000   -0.267   -0.267
#    O1                0.045    0.052    0.869    0.385    0.045    0.045
#    O2                0.187    0.051    3.635    0.000    0.187    0.187
#    O3                0.204    0.048    4.239    0.000    0.204    0.204
#    O4               -0.026    0.050   -0.521    0.603   -0.026   -0.026
#    O5                0.081    0.050    1.596    0.110    0.081    0.081
#    O6                0.136    0.051    2.656    0.008    0.136    0.136
#    O7               -0.034    0.052   -0.660    0.509   -0.034   -0.034
#    O8                0.030    0.050    0.595    0.552    0.030    0.030
#    O9                0.127    0.053    2.424    0.015    0.127    0.127
#    O10               0.126    0.050    2.542    0.011    0.126    0.126
#    O11               0.126    0.051    2.456    0.014    0.126    0.126
#    O12              -0.008    0.052   -0.157    0.875   -0.008   -0.008
#    O13               0.014    0.054    0.252    0.801    0.014    0.014
#    O14               0.118    0.052    2.292    0.022    0.118    0.118
#    O15               0.146    0.049    2.978    0.003    0.146    0.146
#    O16               0.070    0.050    1.389    0.165    0.070    0.070
#    X1                0.716    0.036   20.158    0.000    0.716    0.716
#    X2                0.344    0.048    7.166    0.000    0.344    0.344
#    X3                0.432    0.046    9.372    0.000    0.432    0.432
#    X4                0.696    0.030   23.008    0.000    0.696    0.696
#    X5                0.337    0.048    6.999    0.000    0.337    0.337
#    X6                0.548    0.047   11.741    0.000    0.548    0.548
#    X7                0.389    0.049    7.919    0.000    0.389    0.389
#    X8                0.797    0.023   34.795    0.000    0.797    0.797
#    X9                0.550    0.038   14.288    0.000    0.550    0.550
#    X10               0.440    0.048    9.184    0.000    0.440    0.440
#    X11               0.502    0.043   11.648    0.000    0.502    0.502
#    X12               0.563    0.037   15.189    0.000    0.563    0.563
#    X13               0.748    0.029   26.083    0.000    0.748    0.748
#    X14               0.376    0.046    8.187    0.000    0.376    0.376
#    X15               0.436    0.045    9.617    0.000    0.436    0.436
#    X16               0.641    0.031   20.620    0.000    0.641    0.641


bifactorIndices(CFA_model3)
$ModelLevelIndices
ECV.global           PUC  Omega.global OmegaH.global 
0.2443841     0.8421053     0.8978254     0.4174395 

$FactorLevelIndices
#         ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8674604 0.13210199 0.13253957 0.8519069 0.7555523 0.8599401 0.9272945
#factor2 0.8765483 0.14506614 0.12345172 0.8548130 0.7990429 0.8757369 0.9372008
#factor3 0.7780358 0.12619331 0.22196420 0.8393431 0.7898870 0.8528416 0.9274576
#factor4 0.9447302 0.16608208 0.05526982 0.8672107 0.8650831 0.8946445 0.9476152
#factor5 0.9534080 0.12419576 0.04659202 0.8089543 0.7850749 0.8602615 0.9277829
#factor6 0.2896663 0.06197666 0.71033375 0.8829112 0.1229757 0.7334028 0.8823083
#global  0.2443841 0.24438407 0.24438407 0.8978254 0.4174395 0.9229271 0.9588863

max.print
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")







####### TIPI
### Coelho et al. (2018) https://pdf.sciencedirectassets.com/271782/1-s2.0-S0191886918X00099/1-s2.0-S0191886918303623/main.pdf
### UK sample https://osf.io/9kjxs/?view_only=d532eaf8b5234c8fb1b5eeeb27f0190f

library(haven)
TIPI_Coelho_UK <- read_sav("TIPI_Coelho_UK.sav")
colnames(TIPI_Coelho_UK)
mydata <- as.data.frame(TIPI_Coelho_UK[,1:10])
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factor and 3 components
# Eigenvalue1 = 1.98, eigenvalue2 = .98

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components 
# Eigenvalue 1 = 2.11; eigenvalue 2 = 1.07


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.20, RMSEA=.179, RMSR=.14, TLI=.326

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.21, RMSEA=.211, RMSR=.16, TLI=.265


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.58, RMSEA=.07, RMSR=.02, TLI=.896
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.14  0.22 -0.35 -0.47
#MR2 -0.14  1.00 -0.29 -0.03 -0.03
#MR4  0.22 -0.29  1.00  0.05 -0.11
#MR3 -0.35 -0.03  0.05  1.00  0.18
#MR5 -0.47 -0.03 -0.11  0.18  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.63, RMSEA=.101, RMSR=.02, TLI=.831
#      MR1   MR2   MR4   MR5   MR3
#MR1  1.00 -0.19  0.22  0.41 -0.35
#MR2 -0.19  1.00 -0.28  0.06 -0.06
#MR4  0.22 -0.28  1.00  0.17  0.06
#MR5  0.41  0.06  0.17  1.00 -0.16
#MR3 -0.35 -0.06  0.06 -0.16  1.00

colnames(mydata)

#Analysis lavaan
library(lavaan)

UNImodel= '
 global =~  TIPI01+TIPI02+TIPI03+TIPI04+TIPI05+TIPI06+TIPI07+TIPI08+TIPI09+TIPI10
'
CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.691       0.563
#Tucker-Lewis Index (TLI)                       0.603       0.438
#Robust Comparative Fit Index (CFI)                         0.384
#Robust Tucker-Lewis Index (TLI)                            0.208
#RMSEA                                          0.233       0.228
#Robust RMSEA                                               0.222
#SRMR                                           0.149       0.149

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .208

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.502       0.470
#Tucker-Lewis Index (TLI)                       0.359       0.318
#Robust Comparative Fit Index (CFI)                         0.502
#Robust Tucker-Lewis Index (TLI)                            0.359
#RMSEA                                          0.177       0.165
#Robust RMSEA                                               0.176
#SRMR                                           0.131       0.131

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .117

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


EGAmodel= '
 factor1=~ TIPI01+TIPI06
 factor2=~ TIPI02+TIPI07
 factor3=~ TIPI03+TIPI08
 factor4=~ TIPI04+TIPI09
 factor5=~ TIPI05+TIPI10
'
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.966       0.931
#Tucker-Lewis Index (TLI)                       0.939       0.876
#Robust Comparative Fit Index (CFI)                         0.860
#Robust Tucker-Lewis Index (TLI)                            0.747
#RMSEA                                          0.092       0.107
#Robust RMSEA                                               0.126
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .501

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.000    0.065   -0.002    0.999   -0.000   -0.000
#    factor3           0.087    0.056    1.561    0.118    0.087    0.087
#    factor4          -0.240    0.052   -4.605    0.000   -0.240   -0.240
#    factor5           0.351    0.075    4.688    0.000    0.351    0.351
#  factor2 ~~                                                            
#    factor3          -0.308    0.091   -3.383    0.001   -0.308   -0.308
#    factor4           0.584    0.101    5.806    0.000    0.584    0.584
#    factor5          -0.036    0.060   -0.596    0.551   -0.036   -0.036
#  factor3 ~~                                                            
#    factor4          -0.557    0.052  -10.767    0.000   -0.557   -0.557
#    factor5           0.059    0.055    1.067    0.286    0.059    0.059
#  factor4 ~~                                                            
#    factor5          -0.189    0.054   -3.530    0.000   -0.189   -0.189

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning 
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.908       0.896
#Tucker-Lewis Index (TLI)                       0.835       0.812
#Robust Comparative Fit Index (CFI)                         0.909
#Robust Tucker-Lewis Index (TLI)                            0.836
#RMSEA                                          0.090       0.087
#Robust RMSEA                                               0.089
#SRMR                                           0.064       0.064

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.000    0.065   -0.002    0.999   -0.000   -0.000
#    factor3           0.087    0.056    1.561    0.118    0.087    0.087
#    factor4          -0.240    0.052   -4.605    0.000   -0.240   -0.240
#    factor5           0.351    0.075    4.688    0.000    0.351    0.351
#  factor2 ~~                                                            
#    factor3          -0.308    0.091   -3.383    0.001   -0.308   -0.308
#    factor4           0.584    0.101    5.806    0.000    0.584    0.584
#    factor5          -0.036    0.060   -0.596    0.551   -0.036   -0.036
#  factor3 ~~                                                            
#    factor4          -0.557    0.052  -10.767    0.000   -0.557   -0.557
#    factor5           0.059    0.055    1.067    0.286    0.059    0.059
#  factor4 ~~                                                            
#    factor5          -0.189    0.054   -3.530    0.000   -0.189   -0.189

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .473

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
# also warning message
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.610       0.660
#Tucker-Lewis Index (TLI)                       0.499       0.563
#Robust Comparative Fit Index (CFI)                         0.622
#Robust Tucker-Lewis Index (TLI)                            0.514
#RMSEA                                          0.156       0.132
#Robust RMSEA                                               0.153
#SRMR                                           0.157       0.157

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .544


#bifactor model 
BIFmodel= '
 factor1=~ TIPI01+TIPI06
 factor2=~ TIPI02+TIPI07
 factor3=~ TIPI03+TIPI08
 factor4=~ TIPI04+TIPI09
 factor5=~ TIPI05+TIPI10
 global =~  TIPI01+TIPI02+TIPI03+TIPI04+TIPI05+TIPI06+TIPI07+TIPI08+TIPI09+TIPI10
'

CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.784       0.724
#Tucker-Lewis Index (TLI)                       0.611       0.503
#Robust Comparative Fit Index (CFI)                         0.779
#Robust Tucker-Lewis Index (TLI)                            0.603
#RMSEA                                          0.138       0.141
#Robust RMSEA                                               0.138
#SRMR                                           0.102       0.102

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    TIPI01            1.011    0.163    6.210    0.000    1.011    0.584
#    TIPI06           -1.196    0.441   -2.713    0.007   -1.196   -0.720
#  factor2 =~                                                            
#    TIPI02            1.271    0.475    2.673    0.008    1.271    0.712
#    TIPI07           -0.280    0.128   -2.181    0.029   -0.280   -0.218
#  factor3 =~                                                            
#    TIPI03            0.884    0.081   10.945    0.000    0.884    0.705
#    TIPI08           -0.920    0.241   -3.825    0.000   -0.920   -0.579
#  factor4 =~                                                            
#    TIPI04            1.433    0.350    4.091    0.000    1.433    0.751
#    TIPI09           -1.062    0.190   -5.602    0.000   -1.062   -0.689
#  factor5 =~                                                            
#    TIPI05            0.671    0.485    1.383    0.167    0.671    0.497
#    TIPI10           -0.712    0.861   -0.828    0.408   -0.712   -0.410
#  global =~                                                             
#    TIPI01            1.148    0.579    1.983    0.047    1.148    0.662
#    TIPI02           -0.127    0.653   -0.195    0.846   -0.127   -0.071
#    TIPI03            0.396    0.294    1.344    0.179    0.396    0.316
#    TIPI04           -0.575    0.489   -1.176    0.239   -0.575   -0.301
#    TIPI05            0.757    0.197    3.836    0.000    0.757    0.561
#    TIPI06           -0.413    0.386   -1.071    0.284   -0.413   -0.249
#    TIPI07            0.438    0.138    3.160    0.002    0.438    0.341
#    TIPI08           -0.323    0.516   -0.626    0.531   -0.323   -0.203
#    TIPI09            0.651    0.600    1.086    0.277    0.651    0.423
#    TIPI10           -0.439    0.233   -1.889    0.059   -0.439   -0.253

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .571

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.global           PUC  Omega.global OmegaH.global 
#0.2760743     0.8888889     0.2682120     0.2248400 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS      Omega      OmegaH         H        FD
#factor1 0.6316312 0.16796234 0.3683688 0.22824563 0.022273852 0.6139310 0.7854048
#factor2 0.8205619 0.10856157 0.1794381 0.19300159 0.148672007 0.5189982 0.7211370
#factor3 0.8551518 0.16292512 0.1448482 0.02709341 0.015081831 0.5989857 0.7834612
#factor4 0.7940608 0.20320967 0.2059392 0.02614126 0.005368293 0.6870994 0.8369346
#factor5 0.5231342 0.08126697 0.4768658 0.07849000 0.005865534 0.3466191 0.6105094
#global  0.2760743 0.27607434 0.2760743 0.26821203 0.224839952 0.6647393 0.8028605






####### TIPI
### Coelho et al. (2018) https://pdf.sciencedirectassets.com/271782/1-s2.0-S0191886918X00099/1-s2.0-S0191886918303623/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEH8aCXVzLWVhc3QtMSJHMEUCIQCDS4QpaofSiVW%2BFHG8qHNzWK%2FU3zqoXiqNltIwA8OVNAIgPB6oxQLnptnMSC%2BA3Y8gBiK5cFLxCOMdV4Cp5LY2vEkquwUI2P%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDNxdNdxA9YyotQ5yKCqPBdHdZM%2F%2F%2Fj%2BotQx%2FbPfJs0bKKaRf2yXTOGdMR2u5ZI5NNR7kW%2BrNPGD5xH2JiXTWLQEd%2FCpUy43nVU7JIAl6tD7NcRtj2geMorr06%2FckH%2F7X3YEj08bllzF9%2BQFCG9M%2FMQTpHGlg076fbLZDdujjtJzwURl%2B99a9HVu7%2BaC%2Fdg4OFVbR9Aq8DZi4RRx8t2gySu7eXqtcCloXLZPUE%2FWBcb7jL5S6TT994h%2BQ5Jn5YVRzedD6v%2Fm%2B%2F%2Fipzkr21H%2FcNRFJv%2BkjMK5bTHDVkQXY1uhiOH5YThhTaV4N9wjiqCtyLa%2BhskhuFPj1xrIwYQvoLLny0JFS%2B3H%2FfwimXkRmpIoqa0Hiwh%2Bi%2BwHqgJudy7wkarIB1BKFdSkkBBSa410hD4TX3beUc2uJDh1bqx3Aj5G5G6MtFcpRct3bnkUEy%2FyI9xQhwEwdZGFnrHSDKqalCSPGRVYTGhiFdCSMZ1tdIJwLti6GMNEkqgCOQrb94qxhH%2Fc5wUCI7etYJ7NlQtXT4l0xbQiNML6C1a%2BClFbM6PX2E6rTe4xtZMVXnah4m0YAeNw10LWDUfazZJMmUNJdXsFPT10SZp4o6UjihTaMwGBlYV8V%2BgoDGPqb1GR5YB%2FsyNQJkYdrM%2BQU%2F5cHZPv2PsZeak%2FrAI4WeuEex15ov14609aKGlYh2avVIhsOm7%2FUkaJt9DJJmcg3byukIe99OcyN5Bs0lwzpRTx4jzDK38xpErsRGnpIeG20PIBQSBfBHTc3dnYbYxxUFvuI3BbGoE3qkCPQzKU1EQQne9sIOYYbwSOziWEmbcNUdSSQA3YCT%2BKDSReltHR3dQu54N3%2BsLVo4RcEfQmgEaM26RlEY6dkdcIJ7NKwNQGLLZsfr3Iwl8LRpAY6sQFr73KCLasp8o76%2FM6CmVmeuG2XeMT%2FK3olT61QtTMpOyPAhmwrWxkPLBU%2BDB11ajfpd4oDQE86k2US1xloXkKTYLnLvsNKjNwRa3WBnEsSwLbutuJJ8moqOjlOLse1eELnwvfLWVMMgn8Jze5jxdeJj6TswnGnDgJ%2B3IvE6Vt8MoWejGJTSAORHUzp71tCAzplO0hJ9gt472NOFQ1lrvgCddHmZJ9lBRFS7rfTxvklugo%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230622T152805Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQDYMSG3B%2F20230622%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=be147c12ca6fd5e1aa076585613a9e28e63f3cfb0ddacfa8e481a47d0fb7d5df&hash=92ae8e150cf24a3a7da7bbedac910e6dc8039a3cd4ea11a8a639f8e0c6defaeb&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0191886918303623&tid=spdf-6c91882a-2fe4-4725-83ce-42f6973f495a&sid=6551f5137aa8174ba29ac3f7d2fe0b003b38gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0703520b065759510504&rr=7db584dfcb512e56&cc=be
### Brasilian sample https://osf.io/9kjxs/?view_only=d532eaf8b5234c8fb1b5eeeb27f0190f

library(haven)
TIPI_Coelho_Bra <- read_sav("TIPI_Coelho_Bra.sav")
colnames(TIPI_Coelho_Bra)
mydata <- as.data.frame(TIPI_Coelho_Bra[,11:20])
colnames(mydata)
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue1 = 1.34, eigenvalue2 = .78

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components 
# Eigenvalue 1 = 1.53; eigenvalue 2 = .88


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.13, RMSEA=.142, RMSR=.12, TLI=.24

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.15, RMSEA=.17, RMSR=.14, TLI=.21


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.51, RMSEA=.074, RMSR=.02, TLI=.793
#      MR1   MR2   MR3   MR4   MR5
#MR1  1.00 -0.10  0.08  0.23 -0.11
#MR2 -0.10  1.00 -0.11 -0.07 -0.04
#MR3  0.08 -0.11  1.00  0.05  0.05
#MR4  0.23 -0.07  0.05  1.00  0.06
#MR5 -0.11 -0.04  0.05  0.06  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.56, RMSEA=.108, RMSR=.02, TLI=.68
#      MR1   MR2   MR3   MR4   MR5
#MR1  1.00 -0.07  0.05  0.25 -0.04
#MR2 -0.07  1.00 -0.12 -0.10  0.01
#MR3  0.05 -0.12  1.00  0.08 -0.09
#MR4  0.25 -0.10  0.08  1.00  0.10
#MR5 -0.04  0.01 -0.09  0.10  1.00

colnames(mydata)

#Analysis lavaan
library(lavaan)

UNImodel= '
 global =~  TIPI01+TIPI02+TIPI03+TIPI04+TIPI05+TIPI06+TIPI07+TIPI08+TIPI09+TIPI10
'
CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.573       0.430
#Tucker-Lewis Index (TLI)                       0.452       0.268
#Robust Comparative Fit Index (CFI)                         0.377
#Robust Tucker-Lewis Index (TLI)                            0.199
#RMSEA                                          0.166       0.168
#Robust RMSEA                                               0.173
#SRMR                                           0.124       0.124

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .114

CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.436       0.411
#Tucker-Lewis Index (TLI)                       0.275       0.242
#Robust Comparative Fit Index (CFI)                         0.434
#Robust Tucker-Lewis Index (TLI)                            0.273
#RMSEA                                          0.140       0.137
#Robust RMSEA                                               0.140
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .034

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


EGAmodel= '
 factor1=~ TIPI01+TIPI06
 factor2=~ TIPI02+TIPI07
 factor3=~ TIPI03+TIPI08
 factor4=~ TIPI04+TIPI09
 factor5=~ TIPI05+TIPI10
'

CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.880       0.798
#Tucker-Lewis Index (TLI)                       0.785       0.637
#Robust Comparative Fit Index (CFI)                         0.769
#Robust Tucker-Lewis Index (TLI)                            0.585
#RMSEA                                          0.104       0.118
#Robust RMSEA                                               0.124
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .296

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.216    0.226   -0.953    0.340   -0.216   -0.216
#    factor3           0.074    0.053    1.411    0.158    0.074    0.074
#    factor4          -0.074    0.060   -1.231    0.218   -0.074   -0.074
#    factor5           0.375    0.071    5.268    0.000    0.375    0.375
#  factor2 ~~                                                            
#    factor3          -0.113    0.125   -0.900    0.368   -0.113   -0.113
#    factor4           0.161    0.173    0.933    0.351    0.161    0.161
#    factor5          -0.141    0.155   -0.910    0.363   -0.141   -0.141
#  factor3 ~~                                                            
#    factor4          -0.178    0.079   -2.266    0.023   -0.178   -0.178
#    factor5           0.139    0.077    1.797    0.072    0.139    0.139
#  factor4 ~~                                                            
#    factor5          -0.157    0.079   -2.002    0.045   -0.157   -0.157

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning that did not result in data
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .32

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.591       0.651
#Tucker-Lewis Index (TLI)                       0.474       0.551
#Robust Comparative Fit Index (CFI)                         0.603
#Robust Tucker-Lewis Index (TLI)                            0.490
#RMSEA                                          0.119       0.106
#Robust RMSEA                                               0.117
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .474


#bifactor model 
BIFmodel= '
 factor1=~ TIPI01+TIPI06
 factor2=~ TIPI02+TIPI07
 factor3=~ TIPI03+TIPI08
 factor4=~ TIPI04+TIPI09
 factor5=~ TIPI05+TIPI10
 global =~  TIPI01+TIPI02+TIPI03+TIPI04+TIPI05+TIPI06+TIPI07+TIPI08+TIPI09+TIPI10
'

CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
# did not converge
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.807       0.833
#Tucker-Lewis Index (TLI)                       0.653       0.699
#Robust Comparative Fit Index (CFI)                         0.815
#Robust Tucker-Lewis Index (TLI)                            0.667
#RMSEA                                          0.097       0.087
#Robust RMSEA                                               0.094
#SRMR                                           0.070       0.070

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    TIPI01            1.001    0.151    6.609    0.000    1.001    0.562
#    TIPI06           -1.147    0.072  -15.934    0.000   -1.147   -0.597
#  factor2 =~                                                            
#    TIPI02            1.315    0.067   19.553    0.000    1.315    0.725
#    TIPI07           -0.356    0.102   -3.507    0.000   -0.356   -0.272
#  factor3 =~                                                            
#    TIPI03            1.003    0.089   11.207    0.000    1.003    0.691
#    TIPI08           -0.818    0.066  -12.405    0.000   -0.818   -0.441
#  factor4 =~                                                            
#    TIPI04            1.314    0.075   17.436    0.000    1.314    0.723
#    TIPI09           -1.119    0.058  -19.262    0.000   -1.119   -0.619
#  factor5 =~                                                            
#    TIPI05            0.895    2.403    0.372    0.710    0.895    0.657
#    TIPI10           -0.580    1.432   -0.405    0.686   -0.580   -0.351
#  global =~                                                             
#    TIPI01            1.061    0.158    6.713    0.000    1.061    0.596
#    TIPI02            0.241    0.193    1.250    0.211    0.241    0.133
#    TIPI03            0.416    0.126    3.306    0.001    0.416    0.287
#    TIPI04           -0.314    0.160   -1.961    0.050   -0.314   -0.173
#    TIPI05            0.400    0.099    4.048    0.000    0.400    0.294
#    TIPI06           -0.558    0.197   -2.838    0.005   -0.558   -0.291
#    TIPI07            0.788    0.118    6.655    0.000    0.788    0.601
#    TIPI08            0.055    0.137    0.402    0.688    0.055    0.030
#    TIPI09            0.252    0.167    1.513    0.130    0.252    0.140
#    TIPI10           -0.294    0.122   -2.419    0.016   -0.294   -0.178

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .48

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.global           PUC  Omega.global OmegaH.global 
#0.2389656     0.8888889     0.3063511     0.2595142 

#$FactorLevelIndices
#         ECV_SS    ECV_SG     ECV_GS      Omega      OmegaH         H        FD
#factor1 0.6048682 0.1504177 0.39513177 0.09612334 0.001253817 0.5041332 0.7138805
#factor2 0.6125789 0.1338739 0.38742105 0.42103078 0.116222016 0.5423551 0.7487709
#factor3 0.8898985 0.1501349 0.11010149 0.11526677 0.044185701 0.5355983 0.7453857
#factor4 0.9483414 0.2025210 0.05165862 0.01117518 0.010140264 0.6317635 0.7967356
#factor5 0.8244040 0.1240868 0.17559602 0.07477276 0.065418044 0.4738583 0.6989244
#global  0.2389656 0.2389656 0.23896559 0.30635107 0.259514162 0.5990903 0.7719368









####### OCEAN.20
### Boutilier https://psyarxiv.com/cmtz3/download?format=pdf
### https://osf.io/28b9q

### Not included because overlap with Hope

library(haven)
Boutilier_RRB_Hope_Optimism_Data <- read_sav("Boutilier_RRB_Hope_Optimism_Data.sav")
colnames(Boutilier_RRB_Hope_Optimism_Data)
mydata <- as.data.frame(Boutilier_RRB_Hope_Optimism_Data[,c(6:14,16,18,20,22:29)])
colnames(mydata)
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 components
# Eigenvalue1 = 3.2, eigenvalue2 = 2.59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components 
# Eigenvalue 1 = 3.43; eigenvalue 2 = 2.78


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.197, RMSR=.21, TLI=.19

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.17, RMSEA=.223, RMSR=.23, TLI=.154


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.61, RMSEA=.031, RMSR=.02, TLI=.98
#      MR1  MR3  MR2   MR5   MR4
#MR1  1.00 0.19 0.07 -0.16  0.05
#MR3  0.19 1.00 0.01  0.17  0.11
#MR2  0.07 0.01 1.00  0.14  0.21
#MR5 -0.16 0.17 0.14  1.00 -0.08
#MR4  0.05 0.11 0.21 -0.08  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.66, RMSEA=.06, RMSR=.02, TLI=.938
#      MR1  MR3  MR2   MR5   MR4
#MR1  1.00 0.22 0.07 -0.16  0.05
#MR3  0.22 1.00 0.01  0.16  0.08
#MR2  0.07 0.01 1.00  0.14  0.22
#MR5 -0.16 0.16 0.14  1.00 -0.08
#MR4  0.05 0.08 0.22 -0.08  1.00

colnames(mydata)

#Analysis lavaan
library(lavaan)

UNImodel= '
 global =~  O1+O2+O3+O4+C1+C2+C3+C4+E1+E2+E3+E4+A1+A2+A3+A4+N1+N2+N3+N4
'
CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
# gives warning
CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.212       0.173
#Tucker-Lewis Index (TLI)                       0.119       0.076
#Robust Comparative Fit Index (CFI)                         0.209
#Robust Tucker-Lewis Index (TLI)                            0.116
#RMSEA                                          0.210       0.199
#Robust RMSEA                                               0.209
#SRMR                                           0.212       0.212

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .028

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


EGAmodel= '
 factor1=~ O1+O2+O3+O4
 factor2=~ C1+C2+C3+C4
 factor3=~ E1+E2+E3+E4
 factor4=~ A1+A2+A3+A4
 factor5=~ N1+N2+N3+N4
'

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.939       0.942
#Tucker-Lewis Index (TLI)                       0.927       0.931
#Robust Comparative Fit Index (CFI)                         0.945
#Robust Tucker-Lewis Index (TLI)                            0.935
#RMSEA                                          0.060       0.054
#Robust RMSEA                                               0.057
#SRMR                                           0.065       0.065

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.059    0.080    0.744    0.457    0.059    0.059
#    factor3           0.129    0.082    1.578    0.115    0.129    0.129
#    factor4           0.232    0.079    2.933    0.003    0.232    0.232
#    factor5          -0.113    0.100   -1.131    0.258   -0.113   -0.113
#  factor2 ~~                                                            
#    factor3           0.216    0.086    2.500    0.012    0.216    0.216
#    factor4           0.087    0.078    1.123    0.262    0.087    0.087
#    factor5          -0.145    0.082   -1.760    0.078   -0.145   -0.145
#  factor3 ~~                                                            
#    factor4           0.026    0.088    0.291    0.771    0.026    0.026
#    factor5           0.195    0.092    2.125    0.034    0.195    0.195
#  factor4 ~~                                                            
#    factor5           0.217    0.078    2.786    0.005    0.217    0.217

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .587

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.921       0.924
#Tucker-Lewis Index (TLI)                       0.912       0.915
#Robust Comparative Fit Index (CFI)                         0.928
#Robust Tucker-Lewis Index (TLI)                            0.919
#RMSEA                                          0.066       0.060
#Robust RMSEA                                               0.063
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .588


#bifactor model 
BIFmodel= '
 factor1=~ O1+O2+O3+O4
 factor2=~ C1+C2+C3+C4
 factor3=~ E1+E2+E3+E4
 factor4=~ A1+A2+A3+A4
 factor5=~ N1+N2+N3+N4
 global =~  O1+O2+O3+O4+C1+C2+C3+C4+E1+E2+E3+E4+A1+A2+A3+A4+N1+N2+N3+N4
'

CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.954       0.948
#Tucker-Lewis Index (TLI)                       0.942       0.934
#Robust Comparative Fit Index (CFI)                         0.954
#Robust Tucker-Lewis Index (TLI)                            0.942
#RMSEA                                          0.054       0.053
#Robust RMSEA                                               0.053
#SRMR                                           0.091       0.091

# Because of high values Std.lv below, Std.all taken to calculate average loading
#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    O1                0.731    0.114    6.437    0.000    0.731    0.536
#    O2                1.268    0.100   12.728    0.000    1.268    0.850
#    O3                1.176    0.113   10.386    0.000    1.176    0.725
#    O4                1.238    0.117   10.582    0.000    1.238    0.659
#  factor2 =~                                                            
#    C1                1.344    0.196    6.868    0.000    1.344    0.811
#    C2                1.237    0.173    7.148    0.000    1.237    0.755
#    C3                1.146    0.127    9.018    0.000    1.146    0.718
#    C4                1.234    0.229    5.398    0.000    1.234    0.739
#  factor3 =~                                                            
#    E1                1.320    0.328    4.019    0.000    1.320    0.756
#    E2                1.549    0.263    5.894    0.000    1.549    0.799
#    E3                1.438    0.146    9.859    0.000    1.438    0.811
#    E4                1.595    0.202    7.914    0.000    1.595    0.816
#  factor4 =~                                                            
#    A1                0.793    0.104    7.644    0.000    0.793    0.621
#    A2                1.222    0.097   12.539    0.000    1.222    0.898
#    A3                1.125    0.103   10.969    0.000    1.125    0.759
#    A4                0.953    0.099    9.667    0.000    0.953    0.755
#  factor5 =~                                                            
#    N1                1.350    0.143    9.443    0.000    1.350    0.656
#    N2                1.564    0.299    5.233    0.000    1.564    0.756
#    N3                1.277    0.205    6.235    0.000    1.277    0.634
#    N4                1.192    0.284    4.199    0.000    1.192    0.610
#  global =~                                                             
#    O1                0.124    0.134    0.924    0.355    0.124    0.091
#    O2                0.141    0.164    0.858    0.391    0.141    0.094
#    O3               -0.016    0.268   -0.059    0.953   -0.016   -0.010
#    O4                0.328    0.212    1.546    0.122    0.328    0.174
#    C1                0.641    0.383    1.673    0.094    0.641    0.387
#    C2                0.708    0.306    2.314    0.021    0.708    0.432
#    C3                0.466    0.227    2.058    0.040    0.466    0.292
#    C4                0.774    0.404    1.916    0.055    0.774    0.464
#    E1                0.732    0.666    1.099    0.272    0.732    0.419
#    E2                0.547    0.697    0.785    0.433    0.547    0.282
#    E3               -0.211    0.706   -0.299    0.765   -0.211   -0.119
#    E4               -0.411    0.810   -0.507    0.612   -0.411   -0.210
#    A1               -0.001    0.178   -0.006    0.995   -0.001   -0.001
#    A2                0.251    0.302    0.831    0.406    0.251    0.184
#    A3                0.125    0.258    0.486    0.627    0.125    0.085
#    A4                0.116    0.213    0.545    0.586    0.116    0.092
#    N1               -0.371    0.442   -0.839    0.402   -0.371   -0.180
#    N2               -0.683    0.645   -1.059    0.290   -0.683   -0.330
#    N3               -0.606    0.397   -1.524    0.127   -0.606   -0.301
#    N4               -0.465    0.643   -0.723    0.470   -0.465   -0.238

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .637

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.global           PUC  Omega.global OmegaH.global 
#0.1091324     0.8421053     0.8551177     0.0482618 

#$FactorLevelIndices
#         ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.9763745 0.1607791 0.02362553 0.7970545 0.7845352 0.8296843 0.9110781
#factor2 0.7824958 0.1870545 0.21750423 0.9154658 0.7201627 0.8466100 0.9031866
#factor3 0.8898310 0.2069857 0.11016895 0.8991185 0.8869918 0.8753357 0.9449480
#factor4 0.9792219 0.1909223 0.02077808 0.8525742 0.8407305 0.8820336 0.9411779
#factor5 0.8602807 0.1451259 0.13971929 0.8083307 0.6992854 0.7705844 0.8724578
#global  0.1091324 0.1091324 0.10913241 0.8551177 0.0482618 0.6075646 0.7841816









####### HEXACO Grabovac https://osf.io/quspg

### not included because other dataset used

osf_s4fpj_Grabova <- read.csv("~/Self regulation/osf_s4fpj_Grabova.csv")
colnames(osf_s4fpj_Grabova)
mydata <- as.data.frame(osf_s4fpj_Grabova[,69:128])
colnames(mydata)
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 components
# Eigenvalue1 = 5.18, eigenvalue2 = 3.7

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 components 
# Eigenvalue 1 = 6.35; eigenvalue 2 = 3.99


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.09, RMSEA=.086, RMSR=.13, TLI=.167

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.11, RMSEA=.121, RMSR=.14, TLI=.11


# Give solution with 6 factors
fit5 <- fa(mydata,6)
fit5
diagram(fit5)
# %variance explained=.34, RMSEA=.053, RMSR=.06, TLI=.673
#      MR2   MR3   MR1  MR4   MR5   MR6
#MR2  1.00 -0.04  0.07 0.02 -0.04 -0.09
#MR3 -0.04  1.00 -0.12 0.03 -0.02  0.16
#MR1  0.07 -0.12  1.00 0.08 -0.06 -0.03
#MR4  0.02  0.03  0.08 1.00  0.05  0.08
#MR5 -0.04 -0.02 -0.06 0.05  1.00 -0.04
#MR6 -0.09  0.16 -0.03 0.08 -0.04  1.00

fit6 <- fa(rho,6,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.39, RMSEA=.098, RMSR=.06, TLI=.394
#      MR4   MR1   MR2   MR3   MR5   MR6
#MR4  1.00 -0.14 -0.04 -0.03  0.00  0.16
#MR1 -0.14  1.00  0.13  0.06 -0.09 -0.03
#MR2 -0.04  0.13  1.00  0.03  0.04  0.09
#MR3 -0.03  0.06  0.03  1.00 -0.05 -0.10
#MR5  0.00 -0.09  0.04 -0.05  1.00 -0.04
#MR6  0.16 -0.03  0.09 -0.10 -0.04  1.00

colnames(mydata)

# Lavaan
UNImodel= '
 general =~ h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+
            h21+h22+h23+h24+h25+h26+h27+h28+h29+h30+h31+h32+h33+h34+h35+h36+h37+h38+h39+h40+
            h41+h42+h43+h44+h45+h46+h47+h48+h49+h50+h51+h52+h53+h54+h55+h56+h57+h58+h59+h60
'

library(lavaan)
CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
# gives warning
CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.186       0.164
#Tucker-Lewis Index (TLI)                       0.157       0.135
#Robust Comparative Fit Index (CFI)                         0.180
#Robust Tucker-Lewis Index (TLI)                            0.151
#RMSEA                                          0.095       0.093
#Robust RMSEA                                               0.094
#SRMR                                           0.124       0.124

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .335

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

EGAmodel= '
 honesty =~ h6+h30+h54+h12+h36+h60+h18+h42+h24+h48
 emotion =~ h5+h29+h53+h11+h35+h17+h41+h23+h47+h59
 extrave =~ h4+h28+h52+h10+h34+h58+h16+h40+h22+h46
 agreeab =~ h3+h27+h9+h33+h51+h15+h39+h57+h21+h45
 conscie =~ h2+h26+h8+h32+h14+h38+h50+h20+h44+h56
 opennes =~ h1+h25+h7+h31+h13+h37+h49+h19+h43+h55
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.523       0.519
#Tucker-Lewis Index (TLI)                       0.502       0.497
#Robust Comparative Fit Index (CFI)                         0.526
#Robust Tucker-Lewis Index (TLI)                            0.506
#RMSEA                                          0.073       0.071
#Robust RMSEA                                               0.072
#SRMR                                           0.098       0.098

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  honesty ~~                                                            
#    emotion          -0.060    0.107   -0.556    0.578   -0.060   -0.060
#    extrave          -0.077    0.146   -0.526    0.599   -0.077   -0.077
#    agreeab           0.274    0.132    2.082    0.037    0.274    0.274
#    conscie           0.454    0.097    4.663    0.000    0.454    0.454
#    opennes          -0.101    0.100   -1.014    0.311   -0.101   -0.101
#  emotion ~~                                                            
#    extrave          -0.187    0.142   -1.313    0.189   -0.187   -0.187
#    agreeab          -0.206    0.112   -1.838    0.066   -0.206   -0.206
#    conscie          -0.095    0.102   -0.933    0.351   -0.095   -0.095
#    opennes          -0.147    0.105   -1.408    0.159   -0.147   -0.147
#  extrave ~~                                                            
#    agreeab          -0.113    0.119   -0.953    0.341   -0.113   -0.113
#    conscie           0.245    0.125    1.971    0.049    0.245    0.245
#    opennes          -0.090    0.083   -1.087    0.277   -0.090   -0.090
#  agreeab ~~                                                            
#    conscie           0.027    0.104    0.259    0.796    0.027    0.027
#    opennes          -0.101    0.094   -1.072    0.284   -0.101   -0.101
#  conscie ~~                                                            
#    opennes          -0.010    0.093   -0.104    0.917   -0.010   -0.010


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .273

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.509       0.507
#Tucker-Lewis Index (TLI)                       0.492       0.490
#Robust Comparative Fit Index (CFI)                         0.514
#Robust Tucker-Lewis Index (TLI)                            0.497
#RMSEA                                          0.074       0.071
#Robust RMSEA                                               0.072
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .269


#bifactor model 
BIFmodel= '
 honesty =~ h6+h30+h54+h12+h36+h60+h18+h42+h24+h48
 emotion =~ h5+h29+h53+h11+h35+h17+h41+h23+h47+h59
 extrave =~ h4+h28+h52+h10+h34+h58+h16+h40+h22+h46
 agreeab =~ h3+h27+h9+h33+h51+h15+h39+h57+h21+h45
 conscie =~ h2+h26+h8+h32+h14+h38+h50+h20+h44+h56
 opennes =~ h1+h25+h7+h31+h13+h37+h49+h19+h43+h55
 general =~ h1+h2+h3+h4+h5+h6+h7+h8+h9+h10+h11+h12+h13+h14+h15+h16+h17+h18+h19+h20+
            h21+h22+h23+h24+h25+h26+h27+h28+h29+h30+h31+h32+h33+h34+h35+h36+h37+h38+h39+h40+
            h41+h42+h43+h44+h45+h46+h47+h48+h49+h50+h51+h52+h53+h54+h55+h56+h57+h58+h59+h60
'
CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.603       0.562
#Tucker-Lewis Index (TLI)                       0.574       0.531
#Robust Comparative Fit Index (CFI)                         0.588
#Robust Tucker-Lewis Index (TLI)                            0.558
#RMSEA                                          0.067       0.068
#Robust RMSEA                                               0.068
#SRMR                                           0.092       0.092

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  honesty =~                                                            
#    h6                0.557    0.370    1.506    0.132    0.557    0.408
#    h30              -0.402    0.360   -1.117    0.264   -0.402   -0.368
#    h54               0.306    0.558    0.548    0.584    0.306    0.226
#    h12              -0.534    0.292   -1.830    0.067   -0.534   -0.394
#    h36               0.408    0.306    1.332    0.183    0.408    0.305
#    h60              -0.481    0.261   -1.844    0.065   -0.481   -0.411
#    h18               0.604    0.293    2.059    0.040    0.604    0.520
#    h42              -0.741    0.410   -1.808    0.071   -0.741   -0.554
#    h24              -0.488    0.133   -3.665    0.000   -0.488   -0.410
#    h48              -0.634    0.169   -3.754    0.000   -0.634   -0.549
#  emotion =~                                                            
#    h5                0.233    0.141    1.645    0.100    0.233    0.193
#    h29               0.280    0.157    1.788    0.074    0.280    0.219
#    h53              -0.425    0.129   -3.301    0.001   -0.425   -0.341
#    h11               0.192    0.151    1.272    0.203    0.192    0.170
#    h35              -0.179    0.118   -1.515    0.130   -0.179   -0.164
#    h17               0.738    0.124    5.968    0.000    0.738    0.566
#    h41              -0.597    0.153   -3.897    0.000   -0.597   -0.490
#    h23               0.696    0.110    6.334    0.000    0.696    0.539
#    h47               0.645    0.131    4.918    0.000    0.645    0.537
#    h59              -0.989    0.106   -9.342    0.000   -0.989   -0.776
#  extrave =~                                                            
#    h4                0.454    0.178    2.548    0.011    0.454    0.448
#    h28              -0.610    0.132   -4.604    0.000   -0.610   -0.528
#    h52              -0.495    0.266   -1.861    0.063   -0.495   -0.363
#    h10              -0.669    0.113   -5.907    0.000   -0.669   -0.500
#    h34               0.696    0.128    5.425    0.000    0.696    0.610
#    h58               0.808    0.138    5.859    0.000    0.808    0.626
#    h16               0.463    0.119    3.879    0.000    0.463    0.365
#    h40               0.873    0.100    8.748    0.000    0.873    0.698
#    h22               0.648    0.158    4.091    0.000    0.648    0.605
#    h46              -0.643    0.091   -7.086    0.000   -0.643   -0.553
#  agreeab =~                                                            
#    h3                0.644    0.114    5.655    0.000    0.644    0.508
#    h27               0.508    0.113    4.508    0.000    0.508    0.382
#    h9               -0.769    0.098   -7.856    0.000   -0.769   -0.599
#    h33               0.806    0.094    8.554    0.000    0.806    0.655
#    h51               0.622    0.081    7.662    0.000    0.622    0.602
#    h15              -0.520    0.110   -4.742    0.000   -0.520   -0.421
#    h39               0.343    0.100    3.434    0.001    0.343    0.310
#    h57              -0.441    0.108   -4.100    0.000   -0.441   -0.364
#    h21              -0.592    0.117   -5.068    0.000   -0.592   -0.482
#    h45               0.630    0.099    6.352    0.000    0.630    0.496
#  conscie =~                                                            
#    h2                0.779    0.178    4.384    0.000    0.779    0.640
#    h26              -0.736    0.203   -3.629    0.000   -0.736   -0.579
#    h8                0.366    0.124    2.964    0.003    0.366    0.484
#    h32              -0.374    0.190   -1.973    0.048   -0.374   -0.333
#    h14              -0.457    0.103   -4.436    0.000   -0.457   -0.483
#    h38               0.711    0.100    7.082    0.000    0.711    0.689
#    h50               0.882    0.165    5.335    0.000    0.882    0.637
#    h20              -0.450    0.221   -2.034    0.042   -0.450   -0.389
#    h44              -0.414    0.235   -1.762    0.078   -0.414   -0.382
#    h56              -0.587    0.163   -3.596    0.000   -0.587   -0.516
#  opennes =~                                                            
#    h1                0.737    0.091    8.064    0.000    0.737    0.662
#    h25              -0.695    0.101   -6.868    0.000   -0.695   -0.531
#    h7               -0.351    0.113   -3.115    0.002   -0.351   -0.237
#    h31               0.436    0.092    4.724    0.000    0.436    0.424
#    h13              -0.835    0.109   -7.668    0.000   -0.835   -0.687
#    h37              -0.681    0.101   -6.758    0.000   -0.681   -0.546
#    h49               0.767    0.101    7.617    0.000    0.767    0.661
#    h19               0.394    0.068    5.790    0.000    0.394    0.486
#    h43              -0.437    0.076   -5.748    0.000   -0.437   -0.445
#    h55               0.709    0.112    6.321    0.000    0.709    0.546
#  general =~                                                            
#    h1                0.085    0.108    0.790    0.430    0.085    0.076
#    h2               -0.119    0.447   -0.267    0.789   -0.119   -0.098
#    h3                0.000    0.222    0.000    1.000    0.000    0.000
#    h4               -0.288    0.398   -0.723    0.469   -0.288   -0.284
#    h5                0.234    0.162    1.439    0.150    0.234    0.194
#    h6               -0.290    0.339   -0.856    0.392   -0.290   -0.212
#    h7                0.132    0.171    0.771    0.441    0.132    0.089
#    h8               -0.090    0.330   -0.273    0.785   -0.090   -0.119
#    h9                0.317    0.200    1.582    0.114    0.317    0.247
#    h10               0.224    0.157    1.432    0.152    0.224    0.168
#    h11               0.394    0.307    1.287    0.198    0.394    0.350
#    h12               0.473    0.233    2.030    0.042    0.473    0.349
#    h13               0.256    0.200    1.285    0.199    0.256    0.211
#    h14               0.046    0.129    0.355    0.723    0.046    0.048
#    h15               0.182    0.231    0.786    0.432    0.182    0.147
#    h16               0.083    0.146    0.569    0.569    0.083    0.066
#    h17               0.316    0.144    2.198    0.028    0.316    0.243
#    h18              -0.033    0.139   -0.237    0.812   -0.033   -0.028
#    h19               0.124    0.104    1.191    0.234    0.124    0.153
#    h20               0.582    0.131    4.435    0.000    0.582    0.502
#    h21               0.532    0.181    2.944    0.003    0.532    0.433
#    h22              -0.204    0.418   -0.487    0.626   -0.204   -0.190
#    h23               0.379    0.222    1.706    0.088    0.379    0.293
#    h24               0.362    0.144    2.512    0.012    0.362    0.304
#    h25              -0.078    0.264   -0.298    0.766   -0.078   -0.060
#    h26               0.430    0.324    1.329    0.184    0.430    0.339
#    h27              -0.084    0.292   -0.288    0.773   -0.084   -0.063
#    h28               0.383    0.145    2.633    0.008    0.383    0.332
#    h29               0.240    0.189    1.267    0.205    0.240    0.187
#    h30               0.398    0.237    1.683    0.092    0.398    0.364
#    h31               0.083    0.118    0.709    0.478    0.083    0.081
#    h32               0.388    0.236    1.644    0.100    0.388    0.345
#    h33               0.163    0.144    1.132    0.258    0.163    0.133
#    h34               0.228    0.182    1.254    0.210    0.228    0.200
#    h35              -0.205    0.175   -1.166    0.244   -0.205   -0.187
#    h36              -0.242    0.255   -0.952    0.341   -0.242   -0.181
#    h37               0.181    0.310    0.584    0.559    0.181    0.145
#    h38              -0.152    0.373   -0.408    0.683   -0.152   -0.147
#    h39               0.407    0.125    3.245    0.001    0.407    0.368
#    h40               0.192    0.266    0.721    0.471    0.192    0.153
#    h41              -0.134    0.345   -0.387    0.699   -0.134   -0.110
#    h42               0.259    0.237    1.095    0.274    0.259    0.194
#    h43               0.217    0.263    0.827    0.408    0.217    0.222
#    h44               0.629    0.113    5.581    0.000    0.629    0.581
#    h45               0.065    0.185    0.351    0.725    0.065    0.051
#    h46               0.292    0.215    1.358    0.174    0.292    0.251
#    h47               0.320    0.312    1.024    0.306    0.320    0.267
#    h48               0.310    0.125    2.475    0.013    0.310    0.268
#    h49              -0.059    0.117   -0.503    0.615   -0.059   -0.051
#    h50               0.186    0.338    0.551    0.581    0.186    0.135
#    h51               0.136    0.122    1.115    0.265    0.136    0.132
#    h52               0.853    0.192    4.442    0.000    0.853    0.627
#    h53              -0.204    0.272   -0.749    0.454   -0.204   -0.164
#    h54              -0.220    0.467   -0.471    0.637   -0.220   -0.163
#    h55               0.068    0.127    0.535    0.593    0.068    0.052
#    h56               0.336    0.170    1.978    0.048    0.336    0.296
#    h57               0.497    0.134    3.709    0.000    0.497    0.409
#    h58               0.218    0.180    1.215    0.224    0.218    0.169
#    h59              -0.034    0.138   -0.249    0.803   -0.034   -0.027
#    h60               0.390    0.189    2.059    0.040    0.390    0.333

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .293

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.2006983      0.8474576      0.6469267      0.6040124 

#$FactorLevelIndices
#        ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#honesty 0.7298142 0.09783292 0.27018579 0.2858711 0.142802161 0.7010144 0.8315536
#emotion 0.8044429 0.10765898 0.19555708 0.1472760 0.023182103 0.7722676 0.8777151
#extrave 0.7833534 0.15724965 0.21664659 0.4010711 0.189175925 0.8176397 0.9095932
#agreeab 0.7986951 0.13146057 0.20130491 0.4001000 0.102222370 0.7780436 0.8879181
#conscie 0.7389584 0.14909247 0.26104162 0.3649241 0.005445433 0.8088565 0.8976895
#opennes 0.9453593 0.15600708 0.05464073 0.1210308 0.014127265 0.8190827 0.9071766
#general 0.2006983 0.20069833 0.20069833 0.6469267 0.604012364 0.8190903 0.9035965






####### IPIP
### Silvia (2022) https://osf.io/z694h/ https://www.sciencedirect.com/science/article/pii/S0092656609001317?via%3Dihub
### data https://osf.io/nwa3h

IPIP_Silvia <- read.csv("IPIP_Silvia.csv")
colnames(IPIP_Silvia)
mydata <- as.data.frame(IPIP_Silvia[,107:156])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 6 components
# Eigenvalue1 = 6.78, eigenvalue2 = 4.65

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 components 
# Eigenvalue 1 = 7.72; eigenvalue 2 = 5.3


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.112, RMSR=.16, TLI=.208

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.15, RMSEA=.199, RMSR=.18, TLI=.07


# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.43, RMSEA=.063, RMSR=.05, TLI=.741
#      MR1   MR3   MR2   MR4  MR5
#MR1  1.00 -0.14  0.05  0.20 0.08
#MR3 -0.14  1.00 -0.18 -0.03 0.00
#MR2  0.05 -0.18  1.00  0.01 0.21
#MR4  0.20 -0.03  0.01  1.00 0.08
#MR5  0.08  0.00  0.21  0.08 1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.50, RMSEA=.178, RMSR=.06, TLI=.236
#      MR1   MR3   MR2   MR4  MR5
#MR1  1.00 -0.13  0.06  0.19 0.07
#MR3 -0.13  1.00 -0.18 -0.01 0.00
#MR2  0.06 -0.18  1.00  0.02 0.22
#MR4  0.19 -0.01  0.02  1.00 0.09
#MR5  0.07  0.00  0.22  0.09 1.00

colnames(mydata)

#Analysis lavaan
library(lavaan)

UNImodel= '
 global =~  ipip1+ipip2+ipip3+ipip4+ipip5+ipip6+ipip7+ipip8+ipip9+ipip10+
            ipip11+ipip12+ipip13+ipip14+ipip15+ipip16+ipip17+ipip18+ipip19+ipip20+
            ipip21+ipip22+ipip23+ipip24+ipip25+ipip26+ipip27+ipip28+ipip29+ipip30+
            ipip31+ipip32+ipip33+ipip34+ipip35+ipip36+ipip37+ipip38+ipip39+ipip40+
            ipip41+ipip42+ipip43+ipip44+ipip45+ipip46+ipip47+ipip48+ipip49+ipip50
'

CFA_model1 <- cfa(UNImodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.535       0.363
#Tucker-Lewis Index (TLI)                       0.515       0.336
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.228       0.128
#Robust RMSEA                                                  NA
#SRMR                                           0.187       0.187

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .163

CFA_model1 <- cfa(UNImodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.250       0.232
#Tucker-Lewis Index (TLI)                       0.218       0.199
#Robust Comparative Fit Index (CFI)                         0.248
#Robust Tucker-Lewis Index (TLI)                            0.216
#RMSEA                                          0.120       0.116
#Robust RMSEA                                               0.119
#SRMR                                           0.156       0.156

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .382

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

EGAmodel= '
 factor1=~  ipip1+ipip2+ipip3+ipip4+ipip5+ipip6+ipip7+ipip8+ipip9+ipip10
 factor2=~  ipip11+ipip12+ipip13+ipip14+ipip15+ipip16+ipip17+ipip18+ipip19+ipip20
 factor3=~  ipip21+ipip22+ipip23+ipip24+ipip25+ipip26+ipip27+ipip28+ipip29+ipip30
 factor4=~  ipip31+ipip32+ipip33+ipip34+ipip35+ipip36+ipip37+ipip38+ipip39+ipip40
 factor5=~  ipip41+ipip42+ipip43+ipip44+ipip45+ipip46+ipip47+ipip48+ipip49+ipip50
'
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.914       0.848
#Tucker-Lewis Index (TLI)                       0.909       0.840
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.099       0.063
#Robust RMSEA                                                  NA
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .474

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.142    0.069    2.058    0.040    0.142    0.142
#    factor3           0.080    0.069    1.161    0.246    0.080    0.080
#    factor4          -0.139    0.067   -2.060    0.039   -0.139   -0.139
#    factor5           0.288    0.066    4.385    0.000    0.288    0.288
#  factor2 ~~                                                            
#    factor3           0.370    0.060    6.147    0.000    0.370    0.370
#    factor4          -0.040    0.067   -0.595    0.552   -0.040   -0.040
#    factor5           0.137    0.062    2.211    0.027    0.137    0.137
#  factor3 ~~                                                            
#    factor4          -0.238    0.064   -3.707    0.000   -0.238   -0.238
#    factor5           0.107    0.066    1.623    0.104    0.107    0.107
#  factor4 ~~                                                            
#    factor5          -0.064    0.064   -1.010    0.312   -0.064   -0.064


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.702       0.701
#Tucker-Lewis Index (TLI)                       0.687       0.686
#Robust Comparative Fit Index (CFI)                         0.708
#Robust Tucker-Lewis Index (TLI)                            0.693
#RMSEA                                          0.076       0.073
#Robust RMSEA                                               0.074
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .39

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.121    0.107    1.131    0.258    0.121    0.121
#    factor3           0.071    0.087    0.813    0.416    0.071    0.071
#    factor4          -0.151    0.096   -1.575    0.115   -0.151   -0.151
#    factor5           0.308    0.093    3.303    0.001    0.308    0.308
#  factor2 ~~                                                            
#    factor3           0.347    0.080    4.348    0.000    0.347    0.347
#    factor4           0.014    0.097    0.145    0.884    0.014    0.014
#    factor5           0.107    0.099    1.083    0.279    0.107    0.107
#  factor3 ~~                                                            
#    factor4          -0.232    0.091   -2.567    0.010   -0.232   -0.232
#    factor5           0.099    0.103    0.963    0.336    0.099    0.099
#  factor4 ~~                                                            
#    factor5          -0.092    0.095   -0.960    0.337   -0.092   -0.092

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.693       0.693
#Tucker-Lewis Index (TLI)                       0.680       0.680
#Robust Comparative Fit Index (CFI)                         0.700
#Robust Tucker-Lewis Index (TLI)                            0.687
#RMSEA                                          0.077       0.073
#Robust RMSEA                                               0.075
#SRMR                                           0.112       0.112

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .397


#bifactor model 
BIFmodel= '
 factor1=~  ipip1+ipip2+ipip3+ipip4+ipip5+ipip6+ipip7+ipip8+ipip9+ipip10
 factor2=~  ipip11+ipip12+ipip13+ipip14+ipip15+ipip16+ipip17+ipip18+ipip19+ipip20
 factor3=~  ipip21+ipip22+ipip23+ipip24+ipip25+ipip26+ipip27+ipip28+ipip29+ipip30
 factor4=~  ipip31+ipip32+ipip33+ipip34+ipip35+ipip36+ipip37+ipip38+ipip39+ipip40
 factor5=~  ipip41+ipip42+ipip43+ipip44+ipip45+ipip46+ipip47+ipip48+ipip49+ipip50
 global =~  ipip1+ipip2+ipip3+ipip4+ipip5+ipip6+ipip7+ipip8+ipip9+ipip10+
            ipip11+ipip12+ipip13+ipip14+ipip15+ipip16+ipip17+ipip18+ipip19+ipip20+
            ipip21+ipip22+ipip23+ipip24+ipip25+ipip26+ipip27+ipip28+ipip29+ipip30+
            ipip31+ipip32+ipip33+ipip34+ipip35+ipip36+ipip37+ipip38+ipip39+ipip40+
            ipip41+ipip42+ipip43+ipip44+ipip45+ipip46+ipip47+ipip48+ipip49+ipip50
'

CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.746       0.691
#Tucker-Lewis Index (TLI)                       0.724       0.663
#Robust Comparative Fit Index (CFI)                         0.727
#Robust Tucker-Lewis Index (TLI)                            0.702
#RMSEA                                          0.071       0.075
#Robust RMSEA                                               0.073
#SRMR                                           0.089       0.089

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    ipip1             0.805    0.130    6.193    0.000    0.805    0.696
#    ipip2             0.953    0.079   11.993    0.000    0.953    0.772
#    ipip3             0.814    0.126    6.451    0.000    0.814    0.657
#    ipip4             0.418    0.148    2.817    0.005    0.418    0.481
#    ipip5             0.782    0.100    7.840    0.000    0.782    0.623
#    ipip6             0.628    0.092    6.844    0.000    0.628    0.602
#    ipip7             0.835    0.236    3.543    0.000    0.835    0.644
#    ipip8             0.833    0.132    6.292    0.000    0.833    0.656
#    ipip9             0.746    0.218    3.427    0.001    0.746    0.689
#    ipip10            0.841    0.101    8.303    0.000    0.841    0.665
#  factor2 =~                                                            
#    ipip11            0.182    0.108    1.687    0.092    0.182    0.229
#    ipip12            0.441    0.089    4.983    0.000    0.441    0.613
#    ipip13            0.345    0.139    2.474    0.013    0.345    0.370
#    ipip14            0.547    0.084    6.537    0.000    0.547    0.624
#    ipip15            0.498    0.064    7.802    0.000    0.498    0.635
#    ipip16            0.427    0.105    4.065    0.000    0.427    0.422
#    ipip17            0.537    0.130    4.136    0.000    0.537    0.571
#    ipip18            0.254    0.163    1.554    0.120    0.254    0.280
#    ipip19            0.353    0.072    4.898    0.000    0.353    0.456
#    ipip20            0.406    0.088    4.615    0.000    0.406    0.572
#  factor3 =~                                                            
#    ipip21            0.625    0.111    5.627    0.000    0.625    0.657
#    ipip22            0.715    0.137    5.229    0.000    0.715    0.618
#    ipip23            0.475    0.096    4.958    0.000    0.475    0.487
#    ipip24            0.731    0.077    9.543    0.000    0.731    0.738
#    ipip25            0.775    0.089    8.748    0.000    0.775    0.705
#    ipip26            0.793    0.081    9.809    0.000    0.793    0.677
#    ipip27            0.653    0.110    5.936    0.000    0.653    0.579
#    ipip28            0.711    0.110    6.457    0.000    0.711    0.609
#    ipip29            0.555    0.128    4.342    0.000    0.555    0.526
#    ipip30            0.488    0.074    6.563    0.000    0.488    0.564
#  factor4 =~                                                            
#    ipip31            0.912    0.312    2.922    0.003    0.912    0.702
#    ipip32            0.824    0.352    2.344    0.019    0.824    0.676
#    ipip33            0.567    0.286    1.984    0.047    0.567    0.532
#    ipip34            0.686    0.099    6.922    0.000    0.686    0.630
#    ipip35            0.898    0.143    6.291    0.000    0.898    0.713
#    ipip36            0.307    0.310    0.989    0.323    0.307    0.256
#    ipip37            0.832    0.456    1.825    0.068    0.832    0.661
#    ipip38            0.937    0.456    2.057    0.040    0.937    0.742
#    ipip39            0.838    0.135    6.206    0.000    0.838    0.706
#    ipip40            0.557    0.182    3.064    0.002    0.557    0.531
#  factor5 =~                                                            
#    ipip41            0.844    0.239    3.527    0.000    0.844    0.829
#    ipip42            0.542    0.211    2.566    0.010    0.542    0.448
#    ipip43            0.246    0.816    0.302    0.763    0.246    0.235
#    ipip44            0.468    0.292    1.602    0.109    0.468    0.427
#    ipip45            0.405    0.234    1.731    0.083    0.405    0.435
#    ipip46            0.884    0.286    3.094    0.002    0.884    0.825
#    ipip47            0.381    0.791    0.482    0.630    0.381    0.425
#    ipip48            0.427    0.689    0.620    0.535    0.427    0.427
#    ipip49            0.501    0.979    0.512    0.609    0.501    0.493
#    ipip50            0.588    0.686    0.857    0.391    0.588    0.568
#  global =~                                                             
#    ipip1             0.287    0.262    1.097    0.272    0.287    0.248
#    ipip2             0.230    0.190    1.207    0.227    0.230    0.186
#    ipip3             0.102    0.171    0.596    0.551    0.102    0.083
#    ipip4             0.374    0.294    1.275    0.202    0.374    0.430
#    ipip5             0.231    0.165    1.400    0.161    0.231    0.184
#    ipip6             0.452    0.122    3.717    0.000    0.452    0.434
#    ipip7             0.638    0.456    1.398    0.162    0.638    0.492
#    ipip8             0.448    0.231    1.943    0.052    0.448    0.353
#    ipip9             0.307    0.454    0.677    0.498    0.307    0.283
#    ipip10            0.353    0.170    2.080    0.038    0.353    0.279
#    ipip11            0.329    0.249    1.324    0.186    0.329    0.415
#    ipip12            0.187    0.133    1.404    0.160    0.187    0.260
#    ipip13           -0.075    0.395   -0.189    0.850   -0.075   -0.080
#    ipip14            0.072    0.140    0.516    0.606    0.072    0.082
#    ipip15            0.195    0.135    1.447    0.148    0.195    0.248
#    ipip16            0.115    0.150    0.768    0.443    0.115    0.114
#    ipip17            0.111    0.283    0.393    0.694    0.111    0.118
#    ipip18           -0.048    0.370   -0.131    0.896   -0.048   -0.053
#    ipip19            0.281    0.155    1.810    0.070    0.281    0.362
#    ipip20            0.190    0.190    1.002    0.316    0.190    0.268
#    ipip21            0.189    0.381    0.495    0.620    0.189    0.198
#    ipip22           -0.100    0.319   -0.313    0.754   -0.100   -0.087
#    ipip23            0.285    0.277    1.028    0.304    0.285    0.291
#    ipip24            0.147    0.153    0.963    0.336    0.147    0.149
#    ipip25            0.008    0.188    0.044    0.965    0.008    0.008
#    ipip26            0.154    0.171    0.901    0.367    0.154    0.132
#    ipip27            0.216    0.356    0.608    0.543    0.216    0.192
#    ipip28           -0.121    0.187   -0.650    0.515   -0.121   -0.104
#    ipip29            0.236    0.231    1.019    0.308    0.236    0.223
#    ipip30            0.260    0.132    1.968    0.049    0.260    0.300
#    ipip31           -0.376    0.690   -0.544    0.586   -0.376   -0.289
#    ipip32           -0.306    0.917   -0.334    0.739   -0.306   -0.251
#    ipip33           -0.420    0.932   -0.451    0.652   -0.420   -0.394
#    ipip34            0.001    0.341    0.003    0.997    0.001    0.001
#    ipip35           -0.096    0.256   -0.376    0.707   -0.096   -0.076
#    ipip36           -0.195    1.074   -0.181    0.856   -0.195   -0.162
#    ipip37            0.132    0.683    0.193    0.847    0.132    0.105
#    ipip38            0.155    0.564    0.275    0.783    0.155    0.123
#    ipip39           -0.141    0.300   -0.471    0.637   -0.141   -0.119
#    ipip40           -0.232    0.635   -0.366    0.714   -0.232   -0.222
#    ipip41            0.033    0.604    0.055    0.956    0.033    0.032
#    ipip42            0.107    0.629    0.171    0.864    0.107    0.089
#    ipip43            0.446    0.906    0.493    0.622    0.446    0.426
#    ipip44            0.161    0.576    0.280    0.780    0.161    0.147
#    ipip45            0.306    0.191    1.601    0.109    0.306    0.329
#    ipip46            0.095    0.489    0.194    0.846    0.095    0.088
#    ipip47            0.533    0.873    0.611    0.542    0.533    0.595
#    ipip48            0.351    0.994    0.353    0.724    0.351    0.350
#    ipip49            0.580    1.244    0.467    0.641    0.580    0.571
#    ipip50            0.483    0.861    0.561    0.575    0.483    0.467

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .439

semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#ModelLevelIndices
#ECV.global           PUC  Omega.global OmegaH.global 
#0.1751774     0.8163265     0.8877702     0.2381080 

#$FactorLevelIndices
#         ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8043213 0.2008269 0.19567874 0.9153333 0.7563839 0.8868745 0.9324466
#factor2 0.8190430 0.1167846 0.18095701 0.7869435 0.6952392 0.7855597 0.8839308
#factor3 0.9147624 0.1817773 0.08523757 0.8726218 0.8352611 0.8695201 0.9321944
#factor4 0.9035877 0.1873992 0.09641228 0.8756500 0.8390079 0.8803534 0.9394468
#factor5 0.6846324 0.1380345 0.31536762 0.8618153 0.6308860 0.8638497 0.9233103
#global  0.1751774 0.1751774 0.17517742 0.8877702 0.2381080 0.8193195 0.8874817






####################################
####### HEXACO Gallagher 2023 https://www.tandfonline.com/doi/abs/10.1080/00223891.2022.2153690
### https://osf.io/tq64s/
### study 4 https://osf.io/9pwky

HEXACO_Gallagher <- read.csv("HEXACO_Gallagher.csv")
colnames(HEXACO_Gallagher)
mydata <- as.data.frame(HEXACO_Gallagher[,23:82])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 7 components
# Eigenvalue1 = 10, eigenvalue2 = 4.1

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 components 
# Eigenvalue 1 = 11.85; eigenvalue 2 = 4.99


# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.17, RMSEA=.101, RMSR=.12, TLI=.302

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.129, RMSR=.14, TLI=.247


# Give solution with 6 factors
fit5 <- fa(mydata,6)
fit5
diagram(fit5)
# %variance explained=.43, RMSEA=.06, RMSR=.04, TLI=.752
#      MR1   MR6  MR3   MR4   MR2   MR5
#MR1  1.00  0.24 0.12 -0.11  0.07  0.18
#MR6  0.24  1.00 0.18 -0.12  0.33  0.21
#MR3  0.12  0.18 1.00  0.03  0.12  0.17
#MR4 -0.11 -0.12 0.03  1.00 -0.04 -0.14
#MR2  0.07  0.33 0.12 -0.04  1.00  0.24
#MR5  0.18  0.21 0.17 -0.14  0.24  1.00

fit6 <- fa(rho,6,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.51, RMSEA=.089, RMSR=.05, TLI=.635
#      MR2   MR1  MR3   MR5   MR6   MR4
#MR2  1.00  0.23 0.10  0.17  0.06 -0.11
#MR1  0.23  1.00 0.19  0.22  0.35 -0.13
#MR3  0.10  0.19 1.00  0.20  0.14  0.03
#MR5  0.17  0.22 0.20  1.00  0.25 -0.12
#MR6  0.06  0.35 0.14  0.25  1.00 -0.04
#MR4 -0.11 -0.13 0.03 -0.12 -0.04  1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~ HEX_O_1+HEX_C_2+HEX_A_3+HEX_X_4+HEX_E_5+HEX_HH_6+HEX_O_7+HEX_C_8+  
           HEX_A_9+HEX_X_10+HEX_E_11+HEX_HH_12+HEX_O_13+HEX_C_14+HEX_A_15+HEX_X_16+ 
           HEX_E_17+HEX_HH_18+HEX_O_19+HEX_C_20+HEX_A_21+HEX_X_22+HEX_E_23+HEX_HH_24+
           HEX_O_25+HEX_C_26+HEX_A_27+HEX_X_28+HEX_E_29+HEX_HH_30+HEX_O_31+HEX_C_32+ 
           HEX_A_33+HEX_X_34+HEX_E_35+HEX_HH_36+HEX_O_37+HEX_C_38+HEX_A_39+HEX_X_40+ 
           HEX_E_41+HEX_HH_42+HEX_O_43+HEX_C_44+HEX_A_45+HEX_X_46+HEX_E_47+HEX_HH_48+
           HEX_O_49+HEX_C_50+HEX_A_51+HEX_X_52+HEX_E_53+HEX_HH_54+HEX_O_55+HEX_C_56+ 
           HEX_A_57+HEX_X_58+HEX_E_59+HEX_HH_60 
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.683       0.401
#Tucker-Lewis Index (TLI)                       0.672       0.380
#Robust Comparative Fit Index (CFI)                         0.279
#Robust Tucker-Lewis Index (TLI)                            0.254
#RMSEA                                          0.184       0.115
#Robust RMSEA                                               0.132
#SRMR                                           0.148       0.148

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .215

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.326       0.326
#Tucker-Lewis Index (TLI)                       0.302       0.302
#Robust Comparative Fit Index (CFI)                         0.331
#Robust Tucker-Lewis Index (TLI)                            0.308
#RMSEA                                          0.104       0.095
#Robust RMSEA                                               0.102
#SRMR                                           0.123       0.123

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .130

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

EGAmodel= '
 factor1=~ HEX_O_1+HEX_O_7+HEX_O_13+HEX_O_19+HEX_O_25+HEX_O_31+HEX_O_37+HEX_O_43+HEX_O_49+HEX_O_55
 factor2=~ HEX_C_2+HEX_C_8+HEX_C_14+HEX_C_20+HEX_C_26+HEX_C_32+HEX_C_38+HEX_C_44+HEX_C_50+HEX_C_56
 factor3=~ HEX_A_3+HEX_A_9+HEX_A_15+HEX_A_21+HEX_A_27+HEX_A_33+HEX_A_39+HEX_A_45+HEX_A_51+HEX_A_57
 factor4=~ HEX_X_4+HEX_X_10+HEX_X_16+HEX_X_22+HEX_X_28+HEX_X_34+HEX_X_40+HEX_X_46+HEX_X_52+HEX_X_58
 factor5=~ HEX_E_5+HEX_E_11+HEX_E_17+HEX_E_23+HEX_E_29+HEX_E_35+HEX_E_41+HEX_E_47+HEX_E_53+HEX_E_59
 factor6=~ HEX_HH_6+HEX_HH_12+HEX_HH_18+HEX_HH_24+HEX_HH_30+HEX_HH_36+HEX_HH_42+HEX_HH_48+HEX_HH_54+HEX_HH_60 
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.747
#Tucker-Lewis Index (TLI)                       0.881       0.736
#Robust Comparative Fit Index (CFI)                         0.576
#Robust Tucker-Lewis Index (TLI)                            0.557
#RMSEA                                          0.111       0.075
#Robust RMSEA                                               0.102
#SRMR                                           0.103       0.103

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.360    0.044    8.237    0.000    0.360    0.360
#    factor3           0.319    0.043    7.358    0.000    0.319    0.319
#    factor4           0.273    0.043    6.362    0.000    0.273    0.273
#    factor5          -0.039    0.051   -0.764    0.445   -0.039   -0.039
#    factor6           0.206    0.046    4.468    0.000    0.206    0.206
#  factor2 ~~                                                            
#    factor3           0.370    0.040    9.224    0.000    0.370    0.370
#    factor4           0.459    0.034   13.628    0.000    0.459    0.459
#    factor5          -0.247    0.041   -6.021    0.000   -0.247   -0.247
#    factor6           0.360    0.040    8.931    0.000    0.360    0.360
#  factor3 ~~                                                            
#    factor4           0.527    0.032   16.464    0.000    0.527    0.527
#    factor5          -0.343    0.041   -8.361    0.000   -0.343   -0.343
#    factor6           0.408    0.038   10.808    0.000    0.408    0.408
#  factor4 ~~                                                            
#    factor5          -0.552    0.032  -17.285    0.000   -0.552   -0.552
#    factor6           0.146    0.043    3.382    0.001    0.146    0.146
#  factor5 ~~                                                            
#    factor6          -0.112    0.046   -2.442    0.015   -0.112   -0.112

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .465

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.658       0.662
#Tucker-Lewis Index (TLI)                       0.643       0.647
#Robust Comparative Fit Index (CFI)                         0.669
#Robust Tucker-Lewis Index (TLI)                            0.654
#RMSEA                                          0.074       0.068
#Robust RMSEA                                               0.072
#SRMR                                           0.091       0.091

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .376

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.310    0.068    4.570    0.000    0.310    0.310
#    factor3           0.300    0.060    5.019    0.000    0.300    0.300
#    factor4           0.269    0.057    4.743    0.000    0.269    0.269
#    factor5          -0.003    0.071   -0.040    0.968   -0.003   -0.003
#    factor6           0.157    0.063    2.485    0.013    0.157    0.157
#  factor2 ~~                                                            
#    factor3           0.337    0.065    5.181    0.000    0.337    0.337
#    factor4           0.448    0.054    8.229    0.000    0.448    0.448
#    factor5          -0.264    0.061   -4.341    0.000   -0.264   -0.264
#    factor6           0.366    0.063    5.812    0.000    0.366    0.366
#  factor3 ~~                                                            
#    factor4           0.540    0.047   11.560    0.000    0.540    0.540
#    factor5          -0.314    0.084   -3.748    0.000   -0.314   -0.314
#    factor6           0.424    0.055    7.636    0.000    0.424    0.424
#  factor4 ~~                                                            
#    factor5          -0.519    0.067   -7.702    0.000   -0.519   -0.519
#    factor6           0.203    0.094    2.160    0.031    0.203    0.203
#  factor5 ~~                                                            
#    factor6          -0.125    0.067   -1.856    0.063   -0.125   -0.125

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.623       0.627
#Tucker-Lewis Index (TLI)                       0.610       0.614
#Robust Comparative Fit Index (CFI)                         0.634
#Robust Tucker-Lewis Index (TLI)                            0.621
#RMSEA                                          0.078       0.071
#Robust RMSEA                                               0.076
#SRMR                                           0.144       0.144

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .382


#bifactor model 
BIFmodel= '
 factor1=~ HEX_O_1+HEX_O_7+HEX_O_13+HEX_O_19+HEX_O_25+HEX_O_31+HEX_O_37+HEX_O_43+HEX_O_49+HEX_O_55
 factor2=~ HEX_C_2+HEX_C_8+HEX_C_14+HEX_C_20+HEX_C_26+HEX_C_32+HEX_C_38+HEX_C_44+HEX_C_50+HEX_C_56
 factor3=~ HEX_A_3+HEX_A_9+HEX_A_15+HEX_A_21+HEX_A_27+HEX_A_33+HEX_A_39+HEX_A_45+HEX_A_51+HEX_A_57
 factor4=~ HEX_X_4+HEX_X_10+HEX_X_16+HEX_X_22+HEX_X_28+HEX_X_34+HEX_X_40+HEX_X_46+HEX_X_52+HEX_X_58
 factor5=~ HEX_E_5+HEX_E_11+HEX_E_17+HEX_E_23+HEX_E_29+HEX_E_35+HEX_E_41+HEX_E_47+HEX_E_53+HEX_E_59
 factor6=~ HEX_HH_6+HEX_HH_12+HEX_HH_18+HEX_HH_24+HEX_HH_30+HEX_HH_36+HEX_HH_42+HEX_HH_48+HEX_HH_54+HEX_HH_60 
 general=~ HEX_O_1+HEX_C_2+HEX_A_3+HEX_X_4+HEX_E_5+HEX_HH_6+HEX_O_7+HEX_C_8+  
           HEX_A_9+HEX_X_10+HEX_E_11+HEX_HH_12+HEX_O_13+HEX_C_14+HEX_A_15+HEX_X_16+ 
           HEX_E_17+HEX_HH_18+HEX_O_19+HEX_C_20+HEX_A_21+HEX_X_22+HEX_E_23+HEX_HH_24+
           HEX_O_25+HEX_C_26+HEX_A_27+HEX_X_28+HEX_E_29+HEX_HH_30+HEX_O_31+HEX_C_32+ 
           HEX_A_33+HEX_X_34+HEX_E_35+HEX_HH_36+HEX_O_37+HEX_C_38+HEX_A_39+HEX_X_40+ 
           HEX_E_41+HEX_HH_42+HEX_O_43+HEX_C_44+HEX_A_45+HEX_X_46+HEX_E_47+HEX_HH_48+
           HEX_O_49+HEX_C_50+HEX_A_51+HEX_X_52+HEX_E_53+HEX_HH_54+HEX_O_55+HEX_C_56+ 
           HEX_A_57+HEX_X_58+HEX_E_59+HEX_HH_60 
'
CFA_model4 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.723       0.730
#Tucker-Lewis Index (TLI)                       0.703       0.710
#Robust Comparative Fit Index (CFI)                         0.735
#Robust Tucker-Lewis Index (TLI)                            0.716
#RMSEA                                          0.068       0.061
#Robust RMSEA                                               0.065
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .393

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3466601      0.8474576      0.8994756      0.4835475 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#factor1 0.9222891 0.14039966 0.07771095 0.8488834 0.79645426 0.8608966 0.9305001
#factor2 0.7669478 0.10465338 0.23305219 0.8205408 0.66210204 0.7858576 0.8916498
#factor3 0.6718580 0.11009850 0.32814201 0.8660406 0.57762645 0.8022563 0.9034976
#factor4 0.1646412 0.03781172 0.83535883 0.9155284 0.01615711 0.5175316 0.8298983
#factor5 0.7429609 0.12443191 0.25703914 0.8564133 0.72763292 0.8223570 0.9126686
#factor6 0.9044780 0.13594476 0.09552196 0.8427511 0.81500761 0.8399772 0.9201731
#general 0.3466601 0.34666007 0.34666007 0.8994756 0.48354753 0.9341425 0.9670737






####### BFAS study
####### Zajenkowski & Szymaniak (2021) https://link.springer.com/article/10.1007/s12144-019-0147-1
####### Data https://osf.io/7c6e3

library(haven)
BFAS_Zajenkowski <- read_sav("BFAS_Zajenkowski.sav")
colnames(BFAS_Zajenkowski)
mydata <- BFAS_Zajenkowski[,49:148]
#mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 5


library(psych)
omega(mydata) # alpha = .88, omega T = .9

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 9 components
# Eigenvalue1 = 10.29, eigenvalue2 = 8.75

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 10 factors and 10 components
# Eigenvalue 1 = 11.5; eigenvalue 2 = 10.26

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.10, RMSEA=.078, RMSR=.13, TLI=.197

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.11, RMSEA=.107, RMSR=.15, TLI=.129

# Give solution with 5 factors
fit5 <- fa(mydata,5)
fit5
diagram(fit5)
# %variance explained=.34, RMSEA=.048, RMSR=.05, TLI=.691
#      MR3  MR2   MR1   MR4   MR5
#MR3  1.00 0.02 -0.14 -0.14  0.12
#MR2  0.02 1.00  0.03  0.13  0.15
#MR1 -0.14 0.03  1.00  0.07 -0.14
#MR4 -0.14 0.13  0.07  1.00 -0.02
#MR5  0.12 0.15 -0.14 -0.02  1.00

fit6 <- fa(rho,5,n.obs=nrow(mydata))
fit6
diagram(fit6)
# %variance explained=.39, RMSEA=.083, RMSR=.06, TLI=.471
#      MR3   MR1  MR2   MR4   MR5
#MR3  1.00 -0.14 0.04 -0.12  0.12
#MR1 -0.14  1.00 0.02  0.07 -0.14
#MR2  0.04  0.02 1.00  0.16  0.15
#MR4 -0.12  0.07 0.16  1.00 -0.01
#MR5  0.12 -0.14 0.15 -0.01  1.00


#unifactorial model with Lavaan
UNImodel= '
global =~ IPIP1+IPIP2+IPIP3+IPIP4+IPIP5+IPIP6+IPIP7+IPIP8+IPIP9+IPIP10+
          IPIP11+IPIP12+IPIP13+IPIP14+IPIP15+IPIP16+IPIP17+IPIP18+IPIP19+IPIP20+
          IPIP21+IPIP22+IPIP23+IPIP24+IPIP25+IPIP26+IPIP27+IPIP28+IPIP29+IPIP30+
          IPIP31+IPIP32+IPIP33+IPIP34+IPIP35+IPIP36+IPIP37+IPIP38+IPIP39+IPIP40+
          IPIP41+IPIP42+IPIP43+IPIP44+IPIP45+IPIP46+IPIP47+IPIP48+IPIP49+IPIP50+
          IPIP51+IPIP52+IPIP53+IPIP54+IPIP55+IPIP56+IPIP57+IPIP58+IPIP59+IPIP60+
          IPIP61+IPIP62+IPIP63+IPIP64+IPIP65+IPIP66+IPIP67+IPIP68+IPIP69+IPIP70+
          IPIP71+IPIP72+IPIP73+IPIP74+IPIP75+IPIP76+IPIP77+IPIP78+IPIP79+IPIP80+
          IPIP81+IPIP82+IPIP83+IPIP84+IPIP85+IPIP86+IPIP87+IPIP88+IPIP89+IPIP90+
          IPIP91+IPIP92+IPIP93+IPIP94+IPIP95+IPIP96+IPIP97+IPIP98+IPIP99+IPIP100
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.453       0.310
#Tucker-Lewis Index (TLI)                       0.442       0.296
#Robust Comparative Fit Index (CFI)                         0.155
#Robust Tucker-Lewis Index (TLI)                            0.137
#RMSEA                                          0.171       0.082
#Robust RMSEA                                               0.112
#SRMR                                           0.156       0.156

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .079

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.209       0.201
#Tucker-Lewis Index (TLI)                       0.193       0.184
#Robust Comparative Fit Index (CFI)                         0.210
#Robust Tucker-Lewis Index (TLI)                            0.194
#RMSEA                                          0.083       0.080
#Robust RMSEA                                               0.082
#SRMR                                           0.131       0.131

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .048


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1=~ IPIP1+IPIP6+IPIP11+IPIP16+IPIP21+IPIP26+IPIP31+IPIP36+IPIP41+IPIP46+
           IPIP51+IPIP56+IPIP61+IPIP66+IPIP71+IPIP76+IPIP81+IPIP86+IPIP91+IPIP96
 factor2=~ IPIP2+IPIP7+IPIP12+IPIP17+IPIP22+IPIP27+IPIP32+IPIP37+IPIP42+IPIP47+
           IPIP52+IPIP57+IPIP62+IPIP67+IPIP72+IPIP77+IPIP82+IPIP87+IPIP92+IPIP97
 factor3=~ IPIP3+IPIP8+IPIP13+IPIP18+IPIP23+IPIP28+IPIP33+IPIP38+IPIP43+IPIP48+
           IPIP53+IPIP58+IPIP63+IPIP68+IPIP73+IPIP78+IPIP83+IPIP88+IPIP93+IPIP98
 factor4=~ IPIP4+IPIP9+IPIP14+IPIP19+IPIP24+IPIP29+IPIP34+IPIP39+IPIP44+IPIP49+
           IPIP54+IPIP59+IPIP64+IPIP69+IPIP74+IPIP79+IPIP84+IPIP89+IPIP94+IPIP99
 factor5=~ IPIP5+IPIP10+IPIP15+IPIP20+IPIP25+IPIP30+IPIP35+IPIP40+IPIP45+IPIP50+
           IPIP55+IPIP60+IPIP65+IPIP70+IPIP75+IPIP80+IPIP85+IPIP90+IPIP95+IPIP100
'
library(lavaan)
CFA_model1 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.676       0.574
#Tucker-Lewis Index (TLI)                       0.669       0.565
#Robust Comparative Fit Index (CFI)                         0.376
#Robust Tucker-Lewis Index (TLI)                            0.362
#RMSEA                                          0.131       0.065
#Robust RMSEA                                               0.096
#SRMR                                           0.130       0.130

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .284

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.190    0.042    4.561    0.000    0.190    0.190
#    factor3          -0.265    0.038   -6.971    0.000   -0.265   -0.265
#    factor4          -0.256    0.043   -5.958    0.000   -0.256   -0.256
#    factor5          -0.063    0.045   -1.393    0.164   -0.063   -0.063
#  factor2 ~~                                                            
#    factor3           0.132    0.047    2.801    0.005    0.132    0.132
#    factor4           0.092    0.044    2.074    0.038    0.092    0.092
#    factor5           0.386    0.041    9.356    0.000    0.386    0.386
#  factor3 ~~                                                            
#    factor4           0.296    0.042    7.128    0.000    0.296    0.296
#    factor5           0.299    0.039    7.765    0.000    0.299    0.299
#  factor4 ~~                                                            
#    factor5           0.559    0.031   18.210    0.000    0.559    0.559

CFA_model1 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.459       0.462
#Tucker-Lewis Index (TLI)                       0.446       0.449
#Robust Comparative Fit Index (CFI)                         0.466
#Robust Tucker-Lewis Index (TLI)                            0.454
#RMSEA                                          0.069       0.065
#Robust RMSEA                                               0.068
#SRMR                                           0.114       0.114

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.189    0.067    2.796    0.005    0.189    0.189
#    factor3          -0.268    0.069   -3.895    0.000   -0.268   -0.268
#    factor4          -0.313    0.064   -4.905    0.000   -0.313   -0.313
#    factor5          -0.080    0.120   -0.672    0.501   -0.080   -0.080
#  factor2 ~~                                                            
#    factor3           0.116    0.071    1.626    0.104    0.116    0.116
#    factor4          -0.060    0.078   -0.768    0.443   -0.060   -0.060
#    factor5           0.334    0.104    3.201    0.001    0.334    0.334
#  factor3 ~~                                                            
#    factor4           0.295    0.064    4.603    0.000    0.295    0.295
#    factor5           0.299    0.083    3.624    0.000    0.299    0.299
#  factor4 ~~                                                            
#    factor5           0.565    0.080    7.090    0.000    0.565    0.565

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .201

library(semPlot)
semPaths(CFA_model1,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.444       0.447
#Tucker-Lewis Index (TLI)                       0.433       0.436
#Robust Comparative Fit Index (CFI)                         0.452
#Robust Tucker-Lewis Index (TLI)                            0.440
#RMSEA                                          0.070       0.066
#Robust RMSEA                                               0.069
#SRMR                                           0.127       0.127

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .201


#bifactor model better
BIFmodel= '
 factor1=~ IPIP1+IPIP6+IPIP11+IPIP16+IPIP21+IPIP26+IPIP31+IPIP36+IPIP41+IPIP46+
           IPIP51+IPIP56+IPIP61+IPIP66+IPIP71+IPIP76+IPIP81+IPIP86+IPIP91+IPIP96
 factor2=~ IPIP2+IPIP7+IPIP12+IPIP17+IPIP22+IPIP27+IPIP32+IPIP37+IPIP42+IPIP47+
           IPIP52+IPIP57+IPIP62+IPIP67+IPIP72+IPIP77+IPIP82+IPIP87+IPIP92+IPIP97
 factor3=~ IPIP3+IPIP8+IPIP13+IPIP18+IPIP23+IPIP28+IPIP33+IPIP38+IPIP43+IPIP48+
           IPIP53+IPIP58+IPIP63+IPIP68+IPIP73+IPIP78+IPIP83+IPIP88+IPIP93+IPIP98
 factor4=~ IPIP4+IPIP9+IPIP14+IPIP19+IPIP24+IPIP29+IPIP34+IPIP39+IPIP44+IPIP49+
           IPIP54+IPIP59+IPIP64+IPIP69+IPIP74+IPIP79+IPIP84+IPIP89+IPIP94+IPIP99
 factor5=~ IPIP5+IPIP10+IPIP15+IPIP20+IPIP25+IPIP30+IPIP35+IPIP40+IPIP45+IPIP50+
           IPIP55+IPIP60+IPIP65+IPIP70+IPIP75+IPIP80+IPIP85+IPIP90+IPIP95+IPIP100
 global =~ IPIP1+IPIP2+IPIP3+IPIP4+IPIP5+IPIP6+IPIP7+IPIP8+IPIP9+IPIP10+
           IPIP11+IPIP12+IPIP13+IPIP14+IPIP15+IPIP16+IPIP17+IPIP18+IPIP19+IPIP20+
           IPIP21+IPIP22+IPIP23+IPIP24+IPIP25+IPIP26+IPIP27+IPIP28+IPIP29+IPIP30+
           IPIP31+IPIP32+IPIP33+IPIP34+IPIP35+IPIP36+IPIP37+IPIP38+IPIP39+IPIP40+
           IPIP41+IPIP42+IPIP43+IPIP44+IPIP45+IPIP46+IPIP47+IPIP48+IPIP49+IPIP50+
           IPIP51+IPIP52+IPIP53+IPIP54+IPIP55+IPIP56+IPIP57+IPIP58+IPIP59+IPIP60+
           IPIP61+IPIP62+IPIP63+IPIP64+IPIP65+IPIP66+IPIP67+IPIP68+IPIP69+IPIP70+
           IPIP71+IPIP72+IPIP73+IPIP74+IPIP75+IPIP76+IPIP77+IPIP78+IPIP79+IPIP80+
           IPIP81+IPIP82+IPIP83+IPIP84+IPIP85+IPIP86+IPIP87+IPIP88+IPIP89+IPIP90+
           IPIP91+IPIP92+IPIP93+IPIP94+IPIP95+IPIP96+IPIP97+IPIP98+IPIP99+IPIP100
'
CFA_model3 <- cfa(BIFmodel, data = mydata,estimator="MLR",orthogonal=TRUE,std.lv=TRUE)
options(max.print=9000)
summary(CFA_model3, standardized = TRUE, fit.measures = TRUE)

#Comparative Fit Index (CFI)                    0.566       0.566
#Tucker-Lewis Index (TLI)                       0.548       0.547
#Robust Comparative Fit Index (CFI)                         0.573
#Robust Tucker-Lewis Index (TLI)                            0.555
#RMSEA                                          0.062       0.059
#Robust RMSEA                                               0.061
#SRMR                                           0.109       0.109

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# Median %variance = .271

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    IPIP1             0.357    0.074    4.816    0.000    0.357    0.307
#    IPIP6             0.592    0.079    7.507    0.000    0.592    0.476
#    IPIP11            0.936    0.052   17.907    0.000    0.936    0.742
#    IPIP16            0.624    0.121    5.139    0.000    0.624    0.426
#    IPIP21            0.535    0.086    6.213    0.000    0.535    0.383
#    IPIP26            0.899    0.083   10.870    0.000    0.899    0.679
#    IPIP31            0.897    0.083   10.791    0.000    0.897    0.681
#    IPIP36            0.603    0.117    5.175    0.000    0.603    0.435
#    IPIP41            0.495    0.070    7.111    0.000    0.495    0.417
#    IPIP46            0.545    0.060    9.149    0.000    0.545    0.468
#    IPIP51            0.564    0.097    5.803    0.000    0.564    0.476
#    IPIP56            0.383    0.115    3.321    0.001    0.383    0.282
#    IPIP61            0.534    0.097    5.518    0.000    0.534    0.368
#    IPIP66            0.857    0.068   12.678    0.000    0.857    0.685
#    IPIP71            1.045    0.054   19.526    0.000    1.045    0.803
#    IPIP76            0.459    0.069    6.610    0.000    0.459    0.413
#    IPIP81            0.850    0.089    9.552    0.000    0.850    0.663
#    IPIP86            0.827    0.077   10.697    0.000    0.827    0.670
#    IPIP91            0.567    0.075    7.597    0.000    0.567    0.485
#    IPIP96            0.835    0.082   10.234    0.000    0.835    0.694
#  factor2 =~                                                            
#    IPIP2             0.729    0.054   13.608    0.000    0.729    0.712
#    IPIP7             0.123    0.141    0.873    0.383    0.123    0.079
#    IPIP12            0.353    0.159    2.221    0.026    0.353    0.236
#    IPIP17            0.355    0.101    3.527    0.000    0.355    0.315
#    IPIP22            0.832    0.054   15.449    0.000    0.832    0.749
#    IPIP27            0.391    0.151    2.589    0.010    0.391    0.279
#    IPIP32            0.755    0.173    4.353    0.000    0.755    0.487
#    IPIP37            0.308    0.101    3.064    0.002    0.308    0.273
#    IPIP42            0.594    0.050   11.943    0.000    0.594    0.691
#    IPIP47            0.172    0.144    1.196    0.232    0.172    0.123
#    IPIP52            0.751    0.167    4.505    0.000    0.751    0.482
#    IPIP57            0.204    0.058    3.495    0.000    0.204    0.206
#    IPIP62            0.910    0.044   20.571    0.000    0.910    0.827
#    IPIP67            0.282    0.168    1.682    0.093    0.282    0.195
#    IPIP72            0.625    0.168    3.730    0.000    0.625    0.410
#    IPIP77            0.197    0.069    2.842    0.004    0.197    0.201
#    IPIP82            0.604    0.060   10.113    0.000    0.604    0.667
#    IPIP87            0.175    0.122    1.437    0.151    0.175    0.124
#    IPIP92            0.549    0.158    3.471    0.001    0.549    0.344
#    IPIP97            0.397    0.136    2.927    0.003    0.397    0.265
#  factor3 =~                                                            
#    IPIP3             0.676    0.066   10.202    0.000    0.676    0.538
#    IPIP8             0.404    0.071    5.721    0.000    0.404    0.357
#    IPIP13            0.579    0.052   11.069    0.000    0.579    0.548
#    IPIP18            0.637    0.096    6.646    0.000    0.637    0.415
#    IPIP23            0.641    0.073    8.788    0.000    0.641    0.510
#    IPIP28            0.786    0.060   12.999    0.000    0.786    0.643
#    IPIP33            0.811    0.053   15.170    0.000    0.811    0.715
#    IPIP38            0.126    0.073    1.732    0.083    0.126    0.107
#    IPIP43            0.297    0.079    3.740    0.000    0.297    0.239
#    IPIP48            0.310    0.062    5.006    0.000    0.310    0.294
#    IPIP53            0.711    0.052   13.591    0.000    0.711    0.699
#    IPIP58            0.465    0.093    4.997    0.000    0.465    0.319
#    IPIP63            0.469    0.092    5.116    0.000    0.469    0.380
#    IPIP68            0.495    0.071    6.922    0.000    0.495    0.464
#    IPIP73            0.702    0.050   14.096    0.000    0.702    0.681
#    IPIP78            0.462    0.100    4.639    0.000    0.462    0.317
#    IPIP83            0.459    0.098    4.685    0.000    0.459    0.353
#    IPIP88            0.576    0.060    9.672    0.000    0.576    0.543
#    IPIP93           -0.487    0.074   -6.599    0.000   -0.487   -0.420
#    IPIP98            0.332    0.057    5.834    0.000    0.332    0.355
#  factor4 =~                                                            
#    IPIP4             0.361    0.069    5.256    0.000    0.361    0.356
#    IPIP9             0.662    0.109    6.059    0.000    0.662    0.463
#    IPIP14            0.379    0.115    3.282    0.001    0.379    0.298
#    IPIP19            0.690    0.053   12.986    0.000    0.690    0.657
#    IPIP24            0.431    0.089    4.838    0.000    0.431    0.352
#    IPIP29            0.752    0.072   10.398    0.000    0.752    0.532
#    IPIP34            0.305    0.140    2.178    0.029    0.305    0.227
#    IPIP39            0.737    0.062   11.856    0.000    0.737    0.660
#    IPIP44            0.391    0.067    5.832    0.000    0.391    0.385
#    IPIP49            0.712    0.070   10.181    0.000    0.712    0.545
#    IPIP54            0.309    0.146    2.124    0.034    0.309    0.226
#    IPIP59            0.870    0.084   10.302    0.000    0.870    0.720
#    IPIP64            0.201    0.079    2.538    0.011    0.201    0.181
#    IPIP69            0.714    0.108    6.599    0.000    0.714    0.474
#    IPIP74            0.608    0.141    4.318    0.000    0.608    0.412
#    IPIP79            0.611    0.054   11.257    0.000    0.611    0.628
#    IPIP84            0.286    0.065    4.368    0.000    0.286    0.280
#    IPIP89            0.864    0.093    9.268    0.000    0.864    0.705
#    IPIP94            0.405    0.139    2.919    0.004    0.405    0.290
#    IPIP99            0.669    0.070    9.564    0.000    0.669    0.607
#  factor5 =~                                                            
#    IPIP5             0.637    0.119    5.342    0.000    0.637    0.440
#    IPIP10            0.117    0.134    0.873    0.382    0.117    0.115
#    IPIP15            0.456    0.067    6.794    0.000    0.456    0.450
#    IPIP20            0.398    0.140    2.849    0.004    0.398    0.256
#    IPIP25            0.689    0.106    6.523    0.000    0.689    0.449
#    IPIP30            0.328    0.146    2.244    0.025    0.328    0.274
#    IPIP35            0.564    0.076    7.408    0.000    0.564    0.575
#    IPIP40            0.296    0.117    2.528    0.011    0.296    0.190
#    IPIP45            0.549    0.111    4.953    0.000    0.549    0.361
#    IPIP50            0.366    0.130    2.810    0.005    0.366    0.347
#    IPIP55            0.574    0.082    6.988    0.000    0.574    0.595
#    IPIP60            0.139    0.133    1.045    0.296    0.139    0.091
#    IPIP65            0.567    0.097    5.835    0.000    0.567    0.371
#    IPIP70            0.288    0.130    2.224    0.026    0.288    0.282
#    IPIP75            0.631    0.057   11.153    0.000    0.631    0.703
#    IPIP80            0.462    0.105    4.379    0.000    0.462    0.301
#    IPIP85            0.579    0.048   12.056    0.000    0.579    0.601
#    IPIP90            0.391    0.150    2.599    0.009    0.391    0.362
#    IPIP95            0.520    0.074    7.055    0.000    0.520    0.458
#    IPIP100           0.463    0.148    3.115    0.002    0.463    0.366
#  global =~                                                             
#    IPIP1             0.336    0.063    5.334    0.000    0.336    0.290
#    IPIP2             0.101    0.133    0.759    0.448    0.101    0.099
#    IPIP3             0.137    0.112    1.225    0.221    0.137    0.109
#    IPIP4            -0.149    0.073   -2.032    0.042   -0.149   -0.147
#    IPIP5             0.488    0.196    2.489    0.013    0.488    0.337
#    IPIP6             0.127    0.227    0.557    0.577    0.127    0.102
#    IPIP7             0.473    0.125    3.775    0.000    0.473    0.305
#    IPIP8            -0.024    0.072   -0.336    0.737   -0.024   -0.021
#    IPIP9             0.528    0.165    3.204    0.001    0.528    0.369
#    IPIP10            0.024    0.153    0.155    0.877    0.024    0.023
#    IPIP11            0.044    0.131    0.333    0.739    0.044    0.035
#    IPIP12            0.863    0.077   11.232    0.000    0.863    0.576
#    IPIP13           -0.047    0.145   -0.323    0.747   -0.047   -0.044
#    IPIP14            0.418    0.096    4.369    0.000    0.418    0.330
#    IPIP15           -0.085    0.096   -0.882    0.378   -0.085   -0.084
#    IPIP16            0.632    0.123    5.137    0.000    0.632    0.432
#    IPIP17            0.088    0.172    0.510    0.610    0.088    0.078
#    IPIP18            0.538    0.083    6.493    0.000    0.538    0.351
#    IPIP19           -0.216    0.125   -1.734    0.083   -0.216   -0.206
#    IPIP20            0.588    0.123    4.776    0.000    0.588    0.377
#    IPIP21            0.492    0.110    4.491    0.000    0.492    0.353
#    IPIP22            0.115    0.174    0.661    0.509    0.115    0.104
#    IPIP23            0.298    0.114    2.618    0.009    0.298    0.237
#    IPIP24           -0.071    0.088   -0.803    0.422   -0.071   -0.058
#    IPIP25            0.602    0.148    4.074    0.000    0.602    0.392
#    IPIP26            0.067    0.267    0.251    0.801    0.067    0.051
#    IPIP27            0.722    0.074    9.746    0.000    0.722    0.516
#    IPIP28           -0.073    0.070   -1.043    0.297   -0.073   -0.060
#    IPIP29            0.190    0.239    0.797    0.426    0.190    0.135
#    IPIP30           -0.035    0.151   -0.230    0.818   -0.035   -0.029
#    IPIP31            0.012    0.085    0.137    0.891    0.012    0.009
#    IPIP32            0.858    0.128    6.729    0.000    0.858    0.553
#    IPIP33           -0.075    0.103   -0.727    0.467   -0.075   -0.066
#    IPIP34            0.543    0.090    6.008    0.000    0.543    0.406
#    IPIP35           -0.135    0.116   -1.167    0.243   -0.135   -0.137
#    IPIP36            0.674    0.104    6.466    0.000    0.674    0.486
#    IPIP37            0.114    0.169    0.678    0.498    0.114    0.102
#    IPIP38            0.283    0.062    4.592    0.000    0.283    0.241
#    IPIP39           -0.217    0.175   -1.238    0.216   -0.217   -0.194
#    IPIP40            0.654    0.101    6.500    0.000    0.654    0.421
#    IPIP41            0.341    0.066    5.143    0.000    0.341    0.287
#    IPIP42            0.042    0.104    0.400    0.689    0.042    0.048
#    IPIP43            0.399    0.068    5.873    0.000    0.399    0.321
#    IPIP44           -0.220    0.084   -2.620    0.009   -0.220   -0.217
#    IPIP45            0.559    0.154    3.633    0.000    0.559    0.368
#    IPIP46           -0.022    0.172   -0.128    0.898   -0.022   -0.019
#    IPIP47            0.636    0.084    7.587    0.000    0.636    0.454
#    IPIP48            0.056    0.083    0.682    0.495    0.056    0.054
#    IPIP49            0.155    0.205    0.756    0.449    0.155    0.119
#    IPIP50           -0.007    0.127   -0.054    0.957   -0.007   -0.006
#    IPIP51           -0.202    0.083   -2.446    0.014   -0.202   -0.171
#    IPIP52            0.769    0.164    4.697    0.000    0.769    0.494
#    IPIP53           -0.068    0.122   -0.562    0.574   -0.068   -0.067
#    IPIP54            0.566    0.082    6.875    0.000    0.566    0.413
#    IPIP55           -0.083    0.131   -0.633    0.527   -0.083   -0.086
#    IPIP56            0.594    0.127    4.673    0.000    0.594    0.437
#    IPIP57            0.023    0.057    0.410    0.682    0.023    0.024
#    IPIP58            0.437    0.088    4.968    0.000    0.437    0.299
#    IPIP59           -0.287    0.180   -1.593    0.111   -0.287   -0.237
#    IPIP60            0.591    0.159    3.716    0.000    0.591    0.387
#    IPIP61            0.682    0.089    7.671    0.000    0.682    0.469
#    IPIP62            0.032    0.193    0.165    0.869    0.032    0.029
#    IPIP63            0.311    0.190    1.636    0.102    0.311    0.252
#    IPIP64           -0.076    0.121   -0.631    0.528   -0.076   -0.069
#    IPIP65            0.404    0.114    3.534    0.000    0.404    0.265
#    IPIP66            0.058    0.219    0.264    0.792    0.058    0.046
#    IPIP67            0.620    0.114    5.459    0.000    0.620    0.429
#    IPIP68            0.015    0.098    0.158    0.874    0.015    0.015
#    IPIP69            0.493    0.182    2.710    0.007    0.493    0.327
#    IPIP70           -0.001    0.166   -0.008    0.994   -0.001   -0.001
#    IPIP71           -0.075    0.147   -0.511    0.609   -0.075   -0.058
#    IPIP72            0.869    0.094    9.212    0.000    0.869    0.570
#    IPIP73           -0.081    0.129   -0.631    0.528   -0.081   -0.079
#    IPIP74            0.501    0.125    3.995    0.000    0.501    0.339
#    IPIP75           -0.061    0.096   -0.629    0.529   -0.061   -0.068
#    IPIP76            0.233    0.110    2.116    0.034    0.233    0.209
#    IPIP77            0.009    0.140    0.061    0.951    0.009    0.009
#    IPIP78            0.515    0.079    6.524    0.000    0.515    0.354
#    IPIP79           -0.197    0.104   -1.892    0.059   -0.197   -0.203
#    IPIP80            0.691    0.101    6.838    0.000    0.691    0.450
#    IPIP81           -0.018    0.084   -0.218    0.828   -0.018   -0.014
#    IPIP82            0.055    0.151    0.367    0.713    0.055    0.061
#    IPIP83            0.460    0.163    2.828    0.005    0.460    0.353
#    IPIP84           -0.088    0.089   -0.988    0.323   -0.088   -0.087
#    IPIP85            0.034    0.074    0.453    0.651    0.034    0.035
#    IPIP86            0.098    0.219    0.448    0.654    0.098    0.079
#    IPIP87            0.472    0.122    3.877    0.000    0.472    0.334
#    IPIP88            0.058    0.068    0.852    0.394    0.058    0.054
#    IPIP89           -0.351    0.188   -1.862    0.063   -0.351   -0.286
#    IPIP90           -0.004    0.148   -0.027    0.978   -0.004   -0.004
#    IPIP91            0.104    0.218    0.478    0.632    0.104    0.089
#    IPIP92            0.758    0.168    4.514    0.000    0.758    0.475
#    IPIP93           -0.069    0.197   -0.348    0.728   -0.069   -0.059
#    IPIP94            0.546    0.088    6.172    0.000    0.546    0.390
#    IPIP95           -0.053    0.139   -0.384    0.701   -0.053   -0.047
#    IPIP96            0.122    0.231    0.528    0.597    0.122    0.101
#    IPIP97            0.577    0.108    5.327    0.000    0.577    0.385
#    IPIP98           -0.040    0.059   -0.677    0.498   -0.040   -0.043
#    IPIP99           -0.274    0.165   -1.659    0.097   -0.274   -0.249
#    IPIP100          -0.084    0.107   -0.792    0.428   -0.084   -0.067

semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="global")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model3)
#$ModelLevelIndices
#   ECV.global           PUC  Omega.global OmegaH.global 
#    0.2401981     0.8080808     0.8921347     0.3200880 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8294699 0.2043087 0.1705301 0.9053745 0.8286144 0.9148599 0.9565734
#factor2 0.6147299 0.1342071 0.3852701 0.8699558 0.5641460 0.8852621 0.9392864
#factor3 0.8499943 0.1510483 0.1500057 0.8254441 0.7681472 0.8728739 0.9352933
#factor4 0.7683471 0.1562478 0.2316529 0.8538467 0.8458335 0.8793256 0.9429650
#factor5 0.7387172 0.1139899 0.2612828 0.8055312 0.7251753 0.8259547 0.9112884
#global  0.2401981 0.2401981 0.2401981 0.8921347 0.3200880 0.8966528 0.9488303




