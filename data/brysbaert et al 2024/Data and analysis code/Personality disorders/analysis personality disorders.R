# Analysis Personality disorders


################################################################
### ICD-11 Carnovale et al. (2020) https://psycnet.apa.org/record/2019-57538-001
### https://osf.io/nfe6q

ICD11_Carnovale <- read.csv("ICD11_Carnovale.csv")
colnames(ICD11_Carnovale)
mydata <- as.data.frame(ICD11_Carnovale[,6:65])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 7 components
# Eigenvalue 1 = 6.51; eigenvalue 2 = 4.38

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 components
# Eigenvalue 1 = 7.63; eigenvalue 2 = 5.03

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.11, RMSEA=.081, RMSR=.11, TLI=.48

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.13, RMSEA=.102, RMSR=.13, TLI=.222

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities and response bias

# Give solution with 5 factors
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.31, RMSEA=.046, RMSR=.04, TLI=.764
#      MR1   MR2   MR3   MR5   MR4
#MR1  1.00 -0.13  0.28  0.05 -0.02
#MR2 -0.13  1.00  0.09 -0.25  0.26
#MR3  0.28  0.09  1.00  0.12 -0.12
#MR5  0.05 -0.25  0.12  1.00 -0.01
#MR4 -0.02  0.26 -0.12 -0.01  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.31, RMSEA=.046, RMSR=.04, TLI=.764
#      MR1   MR2   MR3   MR5   MR4
#MR1  1.00 -0.13  0.28  0.05 -0.02
#MR2 -0.13  1.00  0.09 -0.25  0.26
#MR3  0.28  0.09  1.00  0.12 -0.12
#MR5  0.05 -0.25  0.12  1.00 -0.01
#MR4 -0.02  0.26 -0.12 -0.01  1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19+p20+
            p21+p22+p23+p24+p25+p26+p27+p28+p29+p30+p31+p32+p33+p34+p35+p36+p37+p38+p39+p40+
            p41+p42+p43+p44+p45+p46+p47+p48+p49+p50+p51+p52+p53+p54+p55+p56+p57+p58+p59+p60
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.535       0.377
#Tucker-Lewis Index (TLI)                       0.518       0.355
#Robust Comparative Fit Index (CFI)                         0.263
#Robust Tucker-Lewis Index (TLI)                            0.237
#RMSEA                                          0.140       0.087
#Robust RMSEA                                               0.103
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .122

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.292       0.280
#Tucker-Lewis Index (TLI)                       0.268       0.255
#Robust Comparative Fit Index (CFI)                         0.295
#Robust Tucker-Lewis Index (TLI)                            0.270
#RMSEA                                          0.084       0.078
#Robust RMSEA                                               0.082
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .078

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 5 factors based on theoretical analysis
EGAmodel= '
  factor1 =~  p1+p6+p11+p16+p21+p26+p31+p36+p41+p46+p51+p56
  factor2 =~  p2+p7+p12+p17+p22+p27+p32+p37+p42+p47+p52+p57
  factor3 =~  p3+p8+p13+p18+p23+p28+p33+p38+p43+p48+p53+p58
  factor4 =~  p4+p9+p14+p19+p24+p29+p34+p39+p44+p49+p54+p59
  factor5 =~  p5+p10+p15+p20+p25+p30+p35+p40+p45+p50+p55+p60
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.805       0.724
#Tucker-Lewis Index (TLI)                       0.797       0.712
#Robust Comparative Fit Index (CFI)                         0.573
#Robust Tucker-Lewis Index (TLI)                            0.555
#RMSEA                                          0.091       0.058
#Robust RMSEA                                               0.079
#SRMR                                           0.095       0.095

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.392    0.042    9.364    0.000    0.392    0.392
#    factor3           0.397    0.041    9.750    0.000    0.397    0.397
#    factor4           0.049    0.053    0.926    0.354    0.049    0.049
#    factor5           0.124    0.048    2.567    0.010    0.124    0.124
#  factor2 ~~                                                            
#    factor3           0.277    0.043    6.449    0.000    0.277    0.277
#    factor4           0.456    0.041   11.118    0.000    0.456    0.456
#    factor5          -0.553    0.035  -15.716    0.000   -0.553   -0.553
#  factor3 ~~                                                            
#    factor4           0.152    0.047    3.240    0.001    0.152    0.152
#    factor5           0.151    0.045    3.369    0.001    0.151    0.151
#  factor4 ~~                                                            
#    factor5          -0.082    0.049   -1.689    0.091   -0.082   -0.082

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.639       0.650
#Tucker-Lewis Index (TLI)                       0.624       0.635
#Robust Comparative Fit Index (CFI)                         0.653
#Robust Tucker-Lewis Index (TLI)                            0.639
#RMSEA                                          0.060       0.055
#Robust RMSEA                                               0.058
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .217

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.394    0.056    7.100    0.000    0.394    0.394
#    factor3           0.390    0.055    7.042    0.000    0.390    0.390
#    factor4           0.028    0.071    0.388    0.698    0.028    0.028
#    factor5           0.115    0.070    1.656    0.098    0.115    0.115
#  factor2 ~~                                                            
#    factor3           0.234    0.068    3.432    0.001    0.234    0.234
#    factor4           0.418    0.058    7.201    0.000    0.418    0.418
#    factor5          -0.534    0.059   -9.014    0.000   -0.534   -0.534
#  factor3 ~~                                                            
#    factor4           0.078    0.073    1.077    0.282    0.078    0.078
#    factor5           0.207    0.064    3.217    0.001    0.207    0.207
#  factor4 ~~                                                            
#    factor5          -0.077    0.066   -1.169    0.242   -0.077   -0.077

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.597       0.607
#Tucker-Lewis Index (TLI)                       0.583       0.593
#Robust Comparative Fit Index (CFI)                         0.610
#Robust Tucker-Lewis Index (TLI)                            0.597
#RMSEA                                          0.063       0.058
#Robust RMSEA                                               0.061
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .213

# Bifactor model
BIFmodel= '
  factor1 =~  p1+p6+p11+p16+p21+p26+p31+p36+p41+p46+p51+p56
  factor2 =~  p2+p7+p12+p17+p22+p27+p32+p37+p42+p47+p52+p57
  factor3 =~  p3+p8+p13+p18+p23+p28+p33+p38+p43+p48+p53+p58
  factor4 =~  p4+p9+p14+p19+p24+p29+p34+p39+p44+p49+p54+p59
  factor5 =~  p5+p10+p15+p20+p25+p30+p35+p40+p45+p50+p55+p60
 general=~  p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19+p20+
            p21+p22+p23+p24+p25+p26+p27+p28+p29+p30+p31+p32+p33+p34+p35+p36+p37+p38+p39+p40+
            p41+p42+p43+p44+p45+p46+p47+p48+p49+p50+p51+p52+p53+p54+p55+p56+p57+p58+p59+p60
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.690       0.696
#Tucker-Lewis Index (TLI)                       0.668       0.674
#Robust Comparative Fit Index (CFI)                         0.702
#Robust Tucker-Lewis Index (TLI)                            0.681
#RMSEA                                          0.056       0.052
#Robust RMSEA                                               0.054
#SRMR                                           0.095       0.095

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .246

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.19078943     0.81355932     0.79172640     0.07497858 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#factor1 0.9398190 0.2089702 0.06018101 0.8340401 0.78326358 0.8582768 0.9243680
#factor2 0.6937852 0.1655953 0.30621476 0.8349354 0.65258586 0.8262689 0.9048229
#factor3 0.9347188 0.1676449 0.06528122 0.7781904 0.76546150 0.8192943 0.9071707
#factor4 0.8272665 0.1251800 0.17273350 0.7357626 0.61358244 0.7767696 0.8797953
#factor5 0.6808661 0.1418201 0.31913387 0.7741352 0.66921298 0.7882977 0.8995563
#general 0.1907894 0.1907894 0.19078943 0.7917264 0.07497858 0.8014199 0.8945739








################################################################
### Personality functioning LPFS-BF 2.0 https://psycnet.apa.org/record/2023-00705-001
### https://osf.io/dmsxr

Natoli_Pers_Functioning <- read.csv("Natoli_Pers_Functioning.csv")
colnames(Natoli_Pers_Functioning)
mydata <- as.data.frame(Natoli_Pers_Functioning[,6:17])
colnames(mydata)
mydata[mydata == -999] <- NA 
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors 2 components
# Eigenvalue 1 = 5.32; eigenvalue 2 = 0.54

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 6.38; eigenvalue 2 = .58

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.44, RMSEA=.094, RMSR=.06, TLI=.89

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.12, RMSR=.06, TLI=.873

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communities and response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.057, RMSR=.02, TLI=.96
#     MR1  MR2
#MR1 1.00 0.73
#MR2 0.73 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.076, RMSR=.03, TLI=.949
#     MR1  MR2
#MR1 1.00 0.76
#MR2 0.76 1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~   lpfsbf1+lpfsbf2+lpfsbf3+lpfsbf4+lpfsbf5+lpfsbf6+lpfsbf7+
             lpfsbf8+lpfsbf9+lpfsbf10+lpfsbf11+lpfsbf12
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.950
#Tucker-Lewis Index (TLI)                       0.987       0.939
#Robust Comparative Fit Index (CFI)                         0.898
#Robust Tucker-Lewis Index (TLI)                            0.875
#RMSEA                                          0.081       0.112
#Robust RMSEA                                               0.119
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .519

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
  self      =~   lpfsbf1+lpfsbf2+lpfsbf3+lpfsbf4+lpfsbf5+lpfsbf6
  interpers =~   lpfsbf7+lpfsbf8+lpfsbf9+lpfsbf10+lpfsbf11+lpfsbf12
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.997       0.983
#Tucker-Lewis Index (TLI)                       0.996       0.979
#Robust Comparative Fit Index (CFI)                         0.957
#Robust Tucker-Lewis Index (TLI)                            0.947
#RMSEA                                          0.045       0.065
#Robust RMSEA                                               0.078
#SRMR                                           0.033       0.033

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#self ~~                                                               
#  interpers         0.824    0.006  129.971    0.000    0.824    0.824

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .573

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.671       0.538
#Tucker-Lewis Index (TLI)                       0.598       0.436
#Robust Comparative Fit Index (CFI)                         0.855
#Robust Tucker-Lewis Index (TLI)                            0.823
#RMSEA                                          0.446       0.339
#Robust RMSEA                                               0.142
#SRMR                                           0.327       0.327

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .575

# Bifactor model
BIFmodel= '
  self      =~   lpfsbf1+lpfsbf2+lpfsbf3+lpfsbf4+lpfsbf5+lpfsbf6
  interpers =~   lpfsbf7+lpfsbf8+lpfsbf9+lpfsbf10+lpfsbf11+lpfsbf12
 general=~   lpfsbf1+lpfsbf2+lpfsbf3+lpfsbf4+lpfsbf5+lpfsbf6+lpfsbf7+
             lpfsbf8+lpfsbf9+lpfsbf10+lpfsbf11+lpfsbf12
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.992
#Tucker-Lewis Index (TLI)                       0.998       0.988
#Robust Comparative Fit Index (CFI)                         0.976
#Robust Tucker-Lewis Index (TLI)                            0.962
#RMSEA                                          0.030       0.050
#Robust RMSEA                                               0.065
#SRMR                                           0.020       0.020

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .596

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  self =~                                                               
#    lpfsbf1           0.489    0.013   38.706    0.000    0.489    0.489
#    lpfsbf2           0.530    0.012   43.188    0.000    0.530    0.530
#    lpfsbf3           0.442    0.013   34.414    0.000    0.442    0.442
#    lpfsbf4           0.418    0.014   30.945    0.000    0.418    0.418
#    lpfsbf5           0.479    0.012   38.916    0.000    0.479    0.479
#    lpfsbf6           0.396    0.013   29.494    0.000    0.396    0.396
#  interpers =~                                                          
#    lpfsbf7           0.389    0.025   15.497    0.000    0.389    0.389
#    lpfsbf8           0.239    0.021   11.438    0.000    0.239    0.239
#    lpfsbf9           0.296    0.021   14.149    0.000    0.296    0.296
#    lpfsbf10         -0.022    0.019   -1.118    0.264   -0.022   -0.022
#    lpfsbf11         -0.312    0.038   -8.209    0.000   -0.312   -0.312
#    lpfsbf12          0.122    0.020    6.193    0.000    0.122    0.122
#  general =~                                                            
#    lpfsbf1           0.690    0.010   68.656    0.000    0.690    0.690
#    lpfsbf2           0.661    0.010   67.612    0.000    0.661    0.661
#    lpfsbf3           0.695    0.009   74.304    0.000    0.695    0.695
#    lpfsbf4           0.642    0.010   61.842    0.000    0.642    0.642
#    lpfsbf5           0.731    0.009   82.490    0.000    0.731    0.731
#    lpfsbf6           0.641    0.010   64.323    0.000    0.641    0.641
#    lpfsbf7           0.653    0.012   55.037    0.000    0.653    0.653
#    lpfsbf8           0.627    0.011   55.637    0.000    0.627    0.627
#    lpfsbf9           0.719    0.010   71.714    0.000    0.719    0.719
#    lpfsbf10          0.715    0.010   75.021    0.000    0.715    0.715
#    lpfsbf11          0.797    0.009   86.678    0.000    0.797    0.797
#    lpfsbf12          0.740    0.009   83.487    0.000    0.740    0.740

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.7742587      0.5454545      0.9245614      0.8111891 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#self      0.3168207 0.17097812 0.6831793 0.8947899 0.26833290 0.6218118 0.8208032
#interpers 0.1189647 0.05476322 0.8810353 0.8388401 0.02244486 0.3143993 0.7142434
#general   0.7742587 0.77425865 0.7742587 0.9245614 0.81118911 0.9207476 0.9500379








################################################################
### Self-Report Psychopathy Scale (SRP4; Paulhus et al., 2017)
### Knack et al. (2021) https://osf.io/5ce7m
### https://osf.io/gau8f/?view_only=

library(haven)
Psychopathy_Knack <- read_sav("Psychopathy_Knack.sav")
colnames(Psychopathy_Knack)
mydata <- as.data.frame(Psychopathy_Knack[,26:89])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 6 components
# Eigenvalue 1 = 12.22; eigenvalue 2 = 2.68

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 components
# Eigenvalue 1 = 16.69; eigenvalue 2 = 3.2

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.065, RMSR=.08, TLI=.54

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.206, RMSR=.1, TLI=.135

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities and response bias

# Give solution with 4 facets (https://arc.psych.wisc.edu/self-report/self-report-psychopathy-srp-iii/)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.049, RMSR=.05, TLI=.743
# With factor correlations of 
#MR1   MR4   MR3   MR2
#MR1  1.00  0.40  0.43 -0.25
#MR4  0.40  1.00  0.36 -0.10
#MR3  0.43  0.36  1.00 -0.11
#MR2 -0.25 -0.10 -0.11  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.207, RMSR=.06, TLI=.124
#      MR1   MR4   MR3   MR2
#MR1  1.00  0.45  0.40 -0.09
#MR4  0.45  1.00  0.45 -0.10
#MR3  0.40  0.45  1.00 -0.03
#MR2 -0.09 -0.10 -0.03  1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~   SRP_B_1+SRP_B_2+SRP_B_3+SRP_B_4+SRP_B_5+SRP_B_6+SRP_B_7+SRP_B_8+SRP_B_9+SRP_B_10+
             SRP_B_11+SRP_B_12+SRP_B_13+SRP_B_14+SRP_B_15+SRP_B_16+SRP_B_17+SRP_B_18+SRP_B_19+
             SRP_B_20+SRP_B_21+SRP_B_22+SRP_B_23+SRP_B_24+SRP_B_25+SRP_B_26+SRP_B_27+SRP_B_28+
             SRP_B_29+SRP_B_30+SRP_B_31+SRP_B_32+SRP_B_33+SRP_B_34+SRP_B_35+SRP_B_36+SRP_B_37+
             SRP_B_38+SRP_B_39+SRP_B_40+SRP_B_41+SRP_B_42+SRP_B_43+SRP_B_44+SRP_B_45+SRP_B_46+
             SRP_B_47+SRP_B_48+SRP_B_49+SRP_B_50+SRP_B_51+SRP_B_52+SRP_B_53+SRP_B_54+SRP_B_55+
             SRP_B_56+SRP_B_57+SRP_B_58+SRP_B_59+SRP_B_60+SRP_B_61+SRP_B_62+SRP_B_63+SRP_B_64
'

library(lavaan)
# because of warnings with ordered, MLR
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.544       0.555
#Tucker-Lewis Index (TLI)                       0.529       0.540
#Robust Comparative Fit Index (CFI)                         0.559
#Robust Tucker-Lewis Index (TLI)                            0.544
#RMSEA                                          0.069       0.064
#Robust RMSEA                                               0.067
#SRMR                                           0.079       0.079

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .188

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
  manip      =~   SRP_B_3+SRP_B_8+SRP_B_13+SRP_B_16+SRP_B_20+SRP_B_24+SRP_B_27+SRP_B_31+SRP_B_35+
                  SRP_B_38+SRP_B_41+SRP_B_45+SRP_B_50+SRP_B_54+SRP_B_58+SRP_B_61
  callous    =~   SRP_B_2+SRP_B_7+SRP_B_11+SRP_B_15+SRP_B_19+SRP_B_23+SRP_B_26+SRP_B_30+SRP_B_33+
                  SRP_B_37+SRP_B_40+SRP_B_44+SRP_B_48+SRP_B_53+SRP_B_56+SRP_B_60
  erratic    =~   SRP_B_1+SRP_B_4+SRP_B_9+SRP_B_14+SRP_B_17+SRP_B_22+SRP_B_25+SRP_B_28+SRP_B_32+
                  SRP_B_36+SRP_B_39+SRP_B_42+SRP_B_47+SRP_B_51+SRP_B_55+SRP_B_59
  antiscoc   =~   SRP_B_5+SRP_B_6+SRP_B_10+SRP_B_12+SRP_B_18+SRP_B_21+SRP_B_29+SRP_B_34+SRP_B_43+
                  SRP_B_46+SRP_B_49+SRP_B_52+SRP_B_57+SRP_B_62+SRP_B_63+SRP_B_64
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.616       0.627
#Tucker-Lewis Index (TLI)                       0.602       0.613
#Robust Comparative Fit Index (CFI)                         0.632
#Robust Tucker-Lewis Index (TLI)                            0.618
#RMSEA                                          0.064       0.058
#Robust RMSEA                                               0.061
#SRMR                                           0.078       0.078

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  manip ~~                                                              
#    callous           0.869    0.037   23.260    0.000    0.869    0.869
#    erratic           0.765    0.037   20.944    0.000    0.765    0.765
#    antiscoc         -0.670    0.051  -13.074    0.000   -0.670   -0.670
#  callous ~~                                                            
#    erratic           0.696    0.056   12.330    0.000    0.696    0.696
#    antiscoc         -0.610    0.069   -8.802    0.000   -0.610   -0.610
#  erratic ~~                                                            
#    antiscoc         -0.741    0.041  -18.035    0.000   -0.741   -0.741

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .235

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.526       0.532
#Tucker-Lewis Index (TLI)                       0.511       0.517
#Robust Comparative Fit Index (CFI)                         0.539
#Robust Tucker-Lewis Index (TLI)                            0.524
#RMSEA                                          0.070       0.065
#Robust RMSEA                                               0.068
#SRMR                                           0.168       0.168

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .232

# Bifactor model
BIFmodel= '
  manip      =~   SRP_B_3+SRP_B_8+SRP_B_13+SRP_B_16+SRP_B_20+SRP_B_24+SRP_B_27+SRP_B_31+SRP_B_35+
                  SRP_B_38+SRP_B_41+SRP_B_45+SRP_B_50+SRP_B_54+SRP_B_58+SRP_B_61
  callous    =~   SRP_B_2+SRP_B_7+SRP_B_11+SRP_B_15+SRP_B_19+SRP_B_23+SRP_B_26+SRP_B_30+SRP_B_33+
                  SRP_B_37+SRP_B_40+SRP_B_44+SRP_B_48+SRP_B_53+SRP_B_56+SRP_B_60
  erratic    =~   SRP_B_1+SRP_B_4+SRP_B_9+SRP_B_14+SRP_B_17+SRP_B_22+SRP_B_25+SRP_B_28+SRP_B_32+
                  SRP_B_36+SRP_B_39+SRP_B_42+SRP_B_47+SRP_B_51+SRP_B_55+SRP_B_59
  antiscoc   =~   SRP_B_5+SRP_B_6+SRP_B_10+SRP_B_12+SRP_B_18+SRP_B_21+SRP_B_29+SRP_B_34+SRP_B_43+
                  SRP_B_46+SRP_B_49+SRP_B_52+SRP_B_57+SRP_B_62+SRP_B_63+SRP_B_64
  general=~   SRP_B_1+SRP_B_2+SRP_B_3+SRP_B_4+SRP_B_5+SRP_B_6+SRP_B_7+SRP_B_8+SRP_B_9+SRP_B_10+
              SRP_B_11+SRP_B_12+SRP_B_13+SRP_B_14+SRP_B_15+SRP_B_16+SRP_B_17+SRP_B_18+SRP_B_19+
              SRP_B_20+SRP_B_21+SRP_B_22+SRP_B_23+SRP_B_24+SRP_B_25+SRP_B_26+SRP_B_27+SRP_B_28+
              SRP_B_29+SRP_B_30+SRP_B_31+SRP_B_32+SRP_B_33+SRP_B_34+SRP_B_35+SRP_B_36+SRP_B_37+
              SRP_B_38+SRP_B_39+SRP_B_40+SRP_B_41+SRP_B_42+SRP_B_43+SRP_B_44+SRP_B_45+SRP_B_46+
              SRP_B_47+SRP_B_48+SRP_B_49+SRP_B_50+SRP_B_51+SRP_B_52+SRP_B_53+SRP_B_54+SRP_B_55+
              SRP_B_56+SRP_B_57+SRP_B_58+SRP_B_59+SRP_B_60+SRP_B_61+SRP_B_62+SRP_B_63+SRP_B_64
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.683       0.681
#Tucker-Lewis Index (TLI)                       0.662       0.659
#Robust Comparative Fit Index (CFI)                         0.693
#Robust Tucker-Lewis Index (TLI)                            0.672
#RMSEA                                          0.059       0.055
#Robust RMSEA                                               0.057
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .269

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  manip =~                                                              
#    SRP_B_3           0.005    0.090    0.051    0.959    0.005    0.004
#    SRP_B_8           0.465    0.132    3.533    0.000    0.465    0.433
#    SRP_B_13          0.089    0.087    1.025    0.305    0.089    0.084
#    SRP_B_16         -0.047    0.126   -0.374    0.709   -0.047   -0.039
#    SRP_B_20          0.055    0.131    0.418    0.676    0.055    0.062
#    SRP_B_24         -0.086    0.156   -0.555    0.579   -0.086   -0.077
#    SRP_B_27          0.082    0.155    0.527    0.598    0.082    0.082
#    SRP_B_31         -0.115    0.198   -0.581    0.561   -0.115   -0.101
#    SRP_B_35          0.252    0.150    1.685    0.092    0.252    0.254
#    SRP_B_38         -0.060    0.143   -0.417    0.677   -0.060   -0.053
#    SRP_B_41          0.556    0.109    5.079    0.000    0.556    0.471
#    SRP_B_45          0.283    0.105    2.691    0.007    0.283    0.256
#    SRP_B_50          0.255    0.091    2.816    0.005    0.255    0.232
#    SRP_B_54          0.653    0.132    4.962    0.000    0.653    0.590
#    SRP_B_58          0.434    0.104    4.169    0.000    0.434    0.361
#    SRP_B_61         -0.132    0.178   -0.737    0.461   -0.132   -0.112
#  callous =~                                                            
#    SRP_B_2           0.036    0.100    0.362    0.718    0.036    0.033
#    SRP_B_7           0.180    0.094    1.927    0.054    0.180    0.193
#    SRP_B_11         -0.429    0.103   -4.186    0.000   -0.429   -0.401
#    SRP_B_15          0.057    0.126    0.451    0.652    0.057    0.048
#    SRP_B_19         -0.634    0.080   -7.899    0.000   -0.634   -0.673
#    SRP_B_23          0.091    0.133    0.688    0.492    0.091    0.061
#    SRP_B_26         -0.465    0.098   -4.726    0.000   -0.465   -0.473
#    SRP_B_30          0.156    0.101    1.539    0.124    0.156    0.129
#    SRP_B_33          0.224    0.106    2.108    0.035    0.224    0.180
#    SRP_B_37          0.437    0.089    4.925    0.000    0.437    0.408
#    SRP_B_40          0.092    0.146    0.628    0.530    0.092    0.068
#    SRP_B_44         -0.679    0.086   -7.860    0.000   -0.679   -0.696
#    SRP_B_48          0.132    0.089    1.492    0.136    0.132    0.124
#    SRP_B_53          0.248    0.124    2.006    0.045    0.248    0.229
#    SRP_B_56          0.248    0.100    2.475    0.013    0.248    0.276
#    SRP_B_60          0.043    0.106    0.404    0.686    0.043    0.039
#  erratic =~                                                            
#    SRP_B_1           0.180    0.133    1.354    0.176    0.180    0.145
#    SRP_B_4           0.242    0.216    1.123    0.261    0.242    0.155
#    SRP_B_9           0.632    0.142    4.460    0.000    0.632    0.485
#    SRP_B_14         -0.061    0.131   -0.469    0.639   -0.061   -0.055
#    SRP_B_17          0.641    0.095    6.761    0.000    0.641    0.512
#    SRP_B_22         -0.014    0.145   -0.099    0.921   -0.014   -0.012
#    SRP_B_25         -0.395    0.101   -3.915    0.000   -0.395   -0.304
#    SRP_B_28          0.672    0.092    7.335    0.000    0.672    0.555
#    SRP_B_32          0.048    0.071    0.676    0.499    0.048    0.053
#    SRP_B_36         -0.469    0.102   -4.595    0.000   -0.469   -0.326
#    SRP_B_39          0.120    0.121    0.992    0.321    0.120    0.093
#    SRP_B_42          0.093    0.160    0.580    0.562    0.093    0.084
#    SRP_B_47         -0.649    0.081   -8.018    0.000   -0.649   -0.535
#    SRP_B_51          0.009    0.120    0.074    0.941    0.009    0.009
#    SRP_B_55         -0.157    0.131   -1.197    0.231   -0.157   -0.128
#    SRP_B_59         -0.043    0.154   -0.280    0.779   -0.043   -0.037
#  antiscoc =~                                                           
#    SRP_B_5           0.631    0.444    1.423    0.155    0.631    0.506
#    SRP_B_6           0.556    0.404    1.377    0.169    0.556    0.515
#    SRP_B_10         -0.214    0.094   -2.265    0.024   -0.214   -0.187
#    SRP_B_12         -0.196    0.154   -1.272    0.203   -0.196   -0.318
#    SRP_B_18          0.247    0.185    1.332    0.183    0.247    0.251
#    SRP_B_21          0.333    0.095    3.520    0.000    0.333    0.244
#    SRP_B_29         -0.477    0.097   -4.896    0.000   -0.477   -0.441
#    SRP_B_34          0.611    0.248    2.467    0.014    0.611    0.367
#    SRP_B_43         -0.483    0.175   -2.759    0.006   -0.483   -0.342
#    SRP_B_46          0.426    0.232    1.839    0.066    0.426    0.266
#    SRP_B_49         -0.373    0.257   -1.453    0.146   -0.373   -0.415
#    SRP_B_52         -0.100    0.275   -0.365    0.715   -0.100   -0.079
#    SRP_B_57         -0.217    0.105   -2.077    0.038   -0.217   -0.299
#    SRP_B_62         -0.458    0.228   -2.008    0.045   -0.458   -0.333
#    SRP_B_63         -0.184    0.216   -0.851    0.395   -0.184   -0.296
#    SRP_B_64         -0.179    0.225   -0.795    0.427   -0.179   -0.290
#  general =~                                                            
#    SRP_B_1           0.760    0.072   10.492    0.000    0.760    0.613
#    SRP_B_2           0.368    0.069    5.363    0.000    0.368    0.334
#    SRP_B_3           0.548    0.061    9.031    0.000    0.548    0.495
#    SRP_B_4           0.648    0.141    4.610    0.000    0.648    0.414
#    SRP_B_5          -0.351    0.092   -3.833    0.000   -0.351   -0.282
#    SRP_B_6          -0.189    0.082   -2.297    0.022   -0.189   -0.175
#    SRP_B_7           0.499    0.065    7.663    0.000    0.499    0.533
#    SRP_B_8           0.572    0.076    7.500    0.000    0.572    0.533
#    SRP_B_9           0.794    0.104    7.631    0.000    0.794    0.610
#    SRP_B_10          0.755    0.060   12.611    0.000    0.755    0.662
#    SRP_B_11         -0.095    0.089   -1.071    0.284   -0.095   -0.089
#    SRP_B_12          0.202    0.063    3.216    0.001    0.202    0.328
#    SRP_B_13          0.631    0.063    9.969    0.000    0.631    0.595
#    SRP_B_14         -0.185    0.087   -2.125    0.034   -0.185   -0.166
#    SRP_B_15          0.682    0.066   10.303    0.000    0.682    0.576
#    SRP_B_16         -0.551    0.076   -7.226    0.000   -0.551   -0.461
#    SRP_B_17          0.428    0.077    5.583    0.000    0.428    0.342
#    SRP_B_18         -0.131    0.069   -1.905    0.057   -0.131   -0.133
#    SRP_B_19         -0.184    0.088   -2.081    0.037   -0.184   -0.195
#    SRP_B_20          0.484    0.067    7.243    0.000    0.484    0.541
#    SRP_B_21         -0.405    0.084   -4.795    0.000   -0.405   -0.296
#    SRP_B_22         -0.417    0.086   -4.872    0.000   -0.417   -0.358
#    SRP_B_23         -0.323    0.100   -3.244    0.001   -0.323   -0.216
#    SRP_B_24         -0.308    0.090   -3.431    0.001   -0.308   -0.276
#    SRP_B_25         -0.324    0.080   -4.031    0.000   -0.324   -0.250
#    SRP_B_26         -0.190    0.085   -2.245    0.025   -0.190   -0.194
#    SRP_B_27          0.584    0.084    6.959    0.000    0.584    0.586
#    SRP_B_28          0.701    0.085    8.236    0.000    0.701    0.579
#    SRP_B_29          0.496    0.071    6.982    0.000    0.496    0.458
#    SRP_B_30          0.493    0.075    6.551    0.000    0.493    0.407
#    SRP_B_31         -0.419    0.089   -4.704    0.000   -0.419   -0.369
#    SRP_B_32          0.539    0.054   10.044    0.000    0.539    0.593
#    SRP_B_33         -0.020    0.079   -0.258    0.796   -0.020   -0.016
#    SRP_B_34         -0.515    0.148   -3.481    0.000   -0.515   -0.310
#    SRP_B_35          0.555    0.085    6.504    0.000    0.555    0.560
#    SRP_B_36         -0.227    0.095   -2.378    0.017   -0.227   -0.158
#    SRP_B_37          0.545    0.083    6.553    0.000    0.545    0.508
#    SRP_B_38         -0.126    0.080   -1.579    0.114   -0.126   -0.113
#    SRP_B_39          0.644    0.072    8.923    0.000    0.644    0.499
#    SRP_B_40          0.632    0.080    7.922    0.000    0.632    0.469
#    SRP_B_41          0.537    0.074    7.263    0.000    0.537    0.455
#    SRP_B_42          0.538    0.077    6.977    0.000    0.538    0.486
#    SRP_B_43          0.596    0.125    4.777    0.000    0.596    0.422
#    SRP_B_44         -0.113    0.084   -1.352    0.176   -0.113   -0.116
#    SRP_B_45          0.542    0.070    7.802    0.000    0.542    0.491
#    SRP_B_46         -0.428    0.110   -3.903    0.000   -0.428   -0.268
#    SRP_B_47         -0.405    0.090   -4.494    0.000   -0.405   -0.334
#    SRP_B_48          0.364    0.063    5.747    0.000    0.364    0.341
#    SRP_B_49          0.328    0.091    3.603    0.000    0.328    0.366
#    SRP_B_50          0.512    0.056    9.114    0.000    0.512    0.466
#    SRP_B_51          0.583    0.060    9.783    0.000    0.583    0.571
#    SRP_B_52          0.503    0.090    5.571    0.000    0.503    0.396
#    SRP_B_53          0.300    0.084    3.592    0.000    0.300    0.277
#    SRP_B_54          0.517    0.064    8.132    0.000    0.517    0.467
#    SRP_B_55          0.602    0.067    8.920    0.000    0.602    0.488
#    SRP_B_56          0.392    0.085    4.587    0.000    0.392    0.436
#    SRP_B_57          0.390    0.074    5.300    0.000    0.390    0.537
#    SRP_B_58          0.512    0.075    6.818    0.000    0.512    0.426
#    SRP_B_59          0.657    0.068    9.642    0.000    0.657    0.561
#    SRP_B_60          0.668    0.064   10.416    0.000    0.668    0.602
#    SRP_B_61         -0.385    0.089   -4.297    0.000   -0.385   -0.326
#    SRP_B_62          0.712    0.106    6.733    0.000    0.712    0.517
#    SRP_B_63          0.243    0.069    3.512    0.000    0.243    0.392
#    SRP_B_64          0.184    0.068    2.714    0.007    0.184    0.299

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.6551547      0.7619048      0.8359506      0.8097623 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#manip    0.2459624 0.06367697 0.7540376 0.6635542 0.176128992 0.5913324 0.7860552
#callous  0.4338565 0.09829169 0.5661435 0.5307667 0.008192457 0.7302272 0.8569201
#erratic  0.2855474 0.07799133 0.7144526 0.6484930 0.015111370 0.6448512 0.8295300
#antiscoc 0.4344367 0.10488533 0.5655633 0.4395290 0.034613592 0.6892364 0.8261161
#general  0.6551547 0.65515468 0.6551547 0.8359506 0.809762337 0.9402430 0.9587068







################################################################
### Levensohn Self Report Psychopathy Scale  
### Knack et al. (2021) https://osf.io/5ce7m
### https://osf.io/gau8f/?view_only=

mydata <- as.data.frame(Psychopathy_Knack[,3:21])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 5.42; eigenvalue 2 = .92

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 7.06; eigenvalue 2 = 1.1

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.29, RMSEA=.096, RMSR=.08, TLI=.725

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.135, RMSR=.09, TLI=.672

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities with response bias

# Give solution with hypothesized 3 factors
# %variance explained=.39, RMSEA=.072, RMSR=.05, TLI=.846
# With factor correlations of 
#      MR1   MR3   MR2
#MR1  1.00 -0.49  0.45
#MR3 -0.49  1.00 -0.36
#MR2  0.45 -0.36  1.00


fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.11, RMSR=.05, TLI=.782
#      MR1   MR3   MR2
#MR1  1.00 -0.53  0.46
#MR3 -0.53  1.00 -0.37
#MR2  0.46 -0.37  1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~   LSRP_1+LSRP_2+LSRP_3+LSRP_4+LSRP_5+LSRP_7+LSRP_9+LSRP_10+LSRP_11+LSRP_13+
             LSRP_16+LSRP_17+LSRP_18+LSRP_21+LSRP_22+LSRP_23+LSRP_24+LSRP_25+LSRP_26
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.952       0.867
#Tucker-Lewis Index (TLI)                       0.946       0.850
#Robust Comparative Fit Index (CFI)                         0.728
#Robust Tucker-Lewis Index (TLI)                            0.694
#RMSEA                                          0.104       0.114
#Robust RMSEA                                               0.132
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .361

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on theoretical analysis (https://journals.sagepub.com/doi/pdf/10.1177/1073191116637421)
EGAmodel= '
  egocentric =~   LSRP_1+LSRP_3+LSRP_5+LSRP_7+LSRP_9+LSRP_11+LSRP_13+LSRP_17+LSRP_21+LSRP_23
  callous    =~   LSRP_22+LSRP_24+LSRP_25+LSRP_26
  antiscoc   =~   LSRP_2+LSRP_4+LSRP_10+LSRP_16+LSRP_18
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.981       0.938
#Tucker-Lewis Index (TLI)                       0.978       0.928
#Robust Comparative Fit Index (CFI)                         0.823
#Robust Tucker-Lewis Index (TLI)                            0.797
#RMSEA                                          0.067       0.079
#Robust RMSEA                                               0.108
#SRMR                                           0.069       0.069

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#egocentric ~~                                                         
#  callous          -0.723    0.034  -21.451    0.000   -0.723   -0.723
#  antiscoc          0.682    0.036   18.975    0.000    0.682    0.682
#callous ~~                                                            
#  antiscoc         -0.415    0.055   -7.528    0.000   -0.415   -0.415

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .447

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.642       0.537
#Tucker-Lewis Index (TLI)                       0.597       0.479
#Robust Comparative Fit Index (CFI)                         0.721
#Robust Tucker-Lewis Index (TLI)                            0.687
#RMSEA                                          0.286       0.213
#Robust RMSEA                                               0.134
#SRMR                                           0.245       0.245

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .437

# Bifactor model
BIFmodel= '
  egocentric =~   LSRP_1+LSRP_3+LSRP_5+LSRP_7+LSRP_9+LSRP_11+LSRP_13+LSRP_17+LSRP_21+LSRP_23
  callous    =~   LSRP_22+LSRP_24+LSRP_25+LSRP_26
  antiscoc   =~   LSRP_2+LSRP_4+LSRP_10+LSRP_16+LSRP_18
 general=~   LSRP_1+LSRP_2+LSRP_3+LSRP_4+LSRP_5+LSRP_7+LSRP_9+LSRP_10+LSRP_11+LSRP_13+
             LSRP_16+LSRP_17+LSRP_18+LSRP_21+LSRP_22+LSRP_23+LSRP_24+LSRP_25+LSRP_26
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.954
#Tucker-Lewis Index (TLI)                       0.986       0.940
#Robust Comparative Fit Index (CFI)                         0.862
#Robust Tucker-Lewis Index (TLI)                            0.823
#RMSEA                                          0.054       0.072
#Robust RMSEA                                               0.101
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .49

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  egocentric =~                                                         
#    LSRP_1            0.445    0.059    7.481    0.000    0.445    0.445
#    LSRP_3            0.077    0.062    1.246    0.213    0.077    0.077
#    LSRP_5            0.105    0.057    1.843    0.065    0.105    0.105
#    LSRP_7            0.398    0.059    6.768    0.000    0.398    0.398
#    LSRP_9            0.550    0.071    7.744    0.000    0.550    0.550
#    LSRP_11           0.357    0.056    6.423    0.000    0.357    0.357
#    LSRP_13           0.262    0.060    4.379    0.000    0.262    0.262
#    LSRP_17           0.050    0.059    0.838    0.402    0.050    0.050
#    LSRP_21          -0.249    0.075   -3.308    0.001   -0.249   -0.249
#    LSRP_23          -0.072    0.064   -1.128    0.259   -0.072   -0.072
#  callous =~                                                            
#    LSRP_22           0.458    0.056    8.136    0.000    0.458    0.458
#    LSRP_24           0.402    0.060    6.674    0.000    0.402    0.402
#    LSRP_25           0.505    0.050   10.025    0.000    0.505    0.505
#    LSRP_26           0.573    0.052   10.997    0.000    0.573    0.573
#  antiscoc =~                                                           
#    LSRP_2            0.417    0.054    7.700    0.000    0.417    0.417
#    LSRP_4            0.488    0.059    8.332    0.000    0.488    0.488
#    LSRP_10           0.503    0.048   10.379    0.000    0.503    0.503
#    LSRP_16           0.572    0.056   10.159    0.000    0.572    0.572
#    LSRP_18           0.461    0.047    9.723    0.000    0.461    0.461
#  general =~                                                            
#    LSRP_1            0.512    0.048   10.706    0.000    0.512    0.512
#    LSRP_2            0.563    0.043   13.149    0.000    0.563    0.563
#    LSRP_3            0.840    0.025   33.060    0.000    0.840    0.840
#    LSRP_4            0.326    0.051    6.449    0.000    0.326    0.326
#    LSRP_5            0.873    0.021   42.083    0.000    0.873    0.873
#    LSRP_7            0.645    0.041   15.681    0.000    0.645    0.645
#    LSRP_9            0.480    0.052    9.241    0.000    0.480    0.480
#    LSRP_10           0.433    0.046    9.477    0.000    0.433    0.433
#    LSRP_11           0.575    0.040   14.261    0.000    0.575    0.575
#    LSRP_13           0.580    0.041   14.183    0.000    0.580    0.580
#    LSRP_16           0.440    0.049    8.949    0.000    0.440    0.440
#    LSRP_17           0.695    0.030   22.932    0.000    0.695    0.695
#    LSRP_18           0.452    0.048    9.370    0.000    0.452    0.452
#    LSRP_21           0.726    0.034   21.339    0.000    0.726    0.726
#    LSRP_22          -0.474    0.040  -11.710    0.000   -0.474   -0.474
#    LSRP_23           0.823    0.031   26.277    0.000    0.823    0.823
#    LSRP_24          -0.406    0.046   -8.848    0.000   -0.406   -0.406
#    LSRP_25          -0.602    0.038  -16.045    0.000   -0.602   -0.602
#    LSRP_26          -0.524    0.046  -11.334    0.000   -0.524   -0.524

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.6858524      0.6432749      0.8282266      0.6361692 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#egocentric 0.1659162 0.09539289 0.8340838 0.8823555 0.08050587 0.5410290 0.7884047
#callous    0.4818098 0.09676389 0.5181902 0.7384009 0.34724008 0.5632865 0.8017444
#antiscoc   0.5440683 0.12199077 0.4559317 0.7435839 0.39608329 0.6177665 0.8126648
#general    0.6858524 0.68585244 0.6858524 0.8282266 0.63616922 0.9370237 0.9635025







################################################################
###  Five Factor Machiavellianism Inventory (Serbian) https://osf.io/teg8j/download
### Grabovac et al. (2022) https://osf.io/quspg

Machiavellism_Grabovac <- read.csv("Machiavellism_Grabovac.csv")
colnames(Machiavellism_Grabovac)
mydata <- as.data.frame(Machiavellism_Grabovac[,1:52])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 6 components
# Eigenvalue 1 = 6.11; eigenvalue 2 = 4.5

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# does not calculate because min and max are not the same in all columns
# data Pearson used

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.12, RMSEA=.086, RMSR=.13, TLI=.278

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# does not calculate because min and max are not the same in all columns

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities with response bias

# Give solution with hypothesized 5 factors
fit3 <- fa(mydata,5)
fit3
diagram(fit3)
# %variance explained=.34, RMSEA=.054, RMSR=.06, TLI=.712
# With factor correlations of 
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00  0.13  0.19 -0.13  0.14
#MR2  0.13  1.00 -0.13  0.01 -0.13
#MR4  0.19 -0.13  1.00 -0.10  0.05
#MR3 -0.13  0.01 -0.10  1.00  0.10
#MR5  0.14 -0.13  0.05  0.10  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# cannot be calculated

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  ffmi1+ffmi2+ffmi3+ffmi4+ffmi5+ffmi6+ffmi7+ffmi8+ffmi9+ffmi10+ffmi11+ffmi12+
            ffmi13+ffmi14+ffmi15+ffmi16+ffmi17+ffmi18+ffmi19+ffmi20+ffmi21+ffmi22+ffmi23+
            ffmi24+ffmi25+ffmi26+ffmi27+ffmi28+ffmi29+ffmi30+ffmi31+ffmi32+ffmi33+ffmi34+
            ffmi35+ffmi36+ffmi37+ffmi38+ffmi39+ffmi40+ffmi41+ffmi42+ffmi43+ffmi44+ffmi45+
            ffmi46+ffmi47+ffmi48+ffmi49+ffmi50+ffmi51+ffmi52
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
#does not converge
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.294       0.280
#Tucker-Lewis Index (TLI)                       0.265       0.250
#Robust Comparative Fit Index (CFI)                         0.292
#Robust Tucker-Lewis Index (TLI)                            0.263
#RMSEA                                          0.094       0.092
#Robust RMSEA                                               0.094
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .600

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on initial analysis (https://osf.io/teg8j/download)
EGAmodel= '
  agency     =~   ffmi1+ffmi14+ffmi27+ffmi40+ffmi2+ffmi15+ffmi28+ffmi41+ffmi4+ffmi17+ffmi30+ffmi43+
                  ffmi5+ffmi18+ffmi31+ffmi44+ffmi10+ffmi23+ffmi36+ffmi49+ffmi7+ffmi20+ffmi33+ffmi46
  antagonism =~   ffmi3+ffmi16+ffmi29+ffmi42+ffmi8+ffmi21+ffmi34+ffmi47+ffmi11+ffmi24+ffmi37+ffmi50+
                  ffmi12+ffmi25+ffmi38+ffmi51+ffmi13+ffmi26+ffmi39+ffmi52
  planfulness=~   ffmi6+ffmi19+ffmi32+ffmi45+ffmi9+ffmi22+ffmi35+ffmi48
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.458       0.455
#Tucker-Lewis Index (TLI)                       0.435       0.431
#Robust Comparative Fit Index (CFI)                         0.462
#Robust Tucker-Lewis Index (TLI)                            0.439
#RMSEA                                          0.083       0.080
#Robust RMSEA                                               0.082
#SRMR                                           0.115       0.115

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#agency ~~                                                             
#  antagonism        0.086    0.131    0.654    0.513    0.086    0.086
#  planfulness       0.345    0.107    3.225    0.001    0.345    0.345
#antagonism ~~                                                         
#  planfulness      -0.279    0.105   -2.645    0.008   -0.279   -0.279

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .171

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.450       0.448
#Tucker-Lewis Index (TLI)                       0.427       0.425
#Robust Comparative Fit Index (CFI)                         0.454
#Robust Tucker-Lewis Index (TLI)                            0.432
#RMSEA                                          0.083       0.081
#Robust RMSEA                                               0.082
#SRMR                                           0.121       0.121

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .166

# Bifactor model
BIFmodel= '
  agency     =~   ffmi1+ffmi14+ffmi27+ffmi40+ffmi2+ffmi15+ffmi28+ffmi41+ffmi4+ffmi17+ffmi30+ffmi43+
                  ffmi5+ffmi18+ffmi31+ffmi44+ffmi10+ffmi23+ffmi36+ffmi49+ffmi7+ffmi20+ffmi33+ffmi46
  antagonism =~   ffmi3+ffmi16+ffmi29+ffmi42+ffmi8+ffmi21+ffmi34+ffmi47+ffmi11+ffmi24+ffmi37+ffmi50+
                  ffmi12+ffmi25+ffmi38+ffmi51+ffmi13+ffmi26+ffmi39+ffmi52
  planfulness=~   ffmi6+ffmi19+ffmi32+ffmi45+ffmi9+ffmi22+ffmi35+ffmi48
 general=~  ffmi1+ffmi2+ffmi3+ffmi4+ffmi5+ffmi6+ffmi7+ffmi8+ffmi9+ffmi10+ffmi11+ffmi12+
            ffmi13+ffmi14+ffmi15+ffmi16+ffmi17+ffmi18+ffmi19+ffmi20+ffmi21+ffmi22+ffmi23+
            ffmi24+ffmi25+ffmi26+ffmi27+ffmi28+ffmi29+ffmi30+ffmi31+ffmi32+ffmi33+ffmi34+
            ffmi35+ffmi36+ffmi37+ffmi38+ffmi39+ffmi40+ffmi41+ffmi42+ffmi43+ffmi44+ffmi45+
            ffmi46+ffmi47+ffmi48+ffmi49+ffmi50+ffmi51+ffmi52
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.614       0.607
#Tucker-Lewis Index (TLI)                       0.581       0.573
#Robust Comparative Fit Index (CFI)                         0.616
#Robust Tucker-Lewis Index (TLI)                            0.584
#RMSEA                                          0.071       0.069
#Robust RMSEA                                               0.070
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .282

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  agency =~                                                             
#    ffmi1             0.353    0.079    4.485    0.000    0.353    0.328
#    ffmi14           -0.059    0.086   -0.685    0.493   -0.059   -0.050
#    ffmi27           -0.113    0.106   -1.065    0.287   -0.113   -0.090
#    ffmi40            0.030    0.101    0.298    0.766    0.030    0.026
#    ffmi2             0.356    0.095    3.742    0.000    0.356    0.330
#    ffmi15            0.289    0.087    3.333    0.001    0.289    0.325
#    ffmi28            0.569    0.084    6.794    0.000    0.569    0.505
#    ffmi41            0.548    0.090    6.090    0.000    0.548    0.494
#    ffmi4             0.499    0.083    6.005    0.000    0.499    0.464
#    ffmi17            0.247    0.069    3.599    0.000    0.247    0.248
#    ffmi30            0.518    0.082    6.318    0.000    0.518    0.498
#    ffmi43            0.925    0.079   11.782    0.000    0.925    0.696
#    ffmi5             0.702    0.075    9.303    0.000    0.702    0.650
#    ffmi18            0.247    0.091    2.716    0.007    0.247    0.227
#    ffmi31            0.407    0.067    6.043    0.000    0.407    0.502
#    ffmi44            0.935    0.069   13.641    0.000    0.935    0.725
#    ffmi10            0.585    0.076    7.717    0.000    0.585    0.576
#    ffmi23            0.611    0.070    8.667    0.000    0.611    0.634
#    ffmi36            0.742    0.064   11.582    0.000    0.742    0.670
#    ffmi49            0.255    0.106    2.407    0.016    0.255    0.221
#    ffmi7             0.404    0.105    3.839    0.000    0.404    0.382
#    ffmi20            0.488    0.102    4.781    0.000    0.488    0.391
#    ffmi33            0.446    0.114    3.907    0.000    0.446    0.344
#    ffmi46            0.223    0.142    1.571    0.116    0.223    0.164
#  antagonism =~                                                         
#    ffmi3             0.296    0.080    3.709    0.000    0.296    0.409
#    ffmi16            0.431    0.086    4.992    0.000    0.431    0.519
#    ffmi29            0.442    0.076    5.845    0.000    0.442    0.529
#    ffmi42            0.401    0.087    4.620    0.000    0.401    0.466
#    ffmi8             0.227    0.104    2.183    0.029    0.227    0.209
#    ffmi21            0.073    0.117    0.622    0.534    0.073    0.072
#    ffmi34            0.181    0.113    1.599    0.110    0.181    0.164
#    ffmi47            0.163    0.124    1.312    0.190    0.163    0.149
#    ffmi11            0.332    0.129    2.584    0.010    0.332    0.283
#    ffmi24            0.372    0.108    3.443    0.001    0.372    0.322
#    ffmi37            0.296    0.103    2.881    0.004    0.296    0.267
#    ffmi50            0.457    0.116    3.940    0.000    0.457    0.364
#    ffmi12            0.292    0.100    2.916    0.004    0.292    0.241
#    ffmi25            0.605    0.124    4.889    0.000    0.605    0.586
#    ffmi38            0.263    0.091    2.885    0.004    0.263    0.241
#    ffmi51            0.680    0.105    6.472    0.000    0.680    0.660
#    ffmi13           -0.014    0.083   -0.168    0.867   -0.014   -0.016
#    ffmi26            0.285    0.112    2.558    0.011    0.285    0.281
#    ffmi39            0.594    0.124    4.802    0.000    0.594    0.509
#    ffmi52            0.538    0.119    4.505    0.000    0.538    0.499
#  planfulness =~                                                        
#    ffmi6             0.551    0.097    5.676    0.000    0.551    0.509
#    ffmi19            0.295    0.117    2.516    0.012    0.295    0.261
#    ffmi32            0.102    0.104    0.982    0.326    0.102    0.086
#    ffmi45            0.421    0.107    3.919    0.000    0.421    0.388
#    ffmi9             0.641    0.076    8.435    0.000    0.641    0.698
#    ffmi22            0.748    0.075    9.986    0.000    0.748    0.702
#    ffmi35            0.740    0.080    9.219    0.000    0.740    0.631
#    ffmi48            0.695    0.076    9.105    0.000    0.695    0.581
#  general =~                                                            
#    ffmi1             0.150    0.092    1.631    0.103    0.150    0.139
#    ffmi2            -0.065    0.140   -0.463    0.643   -0.065   -0.060
#    ffmi3             0.112    0.102    1.107    0.268    0.112    0.155
#    ffmi4             0.317    0.096    3.289    0.001    0.317    0.295
#    ffmi5             0.092    0.115    0.803    0.422    0.092    0.086
#    ffmi6            -0.222    0.118   -1.883    0.060   -0.222   -0.205
#    ffmi7            -0.154    0.116   -1.328    0.184   -0.154   -0.146
#    ffmi8             0.294    0.109    2.706    0.007    0.294    0.271
#    ffmi9            -0.000    0.113   -0.004    0.997   -0.000   -0.000
#    ffmi10            0.406    0.079    5.164    0.000    0.406    0.400
#    ffmi11            0.494    0.116    4.278    0.000    0.494    0.420
#    ffmi12            0.305    0.109    2.789    0.005    0.305    0.252
#    ffmi13            0.210    0.071    2.955    0.003    0.210    0.239
#    ffmi14            0.706    0.097    7.278    0.000    0.706    0.595
#    ffmi15            0.156    0.101    1.551    0.121    0.156    0.175
#    ffmi16           -0.311    0.098   -3.167    0.002   -0.311   -0.374
#    ffmi17            0.427    0.100    4.256    0.000    0.427    0.428
#    ffmi18            0.307    0.126    2.444    0.015    0.307    0.282
#    ffmi19           -0.243    0.120   -2.031    0.042   -0.243   -0.216
#    ffmi20           -0.197    0.125   -1.581    0.114   -0.197   -0.158
#    ffmi21            0.560    0.090    6.229    0.000    0.560    0.553
#    ffmi22            0.153    0.141    1.081    0.280    0.153    0.144
#    ffmi23           -0.086    0.101   -0.852    0.394   -0.086   -0.089
#    ffmi24            0.245    0.140    1.742    0.082    0.245    0.211
#    ffmi25            0.014    0.150    0.096    0.924    0.014    0.014
#    ffmi26           -0.124    0.102   -1.216    0.224   -0.124   -0.122
#    ffmi27            0.598    0.126    4.730    0.000    0.598    0.474
#    ffmi28            0.207    0.105    1.978    0.048    0.207    0.184
#    ffmi29           -0.184    0.099   -1.862    0.063   -0.184   -0.220
#    ffmi30            0.318    0.102    3.120    0.002    0.318    0.306
#    ffmi31            0.155    0.102    1.526    0.127    0.155    0.191
#    ffmi32           -0.408    0.110   -3.717    0.000   -0.408   -0.346
#    ffmi33           -0.334    0.153   -2.178    0.029   -0.334   -0.257
#    ffmi34            0.534    0.102    5.239    0.000    0.534    0.483
#    ffmi35           -0.096    0.156   -0.613    0.540   -0.096   -0.082
#    ffmi36            0.371    0.106    3.480    0.001    0.371    0.335
#    ffmi37            0.503    0.119    4.217    0.000    0.503    0.455
#    ffmi38            0.313    0.098    3.183    0.001    0.313    0.286
#    ffmi39           -0.264    0.135   -1.951    0.051   -0.264   -0.226
#    ffmi40            0.595    0.108    5.511    0.000    0.595    0.525
#    ffmi41            0.230    0.120    1.913    0.056    0.230    0.207
#    ffmi42            0.399    0.103    3.875    0.000    0.399    0.463
#    ffmi43            0.252    0.135    1.865    0.062    0.252    0.190
#    ffmi44           -0.056    0.130   -0.428    0.668   -0.056   -0.043
#    ffmi45           -0.260    0.102   -2.554    0.011   -0.260   -0.240
#    ffmi46           -0.219    0.163   -1.344    0.179   -0.219   -0.161
#    ffmi47            0.589    0.093    6.329    0.000    0.589    0.537
#    ffmi48           -0.188    0.114   -1.650    0.099   -0.188   -0.158
#    ffmi49            0.290    0.088    3.300    0.001    0.290    0.251
#    ffmi50            0.429    0.114    3.767    0.000    0.429    0.341
#    ffmi51            0.134    0.154    0.866    0.386    0.134    0.130
#    ffmi52           -0.220    0.120   -1.827    0.068   -0.220   -0.204

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.3136153      0.6274510      0.8356320      0.1968341 

#$FactorLevelIndices
#             ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#agency      0.7053929 0.3316437 0.2946071 0.8565810 0.7133793 0.8813421 0.9394664
#antagonism  0.5685940 0.2016625 0.4314060 0.7984257 0.6166564 0.7980241 0.8998202
#planfulness 0.8738521 0.1530785 0.1261479 0.7459134 0.6895771 0.7872419 0.8894029
#general     0.3136153 0.3136153 0.3136153 0.8356320 0.1968341 0.8478001 0.9244244







################################################################
###  Machiavellian Personality Scale (Serbian) https://osf.io/teg8j/download
### Grabovac et al. (2022) https://osf.io/quspg

### not included because too many errors and already dataset from this source

mydata <- as.data.frame(Machiavellism_Grabovac[,53:68])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 3.29; eigenvalue 2 = .89

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 3.91; eigenvalue 2 = 1.07

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.21, RMSEA=.094, RMSR=.1, TLI=.637

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.24, RMSEA=.127, RMSR=.11, TLI=.56

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities with response bias

# Give solution with hypothesized 4 factors
fit3 <- fa(mydata,4)
fit3
diagram(fit3)
# %variance explained=.39, RMSEA=.033, RMSR=.04, TLI=.954
# With factor correlations of 
#     MR1  MR3  MR2  MR4
#MR1 1.00 0.44 0.12 0.19
#MR3 0.44 1.00 0.21 0.12
#MR2 0.12 0.21 1.00 0.15
#MR4 0.19 0.12 0.15 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.071, RMSR=.04, TLI=.861
# With factor correlations of 
#     MR1  MR4  MR2  MR3
#MR1 1.00 0.44 0.10 0.16
#MR4 0.44 1.00 0.23 0.15
#MR2 0.10 0.23 1.00 0.17
#MR3 0.16 0.15 0.17 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  mps1+mps2+mps3+mps4+mps5+mps6+mps7+mps8+mps9+mps10+mps11+mps12+mps13+mps14+mps15+mps16
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.889       0.783
#Tucker-Lewis Index (TLI)                       0.871       0.750
#Robust Comparative Fit Index (CFI)                         0.639
#Robust Tucker-Lewis Index (TLI)                            0.583
#RMSEA                                          0.110       0.113
#Robust RMSEA                                               0.127
#SRMR                                           0.106       0.106

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .297

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 4 hypothezied factors (https://osf.io/quspg)
EGAmodel= '
  amorality  =~   mps4+mps7+mps10+mps12+mps15
  control    =~   mps2+mps5+mps13
  status     =~   mps1+mps3+mps8
  distrust   =~   mps6+mps9+mps11+mps14+mps16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
# CFA_model2 <- cfa(EGAmodel,mydata,estimator='MLR',std.lv=TRUE)
# bad estimate covariance status - distrust
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.915       0.825
#Tucker-Lewis Index (TLI)                       0.896       0.786
#Robust Comparative Fit Index (CFI)                         0.723
#Robust Tucker-Lewis Index (TLI)                            0.661
#RMSEA                                          0.099       0.104
#Robust RMSEA                                               0.115
#SRMR                                           0.098       0.098

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  amorality ~~                                                          
#    control           0.722    0.069   10.406    0.000    0.722    0.722
#    status            0.553    0.088    6.298    0.000    0.553    0.553
#    distrust          0.631    0.078    8.119    0.000    0.631    0.631
#  control ~~                                                            
#    status            0.740    0.067   10.984    0.000    0.740    0.740
#    distrust          0.797    0.066   12.136    0.000    0.797    0.797
#  status ~~                                                             
#    distrust          1.135    0.064   17.726    0.000    1.135    1.135

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .299

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.250       0.187
#Tucker-Lewis Index (TLI)                       0.135       0.062
#Robust Comparative Fit Index (CFI)                         0.426
#Robust Tucker-Lewis Index (TLI)                            0.337
#RMSEA                                          0.286       0.218
#Robust RMSEA                                               0.161
#SRMR                                           0.216       0.216

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .240

# Bifactor model
BIFmodel= '
  amorality  =~   mps4+mps7+mps10+mps12+mps15
  control    =~   mps2+mps5+mps13
  status     =~   mps1+mps3+mps8
  distrust   =~   mps6+mps9+mps11+mps14+mps16
  general=~  mps1+mps2+mps3+mps4+mps5+mps6+mps7+mps8+mps9+mps10+mps11+mps12+mps13+mps14+mps15+mps16
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
# warning
#CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# does not converge either misestimate of item mps8, which will be left out
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .293

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  amorality =~                                                          
#    mps4              0.129    0.095    1.366    0.172    0.129    0.136
#    mps7              0.200    0.122    1.635    0.102    0.200    0.199
#    mps10             0.362    0.109    3.312    0.001    0.362    0.325
#    mps12             0.651    0.129    5.037    0.000    0.651    0.559
#    mps15             0.804    0.184    4.374    0.000    0.804    0.672
#  control =~                                                            
#    mps2              0.339    0.235    1.442    0.149    0.339    0.288
#    mps5              0.921    0.544    1.693    0.090    0.921    0.973
#    mps13             0.152    0.136    1.118    0.264    0.152    0.160
#  status =~                                                             
#    mps1             -0.003    0.009   -0.364    0.716   -0.003   -0.003
#    mps3             -0.011    0.009   -1.304    0.192   -0.011   -0.010
#    mps8              9.494    0.056  170.544    0.000    9.494    9.162
#  distrust =~                                                           
#    mps6              0.148    0.157    0.942    0.346    0.148    0.148
#    mps9              0.271    0.163    1.659    0.097    0.271    0.212
#    mps11            -0.138    0.168   -0.822    0.411   -0.138   -0.151
#    mps14             0.592    0.202    2.930    0.003    0.592    0.470
#    mps16             0.357    0.199    1.791    0.073    0.357    0.429
#  general =~                                                            
#    mps1              0.734    0.089    8.270    0.000    0.734    0.589
#    mps2              0.529    0.097    5.467    0.000    0.529    0.450
#    mps3              0.842    0.086    9.809    0.000    0.842    0.751
#    mps4              0.412    0.085    4.871    0.000    0.412    0.432
#    mps5              0.486    0.089    5.468    0.000    0.486    0.514
#    mps6              0.147    0.086    1.711    0.087    0.147    0.147
#    mps7              0.238    0.091    2.615    0.009    0.238    0.237
#    mps8              0.277    0.111    2.495    0.013    0.277    0.267
#    mps9              0.957    0.080   11.915    0.000    0.957    0.749
#    mps10             0.284    0.093    3.047    0.002    0.284    0.254
#    mps11             0.484    0.083    5.857    0.000    0.484    0.527
#    mps12             0.088    0.104    0.849    0.396    0.088    0.075
#    mps13             0.436    0.089    4.913    0.000    0.436    0.461
#    mps14             0.577    0.097    5.922    0.000    0.577    0.458
#    mps15             0.284    0.106    2.685    0.007    0.284    0.238
#    mps16             0.024    0.085    0.279    0.780    0.024    0.029

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.4150760      0.7714286      0.8465742      0.7148364 

#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#amorality 0.6179480 0.10616703 0.3820520 0.6190849 0.4102827 0.6034586 0.7841833
#control   0.7562073 0.33839580 0.2437927 0.7181382 0.2935341 3.2450298 2.4565897
#status    0.3973672 0.06845743 0.6026328 0.8277382 0.2001994 0.6337789 0.9613315
#distrust  0.3449848 0.07190377 0.6550152 0.6213446 0.1881097 0.4951397 0.7281701
#general   0.4150760 0.41507597 0.4150760 0.8465742 0.7148364 0.8746579 0.9299154







################################################################
### Short dark triad (SD3) scale
### Nielsen et al. (2023) https://journals.sagepub.com/doi/pdf/10.1177/25152459231177713
### Data: https://osf.io/r46uw/?view_only=

Dark_triad_Nielsen <- read.csv("Dark_triad_Nielsen.csv")
colnames(Dark_triad_Nielsen)
mydata <- as.data.frame(Dark_triad_Nielsen[,3:29])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors 3 components
# Eigenvalue 1 = 7.00; eigenvalue 2 = 1.59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 8.21; eigenvalue 2 = 1.85

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.26, RMSEA=.085, RMSR=.08, TLI=.707

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.103, RMSR=.09, TLI=.673

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities and response bias

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities 

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.053, RMSR=.03, TLI=.885
#     MR1  MR2  MR3
#MR1 1.00 0.43 0.58
#MR2 0.43 1.00 0.28
#MR3 0.58 0.28 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.069, RMSR=.04, TLI=.853
#     MR1  MR3  MR2
#MR1 1.00 0.56 0.46
#MR3 0.56 1.00 0.27
#MR2 0.46 0.27 1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~   narc1+narc2+narc3+narc4+narc5+narc6+narc7+narc8+narc9+
             mach1+mach2+mach3+mach4+mach5+mach6+mach7+mach8+mach9+
             psyc1+psyc2+psyc3+psyc4+psyc5+psyc6+psyc7+psyc8+psyc9
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.935       0.801
#Tucker-Lewis Index (TLI)                       0.929       0.785
#Robust Comparative Fit Index (CFI)                         0.700
#Robust Tucker-Lewis Index (TLI)                            0.675
#RMSEA                                          0.102       0.109
#Robust RMSEA                                               0.103
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .286

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.731       0.732
#Tucker-Lewis Index (TLI)                       0.709       0.709
#Robust Comparative Fit Index (CFI)                         0.734
#Robust Tucker-Lewis Index (TLI)                            0.711
#RMSEA                                          0.084       0.077
#Robust RMSEA                                               0.084
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .215

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on theoretical analysis
EGAmodel= '
 narcisism =~   narc1+narc2+narc3+narc4+narc5+narc6+narc7+narc8+narc9
 machiavel =~   mach1+mach2+mach3+mach4+mach5+mach6+mach7+mach8+mach9
 psychotic =~   psyc1+psyc2+psyc3+psyc4+psyc5+psyc6+psyc7+psyc8+psyc9
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.966       0.892
#Tucker-Lewis Index (TLI)                       0.963       0.881
#Robust Comparative Fit Index (CFI)                         0.815
#Robust Tucker-Lewis Index (TLI)                            0.798
#RMSEA                                          0.074       0.081
#Robust RMSEA                                               0.081
#SRMR                                           0.066       0.066

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  narcisism ~~                                                          
#    machiavel         0.575    0.017   34.192    0.000    0.575    0.575
#    psychotic         0.613    0.017   35.947    0.000    0.613    0.613
#  machiavel ~~                                                          
#    psychotic         0.853    0.009   99.686    0.000    0.853    0.853

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .450

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.846       0.846
#Tucker-Lewis Index (TLI)                       0.831       0.832
#Robust Comparative Fit Index (CFI)                         0.848
#Robust Tucker-Lewis Index (TLI)                            0.834
#RMSEA                                          0.064       0.059
#Robust RMSEA                                               0.064
#SRMR                                           0.060       0.060

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  narcisism ~~                                                          
#    machiavel         0.551    0.023   23.718    0.000    0.551    0.551
#    psychotic         0.575    0.024   24.366    0.000    0.575    0.575
#  machiavel ~~                                                          
#    psychotic         0.857    0.012   71.885    0.000    0.857    0.857

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .358

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.726       0.724
#Tucker-Lewis Index (TLI)                       0.703       0.701
#Robust Comparative Fit Index (CFI)                         0.728
#Robust Tucker-Lewis Index (TLI)                            0.705
#RMSEA                                          0.085       0.078
#Robust RMSEA                                               0.085
#SRMR                                           0.189       0.189

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .359

# Bifactor model
BIFmodel= '
 narcisism =~   narc1+narc2+narc3+narc4+narc5+narc6+narc7+narc8+narc9
 machiavel =~   mach1+mach2+mach3+mach4+mach5+mach6+mach7+mach8+mach9
 psychotic =~   psyc1+psyc2+psyc3+psyc4+psyc5+psyc6+psyc7+psyc8+psyc9
 general=~   narc1+narc2+narc3+narc4+narc5+narc6+narc7+narc8+narc9+
             mach1+mach2+mach3+mach4+mach5+mach6+mach7+mach8+mach9+
             psyc1+psyc2+psyc3+psyc4+psyc5+psyc6+psyc7+psyc8+psyc9
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.883       0.881
#Tucker-Lewis Index (TLI)                       0.862       0.860
#Robust Comparative Fit Index (CFI)                         0.885
#Robust Tucker-Lewis Index (TLI)                            0.865
#RMSEA                                          0.058       0.054
#Robust RMSEA                                               0.058
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .353

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6616252      0.6923077      0.9089237      0.7809254 
#
#$FactorLevelIndices
#             ECV_SS     ECV_SG    ECV_GS     Omega       OmegaH         H        FD
#narcisism 0.6334939 0.17846794 0.3665061 0.7820559 0.4993664581 0.6939517 0.8370312
#machiavel 0.3132291 0.12410019 0.6867709 0.8585733 0.2727731267 0.5906426 0.7639020
#psychotic 0.1111721 0.03580671 0.8888279 0.7903563 0.0002473987 0.2796511 0.6508828
#general   0.6616252 0.66162516 0.6616252 0.9089237 0.7809253685 0.9123774 0.9443611











