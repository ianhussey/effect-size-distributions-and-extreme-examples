################################################################
### UPPS-P Impulsive Behavior Scale 
### 
### German version by Wuellhorst et al. (2023) https://doi.org/10.17605/OSF.IO/39DGR
### data validation sample used available at https://osf.io/39dgr/
###
### 1. 

library(readxl)
validation_sample_Wuellhorst <- read_excel("validation_sample_Wuellhorst.xlsx")
colnames(validation_sample_Wuellhorst)
mydata <- as.data.frame(validation_sample_Wuellhorst[,2:60])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
  "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
  "upps16","upps17","upps18","upps19","upps20","upps21","upps22","upps23","upps24","upps25",
  "upps26","upps27","upps28","upps29","upps30","upps31","upps32","upps33","upps34","upps35",
  "upps36","upps37","upps38","upps39","upps40","upps41","upps42","upps43","upps44","upps45",
  "upps46","upps47","upps48","upps49","upps50","upps51","upps52","upps53","upps54","upps55",
  "upps56","upps57","upps58","upps59")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .93, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 12.63; eigenvalue 2 = 5.89

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 5 components
# Eigenvalue 1 = 15.24; eigenvalue 2 = 6.91

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.21, RMSEA=.093, RMSR=.14, TLI=.428

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.177, RMSR=.16, TLI=.206

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.043, RMSR=.04, TLI=.876
#      MR1   MR2   MR3   MR5   MR4
#MR1  1.00  0.21 -0.25  0.48 -0.20
#MR2  0.21  1.00  0.07 -0.09 -0.29
#MR3 -0.25  0.07  1.00 -0.32  0.24
#MR5  0.48 -0.09 -0.32  1.00 -0.08
#MR4 -0.20 -0.29  0.24 -0.08  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.157, RMSR=.04, TLI=.375
#      MR1   MR2   MR3   MR5   MR4
#MR1  1.00  0.21 -0.27  0.49 -0.23
#MR2  0.21  1.00  0.07 -0.11 -0.29
#MR3 -0.27  0.07  1.00 -0.31  0.26
#MR5  0.49 -0.11 -0.31  1.00 -0.07
#MR4 -0.23 -0.29  0.26 -0.07  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20+upps21+upps22+upps23+upps24+
  upps25+upps26+upps27+upps28+upps29+upps30+upps31+upps32+upps33+upps34+upps35+upps36+
  upps37+upps38+upps39+upps40+upps41+upps42+upps43+upps44+upps45+upps46+upps47+upps48+
  upps49+upps50+upps51+upps52+upps53+upps54+upps55+upps56+upps57+upps58+upps59
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.794       0.620
#Tucker-Lewis Index (TLI)                       0.786       0.606
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.175       0.109
#Robust RMSEA                                                  NA
#SRMR                                           0.165       0.165

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .267

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.445       0.444
#Tucker-Lewis Index (TLI)                       0.425       0.424
#Robust Comparative Fit Index (CFI)                         0.449
#Robust Tucker-Lewis Index (TLI)                            0.430
#RMSEA                                          0.098       0.094
#Robust RMSEA                                               0.097
#SRMR                                           0.135       0.135

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .141

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (EGA based) 
EGAmodel= '
 factor1 =~ upps59+upps54+upps49+upps52+upps25+upps45+upps20+upps57+upps58+
            upps12+upps35+upps40+upps30+upps5+upps15+upps10
 factor2 =~ upps1+upps38+upps23+upps41+upps3+upps8+upps31+upps13+upps26+
            upps56+upps18+upps46+upps51+upps36
 factor3 =~ upps50+upps34+upps44+upps29+upps53+upps39+upps2+upps17+upps22+
            upps7+upps24
 factor4 =~ upps14+upps19+upps4+upps27+upps42+upps9+upps37+upps32+upps47
 factor5 =~ upps11+upps16+upps55+upps33+upps48+upps43+upps21+upps6+upps28
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.901
#Tucker-Lewis Index (TLI)                       0.954       0.896
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.081       0.056
#Robust RMSEA                                                  NA
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .490

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.309    0.052   -5.960    0.000   -0.309   -0.309
#    factor3           0.673    0.037   18.147    0.000    0.673    0.673
#    factor4          -0.324    0.052   -6.273    0.000   -0.324   -0.324
#    factor5          -0.415    0.050   -8.285    0.000   -0.415   -0.415
#  factor2 ~~                                                            
#    factor3           0.051    0.058    0.890    0.373    0.051    0.051
#    factor4          -0.032    0.056   -0.574    0.566   -0.032   -0.032
#    factor5           0.403    0.046    8.807    0.000    0.403    0.403
#  factor3 ~~                                                            
#    factor4          -0.538    0.044  -12.297    0.000   -0.538   -0.538
#    factor5          -0.403    0.049   -8.245    0.000   -0.403   -0.403
#  factor4 ~~                                                            
#    factor5           0.468    0.047    9.930    0.000    0.468    0.468

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.801       0.805
#Tucker-Lewis Index (TLI)                       0.792       0.797
#Robust Comparative Fit Index (CFI)                         0.808
#Robust Tucker-Lewis Index (TLI)                            0.800
#RMSEA                                          0.059       0.056
#Robust RMSEA                                               0.057
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .414

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.286    0.062   -4.639    0.000   -0.286   -0.286
#    factor3           0.642    0.050   12.763    0.000    0.642    0.642
#    factor4          -0.285    0.062   -4.608    0.000   -0.285   -0.285
#    factor5          -0.376    0.063   -5.970    0.000   -0.376   -0.376
#  factor2 ~~                                                            
#    factor3           0.045    0.072    0.624    0.533    0.045    0.045
#    factor4          -0.017    0.067   -0.248    0.804   -0.017   -0.017
#    factor5           0.407    0.064    6.396    0.000    0.407    0.407
#  factor3 ~~                                                            
#    factor4          -0.498    0.059   -8.431    0.000   -0.498   -0.498
#    factor5          -0.332    0.073   -4.538    0.000   -0.332   -0.332
#  factor4 ~~                                                            
#    factor5           0.408    0.068    6.031    0.000    0.408    0.408

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.763       0.767
#Tucker-Lewis Index (TLI)                       0.755       0.759
#Robust Comparative Fit Index (CFI)                         0.771
#Robust Tucker-Lewis Index (TLI)                            0.763
#RMSEA                                          0.064       0.061
#Robust RMSEA                                               0.063
#SRMR                                           0.165       0.165

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .415


# Bifactor model with 5 factors
BIFmodel= '
 factor1 =~ upps59+upps54+upps49+upps52+upps25+upps45+upps20+upps57+upps58+
            upps12+upps35+upps40+upps30+upps5+upps15+upps10
 factor2 =~ upps1+upps38+upps23+upps41+upps3+upps8+upps31+upps13+upps26+
            upps56+upps18+upps46+upps51+upps36
 factor3 =~ upps50+upps34+upps44+upps29+upps53+upps39+upps2+upps17+upps22+
            upps7+upps24
 factor4 =~ upps14+upps19+upps4+upps27+upps42+upps9+upps37+upps32+upps47
 factor5 =~ upps11+upps16+upps55+upps33+upps48+upps43+upps21+upps6+upps28
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20+upps21+upps22+upps23+upps24+
  upps25+upps26+upps27+upps28+upps29+upps30+upps31+upps32+upps33+upps34+upps35+upps36+
  upps37+upps38+upps39+upps40+upps41+upps42+upps43+upps44+upps45+upps46+upps47+upps48+
  upps49+upps50+upps51+upps52+upps53+upps54+upps55+upps56+upps57+upps58+upps59
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.844       0.848
#Tucker-Lewis Index (TLI)                       0.832       0.836
#Robust Comparative Fit Index (CFI)                         0.851
#Robust Tucker-Lewis Index (TLI)                            0.840
#RMSEA                                          0.053       0.050
#Robust RMSEA                                               0.051
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3696209      0.8024547      0.8640392      0.2576119 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5865290 0.18767342 0.41347102 0.9435187 0.5437052 0.8861917 0.9311421
#factor2 0.9634211 0.19036303 0.03657893 0.7920607 0.7911824 0.8915856 0.9457640
#factor3 0.2493718 0.04658459 0.75062815 0.7920655 0.1607543 0.6552786 0.8385793
#factor4 0.6214050 0.09630596 0.37859503 0.7274370 0.6135533 0.8179713 0.9180438
#factor5 0.7782010 0.10945211 0.22179899 0.8410107 0.7043922 0.8376299 0.9177806
#general 0.3696209 0.36962089 0.36962089 0.8640392 0.2576119 0.9321762 0.9445607







################################################################
### UPPS-P Impulsive Behavior Scale 
### 
### French version by Billieux et al. (2022)
### data JOPI used available at https://osf.io/2p6sx/
###
### 2. 

library(readxl)
Billieux <- read_excel("Billieux.xlsx")
colnames(Billieux)
mydata <- as.data.frame(Billieux[,7:26])
mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
                      "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
                      "upps16","upps17","upps18","upps19","upps20")
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .84, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 4.30; eigenvalue 2 = 2.50

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 5.05; eigenvalue 2 = 2.91

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.172, RMSR=.17, TLI=.309

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.206, RMSR=.19, TLI=.279

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.041, RMSR=.01, TLI=.961
#      MR2  MR4   MR3  MR1  MR5
#MR2  1.00 0.46 -0.03 0.11 0.10
#MR4  0.46 1.00  0.11 0.29 0.31
#MR3 -0.03 0.11  1.00 0.12 0.37
#MR1  0.11 0.29  0.12 1.00 0.57
#MR5  0.10 0.31  0.37 0.57 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.055, RMSR=.02, TLI=.949
#      MR2  MR4   MR3  MR1  MR5
#MR2  1.00 0.48 -0.03 0.11 0.10
#MR4  0.48 1.00  0.11 0.30 0.30
#MR3 -0.03 0.11  1.00 0.13 0.38
#MR1  0.11 0.30  0.13 1.00 0.57
#MR5  0.10 0.30  0.38 0.57 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.706       0.545
#Tucker-Lewis Index (TLI)                       0.672       0.492
#Robust Comparative Fit Index (CFI)                         0.301
#Robust Tucker-Lewis Index (TLI)                            0.218
#RMSEA                                          0.269       0.222
#Robust RMSEA                                               0.215
#SRMR                                           0.211       0.211

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .348

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.385       0.340
#Tucker-Lewis Index (TLI)                       0.313       0.263
#Robust Comparative Fit Index (CFI)                         0.385
#Robust Tucker-Lewis Index (TLI)                            0.313
#RMSEA                                          0.172       0.159
#Robust RMSEA                                               0.172
#SRMR                                           0.162       0.162

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .192

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1 =~ upps13+upps19+upps1+upps6
 factor2 =~ upps8+upps11+upps16+upps5
 factor3 =~ upps9+upps18+upps3+upps14
 factor4 =~ upps7+upps17+upps12+upps4
 factor5 =~ upps10+upps15+upps2+upps20
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.967
#Tucker-Lewis Index (TLI)                       0.987       0.960
#Robust Comparative Fit Index (CFI)                         0.941
#Robust Tucker-Lewis Index (TLI)                            0.930
#RMSEA                                          0.054       0.062
#Robust RMSEA                                               0.064
#SRMR                                           0.045       0.045

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .620

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.527    0.006   81.311    0.000    0.527    0.527
#    factor3           0.120    0.008   14.233    0.000    0.120    0.120
#    factor4           0.339    0.008   44.337    0.000    0.339    0.339
#    factor5           0.388    0.008   51.179    0.000    0.388    0.388
#  factor2 ~~                                                            
#    factor3          -0.035    0.008   -4.176    0.000   -0.035   -0.035
#    factor4           0.135    0.008   16.257    0.000    0.135    0.135
#    factor5           0.134    0.008   15.956    0.000    0.134    0.134
#  factor3 ~~                                                            
#    factor4           0.179    0.008   21.984    0.000    0.179    0.179
#    factor5           0.407    0.007   54.765    0.000    0.407    0.407
#  factor4 ~~                                                            
#    factor5           0.715    0.005  147.052    0.000    0.715    0.715

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.952       0.950
#Tucker-Lewis Index (TLI)                       0.943       0.940
#Robust Comparative Fit Index (CFI)                         0.952
#Robust Tucker-Lewis Index (TLI)                            0.943
#RMSEA                                          0.049       0.045
#Robust RMSEA                                               0.049
#SRMR                                           0.040       0.040

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .529

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.495    0.009   58.052    0.000    0.495    0.495
#    factor3           0.120    0.010   12.447    0.000    0.120    0.120
#    factor4           0.335    0.009   36.349    0.000    0.335    0.335
#    factor5           0.398    0.009   43.607    0.000    0.398    0.398
#  factor2 ~~                                                            
#    factor3          -0.031    0.009   -3.304    0.001   -0.031   -0.031
#    factor4           0.134    0.009   14.653    0.000    0.134    0.134
#    factor5           0.145    0.009   15.607    0.000    0.145    0.145
#  factor3 ~~                                                            
#    factor4           0.176    0.010   18.347    0.000    0.176    0.176
#    factor5           0.372    0.009   40.429    0.000    0.372    0.372
#  factor4 ~~                                                            
#    factor5           0.723    0.006  112.413    0.000    0.723    0.723

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.849       0.842
#Tucker-Lewis Index (TLI)                       0.831       0.823
#Robust Comparative Fit Index (CFI)                         0.849
#Robust Tucker-Lewis Index (TLI)                            0.832
#RMSEA                                          0.085       0.078
#Robust RMSEA                                               0.085
#SRMR                                           0.160       0.160

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .519


# Bifactor model with 5 factors
BIFmodel= '
 factor1 =~ upps13+upps19+upps1+upps6
 factor2 =~ upps8+upps11+upps16+upps5
 factor3 =~ upps9+upps18+upps3+upps14
 factor4 =~ upps7+upps17+upps12+upps4
 factor5 =~ upps10+upps15+upps2+upps20
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.945       0.942
#Tucker-Lewis Index (TLI)                       0.931       0.927
#Robust Comparative Fit Index (CFI)                         0.946
#Robust Tucker-Lewis Index (TLI)                            0.931
#RMSEA                                          0.054       0.050
#Robust RMSEA                                               0.054
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .507

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3115941      0.8421053      0.9057948      0.5912899 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.7683946 0.14915581 0.23160538 0.8222214 0.6319916 0.7409427 0.8651247
#factor2 0.9446795 0.21718843 0.05532052 0.8744420 0.8266446 0.8645938 0.9313948
#factor3 0.8953831 0.17571473 0.10461687 0.8255991 0.7411428 0.7925436 0.8919036
#factor4 0.4346485 0.08925874 0.56535145 0.8319424 0.3506964 0.5988218 0.7744678
#factor5 0.3273858 0.05708814 0.67261415 0.7716681 0.2518367 0.4537093 0.6588707
#general 0.3115941 0.31159415 0.31159415 0.9057948 0.5912899 0.8501974 0.8894678






################################################################
### UPPS-P Impulsive Behavior Scale 
### 
### Eben et al 2020 https://royalsocietypublishing.org/doi/pdf/10.1098/rsos.200664
### data available at https://osf.io/xt58b
###
### 3. 

library(readxl)
Eben2020 <- read_excel("Eben2020.xlsx")
colnames(Eben2020)
UPPS2 <- Eben2020[,c(3,11,21)]
colnames(UPPS2)
library(tidyr)
mydata <- spread(UPPS2, trial_t, resp)
mydata <- mydata[,2:21]
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
                      "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
                      "upps16","upps17","upps18","upps19","upps20")
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .76, omega T = .81

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 4 components
# Eigenvalue 1 = 3.64; eigenvalue 2 = 1.61

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 4 components
# Eigenvalue 1 = 4.36; eigenvalue 2 = 1.95

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.112, RMSR=.12, TLI=.471

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.22, RMSEA=.144, RMSR=.14, TLI=.409

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.029, RMSR=.03, TLI=.964
#     MR4  MR1  MR5   MR3   MR2
#MR4 1.00 0.39 0.27  0.11  0.35
#MR1 0.39 1.00 0.55  0.08  0.09
#MR5 0.27 0.55 1.00  0.24  0.05
#MR3 0.11 0.08 0.24  1.00 -0.08
#MR2 0.35 0.09 0.05 -0.08  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.063, RMSR=.04, TLI=.886
#     MR4  MR1  MR5   MR2   MR3
#MR4 1.00 0.39 0.28  0.37  0.11
#MR1 0.39 1.00 0.55  0.09  0.08
#MR5 0.28 0.55 1.00  0.05  0.24
#MR2 0.37 0.09 0.05  1.00 -0.08
#MR3 0.11 0.08 0.24 -0.08  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.769       0.622
#Tucker-Lewis Index (TLI)                       0.742       0.577
#Robust Comparative Fit Index (CFI)                         0.486
#Robust Tucker-Lewis Index (TLI)                            0.425
#RMSEA                                          0.146       0.135
#Robust RMSEA                                               0.144
#SRMR                                           0.134       0.134

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .257

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.527       0.518
#Tucker-Lewis Index (TLI)                       0.472       0.461
#Robust Comparative Fit Index (CFI)                         0.530
#Robust Tucker-Lewis Index (TLI)                            0.475
#RMSEA                                          0.113       0.107
#Robust RMSEA                                               0.112
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .206

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (EGA based) 
EGAmodel= '
 factor1 =~ upps4+upps7+upps11+upps1
 factor2 =~ upps19+upps2+upps5+upps12
 factor3 =~ upps10+upps3+upps20+upps17
 factor4 =~ upps6+upps8+upps15+upps13
 factor5 =~ upps14+upps9+upps18+upps16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.948
#Tucker-Lewis Index (TLI)                       0.975       0.938
#Robust Comparative Fit Index (CFI)                         0.901
#Robust Tucker-Lewis Index (TLI)                            0.883
#RMSEA                                          0.046       0.052
#Robust RMSEA                                               0.065
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .417

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.930       0.935
#Tucker-Lewis Index (TLI)                       0.916       0.922
#Robust Comparative Fit Index (CFI)                         0.937
#Robust Tucker-Lewis Index (TLI)                            0.925
#RMSEA                                          0.045       0.041
#Robust RMSEA                                               0.042
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .352

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.567    0.061    9.344    0.000    0.567    0.567
#    factor3           0.188    0.078    2.399    0.016    0.188    0.188
#    factor4           0.267    0.077    3.486    0.000    0.267    0.267
#    factor5          -0.056    0.080   -0.700    0.484   -0.056   -0.056
#  factor2 ~~                                                            
#    factor3           0.388    0.066    5.927    0.000    0.388    0.388
#    factor4           0.428    0.062    6.945    0.000    0.428    0.428
#    factor5           0.136    0.080    1.703    0.089    0.136    0.136
#  factor3 ~~                                                            
#    factor4           0.703    0.055   12.847    0.000    0.703    0.703
#    factor5           0.355    0.074    4.831    0.000    0.355    0.355
#  factor4 ~~                                                            
#    factor5           0.131    0.080    1.649    0.099    0.131    0.131

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.782       0.782
#Tucker-Lewis Index (TLI)                       0.756       0.756
#Robust Comparative Fit Index (CFI)                         0.788
#Robust Tucker-Lewis Index (TLI)                            0.763
#RMSEA                                          0.077       0.072
#Robust RMSEA                                               0.075
#SRMR                                           0.139       0.139

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .393


# Bifactor model with 5 factors
BIFmodel= '
 factor1 =~ upps4+upps7+upps11+upps1
 factor2 =~ upps19+upps2+upps5+upps12
 factor3 =~ upps10+upps3+upps20+upps17
 factor4 =~ upps6+upps8+upps15+upps13
 factor5 =~ upps14+upps9+upps18+upps16
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.930       0.931
#Tucker-Lewis Index (TLI)                       0.911       0.913
#Robust Comparative Fit Index (CFI)                         0.935
#Robust Tucker-Lewis Index (TLI)                            0.917
#RMSEA                                          0.046       0.043
#Robust RMSEA                                               0.044
#SRMR                                           0.063       0.063

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .396

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3707430      0.8421053      0.8469713      0.5893664 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8807394 0.14684679 0.11926058 0.6551214 0.6107922 0.6410939 0.8066073
#factor2 0.6833975 0.15692517 0.31660250 0.7744931 0.5268801 0.6656627 0.8174222
#factor3 0.4316273 0.09885268 0.56837265 0.7675133 0.3306323 0.5162342 0.7055532
#factor4 0.3481691 0.07100978 0.65183095 0.7335894 0.2599113 0.4062681 0.6163909
#factor5 0.9118411 0.15562255 0.08815887 0.6720381 0.6262919 0.6623283 0.8148302
#general 0.3707430 0.37074303 0.37074303 0.8469713 0.5893664 0.8139192 0.8639555







################################################################
### UPPS-P Impulsive Behavior Scale 
### 
### Steiner & Frey (2021) https://psycnet.apa.org/record/2021-33516-001
### data available at https://osf.io/cufjv
###
### 4. 

library(readxl)
Steiner <- read_excel("Steiner.xlsx")
head(Steiner)
UPPS2 <- Steiner[,c(2,3,4)]
head(UPPS2)
library(tidyr)
mydata <- spread(UPPS2, quest_label, rating)
colnames(mydata)
mydata <- mydata[,c(2:46)]
mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
                      "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
                      "upps16","upps17","upps18","upps19","upps20","upps21","upps22","upps23","upps24","upps25",
                      "upps26","upps27","upps28","upps29","upps30","upps31","upps32","upps33","upps34","upps35",
                      "upps36","upps37","upps38","upps39","upps40","upps41","upps42","upps43","upps44","upps45")
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .93, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 11.36; eigenvalue 2 = 6.02

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 13.67; eigenvalue 2 = 7.08

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.123, RMSR=.16, TLI=.401

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.159, RMSR=.19, TLI=.342

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities


# Give solution with 4 factors (EFA based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.048, RMSR=.03, TLI=.909
#      MR1   MR2   MR3   MR4
#MR1  1.00  0.21 -0.30 -0.44
#MR2  0.21  1.00 -0.30  0.05
#MR3 -0.30 -0.30  1.00  0.28
#MR4 -0.44  0.05  0.28  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.083, RMSR=.03, TLI=.821
#      MR1   MR2   MR3   MR4
#MR1  1.00  0.21 -0.30 -0.44
#MR2  0.21  1.00 -0.31  0.06
#MR3 -0.30 -0.31  1.00  0.30
#MR4 -0.44  0.06  0.30  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20+upps21+upps22+upps23+upps24+
  upps25+upps26+upps27+upps28+upps29+upps30+upps31+upps32+upps33+upps34+upps35+upps36+
  upps37+upps38+upps39+upps40+upps41+upps42+upps43+upps44+upps45
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.789       0.628
#Tucker-Lewis Index (TLI)                       0.779       0.611
#Robust Comparative Fit Index (CFI)                         0.360
#Robust Tucker-Lewis Index (TLI)                            0.330
#RMSEA                                          0.234       0.146
#Robust RMSEA                                               0.163
#SRMR                                           0.208       0.208

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .404

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.433       0.426
#Tucker-Lewis Index (TLI)                       0.405       0.399
#Robust Comparative Fit Index (CFI)                         0.436
#Robust Tucker-Lewis Index (TLI)                            0.409
#RMSEA                                          0.124       0.115
#Robust RMSEA                                               0.123
#SRMR                                           0.162       0.162

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .252

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 perseverance =~ upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+
            upps10
 premeditation =~ upps11+upps12+upps13+upps14+upps15+upps16+upps17+upps18+upps19+
            upps20+upps20
 sensseeking =~ upps22+upps23+upps24+upps25+upps26+upps27+upps28+upps29+upps30+
            upps31+upps32+upps33
 urgency =~ upps44+upps35+upps36+upps37+upps38+upps39+upps40+upps41+upps42+
            upps43+upps44+upps45
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.960       0.886
#Tucker-Lewis Index (TLI)                       0.958       0.880
#Robust Comparative Fit Index (CFI)                         0.800
#Robust Tucker-Lewis Index (TLI)                            0.788
#RMSEA                                          0.102       0.083
#Robust RMSEA                                               0.093
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .625

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  perseverance ~~                                                       
#    premeditation     0.447    0.033   13.509    0.000    0.447    0.447
#    sensseeking       0.008    0.040    0.197    0.844    0.008    0.008
#    urgency          -0.624    0.025  -24.800    0.000   -0.624   -0.624
#  premeditation ~~                                                      
#    sensseeking      -0.477    0.031  -15.394    0.000   -0.477   -0.477
#    urgency          -0.512    0.030  -17.144    0.000   -0.512   -0.512
#  sensseeking ~~                                                        
#    urgency           0.288    0.037    7.881    0.000    0.288    0.288

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.849       0.852
#Tucker-Lewis Index (TLI)                       0.841       0.844
#Robust Comparative Fit Index (CFI)                         0.856
#Robust Tucker-Lewis Index (TLI)                            0.847
#RMSEA                                          0.066       0.060
#Robust RMSEA                                               0.064
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .511

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.817       0.820
#Tucker-Lewis Index (TLI)                       0.808       0.811
#Robust Comparative Fit Index (CFI)                         0.824
#Robust Tucker-Lewis Index (TLI)                            0.815
#RMSEA                                          0.072       0.066
#Robust RMSEA                                               0.070
#SRMR                                           0.184       0.184

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .497


# Bifactor model with 5 factors
BIFmodel= '
 perseverance =~ upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+
            upps10
 premeditation =~ upps11+upps12+upps13+upps14+upps15+upps16+upps17+upps18+upps19+
            upps20+upps20
 sensseeking =~ upps22+upps23+upps24+upps25+upps26+upps27+upps28+upps29+upps30+
            upps31+upps32+upps33
 urgency =~ upps44+upps35+upps36+upps37+upps38+upps39+upps40+upps41+upps42+
            upps43+upps44+upps45
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20+upps21+upps22+upps23+upps24+
  upps25+upps26+upps27+upps28+upps29+upps30+upps31+upps32+upps33+upps34+upps35+upps36+
  upps37+upps38+upps39+upps40+upps41+upps42+upps43+upps44+upps45
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.864       0.868
#Tucker-Lewis Index (TLI)                       0.851       0.855
#Robust Comparative Fit Index (CFI)                         0.871
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.062       0.056
#Robust RMSEA                                               0.060
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .507

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4347508      0.7868687      0.8552367      0.1297756 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#perseverance  0.6057875 0.13479422 0.39421247 0.7695330 0.62772567 0.8374037 0.9286351
#premeditation 0.7376418 0.15295183 0.26235817 0.8882885 0.67685575 0.8488587 0.9284864
#sensseeking   0.9112803 0.24337440 0.08871973 0.9194132 0.85362056 0.9139212 0.9583098
#urgency       0.1245172 0.03412877 0.87548276 0.9054308 0.06659812 0.5084756 0.7845618
#general       0.4347508 0.43475078 0.43475078 0.8552367 0.12977561 0.9466040 0.9612212






################################################################
### UPPS Impulsive Behavior Scale 
### 
### Racine et al. (2022) https://osf.io/preprints/psyarxiv/3fbm8
### data available at https://osf.io/mxjn9
###
### 5. 

library(haven)
Racine_EMA_Data <- read_sav("Racine_EMA_Data.sav")
colnames(Racine_EMA_Data)
mydata <- as.data.frame(Racine_EMA_Data[,74:132])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
                      "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
                      "upps16","upps17","upps18","upps19","upps20","upps21","upps22","upps23","upps24","upps25",
                      "upps26","upps27","upps28","upps29","upps30","upps31","upps32","upps33","upps34","upps35",
                      "upps36","upps37","upps38","upps39","upps40","upps41","upps42","upps43","upps44","upps45",
                      "upps46","upps47","upps48","upps49","upps50","upps51","upps52","upps53","upps54","upps55",
                      "upps56","upps57","upps58","upps59")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .94, omega T = .96

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 15 factors and 10 components
# Eigenvalue 1 = 14.93; eigenvalue 2 = 5.70

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 15 factors and 10 components
# Eigenvalue 1 = 17.71; eigenvalue 2 = 6.67

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.131, RMSR=.14, TLI=.326

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.305, RMSR=.16, TLI=.082

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 5 factors (theory-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.091, RMSR=.04, TLI=.672
#      MR1   MR4   MR2   MR3   MR5
#MR1  1.00  0.43  0.23 -0.30 -0.25
#MR4  0.43  1.00  0.04 -0.27 -0.24
#MR2  0.23  0.04  1.00 -0.14  0.09
#MR3 -0.30 -0.27 -0.14  1.00  0.38
#MR5 -0.25 -0.24  0.09  0.38  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.302, RMSR=.05, TLI=.099
#      MR1   MR4   MR2   MR3   MR5
#MR1  1.00  0.45  0.25 -0.31 -0.26
#MR4  0.45  1.00  0.03 -0.25 -0.22
#MR2  0.25  0.03  1.00 -0.14  0.09
#MR3 -0.31 -0.25 -0.14  1.00  0.38
#MR5 -0.26 -0.22  0.09  0.38  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20+upps21+upps22+upps23+upps24+
  upps25+upps26+upps27+upps28+upps29+upps30+upps31+upps32+upps33+upps34+upps35+upps36+
  upps37+upps38+upps39+upps40+upps41+upps42+upps43+upps44+upps45+upps46+upps47+upps48+
  upps49+upps50+upps51+upps52+upps53+upps54+upps55+upps56+upps57+upps58+upps59
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.812       0.611
#Tucker-Lewis Index (TLI)                       0.805       0.597
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.196       0.115
#Robust RMSEA                                                  NA
#SRMR                                           0.167       0.167

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .342

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.349       0.342
#Tucker-Lewis Index (TLI)                       0.326       0.319
#Robust Comparative Fit Index (CFI)                         0.349
#Robust Tucker-Lewis Index (TLI)                            0.326
#RMSEA                                          0.131       0.130
#Robust RMSEA                                               0.131
#SRMR                                           0.136       0.136

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .262

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 negurg =~ upps2+upps7+upps12+upps17+upps22+upps29+upps34+upps39+upps44+
            upps51+upps54+upps58
 premed =~ upps1+upps6+upps11+upps16+upps21+upps28+upps33+upps38+upps43+
            upps48+upps59
 perseve =~ upps4+upps9+upps14+upps19+upps24+upps27+upps32+upps37+upps42+
            upps47
 sensseek =~ upps3+upps8+upps13+upps18+upps23+upps9+upps26+upps31+upps36+
             upps41+upps46+upps52+upps56
 posurg =~ upps5+upps10+upps15+upps20+upps30+upps35+upps40+upps45+upps50
            +upps53+upps55+upps57+upps59
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.947       0.880
#Tucker-Lewis Index (TLI)                       0.945       0.875
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.103       0.065
#Robust RMSEA                                                  NA
#SRMR                                           0.099       0.099

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .597

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.648       0.642
#Tucker-Lewis Index (TLI)                       0.632       0.626
#Robust Comparative Fit Index (CFI)                         0.648
#Robust Tucker-Lewis Index (TLI)                            0.632
#RMSEA                                          0.098       0.097
#Robust RMSEA                                               0.098
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .496

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negurg ~~                                                             
#    premed           -0.524    0.008  -62.383    0.000   -0.524   -0.524
#    perseve          -0.520    0.008  -62.785    0.000   -0.520   -0.520
#    sensseek          0.098    0.010    9.798    0.000    0.098    0.098
#    posurg            0.611    0.007   88.675    0.000    0.611    0.611
#  premed ~~                                                             
#    perseve           0.604    0.008   71.658    0.000    0.604    0.604
#    sensseek         -0.251    0.010  -25.632    0.000   -0.251   -0.251
#    posurg           -0.437    0.009  -46.786    0.000   -0.437   -0.437
#  perseve ~~                                                            
#    sensseek          0.079    0.010    8.134    0.000    0.079    0.079
#    posurg           -0.377    0.010  -38.011    0.000   -0.377   -0.377
#  sensseek ~~                                                           
#    posurg            0.293    0.009   32.061    0.000    0.293    0.293

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.616       0.610
#Tucker-Lewis Index (TLI)                       0.601       0.595
#Robust Comparative Fit Index (CFI)                         0.616
#Robust Tucker-Lewis Index (TLI)                            0.601
#RMSEA                                          0.102       0.101
#Robust RMSEA                                               0.102
#SRMR                                           0.201       0.201

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .490


# Bifactor model with 5 factors
BIFmodel= '
 negurg =~ upps2+upps7+upps12+upps17+upps22+upps29+upps34+upps39+upps44+
            upps51+upps54+upps58
 premed =~ upps1+upps6+upps11+upps16+upps21+upps28+upps33+upps38+upps43+
            upps48+upps59
 perseve =~ upps4+upps9+upps14+upps19+upps24+upps27+upps32+upps37+upps42+
            upps47
 sensseek =~ upps3+upps8+upps13+upps18+upps23+upps9+upps26+upps31+upps36+
             upps41+upps46+upps52+upps56
 posurg =~ upps5+upps10+upps15+upps20+upps30+upps35+upps40+upps45+upps50
            +upps53+upps55+upps57+upps59
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20+upps21+upps22+upps23+upps24+
  upps25+upps26+upps27+upps28+upps29+upps30+upps31+upps32+upps33+upps34+upps35+upps36+
  upps37+upps38+upps39+upps40+upps41+upps42+upps43+upps44+upps45+upps46+upps47+upps48+
  upps49+upps50+upps51+upps52+upps53+upps54+upps55+upps56+upps57+upps58+upps59
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.655       0.647
#Tucker-Lewis Index (TLI)                       0.630       0.621
#Robust Comparative Fit Index (CFI)                         0.655
#Robust Tucker-Lewis Index (TLI)                            0.630
#RMSEA                                          0.097       0.097
#Robust RMSEA                                               0.097
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .490

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  negurg =~                                                             
#    upps2             0.057    0.012    4.894    0.000    0.057    0.058
#    upps7             0.122    0.012   10.488    0.000    0.122    0.124
#    upps12            0.119    0.010   12.322    0.000    0.119    0.126
#    upps17           -0.001    0.016   -0.077    0.938   -0.001   -0.001
#    upps22            0.045    0.015    3.018    0.003    0.045    0.044
#    upps29            0.287    0.013   22.944    0.000    0.287    0.307
#    upps34            0.535    0.016   33.948    0.000    0.535    0.516
#    upps39            0.162    0.010   16.542    0.000    0.162    0.169
#    upps44            0.388    0.017   23.380    0.000    0.388    0.355
#    upps51            0.776    0.014   54.457    0.000    0.776    0.770
#    upps54           -0.166    0.009  -17.950    0.000   -0.166   -0.189
#    upps58            0.186    0.011   17.377    0.000    0.186    0.202
#  premed =~                                                             
#    upps1             0.301    0.008   36.304    0.000    0.301    0.362
#    upps6             0.431    0.008   55.464    0.000    0.431    0.507
#    upps11            0.145    0.010   14.922    0.000    0.145    0.146
#    upps16            0.519    0.006   88.089    0.000    0.519    0.663
#    upps21            0.412    0.008   48.453    0.000    0.412    0.461
#    upps28            0.444    0.008   56.833    0.000    0.444    0.542
#    upps33            0.423    0.007   64.228    0.000    0.423    0.553
#    upps38            0.453    0.007   60.637    0.000    0.453    0.548
#    upps43            0.361    0.007   53.800    0.000    0.361    0.498
#    upps48            0.648    0.006  116.445    0.000    0.648    0.811
#    upps59           -0.064    0.009   -6.806    0.000   -0.064   -0.071
#  perseve =~                                                            
#    upps4             0.520    0.007   72.937    0.000    0.520    0.638
#    upps9            -0.192    0.007  -26.712    0.000   -0.192   -0.232
#    upps14            0.484    0.009   56.067    0.000    0.484    0.568
#    upps19            0.369    0.007   54.595    0.000    0.369    0.460
#    upps24            0.202    0.010   21.148    0.000    0.202    0.209
#    upps27            0.660    0.006  107.435    0.000    0.660    0.786
#    upps32            0.295    0.009   32.207    0.000    0.295    0.312
#    upps37            0.516    0.007   70.832    0.000    0.516    0.629
#    upps42            0.622    0.006  110.206    0.000    0.622    0.723
#    upps47           -0.201    0.008  -24.249    0.000   -0.201   -0.199
#  sensseek =~                                                           
#    upps3             0.575    0.007   80.415    0.000    0.575    0.694
#    upps8             0.432    0.008   50.845    0.000    0.432    0.435
#    upps13            0.343    0.009   38.888    0.000    0.343    0.328
#    upps18            0.685    0.008   86.645    0.000    0.685    0.678
#    upps23            0.608    0.007   89.671    0.000    0.608    0.705
#    upps9            -0.243    0.007  -35.631    0.000   -0.243   -0.293
#    upps26            0.861    0.007  132.173    0.000    0.861    0.761
#    upps31            0.664    0.007   96.470    0.000    0.664    0.755
#    upps36            0.586    0.010   59.156    0.000    0.586    0.513
#    upps41            0.678    0.007  103.341    0.000    0.678    0.755
#    upps46            0.689    0.008   87.309    0.000    0.689    0.673
#    upps52            0.658    0.009   76.066    0.000    0.658    0.610
#    upps56            0.613    0.009   66.072    0.000    0.613    0.540
#  posurg =~                                                             
#    upps5             0.425    0.007   56.786    0.000    0.425    0.497
#    upps10            0.516    0.006   86.561    0.000    0.516    0.675
#    upps15            0.407    0.006   68.279    0.000    0.407    0.571
#    upps20            0.479    0.007   73.318    0.000    0.479    0.582
#    upps30            0.603    0.008   77.974    0.000    0.603    0.717
#    upps35            0.482    0.007   72.072    0.000    0.482    0.590
#    upps40            0.491    0.008   62.826    0.000    0.491    0.654
#    upps45            0.525    0.007   79.724    0.000    0.525    0.593
#    upps50            0.418    0.008   51.134    0.000    0.418    0.458
#    upps53            0.410    0.008   48.949    0.000    0.410    0.463
#    upps55            0.333    0.008   39.177    0.000    0.333    0.384
#    upps57            0.014    0.008    1.732    0.083    0.014    0.017
#    upps59            0.319    0.008   37.557    0.000    0.319    0.351
#  general =~                                                            
#    upps1             0.246    0.008   29.912    0.000    0.246    0.296
#    upps2            -0.764    0.006 -127.561    0.000   -0.764   -0.780
#    upps3            -0.075    0.009   -8.496    0.000   -0.075   -0.091
#    upps4             0.274    0.008   33.618    0.000    0.274    0.336
#    upps5            -0.435    0.008  -56.758    0.000   -0.435   -0.510
#    upps6             0.473    0.008   59.804    0.000    0.473    0.556
#    upps7            -0.638    0.007  -87.747    0.000   -0.638   -0.645
#    upps8            -0.350    0.009  -38.608    0.000   -0.350   -0.352
#    upps9            -0.478    0.008  -63.513    0.000   -0.478   -0.576
#    upps10           -0.425    0.007  -59.263    0.000   -0.425   -0.556
#    upps11            0.262    0.010   25.713    0.000    0.262    0.264
#    upps12           -0.434    0.008  -54.192    0.000   -0.434   -0.462
#    upps13            0.011    0.011    0.937    0.349    0.011    0.010
#    upps14            0.191    0.009   21.439    0.000    0.191    0.225
#    upps15           -0.363    0.007  -55.757    0.000   -0.363   -0.510
#    upps16            0.338    0.009   39.665    0.000    0.338    0.432
#    upps17           -0.794    0.007 -114.602    0.000   -0.794   -0.751
#    upps18            0.019    0.011    1.726    0.084    0.019    0.019
#    upps19           -0.024    0.008   -2.916    0.004   -0.024   -0.030
#    upps20           -0.469    0.008  -59.346    0.000   -0.469   -0.570
#    upps21            0.079    0.010    8.134    0.000    0.079    0.088
#    upps22           -0.705    0.008  -90.634    0.000   -0.705   -0.694
#    upps23           -0.154    0.009  -17.927    0.000   -0.154   -0.178
#    upps24            0.536    0.008   65.890    0.000    0.536    0.555
#    upps25           -0.479    0.009  -54.142    0.000   -0.479   -0.526
#    upps26           -0.200    0.012  -17.381    0.000   -0.200   -0.177
#    upps27            0.400    0.008   51.366    0.000    0.400    0.476
#    upps28            0.410    0.009   47.498    0.000    0.410    0.501
#    upps29           -0.648    0.008  -81.021    0.000   -0.648   -0.694
#    upps30           -0.379    0.008  -46.607    0.000   -0.379   -0.450
#    upps31           -0.065    0.010   -6.685    0.000   -0.065   -0.074
#    upps32            0.451    0.009   48.676    0.000    0.451    0.477
#    upps33            0.417    0.008   54.939    0.000    0.417    0.546
#    upps34           -0.597    0.011  -53.148    0.000   -0.597   -0.576
#    upps35           -0.354    0.008  -43.695    0.000   -0.354   -0.434
#    upps36           -0.092    0.012   -7.994    0.000   -0.092   -0.081
#    upps37            0.360    0.009   40.904    0.000    0.360    0.438
#    upps38            0.288    0.009   31.296    0.000    0.288    0.349
#    upps39           -0.582    0.008  -77.093    0.000   -0.582   -0.605
#    upps40           -0.382    0.007  -51.871    0.000   -0.382   -0.508
#    upps41           -0.132    0.010  -13.827    0.000   -0.132   -0.147
#    upps42            0.414    0.008   48.972    0.000    0.414    0.482
#    upps43            0.189    0.009   21.989    0.000    0.189    0.260
#    upps44           -0.821    0.008 -105.170    0.000   -0.821   -0.751
#    upps45           -0.456    0.009  -52.964    0.000   -0.456   -0.515
#    upps46           -0.111    0.010  -10.784    0.000   -0.111   -0.109
#    upps47           -0.604    0.008  -73.115    0.000   -0.604   -0.597
#    upps48            0.330    0.009   35.282    0.000    0.330    0.413
#    upps49            0.325    0.009   34.750    0.000    0.325    0.391
#    upps50           -0.487    0.008  -57.254    0.000   -0.487   -0.533
#    upps51           -0.555    0.010  -54.472    0.000   -0.555   -0.551
#    upps52            0.007    0.011    0.623    0.533    0.007    0.006
#    upps53           -0.448    0.009  -49.759    0.000   -0.448   -0.506
#    upps54            0.508    0.008   67.099    0.000    0.508    0.580
#    upps55           -0.332    0.008  -39.312    0.000   -0.332   -0.383
#    upps56           -0.286    0.011  -25.316    0.000   -0.286   -0.251
#    upps57           -0.354    0.008  -43.571    0.000   -0.354   -0.419
#    upps58           -0.680    0.008  -89.660    0.000   -0.680   -0.740
#    upps59           -0.316    0.009  -33.471    0.000   -0.316   -0.348

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.4497309      0.8772134      0.3195827 
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#negurg   0.1896689 0.04353229 0.8103311 0.9009547 0.1095031 0.6884196 0.9047109
#premed   0.6128224 0.10247342 0.3608707 0.8529424 0.5873566 0.8366128 0.9312152
#perseve  0.5606744 0.09672947 0.4216139 0.7812862 0.6425402 0.8345957 0.9345224
#sensseek 0.8770015 0.17542364 0.1134179 0.8819412 0.8171183 0.9041178 0.9521815
#posurg   0.5479484 0.13211029 0.4513109 0.9291831 0.4871415 0.8549487 0.9306560
#general  0.4497309 0.44973089 0.4497309 0.8772134 0.3195827 0.9524484 0.9607406






################################################################
### UPPS Impulsive Behavior Scale 
### 
### Flayelle et al. (2022) https://www.sciencedirect.com/science/article/pii/S0736585322001137
### data available at https://osf.io/unjzq/
###
### 6. 

library(readxl)
Flayelle <- read_excel("Flayelle.xlsx")
colnames(Flayelle)
mydata <- as.data.frame(Flayelle[,27:46])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
                      "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
                      "upps16","upps17","upps18","upps19","upps20")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .83, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 4.37; eigenvalue 2 = 2.44

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 5.10; eigenvalue 2 = 2.84

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.177, RMSR=.17, TLI=.305

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.214, RMSR=.20, TLI=.270

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.042, RMSR=.01, TLI=.962
#      MR2  MR1  MR4   MR3  MR5
#MR2  1.00 0.06 0.37 -0.05 0.14
#MR1  0.06 1.00 0.27  0.12 0.57
#MR4  0.37 0.27 1.00  0.12 0.37
#MR3 -0.05 0.12 0.12  1.00 0.33
#MR5  0.14 0.57 0.37  0.33 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.65, RMSEA=.059, RMSR=.02, TLI=.945
#      MR2  MR4  MR1   MR3  MR5
#MR2  1.00 0.39 0.06 -0.06 0.13
#MR4  0.39 1.00 0.28  0.12 0.35
#MR1  0.06 0.28 1.00  0.12 0.56
#MR3 -0.06 0.12 0.12  1.00 0.34
#MR5  0.13 0.35 0.56  0.34 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.715       0.572
#Tucker-Lewis Index (TLI)                       0.682       0.522
#Robust Comparative Fit Index (CFI)                         0.254
#Robust Tucker-Lewis Index (TLI)                            0.166
#RMSEA                                          0.283       0.230
#Robust RMSEA                                               0.229
#SRMR                                           0.221       0.221

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .355

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.385       0.357
#Tucker-Lewis Index (TLI)                       0.312       0.282
#Robust Comparative Fit Index (CFI)                         0.385
#Robust Tucker-Lewis Index (TLI)                            0.313
#RMSEA                                          0.176       0.163
#Robust RMSEA                                               0.176
#SRMR                                           0.165       0.165

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .180

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1 =~ upps13+upps19+upps1+upps6
 factor2 =~ upps8+upps11+upps16+upps5
 factor3 =~ upps9+upps18+upps3+upps14
 factor4 =~ upps7+upps17+upps12+upps4
 factor5 =~ upps10+upps15+upps2+upps20
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.967
#Tucker-Lewis Index (TLI)                       0.987       0.961
#Robust Comparative Fit Index (CFI)                         0.936
#Robust Tucker-Lewis Index (TLI)                            0.924
#RMSEA                                          0.058       0.066
#SRMR                                           0.049       0.049

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .640

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.445    0.015   30.625    0.000    0.445    0.445
#    factor3           0.127    0.017    7.281    0.000    0.127    0.127
#    factor4           0.346    0.015   22.698    0.000    0.346    0.346
#    factor5           0.446    0.015   29.944    0.000    0.446    0.446
#  factor2 ~~                                                            
#    factor3          -0.056    0.017   -3.204    0.001   -0.056   -0.056
#    factor4           0.083    0.017    4.872    0.000    0.083    0.083
#    factor5           0.173    0.017   10.222    0.000    0.173    0.173
#  factor3 ~~                                                            
#    factor4           0.169    0.017   10.186    0.000    0.169    0.169
#    factor5           0.365    0.016   23.181    0.000    0.365    0.365
#  factor4 ~~                                                            
#    factor5           0.703    0.010   70.486    0.000    0.703    0.703

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.949       0.948
#Tucker-Lewis Index (TLI)                       0.940       0.938
#Robust Comparative Fit Index (CFI)                         0.950
#Robust Tucker-Lewis Index (TLI)                            0.941
#RMSEA                                          0.052       0.048
#Robust RMSEA                                               0.052
#SRMR                                           0.043       0.043

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .536

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.406    0.019   21.666    0.000    0.406    0.406
#    factor3           0.129    0.020    6.325    0.000    0.129    0.129
#    factor4           0.343    0.019   18.191    0.000    0.343    0.343
#    factor5           0.457    0.018   25.387    0.000    0.457    0.457
#  factor2 ~~                                                            
#    factor3          -0.053    0.019   -2.712    0.007   -0.053   -0.053
#    factor4           0.081    0.019    4.343    0.000    0.081    0.081
#    factor5           0.180    0.019    9.636    0.000    0.180    0.180
#  factor3 ~~                                                            
#    factor4           0.172    0.020    8.712    0.000    0.172    0.172
#    factor5           0.328    0.019   17.082    0.000    0.328    0.328
#  factor4 ~~                                                            
#    factor5           0.708    0.013   53.268    0.000    0.708    0.708

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.857       0.851
#Tucker-Lewis Index (TLI)                       0.840       0.833
#Robust Comparative Fit Index (CFI)                         0.857
#Robust Tucker-Lewis Index (TLI)                            0.840
#RMSEA                                          0.085       0.079
#Robust RMSEA                                               0.085
#SRMR                                           0.159       0.159

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .540


# Bifactor model with 5 factors
BIFmodel= '
 factor1 =~ upps13+upps19+upps1+upps6
 factor2 =~ upps8+upps11+upps16+upps5
 factor3 =~ upps9+upps18+upps3+upps14
 factor4 =~ upps7+upps17+upps12+upps4
 factor5 =~ upps10+upps15+upps2+upps20
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.954       0.952
#Tucker-Lewis Index (TLI)                       0.942       0.940
#Robust Comparative Fit Index (CFI)                         0.955
#Robust Tucker-Lewis Index (TLI)                            0.943
#RMSEA                                          0.051       0.047
#Robust RMSEA                                               0.051
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .570

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3131273      0.8421053      0.9105059      0.5932409 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.7115898 0.13467816 0.28841022 0.8224277 0.5848955 0.7191184 0.8550491
#factor2 0.9453497 0.22061453 0.05465029 0.8870163 0.8389231 0.8788241 0.9390712
#factor3 0.9145596 0.17179862 0.08544040 0.8204495 0.7520485 0.7935452 0.8926103
#factor4 0.4766346 0.10347088 0.52336544 0.8587498 0.4033321 0.6587066 0.8143272
#factor5 0.3265639 0.05631047 0.67343606 0.7751751 0.2564940 0.4555199 0.6609873
#general 0.3131273 0.31312734 0.31312734 0.9105059 0.5932409 0.8577259 0.8945531







################################################################
### UPPS Impulsive Behavior Scale 
### 
### Zsila et al. (2020) https://link.springer.com/article/10.1007/s12144-017-9773-7
### data from author Beata Bothe 
###
### 7. 

library(haven)
Bothe <- read_sav("Bothe.sav")
colnames(Bothe)
mydata <- as.data.frame(Bothe[,1:20])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("upps1","upps2","upps3","upps4","upps5","upps6","upps7","upps8",
                      "upps9","upps10","upps11","upps12","upps13","upps14","upps15",
                      "upps16","upps17","upps18","upps19","upps20")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .85, omega T = .90

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 4 components
# Eigenvalue 1 = 4.91; eigenvalue 2 = 2.45

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 4 components
# Eigenvalue 1 = 5.74; eigenvalue 2 = 2.88

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.157, RMSR=.15, TLI=.423

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.29, RMSEA=.191, RMSR=.18, TLI=.383

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.042, RMSR=.02, TLI=.959
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.20 -0.43  0.18  0.54
#MR2 -0.20  1.00  0.48  0.10 -0.11
#MR4 -0.43  0.48  1.00 -0.06 -0.21
#MR3  0.18  0.10 -0.06  1.00  0.50
#MR5  0.54 -0.11 -0.21  0.50  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.62, RMSEA=.057, RMSR=.02, TLI=.945
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.21 -0.44  0.18  0.53
#MR2 -0.21  1.00  0.50  0.11 -0.11
#MR4 -0.44  0.50  1.00 -0.06 -0.20
#MR3  0.18  0.11 -0.06  1.00  0.49
#MR5  0.53 -0.11 -0.20  0.49  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.777       0.597
#Tucker-Lewis Index (TLI)                       0.751       0.549
#Robust Comparative Fit Index (CFI)                         0.394
#Robust Tucker-Lewis Index (TLI)                            0.322
#RMSEA                                          0.236       0.209
#Robust RMSEA                                               0.200
#SRMR                                           0.188       0.188

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .38

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.490       0.473
#Tucker-Lewis Index (TLI)                       0.430       0.411
#Robust Comparative Fit Index (CFI)                         0.490
#Robust Tucker-Lewis Index (TLI)                            0.430
#RMSEA                                          0.156       0.144
#Robust RMSEA                                               0.156
#SRMR                                           0.150       0.150

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .203

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1 =~ upps13+upps19+upps1+upps6
 factor2 =~ upps8+upps11+upps16+upps5
 factor3 =~ upps9+upps18+upps3+upps14
 factor4 =~ upps7+upps17+upps12+upps4
 factor5 =~ upps10+upps15+upps2+upps20
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.974       0.929
#Tucker-Lewis Index (TLI)                       0.969       0.915
#Robust Comparative Fit Index (CFI)                         0.900
#Robust Tucker-Lewis Index (TLI)                            0.882
#RMSEA                                          0.084       0.091
#Robust RMSEA                                               0.084
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .570

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.617    0.006   99.297    0.000    0.617    0.617
#    factor3          -0.060    0.009   -6.437    0.000   -0.060   -0.060
#    factor4          -0.450    0.007  -60.007    0.000   -0.450   -0.450
#    factor5          -0.513    0.007  -71.915    0.000   -0.513   -0.513
#  factor2 ~~                                                            
#    factor3           0.100    0.009   10.774    0.000    0.100    0.100
#    factor4          -0.229    0.009  -26.125    0.000   -0.229   -0.229
#    factor5          -0.241    0.009  -27.568    0.000   -0.241   -0.241
#  factor3 ~~                                                            
#    factor4           0.220    0.009   24.806    0.000    0.220    0.220
#    factor5           0.498    0.008   63.931    0.000    0.498    0.498
#  factor4 ~~                                                            
#    factor5           0.859    0.004  231.483    0.000    0.859    0.859

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.920       0.917
#Tucker-Lewis Index (TLI)                       0.906       0.901
#Robust Comparative Fit Index (CFI)                         0.921
#Robust Tucker-Lewis Index (TLI)                            0.906
#RMSEA                                          0.064       0.059
#Robust RMSEA                                               0.063
#SRMR                                           0.063       0.063

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .498

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.561    0.009   61.787    0.000    0.561    0.561
#    factor3          -0.048    0.012   -4.206    0.000   -0.048   -0.048
#    factor4          -0.448    0.009  -51.084    0.000   -0.448   -0.448
#    factor5          -0.526    0.009  -59.692    0.000   -0.526   -0.526
#  factor2 ~~                                                            
#    factor3           0.095    0.011    8.863    0.000    0.095    0.095
#    factor4          -0.217    0.009  -22.931    0.000   -0.217   -0.217
#    factor5          -0.255    0.010  -26.446    0.000   -0.255   -0.255
#  factor3 ~~                                                            
#    factor4           0.192    0.011   17.270    0.000    0.192    0.192
#    factor5           0.410    0.011   36.267    0.000    0.410    0.410
#  factor4 ~~                                                            
#    factor5           0.867    0.005  170.254    0.000    0.867    0.867

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.761       0.751
#Tucker-Lewis Index (TLI)                       0.733       0.722
#Robust Comparative Fit Index (CFI)                         0.761
#Robust Tucker-Lewis Index (TLI)                            0.733
#RMSEA                                          0.107       0.099
#SRMR                                           0.197       0.197

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .501


# Bifactor model with 5 factors
BIFmodel= '
 factor1 =~ upps13+upps19+upps1+upps6
 factor2 =~ upps8+upps11+upps16+upps5
 factor3 =~ upps9+upps18+upps3+upps14
 factor4 =~ upps7+upps17+upps12+upps4
 factor5 =~ upps10+upps15+upps2+upps20
 general =~  upps1+upps2+upps3+upps4+upps5+upps6+upps7+upps8+upps9+upps10+upps11+upps12+
  upps13+upps14+upps15+upps16+upps17+upps18+upps19+upps20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.920       0.915
#Tucker-Lewis Index (TLI)                       0.898       0.893
#Robust Comparative Fit Index (CFI)                         0.920
#Robust Tucker-Lewis Index (TLI)                            0.899
#RMSEA                                          0.066       0.061
#Robust RMSEA                                               0.066
#SRMR                                           0.087       0.087

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .532

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3987054      0.8421053      0.7875835      0.2244581 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.6794566 0.13432538 0.3205434 0.8188166 0.5559262 0.7021283 0.8547353
#factor2 0.8956142 0.18861481 0.1043858 0.8375722 0.7492966 0.8303112 0.9191965
#factor3 0.8793488 0.15725388 0.1206512 0.7827360 0.7026824 0.7565483 0.8730523
#factor4 0.3036267 0.06484458 0.6963733 0.8390542 0.2435282 0.4845505 0.7273966
#factor5 0.2822539 0.05625595 0.7177461 0.7911980 0.1644229 0.4830687 0.7304242
#general 0.3987054 0.39870541 0.3987054 0.7875835 0.2244581 0.8900960 0.9260850





################################################################
### Barratt Impulsiveness Scale (BIS-11)
### 
###  Assumed model = 3 factors based on Kapitány-Fövény et al (2020)
###  https://www.researchgate.net/publication/342639896_The_21-item_Barratt_Impulsiveness_Scale_Revised_BIS-R-21_An_alternative_three-factor_model
###
###




################################################################
### Barratt Impulsiveness Scale (BIS-11)
### 
### Littrell (2020) https://www.sciencedirect.com/science/article/abs/pii/S0191886919306105
### data available at https://osf.io/cja5k
###
### 7. 

### Not included because too many calculation problems

Littrell <- read.csv("Littrell.csv")
colnames(Littrell)
mydata <- as.data.frame(Littrell[,44:73])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .89, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 4 components
# Eigenvalue 1 = 6.93; eigenvalue 2 = 2.87

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 8.71; eigenvalue 2 = 3.33

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.23, RMSEA=.121, RMSR=.13, TLI=.437

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.29, RMSEA=.283, RMSR=.15, TLI=.133

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al., 2020)
library(psych)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.083, RMSR=.07, TLI=.730
#      MR1   MR2   MR3
#MR1  1.00 -0.20  0.32
#MR2 -0.20  1.00 -0.23
#MR3  0.32 -0.23  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.279, RMSR=.08, TLI=.151
#      MR1   MR2   MR3
#MR1  1.00 -0.22  0.36
#MR2 -0.22  1.00 -0.25
#MR3  0.36 -0.25  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15+BIS_16+BIS_17+BIS_18+BIS_19+BIS_20+
            BIS_21+BIS_22+BIS_23+BIS_24+BIS_25+BIS_26+BIS_27+BIS_28+BIS_29+BIS_30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
#warning
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.473       0.443
#Tucker-Lewis Index (TLI)                       0.434       0.402
#Robust Comparative Fit Index (CFI)                         0.471
#Robust Tucker-Lewis Index (TLI)                            0.432
#RMSEA                                          0.126       0.123
#Robust RMSEA                                               0.125
#SRMR                                           0.126       0.126

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .219

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (based on Kapitány-Fövény et al.) 
EGAmodel= '
 factor1 =~ BIS_9+BIS_12+BIS_8+BIS_13+BIS_20+BIS_1+BIS_7+BIS_15+BIS_30+BIS_10+
            BIS_29
 factor2 =~ BIS_19+BIS_17+BIS_14+BIS_18+BIS_2+BIS_4+BIS_5+BIS_6+BIS_3+BIS_27+BIS_23
 factor3 =~ BIS_24+BIS_21+BIS_22+BIS_26+BIS_25+BIS_28+BIS_11+BIS_16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.636       0.630
#Tucker-Lewis Index (TLI)                       0.606       0.600
#Robust Comparative Fit Index (CFI)                         0.641
#Robust Tucker-Lewis Index (TLI)                            0.611
#RMSEA                                          0.105       0.100
#Robust RMSEA                                               0.103
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .340

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#  factor2          -0.488    0.085   -5.724    0.000   -0.488   -0.488
#  factor3          -0.370    0.098   -3.766    0.000   -0.370   -0.370
#factor2 ~~                                                            
#  factor3           0.827    0.061   13.463    0.000    0.827    0.827

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.582       0.571
#Tucker-Lewis Index (TLI)                       0.551       0.539
#Robust Comparative Fit Index (CFI)                         0.585
#Robust Tucker-Lewis Index (TLI)                            0.554
#RMSEA                                          0.112       0.108
#Robust RMSEA                                               0.111
#SRMR                                           0.186       0.186

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .250


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ BIS_12+BIS_8+BIS_9+BIS_13+BIS_20+BIS_1+BIS_7+BIS_15+BIS_30+BIS_10+
            BIS_29
 factor2 =~ BIS_19+BIS_17+BIS_14+BIS_18+BIS_2+BIS_4+BIS_5+BIS_6+BIS_3+BIS_27+BIS_23
 factor3 =~ BIS_28+BIS_24+BIS_21+BIS_22+BIS_26+BIS_25+BIS_11+BIS_16
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15+BIS_16+BIS_17+BIS_18+BIS_19+BIS_20+
            BIS_21+BIS_22+BIS_23+BIS_24+BIS_25+BIS_26+BIS_27+BIS_28+BIS_29+BIS_30
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning and no estimates





################################################################
### Barratt Impulsiveness Scale (BIS-11)
### 
### Haines et al (2020) https://journals.sagepub.com/doi/full/10.1177/2167702620929636
### data available at https://osf.io/4a3qt
###
### 8. 

Haines <- read.csv("Haines.csv")
colnames(Haines)
mydata <- as.data.frame(Haines[,26:55])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .89, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 5 components
# Eigenvalue 1 = 7.24; eigenvalue 2 = 2.23

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 6 components
# Eigenvalue 1 = 8.50; eigenvalue 2 = 2.55

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.24, RMSEA=.111, RMSR=.10, TLI=.523

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.28, RMSEA=.139, RMSR=.12, TLI=.469

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al.)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.08, RMSR=.05, TLI=.75
#      MR2   MR1   MR3
#MR2  1.00 -0.39 -0.23
#MR1 -0.39  1.00  0.29
#MR3 -0.23  0.29  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.108, RMSR=.06, TLI=.681
#      MR1   MR2   MR3
#MR1  1.00 -0.38  0.27
#MR2 -0.38  1.00 -0.24
#MR3  0.27 -0.24  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15+BIS_16+BIS_17+BIS_18+BIS_19+BIS_20+
            BIS_21+BIS_22+BIS_23+BIS_24+BIS_25+BIS_26+BIS_27+BIS_28+BIS_29+BIS_30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.861       0.686
#Tucker-Lewis Index (TLI)                       0.851       0.663
#Robust Comparative Fit Index (CFI)                         0.491
#Robust Tucker-Lewis Index (TLI)                            0.453
#RMSEA                                          0.149       0.135
#Robust RMSEA                                               0.142
#SRMR                                           0.128       0.128

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .336

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.557       0.552
#Tucker-Lewis Index (TLI)                       0.524       0.519
#Robust Comparative Fit Index (CFI)                         0.559
#Robust Tucker-Lewis Index (TLI)                            0.527
#RMSEA                                          0.112       0.104
#Robust RMSEA                                               0.111
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .230

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 factor1 =~ BIS_12+BIS_8+BIS_9+BIS_13+BIS_20+BIS_1+BIS_7+BIS_15+BIS_30+BIS_10+
            BIS_29
 factor2 =~ BIS_19+BIS_17+BIS_14+BIS_18+BIS_2+BIS_4+BIS_5+BIS_6+BIS_3+BIS_27+BIS_23
 factor3 =~ BIS_28+BIS_24+BIS_21+BIS_22+BIS_26+BIS_25+BIS_11+BIS_16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.927       0.824
#Tucker-Lewis Index (TLI)                       0.921       0.809
#Robust Comparative Fit Index (CFI)                         0.639
#Robust Tucker-Lewis Index (TLI)                            0.609
#RMSEA                                          0.108       0.102
#Robust RMSEA                                               0.120
#SRMR                                           0.101       0.101

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .475

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.642    0.024  -27.016    0.000   -0.642   -0.642
#    factor3          -0.398    0.029  -13.524    0.000   -0.398   -0.398
#  factor2 ~~                                                            
#    factor3           0.712    0.019   38.414    0.000    0.712    0.712

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.707       0.706
#Tucker-Lewis Index (TLI)                       0.683       0.682
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.687
#RMSEA                                          0.091       0.085
#Robust RMSEA                                               0.090
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .280

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.606    0.033  -18.119    0.000   -0.606   -0.606
#  factor3          -0.465    0.040  -11.650    0.000   -0.465   -0.465
#factor2 ~~                                                            
#  factor3           0.787    0.026   30.720    0.000    0.787    0.787

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.634       0.628
#Tucker-Lewis Index (TLI)                       0.607       0.601
#Robust Comparative Fit Index (CFI)                         0.637
#Robust Tucker-Lewis Index (TLI)                            0.610
#RMSEA                                          0.102       0.095
#Robust RMSEA                                               0.101
#SRMR                                           0.182       0.182

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .259


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ BIS_12+BIS_8+BIS_9+BIS_13+BIS_20+BIS_1+BIS_7+BIS_15+BIS_30+BIS_10+
            BIS_29
 factor2 =~ BIS_19+BIS_17+BIS_14+BIS_18+BIS_2+BIS_4+BIS_5+BIS_6+BIS_3+BIS_27+BIS_23
 factor3 =~ BIS_28+BIS_24+BIS_21+BIS_22+BIS_26+BIS_25+BIS_11+BIS_16 
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15+BIS_16+BIS_17+BIS_18+BIS_19+BIS_20+
            BIS_21+BIS_22+BIS_23+BIS_24+BIS_25+BIS_26+BIS_27+BIS_28+BIS_29+BIS_30
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.803       0.789
#Tucker-Lewis Index (TLI)                       0.772       0.756
#Robust Comparative Fit Index (CFI)                         0.805
#Robust Tucker-Lewis Index (TLI)                            0.773
#RMSEA                                          0.077       0.074
#Robust RMSEA                                               0.077
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .345

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5221267      0.6827586      0.7651771      0.2641590 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.6047738 0.21307057 0.3952262 0.8678704 0.5592475 0.7787250 0.8779084
#factor2 0.2803114 0.09019054 0.7196886 0.7952617 0.1745687 0.5829661 0.8052024
#factor3 0.5357284 0.17461224 0.4642716 0.8175715 0.1620586 2.2254999 1.7231290
#general 0.5221267 0.52212665 0.5221267 0.7651771 0.2641590 0.9045196 0.9366486








################################################################
### Barratt Impulsiveness Scale (BIS-11)
### 
### Finley et al (2019) https://econtent.hogrefe.com/doi/abs/10.1027/1864-9335/a000376?journalCode=zsp
### data available at https://osf.io/kwbmf/?view_only=
###
### 9. 

library(haven)
Finley <- read_sav("Finley.sav")
colnames(Finley)
# 34 items because BIS-10; items 19, 26, 27 and 29 were dropped (https://homepages.se.edu/cvonbergen/files/2013/01/Factor-Structure-of-the-Barratt-Impulsiveness-Scale.pdf)
mydata <- as.data.frame(Finley[,129:162])
mydata <- mydata[,-c(19,26,27,29)]
#mydata[mydata<=0] <- NA
mydata[mydata>=5] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("BIS_1","BIS_2","BIS_3","BIS_4","BIS_5","BIS_6","BIS_7","BIS_8",
                      "BIS_9","BIS_10","BIS_11","BIS_12","BIS_13","BIS_14","BIS_15",
                      "BIS_16","BIS_17","BIS_18","BIS_19","BIS_20","BIS_21","BIS_22",
                      "BIS_23","BIS_24","BIS_25","BIS_26","BIS_27","BIS_28","BIS_29","BIS_30")
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .80, omega T = .83

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 4.67; eigenvalue 2 = 2.01

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 5.46; eigenvalue 2 = 2.38

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.092, RMSR=.10, TLI=.450

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.18, RMSEA=.124, RMSR=.12, TLI=.358

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 7 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al.)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.28, RMSEA=.062, RMSR=.06, TLI=.747
#      MR2   MR1   MR3
#MR2  1.00 -0.34 -0.12
#MR1 -0.34  1.00  0.25
#MR3 -0.12  0.25  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.097, RMSR=.07, TLI=.606
#      MR2   MR1   MR3
#MR2  1.00 -0.34 -0.12
#MR1 -0.34  1.00  0.25
#MR3 -0.12  0.25  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15+BIS_16+BIS_17+BIS_18+BIS_19+BIS_20+
            BIS_21+BIS_22+BIS_23+BIS_24+BIS_25+BIS_26+BIS_27+BIS_28+BIS_29+BIS_30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.781       0.665
#Tucker-Lewis Index (TLI)                       0.765       0.641
#Robust Comparative Fit Index (CFI)                         0.431
#Robust Tucker-Lewis Index (TLI)                            0.389
#RMSEA                                          0.120       0.099
#Robust RMSEA                                               0.123
#SRMR                                           0.119       0.119

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .131

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.486       0.480
#Tucker-Lewis Index (TLI)                       0.448       0.441
#Robust Comparative Fit Index (CFI)                         0.489
#Robust Tucker-Lewis Index (TLI)                            0.451
#RMSEA                                          0.095       0.091
#Robust RMSEA                                               0.094
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .074

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (based on Kapitány-Fövény) 
EGAmodel= '
 factor1 =~ BIS_12+BIS_8+BIS_9+BIS_13+BIS_20+BIS_1+BIS_7+BIS_15+BIS_30+BIS_10+
            BIS_29
 factor2 =~ BIS_19+BIS_17+BIS_14+BIS_18+BIS_2+BIS_4+BIS_5+BIS_6+BIS_3+BIS_27+BIS_23
 factor3 =~ BIS_28+BIS_24+BIS_21+BIS_22+BIS_26+BIS_25+BIS_11+BIS_16'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.866       0.789
#Tucker-Lewis Index (TLI)                       0.855       0.772
#Robust Comparative Fit Index (CFI)                         0.588
#Robust Tucker-Lewis Index (TLI)                            0.554
#RMSEA                                          0.094       0.079
#Robust RMSEA                                               0.105
#SRMR                                           0.104       0.104

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .297

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.548    0.047  -11.728    0.000   -0.548   -0.548
#    factor3          -0.217    0.061   -3.548    0.000   -0.217   -0.217
#  factor2 ~~                                                            
#    factor3           0.616    0.044   14.166    0.000    0.616    0.616

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.649       0.647
#Tucker-Lewis Index (TLI)                       0.620       0.618
#Robust Comparative Fit Index (CFI)                         0.654
#Robust Tucker-Lewis Index (TLI)                            0.626
#RMSEA                                          0.079       0.075
#Robust RMSEA                                               0.077
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .154

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.516    0.063   -8.149    0.000   -0.516   -0.516
#  factor3          -0.082    0.091   -0.897    0.370   -0.082   -0.082
#factor2 ~~                                                            
#  factor3           0.328    0.102    3.203    0.001    0.328    0.328

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.611       0.611
#Tucker-Lewis Index (TLI)                       0.583       0.583
#Robust Comparative Fit Index (CFI)                         0.617
#Robust Tucker-Lewis Index (TLI)                            0.589
#RMSEA                                          0.082       0.078
#Robust RMSEA                                               0.081
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .144


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ BIS_12+BIS_8+BIS_9+BIS_13+BIS_20+BIS_1+BIS_7+BIS_15+BIS_30+BIS_10+
            BIS_29
 factor2 =~ BIS_19+BIS_17+BIS_14+BIS_18+BIS_2+BIS_4+BIS_5+BIS_6+BIS_3+BIS_27+BIS_23
 factor3 =~ BIS_28+BIS_24+BIS_21+BIS_22+BIS_26+BIS_25+BIS_11+BIS_16 
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15+BIS_16+BIS_17+BIS_18+BIS_19+BIS_20+
            BIS_21+BIS_22+BIS_23+BIS_24+BIS_25+BIS_26+BIS_27+BIS_28+BIS_29+BIS_30
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.769       0.765
#Tucker-Lewis Index (TLI)                       0.732       0.727
#Robust Comparative Fit Index (CFI)                         0.773
#Robust Tucker-Lewis Index (TLI)                            0.736
#RMSEA                                          0.066       0.063
#Robust RMSEA                                               0.065
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .275

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4358946      0.6827586      0.6939282      0.2414379 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.7243605 0.2848826 0.2756395 0.8069519 0.65222673 0.7741739 0.8744949
#factor2 0.4212803 0.1612086 0.5787197 0.7268441 0.17127621 0.6667920 0.8450377
#factor3 0.5267361 0.1180143 0.4732639 0.5663183 0.07623748 0.6725345 0.8573351
#general 0.4358946 0.4358946 0.4358946 0.6939282 0.24143794 0.8300007 0.8997135








################################################################
### Barratt Impulsiveness Scale (BIS-15)
### 
### Meule et al (2019) 
### data available at https://data.mendeley.com/datasets/3x3bzfdmbs/1
###
### 10. 

library(haven)
Meule <- read_sav("Meule.sav")
colnames(Meule)
# 15 items because BIS-15
mydata <- as.data.frame(Meule[,5:19])
#mydata[mydata<=0] <- NA
#mydata[mydata>=5] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("BIS_1","BIS_2","BIS_3","BIS_4","BIS_5","BIS_6","BIS_7","BIS_8",
                      "BIS_9","BIS_10","BIS_11","BIS_12","BIS_13","BIS_14","BIS_15")
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .82, omega T = .87

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 3.69; eigenvalue 2 = 1.24

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 4.24; eigenvalue 2 = 1.42

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.152, RMSR=.13, TLI=.467

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.28, RMSEA=.188, RMSR=.15, TLI=.412

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al.)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.082, RMSR=.04, TLI=.844
#     MR2  MR1  MR3
#MR2 1.00 0.32 0.22
#MR1 0.32 1.00 0.27
#MR3 0.22 0.27 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.116, RMSR=.05, TLI=.775
#     MR3  MR1  MR2
#MR3 1.00 0.33 0.22
#MR1 0.33 1.00 0.27
#MR2 0.22 0.27 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.812       0.660
#Tucker-Lewis Index (TLI)                       0.781       0.603
#Robust Comparative Fit Index (CFI)                         0.483
#Robust Tucker-Lewis Index (TLI)                            0.397
#RMSEA                                          0.183       0.182
#Robust RMSEA                                               0.192
#SRMR                                           0.148       0.148

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .379

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.544       0.515
#Tucker-Lewis Index (TLI)                       0.468       0.435
#Robust Comparative Fit Index (CFI)                         0.544
#Robust Tucker-Lewis Index (TLI)                            0.468
#RMSEA                                          0.153       0.148
#Robust RMSEA                                               0.153
#SRMR                                           0.119       0.119

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .232

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (based on theory; https://www.sciencedirect.com/science/article/pii/S0165178119318864?via%3Dihub) 
EGAmodel= '
 factor1 =~ BIS_15+BIS_8+BIS_5+BIS_1+BIS_7
 factor2 =~ BIS_10+BIS_12+BIS_2+BIS_9+BIS_13
 factor3 =~ BIS_14+BIS_6+BIS_4+BIS_3+BIS_11
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.903
#Tucker-Lewis Index (TLI)                       0.947       0.883
#Robust Comparative Fit Index (CFI)                         0.801
#Robust Tucker-Lewis Index (TLI)                            0.759
#RMSEA                                          0.090       0.098
#Robust RMSEA                                               0.121
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .553

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.479    0.043   11.094    0.000    0.479    0.479
#    factor3           0.287    0.053    5.442    0.000    0.287    0.287
#  factor2 ~~                                                            
#    factor3           0.437    0.046    9.581    0.000    0.437    0.437

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.843       0.839
#Tucker-Lewis Index (TLI)                       0.811       0.805
#Robust Comparative Fit Index (CFI)                         0.845
#Robust Tucker-Lewis Index (TLI)                            0.813
#RMSEA                                          0.091       0.087
#Robust RMSEA                                               0.090
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .446

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.424    0.062    6.850    0.000    0.424    0.424
#    factor3           0.257    0.068    3.791    0.000    0.257    0.257
#  factor2 ~~                                                            
#    factor3           0.403    0.065    6.182    0.000    0.403    0.403

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.795       0.791
#Tucker-Lewis Index (TLI)                       0.761       0.756
#Robust Comparative Fit Index (CFI)                         0.798
#Robust Tucker-Lewis Index (TLI)                            0.764
#RMSEA                                          0.103       0.097
#Robust RMSEA                                               0.102
#SRMR                                           0.153       0.153

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .422


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ BIS_15+BIS_8+BIS_5+BIS_1+BIS_7
 factor2 =~ BIS_10+BIS_12+BIS_2+BIS_9+BIS_13
 factor3 =~ BIS_14+BIS_6+BIS_4+BIS_3+BIS_11
 general =~ BIS_1+BIS_2+BIS_3+BIS_4+BIS_5+BIS_6+BIS_7+BIS_8+BIS_9+BIS_10+
            BIS_11+BIS_12+BIS_13+BIS_14+BIS_15
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.921       0.917
#Tucker-Lewis Index (TLI)                       0.889       0.883
#Robust Comparative Fit Index (CFI)                         0.922
#Robust Tucker-Lewis Index (TLI)                            0.891
#RMSEA                                          0.070       0.067
#Robust RMSEA                                               0.069
#SRMR                                           0.053       0.053

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .468

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3970466      0.7142857      0.8716884      0.5839203 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.6374961 0.2345979 0.3625039 0.8229405 0.5124580 0.7653610 0.8688979
#factor2 0.4827714 0.1768180 0.5172286 0.8183175 0.3783986 0.6884003 0.8223110
#factor3 0.7207572 0.1915374 0.2792428 0.7236189 0.5145604 0.7039683 0.8384095
#general 0.3970466 0.3970466 0.3970466 0.8716884 0.5839203 0.8030056 0.8566269










################################################################
### Barratt Impulsiveness Scale (BIS-11)
### Demetrovics (2013) https://core.ac.uk/download/pdf/35136241.pdf
### 
### data obtained from Demetrovics from a gambling study collected in 2011
###
### 

library(haven)
BIS_Demetrovics1 <- read_sav("BIS_Demetrovics1.sav")
colnames(BIS_Demetrovics1)
mydata <- as.data.frame(BIS_Demetrovics1[,3:32])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .81, omega T = .84

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 4.57; eigenvalue 2 = 2.54

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 5.66; eigenvalue 2 = 3.02

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.15, RMSEA=.102, RMSR=.11, TLI=.384

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.137, RMSR=.13, TLI=.315

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al.)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.28, RMSEA=.072, RMSR=.06, TLI=.692
#      MR1   MR3   MR2
#MR1  1.00 -0.22 -0.16
#MR3 -0.22  1.00  0.36
#MR2 -0.16  0.36  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.109, RMSR=.07, TLI=.565
#      MR2   MR1   MR3
#MR2  1.00 -0.27 -0.08
#MR1 -0.27  1.00  0.27
#MR3 -0.08  0.27  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS1+BIS2+BIS3+BIS4+BIS5+BIS6+BIS7+BIS8+BIS9+BIS10+
            BIS11+BIS12+BIS13+BIS14+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20+
            BIS21+BIS22+BIS23+BIS24+BIS25+BIS26+BIS27+BIS28+BIS29+BIS30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.710       0.618
#Tucker-Lewis Index (TLI)                       0.689       0.590
#Robust Comparative Fit Index (CFI)                         0.382
#Robust Tucker-Lewis Index (TLI)                            0.336
#RMSEA                                          0.136       0.100
#Robust RMSEA                                               0.135
#SRMR                                           0.133       0.133

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .205

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.425       0.401
#Tucker-Lewis Index (TLI)                       0.383       0.357
#Robust Comparative Fit Index (CFI)                         0.425
#Robust Tucker-Lewis Index (TLI)                            0.383
#RMSEA                                          0.103       0.099
#Robust RMSEA                                               0.102
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .146

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (based on Kapitány-Fövény et al. and EFA because order of items did not seem to correspond)
EGAmodel= '
factor1 =~ BIS9+BIS12+BIS8+BIS13+BIS20+BIS1+BIS7+BIS15+BIS30+BIS10+
  BIS29
factor2 =~ BIS19+BIS17+BIS14+BIS18+BIS2+BIS4+BIS5+BIS6+BIS3+BIS27+BIS23
factor3 =~ BIS24+BIS21+BIS22+BIS26+BIS25+BIS28+BIS11+BIS16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.861       0.735
#Tucker-Lewis Index (TLI)                       0.850       0.713
#Robust Comparative Fit Index (CFI)                         0.596
#Robust Tucker-Lewis Index (TLI)                            0.563
#RMSEA                                          0.094       0.084
#Robust RMSEA                                               0.110
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .333

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.357    0.050   -7.161    0.000   -0.357   -0.357
#    factor3          -0.278    0.054   -5.111    0.000   -0.278   -0.278
#  factor2 ~~                                                            
#    factor3           0.691    0.033   21.254    0.000    0.691    0.691

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.657       0.653
#Tucker-Lewis Index (TLI)                       0.629       0.625
#Robust Comparative Fit Index (CFI)                         0.662
#Robust Tucker-Lewis Index (TLI)                            0.634
#RMSEA                                          0.080       0.076
#Robust RMSEA                                               0.079
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .237

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.321    0.086   -3.746    0.000   -0.321   -0.321
#    factor3          -0.275    0.083   -3.323    0.001   -0.275   -0.275
#  factor2 ~~                                                            
#    factor3           0.720    0.046   15.745    0.000    0.720    0.720


library(semPlot)
semPaths(CFA_model2,"std", layout="tree", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.606       0.594
#Tucker-Lewis Index (TLI)                       0.577       0.564
#Robust Comparative Fit Index (CFI)                         0.609
#Robust Tucker-Lewis Index (TLI)                            0.580
#RMSEA                                          0.085       0.081
#Robust RMSEA                                               0.085
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .239


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ BIS9+BIS12+BIS8+BIS13+BIS20+BIS1+BIS7+BIS15+BIS30+BIS10+BIS29
 factor2 =~ BIS19+BIS17+BIS14+BIS18+BIS2+BIS4+BIS5+BIS6+BIS3+BIS27+BIS23
 factor3 =~ BIS24+BIS21+BIS22+BIS26+BIS25+BIS28+BIS11+BIS16
 general =~ BIS1+BIS2+BIS3+BIS4+BIS5+BIS6+BIS7+BIS8+BIS9+BIS10+
            BIS11+BIS12+BIS13+BIS14+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20+
            BIS21+BIS22+BIS23+BIS24+BIS25+BIS26+BIS27+BIS28+BIS29+BIS30
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.733       0.717
#Tucker-Lewis Index (TLI)                       0.690       0.672
#Robust Comparative Fit Index (CFI)                         0.734
#Robust Tucker-Lewis Index (TLI)                            0.692
#RMSEA                                          0.073       0.071
#Robust RMSEA                                               0.072
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .280

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4021181      0.6827586      0.7547049      0.2217674 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8046458 0.3089872 0.1953542 0.8200450 0.6955131 0.7948698 0.8845712
#factor2 0.4605790 0.1551405 0.5394210 0.7567870 0.3410170 0.6459446 0.7898722
#factor3 0.4791342 0.1337541 0.5208658 0.7184913 0.2713925 0.6823190 0.8235645
#general 0.4021181 0.4021181 0.4021181 0.7547049 0.2217674 0.8216903 0.8800858







################################################################
### Barratt Impulsiveness Scale (BIS-11)
### Demetrovics (2013) https://core.ac.uk/download/pdf/35136241.pdf
### 
### data obtained from Demetrovics from a gambling study collected in 2011 (https://core.ac.uk/download/pdf/35136241.pdf)
###
### 

library(haven)
BIS_Demetrovics1 <- read_sav("BIS_Demetrovics1.sav")
colnames(BIS_Demetrovics1)
mydata <- as.data.frame(BIS_Demetrovics1[,3:32])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .81, omega T = .84

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 4.57; eigenvalue 2 = 2.54

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 5.66; eigenvalue 2 = 3.02

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.15, RMSEA=.102, RMSR=.11, TLI=.384

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.137, RMSR=.13, TLI=.315

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al.)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.28, RMSEA=.072, RMSR=.06, TLI=.692
#      MR1   MR3   MR2
#MR1  1.00 -0.22 -0.16
#MR3 -0.22  1.00  0.36
#MR2 -0.16  0.36  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.109, RMSR=.07, TLI=.565
#      MR2   MR1   MR3
#MR2  1.00 -0.27 -0.08
#MR1 -0.27  1.00  0.27
#MR3 -0.08  0.27  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS1+BIS2+BIS3+BIS4+BIS5+BIS6+BIS7+BIS8+BIS9+BIS10+
            BIS11+BIS12+BIS13+BIS14+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20+
            BIS21+BIS22+BIS23+BIS24+BIS25+BIS26+BIS27+BIS28+BIS29+BIS30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.710       0.618
#Tucker-Lewis Index (TLI)                       0.689       0.590
#Robust Comparative Fit Index (CFI)                         0.382
#Robust Tucker-Lewis Index (TLI)                            0.336
#RMSEA                                          0.136       0.100
#Robust RMSEA                                               0.135
#SRMR                                           0.133       0.133

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .205

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.425       0.401
#Tucker-Lewis Index (TLI)                       0.383       0.357
#Robust Comparative Fit Index (CFI)                         0.425
#Robust Tucker-Lewis Index (TLI)                            0.383
#RMSEA                                          0.103       0.099
#Robust RMSEA                                               0.102
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .146

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (based on Kapitány-Fövény et al. and EFA because order of items did not seem to correspond)
EGAmodel= '
factor1 =~ BIS9+BIS12+BIS8+BIS13+BIS20+BIS1+BIS7+BIS15+BIS30+BIS10+
  BIS29
factor2 =~ BIS19+BIS17+BIS14+BIS18+BIS2+BIS4+BIS5+BIS6+BIS3+BIS27+BIS23
factor3 =~ BIS24+BIS21+BIS22+BIS26+BIS25+BIS28+BIS11+BIS16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.861       0.735
#Tucker-Lewis Index (TLI)                       0.850       0.713
#Robust Comparative Fit Index (CFI)                         0.596
#Robust Tucker-Lewis Index (TLI)                            0.563
#RMSEA                                          0.094       0.084
#Robust RMSEA                                               0.110
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .333

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.357    0.050   -7.161    0.000   -0.357   -0.357
#    factor3          -0.278    0.054   -5.111    0.000   -0.278   -0.278
#  factor2 ~~                                                            
#    factor3           0.691    0.033   21.254    0.000    0.691    0.691

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.657       0.653
#Tucker-Lewis Index (TLI)                       0.629       0.625
#Robust Comparative Fit Index (CFI)                         0.662
#Robust Tucker-Lewis Index (TLI)                            0.634
#RMSEA                                          0.080       0.076
#Robust RMSEA                                               0.079
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .237

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.321    0.086   -3.746    0.000   -0.321   -0.321
#    factor3          -0.275    0.083   -3.323    0.001   -0.275   -0.275
#  factor2 ~~                                                            
#    factor3           0.720    0.046   15.745    0.000    0.720    0.720


library(semPlot)
semPaths(CFA_model2,"std", layout="tree", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.606       0.594
#Tucker-Lewis Index (TLI)                       0.577       0.564
#Robust Comparative Fit Index (CFI)                         0.609
#Robust Tucker-Lewis Index (TLI)                            0.580
#RMSEA                                          0.085       0.081
#Robust RMSEA                                               0.085
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .239


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ BIS9+BIS12+BIS8+BIS13+BIS20+BIS1+BIS7+BIS15+BIS30+BIS10+BIS29
 factor2 =~ BIS19+BIS17+BIS14+BIS18+BIS2+BIS4+BIS5+BIS6+BIS3+BIS27+BIS23
 factor3 =~ BIS24+BIS21+BIS22+BIS26+BIS25+BIS28+BIS11+BIS16
 general =~ BIS1+BIS2+BIS3+BIS4+BIS5+BIS6+BIS7+BIS8+BIS9+BIS10+
            BIS11+BIS12+BIS13+BIS14+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20+
            BIS21+BIS22+BIS23+BIS24+BIS25+BIS26+BIS27+BIS28+BIS29+BIS30
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.733       0.717
#Tucker-Lewis Index (TLI)                       0.690       0.672
#Robust Comparative Fit Index (CFI)                         0.734
#Robust Tucker-Lewis Index (TLI)                            0.692
#RMSEA                                          0.073       0.071
#Robust RMSEA                                               0.072
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .280

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4021181      0.6827586      0.7547049      0.2217674 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.8046458 0.3089872 0.1953542 0.8200450 0.6955131 0.7948698 0.8845712
#factor2 0.4605790 0.1551405 0.5394210 0.7567870 0.3410170 0.6459446 0.7898722
#factor3 0.4791342 0.1337541 0.5208658 0.7184913 0.2713925 0.6823190 0.8235645
#general 0.4021181 0.4021181 0.4021181 0.7547049 0.2217674 0.8216903 0.8800858







################################################################
### Barratt Impulsiveness Scale (BIS-R_ 21)
### 21 items
### see Kapitány-Fövény et al. (2020), https://akjournals.com/view/journals/2006/9/2/article-p225.xml
### 
### data obtained from Demetrovics from a 2014 international gaming sample
### https://psycnet.apa.org/record/2018-64911-001
### 

library(haven)
BIS_Demetrovics2 <- read_sav("BIS_Demetrovics2.sav")
colnames(BIS_Demetrovics2)
mydata <- as.data.frame(BIS_Demetrovics2[,10:30])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 4

library(psych)
omega(mydata) # alpha = .81, omega T = .85

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 3.8; eigenvalue 2 = 2.21

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 4.61; eigenvalue 2 = 2.58

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.116, RMSR=.12, TLI=.441

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.22, RMSEA=.142, RMSR=.14, TLI=.41

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities


# Give solution with 3 factors (based on Kapitány-Fövény et al.)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.065, RMSR=.04, TLI=.826
#      MR2   MR1   MR3
#MR2  1.00 -0.24 -0.13
#MR1 -0.24  1.00  0.45
#MR3 -0.13  0.45  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.087, RMSR=.05, TLI=.778
#      MR2   MR1   MR3
#MR2  1.00 -0.27 -0.08
#MR1 -0.27  1.00  0.27
#MR3 -0.08  0.27  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~ BIS1+BIS2+BIS3+BIS4+BIS5+BIS6+BIS7+BIS8+BIS9+BIS10+
            BIS11+BIS12+BIS13+BIS14+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20+BIS21
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.743       0.609
#Tucker-Lewis Index (TLI)                       0.714       0.566
#Robust Comparative Fit Index (CFI)                         0.461
#Robust Tucker-Lewis Index (TLI)                            0.401
#RMSEA                                          0.149       0.135
#Robust RMSEA                                               0.143
#SRMR                                           0.140       0.140

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .232

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.497       0.212
#Tucker-Lewis Index (TLI)                       0.441       0.124
#Robust Comparative Fit Index (CFI)                         0.496
#Robust Tucker-Lewis Index (TLI)                            0.440
#RMSEA                                          0.116       0.134
#Robust RMSEA                                               0.116
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .172

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (based on Kapitány-Fövény et al. and EFA because order of items did not seem to correspond)
EGAmodel= '
 cognitive  =~ BIS1+BIS3+BIS4+BIS5+BIS6+BIS8+BIS9+BIS14+BIS21
 behavioral =~ BIS2+BIS10+BIS11+BIS12+BIS13
 restless   =~ BIS7+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.953       0.905
#Tucker-Lewis Index (TLI)                       0.947       0.893
#Robust Comparative Fit Index (CFI)                         0.810
#Robust Tucker-Lewis Index (TLI)                            0.785
#RMSEA                                          0.064       0.067
#Robust RMSEA                                               0.086
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .387

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  cognitive ~~                                                          
#    behavioral       -0.323    0.011  -30.098    0.000   -0.323   -0.323
#    restless         -0.184    0.011  -16.096    0.000   -0.184   -0.184
#  behavioral ~~                                                         
#    restless          0.702    0.007   95.423    0.000    0.702    0.702


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.850       0.846
#Tucker-Lewis Index (TLI)                       0.831       0.827
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.831
#RMSEA                                          0.064       0.060
#Robust RMSEA                                               0.064
#SRMR                                           0.050       0.050

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .302

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  cognitive ~~                                                          
#    behavioral       -0.301    0.013  -23.409    0.000   -0.301   -0.301
#    restless         -0.169    0.013  -12.534    0.000   -0.169   -0.169
#  behavioral ~~                                                         
#    restless          0.717    0.010   75.183    0.000    0.717    0.717



library(semPlot)
semPaths(CFA_model2,"std", layout="tree", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.774       0.764
#Tucker-Lewis Index (TLI)                       0.749       0.738
#Robust Comparative Fit Index (CFI)                         0.774
#Robust Tucker-Lewis Index (TLI)                            0.749
#RMSEA                                          0.078       0.073
#Robust RMSEA                                               0.078
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .301


# Bifactor model with 3 factors
BIFmodel= '
 cognitive  =~ BIS1+BIS3+BIS4+BIS5+BIS6+BIS8+BIS9+BIS14+BIS21
 behavioral =~ BIS2+BIS10+BIS11+BIS12+BIS13
 restless   =~ BIS7+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20
 general =~ BIS1+BIS2+BIS3+BIS4+BIS5+BIS6+BIS7+BIS8+BIS9+BIS10+
            BIS11+BIS12+BIS13+BIS14+BIS15+BIS16+BIS17+BIS18+BIS19+BIS20+BIS21
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.905       0.903
#Tucker-Lewis Index (TLI)                       0.881       0.878
#Robust Comparative Fit Index (CFI)                         0.905
#Robust Tucker-Lewis Index (TLI)                            0.882
#RMSEA                                          0.053       0.050
#Robust RMSEA                                               0.053
#SRMR                                           0.045       0.045

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .338

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3947068      0.6809524      0.7934772      0.3282246 
#
#$FactorLevelIndices
#              ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#cognitive  0.9424159 0.38285027 0.05758405 0.8198246 0.780322474 0.8173792 0.9028096
#behavioral 0.4901805 0.13485695 0.50981951 0.7712781 0.357349228 0.5803643 0.7617445
#restless   0.2748746 0.08758594 0.72512542 0.7237897 0.002269211 0.4631095 0.7518075
#general    0.3947068 0.39470684 0.39470684 0.7934772 0.328224599 0.7984122 0.8857583








