################################################################
### COVIDSTRESS database
### 
### Blackburn et al. (2022) https://www.nature.com/articles/s41597-022-01383-6
### data available at https://osf.io/36tsd/
###
### Moral Values Scale
###
### Full dataset to determine 2-factor solution
###

Final_COVIDiSTRESS_Vol2_cleaned <- read.csv("Final_COVIDiSTRESS_Vol2_cleaned.csv")
colnames(Final_COVIDiSTRESS_Vol2_cleaned)
table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)

mydata <- Final_COVIDiSTRESS_Vol2_cleaned[,152:162]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
fit1 <- fa(mydata,2)
fit1
diagram(fit1)

EGAmodel= '
 factor1 =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral
 factor2 =~ moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'


################################################################
### COVIDSTRESS database
### 
### UserLanguage = CS (Czech )
### 
###

mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="CS")
mydata <- mydata[,152:162]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .70, omega T = .78

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 2.39; eigenvalue 2 = .98

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 2.7; eigenvalue 2 = 1.07

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.128, RMSR=.11, TLI=.55

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.15, RMSR=.13, TLI=.522

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.078, RMSR=.05, TLI=.831
#      MR1   MR2
#MR1  1.00  0.15
#MR2  0.15  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.103, RMSR=.05, TLI=.772
#      MR1   MR2
#MR1  1.00  0.14
#MR2  0.14  1.00

# Single factor model lavaan
UNImodel= '
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# Warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.645       0.627
#Comparative Fit Index (CFI)                    0.813       0.665
#Tucker-Lewis Index (TLI)                       0.766       0.581
#Robust Comparative Fit Index (CFI)                         0.621
#Robust Tucker-Lewis Index (TLI)                            0.527
#RMSEA                                          0.156       0.163
#Robust RMSEA                                               0.152
#SRMR                                           0.116       0.116

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .271

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.645       0.627
#Tucker-Lewis Index (TLI)                       0.557       0.534
#Robust Comparative Fit Index (CFI)                         0.648
#Robust Tucker-Lewis Index (TLI)                            0.560
#RMSEA                                          0.129       0.120
#Robust RMSEA                                               0.127
#SRMR                                           0.106       0.106

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .202

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.860       0.735
#Tucker-Lewis Index (TLI)                       0.822       0.661
#Robust Comparative Fit Index (CFI)                         0.717
#Robust Tucker-Lewis Index (TLI)                            0.638
#RMSEA                                          0.136       0.147
#Robust RMSEA                                               0.133
#SRMR                                           0.103       0.103

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .291

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          0.527    0.062    8.517    0.000    0.527    0.527

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.743       0.719
#Tucker-Lewis Index (TLI)                       0.671       0.640
#Robust Comparative Fit Index (CFI)                         0.744
#Robust Tucker-Lewis Index (TLI)                            0.673
#RMSEA                                          0.111       0.106
#Robust RMSEA                                               0.110
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .261

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          0.464    0.148    3.130    0.002    0.464    0.464

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.706       0.684
#Tucker-Lewis Index (TLI)                       0.633       0.605
#Robust Comparative Fit Index (CFI)                         0.708
#Robust Tucker-Lewis Index (TLI)                            0.635
#RMSEA                                          0.117       0.111
#Robust RMSEA                                               0.116
#SRMR                                           0.130       0.130

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .265


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral
 factor2 =~ moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.922       0.919
#Tucker-Lewis Index (TLI)                       0.870       0.864
#Robust Comparative Fit Index (CFI)                         0.925
#Robust Tucker-Lewis Index (TLI)                            0.875
#RMSEA                                          0.070       0.065
#Robust RMSEA                                               0.068
#SRMR                                           0.056       0.056

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .290

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    mrl.vls_2_mdnt    0.602    0.075    8.032    0.000    0.602    0.620
#    mrl.vls_8_mdnt    0.788    0.103    7.652    0.000    0.788    0.740
#    mrl.vls_6_mdnt    0.300    0.099    3.035    0.002    0.300    0.328
#    mrl.vls_1_mdnt    0.443    0.182    2.426    0.015    0.443    0.313
#    mrl.vls_7_mdnt    0.615    0.149    4.136    0.000    0.615    0.383
#    mrl.vls_5_mdnt    0.031    0.154    0.204    0.838    0.031    0.018
#  factor2 =~                                                            
#    mrl.vls_11_mdn    0.370    0.603    0.613    0.540    0.370    0.245
#    mrl.vls_10_mdn    3.325    5.271    0.631    0.528    3.325    1.793
#    mrl.vls_9_mdnt    0.141    0.285    0.494    0.621    0.141    0.092
#    mrl.vls_4_mdnt    0.204    0.360    0.567    0.570    0.204    0.129
#    mrl.vls_3_mdnt    0.150    0.268    0.559    0.576    0.150    0.110
#  general =~                                                            
#    mrl.vls_2_mdnt    0.355    0.117    3.045    0.002    0.355    0.366
#    mrl.vls_8_mdnt    0.335    0.125    2.685    0.007    0.335    0.314
#    mrl.vls_6_mdnt    0.303    0.120    2.525    0.012    0.303    0.331
#    mrl.vls_1_mdnt    0.659    0.174    3.789    0.000    0.659    0.466
#    mrl.vls_7_mdnt    0.549    0.173    3.167    0.002    0.549    0.342
#    mrl.vls_5_mdnt    1.001    0.135    7.396    0.000    1.001    0.565
#    mrl.vls_11_mdn    0.524    0.230    2.283    0.022    0.524    0.347
#    mrl.vls_10_mdn   -0.078    0.259   -0.302    0.762   -0.078   -0.042
#    mrl.vls_9_mdnt    0.545    0.122    4.462    0.000    0.545    0.359
#    mrl.vls_4_mdnt    1.046    0.148    7.046    0.000    1.046    0.660
#    mrl.vls_3_mdnt    0.646    0.105    6.141    0.000    0.646    0.476

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.2934469      0.5454545      0.8653253      0.5243485 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5633366 0.1975469 0.4366634 0.7549153 0.3802036 0.6911849 0.8298090
#factor2 0.7838980 0.5090062 0.2161020 0.9193437 0.5830805 3.8552324 1.9244119
#general 0.2934469 0.2934469 0.2934469 0.8653253 0.5243485 0.7236094 0.8583543





################################################################
### COVIDSTRESS database
### 
### UserLanguage = ES-CO (Spanish)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="ES-CO")
mydata <- mydata[,152:162]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .69, omega T = .76

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.05; eigenvalue 2 = 1.05

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 2.91; eigenvalue 2 = 1.31

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.122, RMSR=.11, TLI=.508

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.168, RMSR=.14, TLI=.507

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
# %variance explained=.30, RMSEA=.055, RMSR=.04, TLI=.898
#      MR1   MR2
#MR1  1.00  0.28
#MR2  0.28  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.098, RMSR=.05, TLI=.831
#      MR1   MR2
#MR1  1.00  0.28
#MR2  0.28  1.00

# Single factor model lavaan
UNImodel= '
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.807       0.692
#Tucker-Lewis Index (TLI)                       0.759       0.615
#Robust Comparative Fit Index (CFI)                         0.582
#Robust Tucker-Lewis Index (TLI)                            0.478
#RMSEA                                          0.161       0.166
#Robust RMSEA                                               0.173
#SRMR                                           0.142       0.142

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .311

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.612       0.595
#Tucker-Lewis Index (TLI)                       0.515       0.494
#Robust Comparative Fit Index (CFI)                         0.617
#Robust Tucker-Lewis Index (TLI)                            0.521
#RMSEA                                          0.122       0.108
#Robust RMSEA                                               0.120
#SRMR                                           0.105       0.105

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .174

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.954       0.916
#Tucker-Lewis Index (TLI)                       0.941       0.893
#Robust Comparative Fit Index (CFI)                         0.851
#Robust Tucker-Lewis Index (TLI)                            0.809
#RMSEA                                          0.079       0.088
#Robust RMSEA                                               0.105
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .397

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.342    0.055    6.232    0.000    0.342    0.342

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.882       0.855
#Tucker-Lewis Index (TLI)                       0.849       0.815
#Robust Comparative Fit Index (CFI)                         0.882
#Robust Tucker-Lewis Index (TLI)                            0.850
#RMSEA                                          0.068       0.066
#Robust RMSEA                                               0.067
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .270

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.348    0.078    4.448    0.000    0.348    0.348

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.850       0.813
#Tucker-Lewis Index (TLI)                       0.812       0.766
#Robust Comparative Fit Index (CFI)                         0.849
#Robust Tucker-Lewis Index (TLI)                            0.811
#RMSEA                                          0.076       0.074
#Robust RMSEA                                               0.075
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .275


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral
 factor2 =~ moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.968       0.958
#Tucker-Lewis Index (TLI)                       0.946       0.931
#Robust Comparative Fit Index (CFI)                         0.968
#Robust Tucker-Lewis Index (TLI)                            0.946
#RMSEA                                          0.041       0.040
#Robust RMSEA                                               0.040
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .293

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    mrl.vls_2_mdnt    0.448    0.096    4.674    0.000    0.448    0.497
#    mrl.vls_8_mdnt    0.407    0.084    4.870    0.000    0.407    0.469
#    mrl.vls_6_mdnt    0.339    0.091    3.713    0.000    0.339    0.491
#    mrl.vls_1_mdnt    0.399    0.079    5.078    0.000    0.399    0.325
#    mrl.vls_7_mdnt    0.578    0.099    5.838    0.000    0.578    0.504
#    mrl.vls_5_mdnt    0.563    0.110    5.121    0.000    0.563    0.469
#  factor2 =~                                                            
#    mrl.vls_11_mdn    1.333    0.171    7.813    0.000    1.333    0.617
#    mrl.vls_10_mdn    1.310    0.272    4.809    0.000    1.310    0.598
#    mrl.vls_9_mdnt    0.572    0.167    3.425    0.001    0.572    0.298
#    mrl.vls_4_mdnt    0.187    0.358    0.523    0.601    0.187    0.109
#    mrl.vls_3_mdnt    0.384    0.149    2.581    0.010    0.384    0.226
#  general =~                                                            
#    mrl.vls_2_mdnt    0.035    0.044    0.796    0.426    0.035    0.039
#    mrl.vls_8_mdnt    0.234    0.046    5.113    0.000    0.234    0.270
#    mrl.vls_6_mdnt    0.234    0.059    3.983    0.000    0.234    0.339
#    mrl.vls_1_mdnt    0.169    0.071    2.385    0.017    0.169    0.138
#    mrl.vls_7_mdnt    0.142    0.068    2.087    0.037    0.142    0.124
#    mrl.vls_5_mdnt    0.396    0.093    4.282    0.000    0.396    0.330
#    mrl.vls_11_mdn    0.894    0.216    4.143    0.000    0.894    0.414
#    mrl.vls_10_mdn    0.527    0.294    1.794    0.073    0.527    0.241
#    mrl.vls_9_mdnt    0.590    0.126    4.698    0.000    0.590    0.307
#    mrl.vls_4_mdnt    1.548    0.275    5.620    0.000    1.548    0.900
#    mrl.vls_3_mdnt    0.703    0.091    7.727    0.000    0.703    0.414

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.4289247      0.5454545      0.7647765      0.4043177 
#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.7949070 0.3378311 0.2050930 0.6759391 0.5622456 0.6255077 0.7930612
#factor2 0.4056382 0.2332442 0.5943618 0.7537121 0.2995644 0.5716041 0.7536639
#general 0.4289247 0.4289247 0.4289247 0.7647765 0.4043177 0.8384236 0.9070582





################################################################
### COVIDSTRESS database
### 
### UserLanguage = ET (Estonian)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="ET")
mydata <- mydata[,152:162]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .67, omega T = .74

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 1.94; eigenvalue 2 = .86

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 2.35; eigenvalue 2 = .91

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.116, RMSR=.11, TLI=.498

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.21, RMSEA=.131, RMSR=.12, TLI=.527

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 2 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.28, RMSEA=.081, RMSR=.06, TLI=.753
#      MR1   MR2
#MR1  1.00  0.11
#MR2  0.11  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.01, RMSR=.06, TLI=.723
#      MR1   MR2
#MR1  1.00  0.11
#MR2  0.11  1.00

# Single factor model lavaan
UNImodel= '
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.789       0.633
#Tucker-Lewis Index (TLI)                       0.737       0.542
#Robust Comparative Fit Index (CFI)                         0.630
#Robust Tucker-Lewis Index (TLI)                            0.537
#RMSEA                                          0.131       0.143
#Robust RMSEA                                               0.133
#SRMR                                           0.107       0.107

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .267

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.598       0.567
#Tucker-Lewis Index (TLI)                       0.498       0.459
#Robust Comparative Fit Index (CFI)                         0.596
#Robust Tucker-Lewis Index (TLI)                            0.495
#RMSEA                                          0.118       0.115
#Robust RMSEA                                               0.117
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .233

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.896       0.820
#Tucker-Lewis Index (TLI)                       0.867       0.770
#Robust Comparative Fit Index (CFI)                         0.774
#Robust Tucker-Lewis Index (TLI)                            0.711
#RMSEA                                          0.094       0.101
#Robust RMSEA                                               0.105
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .339

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.352    0.078    4.539    0.000    0.352    0.352

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.770       0.752
#Tucker-Lewis Index (TLI)                       0.705       0.683
#Robust Comparative Fit Index (CFI)                         0.770
#Robust Tucker-Lewis Index (TLI)                            0.705
#RMSEA                                          0.090       0.088
#Robust RMSEA                                               0.089
#SRMR                                           0.079       0.079

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .264

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.324    0.124    2.608    0.009    0.324    0.324

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.744       0.724
#Tucker-Lewis Index (TLI)                       0.680       0.655
#Robust Comparative Fit Index (CFI)                         0.744
#Robust Tucker-Lewis Index (TLI)                            0.679
#RMSEA                                          0.094       0.092
#Robust RMSEA                                               0.093
#SRMR                                           0.095       0.095

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .281


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral
 factor2 =~ moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.881       0.876
#Tucker-Lewis Index (TLI)                       0.802       0.793
#Robust Comparative Fit Index (CFI)                         0.883
#Robust Tucker-Lewis Index (TLI)                            0.805
#RMSEA                                          0.074       0.071
#Robust RMSEA                                               0.073
#SRMR                                           0.070       0.070

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .288

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#   0.003006364    0.545454545   36.603725692    0.468510565 
#
#$FactorLevelIndices
#             ECV_SS      ECV_SG      ECV_GS       Omega      OmegaH            H         FD
#factor1 0.959884754 0.001840036 0.040115246   0.7104111   0.7016864    0.7072453  0.8471314
#factor2 0.997064908 0.995153600 0.002935092 103.3704162 102.2903278 6662.4815899 36.1096991
#general 0.003006364 0.003006364 0.003006364  36.6037257   0.4685106   -2.4143638  1.4729212






################################################################
### COVIDSTRESS database
### 
### UserLanguage = PT (Portugal)
### 
###

table(Final_COVIDiSTRESS_Vol2_cleaned$UserLanguage)
mydata <- subset(Final_COVIDiSTRESS_Vol2_cleaned, UserLanguage =="PT")
mydata <- mydata[,152:162]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 7

library(psych)
omega(mydata) # alpha = .68, omega T = .75

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 1.98; eigenvalue 2 = .93

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 2.44; eigenvalue 2 = 1.23

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.117, RMSR=.11, TLI=.517

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.22, RMSEA=.158, RMSR=.14, TLI=.438

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 2 factors (general EFA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.29, RMSEA=.072, RMSR=.05, TLI=.814
#      MR1   MR2
#MR1  1.00  0.19
#MR2  0.19  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.01, RMSR=.06, TLI=.772
#      MR1   MR2
#MR1  1.00  0.20
#MR2  0.20  1.00

# Single factor model lavaan
UNImodel= '
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.759       0.595
#Tucker-Lewis Index (TLI)                       0.699       0.494
#Robust Comparative Fit Index (CFI)                         0.560
#Robust Tucker-Lewis Index (TLI)                            0.449
#RMSEA                                          0.153       0.164
#Robust RMSEA                                               0.157
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .236

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.620       0.575
#Tucker-Lewis Index (TLI)                       0.525       0.468
#Robust Comparative Fit Index (CFI)                         0.618
#Robust Tucker-Lewis Index (TLI)                            0.523
#RMSEA                                          0.117       0.112
#Robust RMSEA                                               0.116
#SRMR                                           0.102       0.102

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .142

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (general EFA based) 
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.871       0.788
#Tucker-Lewis Index (TLI)                       0.836       0.729
#Robust Comparative Fit Index (CFI)                         0.739
#Robust Tucker-Lewis Index (TLI)                            0.666
#RMSEA                                          0.113       0.120
#Robust RMSEA                                               0.123
#SRMR                                           0.099       0.099

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .300

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.432    0.052    8.383    0.000    0.432    0.432

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.766       0.735
#Tucker-Lewis Index (TLI)                       0.701       0.661
#Robust Comparative Fit Index (CFI)                         0.766
#Robust Tucker-Lewis Index (TLI)                            0.701
#RMSEA                                          0.093       0.090
#Robust RMSEA                                               0.092
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .228

#Covariances:
#                 Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.412    0.114    3.625    0.000    0.412    0.412

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.726       0.701
#Tucker-Lewis Index (TLI)                       0.657       0.626
#Robust Comparative Fit Index (CFI)                         0.727
#Robust Tucker-Lewis Index (TLI)                            0.659
#RMSEA                                          0.099       0.094
#Robust RMSEA                                               0.098
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .240


# Bifactor model with 2 factors
BIFmodel= '
 factor1 =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral
 factor2 =~ moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
 general =~ moral.values_2_midneutral+moral.values_8_midneutral+moral.values_6_midneutral+
            moral.values_1_midneutral+moral.values_7_midneutral+moral.values_5_midneutral+
            moral.values_11_midneutral+moral.values_10_midneutral+moral.values_9_midneutral+
            moral.values_4_midneutral+moral.values_3_midneutral
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.971       0.970
#Tucker-Lewis Index (TLI)                       0.952       0.951
#Robust Comparative Fit Index (CFI)                         0.974
#Robust Tucker-Lewis Index (TLI)                            0.956
#RMSEA                                          0.037       0.034
#Robust RMSEA                                               0.035
#SRMR                                           0.039       0.039

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .304

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5055094      0.5454545      0.7533811      0.4875452 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.7724909 0.3298757 0.2275091 0.6518953 0.54732504 0.6347701 0.7970635
#factor2 0.2873003 0.1646148 0.7126997 0.7336005 0.08335014 0.4566452 0.7174505
#general 0.5055094 0.5055094 0.5055094 0.7533811 0.48754520 0.7780190 0.8879172






