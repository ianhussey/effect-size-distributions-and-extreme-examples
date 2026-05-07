################################################################
### Interoceptive accuracy scale
### 
### 
### 

################################################################
###
### First determine overall FA structure based on all data
###

library(haven)
IAS_Murphy1 <- read_sav("IAS_Murphy1.sav")
colnames(IAS_Murphy1)
mydata <- IAS_Murphy1[,3:23]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
Names <- colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5
mydata1 <- mydata

library(haven)
IAS_Brand1 <- read_sav("IAS_Brand1.sav")
colnames(IAS_Brand1)
mydata <- IAS_Brand1[,13:33]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 5
mydata2 <- mydata

library(haven)
IAS_Brand2 <- read_sav("IAS_Brand2.sav")
colnames(IAS_Brand2)
mydata <- IAS_Brand2[,11:31]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 5
mydata3 <- mydata

library(haven)
IAS_Brand3 <- read_sav("IAS_Brand3.sav")
colnames(IAS_Brand3)
mydata <- IAS_Brand3[,4:24]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 5
mydata4 <- mydata

library(haven)
IAS_Campos <- read_sav("IAS_Campos.sav")
colnames(IAS_Campos)
mydata <- IAS_Campos[,76:96]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 5
mydata5 <- mydata

library(haven)
IAS_Todd <- read_sav("IAS_Todd.sav")
colnames(IAS_Todd)
mydata <- IAS_Todd[,80:100]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 5
mydata6 <- mydata

load("IAS_Gaggero.Rda")
colnames(data)
mydata <- data[,82:102]
mydata[mydata<=0] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- Names
min(mydata)
max(mydata) # response alternatives = 5
mydata7 <- mydata

mydata <- rbind(mydata1,mydata2,mydata3,mydata4,mydata5,mydata6,mydata7)
haven::zap_labels(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 4 components
# Eigenvalue 1 = 5.68; eigenvalue 2 = .70

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 8 factors and 4 components
# Eigenvalue 1 = 7.11; eigenvalue 2 = .78

# Give solution with 3 factors
fit3 <- fa(mydata,3)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.072, RMSR=.05, TLI=.814
# No clear evidence for second factor; so, used the solution from the initial Murphy et al. (2020) study


################################################################
### 
### Murphy et al (2020) https://journals.sagepub.com/doi/full/10.1177/1747021819879826
### data Study 1 https://osf.io/aswv3

mydata <- mydata1
haven::zap_labels(mydata)

min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .88, omega T = .9

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 2 components
# Eigenvalue 1 = 5.57; eigenvalue 2 = .86

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 factors and 2 components
# Eigenvalue 1 = 7.02; eigenvalue 2 = .98

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.27, RMSEA=.082, RMSR=.07, TLI=.744

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.11, RMSR=.08, TLI=.695

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 3 factors (EFA-based on the joint data)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.061, RMSR=.04, TLI=.857
#     MR1  MR3  MR2
#MR1 1.00 0.49 0.35
#MR3 0.49 1.00 0.42
#MR2 0.35 0.42 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.063, RMSR=.04, TLI=.895
#     MR1  MR3  MR2
#MR1 1.00 0.67 0.53
#MR3 0.67 1.00 0.47
#MR2 0.53 0.47 1.00


# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.955       0.868
#Tucker-Lewis Index (TLI)                       0.950       0.854
#Robust Comparative Fit Index (CFI)                         0.736
#Robust Tucker-Lewis Index (TLI)                            0.706
#RMSEA                                          0.089       0.096
#Robust RMSEA                                               0.110
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .377

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.769       0.788
#Tucker-Lewis Index (TLI)                       0.743       0.764
#Robust Comparative Fit Index (CFI)                         0.789
#Robust Tucker-Lewis Index (TLI)                            0.766
#RMSEA                                          0.083       0.067
#Robust RMSEA                                               0.078
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .289

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (general EFA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.975       0.921
#Tucker-Lewis Index (TLI)                       0.972       0.910
#Robust Comparative Fit Index (CFI)                         0.814
#Robust Tucker-Lewis Index (TLI)                            0.790
#RMSEA                                          0.066       0.075
#Robust RMSEA                                               0.093
#SRMR                                           0.063       0.063

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .438

#Covariances:
#                     Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#       factor2           0.783    0.024   31.959    0.000    0.783    0.783
#       factor3           0.717    0.031   23.261    0.000    0.717    0.717
#  factor2 ~~                                                            
#       factor3           0.679    0.030   22.679    0.000    0.679    0.679


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.843       0.863
#Tucker-Lewis Index (TLI)                       0.822       0.845
#Robust Comparative Fit Index (CFI)                         0.865
#Robust Tucker-Lewis Index (TLI)                            0.847
#RMSEA                                          0.069       0.054
#Robust RMSEA                                               0.063
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .330

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.759    0.043   17.500    0.000    0.759    0.759
#    factor3           0.711    0.052   13.619    0.000    0.711    0.711
#  factor2 ~~                                                            
#    factor3           0.649    0.053   12.277    0.000    0.649    0.649

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.693       0.707
#Tucker-Lewis Index (TLI)                       0.659       0.674
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.679
#RMSEA                                          0.096       0.078
#Robust RMSEA                                               0.091
#SRMR                                           0.193       0.193

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .313


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.895
#Tucker-Lewis Index (TLI)                       0.858       0.869
#Robust Comparative Fit Index (CFI)                         0.902
#Robust Tucker-Lewis Index (TLI)                            0.878
#RMSEA                                          0.062       0.050
#Robust RMSEA                                               0.056
#SRMR                                           0.050       0.050

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .312

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.01184079     0.63809524     4.14916441     0.80366093 
#
#$FactorLevelIndices
#            ECV_SS      ECV_SG      ECV_GS      Omega     OmegaH           H         FD
#factor1 0.36459150 0.002066731 0.635408501  0.7845238  0.2264502   0.5742528  0.7752010
#factor2 0.27842294 0.002203496 0.721577062  0.8134169  0.1691717   0.5732525  0.7508938
#factor3 0.99743700 0.983888986 0.002563003 53.5730990 53.0407126 470.3232279 26.7701261
#general 0.01184079 0.011840787 0.011840787  4.1491644  0.8036609   0.8748030  0.9150834








################################################################
### 
### Brand et al (2023) https://www.nature.com/articles/s44271-023-00016-x
### data Sample 1 Giessen https://osf.io/ex684?view_only=0bacd43b0c2c45e7b1c9d1227b33a94d


mydata <- mydata2
haven::zap_labels(mydata)

min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .87, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 5.2; eigenvalue 2 = .88

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 6.86; eigenvalue 2 = .97

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.084, RMSR=.07, TLI=.716

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.113, RMSR=.08, TLI=.678

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 3 factors (EFA-based)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.058, RMSR=.04, TLI=.862
#     MR1  MR3  MR2
#MR1 1.00 0.45 0.39
#MR3 0.45 1.00 0.49
#MR2 0.39 0.49 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.065, RMSR=.04, TLI=.888
#     MR1  MR3  MR2
#MR1 1.00 0.67 0.53
#MR3 0.67 1.00 0.47
#MR2 0.53 0.47 1.00


# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.950       0.879
#Tucker-Lewis Index (TLI)                       0.944       0.866
#Robust Comparative Fit Index (CFI)                         0.738
#Robust Tucker-Lewis Index (TLI)                            0.709
#RMSEA                                          0.081       0.087
#Robust RMSEA                                               0.108
#SRMR                                           0.079       0.079

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .351

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.743       0.756
#Tucker-Lewis Index (TLI)                       0.715       0.729
#Robust Comparative Fit Index (CFI)                         0.757
#Robust Tucker-Lewis Index (TLI)                            0.730
#RMSEA                                          0.085       0.072
#Robust RMSEA                                               0.082
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .258

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (EGA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.972       0.927
#Tucker-Lewis Index (TLI)                       0.969       0.918
#Robust Comparative Fit Index (CFI)                         0.828
#Robust Tucker-Lewis Index (TLI)                            0.806
#RMSEA                                          0.061       0.068
#Robust RMSEA                                               0.088
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .402

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.741    0.028   26.453    0.000    0.741    0.741
#    factor3           0.667    0.033   20.110    0.000    0.667    0.667
#  factor2 ~~                                                            
#    factor3           0.810    0.027   29.950    0.000    0.810    0.810

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.824       0.836
#Tucker-Lewis Index (TLI)                       0.802       0.815
#Robust Comparative Fit Index (CFI)                         0.839
#Robust Tucker-Lewis Index (TLI)                            0.818
#RMSEA                                          0.071       0.060
#Robust RMSEA                                               0.067
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .290

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.702    0.042   16.571    0.000    0.702    0.702
#    factor3           0.620    0.048   12.826    0.000    0.620    0.620
#  factor2 ~~                                                            
#    factor3           0.798    0.042   18.806    0.000    0.798    0.798


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (3 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.665       0.668
#Tucker-Lewis Index (TLI)                       0.627       0.631
#Robust Comparative Fit Index (CFI)                         0.675
#Robust Tucker-Lewis Index (TLI)                            0.639
#RMSEA                                          0.097       0.084
#Robust RMSEA                                               0.094
#SRMR                                           0.185       0.185

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .266


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.885       0.889
#Tucker-Lewis Index (TLI)                       0.857       0.862
#Robust Comparative Fit Index (CFI)                         0.896
#Robust Tucker-Lewis Index (TLI)                            0.870
#RMSEA                                          0.060       0.052
#Robust RMSEA                                               0.057
#SRMR                                           0.049       0.049

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .314

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    Burp              0.280    0.098    2.849    0.004    0.280    0.287
#    Cough             0.601    0.052   11.502    0.000    0.601    0.723
#    Wind              0.152    0.099    1.534    0.125    0.152    0.161
#    Sneeze            0.411    0.043    9.617    0.000    0.411    0.591
#    Vomit             0.233    0.070    3.338    0.001    0.233    0.255
#    Bruise           -0.007    0.086   -0.086    0.932   -0.007   -0.006
#  factor2 =~                                                            
#    Hungry            0.217    0.185    1.176    0.239    0.217    0.232
#    Thirsty           0.207    0.143    1.450    0.147    0.207    0.191
#    Urinate           0.360    0.062    5.791    0.000    0.360    0.516
#    Defecate          0.340    0.058    5.870    0.000    0.340    0.509
#    Breathing        -0.063    0.112   -0.559    0.576   -0.063   -0.066
#    Pain              0.029    0.127    0.226    0.821    0.029    0.039
#    Sex_arousal       0.145    0.092    1.584    0.113    0.145    0.173
#    Temp              0.057    0.080    0.707    0.479    0.057    0.084
#    Taste             0.064    0.071    0.897    0.370    0.064    0.085
#    Heart            -0.175    0.080   -2.191    0.028   -0.175   -0.196
#    Muscles           0.024    0.113    0.210    0.834    0.024    0.031
#  factor3 =~                                                            
#    Itch              0.268    0.105    2.545    0.011    0.268    0.416
#    Tickle            0.630    0.244    2.581    0.010    0.630    0.727
#    Affective_toch    0.121    0.058    2.106    0.035    0.121    0.134
#    Blood_Sugar      -0.033    0.079   -0.415    0.678   -0.033   -0.027
#  general =~                                                            
#    Urinate           0.259    0.044    5.885    0.000    0.259    0.371
#    Defecate          0.291    0.043    6.852    0.000    0.291    0.437
#    Hungry            0.330    0.073    4.543    0.000    0.330    0.353
#    Thirsty           0.427    0.067    6.351    0.000    0.427    0.394
#    Breathing         0.412    0.045    9.229    0.000    0.412    0.435
#    Heart             0.287    0.059    4.878    0.000    0.287    0.322
#    Pain              0.429    0.055    7.828    0.000    0.429    0.580
#    Taste             0.361    0.040    8.912    0.000    0.361    0.481
#    Temp              0.339    0.048    7.059    0.000    0.339    0.503
#    Muscles           0.429    0.051    8.460    0.000    0.429    0.560
#    Vomit             0.463    0.046   10.103    0.000    0.463    0.506
#    Sex_arousal       0.410    0.041    9.937    0.000    0.410    0.487
#    Sneeze            0.288    0.037    7.791    0.000    0.288    0.415
#    Wind              0.507    0.065    7.853    0.000    0.507    0.539
#    Affective_toch    0.478    0.044   10.734    0.000    0.478    0.528
#    Cough             0.401    0.046    8.753    0.000    0.401    0.482
#    Burp              0.522    0.072    7.273    0.000    0.522    0.536
#    Bruise            0.619    0.057   10.777    0.000    0.619    0.532
#    Tickle            0.427    0.053    8.122    0.000    0.427    0.492
#    Blood_Sugar       0.548    0.066    8.335    0.000    0.548    0.443
#    Itch              0.392    0.034   11.653    0.000    0.392    0.610

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6640252      0.6380952      0.8880634      0.8211221 
#$FactorLevelIndices
#        ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.4072377 0.14210094 0.5927623 0.7922965 0.24449499 0.6448539 0.8364880
#factor2 0.2368737 0.09594072 0.7631263 0.7695658 0.07339525 0.4737378 0.7014888
#factor3 0.3980499 0.09793316 0.6019501 0.7279009 0.19420896 0.5740961 0.8213367
#general 0.6640252 0.66402518 0.6640252 0.8880634 0.82112213 0.8682279 0.9174133





################################################################
### 
### Brand et al (2023) https://www.nature.com/articles/s44271-023-00016-x
### data Sample 2 Mainz https://osf.io/ex684?view_only=0bacd43b0c2c45e7b1c9d1227b33a94d

mydata <- mydata3
haven::zap_labels(mydata)

min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .89, omega T = .9

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 4 components
# Eigenvalue 1 = 5.78; eigenvalue 2 = .55

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 3 components
# Eigenvalue 1 = 7.23; eigenvalue 2 = .61

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.28, RMSEA=.073, RMSR=.05, TLI=.8

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.34, RMSEA=.092, RMSR=.06, TLI=.779

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 3 factors (EFA-based)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.34, RMSEA=.062, RMSR=.04, TLI=.858
#     MR1  MR3  MR2
#MR1 1.00 0.59 0.52
#MR3 0.59 1.00 0.48
#MR2 0.52 0.48 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.076, RMSR=.04, TLI=.855
#     MR1  MR3  MR2
#MR1 1.00 0.67 0.53
#MR3 0.67 1.00 0.47
#MR2 0.53 0.47 1.00


# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.964       0.863
#Tucker-Lewis Index (TLI)                       0.960       0.848
#Robust Comparative Fit Index (CFI)                         0.804
#Robust Tucker-Lewis Index (TLI)                            0.782
#RMSEA                                          0.077       0.091
#Robust RMSEA                                               0.092
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .381

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.820       0.825
#Tucker-Lewis Index (TLI)                       0.800       0.805
#Robust Comparative Fit Index (CFI)                         0.824
#Robust Tucker-Lewis Index (TLI)                            0.805
#RMSEA                                          0.073       0.063
#Robust RMSEA                                               0.072
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .268

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (EGA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.976       0.906
#Tucker-Lewis Index (TLI)                       0.973       0.894
#Robust Comparative Fit Index (CFI)                         0.856
#Robust Tucker-Lewis Index (TLI)                            0.838
#RMSEA                                          0.063       0.076
#Robust RMSEA                                               0.079
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .416

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.815    0.012   66.972    0.000    0.815    0.815
#    factor3           0.708    0.017   41.602    0.000    0.708    0.708
#  factor2 ~~                                                            
#    factor3           0.808    0.014   59.474    0.000    0.808    0.808


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.869       0.872
#Tucker-Lewis Index (TLI)                       0.852       0.856
#Robust Comparative Fit Index (CFI)                         0.873
#Robust Tucker-Lewis Index (TLI)                            0.857
#RMSEA                                          0.063       0.054
#Robust RMSEA                                               0.062
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .318

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.819    0.017   47.206    0.000    0.819    0.819
#    factor3           0.680    0.025   26.763    0.000    0.680    0.680
#  factor2 ~~                                                            
#    factor3           0.798    0.023   34.891    0.000    0.798    0.798

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (3 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.685       0.683
#Tucker-Lewis Index (TLI)                       0.650       0.648
#Robust Comparative Fit Index (CFI)                         0.688
#Robust Tucker-Lewis Index (TLI)                            0.653
#RMSEA                                          0.097       0.085
#Robust RMSEA                                               0.096
#SRMR                                           0.203       0.203

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .321


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.922       0.923
#Tucker-Lewis Index (TLI)                       0.902       0.903
#Robust Comparative Fit Index (CFI)                         0.925
#Robust Tucker-Lewis Index (TLI)                            0.907
#RMSEA                                          0.051       0.045
#Robust RMSEA                                               0.050
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .334

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7235284      0.6380952      0.9004254      0.8464274 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.2869019 0.09044078 0.7130981 0.7853373 0.17662997 0.4668592 0.7091775
#factor2 0.2109279 0.09907979 0.7890721 0.8171433 0.06873477 0.5251447 0.7553351
#factor3 0.4043582 0.08695107 0.5956418 0.6936120 0.19541039 0.5215092 0.7802096
#general 0.7235284 0.72352837 0.7235284 0.9004254 0.84642743 0.8856658 0.9288704





################################################################
### 
### Brand et al (2023) https://www.nature.com/articles/s44271-023-00016-x
### data Sample 3 PotVie https://osf.io/3f2h6/?view_only=0bacd43b0c2c45e7b1c9d1227b33a94d

mydata <- mydata4
haven::zap_labels(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .85, omega T = .87

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 4.71; eigenvalue 2 = .94

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 5.76; eigenvalue 2 = 1.09

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.093, RMSR=.08, TLI=.637

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.27, RMSEA=.116, RMSR=.09, TLI=.604

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 3 factors (EFA-based)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.065, RMSR=.05, TLI=.823
#     MR1  MR3  MR2
#MR1 1.00 0.50 0.25
#MR3 0.50 1.00 0.39
#MR2 0.25 0.39 1.00


fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.087, RMSR=.05, TLI=.777
#     MR1  MR3  MR2
#MR1 1.00 0.52 0.25
#MR3 0.52 1.00 0.41
#MR2 0.25 0.41 1.00


# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.890       0.737
#Tucker-Lewis Index (TLI)                       0.877       0.707
#Robust Comparative Fit Index (CFI)                         0.643
#Robust Tucker-Lewis Index (TLI)                            0.603
#RMSEA                                          0.117       0.121
#Robust RMSEA                                               0.117
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .272

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.673       0.681
#Tucker-Lewis Index (TLI)                       0.637       0.646
#Robust Comparative Fit Index (CFI)                         0.681
#Robust Tucker-Lewis Index (TLI)                            0.646
#RMSEA                                          0.094       0.083
#Robust RMSEA                                               0.092
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .195

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors (EGA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.954       0.881
#Tucker-Lewis Index (TLI)                       0.948       0.866
#Robust Comparative Fit Index (CFI)                         0.795
#Robust Tucker-Lewis Index (TLI)                            0.769
#RMSEA                                          0.076       0.082
#Robust RMSEA                                               0.089
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .348

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.699    0.028   25.331    0.000    0.699    0.699
#    factor3           0.528    0.034   15.538    0.000    0.528    0.528
#  factor2 ~~                                                            
#    factor3           0.463    0.036   12.747    0.000    0.463    0.463

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.823       0.832
#Tucker-Lewis Index (TLI)                       0.800       0.810
#Robust Comparative Fit Index (CFI)                         0.833
#Robust Tucker-Lewis Index (TLI)                            0.811
#RMSEA                                          0.070       0.061
#Robust RMSEA                                               0.067
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .255

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.698    0.041   16.917    0.000    0.698    0.698
#    factor3           0.486    0.050    9.812    0.000    0.486    0.486
#  factor2 ~~                                                            
#    factor3           0.360    0.051    6.995    0.000    0.360    0.360


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (2 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.721       0.726
#Tucker-Lewis Index (TLI)                       0.690       0.696
#Robust Comparative Fit Index (CFI)                         0.729
#Robust Tucker-Lewis Index (TLI)                            0.699
#RMSEA                                          0.087       0.077
#Robust RMSEA                                               0.085
#SRMR                                           0.157       0.157

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .257


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.879       0.872
#Tucker-Lewis Index (TLI)                       0.849       0.840
#Robust Comparative Fit Index (CFI)                         0.883
#Robust Tucker-Lewis Index (TLI)                            0.854
#RMSEA                                          0.060       0.056
#Robust RMSEA                                               0.059
#SRMR                                           0.049       0.049

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
#     0.5576113      0.6380952      0.8781142      0.7378012 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3786264 0.1306119 0.6213736 0.7927254 0.2403989 0.5783191 0.7663403
#factor2 0.3562608 0.1452384 0.6437392 0.7843040 0.2350194 0.5694594 0.7277768
#factor3 0.6732549 0.1665384 0.3267451 0.7098347 0.3956528 0.8925070 1.0005872
#general 0.5576113 0.5576113 0.5576113 0.8781142 0.7378012 0.8429439 0.8889301






################################################################
### 
### Campos et al (2023) https://recipp.ipp.pt/handle/10400.22/18657
### data https://github.com/RealityBending/InteroceptionIAS/blob/main/data/Campos2022/Dataset_Test.sav

mydata <- mydata5
haven::zap_labels(mydata)

min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .91, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 6.77; eigenvalue 2 = .97

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 8.74; eigenvalue 2 = 1.03

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.103, RMSR=.08, TLI=.713

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.131, RMSR=.08, TLI=.697

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 3 factors (EFA-based)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.073, RMSR=.04, TLI=.855
#     MR1  MR3  MR2
#MR1 1.00 0.55 0.55
#MR3 0.55 1.00 0.37
#MR2 0.55 0.37 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.102, RMSR=.04, TLI=.815
#     MR1  MR3  MR2
#MR1 1.00 0.63 0.58
#MR3 0.63 1.00 0.44
#MR2 0.58 0.44 1.00


# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.864
#Tucker-Lewis Index (TLI)                       0.952       0.848
#Robust Comparative Fit Index (CFI)                         0.731
#Robust Tucker-Lewis Index (TLI)                            0.701
#RMSEA                                          0.119       0.130
#Robust RMSEA                                               0.132
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .467

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.741       0.750
#Tucker-Lewis Index (TLI)                       0.712       0.722
#Robust Comparative Fit Index (CFI)                         0.754
#Robust Tucker-Lewis Index (TLI)                            0.727
#RMSEA                                          0.104       0.085
#Robust RMSEA                                               0.100
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .303

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (EGA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.929
#Tucker-Lewis Index (TLI)                       0.977       0.919
#Robust Comparative Fit Index (CFI)                         0.830
#Robust Tucker-Lewis Index (TLI)                            0.808
#RMSEA                                          0.083       0.095
#Robust RMSEA                                               0.105
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .518

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.819    0.018   46.568    0.000    0.819    0.819
#    factor3           0.661    0.029   22.627    0.000    0.661    0.661
#  factor2 ~~                                                            
#    factor3           0.644    0.027   23.561    0.000    0.644    0.644

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.843       0.853
#Tucker-Lewis Index (TLI)                       0.822       0.834
#Robust Comparative Fit Index (CFI)                         0.857
#Robust Tucker-Lewis Index (TLI)                            0.838
#RMSEA                                          0.082       0.066
#Robust RMSEA                                               0.077
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .387

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.789    0.034   23.507    0.000    0.789    0.789
#    factor3           0.605    0.045   13.349    0.000    0.605    0.605
#  factor2 ~~                                                            
#    factor3           0.552    0.047   11.759    0.000    0.552    0.552


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (3 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.727       0.731
#Tucker-Lewis Index (TLI)                       0.697       0.701
#Robust Comparative Fit Index (CFI)                         0.739
#Robust Tucker-Lewis Index (TLI)                            0.710
#RMSEA                                          0.107       0.088
#Robust RMSEA                                               0.103
#SRMR                                           0.232       0.232

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .403


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.899       0.906
#Tucker-Lewis Index (TLI)                       0.874       0.882
#Robust Comparative Fit Index (CFI)                         0.912
#Robust Tucker-Lewis Index (TLI)                            0.889
#RMSEA                                          0.069       0.055
#Robust RMSEA                                               0.064
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .415

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6796611      0.6380952      0.9229461      0.8240271 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2374410 0.0780483 0.7625590 0.8424342 0.1262533 0.4831557 0.7490677
#factor2 0.3082713 0.1404976 0.6917287 0.8563082 0.2104608 0.6275457 0.7822976
#factor3 0.4722824 0.1017930 0.5277176 0.7735342 0.3141580 0.6133025 0.8365469
#general 0.6796611 0.6796611 0.6796611 0.9229461 0.8240271 0.9028337 0.9286362






################################################################
### 
### Todd et al (2022) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0277894
### data  https://osf.io/ms354/

mydata <- mydata6
haven::zap_labels(mydata)

min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .89, omega T = .9

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 5.95; eigenvalue 2 = .86

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 7.29; eigenvalue 2 = .97

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.28, RMSEA=.091, RMSR=.07, TLI=.728

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.113, RMSR=.08, TLI=.702

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 3 factors (EFA-based)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.07, RMSR=.04, TLI=.838
#     MR1  MR3  MR2
#MR1 1.00 0.45 0.38
#MR3 0.45 1.00 0.18
#MR2 0.38 0.18 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.092, RMSR=.05, TLI=.803
#     MR1  MR3  MR2
#MR1 1.00 0.53 0.32
#MR3 0.53 1.00 0.16
#MR2 0.32 0.16 1.00

# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
#warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.941       0.796
#Tucker-Lewis Index (TLI)                       0.935       0.773
#Robust Comparative Fit Index (CFI)                         0.737
#Robust Tucker-Lewis Index (TLI)                            0.707
#RMSEA                                          0.108       0.120
#Robust RMSEA                                               0.113
#SRMR                                           0.079       0.079

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .374

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.755       0.763
#Tucker-Lewis Index (TLI)                       0.728       0.736
#Robust Comparative Fit Index (CFI)                         0.763
#Robust Tucker-Lewis Index (TLI)                            0.737
#RMSEA                                          0.091       0.078
#Robust RMSEA                                               0.089
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .296

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (EGA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.968       0.882
#Tucker-Lewis Index (TLI)                       0.963       0.866
#Robust Comparative Fit Index (CFI)                         0.823
#Robust Tucker-Lewis Index (TLI)                            0.801
#RMSEA                                          0.081       0.092
#Robust RMSEA                                               0.093
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .461

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.767    0.018   42.858    0.000    0.767    0.767
#    factor3           0.636    0.028   22.748    0.000    0.636    0.636
#  factor2 ~~                                                            
#    factor3           0.626    0.026   24.075    0.000    0.626    0.626

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.846       0.852
#Tucker-Lewis Index (TLI)                       0.826       0.833
#Robust Comparative Fit Index (CFI)                         0.854
#Robust Tucker-Lewis Index (TLI)                            0.835
#RMSEA                                          0.073       0.062
#Robust RMSEA                                               0.070
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .355

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.776    0.027   28.959    0.000    0.776    0.776
#    factor3           0.618    0.039   15.679    0.000    0.618    0.618
#  factor2 ~~                                                            
#    factor3           0.591    0.040   14.798    0.000    0.591    0.591


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (3 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.717       0.718
#Tucker-Lewis Index (TLI)                       0.685       0.686
#Robust Comparative Fit Index (CFI)                         0.723
#Robust Tucker-Lewis Index (TLI)                            0.693
#RMSEA                                          0.098       0.085
#Robust RMSEA                                               0.096
#SRMR                                           0.200       0.200

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .346


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.886
#Tucker-Lewis Index (TLI)                       0.857       0.857
#Robust Comparative Fit Index (CFI)                         0.892
#Robust Tucker-Lewis Index (TLI)                            0.865
#RMSEA                                          0.066       0.058
#Robust RMSEA                                               0.064
#SRMR                                           0.049       0.049

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .365

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    Burp              0.261    0.113    2.307    0.021    0.261    0.262
#    Cough             0.482    0.061    7.948    0.000    0.482    0.525
#    Wind              0.218    0.106    2.063    0.039    0.218    0.226
#    Sneeze            0.401    0.064    6.255    0.000    0.401    0.464
#    Vomit             0.193    0.060    3.241    0.001    0.193    0.191
#    Bruise            0.053    0.071    0.751    0.453    0.053    0.047
#  factor2 =~                                                            
#    Hungry            0.668    0.062   10.823    0.000    0.668    0.674
#    Thirsty           0.541    0.053   10.216    0.000    0.541    0.529
#    Urinate           0.282    0.081    3.504    0.000    0.282    0.308
#    Defecate          0.211    0.075    2.800    0.005    0.211    0.230
#    Breathing         0.250    0.077    3.256    0.001    0.250    0.266
#    Pain              0.078    0.060    1.286    0.198    0.078    0.093
#    Sex_arousal       0.116    0.046    2.538    0.011    0.116    0.138
#    Temp              0.161    0.063    2.567    0.010    0.161    0.196
#    Taste             0.166    0.056    2.950    0.003    0.166    0.184
#    Heart             0.149    0.085    1.748    0.081    0.149    0.154
#    Muscles           0.052    0.062    0.849    0.396    0.052    0.058
#  factor3 =~                                                            
#    Itch              0.635    0.060   10.532    0.000    0.635    0.631
#    Tickle            0.641    0.060   10.686    0.000    0.641    0.623
#    Affective_toch    0.318    0.054    5.890    0.000    0.318    0.290
#    Blood_Sugar       0.275    0.058    4.743    0.000    0.275    0.231
#  general =~                                                            
#    Urinate           0.459    0.043   10.704    0.000    0.459    0.500
#    Defecate          0.474    0.041   11.634    0.000    0.474    0.518
#    Hungry            0.352    0.052    6.784    0.000    0.352    0.356
#    Thirsty           0.352    0.049    7.156    0.000    0.352    0.344
#    Breathing         0.414    0.041    9.981    0.000    0.414    0.441
#    Heart             0.366    0.047    7.837    0.000    0.366    0.377
#    Pain              0.571    0.037   15.277    0.000    0.571    0.686
#    Taste             0.507    0.037   13.554    0.000    0.507    0.561
#    Temp              0.501    0.037   13.466    0.000    0.501    0.611
#    Muscles           0.605    0.037   16.219    0.000    0.605    0.666
#    Vomit             0.553    0.038   14.544    0.000    0.553    0.546
#    Sex_arousal       0.440    0.032   13.605    0.000    0.440    0.526
#    Sneeze            0.431    0.036   11.869    0.000    0.431    0.500
#    Wind              0.563    0.043   12.957    0.000    0.563    0.583
#    Affective_toch    0.564    0.040   14.262    0.000    0.564    0.515
#    Cough             0.525    0.038   13.658    0.000    0.525    0.572
#    Burp              0.541    0.046   11.637    0.000    0.541    0.545
#    Bruise            0.455    0.049    9.349    0.000    0.455    0.400
#    Tickle            0.463    0.041   11.243    0.000    0.463    0.450
#    Blood_Sugar       0.337    0.051    6.561    0.000    0.337    0.283
#    Itch              0.466    0.043   10.887    0.000    0.466    0.463

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6709941      0.6380952      0.9049973      0.8013466 
#$FactorLevelIndices
#        ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2797283 0.08048589 0.7202717 0.7772857 0.1780108 0.4513516 0.6947926
#factor2 0.2661720 0.13394876 0.7338280 0.8496821 0.1735946 0.6141336 0.7727803
#factor3 0.5481097 0.11457124 0.4518903 0.7244474 0.3757119 0.5910600 0.7987801
#general 0.6709941 0.67099411 0.6709941 0.9049973 0.8013466 0.8871230 0.9203340





################################################################
### 
### Gaggero et al (2021) https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0261126
### data  https://osf.io/5x9sg

mydata <- mydata7
haven::zap_labels(mydata)

min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .88, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 2 components
# Eigenvalue 1 = 5.48; eigenvalue 2 = .78

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 8 factors and 2 components
# Eigenvalue 1 = 6.86; eigenvalue 2 = .9

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.26, RMSEA=.074, RMSR=.06, TLI=.783

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.094, RMSR=.07, TLI=.757

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communites if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities

# Give solution with 3 factors (EFA-based)
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.056, RMSR=.04, TLI=.876
#     MR1  MR3  MR2
#MR1 1.00 0.57 0.33
#MR3 0.57 1.00 0.37
#MR2 0.33 0.37 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.076, RMSR=.04, TLI=.842
#     MR1  MR3  MR2
#MR1 1.00 0.58 0.32
#MR3 0.58 1.00 0.21
#MR2 0.32 0.21 1.00

# Single factor model lavaan
UNImodel= '
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.959       0.855
#Tucker-Lewis Index (TLI)                       0.955       0.839
#Robust Comparative Fit Index (CFI)                         0.788
#Robust Tucker-Lewis Index (TLI)                            0.765
#RMSEA                                          0.079       0.090
#Robust RMSEA                                               0.093
#SRMR                                           0.068       0.068

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .33

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.804       0.816
#Tucker-Lewis Index (TLI)                       0.783       0.796
#Robust Comparative Fit Index (CFI)                         0.816
#Robust Tucker-Lewis Index (TLI)                            0.795
#RMSEA                                          0.074       0.062
#Robust RMSEA                                               0.071
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .263

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors (EGA based) 
EGAmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.910
#Tucker-Lewis Index (TLI)                       0.973       0.898
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.831
#RMSEA                                          0.061       0.072
#Robust RMSEA                                               0.079
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .369

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.864    0.016   54.542    0.000    0.864    0.864
#    factor3           0.658    0.032   20.316    0.000    0.658    0.658
#  factor2 ~~                                                            
#    factor3           0.593    0.032   18.413    0.000    0.593    0.593

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.865       0.876
#Tucker-Lewis Index (TLI)                       0.847       0.860
#Robust Comparative Fit Index (CFI)                         0.877
#Robust Tucker-Lewis Index (TLI)                            0.861
#RMSEA                                          0.062       0.051
#Robust RMSEA                                               0.059
#SRMR                                           0.051       0.051

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .297

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.864    0.016   54.542    0.000    0.864    0.864
#    factor3           0.658    0.032   20.316    0.000    0.658    0.658
#  factor2 ~~                                                            
#    factor3           0.593    0.032   18.413    0.000    0.593    0.593


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (3 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.696       0.701
#Tucker-Lewis Index (TLI)                       0.662       0.667
#Robust Comparative Fit Index (CFI)                         0.705
#Robust Tucker-Lewis Index (TLI)                            0.672
#RMSEA                                          0.092       0.079
#Robust RMSEA                                               0.090
#SRMR                                           0.188       0.188

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .313


# Bifactor model with 3 factors
BIFmodel= '
 factor1 =~ Burp+Cough+Wind+Sneeze+Vomit+Bruise
 factor2 =~ Hungry+Thirsty+Urinate+Defecate+Breathing+Pain+Sex_arousal+Temp+Taste+Heart+Muscles
 factor3 =~ Itch+Tickle+Affective_touch+Blood_Sugar
 general =~ Urinate+Defecate+Hungry+Thirsty+Breathing+Heart+Pain+Taste+Temp+
            Muscles+Vomit+Sex_arousal+Sneeze+Wind+Affective_touch+Cough+Burp+
            Bruise+Tickle+Blood_Sugar+Itch
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.895
#Tucker-Lewis Index (TLI)                       0.867       0.869
#Robust Comparative Fit Index (CFI)                         0.902
#Robust Tucker-Lewis Index (TLI)                            0.877
#RMSEA                                          0.058       0.049
#Robust RMSEA                                               0.055
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .314

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7147285      0.6380952      0.8922346      0.8250932 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2224276 0.07071369 0.7775724 0.7698736 0.1109799 0.3954695 0.6622691
#factor2 0.2192411 0.10615708 0.7807589 0.8179806 0.1036646 0.4974147 0.7148926
#factor3 0.5478112 0.10840072 0.4521888 0.6624730 0.3314570 0.5622984 0.7824257
#general 0.7147285 0.71472851 0.7147285 0.8922346 0.8250932 0.8764305 0.9207254






