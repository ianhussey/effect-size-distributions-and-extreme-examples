################################################################
### Analysis PANAS 20 items
### Spanish study Diaz-Garcia https://zenodo.org/record/1477084

library(haven)
PANAS_Diaz_Garcia <- read_sav("PANAS_Diaz-Garcia.sav")
colnames(PANAS_Diaz_Garcia)
mydata  <- as.data.frame(PANAS_Diaz_Garcia[,12:31])
print(colnames(mydata))
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 6.19; eigenvalue 2 = 2.86

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 7.12; eigenvalue 2 = 3.1

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.31, RMSEA=.168, RMSR=.16, TLI=.472

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.196, RMSR=.18, TLI=.457

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities 

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.11, RMSR=.06, TLI=.772
#      MR1   MR2
#MR1  1.00 -0.27
#MR2 -0.27  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.137, RMSR=.07, TLI=.732
#      MR1   MR2
#MR1  1.00 -0.29
#MR2 -0.29  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ PANAS.R1pre+PANAS.R2pre+PANAS.R3pre+PANAS.R4pre+PANAS.R5pre+PANAS.R6pre+PANAS.R7pre+
           PANAS.R8pre+PANAS.R9pre+PANAS.R10pre+PANAS.R11pre+PANAS.R12pre+PANAS.R13pre+PANAS.R14pre+
           PANAS.R15pre+PANAS.R16pre+PANAS.R17pre+PANAS.R18pre+PANAS.R19pre+PANAS.R20pre
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.848       0.687
#Tucker-Lewis Index (TLI)                       0.830       0.650
#Robust Comparative Fit Index (CFI)                         0.457
#Robust Tucker-Lewis Index (TLI)                            0.393
#RMSEA                                          0.256       0.217
#Robust RMSEA                                               0.209
#SRMR                                           0.201       0.201

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .43

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.547       0.529
#Tucker-Lewis Index (TLI)                       0.494       0.473
#Robust Comparative Fit Index (CFI)                         0.549
#Robust Tucker-Lewis Index (TLI)                            0.496
#RMSEA                                          0.166       0.152
#Robust RMSEA                                               0.165
#SRMR                                           0.166       0.166

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .259

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 pos =~ PANAS.R9pre+PANAS.R19pre+PANAS.R5pre+PANAS.R3pre+PANAS.R1pre+PANAS.R14pre+PANAS.R16pre+
           PANAS.R12pre+PANAS.R10pre+PANAS.R17pre
 neg =~ PANAS.R15pre+PANAS.R18pre+PANAS.R7pre+PANAS.R2pre+
           PANAS.R20pre+PANAS.R4pre+PANAS.R11pre+PANAS.R8pre+PANAS.R6pre+PANAS.R13pre
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.968       0.917
#Tucker-Lewis Index (TLI)                       0.964       0.907
#Robust Comparative Fit Index (CFI)                         0.767
#Robust Tucker-Lewis Index (TLI)                            0.738
#RMSEA                                          0.118       0.112
#Robust RMSEA                                               0.137
#SRMR                                           0.090       0.090

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#pos ~~                                                                
#  neg              -0.328    0.035   -9.356    0.000   -0.328   -0.328

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .55

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.804       0.796
#Tucker-Lewis Index (TLI)                       0.779       0.771
#Robust Comparative Fit Index (CFI)                         0.807
#Robust Tucker-Lewis Index (TLI)                            0.783
#RMSEA                                          0.110       0.100
#Robust RMSEA                                               0.108
#SRMR                                           0.074       0.074

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#pos ~~                                                                
#  neg              -0.315    0.044   -7.186    0.000   -0.315   -0.315

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .439


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.796       0.788
#Tucker-Lewis Index (TLI)                       0.772       0.763
#Robust Comparative Fit Index (CFI)                         0.799
#Robust Tucker-Lewis Index (TLI)                            0.775
#RMSEA                                          0.111       0.102
#Robust RMSEA                                               0.110
#SRMR                                           0.130       0.130

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .438

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 pos =~ PANAS.R9pre+PANAS.R19pre+PANAS.R5pre+PANAS.R3pre+PANAS.R1pre+PANAS.R14pre+PANAS.R16pre+
           PANAS.R12pre+PANAS.R10pre+PANAS.R17pre
 neg =~ PANAS.R15pre+PANAS.R18pre+PANAS.R7pre+PANAS.R2pre+
           PANAS.R20pre+PANAS.R4pre+PANAS.R11pre+PANAS.R8pre+PANAS.R6pre+PANAS.R13pre
 general=~ PANAS.R1pre+PANAS.R2pre+PANAS.R3pre+PANAS.R4pre+PANAS.R5pre+PANAS.R6pre+PANAS.R7pre+
           PANAS.R8pre+PANAS.R9pre+PANAS.R10pre+PANAS.R11pre+PANAS.R12pre+PANAS.R13pre+PANAS.R14pre+
           PANAS.R15pre+PANAS.R16pre+PANAS.R17pre+PANAS.R18pre+PANAS.R19pre+PANAS.R20pre
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.885       0.881
#Tucker-Lewis Index (TLI)                       0.855       0.849
#Robust Comparative Fit Index (CFI)                         0.888
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.089       0.081
#Robust RMSEA                                               0.087
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .504

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4510692      0.5263158      0.8647986      0.2065405 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#pos     0.8777717 0.4435793 0.1222283 0.9157086 0.80635271 0.9031565 0.9509596
#neg     0.2129807 0.1053515 0.7870193 0.8945169 0.02937797 0.6156237 0.8787653
#general 0.4510692 0.4510692 0.4510692 0.8647986 0.20654046 0.8880724 0.9485943









################################################################
### Analysis PANAS 20 items
### Ceccato et al. 2021 https://www.sciencedirect.com/science/article/pii/S2352340921001761

### Not included because dataset already used for other variable

library(readxl)
Ceccato_et_al_2021_dataset <- read_excel("Ceccato et  al_2021_dataset.xlsx")
mydata  <- as.data.frame(Ceccato_et_al_2021_dataset[,33:52])
print(colnames(mydata))
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 components
# Eigenvalue 1 = 4.73; eigenvalue 2 = 2.25

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 7.12; eigenvalue 2 = 3.1

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.24, RMSEA=.145, RMSR=.15, TLI=.512

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.182, RMSR=.18, TLI=.441

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities 

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.77, RMSEA=.112, RMSR=.08, TLI=.71
#      MR1   MR2
#MR1  1.00 -0.17
#MR2 -0.17  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.153, RMSR=.11, TLI=.605
#      MR1   MR2
#MR1  1.00 -0.18
#MR2 -0.18  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ PANAS_1+PANAS_2+PANAS_3+PANAS_4+PANAS_5+PANAS_6+PANAS_7+PANAS_8+PANAS_9+PANAS_10+
           PANAS_11+PANAS_12+PANAS_13+PANAS_14+PANAS_15+PANAS_16+PANAS_17+PANAS_18+PANAS_19+PANAS_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.856       0.735
#Tucker-Lewis Index (TLI)                       0.839       0.704
#Robust Comparative Fit Index (CFI)                         0.473
#Robust Tucker-Lewis Index (TLI)                            0.410
#RMSEA                                          0.214       0.200
#Robust RMSEA                                               0.187
#SRMR                                           0.179       0.179

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .13

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 pos =~ PANAS_9+PANAS_19+PANAS_5+PANAS_3+PANAS_1+PANAS_14+PANAS_16+PANAS_12+PANAS_10+PANAS_17
 neg =~ PANAS_15+PANAS_18+PANAS_7+PANAS_2+PANAS_20+PANAS_4+PANAS_11+PANAS_8+PANAS_6+PANAS_13
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.945       0.907
#Tucker-Lewis Index (TLI)                       0.938       0.896
#Robust Comparative Fit Index (CFI)                         0.646
#Robust Tucker-Lewis Index (TLI)                            0.602
#RMSEA                                          0.133       0.119
#Robust RMSEA                                               0.154
#SRMR                                           0.126       0.126

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#pos ~~                                                                
#        neg              -0.148    0.016   -9.274    0.000   -0.148   -0.148

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .44


CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.939       0.928
#Tucker-Lewis Index (TLI)                       0.932       0.920
#Robust Comparative Fit Index (CFI)                         0.644
#Robust Tucker-Lewis Index (TLI)                            0.602
#RMSEA                                          0.139       0.104
#Robust RMSEA                                               0.154
#SRMR                                           0.133       0.133

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .418

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 pos =~ PANAS_9+PANAS_19+PANAS_5+PANAS_3+PANAS_1+PANAS_14+PANAS_16+PANAS_12+PANAS_10+PANAS_17
 neg =~ PANAS_15+PANAS_18+PANAS_7+PANAS_2+PANAS_20+PANAS_4+PANAS_11+PANAS_8+PANAS_6+PANAS_13
 general=~ PANAS_1+PANAS_2+PANAS_3+PANAS_4+PANAS_5+PANAS_6+PANAS_7+PANAS_8+PANAS_9+PANAS_10+
           PANAS_11+PANAS_12+PANAS_13+PANAS_14+PANAS_15+PANAS_16+PANAS_17+PANAS_18+PANAS_19+PANAS_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .540

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  pos =~                                                                
#    PANAS_9           0.364    0.025   14.416    0.000    0.364    0.364
#    PANAS_19          0.658    0.010   62.801    0.000    0.658    0.658
#    PANAS_5           0.473    0.014   33.454    0.000    0.473    0.473
#    PANAS_3           0.197    0.023    8.403    0.000    0.197    0.197
#    PANAS_1           0.389    0.016   23.596    0.000    0.389    0.389
#    PANAS_14          0.447    0.017   26.098    0.000    0.447    0.447
#    PANAS_16          0.821    0.008  109.119    0.000    0.821    0.821
#    PANAS_12          0.722    0.009   81.424    0.000    0.722    0.722
#    PANAS_10          0.415    0.019   22.207    0.000    0.415    0.415
#    PANAS_17          0.698    0.011   64.373    0.000    0.698    0.698
#  neg =~                                                                
#    PANAS_15          0.902    0.005  193.667    0.000    0.902    0.902
#    PANAS_18          0.839    0.007  122.680    0.000    0.839    0.839
#    PANAS_7           0.701    0.013   53.958    0.000    0.701    0.701
#    PANAS_2           0.719    0.010   68.458    0.000    0.719    0.719
#    PANAS_20          0.732    0.012   60.193    0.000    0.732    0.732
#    PANAS_4           0.708    0.011   65.979    0.000    0.708    0.708
#    PANAS_11          0.849    0.007  119.499    0.000    0.849    0.849
#    PANAS_8           0.520    0.014   36.099    0.000    0.520    0.520
#    PANAS_6           0.393    0.021   18.681    0.000    0.393    0.393
#    PANAS_13          0.456    0.021   21.415    0.000    0.456    0.456
#  general =~                                                            
#    PANAS_1           0.478    0.019   24.776    0.000    0.478    0.478
#    PANAS_2           0.423    0.017   24.213    0.000    0.423    0.423
#    PANAS_3          -0.452    0.020  -22.871    0.000   -0.452   -0.452
#    PANAS_4           0.407    0.017   23.532    0.000    0.407    0.407
#    PANAS_5          -0.306    0.019  -15.912    0.000   -0.306   -0.306
#    PANAS_6          -0.202    0.025   -8.181    0.000   -0.202   -0.202
#    PANAS_7           0.574    0.016   35.476    0.000    0.574    0.574
#    PANAS_8          -0.204    0.020  -10.032    0.000   -0.204   -0.204
#    PANAS_9          -0.762    0.018  -41.750    0.000   -0.762   -0.762
#    PANAS_10         -0.396    0.019  -20.376    0.000   -0.396   -0.396
#    PANAS_11         -0.167    0.020   -8.313    0.000   -0.167   -0.167
#    PANAS_12         -0.036    0.022   -1.646    0.100   -0.036   -0.036
#    PANAS_13         -0.237    0.027   -8.922    0.000   -0.237   -0.237
#    PANAS_14         -0.425    0.019  -22.803    0.000   -0.425   -0.425
#    PANAS_15         -0.042    0.021   -2.007    0.045   -0.042   -0.042
#    PANAS_16         -0.111    0.023   -4.922    0.000   -0.111   -0.111
#    PANAS_17          0.264    0.022   12.113    0.000    0.264    0.264
#    PANAS_18          0.243    0.019   12.473    0.000    0.243    0.243
#    PANAS_19         -0.111    0.021   -5.236    0.000   -0.111   -0.111
#    PANAS_20          0.537    0.017   31.604    0.000    0.537    0.537

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.255685960    0.526315789    0.871980039    0.005617926 

#$FactorLevelIndices
#ECV_SS    ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#pos     0.6631743 0.2839383 0.3368257 0.8126320 0.736913880 0.8567695 0.9304877
#neg     0.8050643 0.4603758 0.1949357 0.8863949 0.856879500 0.9348628 0.9705510
#general 0.2556860 0.2556860 0.2556860 0.8719800 0.005617926 0.8041042 0.9320953





################################################################
### Analysis PANAS 20 items
### Samo et al https://osf.io/k6dtz
### https://osf.io/9pwky

### not included because dataset already used for other variable

bfas_study4_various_tests <- read.csv("~/Self regulation/bfas_study4_various_tests.csv")
colnames(bfas_study4_various_tests)
mydata  <- as.data.frame(bfas_study4_various_tests[,3:22])
colnames(mydata)
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata)

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 components
# Eigenvalue 1 = 8.32; eigenvalue 2 = 2.45

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 9.79; eigenvalue 2 = 2.74

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.42, RMSEA=.174, RMSR=.14, TLI=.568

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.221, RMSR=.16, TLI=.517

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias

# Give solution with 2 factors (theory-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.109, RMSR=.05, TLI=.83
#      MR1   MR2
#MR1  1.00 -0.47
#MR2 -0.47  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.64, RMSEA=.153, RMSR=.05, TLI=.768
#      MR1   MR2
#MR1  1.00 -0.49
#MR2 -0.49  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ PANAS_1+PANAS_2+PANAS_3+PANAS_4+PANAS_5+PANAS_6+PANAS_7+PANAS_8+PANAS_9+PANAS_10+
           PANAS_11+PANAS_12+PANAS_13+PANAS_14+PANAS_15+PANAS_16+PANAS_17+PANAS_18+PANAS_19+PANAS_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.938       0.831
#Tucker-Lewis Index (TLI)                       0.931       0.811
#Robust Comparative Fit Index (CFI)                         0.550
#Robust Tucker-Lewis Index (TLI)                            0.497
#RMSEA                                          0.221       0.189
#Robust RMSEA                                               0.228
#SRMR                                           0.174       0.174

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .551

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 pos =~ PANAS_9+PANAS_19+PANAS_5+PANAS_3+PANAS_1+PANAS_14+PANAS_16+PANAS_12+PANAS_10+PANAS_17
 neg =~ PANAS_15+PANAS_18+PANAS_7+PANAS_2+PANAS_20+PANAS_4+PANAS_11+PANAS_8+PANAS_6+PANAS_13
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.988       0.949
#Tucker-Lewis Index (TLI)                       0.986       0.943
#Robust Comparative Fit Index (CFI)                         0.811
#Robust Tucker-Lewis Index (TLI)                            0.788
#RMSEA                                          0.099       0.104
#Robust RMSEA                                               0.148
#SRMR                                           0.074       0.074

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#pos ~~                                                                
#        neg              -0.544    0.033  -16.648    0.000   -0.544   -0.544

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .665

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.887       0.849
#Tucker-Lewis Index (TLI)                       0.873       0.832
#Robust Comparative Fit Index (CFI)                         0.795
#Robust Tucker-Lewis Index (TLI)                            0.771
#RMSEA                                          0.300       0.178
#Robust RMSEA                                               0.154
#SRMR                                           0.251       0.251

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .665


# Bifactor model
BIFmodel= '
 pos =~ PANAS_9+PANAS_19+PANAS_5+PANAS_3+PANAS_1+PANAS_14+PANAS_16+PANAS_12+PANAS_10+PANAS_17
 neg =~ PANAS_15+PANAS_18+PANAS_7+PANAS_2+PANAS_20+PANAS_4+PANAS_11+PANAS_8+PANAS_6+PANAS_13
 general=~ PANAS_1+PANAS_2+PANAS_3+PANAS_4+PANAS_5+PANAS_6+PANAS_7+PANAS_8+PANAS_9+PANAS_10+
           PANAS_11+PANAS_12+PANAS_13+PANAS_14+PANAS_15+PANAS_16+PANAS_17+PANAS_18+PANAS_19+PANAS_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .699

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  pos =~                                                                
#    PANAS_9           0.445    0.061    7.319    0.000    0.445    0.445
#    PANAS_19          0.172    0.059    2.895    0.004    0.172    0.172
#    PANAS_5           0.153    0.061    2.493    0.013    0.153    0.153
#    PANAS_3           0.606    0.058   10.367    0.000    0.606    0.606
#    PANAS_1           0.098    0.064    1.529    0.126    0.098    0.098
#    PANAS_14          0.349    0.062    5.650    0.000    0.349    0.349
#    PANAS_16          0.037    0.063    0.590    0.555    0.037    0.037
#    PANAS_12         -0.307    0.069   -4.448    0.000   -0.307   -0.307
#    PANAS_10          0.346    0.058    5.985    0.000    0.346    0.346
#    PANAS_17         -0.458    0.085   -5.409    0.000   -0.458   -0.458
#  neg =~                                                                
#    PANAS_15          0.674    0.028   23.676    0.000    0.674    0.674
#    PANAS_18          0.636    0.031   20.701    0.000    0.636    0.636
#    PANAS_7           0.824    0.021   38.898    0.000    0.824    0.824
#    PANAS_2           0.705    0.025   27.989    0.000    0.705    0.705
#    PANAS_20          0.840    0.021   40.294    0.000    0.840    0.840
#    PANAS_4           0.711    0.026   27.867    0.000    0.711    0.711
#    PANAS_11          0.610    0.028   21.898    0.000    0.610    0.610
#    PANAS_8           0.600    0.034   17.405    0.000    0.600    0.600
#    PANAS_6           0.722    0.027   26.666    0.000    0.722    0.722
#    PANAS_13          0.705    0.029   24.389    0.000    0.705    0.705
#  general =~                                                            
#    PANAS_1           0.825    0.019   42.475    0.000    0.825    0.825
#    PANAS_2          -0.571    0.034  -16.743    0.000   -0.571   -0.571
#    PANAS_3           0.640    0.053   12.057    0.000    0.640    0.640
#    PANAS_4          -0.547    0.035  -15.479    0.000   -0.547   -0.547
#    PANAS_5           0.806    0.025   32.169    0.000    0.806    0.806
#    PANAS_6          -0.488    0.042  -11.701    0.000   -0.488   -0.488
#    PANAS_7          -0.396    0.044   -9.072    0.000   -0.396   -0.396
#    PANAS_8          -0.488    0.043  -11.411    0.000   -0.488   -0.488
#    PANAS_9           0.781    0.038   20.608    0.000    0.781    0.781
#    PANAS_10          0.711    0.037   19.289    0.000    0.711    0.711
#    PANAS_11         -0.529    0.035  -15.283    0.000   -0.529   -0.529
#    PANAS_12          0.635    0.042   15.247    0.000    0.635    0.635
#    PANAS_13         -0.505    0.040  -12.693    0.000   -0.505   -0.505
#    PANAS_14          0.736    0.034   21.722    0.000    0.736    0.736
#    PANAS_15         -0.505    0.035  -14.276    0.000   -0.505   -0.505
#    PANAS_16          0.787    0.023   33.848    0.000    0.787    0.787
#    PANAS_17          0.791    0.042   18.860    0.000    0.791    0.791
#    PANAS_18         -0.354    0.043   -8.315    0.000   -0.354   -0.354
#    PANAS_19          0.758    0.026   29.435    0.000    0.758    0.758
#    PANAS_20         -0.380    0.042   -8.945    0.000   -0.380   -0.380

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.5626397      0.5263158      0.8803844      0.1168709 

#$FactorLevelIndices
#        ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#pos     0.1727931 0.0832390 0.8272069 0.9301665 0.05165897 0.6058019 0.9041960
#neg     0.6832709 0.3541213 0.3167291 0.9503348 0.64320227 0.9183973 0.9689883
#general 0.5626397 0.5626397 0.5626397 0.8803844 0.11687095 0.9440799 0.9744856




