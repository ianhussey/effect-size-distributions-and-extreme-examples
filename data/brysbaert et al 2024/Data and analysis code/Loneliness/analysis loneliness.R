################################################################
### Analysis loneliness as example of unifactorial scale
###
### Three studies 
###   1. Lardone et al. (2020) https://peerj.com/articles/10611/
###   2. Weinstein & Nguyen (2020) https://royalsocietypublishing.org/doi/full/10.1098/rsos.200458
###   3. Panayiotou et al. (2023) https://journals.sagepub.com/doi/full/10.1177/10731911221119533
###   4. McDanal et al. (2021) https://archive.org/details/osf-registrations-52ctd-v1


################################################################
###
### First determine likely factor structure on the basis of all datasets

library(readxl)
Lardone <- read_excel("Lardone.xlsx")
colnames(Lardone)
mydata <- as.data.frame(Lardone[,109:128])
colnames(mydata)
mydata1 <- na.omit(mydata)
min(mydata1)
max(mydata1) # 5 response alternatives

Weinstein <- read.csv("Weinstein.csv")
colnames(Weinstein)
mydata <- as.data.frame(Weinstein[,50:69])
colnames(mydata)
mydata2 <- na.omit(mydata)
min(mydata2)
max(mydata2) # 7 response alternatives

Panayiotou <- read_excel("Panayiotou.xlsx")
colnames(Panayiotou)
mydata <- as.data.frame(Panayiotou[,2:21])
colnames(mydata)
mydata3 <- na.omit(mydata)
min(mydata3)
max(mydata3) # 4 response alternatives

McDanal <- read.csv("McDanal.csv")
colnames(McDanal)
mydata <- as.data.frame(McDanal[,30:49])
colnames(mydata)
mydata4 <- na.omit(mydata)
min(mydata4)
max(mydata4) # 4 response alternatives

fit1 <- fa.parallel(mydata1)
fit1
# 3 components
fit2 <- fa.parallel(mydata2)
fit2
# 1 component
fit3 <- fa.parallel(mydata3)
fit3
# 3 components
fit4 <- fa.parallel(mydata4)
fit4
# 3 components

colnames(mydata2) <- colnames(mydata1)
colnames(mydata3) <- colnames(mydata1)
colnames(mydata4) <- colnames(mydata1)
mydata <- rbind(mydata1,mydata2,mydata3,mydata4)

fit1 <- fa.parallel(mydata)
fit1
# 2 components

fit2 <- fa(mydata,2)
fit2
diagram(fit2)

# factor 1 = items 20+19+10+6+16+1+5+9+15
# factor 2 = items 14+11+18+17+2+13+12+7+3+8+4


################################################################
### Lardone et al. (2020) https://peerj.com/articles/10611/
### Data https://osf.io/j9cnq

mydata <- mydata1

library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 4.67; eigenvalue 2 = 1.83

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 5.98; eigenvalue 2 = 2.33

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.23, RMSEA=.106, RMSR=.12, TLI=.607

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.166, RMSR=.15, TLI=.47

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities with response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.34, RMSEA=.073, RMSR=.06, TLI=.809
#     MR1  MR2
#MR1 1.00 -.26
#MR2 -.26 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.136, RMSR=.08, TLI=.643
#     MR1  MR2
#MR1 1.00 -.28
#MR2 -.28 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.873       0.753
#Tucker-Lewis Index (TLI)                       0.858       0.723
#Robust Comparative Fit Index (CFI)                         0.569
#Robust Tucker-Lewis Index (TLI)                            0.518
#RMSEA                                          0.147       0.142
#Robust RMSEA                                               0.161
#SRMR                                           0.144       0.144

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .378

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.648       0.639
#Tucker-Lewis Index (TLI)                       0.607       0.596
#Robust Comparative Fit Index (CFI)                         0.650
#Robust Tucker-Lewis Index (TLI)                            0.608
#RMSEA                                          0.108       0.106
#Robust RMSEA                                               0.108
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .208

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.946       0.890
#Tucker-Lewis Index (TLI)                       0.940       0.876
#Robust Comparative Fit Index (CFI)                         0.756
#Robust Tucker-Lewis Index (TLI)                            0.725
#RMSEA                                          0.096       0.095
#Robust RMSEA                                               0.122
#SRMR                                           0.109       0.109

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#        factor2          -0.433    0.054   -7.983    0.000   -0.433   -0.433

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .446

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.815       0.801
#Tucker-Lewis Index (TLI)                       0.792       0.777
#Robust Comparative Fit Index (CFI)                         0.814
#Robust Tucker-Lewis Index (TLI)                            0.790
#RMSEA                                          0.079       0.079
#Robust RMSEA                                               0.079
#SRMR                                           0.084       0.084

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#        factor2           -0.417    0.098   -4.270    0.000   -0.417   -0.417

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .344

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.796       0.784
#Tucker-Lewis Index (TLI)                       0.772       0.759
#Robust Comparative Fit Index (CFI)                         0.795
#Robust Tucker-Lewis Index (TLI)                            0.771
#RMSEA                                          0.083       0.082
#Robust RMSEA                                               0.082
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .368

# Bifactor model
BIFmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.916       0.909
#Tucker-Lewis Index (TLI)                       0.894       0.885
#Robust Comparative Fit Index (CFI)                         0.915
#Robust Tucker-Lewis Index (TLI)                            0.893
#RMSEA                                          0.056       0.056
#Robust RMSEA                                               0.056
#SRMR                                           0.058       0.058

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .397

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#  0.3578095057   0.5210526316   0.7851876297   0.0006937443 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega       OmegaH         H        FD
#factor1 0.4797876 0.2157189 0.5202124 0.7914766 0.4276001281 0.7314475 0.8496302
#factor2 0.7748583 0.4264716 0.2251417 0.8614041 0.6879183680 0.8341454 0.9052604
#general 0.3578095 0.3578095 0.3578095 0.7851876 0.0006937443 0.8415013 0.9022048






################################################################
### Weinstein & Nguyen 

mydata <- mydata2

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 1 component
# Eigenvalue 1 = 13.09; eigenvalue 2 = .58

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 1 component
# Eigenvalue 1 = 14.39; eigenvalue 2 = .53

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.65, RMSEA=.116, RMSR=.04, TLI=.874

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.72, RMSEA=.129, RMSR=.04, TLI=.872

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communities with response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 1 community

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.094, RMSR=.03, TLI=.917
#     MR1  MR2
#MR1 1.00  .28
#MR2  .28 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.75, RMSEA=.107, RMSR=.03, TLI=.913
#     MR1  MR2
#MR1 1.00  .13
#MR2  .13 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.997       0.971
#Tucker-Lewis Index (TLI)                       0.997       0.967
#Robust Comparative Fit Index (CFI)                         0.888
#Robust Tucker-Lewis Index (TLI)                            0.875
#RMSEA                                          0.082       0.113
#Robust RMSEA                                               0.129
#SRMR                                           0.040       0.040

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .737

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.888       0.890
#Tucker-Lewis Index (TLI)                       0.874       0.877
#Robust Comparative Fit Index (CFI)                         0.893
#Robust Tucker-Lewis Index (TLI)                            0.881
#RMSEA                                          0.116       0.088
#Robust RMSEA                                               0.113
#SRMR                                           0.042       0.042

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .654

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 2 factors
EGAmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered="TRUE",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.819
#Tucker-Lewis Index (TLI)                       1.001       0.797
#Robust Comparative Fit Index (CFI)                         0.996
#Robust Tucker-Lewis Index (TLI)                            0.995
#RMSEA                                          0.000       0.070
#Robust RMSEA                                               0.031
#SRMR                                           0.040       0.040

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#        factor2           1.014    0.003  357.893    0.000    1.014    1.014

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .653

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.891       0.894
#Tucker-Lewis Index (TLI)                       0.878       0.880
#Robust Comparative Fit Index (CFI)                         0.897
#Robust Tucker-Lewis Index (TLI)                            0.884
#RMSEA                                          0.115       0.087
#Robust RMSEA                                               0.111
#SRMR                                           0.042       0.042

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#        factor2           1.017    0.003  383.774    0.000    1.017    1.017

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .650

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.755       0.753
#Tucker-Lewis Index (TLI)                       0.726       0.724
#Robust Comparative Fit Index (CFI)                         0.760
#Robust Tucker-Lewis Index (TLI)                            0.731
#RMSEA                                          0.172       0.132
#Robust RMSEA                                               0.169
#SRMR                                           0.455       0.455

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .646

# Bifactor model
BIFmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.912       0.909
#Tucker-Lewis Index (TLI)                       0.889       0.885
#Robust Comparative Fit Index (CFI)                         0.917
#Robust Tucker-Lewis Index (TLI)                            0.894
#RMSEA                                          0.109       0.085
#Robust RMSEA                                               0.106
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .700

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    UCLALS_20        -0.251    0.094   -2.667    0.008   -0.251   -0.167
#    UCLALS_19         0.525    0.161    3.261    0.001    0.525    0.362
#    UCLALS_10        -0.155    0.050   -3.093    0.002   -0.155   -0.107
#    UCLALS_6          0.091    0.166    0.546    0.585    0.091    0.055
#    UCLALS_16        -0.098    0.067   -1.466    0.143   -0.098   -0.072
#    UCLALS_1         -0.237    0.105   -2.265    0.024   -0.237   -0.147
#    UCLALS_5         -0.261    0.082   -3.174    0.002   -0.261   -0.161
#    UCLALS_9          0.221    0.166    1.332    0.183    0.221    0.139
#    UCLALS_15         0.236    0.099    2.374    0.018    0.236    0.163
#  factor2 =~                                                            
#    UCLALS_14         0.321    0.067    4.789    0.000    0.321    0.247
#    UCLALS_11         0.403    0.080    5.032    0.000    0.403    0.248
#    UCLALS_18         0.270    0.067    4.050    0.000    0.270    0.188
#    UCLALS_17         0.108    0.063    1.733    0.083    0.108    0.078
#    UCLALS_2         -0.258    0.057   -4.523    0.000   -0.258   -0.165
#    UCLALS_13        -0.339    0.065   -5.224    0.000   -0.339   -0.212
#    UCLALS_12        -0.112    0.078   -1.446    0.148   -0.112   -0.074
#    UCLALS_7         -0.185    0.075   -2.473    0.013   -0.185   -0.131
#    UCLALS_3          0.186    0.095    1.952    0.051    0.186    0.109
#    UCLALS_8         -0.053    0.066   -0.806    0.420   -0.053   -0.036
#    UCLALS_4          0.402    0.072    5.568    0.000    0.402    0.269
#  general =~                                                            
#    UCLALS_1          1.128    0.047   24.230    0.000    1.128    0.701
#    UCLALS_2          1.382    0.049   28.171    0.000    1.382    0.886
#    UCLALS_3          1.219    0.048   25.286    0.000    1.219    0.713
#    UCLALS_4          1.211    0.051   23.757    0.000    1.211    0.810
#    UCLALS_5          1.385    0.052   26.497    0.000    1.385    0.855
#    UCLALS_6          1.351    0.046   29.297    0.000    1.351    0.825
#    UCLALS_7          1.133    0.051   22.323    0.000    1.133    0.805
#    UCLALS_8          1.026    0.048   21.306    0.000    1.026    0.691
#    UCLALS_9          1.208    0.048   25.020    0.000    1.208    0.763
#    UCLALS_10         1.280    0.050   25.448    0.000    1.280    0.879
#    UCLALS_11         1.282    0.054   23.810    0.000    1.282    0.791
#    UCLALS_12         1.328    0.053   25.023    0.000    1.328    0.877
#    UCLALS_13         1.404    0.049   28.518    0.000    1.404    0.880
#    UCLALS_14         0.973    0.056   17.472    0.000    0.973    0.749
#    UCLALS_15         1.203    0.056   21.618    0.000    1.203    0.829
#    UCLALS_16         1.094    0.055   19.955    0.000    1.094    0.808
#    UCLALS_17         1.105    0.057   19.535    0.000    1.105    0.790
#    UCLALS_18         1.180    0.053   22.094    0.000    1.180    0.823
#    UCLALS_19         1.173    0.051   23.185    0.000    1.173    0.809
#    UCLALS_20         1.306    0.051   25.855    0.000    1.306    0.865

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
$ModelLevelIndices
ECV.general            PUC  Omega.general OmegaH.general 
0.9550560      0.5210526      0.9765176      0.9754883 

$FactorLevelIndices
ECV_SS     ECV_SG    ECV_GS     Omega       OmegaH         H        FD
factor1 0.04338372 0.01983349 0.9566163 0.9517647 7.699633e-05 0.2279266 0.7243342
factor2 0.04625804 0.02511051 0.9537420 0.9564348 3.322468e-03 0.2662942 0.7465050
general 0.95505600 0.95505600 0.9550560 0.9765176 9.754883e-01 0.9769316 0.9896466






################################################################
### Panayiotou et al. 
### 

mydata <- mydata3

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 8.76; eigenvalue 2 = .69

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 10.48; eigenvalue 2 = .69

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.44, RMSEA=.106, RMSR=.06, TLI=.807

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.133, RMSR=.06, TLI=.875

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities with response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.085, RMSR=.04, TLI=.877
#     MR1  MR2
#MR1 1.00  .71
#MR2  .71 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.125, RMSR=.05, TLI=.812
#     MR1  MR2
#MR1 1.00  .73
#MR2  .73 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.988       0.940
#Tucker-Lewis Index (TLI)                       0.986       0.933
#Robust Comparative Fit Index (CFI)                         0.790
#Robust Tucker-Lewis Index (TLI)                            0.766
#Root Mean Square Error of Approximation:
#RMSEA                                          0.093       0.116
#Robust RMSEA                                               0.140
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .561

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.828       0.822
#Tucker-Lewis Index (TLI)                       0.808       0.801
#Robust Comparative Fit Index (CFI)                         0.828
#Robust Tucker-Lewis Index (TLI)                            0.808
#RMSEA                                          0.106       0.097
#Robust RMSEA                                               0.106
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .463

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.957
#Tucker-Lewis Index (TLI)                       0.990       0.951
#Robust Comparative Fit Index (CFI)                         0.845
#Robust Tucker-Lewis Index (TLI)                            0.826
#RMSEA                                          0.077       0.099
#Robust RMSEA                                               0.120
#SRMR                                           0.053       0.053

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#        factor2           0.875    0.002  398.064    0.000    0.875    0.875

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .599

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.869       0.863
#Tucker-Lewis Index (TLI)                       0.853       0.846
#Robust Comparative Fit Index (CFI)                         0.869
#Robust Tucker-Lewis Index (TLI)                            0.853
#RMSEA                                          0.093       0.085
#Robust RMSEA                                               0.093
#SRMR                                           0.051       0.051

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#        factor2           0.859    0.003  265.690    0.000    0.859    0.859

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .493

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.784       0.775
#Tucker-Lewis Index (TLI)                       0.759       0.748
#Robust Comparative Fit Index (CFI)                         0.784
#Robust Tucker-Lewis Index (TLI)                            0.759
#RMSEA                                          0.118       0.109
#Robust RMSEA                                               0.118
#SRMR                                           0.287       0.287

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .480

# Bifactor model
BIFmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.934       0.928
#Tucker-Lewis Index (TLI)                       0.916       0.909
#Robust Comparative Fit Index (CFI)                         0.934
#Robust Tucker-Lewis Index (TLI)                            0.916
#RMSEA                                          0.070       0.065
#Robust RMSEA                                               0.070
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .513

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.8064852      0.5210526      0.9451741      0.8735521 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#factor1 0.1118814 0.05687672 0.8881186 0.9135948 0.001412039 0.3931387 0.8221474
#factor2 0.2779266 0.13663806 0.7220734 0.8941383 0.225462169 0.6374030 0.8372201
#general 0.8064852 0.80648522 0.8064852 0.9451741 0.873552129 0.9444166 0.9711274






################################################################
### McDanal et al. 
### 

mydata <- mydata4

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 7.19; eigenvalue 2 = .91

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 8.65; eigenvalue 2 = .95

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.36, RMSEA=.099, RMSR=.07, TLI=.773

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.125, RMSR=.07, TLI=.743

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities with response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.082, RMSR=.05, TLI=.846
#     MR1  MR2
#MR1 1.00  .73
#MR2  .73 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.108, RMSR=.05, TLI=.809
#     MR1  MR2
#MR1 1.00  .76
#MR2  .76 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.972       0.887
#Tucker-Lewis Index (TLI)                       0.969       0.874
#Robust Comparative Fit Index (CFI)                         0.755
#Robust Tucker-Lewis Index (TLI)                            0.727
#RMSEA                                          0.097       0.115
#Robust RMSEA                                               0.129
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .442

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.797       0.795
#Tucker-Lewis Index (TLI)                       0.773       0.770
#Robust Comparative Fit Index (CFI)                         0.798
#Robust Tucker-Lewis Index (TLI)                            0.774
#RMSEA                                          0.100       0.092
#Robust RMSEA                                               0.099
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .369

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 3 factors based on Whisman et al
EGAmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.985       0.935
#Tucker-Lewis Index (TLI)                       0.983       0.927
#Robust Comparative Fit Index (CFI)                         0.844
#Robust Tucker-Lewis Index (TLI)                            0.824
#RMSEA                                          0.071       0.088
#Robust RMSEA                                               0.104
#SRMR                                           0.058       0.058

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#       factor2           0.783    0.010   76.696    0.000    0.783    0.783

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .494

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.875       0.873
#Tucker-Lewis Index (TLI)                       0.860       0.857
#Robust Comparative Fit Index (CFI)                         0.876
#Robust Tucker-Lewis Index (TLI)                            0.861
#RMSEA                                          0.078       0.072
#Robust RMSEA                                               0.078
#SRMR                                           0.049       0.049

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
# factor1 ~~                                                            
#       factor2           0.767    0.012   61.547    0.000    0.767    0.767

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .417

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.801       0.796
#Tucker-Lewis Index (TLI)                       0.778       0.771
#Robust Comparative Fit Index (CFI)                         0.802
#Robust Tucker-Lewis Index (TLI)                            0.779
#RMSEA                                          0.099       0.092
#Robust RMSEA                                               0.098
#SRMR                                           0.222       0.222

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .418

# Bifactor model
BIFmodel= '
 factor1 =~  UCLALS_20+UCLALS_19+UCLALS_10+UCLALS_6+UCLALS_16+UCLALS_1+UCLALS_5+
             UCLALS_9+UCLALS_15
 factor2 =~  UCLALS_14+UCLALS_11+UCLALS_18+UCLALS_17+UCLALS_2+UCLALS_13+
             UCLALS_12+UCLALS_7+UCLALS_3+UCLALS_8+UCLALS_4
 general=~  UCLALS_1+UCLALS_2+UCLALS_3+UCLALS_4+UCLALS_5+UCLALS_6+UCLALS_7+UCLALS_8+
            UCLALS_9+UCLALS_10+UCLALS_11+UCLALS_12+UCLALS_13+UCLALS_14+UCLALS_15+
            UCLALS_16+UCLALS_17+UCLALS_18+UCLALS_19+UCLALS_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.929       0.927
#Tucker-Lewis Index (TLI)                       0.910       0.908
#Robust Comparative Fit Index (CFI)                         0.930
#Robust Tucker-Lewis Index (TLI)                            0.911
#RMSEA                                          0.063       0.058
#Robust RMSEA                                               0.062
#SRMR                                           0.040       0.040

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .434

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7177554      0.5210526      0.9261203      0.7907946 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.1283657 0.06398264 0.8716343 0.8805257 0.01067098 0.3842118 0.7813130
#factor2 0.4351665 0.21826197 0.5648335 0.8762985 0.37727884 0.7098314 0.8484413
#general 0.7177554 0.71775538 0.7177554 0.9261203 0.79079459 0.9210905 0.9589432





