################################################################
### Analysis The Questionnaire for Eudaimonic Well-Being
### Ishii et al. (2022) https://bmcpsychology.biomedcentral.com/articles/10.1186/s40359-021-00707-2 
### data https://osf.io/62daf/  
###

library(readxl)
QEWB_Ishii <- read_excel("QEWB_Ishii.xlsx")
colnames(QEWB_Ishii)
mydata  <- as.data.frame(QEWB_Ishii[,5:25])
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 4.31; eigenvalue 2 = 1.15

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components 
# Eigenvalue 1 = 5.01; eigenvalue 2 = 1.37

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.21, RMSEA=.086, RMSR=.08, TLI=.669

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.24, RMSEA=.108, RMSR=.10, TLI=.611

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.045, RMSR=.03, TLI=.909
#     MR1  MR3  MR2
#MR1 1.00 0.42 0.16
#MR3 0.42 1.00 0.05
#MR2 0.16 0.05 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.06, RMSR=.04, TLI=.879
#     MR1  MR3  MR2
#MR1 1.00 0.42 0.17
#MR3 0.42 1.00 0.10
#MR2 0.17 0.10 1.00

# Single factor model lavaan
UNImodel= '
 general=~ QEWB1+QEWB2+QEWB3rev+QEWB4+QEWB5+QEWB6+QEWB7rev+QEWB8+QEWB9+QEWB10+QEWB11rev+QEWB12rev+QEWB13+QEWB14+QEWB15+
           QEWB16rev+QEWB17+QEWB18+QEWB19rev+QEWB20rev+QEWB21
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.908       0.808
#Tucker-Lewis Index (TLI)                       0.898       0.787
#Robust Comparative Fit Index (CFI)                         0.661
#Robust Tucker-Lewis Index (TLI)                            0.623
#RMSEA                                          0.098       0.104
#Robust RMSEA                                               0.107
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .17

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.707       0.708
#Tucker-Lewis Index (TLI)                       0.675       0.675
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.679
#RMSEA                                          0.086       0.078
#Robust RMSEA                                               0.085
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .114

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on EGA analysis
EGAmodel= '
 factor1=~ QEWB8+QEWB4+QEWB1+QEWB15+QEWB18+QEWB13+QEWB14+QEWB17+QEWB10+QEWB5
 factor2=~ QEWB9+QEWB21+QEWB2+QEWB6+QEWB11rev+QEWB16rev
 factor3=~ QEWB20rev+QEWB3rev+QEWB7rev+QEWB12rev+QEWB19rev
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.948       0.890
#Tucker-Lewis Index (TLI)                       0.941       0.876
#Robust Comparative Fit Index (CFI)                         0.796
#Robust Tucker-Lewis Index (TLI)                            0.770
#RMSEA                                          0.075       0.079
#Robust RMSEA                                               0.084
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .26

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.732    0.020   36.285    0.000    0.732    0.732
#  factor3           0.356    0.036    9.902    0.000    0.356    0.356
#factor2 ~~                                                            
#  factor3           0.373    0.030   12.467    0.000    0.373    0.373


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.835       0.836
#Tucker-Lewis Index (TLI)                       0.814       0.815
#Robust Comparative Fit Index (CFI)                         0.840
#Robust Tucker-Lewis Index (TLI)                            0.819
#RMSEA                                          0.065       0.059
#Robust RMSEA                                               0.064
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .305

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.716    0.030   24.028    0.000    0.716    0.716
#    factor3           0.274    0.054    5.035    0.000    0.274    0.274
#  factor2 ~~                                                            
#    factor3           0.253    0.051    4.971    0.000    0.253    0.253



CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.749       0.749
#Tucker-Lewis Index (TLI)                       0.721       0.721
#Robust Comparative Fit Index (CFI)                         0.753
#Robust Tucker-Lewis Index (TLI)                            0.726
#RMSEA                                          0.079       0.072
#Robust RMSEA                                               0.078
#SRMR                                           0.136       0.136

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .277

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1=~ QEWB8+QEWB4+QEWB1+QEWB15+QEWB18+QEWB13+QEWB14+QEWB17+QEWB10+QEWB5
 factor2=~ QEWB9+QEWB21+QEWB2+QEWB6+QEWB11rev+QEWB16rev
 factor3=~ QEWB20rev+QEWB3rev+QEWB7rev+QEWB12rev+QEWB19rev
 general=~ QEWB8+QEWB4+QEWB1+QEWB15+QEWB18+QEWB13+QEWB14+QEWB17+QEWB10+QEWB5+
           QEWB9+QEWB21+QEWB2+QEWB6+QEWB11rev+QEWB16rev+
           QEWB20rev+QEWB3rev+QEWB7rev+QEWB12rev+QEWB19rev
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.897       0.893
#Tucker-Lewis Index (TLI)                       0.871       0.866
#Robust Comparative Fit Index (CFI)                         0.899
#Robust Tucker-Lewis Index (TLI)                            0.874
#RMSEA                                          0.054       0.050
#Robust RMSEA                                               0.053
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .305

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5773209      0.6666667      0.8416109      0.6708112 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4100839 0.1334558 0.5899161 0.7054497 0.2651607 0.5290393 0.7160665
#factor2 0.2276316 0.1054980 0.7723684 0.8554187 0.1398968 0.4688723 0.6990642
#factor3 0.8703015 0.1837254 0.1296985 0.6226970 0.5871586 0.6528384 0.8119629
#general 0.5773209 0.5773209 0.5773209 0.8416109 0.6708112 0.8573061 0.9027410




