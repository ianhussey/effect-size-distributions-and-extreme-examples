######################################################################
### Analysis Domain-Specific Risk-Taking (DOSPERT) Scale
###

################################################################
### Frey et al. (2023) https://link.springer.com/article/10.1007/s11166-022-09398-5
### data https://osf.io/pjt57/

DOSPERT_Frey <- read.csv("DOSPERT_Frey.csv")
colnames(DOSPERT_Frey)
mydata <- DOSPERT_Frey[,c(3:8,21:26,39:44,57:62,75:80)]
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 5.69; eigenvalue 2 = 2.30

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 8 factors and 5 components
# Eigenvalue 1 = 7.37; eigenvalue 2 = 2.88

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.12, RMSR=.11, TLI=.388

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.147, RMSR=.13, TLI=.385

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 5 communities

# Give solution with 5 factors (theory-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.06, RMSR=.03, TLI=.847
#     MR2  MR3  MR1  MR7  MR4  MR6  MR8  MR5  MR9
#MR1  MR4   MR3   MR2  MR5
#MR1  1.00 0.34  0.45 -0.05 0.20
#MR4  0.34 1.00  0.32  0.12 0.36
#MR3  0.45 0.32  1.00 -0.05 0.41
#MR2 -0.05 0.12 -0.05  1.00 0.17
#MR5  0.20 0.36  0.41  0.17 1.00
#
# Average absolute correlation = .247

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.08, RMSR=.04, TLI=.816
#      MR3  MR4   MR1   MR2  MR5
#MR3  1.00 0.37  0.52 -0.09 0.09
#MR4  0.37 1.00  0.41  0.13 0.29
#MR1  0.52 0.41  1.00 -0.07 0.37
#MR2 -0.09 0.13 -0.07  1.00 0.15
#MR5  0.09 0.29  0.37  0.15 1.00
#
# Average absolute correlation = .249

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  ethR_1+ethR_2+ethR_3+ethR_4+ethR_5+ethR_6+
            finR_1+finR_2+finR_3+finR_4+finR_5+finR_6+
            heaR_1+heaR_2+heaR_3+heaR_4+heaR_5+heaR_6+
            recR_1+recR_2+recR_3+recR_4+recR_5+recR_6+
            socR_1+socR_2+socR_3+socR_4+socR_5+socR_6
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.785       0.588
#Tucker-Lewis Index (TLI)                       0.769       0.557
#Robust Comparative Fit Index (CFI)                         0.371
#Robust Tucker-Lewis Index (TLI)                            0.325
#RMSEA                                          0.176       0.160
#Robust RMSEA                                               0.155
#SRMR                                           0.144       0.144

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .261

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.431       0.422
#Tucker-Lewis Index (TLI)                       0.389       0.379
#Robust Comparative Fit Index (CFI)                         0.432
#Robust Tucker-Lewis Index (TLI)                            0.390
#RMSEA                                          0.120       0.110
#Robust RMSEA                                               0.120
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .208

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 59 factors based on the scale construction
EGAmodel= '
 ethical   =~ ethR_1+ethR_2+ethR_3+ethR_4+ethR_5+ethR_6
 financial =~ finR_1+finR_2+finR_3+finR_4+finR_5+finR_6
 health    =~ heaR_1+heaR_2+heaR_3+heaR_4+heaR_5+heaR_6
 recrea    =~ recR_1+recR_2+recR_3+recR_4+recR_5+recR_6
 social    =~ socR_1+socR_2+socR_3+socR_4+socR_5+socR_6
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.949       0.889
#Tucker-Lewis Index (TLI)                       0.944       0.878
#Robust Comparative Fit Index (CFI)                         0.772
#Robust Tucker-Lewis Index (TLI)                            0.749
#RMSEA                                          0.087       0.084
#Robust RMSEA                                               0.094
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .467

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  ethical ~~                                                            
#    financial         0.480    0.017   28.877    0.000    0.480    0.480
#    health            0.719    0.015   48.628    0.000    0.719    0.719
#    recrea            0.384    0.018   21.164    0.000    0.384    0.384
#    social           -0.145    0.020   -7.111    0.000   -0.145   -0.145
#  financial ~~                                                          
#    health            0.484    0.017   28.818    0.000    0.484    0.484
#    recrea            0.475    0.016   29.716    0.000    0.475    0.475
#    social            0.040    0.019    2.083    0.037    0.040    0.040
#  health ~~                                                             
#    recrea            0.551    0.015   36.268    0.000    0.551    0.551
#    social            0.096    0.020    4.679    0.000    0.096    0.096
#  recrea ~~                                                             
#    social            0.150    0.019    8.110    0.000    0.150    0.150
#
# average correlation = .352


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.803       0.796
#Tucker-Lewis Index (TLI)                       0.783       0.775
#Robust Comparative Fit Index (CFI)                         0.804
#Robust Tucker-Lewis Index (TLI)                            0.784
#RMSEA                                          0.072       0.066
#Robust RMSEA                                               0.071
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .327

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  ethical ~~                                                            
#    financial         0.469    0.023   20.295    0.000    0.469    0.469
#    health            0.685    0.022   31.847    0.000    0.685    0.685
#    recrea            0.323    0.022   14.564    0.000    0.323    0.323
#    social           -0.113    0.026   -4.282    0.000   -0.113   -0.113
#  financial ~~                                                          
#    health            0.472    0.022   21.538    0.000    0.472    0.472
#    recrea            0.347    0.022   15.842    0.000    0.347    0.347
#    social           -0.011    0.023   -0.496    0.620   -0.011   -0.011
#  health ~~                                                             
#    recrea            0.483    0.021   22.566    0.000    0.483    0.483
#    social            0.134    0.026    5.045    0.000    0.134    0.134
#  recrea ~~                                                             
#    social            0.137    0.022    6.301    0.000    0.137    0.137
#
# average correlation = .317

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.729       0.719
#Tucker-Lewis Index (TLI)                       0.709       0.698
#Robust Comparative Fit Index (CFI)                         0.730
#Robust Tucker-Lewis Index (TLI)                            0.710
#RMSEA                                          0.083       0.077
#Robust RMSEA                                               0.083
#SRMR                                           0.149       0.149

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .334


# Bifactor model
BIFmodel= '
 ethical   =~ ethR_1+ethR_2+ethR_3+ethR_4+ethR_5+ethR_6
 financial =~ finR_1+finR_2+finR_3+finR_4+finR_5+finR_6
 health    =~ heaR_1+heaR_2+heaR_3+heaR_4+heaR_5+heaR_6
 recrea    =~ recR_1+recR_2+recR_3+recR_4+recR_5+recR_6
 social    =~ socR_1+socR_2+socR_3+socR_4+socR_5+socR_6
 general=~  ethR_1+ethR_2+ethR_3+ethR_4+ethR_5+ethR_6+
            finR_1+finR_2+finR_3+finR_4+finR_5+finR_6+
            heaR_1+heaR_2+heaR_3+heaR_4+heaR_5+heaR_6+
            recR_1+recR_2+recR_3+recR_4+recR_5+recR_6+
            socR_1+socR_2+socR_3+socR_4+socR_5+socR_6
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.859       0.852
#Tucker-Lewis Index (TLI)                       0.836       0.828
#Robust Comparative Fit Index (CFI)                         0.860
#Robust Tucker-Lewis Index (TLI)                            0.837
#RMSEA                                          0.062       0.058
#Robust RMSEA                                               0.062
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .373

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3543157      0.8275862      0.8855617      0.6116455 
#
#$FactorLevelIndices
#             ECV_SS    ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#ethical   0.7536254 0.1259793 0.24637459 0.7627817 0.577693424 0.6991581 0.8461791
#financial 0.2712249 0.0815887 0.72877512 0.8827410 0.008101717 0.5634638 0.8574866
#health    0.7272392 0.1040386 0.27276081 0.7197323 0.526078161 0.6351845 0.8046558
#recrea    0.7645200 0.1674309 0.23548004 0.8263232 0.622228741 0.8390488 0.9308969
#social    0.9805119 0.1666469 0.01948809 0.7642647 0.761614808 0.7802887 0.8846209
#general   0.3543157 0.3543157 0.35431570 0.8855617 0.611645548 0.8916107 0.9528686





