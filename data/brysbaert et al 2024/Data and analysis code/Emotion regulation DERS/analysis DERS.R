######################################################################
### Analysis Difficulties in Emotion Regulation Scale (DERS)
###

################################################################
### Valencia et al. (2024)  https://osf.io/t6adk/
### 

library(haven)
DERS_Valencia <- read_sav("DERS_Valencia.sav")
colnames(DERS_Valencia)
mydata <- DERS_Valencia[,c(8:43)]
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 13.67; eigenvalue 2 = 3.78

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 15.69; eigenvalue 2 = 4.15

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.38, RMSEA=.124, RMSR=.13, TLI=.604

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.155, RMSR=.14, TLI=.548

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 5 communities

# Give solution with 6 factors (theory-based)
fit4 <- fa(mydata,6)
fit4
diagram(fit4)
# %variance explained=.62, RMSEA=.057, RMSR=.02, TLI=.915
#       MR1  MR4  MR3  MR6  MR2  MR5
#MR1 1.00 0.63 0.60 0.59 0.12 0.12
#MR4 0.63 1.00 0.52 0.44 0.06 0.10
#MR3 0.60 0.52 1.00 0.52 0.08 0.19
#MR6 0.59 0.44 0.52 1.00 0.30 0.24
#MR2 0.12 0.06 0.08 0.30 1.00 0.46
#MR5 0.12 0.10 0.19 0.24 0.46 1.00
#
# Average absolute correlation = .331

fit4 <- fa(rho,6,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.097, RMSR=.03, TLI=.822
#     MR1  MR3  MR4  MR6  MR2  MR5
#MR1 1.00 0.61 0.63 0.59 0.15 0.13
#MR3 0.61 1.00 0.53 0.53 0.10 0.22
#MR4 0.63 0.53 1.00 0.44 0.07 0.10
#MR6 0.59 0.53 0.44 1.00 0.32 0.26
#MR2 0.15 0.10 0.07 0.32 1.00 0.47
#MR5 0.13 0.22 0.10 0.26 0.47 1.00
#
# Average absolute correlation = .343

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  ders1+ders2+ders3+ders4+ders5+ders6+
            ders7+ders8+ders9+ders10+ders11+ders12+
            ders13+ders14+ders15+ders16+ders17+ders18+
            ders19+ders20+ders21+ders22+ders23+ders24+
            ders25+ders26+ders27+ders28+ders29+ders30+
            ders31+ders32+ders33+ders34+ders35+ders36
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.929       0.760
#Tucker-Lewis Index (TLI)                       0.924       0.746
#Robust Comparative Fit Index (CFI)                         0.578
#Robust Tucker-Lewis Index (TLI)                            0.552
#RMSEA                                          0.214       0.170
#Robust RMSEA                                               0.160
#SRMR                                           0.144       0.144

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .568

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.626       0.630
#Tucker-Lewis Index (TLI)                       0.603       0.608
#Robust Comparative Fit Index (CFI)                         0.638
#Robust Tucker-Lewis Index (TLI)                            0.616
#RMSEA                                          0.127       0.109
#Robust RMSEA                                               0.124
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .448

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 6 factors based on the scale construction
EGAmodel= '
 nonaccept     =~ ders11+ders12+ders21+ders23+ders25+ders29
 goal_behav    =~ ders13+ders18+ders20+ders26+ders33
 impulse_contr =~ ders3+ders14+ders19+ders24+ders27+ders32
 emo_awareness =~ ders2+ders6+ders8+ders10+ders17+ders34
 emo_regul     =~ ders15+ders16+ders22+ders28+ders30+ders31+ders35+ders36
 emo_clarity   =~ ders1+ders4+ders5+ders7+ders9
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.970       0.882
#Tucker-Lewis Index (TLI)                       0.968       0.872
#Robust Comparative Fit Index (CFI)                         0.773
#Robust Tucker-Lewis Index (TLI)                            0.753
#RMSEA                                          0.140       0.120
#Robust RMSEA                                               0.119
#SRMR                                           0.103       0.103

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .708

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  nonaccept ~~                                                          
#    goal_behav        0.757    0.025   29.849    0.000    0.757    0.757
#    impulse_contr     0.705    0.028   25.016    0.000    0.705    0.705
#    emo_awareness     0.240    0.042    5.655    0.000    0.240    0.240
#    emo_regul         0.922    0.011   87.409    0.000    0.922    0.922
#    emo_clarity       0.676    0.030   22.302    0.000    0.676    0.676
#  goal_behav ~~                                                         
#    impulse_contr     0.684    0.031   22.275    0.000    0.684    0.684
#    emo_awareness     0.205    0.044    4.644    0.000    0.205    0.205
#    emo_regul         0.853    0.017   50.460    0.000    0.853    0.853
#    emo_clarity       0.573    0.036   15.755    0.000    0.573    0.573
#  impulse_contr ~~                                                      
#    emo_awareness     0.259    0.049    5.281    0.000    0.259    0.259
#    emo_regul         0.833    0.019   44.144    0.000    0.833    0.833
#    emo_clarity       0.634    0.034   18.727    0.000    0.634    0.634
#  emo_awareness ~~                                                      
#    emo_regul         0.371    0.041    9.027    0.000    0.371    0.371
#    emo_clarity       0.718    0.025   28.686    0.000    0.718    0.718
#  emo_regul ~~                                                          
#    emo_clarity       0.760    0.026   29.611    0.000    0.760    0.760
#
# average correlation = .613


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.832       0.839
#Tucker-Lewis Index (TLI)                       0.817       0.825
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.832
#RMSEA                                          0.086       0.073
#Robust RMSEA                                               0.082
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .636

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  nonaccept ~~                                                          
#    goal_behav        0.731    0.038   19.468    0.000    0.731    0.731
#    impulse_contr     0.673    0.044   15.451    0.000    0.673    0.673
#    emo_awareness     0.186    0.063    2.959    0.003    0.186    0.186
#    emo_regul         0.903    0.020   44.821    0.000    0.903    0.903
#    emo_clarity       0.658    0.044   14.950    0.000    0.658    0.658
#  goal_behav ~~                                                         
#    impulse_contr     0.633    0.045   13.919    0.000    0.633    0.633
#    emo_awareness     0.130    0.070    1.854    0.064    0.130    0.130
#    emo_regul         0.823    0.027   30.930    0.000    0.823    0.823
#    emo_clarity       0.547    0.050   11.037    0.000    0.547    0.547
#  impulse_contr ~~                                                      
#    emo_awareness     0.165    0.064    2.575    0.010    0.165    0.165
#    emo_regul         0.788    0.034   22.873    0.000    0.788    0.788
#    emo_clarity       0.606    0.050   12.223    0.000    0.606    0.606
#  emo_awareness ~~                                                      
#    emo_regul         0.263    0.066    3.978    0.000    0.263    0.263
#    emo_clarity       0.451    0.074    6.074    0.000    0.451    0.451
#  emo_regul ~~                                                          
#    emo_clarity       0.745    0.038   19.680    0.000    0.745    0.745
#
# average correlation = .553

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.695       0.695
#Tucker-Lewis Index (TLI)                       0.676       0.676
#Robust Comparative Fit Index (CFI)                         0.706
#Robust Tucker-Lewis Index (TLI)                            0.688
#RMSEA                                          0.115       0.099
#Robust RMSEA                                               0.111
#SRMR                                           0.330       0.330

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .613


# Bifactor model
BIFmodel= '
 nonaccept     =~ ders11+ders12+ders21+ders23+ders25+ders29
 goal_behav    =~ ders13+ders18+ders20+ders26+ders33
 impulse_contr =~ ders3+ders14+ders19+ders24+ders27+ders32
 emo_awareness =~ ders2+ders6+ders8+ders10+ders17+ders34
 emo_regul     =~ ders15+ders16+ders22+ders28+ders30+ders31+ders35+ders36
 emo_clarity   =~ ders1+ders4+ders5+ders7+ders9
 general=~  ders1+ders2+ders3+ders4+ders5+ders6+
            ders7+ders8+ders9+ders10+ders11+ders12+
            ders13+ders14+ders15+ders16+ders17+ders18+
            ders19+ders20+ders21+ders22+ders23+ders24+
            ders25+ders26+ders27+ders28+ders29+ders30+
            ders31+ders32+ders33+ders34+ders35+ders36
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.842       0.848
#Tucker-Lewis Index (TLI)                       0.822       0.828
#Robust Comparative Fit Index (CFI)                         0.855
#Robust Tucker-Lewis Index (TLI)                            0.836
#RMSEA                                          0.085       0.072
#Robust RMSEA                                               0.081
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .640

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6481753      0.8523810      0.9639922      0.8816672 
#
#$FactorLevelIndices
#                  ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#nonaccept     0.21259169 0.04222802 0.7874083 0.9216899 0.18281646 0.5112170 0.8115870
#goal_behav    0.34657491 0.05114839 0.6534251 0.8652989 0.29303748 0.5862858 0.8637523
#impulse_contr 0.37409478 0.06105346 0.6259052 0.8673779 0.30906965 0.6359336 0.8768051
#emo_awareness 0.93525660 0.12778645 0.0647434 0.8285942 0.78525246 0.8547428 0.9295496
#emo_regul     0.05541476 0.01294448 0.9445852 0.9115039 0.01384748 0.2283140 0.6841153
#emo_clarity   0.47080708 0.05666390 0.5291929 0.8108117 0.40064139 0.6121674 0.8586685
#general       0.64817530 0.64817530 0.6481753 0.9639922 0.88166718 0.9669274 0.9759064






