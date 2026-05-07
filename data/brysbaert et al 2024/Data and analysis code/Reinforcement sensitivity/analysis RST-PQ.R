################################################################
### Analysis RST-PQ Questionnaire
### Sözer et al Turkey https://dusunenadamdergisi.org/storage/upload/pdfs/1655381096-en.pdf

library(haven)
RST_Sözer <- read_sav("RST_Sözer.sav")
colnames(RST_Sözer)
mydata <- as.data.frame(RST_Sözer[,3:67])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) #4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 7 components
# Eigenvalue 1 = 10.5; eigenvalue 2 = 5.94

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 components
# Eigenvalue 1 = 12.25; eigenvalue 2 = 7.11

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.081, RMSR=.12, TLI=.392

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.112, RMSR=.14, TLI=.293

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities with response bias

# Give solution with 6 factors
fit4 <- fa(mydata,6)
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.042, RMSR=.04, TLI=.834
#      MR1   MR3   MR2   MR4   MR5   MR6
#MR1  1.00 -0.24 -0.07  0.36  0.01  0.20
#MR3 -0.24  1.00  0.23 -0.11  0.15 -0.08
#MR2 -0.07  0.23  1.00 -0.02  0.28  0.01
#MR4  0.36 -0.11 -0.02  1.00  0.14 -0.01
#MR5  0.01  0.15  0.28  0.14  1.00 -0.06
#MR6  0.20 -0.08  0.01 -0.01 -0.06  1.00

fit4 <- fa(rho,6,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.08, RMSR=.05, TLI=.643
#      MR1   MR3   MR2   MR4   MR5   MR6
#MR1  1.00 -0.23 -0.06  0.36  0.03  0.18
#MR3 -0.23  1.00  0.22 -0.11  0.17 -0.08
#MR2 -0.06  0.22  1.00 -0.03  0.28 -0.01
#MR4  0.36 -0.11 -0.03  1.00  0.13 -0.02
#MR5  0.03  0.17  0.28  0.13  1.00 -0.08
#MR6  0.18 -0.08 -0.01 -0.02 -0.08  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ RST_PQ1+RST_PQ2+RST_PQ3+RST_PQ4+RST_PQ5+RST_PQ6+RST_PQ7+RST_PQ8+RST_PQ9+RST_PQ10+RST_PQ11+
           RST_PQ12+RST_PQ13+RST_PQ14+RST_PQ15+RST_PQ16+RST_PQ17+RST_PQ18+RST_PQ19+RST_PQ20+RST_PQ21+
           RST_PQ22+RST_PQ23+RST_PQ24+RST_PQ25+RST_PQ26+RST_PQ27+RST_PQ28+RST_PQ29+RST_PQ30+RST_PQ31+
           RST_PQ32+RST_PQ33+RST_PQ34+RST_PQ35+RST_PQ36+RST_PQ37+RST_PQ38+RST_PQ39+RST_PQ40+RST_PQ41+
           RST_PQ42+RST_PQ43+RST_PQ44+RST_PQ45+RST_PQ46+RST_PQ47+RST_PQ48+RST_PQ49+RST_PQ50+RST_PQ51+
           RST_PQ52+RST_PQ53+RST_PQ54+RST_PQ55+RST_PQ56+RST_PQ57+RST_PQ58+RST_PQ59+RST_PQ60+RST_PQ61+
           RST_PQ62+RST_PQ63+RST_PQ64+RST_PQ65
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.720       0.539
#Tucker-Lewis Index (TLI)                       0.711       0.524
#Robust Comparative Fit Index (CFI)                         0.336
#Robust Tucker-Lewis Index (TLI)                            0.314
#RMSEA                                          0.142       0.091
#Robust RMSEA                                               0.113
#SRMR                                           0.145       0.145

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .138

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.408       0.408
#Tucker-Lewis Index (TLI)                       0.389       0.389
#Robust Comparative Fit Index (CFI)                         0.411
#Robust Tucker-Lewis Index (TLI)                            0.392
#RMSEA                                          0.084       0.082
#Robust RMSEA                                               0.084
#SRMR                                           0.118       0.118

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .085

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 6 factors based on theoretical analysis
EGAmodel= '
 FFFs =~ RST_PQ9+RST_PQ19+RST_PQ39+RST_PQ45+RST_PQ46+RST_PQ48+RST_PQ52+RST_PQ58+RST_PQ59+RST_PQ62
 BIS  =~ RST_PQ1+RST_PQ2+RST_PQ6+RST_PQ7+RST_PQ10+RST_PQ17+RST_PQ18+RST_PQ21+RST_PQ29+RST_PQ33+RST_PQ34+
         RST_PQ42+RST_PQ43+RST_PQ47+RST_PQ49+RST_PQ50+RST_PQ55+RST_PQ56+RST_PQ57+RST_PQ60+RST_PQ61+
         RST_PQ63+RST_PQ64
 RI   =~ RST_PQ11+RST_PQ13+RST_PQ14+RST_PQ15+RST_PQ26+RST_PQ32+RST_PQ35
 GDP  =~ RST_PQ5+RST_PQ12+RST_PQ20+RST_PQ31+RST_PQ41+RST_PQ54+RST_PQ65
 RR   =~ RST_PQ3+RST_PQ4+RST_PQ8+RST_PQ16+RST_PQ23+RST_PQ24+RST_PQ25+RST_PQ30+RST_PQ36+RST_PQ37
 Imp  =~ RST_PQ22+RST_PQ27+RST_PQ28+RST_PQ38+RST_PQ40+RST_PQ44+RST_PQ51+RST_PQ53
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.929       0.867
#Tucker-Lewis Index (TLI)                       0.927       0.862
#Robust Comparative Fit Index (CFI)                         0.672
#Robust Tucker-Lewis Index (TLI)                            0.659
#RMSEA                                          0.072       0.049
#Robust RMSEA                                               0.080
#SRMR                                           0.088       0.088

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  FFFs ~~                                                               
#    BIS               0.511    0.040   12.909    0.000    0.511    0.511
#    RI               -0.184    0.055   -3.330    0.001   -0.184   -0.184
#    GDP              -0.118    0.054   -2.172    0.030   -0.118   -0.118
#    RR                0.072    0.059    1.216    0.224    0.072    0.072
#    Imp               0.264    0.059    4.466    0.000    0.264    0.264
#  BIS ~~                                                                
#    RI               -0.231    0.048   -4.833    0.000   -0.231   -0.231
#    GDP              -0.330    0.046   -7.245    0.000   -0.330   -0.330
#    RR               -0.108    0.053   -2.051    0.040   -0.108   -0.108
#    Imp               0.318    0.054    5.905    0.000    0.318    0.318
#  RI ~~                                                                 
#    GDP               0.575    0.037   15.556    0.000    0.575    0.575
#    RR                0.624    0.038   16.522    0.000    0.624    0.624
#    Imp               0.432    0.052    8.325    0.000    0.432    0.432
#  GDP ~~                                                                
#    RR                0.425    0.047    9.095    0.000    0.425    0.425
#    Imp               0.038    0.056    0.680    0.497    0.038    0.038
#  RR ~~                                                                 
#    Imp               0.581    0.048   12.128    0.000    0.581    0.581

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .404

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.756       0.761
#Tucker-Lewis Index (TLI)                       0.746       0.751
#Robust Comparative Fit Index (CFI)                         0.763
#Robust Tucker-Lewis Index (TLI)                            0.753
#RMSEA                                          0.055       0.052
#Robust RMSEA                                               0.053
##SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .335

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  FFFs ~~                                                               
#    BIS               0.499    0.048   10.378    0.000    0.499    0.499
#    RI               -0.177    0.062   -2.827    0.005   -0.177   -0.177
#    GDP              -0.116    0.060   -1.933    0.053   -0.116   -0.116
#    RR                0.069    0.069    1.005    0.315    0.069    0.069
#    Imp               0.223    0.074    3.013    0.003    0.223    0.223
#  BIS ~~                                                                
#    RI               -0.235    0.055   -4.253    0.000   -0.235   -0.235
#    GDP              -0.330    0.051   -6.413    0.000   -0.330   -0.330
#    RR               -0.116    0.067   -1.727    0.084   -0.116   -0.116
#    Imp               0.277    0.067    4.116    0.000    0.277    0.277
#  RI ~~                                                                 
#    GDP               0.549    0.046   12.013    0.000    0.549    0.549
#    RR                0.608    0.049   12.322    0.000    0.608    0.608
#    Imp               0.428    0.064    6.693    0.000    0.428    0.428
#  GDP ~~                                                                
#    RR                0.405    0.057    7.042    0.000    0.405    0.405
#    Imp              -0.005    0.072   -0.076    0.940   -0.005   -0.005
#  RR ~~                                                                 
#    Imp               0.551    0.066    8.416    0.000    0.551    0.551


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.706       0.711
#Tucker-Lewis Index (TLI)                       0.696       0.701
#Robust Comparative Fit Index (CFI)                         0.712
#Robust Tucker-Lewis Index (TLI)                            0.703
#RMSEA                                          0.060       0.057
#Robust RMSEA                                               0.059
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .322

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 FFFs =~ RST_PQ9+RST_PQ19+RST_PQ39+RST_PQ45+RST_PQ46+RST_PQ48+RST_PQ52+RST_PQ58+RST_PQ59+RST_PQ62
 BIS  =~ RST_PQ1+RST_PQ2+RST_PQ6+RST_PQ7+RST_PQ10+RST_PQ17+RST_PQ18+RST_PQ21+RST_PQ29+RST_PQ33+RST_PQ34+
         RST_PQ42+RST_PQ43+RST_PQ47+RST_PQ49+RST_PQ50+RST_PQ55+RST_PQ56+RST_PQ57+RST_PQ60+RST_PQ61+
         RST_PQ63+RST_PQ64
 RI   =~ RST_PQ11+RST_PQ13+RST_PQ14+RST_PQ15+RST_PQ26+RST_PQ32+RST_PQ35
 GDP  =~ RST_PQ5+RST_PQ12+RST_PQ20+RST_PQ31+RST_PQ41+RST_PQ54+RST_PQ65
 RR   =~ RST_PQ3+RST_PQ4+RST_PQ8+RST_PQ16+RST_PQ23+RST_PQ24+RST_PQ25+RST_PQ30+RST_PQ36+RST_PQ37
 Imp  =~ RST_PQ22+RST_PQ27+RST_PQ28+RST_PQ38+RST_PQ40+RST_PQ44+RST_PQ51+RST_PQ53
 general=~ RST_PQ1+RST_PQ2+RST_PQ3+RST_PQ4+RST_PQ5+RST_PQ6+RST_PQ7+RST_PQ8+RST_PQ9+RST_PQ10+RST_PQ11+
           RST_PQ12+RST_PQ13+RST_PQ14+RST_PQ15+RST_PQ16+RST_PQ17+RST_PQ18+RST_PQ19+RST_PQ20+RST_PQ21+
           RST_PQ22+RST_PQ23+RST_PQ24+RST_PQ25+RST_PQ26+RST_PQ27+RST_PQ28+RST_PQ29+RST_PQ30+RST_PQ31+
           RST_PQ32+RST_PQ33+RST_PQ34+RST_PQ35+RST_PQ36+RST_PQ37+RST_PQ38+RST_PQ39+RST_PQ40+RST_PQ41+
           RST_PQ42+RST_PQ43+RST_PQ44+RST_PQ45+RST_PQ46+RST_PQ47+RST_PQ48+RST_PQ49+RST_PQ50+RST_PQ51+
           RST_PQ52+RST_PQ53+RST_PQ54+RST_PQ55+RST_PQ56+RST_PQ57+RST_PQ58+RST_PQ59+RST_PQ60+RST_PQ61+
           RST_PQ62+RST_PQ63+RST_PQ64+RST_PQ65
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.765       0.766
#Tucker-Lewis Index (TLI)                       0.749       0.750
#Robust Comparative Fit Index (CFI)                         0.770
#Robust Tucker-Lewis Index (TLI)                            0.754
#RMSEA                                          0.054       0.052
#Robust RMSEA                                               0.053
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .336

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4034948      0.8014423      0.8750954      0.5143495 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#FFFs    0.6933021 0.08690116 0.30669788 0.7989295 0.5670337 0.7221188 0.8535389
#BIS     0.1993458 0.07885244 0.80065424 0.9282649 0.1188740 0.7193196 0.8459888
#RI      0.8980817 0.10042714 0.10191832 0.7960980 0.7386377 0.7947327 0.8958291
#GDP     0.8688597 0.13666802 0.13114028 0.8789344 0.7813226 0.8719163 0.9395573
#RR      0.9406754 0.11782517 0.05932461 0.7948018 0.7880421 0.7948930 0.8937148
#Imp     0.8950366 0.07583124 0.10496342 0.6994729 0.6484558 0.7187216 0.8494242
#general 0.4034948 0.40349483 0.40349483 0.8750954 0.5143495 0.9365354 0.9543579









################################################################
### Analysis RST-PQ Questionnaire
### Satchell et al https://psyarxiv.com/zcjnk/download?format=pdf


Satchell_RST.PQ <- read.csv("Satchell_RST-PQ.csv")
colnames(Satchell_RST.PQ)
mydata <- as.data.frame(Satchell_RST.PQ[,5:77])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 components
# Eigenvalue 1 = 12.24; eigenvalue 2 = 7.6

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 components
# Eigenvalue 1 = 14.15; eigenvalue 2 = 8.91

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.17, RMSEA=.081, RMSR=.13, TLI=.377

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.208, RMSR=.16, TLI=.086

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 7 communities with response bias

# Give solution with 7 factors
fit4 <- fa(mydata,7)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.042, RMSR=.04, TLI=.827
#      MR1   MR2  MR6   MR4   MR3   MR5   MR7
#MR1  1.00 -0.17 0.25 -0.03  0.06  0.36 -0.12
#MR2 -0.17  1.00 0.21  0.22  0.15 -0.03  0.17
#MR6  0.25  0.21 1.00  0.22  0.08  0.21  0.11
#MR4 -0.03  0.22 0.22  1.00  0.20  0.02  0.05
#MR3  0.06  0.15 0.08  0.20  1.00  0.02 -0.03
#MR5  0.36 -0.03 0.21  0.02  0.02  1.00 -0.02
#MR7 -0.12  0.17 0.11  0.05 -0.03 -0.02  1.00

fit4 <- fa(rho,7,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.206, RMSR=.05, TLI=.086
#      MR1   MR2  MR6   MR5  MR3   MR4   MR7
#MR1  1.00 -0.17 0.24  0.36 0.09 -0.03 -0.13
#MR2 -0.17  1.00 0.22 -0.04 0.13  0.21  0.19
#MR6  0.24  0.22 1.00  0.19 0.06  0.22  0.13
#MR5  0.36 -0.04 0.19  1.00 0.02  0.01 -0.04
#MR3  0.09  0.13 0.06  0.02 1.00  0.18  0.04
#MR4 -0.03  0.21 0.22  0.01 0.18  1.00  0.08
#MR7 -0.13  0.19 0.13 -0.04 0.04  0.08  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ RST1_BIS+RST2_BIS+RST3_BASRR+RST4_BASRR+RST5_BASGDP+RST6_DF+     
           RST7_BIS+RST8_BIS+RST9_BASRR+RST10_FFA+RST11_BIS+RST12_BASRI+ 
           RST13_BASGDP+RST14_DF+RST15_BASRI+RST16_BASRI+RST17_BASRI+RST18_BASRR+ 
           RST19_DF+RST20_BIS+RST21_BIS+RST22_FFA+RST23_BASGDP+RST24_DF+    
           RST25_BIS+RST26_BASImp+RST27_BASRR+RST28_BASRR+RST29_BASRR+RST30_BASRI+ 
           RST31_DF+RST32_BASImp+RST33_BASImp+RST34_BIS+RST35_BASRR+RST36_BASGDP+
           RST37_BASRI+RST38_BIS+RST39_BIS+RST40_DF+RST41_BASRI+RST42_BASRR+ 
           RST43_BASRR+RST44_BASImp+RST45_DF+RST46_DF+RST47_FFA+RST48_BASImp+
           RST49_BASGDP+RST50_BIS+RST51_BIS+RST52_BASImp+RST53_FFA+RST54_FFA+   
           RST55_BIS+RST56_FFA+RST57_BIS+RST58_BIS+RST59_BASImp+RST60_FFA+   
           RST61_BASImp+RST62_BASGDP+RST63_BIS+RST64_BIS+RST65_BIS+RST66_FFA+   
           RST67_FFA+RST68_BIS+RST69_BIS+RST70_FFA+RST71_BIS+RST72_BIS+RST73_BASGDP
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.690       0.510
#Tucker-Lewis Index (TLI)                       0.681       0.496
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.162       0.089
#Robust RMSEA                                                  NA
#SRMR                                           0.157       0.157

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .100

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.379       0.377
#Tucker-Lewis Index (TLI)                       0.362       0.359
#Robust Comparative Fit Index (CFI)                         0.382
#Robust Tucker-Lewis Index (TLI)                            0.364
#RMSEA                                          0.090       0.088
#Robust RMSEA                                               0.089
#SRMR                                           0.132       0.132

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .068

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 7 factors based on theoretical analysis
EGAmodel= '
 BIS  =~   RST1_BIS+RST2_BIS+RST7_BIS+RST8_BIS+RST11_BIS+RST20_BIS+RST21_BIS+RST25_BIS+RST34_BIS+
           RST38_BIS+RST39_BIS+RST50_BIS+RST51_BIS+RST55_BIS+RST57_BIS+RST58_BIS+RST63_BIS+
           RST64_BIS+RST65_BIS+RST68_BIS+RST69_BIS+RST71_BIS+RST72_BIS
 BASRR =~  RST3_BASRR+RST4_BASRR+RST9_BASRR+RST18_BASRR+RST27_BASRR+RST28_BASRR+RST29_BASRR+
           RST35_BASRR+RST42_BASRR+RST43_BASRR
 BASGDP =~ RST5_BASGDP+RST13_BASGDP+RST23_BASGDP+RST36_BASGDP+RST49_BASGDP+RST62_BASGDP+RST73_BASGDP
 BASDF  =~ RST6_DF+RST14_DF+RST19_DF+RST24_DF+RST31_DF+RST40_DF+RST45_DF+RST46_DF
 BASFFA =~ RST10_FFA+RST22_FFA+RST47_FFA+RST53_FFA+RST54_FFA+RST56_FFA+RST60_FFA+   
           RST66_FFA+RST67_FFA+RST70_FFA
 BASRI  =~ RST12_BASRI+RST15_BASRI+RST16_BASRI+RST17_BASRI+RST30_BASRI+RST37_BASRI+RST41_BASRI
 BASImp =~ RST26_BASImp+RST32_BASImp+RST33_BASImp+RST44_BASImp+RST48_BASImp+
           RST52_BASImp+RST59_BASImp+RST61_BASImp
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.911       0.834
#Tucker-Lewis Index (TLI)                       0.907       0.828
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.087       0.052
#Robust RMSEA                                                  NA
#SRMR                                           0.109       0.109

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  BIS ~~                                                                
#    BASRR             0.295    0.059    5.011    0.000    0.295    0.295
#    BASGDP           -0.172    0.062   -2.764    0.006   -0.172   -0.172
#    BASDF             0.108    0.066    1.640    0.101    0.108    0.108
#    BASFFA            0.529    0.052   10.227    0.000    0.529    0.529
#    BASRI            -0.214    0.064   -3.332    0.001   -0.214   -0.214
#    BASImp            0.321    0.069    4.669    0.000    0.321    0.321
#  BASRR ~~                                                              
#    BASGDP            0.391    0.060    6.525    0.000    0.391    0.391
#    BASDF             0.348    0.065    5.376    0.000    0.348    0.348
#    BASFFA            0.278    0.068    4.084    0.000    0.278    0.278
#    BASRI             0.471    0.056    8.369    0.000    0.471    0.471
#    BASImp            0.628    0.051   12.357    0.000    0.628    0.628
#  BASGDP ~~                                                             
#    BASDF             0.237    0.070    3.367    0.001    0.237    0.237
#    BASFFA           -0.074    0.068   -1.090    0.276   -0.074   -0.074
#    BASRI             0.603    0.051   11.733    0.000    0.603    0.603
#    BASImp            0.233    0.071    3.293    0.001    0.233    0.233
#  BASDF ~~                                                              
#    BASFFA            0.044    0.071    0.620    0.535    0.044    0.044
#    BASRI             0.209    0.068    3.075    0.002    0.209    0.209
#    BASImp            0.534    0.060    8.931    0.000    0.534    0.534
#  BASFFA ~~                                                             
#    BASRI            -0.154    0.064   -2.412    0.016   -0.154   -0.154
#    BASImp            0.324    0.073    4.432    0.000    0.324    0.324
#  BASRI ~~                                                              
#    BASImp            0.546    0.059    9.207    0.000    0.546    0.546


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.703       0.705
#Tucker-Lewis Index (TLI)                       0.692       0.694
#Robust Comparative Fit Index (CFI)                         0.708
#Robust Tucker-Lewis Index (TLI)                            0.698
#RMSEA                                          0.062       0.060
#Robust RMSEA                                               0.061
#SRMR                                           0.095       0.095

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  BIS ~~                                                                
#    BASRR             0.284    0.087    3.255    0.001    0.284    0.284
#    BASGDP           -0.161    0.072   -2.250    0.024   -0.161   -0.161
#    BASDF             0.107    0.082    1.309    0.190    0.107    0.107
#    BASFFA            0.520    0.066    7.844    0.000    0.520    0.520
#    BASRI            -0.200    0.081   -2.463    0.014   -0.200   -0.200
#    BASImp            0.277    0.099    2.782    0.005    0.277    0.277
#  BASRR ~~                                                              
#    BASGDP            0.392    0.078    5.054    0.000    0.392    0.392
#    BASDF             0.322    0.072    4.487    0.000    0.322    0.322
#    BASFFA            0.299    0.084    3.567    0.000    0.299    0.299
#    BASRI             0.454    0.080    5.653    0.000    0.454    0.454
#    BASImp            0.596    0.078    7.689    0.000    0.596    0.596
#  BASGDP ~~                                                             
#    BASDF             0.223    0.080    2.771    0.006    0.223    0.223
#    BASFFA           -0.074    0.086   -0.851    0.395   -0.074   -0.074
#    BASRI             0.566    0.079    7.153    0.000    0.566    0.566
#    BASImp            0.200    0.091    2.207    0.027    0.200    0.200
#  BASDF ~~                                                              
#    BASFFA            0.038    0.089    0.432    0.666    0.038    0.038
#    BASRI             0.183    0.086    2.139    0.032    0.183    0.183
#    BASImp            0.488    0.078    6.281    0.000    0.488    0.488
#  BASFFA ~~                                                             
#    BASRI            -0.157    0.084   -1.870    0.061   -0.157   -0.157
#    BASImp            0.250    0.112    2.230    0.026    0.250    0.250
#  BASRI ~~                                                              
#    BASImp            0.530    0.098    5.408    0.000    0.530    0.530

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .326


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.662       0.665
#Tucker-Lewis Index (TLI)                       0.652       0.655
#Robust Comparative Fit Index (CFI)                         0.667
#Robust Tucker-Lewis Index (TLI)                            0.658
#RMSEA                                          0.066       0.064
#Robust RMSEA                                               0.065
#SRMR                                           0.139       0.139

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .336

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 BIS  =~   RST1_BIS+RST2_BIS+RST7_BIS+RST8_BIS+RST11_BIS+RST20_BIS+RST21_BIS+RST25_BIS+RST34_BIS+
           RST38_BIS+RST39_BIS+RST50_BIS+RST51_BIS+RST55_BIS+RST57_BIS+RST58_BIS+RST63_BIS+
           RST64_BIS+RST65_BIS+RST68_BIS+RST69_BIS+RST71_BIS+RST72_BIS
 BASRR =~  RST3_BASRR+RST4_BASRR+RST9_BASRR+RST18_BASRR+RST27_BASRR+RST28_BASRR+RST29_BASRR+
           RST35_BASRR+RST42_BASRR+RST43_BASRR
 BASGDP =~ RST5_BASGDP+RST13_BASGDP+RST23_BASGDP+RST36_BASGDP+RST49_BASGDP+RST62_BASGDP+RST73_BASGDP
 BASDF  =~ RST6_DF+RST14_DF+RST19_DF+RST24_DF+RST31_DF+RST40_DF+RST45_DF+RST46_DF
 BASFFA =~ RST10_FFA+RST22_FFA+RST47_FFA+RST53_FFA+RST54_FFA+RST56_FFA+RST60_FFA+   
           RST66_FFA+RST67_FFA+RST70_FFA
 BASRI  =~ RST12_BASRI+RST15_BASRI+RST16_BASRI+RST17_BASRI+RST30_BASRI+RST37_BASRI+RST41_BASRI
 BASImp =~ RST26_BASImp+RST32_BASImp+RST33_BASImp+RST44_BASImp+RST48_BASImp+
           RST52_BASImp+RST59_BASImp+RST61_BASImp
 general=~ RST1_BIS+RST2_BIS+RST3_BASRR+RST4_BASRR+RST5_BASGDP+RST6_DF+     
           RST7_BIS+RST8_BIS+RST9_BASRR+RST10_FFA+RST11_BIS+RST12_BASRI+ 
           RST13_BASGDP+RST14_DF+RST15_BASRI+RST16_BASRI+RST17_BASRI+RST18_BASRR+ 
           RST19_DF+RST20_BIS+RST21_BIS+RST22_FFA+RST23_BASGDP+RST24_DF+    
           RST25_BIS+RST26_BASImp+RST27_BASRR+RST28_BASRR+RST29_BASRR+RST30_BASRI+ 
           RST31_DF+RST32_BASImp+RST33_BASImp+RST34_BIS+RST35_BASRR+RST36_BASGDP+
           RST37_BASRI+RST38_BIS+RST39_BIS+RST40_DF+RST41_BASRI+RST42_BASRR+ 
           RST43_BASRR+RST44_BASImp+RST45_DF+RST46_DF+RST47_FFA+RST48_BASImp+
           RST49_BASGDP+RST50_BIS+RST51_BIS+RST52_BASImp+RST53_FFA+RST54_FFA+   
           RST55_BIS+RST56_FFA+RST57_BIS+RST58_BIS+RST59_BASImp+RST60_FFA+   
           RST61_BASImp+RST62_BASGDP+RST63_BIS+RST64_BIS+RST65_BIS+RST66_FFA+   
           RST67_FFA+RST68_BIS+RST69_BIS+RST70_FFA+RST71_BIS+RST72_BIS+RST73_BASGDP
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE,estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.724       0.719
#Tucker-Lewis Index (TLI)                       0.707       0.702
#Robust Comparative Fit Index (CFI)                         0.725
#Robust Tucker-Lewis Index (TLI)                            0.709
#RMSEA                                          0.061       0.060
#Robust RMSEA                                               0.060
#SRMR                                           0.103       0.103

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .376

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  BIS =~                                                                
#    RST1_BIS          0.322    0.182    1.772    0.076    0.322    0.391
#    RST2_BIS          0.399    0.118    3.371    0.001    0.399    0.437
#    RST7_BIS          0.517    0.150    3.447    0.001    0.517    0.506
#    RST8_BIS          0.293    0.161    1.819    0.069    0.293    0.285
#    RST11_BIS         0.317    0.149    2.128    0.033    0.317    0.299
#    RST20_BIS         0.073    0.238    0.307    0.758    0.073    0.078
#    RST21_BIS         0.037    0.217    0.172    0.863    0.037    0.037
#    RST25_BIS         0.552    0.165    3.354    0.001    0.552    0.554
#    RST34_BIS         0.243    0.196    1.240    0.215    0.243    0.239
#    RST38_BIS         0.206    0.299    0.690    0.490    0.206    0.210
#    RST39_BIS         0.302    0.292    1.033    0.301    0.302    0.269
#    RST50_BIS        -0.087    0.383   -0.228    0.820   -0.087   -0.084
#    RST51_BIS         0.199    0.315    0.633    0.527    0.199    0.201
#    RST55_BIS         0.341    0.315    1.082    0.279    0.341    0.310
#    RST57_BIS         0.074    0.163    0.453    0.651    0.074    0.090
#    RST58_BIS         0.280    0.273    1.024    0.306    0.280    0.294
#    RST63_BIS        -0.251    0.194   -1.294    0.196   -0.251   -0.255
#    RST64_BIS         0.174    0.285    0.612    0.541    0.174    0.171
#    RST65_BIS         0.371    0.185    2.000    0.045    0.371    0.335
#    RST68_BIS         0.174    0.180    0.965    0.335    0.174    0.173
#    RST69_BIS         0.338    0.208    1.626    0.104    0.338    0.326
#    RST71_BIS        -0.120    0.269   -0.445    0.656   -0.120   -0.120
#    RST72_BIS        -0.200    0.200   -0.997    0.319   -0.200   -0.219
#  BASRR =~                                                              
#    RST3_BASRR        0.345    0.056    6.172    0.000    0.345    0.450
#    RST4_BASRR        0.363    0.069    5.267    0.000    0.363    0.423
#    RST9_BASRR        0.453    0.051    8.870    0.000    0.453    0.515
#    RST18_BASRR       0.543    0.057    9.492    0.000    0.543    0.661
#    RST27_BASRR       0.422    0.051    8.220    0.000    0.422    0.505
#    RST28_BASRR       0.513    0.067    7.609    0.000    0.513    0.527
#    RST29_BASRR       0.478    0.070    6.788    0.000    0.478    0.544
#    RST35_BASRR       0.528    0.058    9.030    0.000    0.528    0.610
#    RST42_BASRR       0.548    0.062    8.818    0.000    0.548    0.600
#    RST43_BASRR       0.609    0.055   11.167    0.000    0.609    0.710
#  BASGDP =~                                                             
#    RST5_BASGDP       0.584    0.045   12.941    0.000    0.584    0.749
#    RST13_BASGDP      0.585    0.051   11.531    0.000    0.585    0.710
#    RST23_BASGDP      0.505    0.050   10.149    0.000    0.505    0.644
#    RST36_BASGDP      0.652    0.056   11.668    0.000    0.652    0.737
#    RST49_BASGDP      0.708    0.043   16.465    0.000    0.708    0.819
#    RST62_BASGDP      0.297    0.058    5.130    0.000    0.297    0.375
#    RST73_BASGDP      0.662    0.047   14.084    0.000    0.662    0.742
#  BASDF =~                                                              
#    RST6_DF           0.647    0.048   13.432    0.000    0.647    0.715
#    RST14_DF          0.498    0.060    8.315    0.000    0.498    0.551
#    RST19_DF          0.376    0.054    6.910    0.000    0.376    0.493
#    RST24_DF          0.773    0.042   18.396    0.000    0.773    0.846
#    RST31_DF          0.519    0.053    9.853    0.000    0.519    0.601
#    RST40_DF          0.638    0.058   10.913    0.000    0.638    0.634
#    RST45_DF          0.410    0.058    7.083    0.000    0.410    0.424
#    RST46_DF          0.341    0.046    7.372    0.000    0.341    0.511
#  BASFFA =~                                                             
#    RST10_FFA         0.530    0.105    5.048    0.000    0.530    0.534
#    RST22_FFA         0.333    0.097    3.423    0.001    0.333    0.354
#    RST47_FFA         0.358    0.098    3.658    0.000    0.358    0.342
#    RST53_FFA         0.427    0.112    3.804    0.000    0.427    0.423
#    RST54_FFA         0.414    0.103    4.036    0.000    0.414    0.374
#    RST56_FFA         0.380    0.109    3.493    0.000    0.380    0.372
#    RST60_FFA         0.480    0.085    5.668    0.000    0.480    0.470
#    RST66_FFA         0.672    0.140    4.796    0.000    0.672    0.563
#    RST67_FFA         0.385    0.095    4.039    0.000    0.385    0.337
#    RST70_FFA         0.374    0.082    4.548    0.000    0.374    0.416
#  BASRI =~                                                              
#    RST12_BASRI       0.395    0.074    5.324    0.000    0.395    0.411
#    RST15_BASRI       0.383    0.075    5.119    0.000    0.383    0.407
#    RST16_BASRI       0.694    0.049   14.235    0.000    0.694    0.779
#    RST17_BASRI       0.496    0.057    8.758    0.000    0.496    0.563
#    RST30_BASRI       0.344    0.068    5.073    0.000    0.344    0.371
#    RST37_BASRI       0.692    0.057   12.053    0.000    0.692    0.747
#    RST41_BASRI       0.472    0.055    8.628    0.000    0.472    0.578
#  BASImp =~                                                             
#    RST26_BASImp      0.430    0.077    5.604    0.000    0.430    0.456
#    RST32_BASImp      0.475    0.084    5.672    0.000    0.475    0.464
#    RST33_BASImp      0.552    0.069    7.949    0.000    0.552    0.641
#    RST44_BASImp      0.618    0.060   10.314    0.000    0.618    0.705
#    RST48_BASImp      0.391    0.077    5.095    0.000    0.391    0.413
#    RST52_BASImp      0.537    0.078    6.868    0.000    0.537    0.477
#    RST59_BASImp      0.385    0.074    5.173    0.000    0.385    0.392
#    RST61_BASImp      0.324    0.076    4.245    0.000    0.324    0.366
#  general =~                                                            
#    RST1_BIS          0.454    0.105    4.300    0.000    0.454    0.551
#    RST2_BIS          0.582    0.101    5.778    0.000    0.582    0.638
#    RST3_BASRR       -0.137    0.079   -1.738    0.082   -0.137   -0.178
#    RST4_BASRR        0.299    0.069    4.354    0.000    0.299    0.348
#    RST5_BASGDP      -0.151    0.071   -2.116    0.034   -0.151   -0.194
#    RST6_DF           0.106    0.075    1.424    0.154    0.106    0.117
#    RST7_BIS          0.569    0.121    4.697    0.000    0.569    0.557
#    RST8_BIS          0.411    0.109    3.776    0.000    0.411    0.400
#    RST9_BASRR       -0.127    0.077   -1.645    0.100   -0.127   -0.144
#    RST10_FFA         0.213    0.089    2.390    0.017    0.213    0.215
#    RST11_BIS         0.419    0.101    4.162    0.000    0.419    0.395
#    RST12_BASRI      -0.291    0.090   -3.236    0.001   -0.291   -0.303
#    RST13_BASGDP     -0.142    0.091   -1.566    0.117   -0.142   -0.172
#    RST14_DF          0.121    0.071    1.709    0.087    0.121    0.134
#    RST15_BASRI      -0.163    0.087   -1.884    0.060   -0.163   -0.174
#    RST16_BASRI      -0.178    0.068   -2.601    0.009   -0.178   -0.199
#    RST17_BASRI       0.121    0.079    1.541    0.123    0.121    0.138
#    RST18_BASRR       0.070    0.086    0.812    0.417    0.070    0.085
#    RST19_DF         -0.093    0.080   -1.170    0.242   -0.093   -0.123
#    RST20_BIS         0.510    0.083    6.151    0.000    0.510    0.542
#    RST21_BIS         0.743    0.054   13.698    0.000    0.743    0.742
#    RST22_FFA         0.159    0.071    2.255    0.024    0.159    0.169
#    RST23_BASGDP     -0.128    0.077   -1.656    0.098   -0.128   -0.163
#    RST24_DF          0.021    0.082    0.263    0.793    0.021    0.023
#    RST25_BIS         0.559    0.125    4.458    0.000    0.559    0.561
#    RST26_BASImp      0.220    0.079    2.788    0.005    0.220    0.233
#    RST27_BASRR       0.201    0.067    2.987    0.003    0.201    0.240
#    RST28_BASRR       0.150    0.091    1.649    0.099    0.150    0.154
#    RST29_BASRR       0.285    0.077    3.714    0.000    0.285    0.324
#    RST30_BASRI      -0.021    0.096   -0.214    0.831   -0.021   -0.022
#    RST31_DF         -0.186    0.084   -2.219    0.026   -0.186   -0.216
#    RST32_BASImp      0.281    0.084    3.340    0.001    0.281    0.274
#    RST33_BASImp      0.108    0.074    1.458    0.145    0.108    0.125
#    RST34_BIS         0.676    0.091    7.413    0.000    0.676    0.664
#    RST35_BASRR       0.332    0.066    5.007    0.000    0.332    0.383
#    RST36_BASGDP     -0.084    0.099   -0.847    0.397   -0.084   -0.095
#    RST37_BASRI      -0.119    0.088   -1.353    0.176   -0.119   -0.129
#    RST38_BIS         0.622    0.111    5.621    0.000    0.622    0.634
#    RST39_BIS         0.784    0.119    6.615    0.000    0.784    0.700
#    RST40_DF          0.017    0.084    0.203    0.839    0.017    0.017
#    RST41_BASRI      -0.118    0.069   -1.716    0.086   -0.118   -0.144
#    RST42_BASRR       0.089    0.081    1.105    0.269    0.089    0.098
#    RST43_BASRR       0.117    0.079    1.472    0.141    0.117    0.136
#    RST44_BASImp     -0.006    0.078   -0.078    0.938   -0.006   -0.007
#    RST45_DF          0.321    0.074    4.357    0.000    0.321    0.332
#    RST46_DF          0.007    0.050    0.145    0.885    0.007    0.011
#    RST47_FFA         0.353    0.080    4.424    0.000    0.353    0.337
#    RST48_BASImp      0.321    0.073    4.405    0.000    0.321    0.339
#    RST49_BASGDP     -0.144    0.099   -1.445    0.148   -0.144   -0.166
#    RST50_BIS         0.646    0.092    7.034    0.000    0.646    0.623
#    RST51_BIS         0.576    0.115    4.985    0.000    0.576    0.582
#    RST52_BASImp     -0.162    0.088   -1.831    0.067   -0.162   -0.144
#    RST53_FFA         0.278    0.088    3.162    0.002    0.278    0.275
#    RST54_FFA         0.191    0.079    2.414    0.016    0.191    0.172
#    RST55_BIS         0.850    0.119    7.117    0.000    0.850    0.773
#    RST56_FFA         0.316    0.082    3.856    0.000    0.316    0.309
#    RST57_BIS         0.459    0.067    6.880    0.000    0.459    0.557
#    RST58_BIS         0.699    0.114    6.127    0.000    0.699    0.736
#    RST59_BASImp      0.043    0.081    0.529    0.597    0.043    0.044
#    RST60_FFA         0.339    0.074    4.603    0.000    0.339    0.332
#    RST61_BASImp      0.153    0.083    1.849    0.064    0.153    0.173
#    RST62_BASGDP      0.211    0.068    3.107    0.002    0.211    0.267
#    RST63_BIS         0.647    0.151    4.295    0.000    0.647    0.658
#    RST64_BIS         0.777    0.085    9.123    0.000    0.777    0.764
#    RST65_BIS         0.596    0.111    5.353    0.000    0.596    0.539
#    RST66_FFA         0.293    0.106    2.770    0.006    0.293    0.245
#    RST67_FFA         0.162    0.083    1.953    0.051    0.162    0.142
#    RST68_BIS         0.571    0.080    7.143    0.000    0.571    0.567
#    RST69_BIS         0.789    0.094    8.407    0.000    0.789    0.761
#    RST70_FFA         0.516    0.087    5.923    0.000    0.516    0.574
#    RST71_BIS         0.599    0.079    7.555    0.000    0.599    0.602
#    RST72_BIS         0.521    0.145    3.582    0.000    0.521    0.572
#    RST73_BASGDP     -0.139    0.080   -1.729    0.084   -0.139   -0.156

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.3923127      0.8321918      0.9148141      0.6359782 

#$FactorLevelIndices
#        ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#BIS     0.1760566 0.06594640 0.82394343 0.9474162 0.08835573 0.6952509 0.8409738
#BASRR   0.8530640 0.10916553 0.14693599 0.8388197 0.78536524 0.8313156 0.9157622
#BASGDP  0.9373917 0.11729945 0.06260833 0.8727428 0.85544841 0.8838003 0.9424990
#BASDF   0.9357058 0.10327748 0.06429422 0.8261900 0.82301664 0.8597226 0.9285338
#BASFFA  0.6656124 0.06277939 0.33438757 0.7759256 0.53950084 0.6972552 0.8419067
#BASRI   0.9128817 0.07931296 0.08711825 0.7758727 0.74128375 0.8123388 0.9054875
#BASImp  0.8657296 0.06990612 0.13427044 0.7429423 0.69408087 0.7535749 0.8710062
#general 0.3923127 0.39231266 0.39231266 0.9148141 0.63597818 0.9499062 0.9650912







################################################################
### Nuel et al https://online.ucpress.edu/collabra/article/8/1/34197/163623/The-Virtual-Reality-of-Social-Approach-Avoidance

#Experiment 1
Nuel_RST.PQ <- read.delim("Nuel_RST-PQ.txt")
colnames(Nuel_RST.PQ)
mydata <- as.data.frame(Nuel_RST.PQ[,2:55])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 4 response alternatives

#Experiment 2
Nuel_RST.PQ2 <- read.delim("Nuel_RST-PQ2.txt")
mydata2 <- as.data.frame(Nuel_RST.PQ2[,2:55])
colnames(mydata2)
mydata2 <- na.omit(mydata2)
min(mydata2)
max(mydata2)

mydata <- rbind(mydata,mydata2)
library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 8 factors and 6 components
# Eigenvalue 1 = 7.27; eigenvalue 2 = 4.30

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 components
# Eigenvalue 1 = 8.53; eigenvalue 2 = 5.23

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.13, RMSEA=.086, RMSR=.12, TLI=.332

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.16, RMSEA=.118, RMSR=.14, TLI=.24

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities with response bias

# Give solution with 6 factors
fit4 <- fa(mydata,6)
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.042, RMSR=.04, TLI=.841
#      MR1   MR2   MR6  MR4   MR3   MR5
#MR1  1.00 -0.20 -0.22 0.01  0.00  0.25
#MR2 -0.20  1.00  0.32 0.21 -0.03  0.01
#MR6 -0.22  0.32  1.00 0.12  0.22 -0.14
#MR4  0.01  0.21  0.12 1.00  0.11  0.17
#MR3  0.00 -0.03  0.22 0.11  1.00 -0.02
#MR5  0.25  0.01 -0.14 0.17 -0.02  1.00

fit4 <- fa(rho,6,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.079, RMSR=.04, TLI=.654
#      MR1   MR4   MR6  MR2   MR3   MR5
#MR1  1.00 -0.20 -0.22 0.00  0.00  0.25
#MR4 -0.20  1.00  0.31 0.24 -0.03  0.00
#MR6 -0.22  0.31  1.00 0.13  0.23 -0.14
#MR2  0.00  0.24  0.13 1.00  0.11  0.17
#MR3  0.00 -0.03  0.23 0.11  1.00 -0.03
#MR5  0.25  0.00 -0.14 0.17 -0.03  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ BIS1+GDP2+FFFS3+BIS4+RWR5+GDP6+BIS7+RWI8+BIS9+RWI10+FFFS11+
           BIS12+GDP13+RWR14+IMP15+GDP16+FFFS17+GDP18+RWI19+BIS20+IMP21+RWR22+ 
           IMP23+GDP24+FFFS25+BIS26+RWR27+FFFS28+RWI29+BIS30+RWI31+BIS32+IMP33+ 
           RWR34+FFFS35+IMP36+GDP37+BIS38+RWI39+BIS40+FFFS41+BIS42+IMP43+RWR44+ 
           FFFS45+BIS46+RWI47+FFFS48+GDP49+RWI50+BIS51+BIS52+FFFS53+GDP54 
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.656       0.488
#Tucker-Lewis Index (TLI)                       0.643       0.468
#Robust Comparative Fit Index (CFI)                         0.304
#Robust Tucker-Lewis Index (TLI)                            0.277
#RMSEA                                          0.132       0.091
#Robust RMSEA                                               0.117
#SRMR                                           0.140       0.140

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .131

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.360       0.355
#Tucker-Lewis Index (TLI)                       0.335       0.329
#Robust Comparative Fit Index (CFI)                         0.361
#Robust Tucker-Lewis Index (TLI)                            0.336
#RMSEA                                          0.088       0.086
#Robust RMSEA                                               0.088
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .074

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 6 factors based on theoretical analysis
EGAmodel= '
 BIS    =~ BIS1+BIS4+BIS7+BIS9+BIS12+BIS20+BIS26+BIS30+BIS32+BIS38+BIS40+BIS42+BIS46+BIS51+BIS52
 GDP    =~ GDP2+GDP6+GDP13+GDP16+GDP18+GDP24+GDP37+GDP49+GDP54 
 FFFS   =~ FFFS3+FFFS11+FFFS17+FFFS25+FFFS28+FFFS35+FFFS41+FFFS45+FFFS48+FFFS53
 RWR    =~ RWR5+RWR14+RWR22+RWR27+RWR34+RWR44 
 RWI    =~ RWI8+RWI10+RWI19+RWI29+RWI31+RWI39+RWI47+RWI50
 IMP    =~ IMP15+IMP21+IMP23+IMP33+IMP36+IMP43 
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.921       0.862
#Tucker-Lewis Index (TLI)                       0.917       0.855
#Robust Comparative Fit Index (CFI)                         0.703
#Robust Tucker-Lewis Index (TLI)                            0.688
#RMSEA                                          0.064       0.048
#Robust RMSEA                                               0.077
#SRMR                                           0.085       0.085

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  BIS ~~                                                                
#    GDP              -0.273    0.051   -5.384    0.000   -0.273   -0.273
#    FFFS              0.382    0.048    7.992    0.000    0.382    0.382
#    RWR              -0.093    0.057   -1.633    0.102   -0.093   -0.093
#    RWI              -0.292    0.051   -5.777    0.000   -0.292   -0.292
#    IMP               0.191    0.060    3.206    0.001    0.191    0.191
#  GDP ~~                                                                
#    FFFS             -0.045    0.055   -0.813    0.416   -0.045   -0.045
#    RWR               0.453    0.051    8.930    0.000    0.453    0.453
#    RWI               0.535    0.039   13.568    0.000    0.535    0.535
#    IMP              -0.059    0.057   -1.048    0.294   -0.059   -0.059
#  FFFS ~~                                                               
#    RWR               0.256    0.058    4.406    0.000    0.256    0.256
#    RWI              -0.242    0.054   -4.510    0.000   -0.242   -0.242
#    IMP               0.096    0.060    1.592    0.111    0.096    0.096
#  RWR ~~                                                                
#    RWI               0.387    0.051    7.580    0.000    0.387    0.387
#    IMP               0.316    0.059    5.353    0.000    0.316    0.316
#  RWI ~~                                                                
#    IMP               0.340    0.053    6.475    0.000    0.340    0.340

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .436

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.775       0.778
#Tucker-Lewis Index (TLI)                       0.763       0.766
#Robust Comparative Fit Index (CFI)                         0.780
#Robust Tucker-Lewis Index (TLI)                            0.768
#RMSEA                                          0.053       0.051
#Robust RMSEA                                               0.052
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .357

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  BIS ~~                                                                
#    GDP              -0.265    0.057   -4.690    0.000   -0.265   -0.265
#    FFFS              0.375    0.058    6.477    0.000    0.375    0.375
#    RWR              -0.078    0.062   -1.260    0.208   -0.078   -0.078
#    RWI              -0.298    0.054   -5.545    0.000   -0.298   -0.298
#    IMP               0.140    0.066    2.115    0.034    0.140    0.140
#  GDP ~~                                                                
#    FFFS             -0.035    0.065   -0.539    0.590   -0.035   -0.035
#    RWR               0.415    0.061    6.787    0.000    0.415    0.415
#    RWI               0.512    0.055    9.244    0.000    0.512    0.512
#    IMP              -0.050    0.067   -0.742    0.458   -0.050   -0.050
#  FFFS ~~                                                               
#    RWR               0.252    0.057    4.453    0.000    0.252    0.252
#    RWI              -0.236    0.061   -3.878    0.000   -0.236   -0.236
#    IMP               0.046    0.074    0.616    0.538    0.046    0.046
#  RWR ~~                                                                
#    RWI               0.360    0.059    6.089    0.000    0.360    0.360
#    IMP               0.263    0.059    4.454    0.000    0.263    0.263
#  RWI ~~                                                                
#    IMP               0.356    0.067    5.351    0.000    0.356    0.356


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.728       0.731
#Tucker-Lewis Index (TLI)                       0.718       0.720
#Robust Comparative Fit Index (CFI)                         0.733
#Robust Tucker-Lewis Index (TLI)                            0.722
#RMSEA                                          0.058       0.056
#Robust RMSEA                                               0.057
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .354

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 BIS    =~ BIS1+BIS4+BIS7+BIS9+BIS12+BIS20+BIS26+BIS30+BIS32+BIS38+BIS40+BIS42+BIS46+BIS51+BIS52
 GDP    =~ GDP2+GDP6+GDP13+GDP16+GDP18+GDP24+GDP37+GDP49+GDP54 
 FFFS   =~ FFFS3+FFFS11+FFFS17+FFFS25+FFFS28+FFFS35+FFFS41+FFFS45+FFFS48+FFFS53
 RWR    =~ RWR5+RWR14+RWR22+RWR27+RWR34+RWR44 
 RWI    =~ RWI8+RWI10+RWI19+RWI29+RWI31+RWI39+RWI47+RWI50
 IMP    =~ IMP15+IMP21+IMP23+IMP33+IMP36+IMP43 
 general=~ BIS1+GDP2+FFFS3+BIS4+RWR5+GDP6+BIS7+RWI8+BIS9+RWI10+FFFS11+
           BIS12+GDP13+RWR14+IMP15+GDP16+FFFS17+GDP18+RWI19+BIS20+IMP21+RWR22+ 
           IMP23+GDP24+FFFS25+BIS26+RWR27+FFFS28+RWI29+BIS30+RWI31+BIS32+IMP33+ 
           RWR34+FFFS35+IMP36+GDP37+BIS38+RWI39+BIS40+FFFS41+BIS42+IMP43+RWR44+ 
           FFFS45+BIS46+RWI47+FFFS48+GDP49+RWI50+BIS51+BIS52+FFFS53+GDP54 
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.792       0.784
#Tucker-Lewis Index (TLI)                       0.776       0.767
#Robust Comparative Fit Index (CFI)                         0.792
#Robust Tucker-Lewis Index (TLI)                            0.775
#RMSEA                                          0.051       0.051
#Robust RMSEA                                               0.051
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .372

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.28352267     0.82948987     0.79444732     0.08456246 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#BIS     0.8032133 0.23488051 0.19678666 0.8840836 0.72098072 0.8778294 0.9281599
#GDP     0.5075775 0.10815785 0.49242251 0.8584035 0.44086565 0.7463807 0.8464095
#FFFS    0.9107215 0.11493315 0.08927851 0.7301277 0.67905921 0.7546939 0.8680349
#RWR     0.8332000 0.09599409 0.16680001 0.7458639 0.63259460 0.7584873 0.8761457
#RWI     0.4381482 0.06716571 0.56185181 0.7904876 0.33633307 0.6281737 0.7861255
#IMP     0.9555470 0.09534603 0.04445296 0.7072135 0.70493196 0.7415462 0.8645730
#general 0.2835227 0.28352267 0.28352267 0.7944473 0.08456246 0.8786838 0.8991267











