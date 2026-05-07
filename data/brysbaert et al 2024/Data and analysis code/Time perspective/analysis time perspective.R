################################################################
### Analysis Questionnaires Zimbardo time perspective questionnaire
### Polish study with short form (20 questions addressing 4 dimensions)
### Bodecka
### Research on COVID situation in Poland https://osf.io/kney3

library(haven)
Time_perspective_Bodecka <- read_sav("Time_perspective_Bodecka.sav")
colnames(Time_perspective_Bodecka)
mydata  <- as.data.frame(Time_perspective_Bodecka[,41:60])
print(colnames(mydata))
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 3.38; eigenvalue 2 = 2.36

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 3.95; eigenvalue 2 = 2.62

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.17, RMSEA=.159, RMSR=.17, TLI=.257

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.185, RMSR=.19, TLI=.248

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities 

# Give solution with 4 factors
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.063, RMSR=.03, TLI=.882
#     MR3   MR2   MR4   MR1
#MR3  1.00  0.11 -0.09 -0.20
#MR2  0.11  1.00 -0.04  0.20
#MR4 -0.09 -0.04  1.00  0.39
#MR1 -0.20  0.20  0.39  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.53, RMSEA=.086, RMSR=.04, TLI=.838
#       MR2   MR4   MR3   MR1
#MR2  1.00 -0.08  0.10 -0.18
#MR4 -0.08  1.00 -0.04  0.42
#MR3  0.10 -0.04  1.00  0.20
#MR1 -0.18  0.42  0.20  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ ps_ztpi_r1+ps_ztpi_r2+ps_ztpi_r3+ps_ztpi_r4+ps_ztpi_r5+ps_ztpi_r6+ps_ztpi_r7+ps_ztpi_r8+ps_ztpi_r9+ps_ztpi_r10+
            ps_ztpi_r11+ps_ztpi_r12+ps_ztpi_r13+ps_ztpi_r14+ps_ztpi_r15+ps_ztpi_r16+ps_ztpi_r17+ps_ztpi_r18+ps_ztpi_r19+ps_ztpi_r20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.526       0.285
#Tucker-Lewis Index (TLI)                       0.470       0.201
#Robust Comparative Fit Index (CFI)                         0.291
#Robust Tucker-Lewis Index (TLI)                            0.207
#RMSEA                                          0.288       0.239
#Robust RMSEA                                               0.192
#SRMR                                           0.192       0.192

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .223

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.337       0.285
#Tucker-Lewis Index (TLI)                       0.259       0.201
#Robust Comparative Fit Index (CFI)                         0.336
#Robust Tucker-Lewis Index (TLI)                            0.258
#RMSEA                                          0.160       0.143
#Robust RMSEA                                               0.159
#SRMR                                           0.161       0.161

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .156

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
 factor1 =~ ps_ztpi_r19+ps_ztpi_r2+ps_ztpi_r20+ps_ztpi_r13+ps_ztpi_r12
 factor2 =~ ps_ztpi_r14+ps_ztpi_r9+ps_ztpi_r7+ps_ztpi_r17+ps_ztpi_r5
 factor3 =~ ps_ztpi_r15+ps_ztpi_r10+ps_ztpi_r11+ps_ztpi_r4+ps_ztpi_r16
 factor4 =~ ps_ztpi_r8+ps_ztpi_r1+ps_ztpi_r3+ps_ztpi_r6+ps_ztpi_r18
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.904
#Tucker-Lewis Index (TLI)                       0.949       0.888
#Robust Comparative Fit Index (CFI)                         0.841
#Robust Tucker-Lewis Index (TLI)                            0.816
#RMSEA                                          0.090       0.089
#Robust RMSEA                                               0.093
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .546

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.087    0.038   -2.264    0.024   -0.087   -0.087
#    factor3           0.158    0.038    4.149    0.000    0.158    0.158
#    factor4          -0.295    0.035   -8.471    0.000   -0.295   -0.295
#  factor2 ~~                                                            
#    factor3          -0.025    0.033   -0.761    0.447   -0.025   -0.025
#    factor4           0.563    0.031   17.940    0.000    0.563    0.563
#  factor3 ~~                                                            
#    factor4           0.246    0.037    6.589    0.000    0.246    0.246

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.876       0.878
#Tucker-Lewis Index (TLI)                       0.856       0.858
#Robust Comparative Fit Index (CFI)                         0.884
#Robust Tucker-Lewis Index (TLI)                            0.865
#RMSEA                                          0.070       0.060
#Robust RMSEA                                               0.068
#SRMR                                           0.064       0.064

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2          -0.087    0.038   -2.264    0.024   -0.087   -0.087
#  factor3           0.158    0.038    4.149    0.000    0.158    0.158
#  factor4          -0.295    0.035   -8.471    0.000   -0.295   -0.295
#factor2 ~~                                                            
#  factor3          -0.025    0.033   -0.761    0.447   -0.025   -0.025
#  factor4           0.563    0.031   17.940    0.000    0.563    0.563
#factor3 ~~                                                            
#  factor4           0.246    0.037    6.589    0.000    0.246    0.246

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .451


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.826       0.826
#Tucker-Lewis Index (TLI)                       0.806       0.806
#Robust Comparative Fit Index (CFI)                         0.834
#Robust Tucker-Lewis Index (TLI)                            0.814
#RMSEA                                          0.082       0.070
#Robust RMSEA                                               0.080
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .45

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ ps_ztpi_r19+ps_ztpi_r2+ps_ztpi_r20+ps_ztpi_r13+ps_ztpi_r12
 factor2 =~ ps_ztpi_r14+ps_ztpi_r9+ps_ztpi_r7+ps_ztpi_r17+ps_ztpi_r5
 factor3 =~ ps_ztpi_r15+ps_ztpi_r10+ps_ztpi_r11+ps_ztpi_r4+ps_ztpi_r16
 factor4 =~ ps_ztpi_r8+ps_ztpi_r1+ps_ztpi_r3+ps_ztpi_r6+ps_ztpi_r18
 general=~ ps_ztpi_r1+ps_ztpi_r2+ps_ztpi_r3+ps_ztpi_r4+ps_ztpi_r5+ps_ztpi_r6+ps_ztpi_r7+ps_ztpi_r8+ps_ztpi_r9+ps_ztpi_r10+
            ps_ztpi_r11+ps_ztpi_r12+ps_ztpi_r13+ps_ztpi_r14+ps_ztpi_r15+ps_ztpi_r16+ps_ztpi_r17+ps_ztpi_r18+ps_ztpi_r19+ps_ztpi_r20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.897       0.898
#Tucker-Lewis Index (TLI)                       0.870       0.871
#Robust Comparative Fit Index (CFI)                         0.904
#Robust Tucker-Lewis Index (TLI)                            0.879
#RMSEA                                          0.067       0.057
#Robust RMSEA                                               0.064
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .464

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    ps_ztpi_r19       0.849    0.043   19.551    0.000    0.849    0.736
#    ps_ztpi_r2        0.713    0.041   17.559    0.000    0.713    0.679
#    ps_ztpi_r20       0.733    0.038   19.144    0.000    0.733    0.707
#    ps_ztpi_r13       0.682    0.038   17.822    0.000    0.682    0.661
#    ps_ztpi_r12       0.657    0.047   14.121    0.000    0.657    0.591
#  factor2 =~                                                            
#    ps_ztpi_r14       0.659    0.037   17.996    0.000    0.659    0.722
#    ps_ztpi_r9        0.517    0.038   13.731    0.000    0.517    0.647
#    ps_ztpi_r7        0.561    0.042   13.329    0.000    0.561    0.604
#    ps_ztpi_r17       0.450    0.043   10.508    0.000    0.450    0.496
#    ps_ztpi_r5        0.276    0.042    6.559    0.000    0.276    0.326
#  factor3 =~                                                            
#    ps_ztpi_r15       0.647    0.057   11.277    0.000    0.647    0.609
#    ps_ztpi_r10       0.728    0.042   17.204    0.000    0.728    0.757
#    ps_ztpi_r11       0.601    0.053   11.245    0.000    0.601    0.617
#    ps_ztpi_r4        0.661    0.045   14.685    0.000    0.661    0.698
#    ps_ztpi_r16       0.543    0.038   14.314    0.000    0.543    0.567
#  factor4 =~                                                            
#    ps_ztpi_r8       -0.072    0.099   -0.732    0.464   -0.072   -0.080
#    ps_ztpi_r1       -0.272    0.111   -2.445    0.014   -0.272   -0.299
#    ps_ztpi_r3        0.030    0.109    0.273    0.785    0.030    0.031
#    ps_ztpi_r6        0.383    0.182    2.101    0.036    0.383    0.377
#    ps_ztpi_r18      -0.318    0.097   -3.284    0.001   -0.318   -0.326
#  general =~                                                            
#    ps_ztpi_r1        0.518    0.059    8.742    0.000    0.518    0.568
#    ps_ztpi_r2       -0.082    0.065   -1.266    0.206   -0.082   -0.078
#    ps_ztpi_r3        0.573    0.038   14.876    0.000    0.573    0.592
#    ps_ztpi_r4        0.112    0.045    2.512    0.012    0.112    0.118
#    ps_ztpi_r5        0.348    0.040    8.614    0.000    0.348    0.410
#    ps_ztpi_r6        0.785    0.065   12.132    0.000    0.785    0.772
#    ps_ztpi_r7        0.326    0.046    7.100    0.000    0.326    0.350
#    ps_ztpi_r8        0.690    0.038   18.151    0.000    0.690    0.763
#    ps_ztpi_r9        0.289    0.042    6.867    0.000    0.289    0.361
#    ps_ztpi_r10       0.037    0.045    0.816    0.414    0.037    0.038
#    ps_ztpi_r11       0.173    0.048    3.638    0.000    0.173    0.178
#    ps_ztpi_r12      -0.338    0.064   -5.269    0.000   -0.338   -0.304
#    ps_ztpi_r13      -0.162    0.056   -2.915    0.004   -0.162   -0.157
#    ps_ztpi_r14       0.276    0.044    6.298    0.000    0.276    0.303
#    ps_ztpi_r15       0.251    0.053    4.765    0.000    0.251    0.236
#    ps_ztpi_r16       0.135    0.047    2.870    0.004    0.135    0.142
#    ps_ztpi_r17       0.217    0.045    4.786    0.000    0.217    0.239
#    ps_ztpi_r18       0.431    0.072    5.979    0.000    0.431    0.443
#    ps_ztpi_r19      -0.374    0.059   -6.291    0.000   -0.374   -0.324
#    ps_ztpi_r20      -0.194    0.053   -3.627    0.000   -0.194   -0.187
#

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.3184896      0.7894737      0.8247272      0.3300914 

#$FactorLevelIndices
#         ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#factor1 0.8968785 0.24290629 0.10312146 0.8362304 0.762396684 0.8132190 0.9055791
#factor2 0.7440412 0.17574407 0.25595884 0.7922695 0.585016498 0.7388522 0.8666199
#factor3 0.9455599 0.22624864 0.05444015 0.8011280 0.764409593 0.7987390 0.8946456
#factor4 0.1443166 0.03661139 0.85568335 0.7919835 0.007045301 0.2806586 0.6724828
#general 0.3184896 0.31848960 0.31848960 0.8247272 0.330091398 0.8384868 0.9199046



################################################################
### Analysis Questionnaires Zimbardo time perspective questionnaire
### Turkish study Akirmak https://osf.io/k3nvh

ZTPI_Akirmak <- read.csv("ZTPI_Akirmak.csv")
colnames(ZTPI_Akirmak)
mydata  <- as.data.frame(ZTPI_Akirmak[,19:74])
print(colnames(mydata))
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 6 components
# Eigenvalue 1 = 6.34; eigenvalue 2 = 4.19

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 components
# Eigenvalue 1 = 7.06; eigenvalue 2 = 4.73

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.11, RMSEA=.081, RMSR=.12, TLI=.275

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.13, RMSEA=.115, RMSR=.14, TLI=.17

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities 

# Give solution with 5 factors
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.047, RMSR=.05, TLI=.75
#      MR1   MR3   MR2   MR4   MR5
#MR1  1.00 -0.13 -0.08 -0.19  0.17
#MR3 -0.13  1.00 -0.11  0.12 -0.22
#MR2 -0.08 -0.11  1.00  0.12  0.11
#MR4 -0.19  0.12  0.12  1.00  0.06
#MR5  0.17 -0.22  0.11  0.06  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.09, RMSR=.06, TLI=.483
#      MR1   MR3   MR2   MR4   MR5
#MR1  1.00 -0.12 -0.08 -0.16  0.17
#MR3 -0.12  1.00 -0.09  0.14 -0.22
#MR2 -0.08 -0.09  1.00  0.14  0.10
#MR4 -0.16  0.14  0.14  1.00  0.05
#MR5  0.17 -0.22  0.10  0.05  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ Q15_1+Q15_2+Q15_3+Q15_4+Q15_5+Q15_6+Q15_7+Q15_8+Q15_9+Q15_10+Q15_11+Q15_12+Q15_13+
           Q15_14+Q15_15+Q15_16+Q15_17+Q15_18+Q15_19+Q15_20+Q15_21+Q15_22+Q15_23+Q15_24+Q15_25+
           Q15_26+Q15_27+Q15_28+Q15_29+Q15_30+Q15_31+Q15_32+Q15_33+Q15_34+Q15_35+Q15_36+Q15_37+
           Q15_38+Q15_39+Q15_40+Q15_41+Q15_42+Q15_43+Q15_44+Q15_45+Q15_46+Q15_47+Q15_48+Q15_49+
           Q15_50+Q15_51+Q15_52+Q15_53+Q15_54+Q15_55+Q15_56
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.563       0.411
#Tucker-Lewis Index (TLI)                       0.546       0.388
#Robust Comparative Fit Index (CFI)                         0.221
#Robust Tucker-Lewis Index (TLI)                            0.192
#RMSEA                                          0.142       0.082
#Robust RMSEA                                               0.122
#SRMR                                           0.136       0.136

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .096

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.296       0.288
#Tucker-Lewis Index (TLI)                       0.270       0.262
#Robust Comparative Fit Index (CFI)                         0.298
#Robust Tucker-Lewis Index (TLI)                            0.272
#RMSEA                                          0.088       0.084
#Robust RMSEA                                               0.087
#SRMR                                           0.119       0.119

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .04

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors based on theoretical analysis
EGAmodel= '
 pastneg =~ Q15_4+Q15_5+Q15_16+Q15_22+Q15_27+Q15_33+Q15_34+Q15_36+Q15_50+Q15_54
 predhed =~ Q15_1+Q15_8+Q15_12+Q15_17+Q15_19+Q15_23+Q15_26+Q15_28+Q15_31+Q15_32+Q15_42+Q15_44+Q15_46+Q15_48+Q15_55
 F       =~ Q15_6+Q15_9+Q15_10+Q15_13+Q15_18+Q15_21+Q15_24+Q15_30+Q15_40+Q15_43+Q15_45+Q15_51+Q15_56
 PP      =~ Q15_2+Q15_7+Q15_11+Q15_15+Q15_20+Q15_25+Q15_29+Q15_41+Q15_49
 PF      =~ Q15_3+Q15_14+Q15_35+Q15_37+Q15_38+Q15_39+Q15_47+Q15_52+Q15_53
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.755       0.665
#Tucker-Lewis Index (TLI)                       0.744       0.650
#Robust Comparative Fit Index (CFI)                         0.431
#Robust Tucker-Lewis Index (TLI)                            0.405
#RMSEA                                          0.107       0.062
#Robust RMSEA                                               0.105
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .254

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  pastneg ~~                                                            
#    predhed           0.056    0.063    0.882    0.378    0.056    0.056
#    F                -0.201    0.058   -3.485    0.000   -0.201   -0.201
#    PP               -0.637    0.041  -15.529    0.000   -0.637   -0.637
#    PF                0.427    0.054    7.954    0.000    0.427    0.427
#  predhed ~~                                                            
#    F                -0.270    0.061   -4.444    0.000   -0.270   -0.270
#    PP                0.125    0.072    1.727    0.084    0.125    0.125
#    PF                0.309    0.064    4.856    0.000    0.309    0.309
#  F ~~                                                                  
#    PP                0.403    0.060    6.712    0.000    0.403    0.403
#    PF               -0.426    0.057   -7.464    0.000   -0.426   -0.426
#  PP ~~                                                                 
#    PF               -0.197    0.072   -2.748    0.006   -0.197   -0.197

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.565       0.570
#Tucker-Lewis Index (TLI)                       0.545       0.551
#Robust Comparative Fit Index (CFI)                         0.576
#Robust Tucker-Lewis Index (TLI)                            0.557
#RMSEA                                          0.069       0.066
#Robust RMSEA                                               0.067
#SRMR                                           0.101       0.101

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#pastneg ~~                                                            
#        predhed          -0.025    0.104   -0.244    0.807   -0.025   -0.025
#        F                -0.176    0.088   -2.002    0.045   -0.176   -0.176
#        PP               -0.553    0.120   -4.609    0.000   -0.553   -0.553
#        PF                0.381    0.087    4.358    0.000    0.381    0.381
#predhed ~~                                                            
#        F                -0.220    0.145   -1.521    0.128   -0.220   -0.220
#        PP                0.189    0.122    1.555    0.120    0.189    0.189
#        PF                0.261    0.109    2.400    0.016    0.261    0.261
#F ~~                                                                  
#        PP                0.359    0.108    3.332    0.001    0.359    0.359
#        PF               -0.415    0.088   -4.691    0.000   -0.415   -0.415
#PP ~~                                                                 
#        PF               -0.190    0.113   -1.682    0.093   -0.190   -0.190

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .23


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.536       0.543
#Tucker-Lewis Index (TLI)                       0.518       0.525
#Robust Comparative Fit Index (CFI)                         0.547
#Robust Tucker-Lewis Index (TLI)                            0.530
#RMSEA                                          0.072       0.067
#Robust RMSEA                                               0.069
#SRMR                                           0.120       0.120

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .25

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 pastneg =~ Q15_4+Q15_5+Q15_16+Q15_22+Q15_27+Q15_33+Q15_34+Q15_36+Q15_50+Q15_54
 predhed =~ Q15_1+Q15_8+Q15_12+Q15_17+Q15_19+Q15_23+Q15_26+Q15_28+Q15_31+Q15_32+Q15_42+Q15_44+Q15_46+Q15_48+Q15_55
 F       =~ Q15_6+Q15_9+Q15_10+Q15_13+Q15_18+Q15_21+Q15_24+Q15_30+Q15_40+Q15_43+Q15_45+Q15_51+Q15_56
 PP      =~ Q15_2+Q15_7+Q15_11+Q15_15+Q15_20+Q15_25+Q15_29+Q15_41+Q15_49
 PF      =~ Q15_3+Q15_14+Q15_35+Q15_37+Q15_38+Q15_39+Q15_47+Q15_52+Q15_53
 general=~ Q15_1+Q15_2+Q15_3+Q15_4+Q15_5+Q15_6+Q15_7+Q15_8+Q15_9+Q15_10+Q15_11+Q15_12+Q15_13+
           Q15_14+Q15_15+Q15_16+Q15_17+Q15_18+Q15_19+Q15_20+Q15_21+Q15_22+Q15_23+Q15_24+Q15_25+
           Q15_26+Q15_27+Q15_28+Q15_29+Q15_30+Q15_31+Q15_32+Q15_33+Q15_34+Q15_35+Q15_36+Q15_37+
           Q15_38+Q15_39+Q15_40+Q15_41+Q15_42+Q15_43+Q15_44+Q15_45+Q15_46+Q15_47+Q15_48+Q15_49+
           Q15_50+Q15_51+Q15_52+Q15_53+Q15_54+Q15_55+Q15_56
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.634       0.636
#Tucker-Lewis Index (TLI)                       0.605       0.608
#Robust Comparative Fit Index (CFI)                         0.644
#Robust Tucker-Lewis Index (TLI)                            0.616
#RMSEA                                          0.065       0.061
#Robust RMSEA                                               0.063
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .292

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  pastneg =~                                                            
#    Q15_4            -0.395    0.388   -1.018    0.309   -0.395   -0.365
#    Q15_5            -0.330    0.115   -2.873    0.004   -0.330   -0.317
#    Q15_16            0.260    0.307    0.847    0.397    0.260    0.214
#    Q15_22           -0.128    0.241   -0.531    0.595   -0.128   -0.109
#    Q15_27           -0.433    0.323   -1.341    0.180   -0.433   -0.338
#    Q15_33           -0.222    0.113   -1.972    0.049   -0.222   -0.235
#    Q15_34            0.174    0.339    0.515    0.607    0.174    0.140
#    Q15_36            0.125    0.293    0.428    0.669    0.125    0.103
#    Q15_50            0.261    0.387    0.674    0.501    0.261    0.212
#    Q15_54           -0.357    0.251   -1.424    0.155   -0.357   -0.334
#  predhed =~                                                            
#    Q15_1             0.210    0.064    3.261    0.001    0.210    0.254
#    Q15_8             0.359    0.087    4.125    0.000    0.359    0.369
#    Q15_12            0.136    0.076    1.788    0.074    0.136    0.149
#    Q15_17            0.553    0.079    7.029    0.000    0.553    0.562
#    Q15_19            0.550    0.084    6.509    0.000    0.550    0.488
#    Q15_23            0.232    0.082    2.833    0.005    0.232    0.223
#    Q15_26            0.475    0.055    8.623    0.000    0.475    0.606
#    Q15_28            0.381    0.076    4.990    0.000    0.381    0.382
#    Q15_31            0.596    0.107    5.557    0.000    0.596    0.587
#    Q15_32            0.412    0.068    6.044    0.000    0.412    0.481
#    Q15_42            0.674    0.110    6.108    0.000    0.674    0.624
#    Q15_44            0.331    0.097    3.414    0.001    0.331    0.308
#    Q15_46            0.620    0.081    7.683    0.000    0.620    0.589
#    Q15_48            0.415    0.085    4.874    0.000    0.415    0.377
#    Q15_55            0.197    0.063    3.151    0.002    0.197    0.261
#  F =~                                                                  
#    Q15_6             0.547    0.086    6.378    0.000    0.547    0.485
#    Q15_9            -0.580    0.078   -7.410    0.000   -0.580   -0.605
#    Q15_10            0.415    0.072    5.728    0.000    0.415    0.459
#    Q15_13            0.659    0.069    9.564    0.000    0.659    0.593
#    Q15_18            0.531    0.080    6.644    0.000    0.531    0.515
#    Q15_21            0.397    0.075    5.319    0.000    0.397    0.478
#    Q15_24           -0.587    0.077   -7.599    0.000   -0.587   -0.575
#    Q15_30            0.325    0.080    4.066    0.000    0.325    0.358
#    Q15_40            0.446    0.086    5.186    0.000    0.446    0.477
#    Q15_43            0.595    0.090    6.638    0.000    0.595    0.493
#    Q15_45            0.554    0.075    7.380    0.000    0.554    0.545
#    Q15_51            0.232    0.076    3.038    0.002    0.232    0.233
#    Q15_56            0.034    0.094    0.367    0.714    0.034    0.035
#  PP =~                                                                 
#    Q15_2             0.512    0.091    5.650    0.000    0.512    0.536
#    Q15_7             0.437    0.086    5.072    0.000    0.437    0.413
#    Q15_11            0.440    0.080    5.473    0.000    0.440    0.437
#    Q15_15            0.558    0.082    6.790    0.000    0.558    0.601
#    Q15_20            0.485    0.066    7.394    0.000    0.485    0.530
#    Q15_25           -0.257    0.076   -3.362    0.001   -0.257   -0.211
#    Q15_29            0.716    0.108    6.597    0.000    0.716    0.590
#    Q15_41           -0.386    0.096   -4.041    0.000   -0.386   -0.353
#    Q15_49            0.579    0.098    5.887    0.000    0.579    0.477
#  PF =~                                                                 
#    Q15_3             0.502    0.117    4.298    0.000    0.502    0.433
#    Q15_14            0.431    0.077    5.602    0.000    0.431    0.475
#    Q15_35            0.384    0.092    4.153    0.000    0.384    0.391
#    Q15_37            0.531    0.113    4.709    0.000    0.531    0.478
#    Q15_38            0.563    0.101    5.598    0.000    0.563    0.531
#    Q15_39            0.585    0.085    6.846    0.000    0.585    0.585
#    Q15_47            0.241    0.115    2.092    0.036    0.241    0.196
#    Q15_52            0.280    0.117    2.400    0.016    0.280    0.251
#    Q15_53            0.437    0.114    3.825    0.000    0.437    0.397
#  general =~                                                            
#    Q15_1             0.129    0.072    1.784    0.074    0.129    0.156
#    Q15_2             0.104    0.067    1.538    0.124    0.104    0.109
#    Q15_3             0.012    0.105    0.111    0.912    0.012    0.010
#    Q15_4            -0.654    0.124   -5.289    0.000   -0.654   -0.605
#    Q15_5            -0.264    0.148   -1.780    0.075   -0.264   -0.254
#    Q15_6            -0.074    0.107   -0.690    0.490   -0.074   -0.066
#    Q15_7             0.296    0.109    2.717    0.007    0.296    0.280
#    Q15_8            -0.096    0.083   -1.148    0.251   -0.096   -0.098
#    Q15_9            -0.088    0.082   -1.080    0.280   -0.088   -0.092
#    Q15_10            0.257    0.075    3.434    0.001    0.257    0.285
#    Q15_11            0.438    0.069    6.365    0.000    0.438    0.434
#    Q15_12           -0.045    0.074   -0.603    0.547   -0.045   -0.049
#    Q15_13            0.109    0.104    1.046    0.296    0.109    0.098
#    Q15_14           -0.105    0.072   -1.462    0.144   -0.105   -0.116
#    Q15_15            0.017    0.082    0.208    0.835    0.017    0.018
#    Q15_16           -0.875    0.109   -7.994    0.000   -0.875   -0.721
#    Q15_17            0.253    0.079    3.188    0.001    0.253    0.257
#    Q15_18            0.031    0.090    0.343    0.731    0.031    0.030
#    Q15_19            0.044    0.096    0.457    0.648    0.044    0.039
#    Q15_20            0.264    0.081    3.257    0.001    0.264    0.288
#    Q15_21            0.208    0.066    3.153    0.002    0.208    0.251
#    Q15_22           -0.670    0.087   -7.724    0.000   -0.670   -0.569
#    Q15_23           -0.230    0.086   -2.666    0.008   -0.230   -0.222
#    Q15_24           -0.036    0.115   -0.315    0.752   -0.036   -0.036
#    Q15_25           -0.857    0.065  -13.193    0.000   -0.857   -0.706
#    Q15_26            0.026    0.087    0.293    0.770    0.026    0.033
#    Q15_27           -0.793    0.148   -5.351    0.000   -0.793   -0.620
#    Q15_28            0.002    0.096    0.020    0.984    0.002    0.002
#    Q15_29           -0.121    0.098   -1.237    0.216   -0.121   -0.100
#    Q15_30            0.075    0.085    0.881    0.379    0.075    0.083
#    Q15_31           -0.084    0.085   -0.985    0.325   -0.084   -0.083
#    Q15_32            0.158    0.080    1.971    0.049    0.158    0.184
#    Q15_33           -0.354    0.113   -3.135    0.002   -0.354   -0.374
#    Q15_34           -1.027    0.074  -13.881    0.000   -1.027   -0.823
#    Q15_35           -0.155    0.101   -1.534    0.125   -0.155   -0.157
#    Q15_36           -0.732    0.089   -8.206    0.000   -0.732   -0.602
#    Q15_37           -0.385    0.092   -4.204    0.000   -0.385   -0.346
#    Q15_38           -0.305    0.113   -2.696    0.007   -0.305   -0.288
#    Q15_39           -0.154    0.076   -2.025    0.043   -0.154   -0.154
#    Q15_40            0.298    0.085    3.512    0.000    0.298    0.318
#    Q15_41           -0.234    0.091   -2.569    0.010   -0.234   -0.214
#    Q15_42           -0.085    0.089   -0.957    0.339   -0.085   -0.079
#    Q15_43            0.086    0.096    0.901    0.367    0.086    0.071
#    Q15_44           -0.161    0.105   -1.535    0.125   -0.161   -0.150
#    Q15_45            0.221    0.083    2.651    0.008    0.221    0.218
#    Q15_46           -0.011    0.101   -0.104    0.917   -0.011   -0.010
#    Q15_47           -0.297    0.103   -2.891    0.004   -0.297   -0.241
#    Q15_48            0.039    0.104    0.372    0.710    0.039    0.035
#    Q15_49            0.116    0.097    1.190    0.234    0.116    0.095
#    Q15_50           -1.009    0.093  -10.845    0.000   -1.009   -0.819
#    Q15_51           -0.146    0.086   -1.698    0.089   -0.146   -0.147
#    Q15_52           -0.062    0.088   -0.708    0.479   -0.062   -0.056
#    Q15_53           -0.084    0.107   -0.788    0.431   -0.084   -0.076
#    Q15_54           -0.504    0.143   -3.530    0.000   -0.504   -0.472
#    Q15_55           -0.034    0.059   -0.578    0.563   -0.034   -0.045
#    Q15_56            0.005    0.080    0.059    0.953    0.005    0.005

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.3517601      0.8051948      0.7370509      0.2426714 

#$FactorLevelIndices
#        ECV_SS    ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#pastneg 0.1480052 0.0409614 0.85199481 0.8629251 0.02585721 0.4175897 0.7404816
#predhed 0.9286643 0.1869006 0.07133568 0.7683864 0.76836862 0.8064728 0.8997836
#F       0.8932570 0.1856991 0.10674295 0.5766601 0.53140066 0.8034326 0.8979462
#PP      0.6875395 0.1286491 0.31246045 0.6027455 0.59999144 0.7381464 0.8633997
#PF      0.8350044 0.1060297 0.16499565 0.6958920 0.60768311 0.6880821 0.8338963
#general 0.3517601 0.3517601 0.35176015 0.7370509 0.24267139 0.9125692 0.9580311







################################################################
### Analysis Questionnaires Zimbardo time perspective questionnaire
### Ceccato et al. 2021 https://www.sciencedirect.com/science/article/pii/S2352340921001761

library(readxl)
ZTPI_Ceccato <- read_excel("ZTPI_Ceccato.xlsx")
colnames(ZTPI_Ceccato)
mydata  <- as.data.frame(ZTPI_Ceccato[,8:25])
print(colnames(mydata))
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 2.96; eigenvalue 2 = 1.34

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 3.47; eigenvalue 2 = 1.59

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.111, RMSR=.11, TLI=.461

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.135, RMSR=.13, TLI=.423

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities 

# Give solution with 5 factors
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.04, RMSR=.02, TLI=.931
#      MR1   MR4  MR3   MR2   MR5
#MR1  1.00 -0.23 0.25  0.05  0.31
#MR4 -0.23  1.00 0.06  0.28 -0.28
#MR3  0.25  0.06 1.00  0.13  0.08
#MR2  0.05  0.28 0.13  1.00 -0.05
#MR5  0.31 -0.28 0.08 -0.05  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.052, RMSR=.02, TLI=.916
#      MR1   MR5   MR2  MR3   MR4
#MR1  1.00 -0.23  0.06 0.24  0.30
#MR5 -0.23  1.00  0.30 0.05 -0.32
#MR2  0.06  0.30  1.00 0.10 -0.08
#MR3  0.24  0.05  0.10 1.00  0.09
#MR4  0.30 -0.32 -0.08 0.09  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ ZTPI_1+ZTPI_2+ZTPI_3+ZTPI_4+ZTPI_5+ZTPI_6+ZTPI_7+ZTPI_8+ZTPI_9+
           ZTPI_10+ZTPI_11+ZTPI_12+ZTPI_13+ZTPI_14+ZTPI_15+ZTPI_16+ZTPI_17+ZTPI_18
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.707       0.524
#Tucker-Lewis Index (TLI)                       0.668       0.461
#Robust Comparative Fit Index (CFI)                         0.490
#Robust Tucker-Lewis Index (TLI)                            0.422
#RMSEA                                          0.144       0.147
#Robust RMSEA                                               0.135
#SRMR                                           0.120       0.120

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .163

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.525       0.517
#Tucker-Lewis Index (TLI)                       0.462       0.453
#Robust Comparative Fit Index (CFI)                         0.526
#Robust Tucker-Lewis Index (TLI)                            0.463
#RMSEA                                          0.111       0.103
#Robust RMSEA                                               0.111
#SRMR                                           0.101       0.101

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .102

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors based on theoretical analysis
EGAmodel= '
 factor1 =~ ZTPI_2+ZTPI_17+ZTPI_10+ZTPI_9+ZTPI_15+ZTPI_4
 factor2 =~ ZTPI_12+ZTPI_5+ZTPI_14
 factor3 =~ ZTPI_1+ZTPI_3+ZTPI_7
 factor4 =~ ZTPI_13+ZTPI_8+ZTPI_16
 factor5 =~ ZTPI_11+ZTPI_6+ZTPI_18
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.896       0.809
#Tucker-Lewis Index (TLI)                       0.872       0.767
#Robust Comparative Fit Index (CFI)                         0.787
#Robust Tucker-Lewis Index (TLI)                            0.740
#RMSEA                                          0.090       0.097
#Robust RMSEA                                               0.091
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .426

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.376    0.018  -20.478    0.000   -0.376   -0.376
#    factor3           0.082    0.020    4.029    0.000    0.082    0.082
#    factor4           0.344    0.019   18.270    0.000    0.344    0.344
#    factor5           0.665    0.016   41.468    0.000    0.665    0.665
#  factor2 ~~                                                            
#    factor3           0.352    0.021   16.487    0.000    0.352    0.352
#    factor4           0.022    0.022    1.022    0.307    0.022    0.022
#    factor5          -0.364    0.021  -17.417    0.000   -0.364   -0.364
#  factor3 ~~                                                            
#    factor4           0.235    0.023   10.436    0.000    0.235    0.235
#    factor5           0.012    0.024    0.515    0.607    0.012    0.012
#  factor4 ~~                                                            
#    factor5           0.254    0.022   11.328    0.000    0.254    0.254

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.823       0.817
#Tucker-Lewis Index (TLI)                       0.783       0.776
#Robust Comparative Fit Index (CFI)                         0.824
#Robust Tucker-Lewis Index (TLI)                            0.784
#RMSEA                                          0.071       0.066
#Robust RMSEA                                               0.071
#SRMR                                           0.067       0.067

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2          -0.376    0.018  -20.478    0.000   -0.376   -0.376
#    factor3           0.082    0.020    4.029    0.000    0.082    0.082
#    factor4           0.344    0.019   18.270    0.000    0.344    0.344
#    factor5           0.665    0.016   41.468    0.000    0.665    0.665
#  factor2 ~~                                                            
#    factor3           0.352    0.021   16.487    0.000    0.352    0.352
#    factor4           0.022    0.022    1.022    0.307    0.022    0.022
#    factor5          -0.364    0.021  -17.417    0.000   -0.364   -0.364
#  factor3 ~~                                                            
#    factor4           0.235    0.023   10.436    0.000    0.235    0.235
#    factor5           0.012    0.024    0.515    0.607    0.012    0.012
#  factor4 ~~                                                            
#    factor5           0.254    0.022   11.328    0.000    0.254    0.254

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .348

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.717       0.713
#Tucker-Lewis Index (TLI)                       0.679       0.675
#Robust Comparative Fit Index (CFI)                         0.718
#Robust Tucker-Lewis Index (TLI)                            0.680
#RMSEA                                          0.086       0.079
#Robust RMSEA                                               0.086
#SRMR                                           0.119       0.119

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .377

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ ZTPI_2+ZTPI_17+ZTPI_10+ZTPI_9+ZTPI_15+ZTPI_4
 factor2 =~ ZTPI_12+ZTPI_5+ZTPI_14
 factor3 =~ ZTPI_1+ZTPI_3+ZTPI_7
 factor4 =~ ZTPI_13+ZTPI_8+ZTPI_16
 factor5 =~ ZTPI_11+ZTPI_6+ZTPI_18
 general=~ ZTPI_1+ZTPI_2+ZTPI_3+ZTPI_4+ZTPI_5+ZTPI_6+ZTPI_7+ZTPI_8+ZTPI_9+
           ZTPI_10+ZTPI_11+ZTPI_12+ZTPI_13+ZTPI_14+ZTPI_15+ZTPI_16+ZTPI_17+ZTPI_18
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.871       0.865
#Tucker-Lewis Index (TLI)                       0.831       0.824
#Robust Comparative Fit Index (CFI)                         0.872
#Robust Tucker-Lewis Index (TLI)                            0.832
#RMSEA                                          0.062       0.058
#Robust RMSEA                                               0.062
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .384

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    ZTPI_2            0.712    0.113    6.316    0.000    0.712    0.578
#    ZTPI_17           0.409    0.038   10.880    0.000    0.409    0.319
#    ZTPI_10           0.639    0.081    7.866    0.000    0.639    0.488
#    ZTPI_9            0.199    0.049    4.091    0.000    0.199    0.163
#    ZTPI_15           0.103    0.063    1.620    0.105    0.103    0.086
#    ZTPI_4            0.119    0.044    2.714    0.007    0.119    0.102
#  factor2 =~                                                            
#    ZTPI_12           0.687    0.027   25.348    0.000    0.687    0.675
#    ZTPI_5            0.595    0.026   23.275    0.000    0.595    0.586
#    ZTPI_14           0.473    0.025   19.217    0.000    0.473    0.411
#  factor3 =~                                                            
#    ZTPI_1            0.643    0.022   29.304    0.000    0.643    0.630
#    ZTPI_3            0.695    0.023   29.786    0.000    0.695    0.603
#    ZTPI_7            0.688    0.023   29.851    0.000    0.688    0.592
#  factor4 =~                                                            
#    ZTPI_13           0.696    0.027   25.765    0.000    0.696    0.597
#    ZTPI_8            0.705    0.028   25.512    0.000    0.705    0.621
#    ZTPI_16           0.566    0.025   22.583    0.000    0.566    0.475
#  factor5 =~                                                            
#    ZTPI_11           1.130    0.229    4.946    0.000    1.130    0.978
#    ZTPI_6            0.251    0.057    4.423    0.000    0.251    0.263
#    ZTPI_18           0.175    0.044    3.994    0.000    0.175    0.153
#  general =~                                                            
#    ZTPI_1            0.130    0.026    4.983    0.000    0.130    0.127
#    ZTPI_2           -0.484    0.049   -9.805    0.000   -0.484   -0.393
#    ZTPI_3            0.161    0.026    6.219    0.000    0.161    0.140
#    ZTPI_4           -0.488    0.027  -17.765    0.000   -0.488   -0.421
#    ZTPI_5            0.245    0.023   10.507    0.000    0.245    0.242
#    ZTPI_6           -0.381    0.025  -15.494    0.000   -0.381   -0.400
#    ZTPI_7           -0.215    0.030   -7.152    0.000   -0.215   -0.185
#    ZTPI_8           -0.125    0.027   -4.549    0.000   -0.125   -0.110
#    ZTPI_9           -0.741    0.027  -27.103    0.000   -0.741   -0.608
#    ZTPI_10          -0.708    0.034  -20.587    0.000   -0.708   -0.541
#    ZTPI_11          -0.311    0.025  -12.184    0.000   -0.311   -0.269
#    ZTPI_12           0.324    0.023   14.219    0.000    0.324    0.318
#    ZTPI_13          -0.146    0.025   -5.803    0.000   -0.146   -0.125
#    ZTPI_14           0.225    0.024    9.429    0.000    0.225    0.196
#    ZTPI_15          -0.728    0.031  -23.835    0.000   -0.728   -0.610
#    ZTPI_16          -0.392    0.024  -16.324    0.000   -0.392   -0.329
#    ZTPI_17          -0.706    0.029  -24.749    0.000   -0.706   -0.550
#    ZTPI_18          -0.616    0.023  -26.887    0.000   -0.616   -0.539

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.3502257      0.8235294      0.7415763      0.4012408 

#$FactorLevelIndices
#ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3009853 0.09695594 0.69901475 0.7793552 0.1840673 0.4930731 0.6811290
#factor2 0.8303065 0.13075443 0.16969354 0.6475677 0.5377969 0.6101701 0.7888281
#factor3 0.9406393 0.14992427 0.05936066 0.6471341 0.6458119 0.6387601 0.8043034
#factor4 0.8766339 0.13054871 0.12336607 0.6266068 0.5639866 0.5953206 0.7733375
#factor5 0.6676039 0.14159092 0.33239607 0.7040750 0.4022821 0.9566468 1.0052621
#general 0.3502257 0.35022573 0.35022573 0.7415763 0.4012408 0.7795833 0.8557236



