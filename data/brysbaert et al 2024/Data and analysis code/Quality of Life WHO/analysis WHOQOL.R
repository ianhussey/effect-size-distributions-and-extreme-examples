######################################################################
### Analysis WHO dataset quality of life
### Particularly interesting because the dataset has been
### thoroughly investigated by Rogers (2021, 2024)
### https://www.researchgate.net/publication/356264078_Best_Practices_for_Your_Exploratory_Factor_Analysis_A_Factor_Tutorial
### https://www.researchgate.net/publication/375605644_Best_Practices_for_your_Confirmatory_Factor_Analysis_A_JASP_and_lavaan_Tutorial
### data : https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/RCX8FF


Rogers_WHOQOL_Data <- read.delim("Rogers_WHOQOL_Data.dat", header=FALSE)
colnames(Rogers_WHOQOL_Data)
mydata <- Rogers_WHOQOL_Data
mydata <- na.omit(mydata)
colnames(mydata) <- c('Q3_F','Q4_F','Q10_F','Q15_F','Q16_F','Q17_F','Q18_F','Q5_P',
                      'Q6_P','Q7_P','Q11_P','Q19_P','Q26_P','Q20_S','Q21_S','Q22_S',
                      'Q8_A','Q9_A','Q12_A','Q13_A','Q14_A','Q23_A','Q24_A','Q25_A')
min(mydata)
max(mydata) # 5 response alternatives

# parallel EFA analysis
library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 7.51; eigenvalue 2 = 1.18

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 3 components
# Eigenvalue 1 = 8.62; eigenvalue 2 = 1.33

library(EGAnet)
library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 4 communities



# Give solution with 4 factors
fit3 <- fa(mydata,4)
fit3
diagram(fit3)
# %variance explained=.43, RMSEA=.057, RMSR=.03, TLI=.901
#     MR1  MR4  MR2  MR3
#MR1 1.00 0.48 0.38 0.35
#MR4 0.48 1.00 0.44 0.14
#MR2 0.38 0.44 1.00 0.22
#MR3 0.35 0.14 0.22 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.077, RMSR=.03, TLI=.863
#    MR1  MR4  MR2  MR3
#MR1 1.00 0.54 0.34 0.29
#MR4 0.54 1.00 0.43 0.14
#MR2 0.34 0.43 1.00 0.13
#MR3 0.29 0.14 0.13 1.00

#Compare to solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit4)
# %variance explained=.31, RMSEA=.095, RMSR=.07, TLI=.729

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.117, RMSR=.08, TLI=.686


# Single factor model lavaan
UNImodel= '
 general=~  Q3_F+Q4_F+Q10_F+Q15_F+Q16_F+Q17_F+Q18_F+Q5_P+Q6_P+Q7_P+Q11_P+
            Q19_P+Q26_P+Q20_S+Q21_S+Q22_S+Q8_A+Q9_A+Q12_A+Q13_A+Q14_A+
            Q23_A+Q24_A+Q25_A
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.959       0.863
#Tucker-Lewis Index (TLI)                       0.955       0.850
#Robust Comparative Fit Index (CFI)                         0.715
#Robust Tucker-Lewis Index (TLI)                            0.688
#RMSEA                                          0.105       0.115
#Robust RMSEA                                               0.117
#SRMR                                           0.084       0.084

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .358

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.755       0.756
#Tucker-Lewis Index (TLI)                       0.732       0.732
#Robust Comparative Fit Index (CFI)                         0.759
#Robust Tucker-Lewis Index (TLI)                            0.737
#RMSEA                                          0.095       0.084
#Robust RMSEA                                               0.094
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .277

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model based on theory
EGAmodel= '
 physical  =~ Q3_F+Q4_F+Q10_F+Q15_F+Q16_F+Q17_F+Q18_F
 psycholog =~ Q5_P+Q6_P+Q7_P+Q11_P+Q19_P+Q26_P
 social    =~ Q20_S+Q21_S+Q22_S
 environm  =~ Q8_A+Q9_A+Q12_A+Q13_A+Q14_A+Q23_A+Q24_A+Q25_A
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.981       0.929
#Tucker-Lewis Index (TLI)                       0.979       0.921
#Robust Comparative Fit Index (CFI)                         0.838
#Robust Tucker-Lewis Index (TLI)                            0.819
#RMSEA                                          0.073       0.084
#Robust RMSEA                                               0.090
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .428

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  physical ~~                                                           
#    psycholog         0.902    0.010   89.426    0.000    0.902    0.902
#    social            0.631    0.024   26.288    0.000    0.631    0.631
#    environm          0.624    0.022   28.548    0.000    0.624    0.624
#  psycholog ~~                                                          
#    social            0.771    0.019   41.007    0.000    0.771    0.771
#    environm          0.725    0.021   34.935    0.000    0.725    0.725
#  social ~~                                                             
#    environm          0.550    0.028   19.550    0.000    0.550    0.550


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.865       0.867
#Tucker-Lewis Index (TLI)                       0.849       0.850
#Robust Comparative Fit Index (CFI)                         0.870
#Robust Tucker-Lewis Index (TLI)                            0.854
#RMSEA                                          0.071       0.063
#Robust RMSEA                                               0.070
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .367

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  physical ~~                                                           
#    psycholog         0.897    0.017   52.797    0.000    0.897    0.897
#    social            0.618    0.032   19.246    0.000    0.618    0.618
#    environm          0.586    0.032   18.284    0.000    0.586    0.586
#  psycholog ~~                                                          
#    social            0.762    0.026   29.492    0.000    0.762    0.762
#    environm          0.703    0.031   22.669    0.000    0.703    0.703
#  social ~~                                                             
#    environm          0.552    0.037   14.878    0.000    0.552    0.552

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.674       0.668
#Tucker-Lewis Index (TLI)                       0.643       0.637
#Robust Comparative Fit Index (CFI)                         0.677
#Robust Tucker-Lewis Index (TLI)                            0.647
#RMSEA                                          0.109       0.098
#Robust RMSEA                                               0.108
#SRMR                                           0.247       0.247

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .370


# Bifactor model
BIFmodel= '
 physical  =~ Q3_F+Q4_F+Q10_F+Q15_F+Q16_F+Q17_F+Q18_F
 psycholog =~ Q5_P+Q6_P+Q7_P+Q11_P+Q19_P+Q26_P
 social    =~ Q20_S+Q21_S+Q22_S
 environm  =~ Q8_A+Q9_A+Q12_A+Q13_A+Q14_A+Q23_A+Q24_A+Q25_A
 general=~  Q3_F+Q4_F+Q10_F+Q15_F+Q16_F+Q17_F+Q18_F+Q5_P+Q6_P+Q7_P+Q11_P+
            Q19_P+Q26_P+Q20_S+Q21_S+Q22_S+Q8_A+Q9_A+Q12_A+Q13_A+Q14_A+
            Q23_A+Q24_A+Q25_A
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.894
#Tucker-Lewis Index (TLI)                       0.872       0.872
#Robust Comparative Fit Index (CFI)                         0.899
#Robust Tucker-Lewis Index (TLI)                            0.877
#RMSEA                                          0.065       0.058
#Robust RMSEA                                               0.064
#SRMR                                           0.051       0.051

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .401

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6463732      0.7572464      0.9263957      0.8200990 
#
#$FactorLevelIndices
#             ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#physical  0.2608112 0.08989877 0.7391888 0.8606759 0.14205179 0.6459207 0.8370651
#psycholog 0.1904942 0.05254659 0.8095058 0.8380079 0.09996527 0.4666101 0.7747746
#social    0.4996438 0.06579307 0.5003562 0.7135450 0.34185610 0.5004686 0.7650618
#environm  0.5867463 0.14538833 0.4132537 0.7916222 0.47002707 0.6620526 0.8223629
#general   0.6463732 0.64637324 0.6463732 0.9263957 0.82009903 0.9304690 0.9568642





### Analysis group  gambling

mydata <- mydata2

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 60.08; eigenvalue 2 = 1.85

rho <- polychoric(mydata)$rho
# Warning message:
# In cor.smooth(mat) : Matrix was not positive definite, smoothing was done
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 65.72; eigenvalue 2 = 1.58

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.67, RMSEA=.066, RMSR=.04, TLI=.832

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.73, RMSEA=.189, RMSR=.03, TLI=.418

# Alternative for psych smoothing
# library(EGAnet)
# rho2 <- polychoric.matrix(mydata)
# library(fungible)
# rho <- corSmooth(rho2, eps = 1E8 * .Machine$double.eps)
# fit4 <- fa(rho,1,n.obs=nrow(mydata))
# fit4
# diagram(fit4)
# %variance explained=.66, RMSEA=.173, RMSR=.05, TLI=..421

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 4 communities

# Give solution with 2 factors (based on overall analysis)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.66, RMSEA=.061, RMSR=.03, TLI=.857
#     MR1  MR2
#MR1 1.00 0.19
#MR2 0.19 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.75, RMSEA=.189, RMSR=.03, TLI=.417
#     MR1  MR2
#MR1 1.00 0.20
#MR2 0.20 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  COMP1_1+COMP1_2+COMP1_3+COMP1_4+COMP1_5+COMP1_6+COMP1_7+COMP1_8+COMP1_9+
            COMP1_10+COMP1_11+COMP1_12+COMP1_13+COMP1_14+COMP1_15+COMP1_16+COMP1_17+COMP1_18+
            COMP1_19+COMP1_20+COMP1_21+COMP1_22+COMP1_23+COMP1_24+COMP1_25+COMP1_26+COMP1_27+
            COMP1_28+COMP1_29+COMP1_30+COMP1_31+COMP1_32+COMP1_33+COMP1_34+COMP1_35+COMP1_36+
            COMP1_37+COMP1_38+COMP1_39+COMP1_40+COMP1_41+COMP1_42+COMP1_43+COMP1_44+COMP1_45+
            COMP1_46+COMP1_47+COMP1_48+COMP1_49+COMP1_50+COMP1_51+COMP1_52+COMP1_53+COMP1_54+
            COMP1_55+COMP1_56+COMP1_57+COMP1_58+COMP1_59+COMP1_60+COMP1_61+COMP1_62+COMP1_63+
            COMP1_64+COMP1_65+COMP1_66+COMP1_67+COMP1_68+COMP1_69+COMP1_70+COMP1_71+COMP1_72+
            COMP1_73+COMP1_74+COMP1_75+COMP1_76+COMP1_77+COMP1_78+COMP1_79+COMP1_80+COMP1_81+
            COMP1_82+COMP1_83+COMP1_84+COMP1_85+COMP1_86+COMP1_87+COMP1_88+COMP1_89+COMP1_90
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# Warning message: The variance-covariance matrix of the estimated parameters (vcov)
# does not appear to be positive definite!
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.984
#Tucker-Lewis Index (TLI)                       0.999       0.984
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.039       0.044
#Robust RMSEA                                                  NA
#SRMR                                           0.035       0.035

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .772


CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.826       0.847
#Tucker-Lewis Index (TLI)                       0.822       0.843
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.846
#RMSEA                                          0.073       0.059
#Robust RMSEA                                               0.066
#SRMR                                           0.037       0.037

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .701

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 2 factors based on the joint analysis of the two studies
EGAmodel= '
 factor1 =~ COMP1_43+COMP1_45+COMP1_51+COMP1_78+COMP1_79+COMP1_40+COMP1_47+COMP1_74+
            COMP1_56+COMP1_14+COMP1_16+COMP1_90+COMP1_44+COMP1_57+COMP1_80+COMP1_12+
            COMP1_55+COMP1_73+COMP1_46+COMP1_50+COMP1_61+COMP1_82+COMP1_49+COMP1_67+
            COMP1_24+COMP1_38+COMP1_83+COMP1_64+COMP1_26+COMP1_52+COMP1_54+COMP1_76+
            COMP1_9+COMP1_60+COMP1_13+COMP1_62+COMP1_71+COMP1_20+COMP1_34+COMP1_39+
            COMP1_48+COMP1_15+COMP1_25+COMP1_18+COMP1_84+COMP1_85+COMP1_35+COMP1_10+
            COMP1_86+COMP1_28+COMP1_58+COMP1_23+COMP1_89+COMP1_66+COMP1_68+COMP1_8+
            COMP1_32+COMP1_36+COMP1_65+COMP1_77+COMP1_88+COMP1_7+COMP1_19+COMP1_72+
            COMP1_37+COMP1_42+COMP1_75+COMP1_6+COMP1_63+COMP1_21
 factor2 =~ COMP1_17+COMP1_53+COMP1_87+COMP1_22+COMP1_70+COMP1_81+COMP1_11+COMP1_29+
            COMP1_30+COMP1_59+COMP1_41+COMP1_1+COMP1_69+COMP1_27+COMP1_2+COMP1_33+
            COMP1_31+COMP1_5+COMP1_4+COMP1_3
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning The variance-covariance matrix of the estimated parameters (vcov)
# does not appear to be positive definite!
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.985
#Tucker-Lewis Index (TLI)                       1.000       0.985
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.035       0.043
#Robust RMSEA                                                  NA
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .783

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.963    0.004  220.931    0.000    0.963    0.963


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.838       0.859
#Tucker-Lewis Index (TLI)                       0.834       0.856
#Robust Comparative Fit Index (CFI)                         0.862
#Robust Tucker-Lewis Index (TLI)                            0.859
#RMSEA                                          0.070       0.057
#Robust RMSEA                                               0.064
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .718

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.947    0.008  112.902    0.000    0.947    0.947

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.821       0.841
#Tucker-Lewis Index (TLI)                       0.817       0.837
#Robust Comparative Fit Index (CFI)                         0.844
#Robust Tucker-Lewis Index (TLI)                            0.841
#RMSEA                                          0.074       0.060
#Robust RMSEA                                               0.068
#SRMR                                           0.362       0.362

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .716


# Bifactor model
BIFmodel= '
 factor1 =~ COMP1_43+COMP1_45+COMP1_51+COMP1_78+COMP1_79+COMP1_40+COMP1_47+COMP1_74+
            COMP1_56+COMP1_14+COMP1_16+COMP1_90+COMP1_44+COMP1_57+COMP1_80+COMP1_12+
            COMP1_55+COMP1_73+COMP1_46+COMP1_50+COMP1_61+COMP1_82+COMP1_49+COMP1_67+
            COMP1_24+COMP1_38+COMP1_83+COMP1_64+COMP1_26+COMP1_52+COMP1_54+COMP1_76+
            COMP1_9+COMP1_60+COMP1_13+COMP1_62+COMP1_71+COMP1_20+COMP1_34+COMP1_39+
            COMP1_48+COMP1_15+COMP1_25+COMP1_18+COMP1_84+COMP1_85+COMP1_35+COMP1_10+
            COMP1_86+COMP1_28+COMP1_58+COMP1_23+COMP1_89+COMP1_66+COMP1_68+COMP1_8+
            COMP1_32+COMP1_36+COMP1_65+COMP1_77+COMP1_88+COMP1_7+COMP1_19+COMP1_72+
            COMP1_37+COMP1_42+COMP1_75+COMP1_6+COMP1_63+COMP1_21
 factor2 =~ COMP1_17+COMP1_53+COMP1_87+COMP1_22+COMP1_70+COMP1_81+COMP1_11+COMP1_29+
            COMP1_30+COMP1_59+COMP1_41+COMP1_1+COMP1_69+COMP1_27+COMP1_2+COMP1_33+
            COMP1_31+COMP1_5+COMP1_4+COMP1_3
 general=~  COMP1_1+COMP1_2+COMP1_3+COMP1_4+COMP1_5+COMP1_6+COMP1_7+COMP1_8+COMP1_9+
            COMP1_10+COMP1_11+COMP1_12+COMP1_13+COMP1_14+COMP1_15+COMP1_16+COMP1_17+COMP1_18+
            COMP1_19+COMP1_20+COMP1_21+COMP1_22+COMP1_23+COMP1_24+COMP1_25+COMP1_26+COMP1_27+
            COMP1_28+COMP1_29+COMP1_30+COMP1_31+COMP1_32+COMP1_33+COMP1_34+COMP1_35+COMP1_36+
            COMP1_37+COMP1_38+COMP1_39+COMP1_40+COMP1_41+COMP1_42+COMP1_43+COMP1_44+COMP1_45+
            COMP1_46+COMP1_47+COMP1_48+COMP1_49+COMP1_50+COMP1_51+COMP1_52+COMP1_53+COMP1_54+
            COMP1_55+COMP1_56+COMP1_57+COMP1_58+COMP1_59+COMP1_60+COMP1_61+COMP1_62+COMP1_63+
            COMP1_64+COMP1_65+COMP1_66+COMP1_67+COMP1_68+COMP1_69+COMP1_70+COMP1_71+COMP1_72+
            COMP1_73+COMP1_74+COMP1_75+COMP1_76+COMP1_77+COMP1_78+COMP1_79+COMP1_80+COMP1_81+
            COMP1_82+COMP1_83+COMP1_84+COMP1_85+COMP1_86+COMP1_87+COMP1_88+COMP1_89+COMP1_90
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.859       0.871
#Tucker-Lewis Index (TLI)                       0.852       0.865
#Robust Comparative Fit Index (CFI)                         0.878
#Robust Tucker-Lewis Index (TLI)                            0.872
#RMSEA                                          0.066       0.055
#Robust RMSEA                                               0.060
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .728

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.9020227      0.3495630      0.9948648      0.9413746 
#
#$FactorLevelIndices
#            ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.09598964 0.07662666 0.9040104 0.9941543 0.07836670 0.8452991 0.9179902
#factor2 0.10584342 0.02135067 0.8941566 0.9699082 0.06890996 0.6099802 0.8632205
#general 0.90202266 0.90202266 0.9020227 0.9948648 0.94137464 0.9941671 0.9919201






######################################################################
### Analysis WHO dataset quality of life
### 
### McConachie et al. (2018) https://psycnet.apa.org/record/2017-53987-001
### 

library(readxl)
WHOQOL_McConachie <- read_excel("WHOQOL_McConachie.xlsx")
colnames(WHOQOL_McConachie)
mydata <- WHOQOL_McConachie[,3:26]
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

# parallel EFA analysis
library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 7.26; eigenvalue 2 = 1.01

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 4 components
# Eigenvalue 1 = 8.06; eigenvalue 2 = 1.09

library(EGAnet)
library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 4 communities


#Compute to solution with 1 factor
fit1 <- fa(mydata,1)
fit1
diagram(fit1)
# %variance explained=.30, RMSEA=.097, RMSR=.08, TLI=.705

fit2 <- fa(rho,1,n.obs=nrow(mydata))
fit2
diagram(fit2)
# %variance explained=.34, RMSEA=.118, RMSR=.09, TLI=.653

# Give solution with 4 factors
fit3 <- fa(mydata,4)
fit3
diagram(fit3)
# %variance explained=.44, RMSEA=.066, RMSR=.04, TLI=.859
#     MR1  MR4  MR3  MR2
#MR1 1.00 0.66 0.23 0.03
#MR4 0.66 1.00 0.32 0.02
#MR3 0.23 0.32 1.00 0.04
#MR2 0.03 0.02 0.04 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.09, RMSR=.04, TLI=.796
#     MR1  MR4  MR3  MR2
#MR1 1.00 0.66 0.25 0.04
#MR4 0.66 1.00 0.33 0.02
#MR3 0.25 0.33 1.00 0.04
#MR2 0.04 0.02 0.04 1.00

# Single factor model lavaan
UNImodel= '
 general=~  Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q11+Q12+Q13+
            Q14+Q15+Q16+Q17+Q18+Q19+Q20+Q21+Q22+Q23+Q24+Q25+Q26
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.942       0.820
#Tucker-Lewis Index (TLI)                       0.936       0.802
#Robust Comparative Fit Index (CFI)                         0.695
#Robust Tucker-Lewis Index (TLI)                            0.666
#RMSEA                                          0.110       0.114
#Robust RMSEA                                               0.119
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .382

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.730       0.733
#Tucker-Lewis Index (TLI)                       0.704       0.707
#Robust Comparative Fit Index (CFI)                         0.735
#Robust Tucker-Lewis Index (TLI)                            0.710
#RMSEA                                          0.098       0.092
#Robust RMSEA                                               0.097
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .327

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model based on theory
EGAmodel= '
 physical  =~ Q3+Q4+Q10+Q15+Q16+Q17+Q18
 psycholog =~ Q5+Q6+Q7+Q11+Q19+Q26
 social    =~ Q20+Q21+Q22
 environm  =~ Q8+Q9+Q12+Q13+Q14+Q23+Q24+Q25
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.967       0.886
#Tucker-Lewis Index (TLI)                       0.963       0.872
#Robust Comparative Fit Index (CFI)                         0.800
#Robust Tucker-Lewis Index (TLI)                            0.776
#RMSEA                                          0.084       0.092
#Robust RMSEA                                               0.097
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .479

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  physical ~~                                                           
#    psycholog         0.784    0.028   28.034    0.000    0.784    0.784
#    social            0.412    0.054    7.628    0.000    0.412    0.412
#    environm          0.824    0.027   30.033    0.000    0.824    0.824
#  psycholog ~~                                                          
#    social            0.602    0.046   13.182    0.000    0.602    0.602
#    environm          0.700    0.033   20.947    0.000    0.700    0.700
#  social ~~                                                             
#    environm          0.520    0.046   11.252    0.000    0.520    0.520

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.845       0.847
#Tucker-Lewis Index (TLI)                       0.827       0.829
#Robust Comparative Fit Index (CFI)                         0.851
#Robust Tucker-Lewis Index (TLI)                            0.832
#RMSEA                                          0.075       0.071
#Robust RMSEA                                               0.074
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .406

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  physical ~~                                                           
#    psycholog         0.758    0.040   18.852    0.000    0.758    0.758
#    social            0.308    0.070    4.411    0.000    0.308    0.308
#    environm          0.795    0.038   21.038    0.000    0.795    0.795
#  psycholog ~~                                                          
#    social            0.543    0.064    8.438    0.000    0.543    0.543
#    environm          0.689    0.043   16.071    0.000    0.689    0.689
#  social ~~                                                             
#    environm          0.409    0.066    6.214    0.000    0.409    0.409

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
# warning
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.685       0.684
#Tucker-Lewis Index (TLI)                       0.655       0.654
#Robust Comparative Fit Index (CFI)                         0.689
#Robust Tucker-Lewis Index (TLI)                            0.660
#RMSEA                                          0.106       0.101
#Robust RMSEA                                               0.105
#SRMR                                           0.238       0.238

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .398


# Bifactor model
BIFmodel= '
 physical  =~ Q3+Q4+Q10+Q15+Q16+Q17+Q18
 psycholog =~ Q5+Q6+Q7+Q11+Q19+Q26
 social    =~ Q20+Q21+Q22
 environm  =~ Q8+Q9+Q12+Q13+Q14+Q23+Q24+Q25
 general=~  Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q11+Q12+Q13+
            Q14+Q15+Q16+Q17+Q18+Q19+Q20+Q21+Q22+Q23+Q24+Q25+Q26
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.898       0.899
#Tucker-Lewis Index (TLI)                       0.876       0.878
#Robust Comparative Fit Index (CFI)                         0.902
#Robust Tucker-Lewis Index (TLI)                            0.882
#RMSEA                                          0.064       0.060
#Robust RMSEA                                               0.062
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .418

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5948573      0.7572464      0.9196447      0.8146909 
#
#$FactorLevelIndices
#             ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#physical  0.3412863 0.11938134 0.6587137 0.8331159 0.1244052 0.9734921 0.9867907
#psycholog 0.3670833 0.07812090 0.6329167 0.7318010 0.2243069 0.5496815 0.7908510
#social    0.7955083 0.10953771 0.2044917 0.7229988 0.5593244 0.8839808 0.9948550
#environm  0.3273459 0.09810274 0.6726541 0.8501371 0.2728309 0.5641005 0.7591919
#general   0.5948573 0.59485731 0.5948573 0.9196447 0.8146909 0.9220038 0.9501710







### Analysis group  gambling

mydata <- mydata2

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 60.08; eigenvalue 2 = 1.85

rho <- polychoric(mydata)$rho
# Warning message:
# In cor.smooth(mat) : Matrix was not positive definite, smoothing was done
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 65.72; eigenvalue 2 = 1.58

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.67, RMSEA=.066, RMSR=.04, TLI=.832

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.73, RMSEA=.189, RMSR=.03, TLI=.418

# Alternative for psych smoothing
# library(EGAnet)
# rho2 <- polychoric.matrix(mydata)
# library(fungible)
# rho <- corSmooth(rho2, eps = 1E8 * .Machine$double.eps)
# fit4 <- fa(rho,1,n.obs=nrow(mydata))
# fit4
# diagram(fit4)
# %variance explained=.66, RMSEA=.173, RMSR=.05, TLI=..421

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 4 communities

# Give solution with 2 factors (based on overall analysis)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.66, RMSEA=.061, RMSR=.03, TLI=.857
#     MR1  MR2
#MR1 1.00 0.19
#MR2 0.19 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.75, RMSEA=.189, RMSR=.03, TLI=.417
#     MR1  MR2
#MR1 1.00 0.20
#MR2 0.20 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  COMP1_1+COMP1_2+COMP1_3+COMP1_4+COMP1_5+COMP1_6+COMP1_7+COMP1_8+COMP1_9+
            COMP1_10+COMP1_11+COMP1_12+COMP1_13+COMP1_14+COMP1_15+COMP1_16+COMP1_17+COMP1_18+
            COMP1_19+COMP1_20+COMP1_21+COMP1_22+COMP1_23+COMP1_24+COMP1_25+COMP1_26+COMP1_27+
            COMP1_28+COMP1_29+COMP1_30+COMP1_31+COMP1_32+COMP1_33+COMP1_34+COMP1_35+COMP1_36+
            COMP1_37+COMP1_38+COMP1_39+COMP1_40+COMP1_41+COMP1_42+COMP1_43+COMP1_44+COMP1_45+
            COMP1_46+COMP1_47+COMP1_48+COMP1_49+COMP1_50+COMP1_51+COMP1_52+COMP1_53+COMP1_54+
            COMP1_55+COMP1_56+COMP1_57+COMP1_58+COMP1_59+COMP1_60+COMP1_61+COMP1_62+COMP1_63+
            COMP1_64+COMP1_65+COMP1_66+COMP1_67+COMP1_68+COMP1_69+COMP1_70+COMP1_71+COMP1_72+
            COMP1_73+COMP1_74+COMP1_75+COMP1_76+COMP1_77+COMP1_78+COMP1_79+COMP1_80+COMP1_81+
            COMP1_82+COMP1_83+COMP1_84+COMP1_85+COMP1_86+COMP1_87+COMP1_88+COMP1_89+COMP1_90
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# Warning message: The variance-covariance matrix of the estimated parameters (vcov)
# does not appear to be positive definite!
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.984
#Tucker-Lewis Index (TLI)                       0.999       0.984
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.039       0.044
#Robust RMSEA                                                  NA
#SRMR                                           0.035       0.035

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .772


CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.826       0.847
#Tucker-Lewis Index (TLI)                       0.822       0.843
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.846
#RMSEA                                          0.073       0.059
#Robust RMSEA                                               0.066
#SRMR                                           0.037       0.037

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .701

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 2 factors based on the joint analysis of the two studies
EGAmodel= '
 factor1 =~ COMP1_43+COMP1_45+COMP1_51+COMP1_78+COMP1_79+COMP1_40+COMP1_47+COMP1_74+
            COMP1_56+COMP1_14+COMP1_16+COMP1_90+COMP1_44+COMP1_57+COMP1_80+COMP1_12+
            COMP1_55+COMP1_73+COMP1_46+COMP1_50+COMP1_61+COMP1_82+COMP1_49+COMP1_67+
            COMP1_24+COMP1_38+COMP1_83+COMP1_64+COMP1_26+COMP1_52+COMP1_54+COMP1_76+
            COMP1_9+COMP1_60+COMP1_13+COMP1_62+COMP1_71+COMP1_20+COMP1_34+COMP1_39+
            COMP1_48+COMP1_15+COMP1_25+COMP1_18+COMP1_84+COMP1_85+COMP1_35+COMP1_10+
            COMP1_86+COMP1_28+COMP1_58+COMP1_23+COMP1_89+COMP1_66+COMP1_68+COMP1_8+
            COMP1_32+COMP1_36+COMP1_65+COMP1_77+COMP1_88+COMP1_7+COMP1_19+COMP1_72+
            COMP1_37+COMP1_42+COMP1_75+COMP1_6+COMP1_63+COMP1_21
 factor2 =~ COMP1_17+COMP1_53+COMP1_87+COMP1_22+COMP1_70+COMP1_81+COMP1_11+COMP1_29+
            COMP1_30+COMP1_59+COMP1_41+COMP1_1+COMP1_69+COMP1_27+COMP1_2+COMP1_33+
            COMP1_31+COMP1_5+COMP1_4+COMP1_3
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning The variance-covariance matrix of the estimated parameters (vcov)
# does not appear to be positive definite!
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.985
#Tucker-Lewis Index (TLI)                       1.000       0.985
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.035       0.043
#Robust RMSEA                                                  NA
#SRMR                                           0.034       0.034

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .783

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.963    0.004  220.931    0.000    0.963    0.963


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.838       0.859
#Tucker-Lewis Index (TLI)                       0.834       0.856
#Robust Comparative Fit Index (CFI)                         0.862
#Robust Tucker-Lewis Index (TLI)                            0.859
#RMSEA                                          0.070       0.057
#Robust RMSEA                                               0.064
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .718

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.947    0.008  112.902    0.000    0.947    0.947

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.821       0.841
#Tucker-Lewis Index (TLI)                       0.817       0.837
#Robust Comparative Fit Index (CFI)                         0.844
#Robust Tucker-Lewis Index (TLI)                            0.841
#RMSEA                                          0.074       0.060
#Robust RMSEA                                               0.068
#SRMR                                           0.362       0.362

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .716


# Bifactor model
BIFmodel= '
 factor1 =~ COMP1_43+COMP1_45+COMP1_51+COMP1_78+COMP1_79+COMP1_40+COMP1_47+COMP1_74+
            COMP1_56+COMP1_14+COMP1_16+COMP1_90+COMP1_44+COMP1_57+COMP1_80+COMP1_12+
            COMP1_55+COMP1_73+COMP1_46+COMP1_50+COMP1_61+COMP1_82+COMP1_49+COMP1_67+
            COMP1_24+COMP1_38+COMP1_83+COMP1_64+COMP1_26+COMP1_52+COMP1_54+COMP1_76+
            COMP1_9+COMP1_60+COMP1_13+COMP1_62+COMP1_71+COMP1_20+COMP1_34+COMP1_39+
            COMP1_48+COMP1_15+COMP1_25+COMP1_18+COMP1_84+COMP1_85+COMP1_35+COMP1_10+
            COMP1_86+COMP1_28+COMP1_58+COMP1_23+COMP1_89+COMP1_66+COMP1_68+COMP1_8+
            COMP1_32+COMP1_36+COMP1_65+COMP1_77+COMP1_88+COMP1_7+COMP1_19+COMP1_72+
            COMP1_37+COMP1_42+COMP1_75+COMP1_6+COMP1_63+COMP1_21
 factor2 =~ COMP1_17+COMP1_53+COMP1_87+COMP1_22+COMP1_70+COMP1_81+COMP1_11+COMP1_29+
            COMP1_30+COMP1_59+COMP1_41+COMP1_1+COMP1_69+COMP1_27+COMP1_2+COMP1_33+
            COMP1_31+COMP1_5+COMP1_4+COMP1_3
 general=~  COMP1_1+COMP1_2+COMP1_3+COMP1_4+COMP1_5+COMP1_6+COMP1_7+COMP1_8+COMP1_9+
            COMP1_10+COMP1_11+COMP1_12+COMP1_13+COMP1_14+COMP1_15+COMP1_16+COMP1_17+COMP1_18+
            COMP1_19+COMP1_20+COMP1_21+COMP1_22+COMP1_23+COMP1_24+COMP1_25+COMP1_26+COMP1_27+
            COMP1_28+COMP1_29+COMP1_30+COMP1_31+COMP1_32+COMP1_33+COMP1_34+COMP1_35+COMP1_36+
            COMP1_37+COMP1_38+COMP1_39+COMP1_40+COMP1_41+COMP1_42+COMP1_43+COMP1_44+COMP1_45+
            COMP1_46+COMP1_47+COMP1_48+COMP1_49+COMP1_50+COMP1_51+COMP1_52+COMP1_53+COMP1_54+
            COMP1_55+COMP1_56+COMP1_57+COMP1_58+COMP1_59+COMP1_60+COMP1_61+COMP1_62+COMP1_63+
            COMP1_64+COMP1_65+COMP1_66+COMP1_67+COMP1_68+COMP1_69+COMP1_70+COMP1_71+COMP1_72+
            COMP1_73+COMP1_74+COMP1_75+COMP1_76+COMP1_77+COMP1_78+COMP1_79+COMP1_80+COMP1_81+
            COMP1_82+COMP1_83+COMP1_84+COMP1_85+COMP1_86+COMP1_87+COMP1_88+COMP1_89+COMP1_90
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.859       0.871
#Tucker-Lewis Index (TLI)                       0.852       0.865
#Robust Comparative Fit Index (CFI)                         0.878
#Robust Tucker-Lewis Index (TLI)                            0.872
#RMSEA                                          0.066       0.055
#Robust RMSEA                                               0.060
#SRMR                                           0.036       0.036

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .728

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.9020227      0.3495630      0.9948648      0.9413746 
#
#$FactorLevelIndices
#            ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.09598964 0.07662666 0.9040104 0.9941543 0.07836670 0.8452991 0.9179902
#factor2 0.10584342 0.02135067 0.8941566 0.9699082 0.06890996 0.6099802 0.8632205
#general 0.90202266 0.90202266 0.9020227 0.9948648 0.94137464 0.9941671 0.9919201






