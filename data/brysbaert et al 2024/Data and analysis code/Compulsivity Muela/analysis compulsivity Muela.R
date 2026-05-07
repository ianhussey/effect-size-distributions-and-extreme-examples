######################################################################
### Analysis symptom checklist to test the usefulness of the criteria
###

################################################################
### Article on compulsivity in video gaming and gambling
### Muela et al (2023) https://link.springer.com/article/10.1186/s40359-023-01439-1
### Data https://osf.io/xdfmw/?view_only=9831ce7702c34347ac67b45719ddf643


### First combine the two datasets to get an overall EFA model to be used in lavaan

### Video gaming

library(readxl)
compulsivity_videogames_Muela <- read_excel("compulsivity_videogames_Muela.xlsx")
colnames(compulsivity_videogames_Muela)
mydata1 <- compulsivity_videogames_Muela[,10:99]
colnames(mydata1)
mydata1 <- na.omit(mydata1)
min(mydata1)
max(mydata1) # 5 response alternatives

compulsivity_gambling_Muela <- read_excel("compulsivity_gambling_Muela.xlsx")
colnames(compulsivity_gambling_Muela)
mydata2 <- compulsivity_gambling_Muela[,10:99]
colnames(mydata2)
mydata2 <- na.omit(mydata2)
min(mydata2)
max(mydata2) # 5 response alternatives

mydata <- rbind(mydata1,mydata2)

library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 56.41; eigenvalue 2 = 2.74

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 62.57; eigenvalue 2 = 2.34

# Give solution with 3 factors
fit3 <- fa(mydata,3)
fit3
diagram(fit3)
# last factor has only three items; so reduce to 2
fit3 <- fa(mydata,2)
fit3
diagram(fit3)
# first factor has 70 items second 20


### Analysis group video gaming

mydata <- mydata1

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 52.96; eigenvalue 2 = 3.61

rho <- polychoric(mydata)$rho
# Warning message:
# In cor.smooth(mat) : Matrix was not positive definite, smoothing was done
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 59.49; eigenvalue 2 = 3.25

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.59, RMSEA=.072, RMSR=.06, TLI=.772

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.66, RMSEA=.193, RMSR=.05, TLI=.366

# Alternative for psych smoothing
# library(EGAnet)
# rho2 <- polychoric.matrix(mydata)
# library(fungible)
# rho <- corSmooth(rho2, eps = 1E8 * .Machine$double.eps)
# fit4 <- fa(rho,1,n.obs=nrow(mydata))
# fit4
# diagram(fit4)
# %variance explained=.66, RMSEA=.173, RMSR=.05, TLI=..421

library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 6 communities

# Give solution with 2 factors (based on overall analysis)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.062, RMSR=.04, TLI=.832
#     MR1  MR2
#MR1 1.00 0.68
#MR2 0.68 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.70, RMSEA=.191, RMSR=.04, TLI=.378
#     MR1  MR2
#MR1 1.00 0.56
#MR2 0.56 1.00

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
#Comparative Fit Index (CFI)                    0.997       0.968
#Tucker-Lewis Index (TLI)                       0.997       0.967
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.068       0.054
#Robust RMSEA                                                  NA
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .714


CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.768       0.786
#Tucker-Lewis Index (TLI)                       0.763       0.781
#Robust Comparative Fit Index (CFI)                         0.789
#Robust Tucker-Lewis Index (TLI)                            0.784
#RMSEA                                          0.079       0.066
#Robust RMSEA                                               0.074
#SRMR                                           0.056       0.056

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .620

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
#Comparative Fit Index (CFI)                    0.998       0.973
#Tucker-Lewis Index (TLI)                       0.998       0.973
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.058       0.049
#Robust RMSEA                                                  NA
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .730

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.893    0.011   78.482    0.000    0.893    0.893


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.799       0.818
#Tucker-Lewis Index (TLI)                       0.794       0.813
#Robust Comparative Fit Index (CFI)                         0.821
#Robust Tucker-Lewis Index (TLI)                            0.817
#RMSEA                                          0.073       0.061
#Robust RMSEA                                               0.068
#SRMR                                           0.053       0.053

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .636

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#        factor2           0.852    0.018   47.862    0.000    0.852    0.852

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.788       0.806
#Tucker-Lewis Index (TLI)                       0.783       0.801
#Robust Comparative Fit Index (CFI)                         0.809
#Robust Tucker-Lewis Index (TLI)                            0.805
#RMSEA                                          0.075       0.063
#Robust RMSEA                                               0.070
#SRMR                                           0.291       0.291

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .632


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
#Comparative Fit Index (CFI)                    0.837       0.850
#Tucker-Lewis Index (TLI)                       0.829       0.843
#Robust Comparative Fit Index (CFI)                         0.856
#Robust Tucker-Lewis Index (TLI)                            0.849
#RMSEA                                          0.067       0.056
#Robust RMSEA                                               0.062
#SRMR                                           0.040       0.040

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .662

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7838740      0.3495630      0.9931306      0.8494879 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.2311382 0.18618996 0.7688618 0.9925390 0.20842387 0.9328535 0.9538172
#factor2 0.1539409 0.02993609 0.8460591 0.9579432 0.08736435 0.7030314 0.8775085
#general 0.7838740 0.78387395 0.7838740 0.9931306 0.84948794 0.9902131 0.9867312





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



