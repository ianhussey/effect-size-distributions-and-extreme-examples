################################################################
### Analysis Questionnaires for Obsessive-compulsive disorders


################################################################
### Ozcanli et al 2020 https://psychologicabelgica.com/articles/10.5334/pb.537#data-accessibility-statements
library(haven)
Obsession_Ozcanli <- read_sav("Obsession_Ozcanli.sav")
colnames(Obsession_Ozcanli)
mydata  <- as.data.frame(Obsession_Ozcanli[,2:51])
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 components
# Eigenvalue 1 = 14.80; eigenvalue 2 = 3.5

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 4 components
# Eigenvalue 1 = 19.02; eigenvalue 2 = 3.84

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.30, RMSEA=.091, RMSR=.09, TLI=.577

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.133, RMSR=.10, TLI=.476

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities but one very small (4 questions)

# Give solution with 3 factors 4th factor only had 2 significant items
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.064, RMSR=.05, TLI=.79
#     MR1  MR3  MR2
#MR1 1.00 0.45 0.36
#MR3 0.45 1.00 0.43
#MR2 0.36 0.43 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.111, RMSR=.05, TLI=.635
#     MR1  MR2  MR3
#MR1 1.00 0.39 0.45
#MR2 0.39 1.00 0.40
#MR3 0.45 0.40 1.00

# Single factor model lavaan
UNImodel= '
 general=~ LOI_F1+LOI_F2+LOI_F3+LOI_F4+LOI_F5+LOI_F6+LOI_F7+LOI_F8+LOI_F9+LOI_F10+LOI_F11+LOI_F12+LOI_F13+LOI_F14+
   LOI_F15+LOI_F16+LOI_F17+LOI_F18+LOI_F19+LOI_F20+LOI_F21+LOI_F22+LOI_F23+LOI_F24+LOI_F25+LOI_F26+LOI_F27+LOI_F28+
   LOI_F29+LOI_F30+LOI_F31+LOI_F32+LOI_F33+LOI_F34+LOI_F35+LOI_F36+LOI_F37+LOI_F38+LOI_F39+LOI_F40+LOI_F41+LOI_F42+
   LOI_F43+LOI_F44+LOI_F45+LOI_F46+LOI_F47+LOI_F48+LOI_F49+LOI_F50
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.936       0.819
#Tucker-Lewis Index (TLI)                       0.934       0.811
#Robust Comparative Fit Index (CFI)                         0.558
#Robust Tucker-Lewis Index (TLI)                            0.539
#RMSEA                                          0.103       0.081
#Robust RMSEA                                               0.126
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .42

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on EFA analysis
EGAmodel= '
 factor1=~ LOI_F36+LOI_F26+LOI_F12+LOI_F44+LOI_F9+LOI_F43+LOI_F22+LOI_F47+LOI_F31+LOI_F48+LOI_F14+
   LOI_F33+LOI_F19+LOI_F29+LOI_F38+LOI_F6+LOI_F17+LOI_F11+LOI_F23+LOI_F46+LOI_F35+LOI_F41
 factor2=~ LOI_F30+LOI_F40+LOI_F49+LOI_F50+LOI_F37+LOI_F25+LOI_F3+LOI_F27+LOI_F39+LOI_F34+LOI_F4+
   LOI_F5+LOI_F1+LOI_F15+LOI_F8+LOI_F24+LOI_F21+LOI_F28
 factor3=~ LOI_F13+LOI_F7+LOI_F20+LOI_F10+LOI_F2+LOI_F18+LOI_F32+LOI_F42+LOI_F45+LOI_F16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.972       0.910
#Tucker-Lewis Index (TLI)                       0.970       0.906
#Robust Comparative Fit Index (CFI)                         0.700
#Robust Tucker-Lewis Index (TLI)                            0.686
#RMSEA                                          0.069       0.057
#Robust RMSEA                                               0.104
#SRMR                                           0.080       0.080

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.746    0.022   34.041    0.000    0.746    0.746
#  factor3           0.581    0.031   18.827    0.000    0.581    0.581
#factor2 ~~                                                            
#  factor3           0.729    0.023   32.047    0.000    0.729    0.729


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .51


CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.597       0.458
#Tucker-Lewis Index (TLI)                       0.580       0.435
#Robust Comparative Fit Index (CFI)                         0.675
#Robust Tucker-Lewis Index (TLI)                            0.662
#RMSEA                                          0.260       0.140
#Robust RMSEA                                               0.108
#SRMR                                           0.270       0.270

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .49

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1=~ LOI_F36+LOI_F26+LOI_F12+LOI_F44+LOI_F9+LOI_F43+LOI_F22+LOI_F47+LOI_F31+LOI_F48+LOI_F14+
   LOI_F33+LOI_F19+LOI_F29+LOI_F38+LOI_F6+LOI_F17+LOI_F11+LOI_F23+LOI_F46+LOI_F35+LOI_F41
 factor2=~ LOI_F30+LOI_F40+LOI_F49+LOI_F50+LOI_F37+LOI_F25+LOI_F3+LOI_F27+LOI_F39+LOI_F34+LOI_F4+
   LOI_F5+LOI_F1+LOI_F15+LOI_F8+LOI_F24+LOI_F21+LOI_F28
 factor3=~ LOI_F13+LOI_F7+LOI_F20+LOI_F10+LOI_F2+LOI_F18+LOI_F32+LOI_F42+LOI_F45+LOI_F16
 general=~ LOI_F36+LOI_F26+LOI_F12+LOI_F44+LOI_F9+LOI_F43+LOI_F22+LOI_F47+LOI_F31+LOI_F48+LOI_F14+
           LOI_F33+LOI_F19+LOI_F29+LOI_F38+LOI_F6+LOI_F17+LOI_F11+LOI_F23+LOI_F46+LOI_F35+LOI_F41+
           LOI_F30+LOI_F40+LOI_F49+LOI_F50+LOI_F37+LOI_F25+LOI_F3+LOI_F27+LOI_F39+LOI_F34+LOI_F4+
           LOI_F5+LOI_F1+LOI_F15+LOI_F8+LOI_F24+LOI_F21+LOI_F28+
           LOI_F13+LOI_F7+LOI_F20+LOI_F10+LOI_F2+LOI_F18+LOI_F32+LOI_F42+LOI_F45+LOI_F16
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.987       0.943
#Tucker-Lewis Index (TLI)                       0.986       0.938
#Robust Comparative Fit Index (CFI)                         0.751
#Robust Tucker-Lewis Index (TLI)                            0.729
#RMSEA                                          0.047       0.046
#Robust RMSEA                                               0.097
#SRMR                                           0.065       0.065

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .55

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    LOI_F36           0.604    0.039   15.546    0.000    0.604    0.604
#    LOI_F26           0.632    0.037   17.224    0.000    0.632    0.632
#    LOI_F12           0.752    0.044   16.966    0.000    0.752    0.752
#    LOI_F44           0.525    0.046   11.439    0.000    0.525    0.525
#    LOI_F9            0.620    0.034   18.168    0.000    0.620    0.620
#    LOI_F43           0.532    0.038   14.061    0.000    0.532    0.532
#    LOI_F22           0.557    0.039   14.270    0.000    0.557    0.557
#    LOI_F47           0.491    0.039   12.448    0.000    0.491    0.491
#    LOI_F31           0.502    0.041   12.268    0.000    0.502    0.502
#    LOI_F48           0.637    0.062   10.313    0.000    0.637    0.637
#    LOI_F14           0.449    0.045    9.979    0.000    0.449    0.449
#    LOI_F33           0.363    0.050    7.261    0.000    0.363    0.363
#    LOI_F19           0.470    0.043   10.959    0.000    0.470    0.470
#    LOI_F29           0.350    0.067    5.195    0.000    0.350    0.350
#    LOI_F38           0.390    0.061    6.369    0.000    0.390    0.390
#    LOI_F6            0.272    0.043    6.289    0.000    0.272    0.272
#    LOI_F17           0.389    0.044    8.751    0.000    0.389    0.389
#    LOI_F11           0.408    0.050    8.096    0.000    0.408    0.408
#    LOI_F23           0.369    0.052    7.048    0.000    0.369    0.369
#    LOI_F46           0.151    0.058    2.604    0.009    0.151    0.151
#    LOI_F35           0.282    0.049    5.707    0.000    0.282    0.282
#    LOI_F41           0.202    0.052    3.876    0.000    0.202    0.202
#  factor2 =~                                                            
#    LOI_F30           0.398    0.043    9.266    0.000    0.398    0.398
#    LOI_F40           0.437    0.044    9.831    0.000    0.437    0.437
#    LOI_F49           0.188    0.045    4.165    0.000    0.188    0.188
#    LOI_F50           0.710    0.036   19.646    0.000    0.710    0.710
#    LOI_F37           0.181    0.046    3.932    0.000    0.181    0.181
#    LOI_F25           0.198    0.048    4.104    0.000    0.198    0.198
#    LOI_F3            0.417    0.046    9.069    0.000    0.417    0.417
#    LOI_F27           0.097    0.044    2.193    0.028    0.097    0.097
#    LOI_F39           0.170    0.048    3.514    0.000    0.170    0.170
#    LOI_F34           0.137    0.054    2.548    0.011    0.137    0.137
#    LOI_F4            0.232    0.047    4.914    0.000    0.232    0.232
#    LOI_F5            0.169    0.048    3.528    0.000    0.169    0.169
#    LOI_F1            0.675    0.037   18.480    0.000    0.675    0.675
#    LOI_F15           0.047    0.052    0.903    0.367    0.047    0.047
#    LOI_F8           -0.076    0.058   -1.295    0.195   -0.076   -0.076
#    LOI_F24           0.029    0.053    0.540    0.589    0.029    0.029
#    LOI_F21          -0.001    0.053   -0.010    0.992   -0.001   -0.001
#    LOI_F28          -0.020    0.050   -0.408    0.683   -0.020   -0.020
#  factor3 =~                                                            
#    LOI_F13           0.691    0.028   24.346    0.000    0.691    0.691
#    LOI_F7            0.613    0.034   18.270    0.000    0.613    0.613
#    LOI_F20           0.556    0.035   16.115    0.000    0.556    0.556
#    LOI_F10           0.647    0.032   20.269    0.000    0.647    0.647
#    LOI_F2            0.660    0.034   19.468    0.000    0.660    0.660
#    LOI_F18           0.568    0.036   15.956    0.000    0.568    0.568
#    LOI_F32           0.444    0.035   12.737    0.000    0.444    0.444
#    LOI_F42           0.273    0.039    7.043    0.000    0.273    0.273
#    LOI_F45           0.337    0.046    7.314    0.000    0.337    0.337
#    LOI_F16           0.199    0.044    4.528    0.000    0.199    0.199
#  general =~                                                            
#    LOI_F36           0.585    0.041   14.194    0.000    0.585    0.585
#    LOI_F26           0.545    0.042   12.998    0.000    0.545    0.545
#    LOI_F12           0.381    0.057    6.667    0.000    0.381    0.381
#    LOI_F44           0.626    0.042   14.860    0.000    0.626    0.626
#    LOI_F9            0.443    0.042   10.659    0.000    0.443    0.443
#    LOI_F43           0.617    0.035   17.378    0.000    0.617    0.617
#    LOI_F22           0.615    0.036   17.249    0.000    0.615    0.615
#    LOI_F47           0.580    0.037   15.652    0.000    0.580    0.580
#    LOI_F31           0.558    0.040   13.973    0.000    0.558    0.558
#    LOI_F48           0.483    0.061    7.908    0.000    0.483    0.483
#    LOI_F14           0.658    0.038   17.460    0.000    0.658    0.658
#    LOI_F33           0.661    0.037   17.730    0.000    0.661    0.661
#    LOI_F19           0.595    0.038   15.569    0.000    0.595    0.595
#    LOI_F29           0.605    0.045   13.510    0.000    0.605    0.605
#    LOI_F38           0.541    0.050   10.921    0.000    0.541    0.541
#    LOI_F6            0.689    0.032   21.509    0.000    0.689    0.689
#    LOI_F17           0.537    0.041   12.990    0.000    0.537    0.537
#    LOI_F11           0.515    0.047   11.003    0.000    0.515    0.515
#    LOI_F23           0.429    0.048    8.921    0.000    0.429    0.429
#    LOI_F46           0.656    0.044   15.028    0.000    0.656    0.656
#    LOI_F35           0.484    0.046   10.490    0.000    0.484    0.484
#    LOI_F41           0.527    0.043   12.314    0.000    0.527    0.527
#    LOI_F30           0.566    0.036   15.799    0.000    0.566    0.566
#    LOI_F40           0.522    0.039   13.309    0.000    0.522    0.522
#    LOI_F49           0.774    0.027   28.854    0.000    0.774    0.774
#    LOI_F50           0.416    0.042    9.982    0.000    0.416    0.416
#    LOI_F37           0.685    0.030   22.827    0.000    0.685    0.685
#    LOI_F25           0.602    0.034   17.654    0.000    0.602    0.602
#    LOI_F3            0.352    0.045    7.849    0.000    0.352    0.352
#    LOI_F27           0.768    0.027   28.101    0.000    0.768    0.768
#    LOI_F39           0.595    0.035   17.119    0.000    0.595    0.595
#    LOI_F34           0.638    0.037   17.267    0.000    0.638    0.638
#    LOI_F4            0.573    0.038   15.005    0.000    0.573    0.573
#    LOI_F5            0.552    0.037   15.071    0.000    0.552    0.552
#    LOI_F1            0.118    0.049    2.415    0.016    0.118    0.118
#    LOI_F15           0.657    0.034   19.611    0.000    0.657    0.657
#    LOI_F8            0.550    0.040   13.694    0.000    0.550    0.550
#    LOI_F24           0.739    0.029   25.467    0.000    0.739    0.739
#    LOI_F21           0.676    0.033   20.705    0.000    0.676    0.676
#    LOI_F28           0.592    0.038   15.640    0.000    0.592    0.592
#    LOI_F13           0.542    0.036   14.868    0.000    0.542    0.542
#    LOI_F7            0.531    0.039   13.686    0.000    0.531    0.531
#    LOI_F20           0.631    0.034   18.722    0.000    0.631    0.631
#    LOI_F10           0.488    0.040   12.347    0.000    0.488    0.488
#    LOI_F2            0.346    0.045    7.667    0.000    0.346    0.346
#    LOI_F18           0.541    0.038   14.168    0.000    0.541    0.541
#    LOI_F32           0.661    0.032   20.509    0.000    0.661    0.661
#    LOI_F42           0.718    0.030   23.785    0.000    0.718    0.718
#    LOI_F45           0.616    0.037   16.702    0.000    0.616    0.616
#    LOI_F16           0.682    0.032   21.190    0.000    0.682    0.682


library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.6402095      0.6497959      0.9636754      0.8242545 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.4140076 0.18892796 0.5859924 0.9360846 0.3222184 0.8823608 0.9343762
#factor2 0.2135077 0.06595158 0.7864923 0.8788701 0.1089793 0.7335202 0.8673315
#factor3 0.4468764 0.10491098 0.5531236 0.8974977 0.3619923 0.8158277 0.9206038
#general 0.6402095 0.64020949 0.6402095 0.9636754 0.8242545 0.9656269 0.9689135





################################################################
### Analysis Questionnaires for Obsessive-compulsive disorders
### OCI-R González https://www.researchgate.net/publication/277266486_Validacion_del_inventario_de_obsesiones_y_compulsiones_revisado_OCI-R_para_su_uso_en_poblacion_adolescente_espanola

OCD_Piqueras_Gonzalez <- read.csv("~/Self regulation/Second round of analyses/OCD/OCD_Gonzalez.csv")
colnames(OCD_Gonzalez)
mydata  <- as.data.frame(OCD_Gonzalez[,6:23])
mydata[mydata == 987] <- NA
mydata[mydata == 33] <- 3
mydata <- na.omit(mydata)
min(mydata)
max(mydata) #response alternatives = 5

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 2 components
# Eigenvalue 1 = 3.89; eigenvalue 2 = 0.78

rho <- polychoric(mydata)$rho
min(mydata)
max(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 3 components
# Eigenvalue 1 = 4.73; eigenvalue 2 = .85

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.083, RMSR=.07, TLI=.694

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.108, RMSR=.08, TLI=.644

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities

# Give solution with 5 factors 
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.044, RMSR=.03, TLI=.915
#     MR1  MR5  MR4  MR2  MR3
#MR1 1.00 0.42 0.43 0.13 0.30
#MR5 0.42 1.00 0.39 0.44 0.41
#MR4 0.43 0.39 1.00 0.27 0.32
#MR2 0.13 0.44 0.27 1.00 0.20
#MR3 0.30 0.41 0.32 0.20 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.071, RMSR=.03, TLI=.847
#     MR5  MR1  MR2  MR4  MR3
#MR5 1.00 0.45 0.41 0.38 0.48
#MR1 0.45 1.00 0.16 0.42 0.37
#MR2 0.41 0.16 1.00 0.30 0.26
#MR4 0.38 0.42 0.30 1.00 0.35
#MR3 0.48 0.37 0.26 0.35 1.00


# Single factor model lavaan
UNImodel= '
 general=~ OCIR1+OCIR2+OCIR3+OCIR4+OCIR5+OCIR6+OCIR7+OCIR8+OCIR9+
            OCIR10+OCIR11+OCIR12+OCIR13+OCIR14+OCIR15+OCIR16+OCIR17+OCIR18
'
library(lavaan)
# Not all ordered models can be run without warning, hence MLR
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.730       0.744
#Tucker-Lewis Index (TLI)                       0.694       0.710
#Robust Comparative Fit Index (CFI)                         0.743
#Robust Tucker-Lewis Index (TLI)                            0.709
#RMSEA                                          0.084       0.073
#Robust RMSEA                                               0.081
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .21

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 6 factors based on theoretical analysis
EGAmodel= '
 hoarding =~ OCIR1+OCIR7+OCIR13
 washing  =~ OCIR5+OCIR11+OCIR17
 obsessing=~ OCIR6+OCIR12+OCIR18
 ordering =~ OCIR3+OCIR9+OCIR15
 checking =~ OCIR2+OCIR8+OCIR14
 neutralising=~ OCIR4+OCIR10+OCIR16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.895
#Tucker-Lewis Index (TLI)                       0.852       0.866
#Robust Comparative Fit Index (CFI)                         0.896
#Robust Tucker-Lewis Index (TLI)                            0.868
#RMSEA                                          0.059       0.050
#Robust RMSEA                                               0.055
#SRMR                                           0.051       0.051

#Covariances:
#   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#hoarding ~~                                                           
#   washing           0.496    0.063    7.893    0.000    0.496    0.496
#   obsessing         0.749    0.051   14.617    0.000    0.749    0.749
#   ordering          0.715    0.061   11.768    0.000    0.715    0.715
#   checking          0.647    0.054   11.915    0.000    0.647    0.647
#   neutralising      0.572    0.069    8.314    0.000    0.572    0.572
#washing ~~                                                            
#   obsessing         0.475    0.053    8.899    0.000    0.475    0.475
#   ordering          0.749    0.057   13.071    0.000    0.749    0.749
#   checking          0.581    0.055   10.552    0.000    0.581    0.581
#   neutralising      0.723    0.056   12.883    0.000    0.723    0.723
#obsessing ~~                                                          
#   ordering          0.510    0.057    8.978    0.000    0.510    0.510
#   checking          0.528    0.051   10.358    0.000    0.528    0.528
#   neutralising      0.556    0.054   10.290    0.000    0.556    0.556
#ordering ~~                                                           
#   checking          0.692    0.055   12.479    0.000    0.692    0.692
#   neutralising      0.673    0.060   11.202    0.000    0.673    0.673
#checking ~~                                                           
#   neutralising      0.715    0.052   13.878    0.000    0.715    0.715

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .295


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.539       0.538
#Tucker-Lewis Index (TLI)                       0.478       0.477
#Robust Comparative Fit Index (CFI)                         0.546
#Robust Tucker-Lewis Index (TLI)                            0.486
#RMSEA                                          0.110       0.098
#Robust RMSEA                                               0.108
#SRMR                                           0.186       0.186

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .321

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 hoarding =~ OCIR1+OCIR7+OCIR13
 washing  =~ OCIR5+OCIR11+OCIR17
 obsessing=~ OCIR6+OCIR12+OCIR18
 ordering =~ OCIR3+OCIR9+OCIR15
 checking =~ OCIR2+OCIR8+OCIR14
 neutralising=~ OCIR4+OCIR10+OCIR16
 general=~ OCIR1+OCIR2+OCIR3+OCIR4+OCIR5+OCIR6+OCIR7+OCIR8+OCIR9+
            OCIR10+OCIR11+OCIR12+OCIR13+OCIR14+OCIR15+OCIR16+OCIR17+OCIR18
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.889       0.900
#Tucker-Lewis Index (TLI)                       0.855       0.870
#Robust Comparative Fit Index (CFI)                         0.902
#Robust Tucker-Lewis Index (TLI)                            0.871
#RMSEA                                          0.058       0.049
#Robust RMSEA                                               0.054
#SRMR                                           0.050       0.050

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .316

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  hoarding =~                                                           
#    OCIR1             0.134    0.365    0.367    0.713    0.134    0.134
#    OCIR7             1.383    3.725    0.371    0.710    1.383    1.383
#    OCIR13            0.053    0.151    0.351    0.726    0.053    0.053
#  washing =~                                                            
#    OCIR5             0.158    0.072    2.204    0.028    0.158    0.158
#    OCIR11            0.473    0.164    2.883    0.004    0.473    0.473
#    OCIR17            0.574    0.196    2.923    0.003    0.574    0.574
#  obsessing =~                                                          
#    OCIR6             0.241    0.068    3.530    0.000    0.241    0.241
#    OCIR12            1.025    0.217    4.734    0.000    1.025    1.025
#    OCIR18            0.352    0.088    4.002    0.000    0.352    0.352
#  ordering =~                                                           
#    OCIR3             0.448    0.162    2.761    0.006    0.448    0.448
#    OCIR9             0.270    0.104    2.606    0.009    0.270    0.270
#    OCIR15            0.231    0.091    2.548    0.011    0.231    0.231
#  checking =~                                                           
#    OCIR2             0.148    0.090    1.642    0.101    0.148    0.148
#    OCIR8             0.885    0.416    2.127    0.033    0.885    0.885
#    OCIR14            0.285    0.147    1.945    0.052    0.285    0.285
#  neutralising =~                                                       
#    OCIR4             0.217    0.089    2.448    0.014    0.217    0.217
#    OCIR10            0.519    0.170    3.055    0.002    0.519    0.519
#    OCIR16            0.335    0.119    2.810    0.005    0.335    0.335
#  general =~                                                            
#    OCIR1             0.412    0.041    9.949    0.000    0.412    0.412
#    OCIR2             0.477    0.043   11.211    0.000    0.477    0.477
#    OCIR3             0.410    0.043    9.426    0.000    0.410    0.410
#    OCIR4             0.532    0.045   11.938    0.000    0.532    0.532
#    OCIR5             0.477    0.044   10.736    0.000    0.477    0.477
#    OCIR6             0.435    0.042   10.324    0.000    0.435    0.435
#    OCIR7             0.404    0.045    8.967    0.000    0.404    0.404
#    OCIR8             0.555    0.038   14.467    0.000    0.555    0.555
#    OCIR9             0.487    0.042   11.693    0.000    0.487    0.487
#    OCIR10            0.564    0.043   13.156    0.000    0.564    0.564
#    OCIR11            0.501    0.040   12.429    0.000    0.501    0.501
#    OCIR12            0.562    0.036   15.498    0.000    0.562    0.562
#    OCIR13            0.498    0.040   12.456    0.000    0.498    0.498
#    OCIR14            0.579    0.039   14.872    0.000    0.579    0.579
#    OCIR15            0.572    0.035   16.537    0.000    0.572    0.572
#    OCIR16            0.504    0.045   11.114    0.000    0.504    0.504
#    OCIR17            0.494    0.042   11.687    0.000    0.494    0.494
#    OCIR18            0.484    0.041   11.740    0.000    0.484    0.484


library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.4550544      0.8823529      0.8767101      0.7772854 

#$FactorLevelIndices
#              ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#hoarding     0.7692191 0.19561532 0.2307809 0.6460590 0.3311660 1.9319344 1.5084552
#washing      0.4446147 0.05847310 0.5553853 0.6526007 0.2661961 0.4460719 0.7057511
#obsessing    0.6250844 0.12468395 0.3749156 0.7440431 0.4355008 1.0517807 1.2375386
#ordering     0.3089925 0.03310948 0.6910075 0.5738980 0.1623186 0.2789291 0.5473813
#checking     0.5046273 0.08969719 0.4953727 0.7441196 0.3008700 0.7890620 1.0261371
#neutralising 0.3340710 0.04336654 0.6659290 0.6477123 0.2006107 0.3526760 0.6384565
#general      0.4550544 0.45505443 0.4550544 0.8767101 0.7772854 0.8596869 0.9048384


# model gives warning for negative variances; loadings not well estimated;
# model with 5 factors on the basis of EGA was used instead but same warning

BIFmodel= '
 factor1 =~ OCIR1+OCIR7+OCIR3+OCIR9
 factor2 =~ OCIR10+OCIR16+OCIR4+OCIR5
 factor3 =~ OCIR6+OCIR12+OCIR13+OCIR18
 factor4 =~ OCIR2+OCIR8+OCIR14
 factor5 =~ OCIR11+OCIR15+OCIR17
 general=~ OCIR1+OCIR2+OCIR3+OCIR4+OCIR5+OCIR6+OCIR7+OCIR8+OCIR9+
            OCIR10+OCIR11+OCIR12+OCIR13+OCIR14+OCIR15+OCIR16+OCIR17+OCIR18
'
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)








################################################################
### Analysis Questionnaires for Obsessive-compulsive disorders
### OCI-R Ignatova et al https://bulgarian-journal-of-psychiatry.bg/wp-content/uploads/2023/02/BSP-1-2023.pdf#page=12

library(readxl)
OCI_R_dataset_Ignatova <- read_excel("OCI-R dataset Ignatova.xlsx")
colnames(OCI_R_dataset_Ignatova)
mydata  <- as.data.frame(OCI_R_dataset_Ignatova[,2:19])
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata) #response alternatives = 5

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 2 components
# Eigenvalue 1 = 5.46; eigenvalue 2 = 1.27

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 2 components
# Eigenvalue 1 = 6.91; eigenvalue 2 = 1.42

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.30, RMSEA=.144, RMSR=.11, TLI=.55

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.212, RMSR=.12, TLI=.447

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities

# Give solution with 6 factors 
fit4 <- fa(mydata,6)
fit4
diagram(fit4)
# %variance explained=.58, RMSEA=.054, RMSR=.03, TLI=.936
#     MR1  MR2  MR5  MR4  MR3  MR6
#MR1 1.00 0.30 0.51 0.45 0.34 0.35
#MR2 0.30 1.00 0.29 0.15 0.20 0.45
#MR5 0.51 0.29 1.00 0.36 0.32 0.34
#MR4 0.45 0.15 0.36 1.00 0.30 0.25
#MR3 0.34 0.20 0.32 0.30 1.00 0.26
#MR6 0.35 0.45 0.34 0.25 0.26 1.00

fit4 <- fa(rho,6,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.163, RMSR=.03, TLI=.666
#     MR1  MR2  MR4  MR5  MR6  MR3
#MR1 1.00 0.29 0.43 0.35 0.46 0.31
#MR2 0.29 1.00 0.19 0.48 0.27 0.19
#MR4 0.43 0.19 1.00 0.33 0.37 0.34
#MR5 0.35 0.48 0.33 1.00 0.38 0.30
#MR6 0.46 0.27 0.37 0.38 1.00 0.29
#MR3 0.31 0.19 0.34 0.30 0.29 1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~ OC1+OC2+OC3+OC4+OC5+OC6+OC7+OC8+OC9+
            OC10+OC11+OC12+OC13+OC14+OC15+OC16+OC17+OC18
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# did not converge
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.606       0.607
#Tucker-Lewis Index (TLI)                       0.554       0.555
#Robust Comparative Fit Index (CFI)                         0.619
#Robust Tucker-Lewis Index (TLI)                            0.568
#RMSEA                                          0.147       0.125
#Robust RMSEA                                               0.142
#SRMR                                           0.107       0.107

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .226

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 6 factors based on theoretical analysis
EGAmodel= '
 hoarding =~ OC1+OC7+OC13
 washing  =~ OC5+OC11+OC17
 obsessing=~ OC6+OC12+OC18
 ordering =~ OC3+OC9+OC15
 checking =~ OC2+OC8+OC14
 neutralising=~ OC4+OC10+OC16
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.938       0.955
#Tucker-Lewis Index (TLI)                       0.921       0.943
#Robust Comparative Fit Index (CFI)                         0.957
#Robust Tucker-Lewis Index (TLI)                            0.946
#RMSEA                                          0.062       0.045
#Robust RMSEA                                               0.051
#SRMR                                           0.060       0.060

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  hoarding ~~                                                           
#    washing           0.444    0.100    4.460    0.000    0.444    0.444
#    obsessing         0.542    0.088    6.129    0.000    0.542    0.542
#    ordering          0.425    0.087    4.908    0.000    0.425    0.425
#    checking          0.401    0.096    4.159    0.000    0.401    0.401
#    neutralising      0.488    0.095    5.128    0.000    0.488    0.488
#  washing ~~                                                            
#    obsessing         0.314    0.087    3.619    0.000    0.314    0.314
#    ordering          0.516    0.075    6.833    0.000    0.516    0.516
#    checking          0.495    0.089    5.547    0.000    0.495    0.495
#    neutralising      0.579    0.112    5.160    0.000    0.579    0.579
#  obsessing ~~                                                          
#    ordering          0.339    0.086    3.934    0.000    0.339    0.339
#    checking          0.269    0.099    2.707    0.007    0.269    0.269
#    neutralising      0.451    0.094    4.811    0.000    0.451    0.451
#  ordering ~~                                                           
#    checking          0.570    0.074    7.681    0.000    0.570    0.570
#    neutralising      0.681    0.061   11.243    0.000    0.681    0.681
#  checking ~~                                                           
#    neutralising      0.555    0.102    5.434    0.000    0.555    0.555

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .551

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.754       0.761
#Tucker-Lewis Index (TLI)                       0.721       0.729
#Robust Comparative Fit Index (CFI)                         0.771
#Robust Tucker-Lewis Index (TLI)                            0.740
#RMSEA                                          0.116       0.097
#Robust RMSEA                                               0.110
#SRMR                                           0.246       0.246

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .541


# Bifactor model
BIFmodel= '
 hoarding =~ OC1+OC7+OC13
 washing  =~ OC5+OC11+OC17
 obsessing=~ OC6+OC12+OC18
 ordering =~ OC3+OC9+OC15
 checking =~ OC2+OC8+OC14
 neutralising=~ OC4+OC10+OC16
 general=~ OC1+OC2+OC3+OC4+OC5+OC6+OC7+OC8+OC9+
            OC10+OC11+OC12+OC13+OC14+OC15+OC16+OC17+OC18
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.948       0.963
#Tucker-Lewis Index (TLI)                       0.932       0.952
#Robust Comparative Fit Index (CFI)                         0.966
#Robust Tucker-Lewis Index (TLI)                            0.955
#RMSEA                                          0.057       0.041
#Robust RMSEA                                               0.046
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .555

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  hoarding =~                                                           
#    OC1               0.298    0.094    3.163    0.002    0.298    0.285
#    OC7               0.525    0.151    3.471    0.001    0.525    0.560
#    OC13              0.707    0.178    3.980    0.000    0.707    0.586
#  washing =~                                                            
#    OC5               0.256    0.106    2.412    0.016    0.256    0.298
#    OC11              0.608    0.130    4.684    0.000    0.608    0.586
#    OC17              0.500    0.112    4.469    0.000    0.500    0.569
#  obsessing =~                                                          
#    OC6               0.663    0.090    7.336    0.000    0.663    0.577
#    OC12              0.938    0.091   10.358    0.000    0.938    0.829
#    OC18              0.772    0.105    7.319    0.000    0.772    0.664
#  ordering =~                                                           
#    OC3               0.663    0.095    6.993    0.000    0.663    0.564
#    OC9               0.542    0.084    6.435    0.000    0.542    0.450
#    OC15              0.700    0.098    7.146    0.000    0.700    0.583
#  checking =~                                                           
#    OC2               0.276    0.091    3.033    0.002    0.276    0.238
#    OC8               0.693    0.154    4.484    0.000    0.693    0.682
#    OC14              0.570    0.155    3.675    0.000    0.570    0.586
#  neutralising =~                                                       
#    OC4               0.235    0.087    2.704    0.007    0.235    0.271
#    OC10              0.301    0.093    3.249    0.001    0.301    0.473
#    OC16              0.484    0.136    3.566    0.000    0.484    0.493
#  general =~                                                            
#    OC1               0.430    0.091    4.743    0.000    0.430    0.411
#    OC2               0.616    0.083    7.400    0.000    0.616    0.532
#    OC3               0.756    0.083    9.073    0.000    0.756    0.642
#    OC4               0.622    0.083    7.446    0.000    0.622    0.716
#    OC5               0.376    0.087    4.294    0.000    0.376    0.438
#    OC6               0.511    0.100    5.100    0.000    0.511    0.445
#    OC7               0.378    0.084    4.514    0.000    0.378    0.403
#    OC8               0.537    0.095    5.642    0.000    0.537    0.529
#    OC9               0.809    0.077   10.549    0.000    0.809    0.672
#    OC10              0.391    0.082    4.760    0.000    0.391    0.615
#    OC11              0.482    0.099    4.884    0.000    0.482    0.464
#    OC12              0.463    0.102    4.565    0.000    0.463    0.409
#    OC13              0.506    0.104    4.889    0.000    0.506    0.419
#    OC14              0.402    0.096    4.187    0.000    0.402    0.414
#    OC15              0.835    0.087    9.613    0.000    0.835    0.695
#    OC16              0.442    0.102    4.341    0.000    0.442    0.450
#    OC17              0.421    0.092    4.570    0.000    0.421    0.478
#    OC18              0.567    0.106    5.367    0.000    0.567    0.488

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.4853863      0.8823529      0.9270283      0.7897171 

#$FactorLevelIndices
#              ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#hoarding     0.5926635 0.07275144 0.4073365 0.6702585 0.3845153 0.5162934 0.7414267
#washing      0.5426478 0.07445822 0.4573522 0.7140387 0.3750299 0.5231664 0.7529717
#obsessing    0.7074914 0.14403311 0.2925086 0.8667842 0.6100864 0.7769574 0.9095831
#ordering     0.3897476 0.08480566 0.6102524 0.8924866 0.3455188 0.5524051 0.7871874
#checking     0.5408697 0.08531174 0.4591303 0.7603984 0.3880983 0.5923029 0.8139961
#neutralising 0.3305022 0.05325353 0.6694978 0.7749388 0.2520836 0.4076753 0.6660788
#general      0.4853863 0.48538630 0.4853863 0.9270283 0.7897171 0.8829943 0.9055806





################################################################
### Muela et al 2023 https://bmcpsychology.biomedcentral.com/articles/10.1186/s40359-023-01439-1
### Gambling and video gaming
### data https://osf.io/ka2jp?view_only=9831ce7702c34347ac67b45719ddf643

### not included because not clear whether unifactorial or hierarchical model

