################################################################
### Beliefs scale Barnby et al. (2019) https://peerj.com/articles/6819/
### Confirmatory study, agreement https://osf.io/hdkw7

Beliefs_Barnby <- read.csv("Beliefs_Barnby.csv")
colnames(Beliefs_Barnby)
mydata = Beliefs_Barnby[,seq(10, ncol(Beliefs_Barnby), 3)]
mydata <- mydata[,1:50]
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 11 response alternatives


library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 9.56; eigenvalue 2 = 3.08

rho <- polychoric(mydata)$rho
# no results because more than 8 categories
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Eigenvalue 1 = 11.6; eigenvalue 2 = 4.21

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.086, RMSR=.1, TLI=.79

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.121, RMSR=.13, TLI=.378

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities and response bias

# Give solution with 5 factors
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.044, RMSR=.04, TLI=.86
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00  0.35 -0.32  0.01  0.22
#MR2  0.35  1.00 -0.15 -0.03  0.20
#MR4 -0.32 -0.15  1.00  0.04 -0.17
#MR3  0.01 -0.03  0.04  1.00 -0.05
#MR5  0.22  0.20 -0.17 -0.05  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.078, RMSR=.04, TLI=.736
#      MR3   MR1   MR2   MR4   MR5
#MR3  1.00  0.40  0.10 -0.18 -0.12
#MR1  0.40  1.00  0.04 -0.31 -0.20
#MR2  0.10  0.04  1.00 -0.08 -0.06
#MR4 -0.18 -0.31 -0.08  1.00  0.18
#MR5 -0.12 -0.20 -0.06  0.18  1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  S_G1_A_1+S_G2_A_1+S_G3_A_1+S_G4_A_1+S_G5_A_1+S_S1_A_1+S_S2_A_1+S_S3_A_1+
            S_S4_A_1+S_S5_A_1+R_G1_A_1+R_G2_A_1+R_G3_A_1+R_G4_A_1+R_G5_A_1+R_S1_A_1+  
            R_S2_A_1+R_S3_A_1+R_S4_A_1+R_S5_A_1+P_G1_A_1+P_G2_A_1+P_G3_A_1+P_G4_A_1+  
            P_G5_A_1+P_S1_A_1+P_S2_A_1+P_S3_A_1+P_S4_A_1+P_S5_A_1+Pol_S1_A_1+Pol_S2_A_1+
            Pol_S3_A_1+Pol_S4_A_1+Pol_S5_A_1+Pol_G1_A_1+Pol_G2_A_1+Pol_G3_A_1+Pol_G4_A_1+Pol_G5_A_1+
            M_G1_A_1+M_G2_A_1+M_G3_A_1+M_G4_A_1+M_G5_r_A_1+M_S1_A_1+M_S2_r_A_1+M_S3_A_1+  
            M_S4_r_A_1+M_S5_A_1 
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.872       0.764
#Tucker-Lewis Index (TLI)                       0.867       0.754
#Robust Comparative Fit Index (CFI)                         0.428
#Robust Tucker-Lewis Index (TLI)                            0.404
#RMSEA                                          0.132       0.098
#Robust RMSEA                                               0.120
#SRMR                                           0.126       0.126

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .200

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.502       0.506
#Tucker-Lewis Index (TLI)                       0.481       0.485
#Robust Comparative Fit Index (CFI)                         0.512
#Robust Tucker-Lewis Index (TLI)                            0.491
#RMSEA                                          0.088       0.080
#Robust RMSEA                                               0.086
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .167

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 5 factors put into the inventory
EGAmodel= '
 science=~  S_G1_A_1+S_G2_A_1+S_G3_A_1+S_G4_A_1+S_G5_A_1+S_S1_A_1+S_S2_A_1+S_S3_A_1+
            S_S4_A_1+S_S5_A_1
 religio=~  R_G1_A_1+R_G2_A_1+R_G3_A_1+R_G4_A_1+R_G5_A_1+R_S1_A_1+  
            R_S2_A_1+R_S3_A_1+R_S4_A_1+R_S5_A_1
 paranor=~  P_G1_A_1+P_G2_A_1+P_G3_A_1+P_G4_A_1+  
            P_G5_A_1+P_S1_A_1+P_S2_A_1+P_S3_A_1+P_S4_A_1+P_S5_A_1
 politic=~  Pol_S1_A_1+Pol_S2_A_1+
            Pol_S3_A_1+Pol_S4_A_1+Pol_S5_A_1+Pol_G1_A_1+Pol_G2_A_1+Pol_G3_A_1+Pol_G4_A_1+Pol_G5_A_1
 moralit=~  M_G1_A_1+M_G2_A_1+M_G3_A_1+M_G4_A_1+M_G5_r_A_1+M_S1_A_1+M_S2_r_A_1+M_S3_A_1+  
            M_S4_r_A_1+M_S5_A_1 
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.938       0.886
#Tucker-Lewis Index (TLI)                       0.935       0.880
#Robust Comparative Fit Index (CFI)                         0.687
#Robust Tucker-Lewis Index (TLI)                            0.671
#RMSEA                                          0.092       0.069
#Robust RMSEA                                               0.089
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .448

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  science ~~                                                            
#    religio           0.657    0.029   22.303    0.000    0.657    0.657
#    paranor           0.609    0.030   20.264    0.000    0.609    0.609
#    politic          -0.663    0.032  -20.674    0.000   -0.663   -0.663
#    moralit           0.065    0.051    1.275    0.202    0.065    0.065
#  religio ~~                                                            
#    paranor           0.601    0.027   22.043    0.000    0.601    0.601
#    politic          -0.554    0.031  -17.887    0.000   -0.554   -0.554
#    moralit          -0.027    0.048   -0.558    0.577   -0.027   -0.027
#  paranor ~~                                                            
#    politic          -0.406    0.041   -9.815    0.000   -0.406   -0.406
#    moralit           0.195    0.049    3.942    0.000    0.195    0.195
#  politic ~~                                                            
#    moralit          -0.281    0.047   -6.030    0.000   -0.281   -0.281

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.735       0.740
#Tucker-Lewis Index (TLI)                       0.722       0.727
#Robust Comparative Fit Index (CFI)                         0.748
#Robust Tucker-Lewis Index (TLI)                            0.735
#RMSEA                                          0.064       0.058
#Robust RMSEA                                               0.062
#SRMR                                           0.082       0.082

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  science ~~                                                            
#    religio           0.645    0.044   14.777    0.000    0.645    0.645
#    paranor           0.463    0.066    7.067    0.000    0.463    0.463
#    politic          -0.742    0.055  -13.461    0.000   -0.742   -0.742
#    moralit          -0.022    0.084   -0.263    0.793   -0.022   -0.022
#  religio ~~                                                            
#    paranor           0.483    0.050    9.664    0.000    0.483    0.483
#    politic          -0.597    0.065   -9.196    0.000   -0.597   -0.597
#    moralit          -0.072    0.070   -1.033    0.302   -0.072   -0.072
#  paranor ~~                                                            
#    politic          -0.335    0.078   -4.317    0.000   -0.335   -0.335
#    moralit           0.143    0.071    2.016    0.044    0.143    0.143
#  politic ~~                                                            
#    moralit          -0.151    0.139   -1.085    0.278   -0.151   -0.151

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .287

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.681       0.677
#Tucker-Lewis Index (TLI)                       0.667       0.663
#Robust Comparative Fit Index (CFI)                         0.690
#Robust Tucker-Lewis Index (TLI)                            0.677
#RMSEA                                          0.070       0.064
#Robust RMSEA                                               0.068
#SRMR                                           0.150       0.150

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .263

# Bifactor model
BIFmodel= '
 science=~  S_G1_A_1+S_G2_A_1+S_G3_A_1+S_G4_A_1+S_G5_A_1+S_S1_A_1+S_S2_A_1+S_S3_A_1+
            S_S4_A_1+S_S5_A_1
 religio=~  R_G1_A_1+R_G2_A_1+R_G3_A_1+R_G4_A_1+R_G5_A_1+R_S1_A_1+  
            R_S2_A_1+R_S3_A_1+R_S4_A_1+R_S5_A_1
 paranor=~  P_G1_A_1+P_G2_A_1+P_G3_A_1+P_G4_A_1+  
            P_G5_A_1+P_S1_A_1+P_S2_A_1+P_S3_A_1+P_S4_A_1+P_S5_A_1
 politic=~  Pol_S1_A_1+Pol_S2_A_1+
            Pol_S3_A_1+Pol_S4_A_1+Pol_S5_A_1+Pol_G1_A_1+Pol_G2_A_1+Pol_G3_A_1+Pol_G4_A_1+Pol_G5_A_1
 moralit=~  M_G1_A_1+M_G2_A_1+M_G3_A_1+M_G4_A_1+M_G5_r_A_1+M_S1_A_1+M_S2_r_A_1+M_S3_A_1+  
            M_S4_r_A_1+M_S5_A_1 
 general=~  S_G1_A_1+S_G2_A_1+S_G3_A_1+S_G4_A_1+S_G5_A_1+S_S1_A_1+S_S2_A_1+S_S3_A_1+
            S_S4_A_1+S_S5_A_1+R_G1_A_1+R_G2_A_1+R_G3_A_1+R_G4_A_1+R_G5_A_1+R_S1_A_1+  
            R_S2_A_1+R_S3_A_1+R_S4_A_1+R_S5_A_1+P_G1_A_1+P_G2_A_1+P_G3_A_1+P_G4_A_1+  
            P_G5_A_1+P_S1_A_1+P_S2_A_1+P_S3_A_1+P_S4_A_1+P_S5_A_1+Pol_S1_A_1+Pol_S2_A_1+
            Pol_S3_A_1+Pol_S4_A_1+Pol_S5_A_1+Pol_G1_A_1+Pol_G2_A_1+Pol_G3_A_1+Pol_G4_A_1+Pol_G5_A_1+
            M_G1_A_1+M_G2_A_1+M_G3_A_1+M_G4_A_1+M_G5_r_A_1+M_S1_A_1+M_S2_r_A_1+M_S3_A_1+  
            M_S4_r_A_1+M_S5_A_1 
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.782       0.787
#Tucker-Lewis Index (TLI)                       0.762       0.768
#Robust Comparative Fit Index (CFI)                         0.794
#Robust Tucker-Lewis Index (TLI)                            0.776
#RMSEA                                          0.059       0.053
#Robust RMSEA                                               0.057
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .28

# loadings based on Std.all because Std.lv is larger than 1
#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  science =~                                                            
#    S_G1_A_1          0.566    0.200    2.826    0.005    0.566    0.192
#    S_G2_A_1         -1.124    0.139   -8.092    0.000   -1.124   -0.462
#    S_G3_A_1         -1.128    0.149   -7.546    0.000   -1.128   -0.433
#    S_G4_A_1         -1.508    0.148  -10.195    0.000   -1.508   -0.558
#    S_G5_A_1         -1.070    0.119   -8.975    0.000   -1.070   -0.567
#    S_S1_A_1         -0.417    0.086   -4.847    0.000   -0.417   -0.289
#    S_S2_A_1         -0.653    0.175   -3.724    0.000   -0.653   -0.233
#    S_S3_A_1         -0.718    0.133   -5.411    0.000   -0.718   -0.362
#    S_S4_A_1          0.595    0.162    3.681    0.000    0.595    0.213
#    S_S5_A_1         -1.222    0.156   -7.853    0.000   -1.222   -0.442
#  religio =~                                                            
#    R_G1_A_1         -1.457    0.523   -2.786    0.005   -1.457   -0.427
#    R_G2_A_1          0.603    0.523    1.153    0.249    0.603    0.201
#    R_G3_A_1         -1.386    0.526   -2.635    0.008   -1.386   -0.410
#    R_G4_A_1         -0.469    0.580   -0.809    0.419   -0.469   -0.138
#    R_G5_A_1         -0.660    0.429   -1.539    0.124   -0.660   -0.205
#    R_S1_A_1          0.977    0.483    2.022    0.043    0.977    0.326
#    R_S2_A_1          0.490    0.535    0.916    0.360    0.490    0.165
#    R_S3_A_1          0.683    0.293    2.330    0.020    0.683    0.246
#    R_S4_A_1          0.886    0.400    2.213    0.027    0.886    0.313
#    R_S5_A_1          0.562    0.327    1.718    0.086    0.562    0.170
#  paranor =~                                                            
#    P_G1_A_1          1.770    0.225    7.851    0.000    1.770    0.588
#    P_G2_A_1          1.547    0.138   11.227    0.000    1.547    0.608
#    P_G3_A_1          1.286    0.140    9.167    0.000    1.286    0.518
#    P_G4_A_1          1.408    0.144    9.761    0.000    1.408    0.549
#    P_G5_A_1          1.177    0.172    6.846    0.000    1.177    0.394
#    P_S1_A_1          1.808    0.145   12.474    0.000    1.808    0.633
#    P_S2_A_1          1.018    0.152    6.681    0.000    1.018    0.381
#    P_S3_A_1          1.578    0.115   13.762    0.000    1.578    0.693
#    P_S4_A_1          1.413    0.127   11.104    0.000    1.413    0.565
#    P_S5_A_1          1.821    0.134   13.574    0.000    1.821    0.656
#  politic =~                                                            
#    Pol_S1_A_1        0.877    0.211    4.157    0.000    0.877    0.307
#    Pol_S2_A_1        1.261    0.217    5.818    0.000    1.261    0.545
#    Pol_S3_A_1        0.976    0.178    5.471    0.000    0.976    0.375
#    Pol_S4_A_1        0.610    0.154    3.971    0.000    0.610    0.427
#    Pol_S5_A_1       -0.564    0.332   -1.697    0.090   -0.564   -0.174
#    Pol_G1_A_1        0.957    0.179    5.331    0.000    0.957    0.484
#    Pol_G2_A_1        0.671    0.151    4.448    0.000    0.671    0.292
#    Pol_G3_A_1       -0.479    0.180   -2.665    0.008   -0.479   -0.194
#    Pol_G4_A_1        0.672    0.370    1.816    0.069    0.672    0.265
#    Pol_G5_A_1       -0.659    0.449   -1.468    0.142   -0.659   -0.229
#  moralit =~                                                            
#    M_G1_A_1          0.667    0.139    4.787    0.000    0.667    0.538
#    M_G2_A_1          0.943    0.172    5.473    0.000    0.943    0.519
#    M_G3_A_1          1.076    0.155    6.953    0.000    1.076    0.518
#    M_G4_A_1          0.514    0.143    3.598    0.000    0.514    0.191
#    M_G5_r_A_1       -0.325    0.114   -2.848    0.004   -0.325   -0.203
#    M_S1_A_1          1.246    0.158    7.861    0.000    1.246    0.493
#    M_S2_r_A_1       -0.684    0.173   -3.954    0.000   -0.684   -0.249
#    M_S3_A_1          0.587    0.112    5.227    0.000    0.587    0.440
#    M_S4_r_A_1       -0.507    0.196   -2.581    0.010   -0.507   -0.167
#    M_S5_A_1          1.261    0.177    7.129    0.000    1.261    0.458
#  general =~                                                            
#    S_G1_A_1          0.783    0.137    5.697    0.000    0.783    0.266
#    S_G2_A_1         -0.824    0.118   -6.993    0.000   -0.824   -0.339
#    S_G3_A_1         -0.475    0.127   -3.744    0.000   -0.475   -0.182
#    S_G4_A_1         -1.003    0.125   -8.042    0.000   -1.003   -0.371
#    S_G5_A_1         -0.706    0.098   -7.217    0.000   -0.706   -0.374
#    S_S1_A_1         -0.075    0.065   -1.152    0.249   -0.075   -0.052
#    S_S2_A_1         -0.747    0.146   -5.113    0.000   -0.747   -0.267
#    S_S3_A_1         -0.234    0.098   -2.379    0.017   -0.234   -0.118
#    S_S4_A_1          1.101    0.148    7.461    0.000    1.101    0.395
#    S_S5_A_1         -1.559    0.189   -8.241    0.000   -1.559   -0.563
#    R_G1_A_1          2.566    0.337    7.625    0.000    2.566    0.751
#    R_G2_A_1          2.576    0.194   13.268    0.000    2.576    0.858
#    R_G3_A_1          2.358    0.275    8.560    0.000    2.358    0.697
#    R_G4_A_1         -2.670    0.175  -15.279    0.000   -2.670   -0.788
#    R_G5_A_1          2.200    0.184   11.949    0.000    2.200    0.682
#    R_S1_A_1          2.345    0.254    9.242    0.000    2.345    0.782
#    R_S2_A_1          2.451    0.154   15.907    0.000    2.451    0.825
#    R_S3_A_1          1.724    0.204    8.459    0.000    1.724    0.620
#    R_S4_A_1          2.201    0.233    9.466    0.000    2.201    0.777
#    R_S5_A_1          1.620    0.186    8.708    0.000    1.620    0.491
#    P_G1_A_1          1.303    0.328    3.973    0.000    1.303    0.433
#    P_G2_A_1          0.992    0.192    5.163    0.000    0.992    0.390
#    P_G3_A_1          0.915    0.157    5.821    0.000    0.915    0.368
#    P_G4_A_1          1.171    0.161    7.264    0.000    1.171    0.456
#    P_G5_A_1          1.207    0.163    7.397    0.000    1.207    0.404
#    P_S1_A_1          1.026    0.205    5.005    0.000    1.026    0.359
#    P_S2_A_1          0.603    0.162    3.726    0.000    0.603    0.226
#    P_S3_A_1          0.779    0.178    4.369    0.000    0.779    0.342
#    P_S4_A_1          0.931    0.149    6.254    0.000    0.931    0.373
#    P_S5_A_1          1.008    0.221    4.558    0.000    1.008    0.363
#    Pol_S1_A_1       -1.127    0.165   -6.819    0.000   -1.127   -0.395
#    Pol_S2_A_1       -0.099    0.113   -0.876    0.381   -0.099   -0.043
#    Pol_S3_A_1       -1.169    0.246   -4.753    0.000   -1.169   -0.449
#    Pol_S4_A_1        0.033    0.064    0.523    0.601    0.033    0.023
#    Pol_S5_A_1        1.253    0.144    8.728    0.000    1.253    0.387
#    Pol_G1_A_1        0.000    0.091    0.005    0.996    0.000    0.000
#    Pol_G2_A_1       -0.188    0.110   -1.712    0.087   -0.188   -0.082
#    Pol_G3_A_1       -0.122    0.131   -0.932    0.351   -0.122   -0.050
#    Pol_G4_A_1       -0.293    0.160   -1.835    0.067   -0.293   -0.115
#    Pol_G5_A_1        0.645    0.186    3.475    0.001    0.645    0.224
#    M_G1_A_1          0.020    0.050    0.391    0.696    0.020    0.016
#    M_G2_A_1         -0.042    0.079   -0.530    0.596   -0.042   -0.023
#    M_G3_A_1          0.099    0.103    0.962    0.336    0.099    0.048
#    M_G4_A_1         -0.435    0.135   -3.225    0.001   -0.435   -0.162
#    M_G5_r_A_1       -0.121    0.086   -1.412    0.158   -0.121   -0.076
#    M_S1_A_1         -0.297    0.111   -2.685    0.007   -0.297   -0.118
#    M_S2_r_A_1        0.389    0.126    3.075    0.002    0.389    0.141
#    M_S3_A_1          0.048    0.057    0.855    0.393    0.048    0.036
#    M_S4_r_A_1        0.253    0.149    1.700    0.089    0.253    0.083
#    M_S5_A_1         -0.395    0.126   -3.124    0.002   -0.395   -0.143

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.5033857      0.8163265      0.7540430      0.3763234 

#$FactorLevelIndices
#        ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#science 0.5989012 0.09302828 0.40109880 0.6037125 0.464980534 0.6726081 0.8373685
#religio 0.1251526 0.04545537 0.87484740 0.8944758 0.001595108 0.4676771 0.8393197
#paranor 0.6945722 0.18941990 0.30542783 0.8933358 0.619397350 0.8356033 0.9234355
#politic 0.6774279 0.07202356 0.32257207 0.3617176 0.342332249 0.6017179 0.7820289
#moralit 0.9438864 0.09668722 0.05611358 0.4394255 0.436790494 0.6820887 0.8268743
#general 0.5033857 0.50338567 0.50338567 0.7540430 0.376323375 0.9459016 0.9741121



