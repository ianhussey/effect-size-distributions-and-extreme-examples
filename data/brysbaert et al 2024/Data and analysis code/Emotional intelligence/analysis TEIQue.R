################################################################
### Analysis TEIQue Emotional intelligence and other EI tests
### Radeva https://osf.io/snc56/?view_only=

library(readxl)
TEIQue_Radeva <- read_excel("TEIQue_Radeva.xlsx")
colnames(TEIQue_Radeva <- read_excel("TEIQue_Radeva.xlsx"))
mydata <- as.data.frame(TEIQue_Radeva <- read_excel("TEIQue_Radeva.xlsx")[,2:31])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 components
# Eigenvalue 1 = 6.79; eigenvalue 2 = 1.44

rho <- polychoric(mydata)$rho
# does not calculate, because not all question had answers from 1-7
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 3 components
# Eigenvalue 1 = 7.84; eigenvalue 2 = 1.61

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.23, RMSEA=.087, RMSR=.1, TLI=.591

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.117, RMSR=.11, TLI=.499

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities 

# Give solution with 4 factors
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.061, RMSR=.06, TLI=.794
#      MR1   MR4   MR2   MR3
#MR1  1.00 -0.45 -0.34  0.15
#MR4 -0.45  1.00  0.26 -0.12
#MR2 -0.34  0.26  1.00 -0.09
#MR3  0.15 -0.12 -0.09  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.097, RMSR=.06, TLI=.648
#      MR1   MR4   MR2   MR3
#MR1  1.00 -0.49 -0.39  0.15
#MR4 -0.49  1.00  0.23 -0.15
#MR2 -0.39  0.23  1.00 -0.10
#MR3  0.15 -0.15 -0.10  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.915       0.803
#Tucker-Lewis Index (TLI)                       0.909       0.789
#Robust Comparative Fit Index (CFI)                         0.542
#Robust Tucker-Lewis Index (TLI)                            0.508
#RMSEA                                          0.114       0.096
#Robust RMSEA                                               0.128
#SRMR                                           0.106       0.106

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .326

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.599       0.601
#Tucker-Lewis Index (TLI)                       0.570       0.572
#Robust Comparative Fit Index (CFI)                         0.611
#Robust Tucker-Lewis Index (TLI)                            0.582
#RMSEA                                          0.096       0.090
#Robust RMSEA                                               0.093
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .242

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
 wellbe =~ response5+response20+response9+response24+response12+response27
 selfco =~ response4+response19+response7+response22+response15+response30
 emotio =~ response1+response16+response2+response17+response8+response23+response13+response28
 socia  =~ response6+response21+response10+response25+response11+response26+
           response3+response18+response14+response29
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.939       0.847
#Tucker-Lewis Index (TLI)                       0.933       0.834
#Robust Comparative Fit Index (CFI)                         0.620
#Robust Tucker-Lewis Index (TLI)                            0.585
#RMSEA                                          0.097       0.085
#Robust RMSEA                                               0.117
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .379

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  wellbe ~~                                                             
#    selfco            0.683    0.058   11.738    0.000    0.683    0.683
#    emotio           -0.667    0.057  -11.743    0.000   -0.667   -0.667
#    socia            -0.810    0.026  -30.968    0.000   -0.810   -0.810
#  selfco ~~                                                             
#    emotio           -0.496    0.074   -6.689    0.000   -0.496   -0.496
#    socia            -0.698    0.047  -14.724    0.000   -0.698   -0.698
#  emotio ~~                                                             
#    socia             0.677    0.056   12.137    0.000    0.677    0.677

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.701       0.700
#Tucker-Lewis Index (TLI)                       0.674       0.673
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.685
#RMSEA                                          0.084       0.079
#Robust RMSEA                                               0.081
#SRMR                                           0.092       0.092

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#wellbe ~~                                                             
#  selfco            0.649    0.087    7.484    0.000    0.649    0.649
#  emotio           -0.616    0.098   -6.274    0.000   -0.616   -0.616
#  socia            -0.832    0.084   -9.919    0.000   -0.832   -0.832
#selfco ~~                                                             
#  emotio           -0.345    0.127   -2.709    0.007   -0.345   -0.345
#  socia            -0.668    0.078   -8.541    0.000   -0.668   -0.668
#emotio ~~                                                             
#  socia             0.608    0.090    6.782    0.000    0.608    0.608

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .287

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.566       0.558
#Tucker-Lewis Index (TLI)                       0.534       0.525
#Robust Comparative Fit Index (CFI)                         0.572
#Robust Tucker-Lewis Index (TLI)                            0.541
#RMSEA                                          0.100       0.095
#Robust RMSEA                                               0.098
#SRMR                                           0.193       0.193

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .296

# Bifactor model
BIFmodel= '
 wellbe =~ response5+response20+response9+response24+response12+response27
 selfco =~ response4+response19+response7+response22+response15+response30
 emotio =~ response1+response16+response2+response17+response8+response23+response13+response28
 socia  =~ response6+response21+response10+response25+response11+response26+
           response3+response18+response14+response29
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.779       0.780
#Tucker-Lewis Index (TLI)                       0.744       0.744
#Robust Comparative Fit Index (CFI)                         0.789
#Robust Tucker-Lewis Index (TLI)                            0.755
#RMSEA                                          0.074       0.069
#Robust RMSEA                                               0.071
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .375

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  wellbe =~                                                             
#    response5         0.846    0.141    5.980    0.000    0.846    0.603
#    response20       -0.871    0.149   -5.851    0.000   -0.871   -0.649
#    response9        -0.227    0.170   -1.336    0.182   -0.227   -0.216
#    response24       -0.162    0.178   -0.913    0.361   -0.162   -0.135
#    response12        0.619    0.176    3.519    0.000    0.619    0.367
#    response27       -0.232    0.217   -1.070    0.285   -0.232   -0.177
#  selfco =~                                                             
#    response4         1.045    0.206    5.073    0.000    1.045    0.621
#    response19       -0.937    0.209   -4.495    0.000   -0.937   -0.608
#    response7         0.564    0.197    2.858    0.004    0.564    0.345
#    response22        0.412    0.203    2.032    0.042    0.412    0.263
#    response15       -0.459    0.145   -3.162    0.002   -0.459   -0.312
#    response30       -0.181    0.199   -0.909    0.363   -0.181   -0.104
#  emotio =~                                                             
#    response1         0.661    0.192    3.452    0.001    0.661    0.398
#    response16       -1.204    0.170   -7.064    0.000   -1.204   -0.682
#    response2         0.040    0.155    0.257    0.797    0.040    0.033
#    response17        0.185    0.192    0.963    0.336    0.185    0.131
#    response8        -0.474    0.172   -2.754    0.006   -0.474   -0.331
#    response23        0.200    0.202    0.988    0.323    0.200    0.123
#    response13       -0.652    0.191   -3.408    0.001   -0.652   -0.515
#    response28       -0.945    0.153   -6.188    0.000   -0.945   -0.582
#  socia =~                                                              
#    response6         0.122    0.152    0.806    0.420    0.122    0.106
#    response21        0.323    0.225    1.435    0.151    0.323    0.232
#    response10       -0.090    0.187   -0.480    0.631   -0.090   -0.057
#    response25       -0.112    0.264   -0.425    0.671   -0.112   -0.063
#    response11        0.295    0.246    1.197    0.231    0.295    0.239
#    response26       -0.289    0.260   -1.113    0.266   -0.289   -0.207
#    response3        -0.769    0.188   -4.083    0.000   -0.769   -0.616
#    response18        0.944    0.252    3.739    0.000    0.944    0.606
#    response14        0.264    0.171    1.546    0.122    0.264    0.171
#    response29       -0.105    0.117   -0.902    0.367   -0.105   -0.089
#  general =~                                                            
#    response1         0.506    0.170    2.980    0.003    0.506    0.305
#    response2        -0.260    0.134   -1.943    0.052   -0.260   -0.214
#    response3         0.607    0.150    4.048    0.000    0.607    0.487
#    response4        -0.773    0.144   -5.375    0.000   -0.773   -0.459
#    response5        -0.715    0.137   -5.217    0.000   -0.715   -0.510
#    response6         0.671    0.093    7.249    0.000    0.671    0.581
#    response7        -0.400    0.155   -2.589    0.010   -0.400   -0.245
#    response8        -0.773    0.114   -6.761    0.000   -0.773   -0.540
#    response9         0.680    0.089    7.623    0.000    0.680    0.648
#    response10       -0.635    0.135   -4.689    0.000   -0.635   -0.401
#    response11        0.736    0.119    6.183    0.000    0.736    0.596
#    response12       -0.229    0.189   -1.209    0.227   -0.229   -0.136
#    response13       -0.226    0.114   -1.985    0.047   -0.226   -0.179
#    response14       -0.836    0.119   -7.052    0.000   -0.836   -0.541
#    response15        0.852    0.113    7.511    0.000    0.852    0.580
#    response16       -0.479    0.160   -2.998    0.003   -0.479   -0.271
#    response17        0.152    0.134    1.132    0.258    0.152    0.108
#    response18       -0.861    0.183   -4.701    0.000   -0.861   -0.553
#    response19        0.698    0.124    5.651    0.000    0.698    0.453
#    response20        0.811    0.112    7.274    0.000    0.811    0.604
#    response21        0.835    0.118    7.100    0.000    0.835    0.599
#    response22       -0.480    0.161   -2.980    0.003   -0.480   -0.306
#    response23       -0.176    0.159   -1.101    0.271   -0.176   -0.108
#    response24        0.891    0.084   10.592    0.000    0.891    0.744
#    response25       -0.401    0.189   -2.128    0.033   -0.401   -0.225
#    response26       -0.553    0.141   -3.934    0.000   -0.553   -0.396
#    response27        0.798    0.111    7.204    0.000    0.798    0.608
#    response28       -0.733    0.159   -4.595    0.000   -0.733   -0.452
#    response29        0.626    0.092    6.843    0.000    0.626    0.529
#    response30        0.596    0.181    3.293    0.001    0.596    0.344

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.5895472      0.7632184      0.2037233      0.1122320 

#$FactorLevelIndices
#        ECV_SS     ECV_SG    ECV_GS      Omega      OmegaH         H        FD
#wellbe  0.3382038 0.09487259 0.6617962 0.56388440 0.006218344 0.6083904 0.8342080
#selfco  0.5063362 0.09819065 0.4936638 0.04288430 0.010087624 0.6063430 0.8072184
#emotio  0.6425100 0.12805977 0.3574900 0.39658139 0.208738324 0.6761731 0.8349838
#socia   0.2738117 0.08932981 0.7261883 0.07931499 0.014665893 0.5852759 0.8325685
#general 0.5895472 0.58954717 0.5895472 0.20372333 0.112231989 0.9037457 0.9424180









################################################################
### Pérez Díaz, Pablo (2019), “Datasets Spanish-Chilean-TEIQUE-SF”
### https://data.mendeley.com/datasets/gstrr47vpr/3

TEIQue_Perez.Diaz <- read.csv("TEIQue_Perez-Diaz.csv")
colnames(TEIQue_Perez.Diaz)
mydata <- as.data.frame(TEIQue_Perez.Diaz[,1:30])
colnames(mydata) <- c("response1","response2","response3","response4","response5","response6",
                      "response7","response8","response9","response10","response11","response12",
                      "response13","response14","response15","response16","response17","response18",
                      "response19","response20","response21","response22","response23","response24",
                      "response25","response26","response27","response28","response29","response30" )
  
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 6.53; eigenvalue 2 = 2.10

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 7.88; eigenvalue 2 = 2.19

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.08, RMSR=.09, TLI=.639

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.099, RMSR=.1, TLI=.601

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities 

# Give solution with 4 factors
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.041, RMSR=.04, TLI=.905
#      MR1  MR2  MR3   MR4
#MR1  1.00 0.40 0.23 -0.14
#MR2  0.40 1.00 0.25  0.14
#MR3  0.23 0.25 1.00  0.06
#MR4 -0.14 0.14 0.06  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.065, RMSR=.04, TLI=.83
#     MR1   MR2   MR3   MR4
#MR1 1.00  0.42  0.25  0.18
#MR2 0.42  1.00  0.29 -0.07
#MR3 0.25  0.29  1.00 -0.01
#MR4 0.18 -0.07 -0.01  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.913       0.805
#Tucker-Lewis Index (TLI)                       0.906       0.790
#Robust Comparative Fit Index (CFI)                         0.631
#Robust Tucker-Lewis Index (TLI)                            0.604
#RMSEA                                          0.105       0.092
#Robust RMSEA                                               0.102
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .256

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.660       0.676
#Tucker-Lewis Index (TLI)                       0.635       0.652
#Robust Comparative Fit Index (CFI)                         0.683
#Robust Tucker-Lewis Index (TLI)                            0.659
#RMSEA                                          0.082       0.070
#Robust RMSEA                                               0.078
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .195

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
 wellbe =~ response5+response20+response9+response24+response12+response27
 selfco =~ response4+response19+response7+response22+response15+response30
 emotio =~ response1+response16+response2+response17+response8+response23+response13+response28
 socia  =~ response6+response21+response10+response25+response11+response26+
           response3+response18+response14+response29
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.918       0.816
#Tucker-Lewis Index (TLI)                       0.911       0.799
#Robust Comparative Fit Index (CFI)                         0.658
#Robust Tucker-Lewis Index (TLI)                            0.627
#RMSEA                                          0.103       0.090
#Robust RMSEA                                               0.099
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .274

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  wellbe ~~                                                             
#    selfco            0.728    0.039   18.881    0.000    0.728    0.728
#    emotio            0.804    0.031   25.895    0.000    0.804    0.804
#    socia             0.943    0.016   60.475    0.000    0.943    0.943
#  selfco ~~                                                             
#    emotio            0.871    0.030   28.586    0.000    0.871    0.871
#    socia             0.898    0.033   26.823    0.000    0.898    0.898
#  emotio ~~                                                             
#    socia             0.971    0.024   40.549    0.000    0.971    0.971


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.694       0.710
#Tucker-Lewis Index (TLI)                       0.666       0.684
#Robust Comparative Fit Index (CFI)                         0.717
#Robust Tucker-Lewis Index (TLI)                            0.691
#RMSEA                                          0.078       0.066
#Robust RMSEA                                               0.074
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .221

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  wellbe ~~                                                             
#    selfco            0.654    0.078    8.385    0.000    0.654    0.654
#    emotio            0.715    0.063   11.294    0.000    0.715    0.715
#    socia             0.926    0.033   28.071    0.000    0.926    0.926
#  selfco ~~                                                             
#    emotio            0.872    0.062   13.981    0.000    0.872    0.872
#    socia             0.874    0.062   14.198    0.000    0.874    0.874
#  emotio ~~                                                             
#    socia             0.953    0.042   22.503    0.000    0.953    0.953

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.467       0.467
#Tucker-Lewis Index (TLI)                       0.428       0.428
#Robust Comparative Fit Index (CFI)                         0.480
#Robust Tucker-Lewis Index (TLI)                            0.442
#RMSEA                                          0.103       0.089
#Robust RMSEA                                               0.099
#SRMR                                           0.197       0.197

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .221

# Bifactor model
BIFmodel= '
 wellbe =~ response5+response20+response9+response24+response12+response27
 selfco =~ response4+response19+response7+response22+response15+response30
 emotio =~ response1+response16+response2+response17+response8+response23+response13+response28
 socia  =~ response6+response21+response10+response25+response11+response26+
           response3+response18+response14+response29
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .325

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  wellbe =~                                                             
#    response5         0.152    0.323    0.471    0.638    0.152    0.074
#    response20       -0.716    0.251   -2.846    0.004   -0.716   -0.431
#    response9        -0.687    0.155   -4.436    0.000   -0.687   -0.490
#    response24       -0.760    0.173   -4.406    0.000   -0.760   -0.524
#    response12       -0.139    0.249   -0.560    0.575   -0.139   -0.073
#    response27       -0.870    0.235   -3.707    0.000   -0.870   -0.554
#  selfco =~                                                             
#    response4        -0.033    1.135   -0.029    0.977   -0.033   -0.017
#    response19       -0.472    1.306   -0.362    0.718   -0.472   -0.296
#    response7         0.869    1.301    0.668    0.504    0.869    0.451
#    response22        0.854    1.200    0.712    0.476    0.854    0.452
#    response15       -0.371    1.558   -0.238    0.812   -0.371   -0.208
#    response30       -0.371    0.831   -0.447    0.655   -0.371   -0.201
#  emotio =~                                                             
#    response1        -0.758    0.384   -1.977    0.048   -0.758   -0.408
#    response16       -0.009    0.463   -0.019    0.985   -0.009   -0.004
#    response2         0.322    0.337    0.956    0.339    0.322    0.162
#    response17       -0.325    0.350   -0.929    0.353   -0.325   -0.187
#    response8         0.591    0.339    1.744    0.081    0.591    0.286
#    response23       -0.683    0.222   -3.074    0.002   -0.683   -0.410
#    response13        0.041    0.415    0.099    0.921    0.041    0.022
#    response28        0.355    0.448    0.792    0.428    0.355    0.176
#  socia =~                                                              
#    response6        -0.386    0.314   -1.229    0.219   -0.386   -0.252
#    response21       -0.407    0.298   -1.368    0.171   -0.407   -0.290
#    response10        0.947    0.269    3.525    0.000    0.947    0.465
#    response25        0.970    0.428    2.268    0.023    0.970    0.496
#    response11       -0.196    0.223   -0.879    0.379   -0.196   -0.118
#    response26        0.530    0.478    1.109    0.268    0.530    0.287
#    response3        -0.391    0.272   -1.441    0.150   -0.391   -0.269
#    response18       -0.049    0.315   -0.155    0.877   -0.049   -0.025
#    response14        0.153    0.315    0.484    0.628    0.153    0.079
#    response29       -0.239    0.374   -0.641    0.522   -0.239   -0.147
#  general =~                                                            
#    response1         0.802    0.205    3.904    0.000    0.802    0.431
#    response2         0.604    0.180    3.359    0.001    0.604    0.304
#    response3         0.645    0.198    3.264    0.001    0.645    0.444
#    response4         0.909    0.178    5.117    0.000    0.909    0.467
#    response5         1.138    0.115    9.909    0.000    1.138    0.556
#    response6         0.584    0.215    2.719    0.007    0.584    0.381
#    response7         0.979    0.325    3.014    0.003    0.979    0.508
#    response8         1.077    0.243    4.436    0.000    1.077    0.520
#    response9         0.616    0.190    3.249    0.001    0.616    0.439
#    response10        0.898    0.334    2.686    0.007    0.898    0.441
#    response11        0.153    0.144    1.060    0.289    0.153    0.092
#    response12        1.167    0.090   12.925    0.000    1.167    0.615
#    response13        1.002    0.104    9.673    0.000    1.002    0.542
#    response14        1.127    0.099   11.345    0.000    1.127    0.585
#    response15        0.650    0.128    5.078    0.000    0.650    0.364
#    response16        0.864    0.130    6.659    0.000    0.864    0.419
#    response17        0.670    0.175    3.839    0.000    0.670    0.385
#    response18        1.205    0.094   12.800    0.000    1.205    0.629
#    response19        0.532    0.140    3.804    0.000    0.532    0.333
#    response20        0.938    0.212    4.427    0.000    0.938    0.564
#    response21        0.765    0.161    4.747    0.000    0.765    0.544
#    response22        0.740    0.321    2.307    0.021    0.740    0.391
#    response23        0.359    0.224    1.600    0.110    0.359    0.215
#    response24        0.786    0.204    3.857    0.000    0.786    0.542
#    response25        0.272    0.311    0.875    0.382    0.272    0.139
#    response26        0.749    0.144    5.187    0.000    0.749    0.405
#    response27        0.842    0.219    3.841    0.000    0.842    0.536
#    response28        1.196    0.153    7.807    0.000    1.196    0.593
#    response29        0.910    0.146    6.248    0.000    0.910    0.559
#    response30        0.592    0.198    2.987    0.003    0.592    0.320

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.6871222      0.7632184      0.8969846      0.8760572 

#$FactorLevelIndices
#       ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#wellbe  0.3637555 0.10931650 0.6362445 0.8197348 0.224486541 0.5789116 0.7971290
#selfco  0.3723311 0.06220231 0.6276689 0.5626326 0.003221998 0.4100837 0.6856938
#emotio  0.2452620 0.05463696 0.7547380 0.6648823 0.007443399 0.3693943 0.6469816
#socia   0.2800681 0.08672207 0.7199319 0.7149288 0.002022353 0.4930864 0.7362052
#general 0.6871222 0.68712215 0.6871222 0.8969846 0.876057184 0.8980256 0.9463044












################################################################
### Perazzo, Matheus (2020), “Dataset Brazilian TEIQue-SF”
### https://data.mendeley.com/datasets/jjzrpr659v/3

library(haven)
TEIQue_Perazzo <- read_sav("TEIQue_Perazzo.sav")
colnames(TEIQue_Perazzo)
mydata <- as.data.frame(TEIQue_Perazzo[,4:33])
colnames(mydata)
colnames(mydata) <- c("response1","response2","response3","response4","response5","response6",
                      "response7","response8","response9","response10","response11","response12",
                      "response13","response14","response15","response16","response17","response18",
                      "response19","response20","response21","response22","response23","response24",
                      "response25","response26","response27","response28","response29","response30" )

mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 5 components
# Eigenvalue 1 = 6.91; eigenvalue 2 = 1.3

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components
# Eigenvalue 1 = 7.82; eigenvalue 2 = 1.47

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.23, RMSEA=.093, RMSR=.08, TLI=.88

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.109, RMSR=.09, TLI=.565

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities 

# Give solution with 4 factors
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.071, RMSR=.05, TLI=.762
#     MR1  MR3  MR2  MR4
#MR1 1.00 0.40 0.31 0.38
#MR3 0.40 1.00 0.32 0.18
#MR2 0.31 0.32 1.00 0.22
#MR4 0.38 0.18 0.22 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.4, RMSEA=.087, RMSR=.05, TLI=.717
#     MR4  MR1  MR2  MR3
#MR4 1.00 0.42 0.31 0.34
#MR1 0.42 1.00 0.33 0.18
#MR2 0.31 0.33 1.00 0.16
#MR3 0.34 0.18 0.16 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.911       0.785
#Tucker-Lewis Index (TLI)                       0.904       0.769
#Robust Comparative Fit Index (CFI)                         0.585
#Robust Tucker-Lewis Index (TLI)                            0.554
#RMSEA                                          0.113       0.106
#Robust RMSEA                                               0.112
#SRMR                                           0.092       0.092

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .198

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.623       0.622
#Tucker-Lewis Index (TLI)                       0.595       0.594
#Robust Comparative Fit Index (CFI)                         0.629
#Robust Tucker-Lewis Index (TLI)                            0.601
#RMSEA                                          0.095       0.088
#Robust RMSEA                                               0.093
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .174

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
 wellbe =~ response5+response20+response9+response24+response12+response27
 selfco =~ response4+response19+response7+response22+response15+response30
 emotio =~ response1+response16+response2+response17+response8+response23+response13+response28
 socia  =~ response6+response21+response10+response25+response11+response26+
           response3+response18+response14+response29
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.938       0.844
#Tucker-Lewis Index (TLI)                       0.932       0.830
#Robust Comparative Fit Index (CFI)                         0.658
#Robust Tucker-Lewis Index (TLI)                            0.627
#RMSEA                                          0.095       0.091
#Robust RMSEA                                               0.103
#SRMR                                           0.083       0.083

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#wellbe ~~                                                             
# selfco            0.454    0.036   12.565    0.000    0.454    0.454
# emotio            0.702    0.030   23.359    0.000    0.702    0.702
# socia             0.812    0.018   44.650    0.000    0.812    0.812
#selfco ~~                                                             
# emotio            0.552    0.042   13.219    0.000    0.552    0.552
# socia             0.611    0.034   17.901    0.000    0.611    0.611
#emotio ~~                                                             
# socia             0.787    0.027   28.926    0.000    0.787    0.787

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .356

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.705       0.703
#Tucker-Lewis Index (TLI)                       0.679       0.677
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.685
#RMSEA                                          0.084       0.078
#Robust RMSEA                                               0.083
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .234

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  wellbe ~~                                                             
#    selfco            0.445    0.056    7.918    0.000    0.445    0.445
#    emotio            0.724    0.050   14.417    0.000    0.724    0.724
#    socia             0.857    0.025   34.455    0.000    0.857    0.857
#  selfco ~~                                                             
#    emotio            0.493    0.070    7.046    0.000    0.493    0.493
#    socia             0.611    0.049   12.512    0.000    0.611    0.611
#  emotio ~~                                                             
#    socia             0.764    0.056   13.650    0.000    0.764    0.764

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.566       0.555
#Tucker-Lewis Index (TLI)                       0.534       0.522
#Robust Comparative Fit Index (CFI)                         0.570
#Robust Tucker-Lewis Index (TLI)                            0.538
#RMSEA                                          0.101       0.095
#Robust RMSEA                                               0.100
#SRMR                                           0.190       0.190

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .270

# Bifactor model
BIFmodel= '
 wellbe =~ response5+response20+response9+response24+response12+response27
 selfco =~ response4+response19+response7+response22+response15+response30
 emotio =~ response1+response16+response2+response17+response8+response23+response13+response28
 socia  =~ response6+response21+response10+response25+response11+response26+
           response3+response18+response14+response29
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.807       0.797
#Tucker-Lewis Index (TLI)                       0.776       0.765
#Robust Comparative Fit Index (CFI)                         0.810
#Robust Tucker-Lewis Index (TLI)                            0.780
#RMSEA                                          0.070       0.067
#Robust RMSEA                                               0.069
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .334

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6006141      0.7632184      0.9069534      0.8262409 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#wellbe  0.3475008 0.12767362 0.6524992 0.9039253 0.1276372 1.2259580 1.5789120
#selfco  0.6635994 0.11512330 0.3364006 0.6783398 0.4027911 0.6912097 0.8489577
#emotio  0.4966066 0.08710665 0.5033934 0.6630365 0.2973189 0.5630880 0.7577351
#socia   0.2449079 0.06948231 0.7550921 0.7797231 0.1089522 0.5080255 0.7349389
#general 0.6006141 0.60061412 0.6006141 0.9069534 0.8262409 0.9134500 0.9502532












################################################################
### Rua https://osf.io/4shg5/?view_only=
### Test Schutte et al with 33 questions
### CFA is based on A nested version of the Saklofske et al. (2003) four-factor model
### see https://www.researchgate.net/publication/223950298_An_examination_of_the_factor_structure_of_the_Schutte_Self-Report_Emotional_Intelligence_SSREI_Scale_via_confirmatory_factor_analysis

library(readxl)
AES_Rua <- read_excel("AES_Rua.xlsx")
colnames(AES_Rua)
mydata <- as.data.frame(AES_Rua[,5:37])
colnames(mydata)
colnames(mydata) <- c("response1","response2","response3","response4","response5","response6",
                      "response7","response8","response9","response10","response11","response12",
                      "response13","response14","response15","response16","response17","response18",
                      "response19","response20","response21","response22","response23","response24",
                      "response25","response26","response27","response28","response29","response30",
                      "response31","response32","response33" )

mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3factors and 3 components
# Eigenvalue 1 = 9.33; eigenvalue 2 = 2.0

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components
# Eigenvalue 1 = 10.7; eigenvalue 2 = 2.13

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.28, RMSEA=.093, RMSR=.08, TLI=.88

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.257, RMSR=.11, TLI=.165

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities and response bias

# Give solution with 4 factors
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.054, RMSR=.06, TLI=.853
#     MR1  MR3  MR2  MR4
#MR1 1.00 0.42 0.27 0.36
#MR3 0.42 1.00 0.38 0.28
#MR2 0.27 0.38 1.00 0.12
#MR4 0.36 0.28 0.12 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.267, RMSR=.07, TLI=.077
#     MR4  MR1  MR2  MR3
#MR1 1.00 0.49 0.31 0.18
#MR3 0.49 1.00 0.38 0.14
#MR2 0.31 0.38 1.00 0.04
#MR4 0.18 0.14 0.04 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30+response31+response32+response33
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.950       0.857
#Tucker-Lewis Index (TLI)                       0.946       0.847
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.101       0.082
#Robust RMSEA                                                  NA
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .353

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.634       0.627
#Tucker-Lewis Index (TLI)                       0.609       0.603
#Robust Comparative Fit Index (CFI)                         0.639
#Robust Tucker-Lewis Index (TLI)                            0.615
#RMSEA                                          0.100       0.096
#Robust RMSEA                                               0.098
#SRMR                                           0.099       0.099

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .272

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors based on theoretical analysis
EGAmodel= '
 optimi =~ response22+response9+response12+response23+response10+response28+
           response14+response2+response3+response1
 apprais=~ response29+response18+response25+response32+response5+response15
 utiliz =~ response27+response20+response7+response17+response33
 social =~ response30+response11+response4+response24+response13+response16+
           response31
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.969       0.894
#Tucker-Lewis Index (TLI)                       0.966       0.883
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.085       0.080
#Robust RMSEA                                                  NA
#SRMR                                           0.103       0.103

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .458

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  optimi ~~                                                             
#    apprais           0.837    0.036   23.501    0.000    0.837    0.837
#    utiliz            0.727    0.050   14.480    0.000    0.727    0.727
#    social            0.934    0.028   32.917    0.000    0.934    0.934
#  apprais ~~                                                            
#    utiliz            0.551    0.071    7.749    0.000    0.551    0.551
#    social            0.722    0.052   13.976    0.000    0.722    0.722
#  utiliz ~~                                                             
#    social            0.686    0.066   10.347    0.000    0.686    0.686

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.766       0.766
#Tucker-Lewis Index (TLI)                       0.743       0.743
#Robust Comparative Fit Index (CFI)                         0.776
#Robust Tucker-Lewis Index (TLI)                            0.754
#RMSEA                                          0.086       0.080
#Robust RMSEA                                               0.082
#SRMR                                           0.097       0.097

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#optimi ~~                                                             
# apprais           0.780    0.072   10.801    0.000    0.780    0.780
# utiliz            0.686    0.127    5.401    0.000    0.686    0.686
# social            0.927    0.070   13.278    0.000    0.927    0.927
#apprais ~~                                                            
# utiliz            0.416    0.149    2.794    0.005    0.416    0.416
# social            0.705    0.082    8.593    0.000    0.705    0.705
#utiliz ~~                                                             
# social            0.644    0.157    4.088    0.000    0.644    0.644

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .378

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.614       0.611
#Tucker-Lewis Index (TLI)                       0.583       0.580
#Robust Comparative Fit Index (CFI)                         0.623
#Robust Tucker-Lewis Index (TLI)                            0.593
#RMSEA                                          0.109       0.102
#Robust RMSEA                                               0.106
#SRMR                                           0.245       0.245

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .380

# Bifactor model
BIFmodel= '
 optimi =~ response22+response9+response12+response23+response10+response28+
           response14+response2+response3+response1
 apprais=~ response29+response18+response25+response32+response5+response15
 utiliz =~ response27+response20+response7+response17+response33
 social =~ response30+response11+response4+response24+response13+response16+
           response31
 general=~  response1+response2+response3+response4+response5+response6+response7+response8+
            response9+response10+response11+response12+response13+response14+response15+
            response16+response17+response18+response19+response20+response21+response22+
            response23+response24+response25+response26+response27+response28+response29+
            response30+response31+response32+response33
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.759       0.734
#Tucker-Lewis Index (TLI)                       0.728       0.700
#Robust Comparative Fit Index (CFI)                         0.753
#Robust Tucker-Lewis Index (TLI)                            0.721
#RMSEA                                          0.083       0.083
#Robust RMSEA                                               0.083
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .418

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  optimi =~                                                             
#    response22       -0.156    0.156   -0.998    0.318   -0.156   -0.149
#    response9        -0.071    0.167   -0.424    0.671   -0.071   -0.068
#    response12       -0.109    0.162   -0.673    0.501   -0.109   -0.099
#    response23        0.341    0.223    1.531    0.126    0.341    0.333
#    response10        0.337    0.147    2.295    0.022    0.337    0.362
#    response28       -0.477    0.153   -3.124    0.002   -0.477   -0.441
#    response14        0.439    0.211    2.084    0.037    0.439    0.476
#    response2         0.271    0.153    1.774    0.076    0.271    0.300
#    response3         0.396    0.242    1.633    0.102    0.396    0.444
#    response1         0.071    0.140    0.508    0.612    0.071    0.075
#  apprais =~                                                            
#    response29        0.624    0.168    3.719    0.000    0.624    0.580
#    response18        0.513    0.144    3.554    0.000    0.513    0.510
#    response25        0.241    0.171    1.413    0.158    0.241    0.245
#    response32        0.268    0.124    2.155    0.031    0.268    0.256
#    response5        -0.036    0.225   -0.160    0.873   -0.036   -0.030
#    response15        0.483    0.188    2.576    0.010    0.483    0.437
#  utiliz =~                                                             
#    response27        0.347    0.115    3.011    0.003    0.347    0.355
#    response20        0.670    0.091    7.368    0.000    0.670    0.737
#    response7         0.364    0.111    3.281    0.001    0.364    0.381
#    response17        0.652    0.116    5.601    0.000    0.652    0.715
#    response33        0.150    0.119    1.259    0.208    0.150    0.139
#  social =~                                                             
#    response30        0.722    0.174    4.155    0.000    0.722    0.684
#    response11        0.018    0.156    0.114    0.910    0.018    0.015
#    response4         0.259    0.226    1.145    0.252    0.259    0.253
#    response24        0.388    0.159    2.438    0.015    0.388    0.413
#    response13       -0.191    0.191   -1.001    0.317   -0.191   -0.149
#    response16        0.031    0.183    0.171    0.865    0.031    0.031
#    response31        0.405    0.145    2.799    0.005    0.405    0.420
#  general =~                                                            
#    response1         0.458    0.118    3.865    0.000    0.458    0.483
#    response2         0.652    0.124    5.258    0.000    0.652    0.720
#    response3         0.373    0.178    2.097    0.036    0.373    0.419
#    response4         0.225    0.142    1.589    0.112    0.225    0.220
#    response5        -0.330    0.177   -1.868    0.062   -0.330   -0.275
#    response6         0.268    0.187    1.430    0.153    0.268    0.247
#    response7         0.508    0.137    3.701    0.000    0.508    0.532
#    response8         0.459    0.148    3.100    0.002    0.459    0.428
#    response9         0.724    0.115    6.324    0.000    0.724    0.696
#    response10        0.405    0.168    2.409    0.016    0.405    0.435
#    response11        0.693    0.120    5.784    0.000    0.693    0.570
#    response12        0.804    0.107    7.504    0.000    0.804    0.732
#    response13        0.556    0.132    4.205    0.000    0.556    0.435
#    response14        0.569    0.152    3.736    0.000    0.569    0.617
#    response15        0.519    0.160    3.240    0.001    0.519    0.469
#    response16        0.487    0.157    3.099    0.002    0.487    0.476
#    response17        0.382    0.152    2.517    0.012    0.382    0.419
#    response18        0.538    0.117    4.616    0.000    0.538    0.535
#    response19        0.733    0.106    6.906    0.000    0.733    0.693
#    response20        0.378    0.160    2.365    0.018    0.378    0.415
#    response21        0.666    0.109    6.114    0.000    0.666    0.650
#    response22        0.844    0.109    7.769    0.000    0.844    0.806
#    response23        0.689    0.155    4.446    0.000    0.689    0.671
#    response24        0.444    0.163    2.723    0.006    0.444    0.473
#    response25        0.638    0.105    6.064    0.000    0.638    0.648
#    response26        0.236    0.139    1.697    0.090    0.236    0.213
#    response27        0.576    0.123    4.670    0.000    0.576    0.589
#    response28       -0.409    0.173   -2.360    0.018   -0.409   -0.378
#    response29        0.282    0.131    2.147    0.032    0.282    0.262
#    response30        0.491    0.139    3.546    0.000    0.491    0.466
#    response31        0.485    0.140    3.451    0.001    0.485    0.502
#    response32        0.621    0.095    6.519    0.000    0.621    0.594
#    response33       -0.230    0.144   -1.601    0.109   -0.230   -0.214

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.6814135      0.8276515      0.9186264      0.8627728 

#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#optimi  0.2087041 0.07609908 0.7912959 0.8448647 0.04483932 0.5448429 0.7937114
#apprais 0.3910390 0.07013512 0.6089610 0.7103715 0.31595152 0.5515812 0.7730805
#utiliz  0.5676569 0.10316610 0.4323431 0.7625917 0.48882346 0.7198074 0.8838785
#social  0.3785311 0.06918620 0.6214689 0.7324469 0.16079924 0.5818380 0.8151992
#general 0.6814135 0.68141350 0.6814135 0.9186264 0.86277284 0.9368710 0.9638545










