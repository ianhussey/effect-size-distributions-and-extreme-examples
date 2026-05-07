################################################################
### Analysis BISBAS Questionnaire
### Try to find a common factor structure via analysis of different studies

BISBAS_Rutten <- read.csv("BISBAS_Rutten.csv", sep=";")
colnames(BISBAS_Rutten)
mydata  <- as.data.frame(BISBAS_Rutten[,257:276])
mydata <- na.omit(mydata)
print(colnames(mydata))
colnames(mydata) <- c("B2","B3","B4","B5","B7","B8","B9","B10","B12","B13","B14","B15","B16","B18","B19","B20","B21","B22","B23","B24")
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa(mydata,3)
fit1
diagram(fit1)

library(haven)
BISBASdata_Dierickx <- read_sav("BISBASdata_Dierickx.sav")
colnames(BISBASdata_Dierickx)
mydata2  <- as.data.frame(BISBASdata_Dierickx[,c(3:6,8:11,13:17,19:25)])
mydata2 <- na.omit(mydata2)
colnames(mydata2) <- c("B2","B3","B4","B5","B7","B8","B9","B10","B12","B13","B14","B15","B16","B18","B19","B20","B21","B22","B23","B24")

fit2 <- fa(mydata2,3)
fit2
diagram(fit2)


library(haven)
BISBAS_Khaliq <- read_sav("BISBAS_Khaliq.sav")
colnames(BISBAS_Khaliq)
mydata3  <- as.data.frame(BISBAS_Khaliq[,298:321])
mydata3 <- na.omit(mydata3)
mydata3 <- mydata3[-c(1,6,11,17)]
colnames(mydata3) <- c("B2","B3","B4","B5","B7","B8","B9","B10","B12","B13","B14","B15","B16","B18","B19","B20","B21","B22","B23","B24")

fit3 <- fa(mydata3,3)
fit3
diagram(fit3)


BISBAS_Molenda <- read_sav("BISBAS_Molenda.sav")
colnames(BISBAS_Molenda)
mydata4  <- as.data.frame(BISBAS_Molenda[,34:57])
mydata4 <- na.omit(mydata4)
mydata4 <- mydata4[-c(1,6,11,17)]
colnames(mydata4) <- c("B2","B3","B4","B5","B7","B8","B9","B10","B12","B13","B14","B15","B16","B18","B19","B20","B21","B22","B23","B24")

fit4 <- fa(mydata4,3)
fit4
diagram(fit4)


BISBAS_Jonker <- read_sav("BISBAS_Jonker.sav")
colnames(BISBAS_Jonker)
mydata5  <- as.data.frame(BISBAS_Jonker[,23:46])
mydata5 <- na.omit(mydata5)
mydata5 <- mydata5[-c(1,6,11,17)]
colnames(mydata5) <- c("B2","B3","B4","B5","B7","B8","B9","B10","B12","B13","B14","B15","B16","B18","B19","B20","B21","B22","B23","B24")

fit5 <- fa(mydata5,3)
fit5
diagram(fit5)


library(readxl)
BISBAS_Binter <- read_excel("BISBAS_Binter.xlsx")
colnames(BISBAS_Binter)
mydata6  <- as.data.frame(BISBAS_Binter[,10:29])
colnames(mydata6)
colnames(mydata6) <- c("B3","B9","B12","B21","B5","B10","B15","B20","B4","B7","B14","B18","B23","B2","B8","B13","B16","B19","B22","B24")

fit6 <- fa(mydata6,3)
fit6
diagram(fit6)


mydata_all <- rbind(mydata,mydata2,mydata3,mydata4,mydata5,mydata6)
fit7 <- fa.parallel(mydata_all)
plot(fit7,show.legend=TRUE,fa="both")
fit7

fit7 <- fa(mydata_all,3)
fit7
diagram(fit7)


################################################################
### Analysis BISBAS Questionnaire
### Rutten et al Leuven https://osf.io/gawf5

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 3 components
# Eigenvalue 1 = 3.51; eigenvalue 2 = 2.43

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 factors and 3 components
# Eigenvalue 1 = 4.47; eigenvalue 2 = 2.99

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.121, RMSR=.15, TLI=.387

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.22, RMSEA=.197, RMSR=.19, TLI=.24

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities 

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.066, RMSR=.06, TLI=.81
#     MR1  MR2  MR3
#MR1 1.00 0.04 0.12
#MR2 0.04 1.00 0.26
#MR3 0.12 0.26 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.16, RMSR=.07, TLI=.487
#     MR1  MR2  MR3
#MR1 1.00 0.06 0.13
#MR2 0.06 1.00 0.28
#MR3 0.13 0.28 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'
library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.731       0.583
#Tucker-Lewis Index (TLI)                       0.700       0.533
#Robust Comparative Fit Index (CFI)                         0.384
#Robust Tucker-Lewis Index (TLI)                            0.311
#RMSEA                                          0.180       0.148
#Robust RMSEA                                               0.196
#SRMR                                           0.184       0.184

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .218

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.462       0.435
#Tucker-Lewis Index (TLI)                       0.398       0.369
#Robust Comparative Fit Index (CFI)                         0.457
#Robust Tucker-Lewis Index (TLI)                            0.393
#RMSEA                                          0.127       0.125
#Robust RMSEA                                               0.127
#SRMR                                           0.149       0.149

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .056

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.910       0.853
#Tucker-Lewis Index (TLI)                       0.898       0.835
#Robust Comparative Fit Index (CFI)                         0.642
#Robust Tucker-Lewis Index (TLI)                            0.597
#RMSEA                                          0.105       0.088
#Robust RMSEA                                               0.150
#SRMR                                           0.131       0.131

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas              -0.190    0.088   -2.168    0.030   -0.190   -0.190

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .372

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.730       0.723
#Tucker-Lewis Index (TLI)                       0.696       0.688
#Robust Comparative Fit Index (CFI)                         0.732
#Robust Tucker-Lewis Index (TLI)                            0.699
#RMSEA                                          0.091       0.088
#Robust RMSEA                                               0.089
#SRMR                                           0.106       0.106

#Covariances:
#        Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas              -0.158    0.127   -1.242    0.214   -0.158   -0.158

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .326


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.728       0.721
#Tucker-Lewis Index (TLI)                       0.696       0.689
#Robust Comparative Fit Index (CFI)                         0.731
#Robust Tucker-Lewis Index (TLI)                            0.699
#RMSEA                                          0.091       0.088
#Robust RMSEA                                               0.089
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .285

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.841       0.822
#Tucker-Lewis Index (TLI)                       0.799       0.774
#Robust Comparative Fit Index (CFI)                         0.834
#Robust Tucker-Lewis Index (TLI)                            0.790
#RMSEA                                          0.074       0.075
#Robust RMSEA                                               0.074
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .304

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  bis =~                                                                
#    B2                0.593    0.140    4.227    0.000    0.593    0.728
#    B8                0.146    0.149    0.982    0.326    0.146    0.197
#    B13               0.033    0.186    0.176    0.860    0.033    0.041
#    B16              -0.396    0.168   -2.353    0.019   -0.396   -0.480
#    B19              -0.023    0.144   -0.158    0.874   -0.023   -0.031
#    B22               0.334    0.161    2.076    0.038    0.334    0.427
#    B24              -0.090    0.192   -0.466    0.641   -0.090   -0.129
#  bas =~                                                                
#    B3                0.387    0.112    3.453    0.001    0.387    0.462
#    B9                0.341    0.117    2.928    0.003    0.341    0.434
#    B12               0.501    0.100    5.017    0.000    0.501    0.700
#    B21               0.223    0.085    2.634    0.008    0.223    0.340
#    B5                0.346    0.069    5.024    0.000    0.346    0.561
#    B10               0.252    0.090    2.792    0.005    0.252    0.289
#    B15               0.276    0.130    2.115    0.034    0.276    0.289
#    B20               0.315    0.115    2.747    0.006    0.315    0.377
#    B4                0.279    0.074    3.793    0.000    0.279    0.424
#    B7                0.289    0.053    5.442    0.000    0.289    0.467
#    B14               0.366    0.066    5.541    0.000    0.366    0.584
#    B18               0.187    0.067    2.792    0.005    0.187    0.340
#    B23               0.366    0.075    4.846    0.000    0.366    0.515
#  general =~                                                            
#    B2                0.509    0.146    3.500    0.000    0.509    0.626
#    B3               -0.048    0.097   -0.498    0.618   -0.048   -0.058
#    B4               -0.192    0.083   -2.300    0.021   -0.192   -0.291
#    B5                0.046    0.087    0.523    0.601    0.046    0.074
#    B7               -0.162    0.075   -2.156    0.031   -0.162   -0.262
#    B8               -0.596    0.164   -3.644    0.000   -0.596   -0.804
#    B9               -0.159    0.094   -1.699    0.089   -0.159   -0.203
#    B10               0.239    0.099    2.417    0.016    0.239    0.274
#    B12               0.001    0.083    0.013    0.990    0.001    0.001
#    B13              -0.497    0.100   -4.954    0.000   -0.497   -0.627
#    B14              -0.125    0.077   -1.616    0.106   -0.125   -0.199
#    B15              -0.122    0.123   -0.991    0.321   -0.122   -0.128
#    B16              -0.471    0.119   -3.968    0.000   -0.471   -0.572
#    B18              -0.042    0.075   -0.562    0.574   -0.042   -0.077
#    B19              -0.442    0.098   -4.509    0.000   -0.442   -0.599
#    B20               0.142    0.147    0.966    0.334    0.142    0.169
#    B21               0.077    0.086    0.898    0.369    0.077    0.118
#    B22               0.391    0.106    3.681    0.000    0.391    0.500
#    B23              -0.096    0.100   -0.958    0.338   -0.096   -0.135
#    B24              -0.470    0.123   -3.811    0.000   -0.470   -0.675

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.4620601      0.4789474      0.7642932      0.1487028 

#$FactorLevelIndices
#        ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#bis     0.2618863 0.1435590 0.7381137 0.6207758 0.06796232 0.6312662 0.8994678
#bas     0.8728582 0.3943809 0.1271418 0.7751759 0.76350903 0.7963078 0.8943098
#general 0.4620601 0.4620601 0.4620601 0.7642932 0.14870281 0.8521770 0.9260919





################################################################
### Dierickx et al VUB https://osf.io/fv5qg?view_only=d7e4031c26144e729ed1d4084a7acae5

min(mydata2)
max(mydata2) # 4 response alternatives

mydata <- mydata2

library(psych)
fit1 <- fa.parallel(mydata2)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 components
# Eigenvalue 1 = 2.84; eigenvalue 2 = 2.45

rho <- polychoric(mydata2)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata2))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 3.42; eigenvalue 2 = 3.01

# Give solution with 1 factor
fit3 <- fa(mydata2,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.123, RMSR=.14, TLI=.305

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata2))
fit4
diagram(fit4)
# %variance explained=.17, RMSEA=.159, RMSR=.18, TLI=.26

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata2, algorithm="louvain")
EGALV
# suggests 3 communities and response bias

# Give solution with 3 factors
fit4 <- fa(mydata2,3)
fit4
diagram(fit4)
# %variance explained=.33, RMSEA=.057, RMSR=.04, TLI=.849
#      MR2   MR1   MR3
#MR2  1.00 -0.02  0.12
#MR1 -0.02  1.00 -0.19
#MR3  0.12 -0.19  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata2))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.093, RMSR=.05, TLI=.747
#      MR2   MR1   MR3
#MR2  1.00 -0.01 -0.05
#MR1 -0.01  1.00  0.18
#MR3 -0.05  0.18  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata2,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.598       0.445
#Tucker-Lewis Index (TLI)                       0.550       0.380
#Robust Comparative Fit Index (CFI)                         0.331
#Robust Tucker-Lewis Index (TLI)                            0.252
#RMSEA                                          0.177       0.152
#Robust RMSEA                                               0.163
#SRMR                                           0.175       0.175

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .239

CFA_model1 <- cfa(UNImodel, data = mydata2,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.413       0.380
#Tucker-Lewis Index (TLI)                       0.344       0.307
#Robust Comparative Fit Index (CFI)                         0.412
#Robust Tucker-Lewis Index (TLI)                            0.343
#RMSEA                                          0.121       0.117
#Robust RMSEA                                               0.121
#SRMR                                           0.141       0.141

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .274

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata2,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.888       0.822
#Tucker-Lewis Index (TLI)                       0.874       0.800
#Robust Comparative Fit Index (CFI)                         0.700
#Robust Tucker-Lewis Index (TLI)                            0.663
#RMSEA                                          0.094       0.086
#Robust RMSEA                                               0.110
#SRMR                                           0.108       0.108

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas               0.080    0.056    1.415    0.157    0.080    0.080


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .380

CFA_model2 <- cfa(EGAmodel,mydata2,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.761       0.762
#Tucker-Lewis Index (TLI)                       0.731       0.733
#Robust Comparative Fit Index (CFI)                         0.768
#Robust Tucker-Lewis Index (TLI)                            0.739
#RMSEA                                          0.078       0.073
#Robust RMSEA                                               0.076
#SRMR                                           0.083       0.083

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas               0.062    0.078    0.794    0.427    0.062    0.062


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .281


CFA_model3 <- cfa(EGAmodel,mydata2,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.761       0.763
#Tucker-Lewis Index (TLI)                       0.733       0.735
#Robust Comparative Fit Index (CFI)                         0.768
#Robust Tucker-Lewis Index (TLI)                            0.741
#RMSEA                                          0.078       0.072
#Robust RMSEA                                               0.076
#SRMR                                           0.084       0.084

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .280

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata2,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.847       0.842
#Tucker-Lewis Index (TLI)                       0.806       0.800
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.810
#RMSEA                                          0.066       0.063
#Robust RMSEA                                               0.065
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .369

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#   0.405046847    0.478947368    0.721780505    0.005613223 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#bis     0.2965946 0.1512186 0.70340542 0.6079722 0.473081955 0.5516214 0.7478747
#bas     0.9053027 0.4437346 0.09469731 0.7585821 0.749156508 0.8000484 0.8964143
#general 0.4050468 0.4050468 0.40504685 0.7217805 0.005613223 0.7933818 0.8801262







################################################################
### Khaliq et al. https://osf.io/6kds9/

min(mydata3)
max(mydata3) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata3)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 3.76; eigenvalue 2 = 2.38

rho <- polychoric(mydata3)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata3))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 4.67; eigenvalue 2 = 2.95

# Give solution with 1 factor
fit3 <- fa(mydata3,1)
fit3
diagram(fit3)
# %variance explained=.19, RMSEA=.126, RMSR=.15, TLI=.418

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata3))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.173, RMSR=.18, TLI=.337

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata3, algorithm="louvain")
EGALV
# suggests 4 communities and response bias

# Give solution with 3 factors
fit4 <- fa(mydata3,3)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.068, RMSR=.05, TLI=.728
#     MR2  MR1  MR3
#MR2 1.00 0.01 0.01
#MR1 0.01 1.00 0.49
#MR3 0.01 0.49 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata3))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.118, RMSR=.06, TLI=.691
#     MR2  MR1  MR3
#MR2 1.00 0.02 0.04
#MR1 0.02 1.00 0.51
#MR3 0.04 0.51 1.00


# Single factor model lavaan
UNImodel= '
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata3,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.648       0.394
#Tucker-Lewis Index (TLI)                       0.606       0.323
#Robust Comparative Fit Index (CFI)                         0.429
#Robust Tucker-Lewis Index (TLI)                            0.362
#RMSEA                                          0.198       0.179
#Robust RMSEA                                               0.172
#SRMR                                           0.175       0.175

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .283

CFA_model1 <- cfa(UNImodel, data = mydata3,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.479       0.464
#Tucker-Lewis Index (TLI)                       0.417       0.401
#Robust Comparative Fit Index (CFI)                         0.480
#Robust Tucker-Lewis Index (TLI)                            0.419
#RMSEA                                          0.128       0.123
#Robust RMSEA                                               0.127
#SRMR                                           0.139       0.139

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .181

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata3,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.871       0.777
#Tucker-Lewis Index (TLI)                       0.855       0.750
#Robust Comparative Fit Index (CFI)                         0.676
#Robust Tucker-Lewis Index (TLI)                            0.636
#RMSEA                                          0.120       0.109
#Robust RMSEA                                               0.130
#SRMR                                           0.126       0.126

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas               0.003    0.063    0.054    0.957    0.003    0.003

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .452

CFA_model2 <- cfa(EGAmodel,mydata3,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.732       0.728
#Tucker-Lewis Index (TLI)                       0.698       0.695
#Robust Comparative Fit Index (CFI)                         0.736
#Robust Tucker-Lewis Index (TLI)                            0.704
#RMSEA                                          0.092       0.088
#Robust RMSEA                                               0.091
#SRMR                                           0.101       0.101

#Covariances:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas                0.017    0.090    0.195    0.846    0.017    0.017

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .335


CFA_model3 <- cfa(EGAmodel,mydata3,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.732       0.730
#Tucker-Lewis Index (TLI)                       0.701       0.698
#Robust Comparative Fit Index (CFI)                         0.737
#Robust Tucker-Lewis Index (TLI)                            0.706
#RMSEA                                          0.092       0.087
#Robust RMSEA                                               0.090
#SRMR                                           0.100       0.100

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .336

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata3,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.840       0.800
#Tucker-Lewis Index (TLI)                       0.797       0.747
#Robust Comparative Fit Index (CFI)                         0.831
#Robust Tucker-Lewis Index (TLI)                            0.786
#RMSEA                                          0.076       0.080
#Robust RMSEA                                               0.077
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .407

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5071841      0.4789474      0.8279309      0.5631701 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#bis     0.9186914 0.3513901 0.08130856 0.8043028 0.80368445 0.8033395 0.8992990
#bas     0.2290258 0.1414258 0.77097417 0.8422324 0.04160273 0.5964280 0.8035274
#general 0.5071841 0.5071841 0.50718412 0.8279309 0.56317010 0.8358281 0.9122165







################################################################
### Molenda et al. https://cyberpsychology.eu/article/view/15360

min(mydata4)
max(mydata4) # 4 response alternatives

library(psych)
fit1 <- fa.parallel(mydata4)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 3 components
# Eigenvalue 1 = 5.07; eigenvalue 2 = 2.53

rho <- polychoric(mydata4)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata4))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 6.01; eigenvalue 2 = 2.97

# Give solution with 1 factor
fit3 <- fa(mydata4,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.135, RMSR=.14, TLI=.506

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata4))
fit4
diagram(fit4)
# %variance explained=.30, RMSEA=.166, RMSR=.17, TLI=.468

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata4, algorithm="louvain")
EGALV
# suggests 4 communities and response bias

# Give solution with 3 factors
fit4 <- fa(mydata4,3)
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.059, RMSR=.03, TLI=.906
#      MR2  MR1   MR3
#MR2  1.00 0.25 -0.01
#MR1  0.25 1.00  0.46
#MR3 -0.01 0.46  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata4))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.081, RMSR=.04, TLI=.873
#      MR2  MR1   MR3
#MR2  1.00 0.25 -0.01
#MR1  0.25 1.00  0.46
#MR3 -0.01 0.46  1.00

# Single factor model lavaan
UNImodel= '
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata4,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.787       0.282
#Tucker-Lewis Index (TLI)                       0.762       0.198
#Robust Comparative Fit Index (CFI)                         0.519
#Robust Tucker-Lewis Index (TLI)                            0.463
#RMSEA                                          0.216       0.211
#Robust RMSEA                                               0.168
#SRMR                                           0.165       0.165

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .350

CFA_model1 <- cfa(UNImodel, data = mydata4,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.560       0.556
#Tucker-Lewis Index (TLI)                       0.509       0.503
#Robust Comparative Fit Index (CFI)                         0.564
#Robust Tucker-Lewis Index (TLI)                            0.513
#RMSEA                                          0.135       0.115
#Robust RMSEA                                               0.134
#SRMR                                           0.138       0.138

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .261

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata4,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.818       0.341
#Tucker-Lewis Index (TLI)                       0.795       0.259
#Robust Comparative Fit Index (CFI)                         0.622
#Robust Tucker-Lewis Index (TLI)                            0.575
#RMSEA                                          0.200       0.202
#Robust RMSEA                                               0.149
#SRMR                                           0.153       0.153


#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas              -0.627    0.023  -27.441    0.000   -0.627   -0.627

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .407

CFA_model2 <- cfa(EGAmodel,mydata4,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.668       0.667
#Tucker-Lewis Index (TLI)                       0.627       0.626
#Robust Comparative Fit Index (CFI)                         0.673
#Robust Tucker-Lewis Index (TLI)                            0.633
#RMSEA                                          0.118       0.100
#Robust RMSEA                                               0.116
#SRMR                                           0.135       0.135

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas              -0.542    0.054  -10.096    0.000   -0.542   -0.542

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .334


CFA_model3 <- cfa(EGAmodel,mydata4,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.639       0.638
#Tucker-Lewis Index (TLI)                       0.596       0.596
#Robust Comparative Fit Index (CFI)                         0.643
#Robust Tucker-Lewis Index (TLI)                            0.601
#RMSEA                                          0.122       0.104
#Robust RMSEA                                               0.121
#SRMR                                           0.172       0.172

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .340

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
# Not done at the moment, because poor fit of hypothesized model
BIFmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata4,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.855       0.857
#Tucker-Lewis Index (TLI)                       0.816       0.819
#Robust Comparative Fit Index (CFI)                         0.860
#Robust Tucker-Lewis Index (TLI)                            0.823
#RMSEA                                          0.083       0.070
#Robust RMSEA                                               0.081
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .414

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5947459      0.4789474      0.8957320      0.8628916 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#bis     0.6824125 0.2492534 0.3175875 0.7177049 0.233976587 0.7518250 0.8844956
#bas     0.2457684 0.1560007 0.7542316 0.8720560 0.006595933 0.6406828 0.8543429
#general 0.5947459 0.5947459 0.5947459 0.8957320 0.862891564 0.8781333 0.9446563








################################################################
### Jonker et al. https://osf.io/j5x6h/?view_only=

min(mydata5)
max(mydata5) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata5)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 3.18; eigenvalue 2 = 2.55

rho <- polychoric(mydata5)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata6))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 4.04; eigenvalue 2 = 3.11

# Give solution with 1 factor
fit3 <- fa(mydata5,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.124, RMSR=.15, TLI=.365

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata5))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.161, RMSR=.18, TLI=.332

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata5, algorithm="louvain")
EGALV
# suggests 4 communities and response bias

# Give solution with 3 factors
fit4 <- fa(mydata5,3)
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.05, RMSR=.04, TLI=.897
#      MR2  MR1   MR3
#MR2  1.00 0.05 -0.06
#MR1  0.05 1.00  0.39
#MR3 -0.06 0.39  1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata5))
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.084, RMSR=.05, TLI=.817
#      MR2  MR1   MR3
#MR2  1.00 0.07 -0.03
#MR1  0.07 1.00  0.41
#MR3 -0.03 0.41  1.00

# Single factor model lavaan
UNImodel= '
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata5,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.559       0.366
#Tucker-Lewis Index (TLI)                       0.507       0.292
#Robust Comparative Fit Index (CFI)                         0.322
#Robust Tucker-Lewis Index (TLI)                            0.242
#RMSEA                                          0.198       0.177
#Robust RMSEA                                               0.171
#SRMR                                           0.196       0.196

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .187

CFA_model1 <- cfa(UNImodel, data = mydata5,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.385       0.346
#Tucker-Lewis Index (TLI)                       0.313       0.269
#Robust Comparative Fit Index (CFI)                         0.383
#Robust Tucker-Lewis Index (TLI)                            0.310
#RMSEA                                          0.131       0.129
#Robust RMSEA                                               0.130
#SRMR                                           0.158       0.158

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .039

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata5,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.875       0.803
#Tucker-Lewis Index (TLI)                       0.859       0.779
#Robust Comparative Fit Index (CFI)                         0.724
#Robust Tucker-Lewis Index (TLI)                            0.689
#RMSEA                                          0.106       0.099
#Robust RMSEA                                               0.110
#SRMR                                           0.113       0.113

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas              -0.001    0.045   -0.033    0.974   -0.001   -0.001

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .401

CFA_model2 <- cfa(EGAmodel,mydata5,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.763       0.760
#Tucker-Lewis Index (TLI)                       0.734       0.730
#Robust Comparative Fit Index (CFI)                         0.766
#Robust Tucker-Lewis Index (TLI)                            0.736
#RMSEA                                          0.081       0.078
#Robust RMSEA                                               0.081
#SRMR                                           0.090       0.090

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas               0.014    0.058    0.242    0.808    0.014    0.014

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .313


CFA_model3 <- cfa(EGAmodel,mydata5,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.764       0.761
#Tucker-Lewis Index (TLI)                       0.736       0.732
#Robust Comparative Fit Index (CFI)                         0.766
#Robust Tucker-Lewis Index (TLI)                            0.738
#RMSEA                                          0.081       0.078
#Robust RMSEA                                               0.080
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .313

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata5,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.862       0.852
#Tucker-Lewis Index (TLI)                       0.826       0.813
#Robust Comparative Fit Index (CFI)                         0.862
#Robust Tucker-Lewis Index (TLI)                            0.826
#RMSEA                                          0.066       0.065
#Robust RMSEA                                               0.066
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .327

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4049509      0.4789474      0.7945306      0.4240129 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#bis     0.9550776 0.3818335 0.04492241 0.8039946 0.7948256 0.8123607 0.9020111
#bas     0.3552368 0.2132156 0.64476316 0.8121499 0.1352461 0.7012272 0.8412875
#general 0.4049509 0.4049509 0.40495090 0.7945306 0.4240129 0.7796210 0.8720290








################################################################
### Analysis BISBAS Questionnaire
### Binter et al https://osf.io/98f4x

min(mydata6)
max(mydata6) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata6)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5factors and 3 components
# Eigenvalue 1 = 3.34; eigenvalue 2 = 2.29

rho <- polychoric(mydata6)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata6))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 4.06; eigenvalue 2 = 2.76

# Give solution with 1 factor
fit3 <- fa(mydata6,1)
fit3
diagram(fit3)
# %variance explained=.17, RMSEA=.122, RMSR=.14, TLI=.397

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata6))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.153, RMSR=.17, TLI=.361

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata6, algorithm="louvain")
EGALV
# suggests 3 communities and response bias

# Give solution with 3 factors
fit4 <- fa(mydata6,3)
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.059, RMSR=.04, TLI=.86
#     MR2  MR1  MR3
#MR2 1.00 0.03 0.02
#MR1 0.03 1.00 0.40
#MR3 0.02 0.40 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata6))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.086, RMSR=.05, TLI=.796
#     MR2  MR1  MR3
#MR2 1.00 0.05 0.01
#MR1 0.05 1.00 0.39
#MR3 0.01 0.39 1.00

# Single factor model lavaan
UNImodel= '
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata6,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.613       0.412
#Tucker-Lewis Index (TLI)                       0.568       0.343
#Robust Comparative Fit Index (CFI)                         0.400
#Robust Tucker-Lewis Index (TLI)                            0.330
#RMSEA                                          0.191       0.171
#Robust RMSEA                                               0.157
#SRMR                                           0.167       0.167

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .194


CFA_model1 <- cfa(UNImodel, data = mydata6,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.347       0.000
#Tucker-Lewis Index (TLI)                       0.271      -0.194
#Robust Comparative Fit Index (CFI)                         0.332
#Robust Tucker-Lewis Index (TLI)                            0.253
#RMSEA                                          0.135       0.162
#Robust RMSEA                                               0.136
#SRMR                                           0.161       0.161

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .078

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata6,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.894       0.814
#Tucker-Lewis Index (TLI)                       0.880       0.791
#Robust Comparative Fit Index (CFI)                         0.726
#Robust Tucker-Lewis Index (TLI)                            0.692
#RMSEA                                          0.101       0.096
#Robust RMSEA                                               0.106
#SRMR                                           0.098       0.098

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#       bas               0.047    0.040    1.189    0.235    0.047    0.047

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .358

CFA_model2 <- cfa(EGAmodel,mydata6,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.765       0.763
#Tucker-Lewis Index (TLI)                       0.736       0.733
#Robust Comparative Fit Index (CFI)                         0.768
#Robust Tucker-Lewis Index (TLI)                            0.739
#RMSEA                                          0.081       0.077
#Robust RMSEA                                               0.080
#SRMR                                           0.080       0.080

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#       bas               0.059    0.053    1.105    0.269    0.059    0.059

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .267


CFA_model3 <- cfa(EGAmodel,mydata6,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.765       0.763
#Tucker-Lewis Index (TLI)                       0.737       0.735
#Robust Comparative Fit Index (CFI)                         0.768
#Robust Tucker-Lewis Index (TLI)                            0.740
#RMSEA                                          0.081       0.076
#Robust RMSEA                                               0.080
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .267

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 bis =~ B2+B8+B13+B16+B19+B22+B24
 bas =~ B3+B9+B12+B21+B5+B10+B15+B20+B4+B7+B14+B18+B23
 general=~ B2+B3+B4+B5+B7+B8+B9+B10+B12+B13+B14+B15+B16+B18+B19+B20+B21+B22+B23+B24
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata6,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.857       0.849
#Tucker-Lewis Index (TLI)                       0.819       0.808
#Robust Comparative Fit Index (CFI)                         0.858
#Robust Tucker-Lewis Index (TLI)                            0.821
#RMSEA                                          0.067       0.065
#Robust RMSEA                                               0.067
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .367

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  bis =~                                                                
#    B2                0.374    0.031   12.036    0.000    0.374    0.434
#    B8                0.433    0.028   15.387    0.000    0.433    0.552
#    B13               0.507    0.026   19.706    0.000    0.507    0.655
#    B16               0.499    0.025   19.614    0.000    0.499    0.673
#    B19               0.451    0.025   18.164    0.000    0.451    0.656
#    B22               0.424    0.031   13.841    0.000    0.424    0.506
#    B24               0.510    0.027   19.083    0.000    0.510    0.639
#  bas =~                                                                
#    B3                0.347    0.056    6.201    0.000    0.347    0.540
#    B9                0.370    0.067    5.486    0.000    0.370    0.533
#    B12               0.377    0.098    3.835    0.000    0.377    0.528
#    B21               0.032    0.116    0.272    0.785    0.032    0.038
#    B5                0.045    0.102    0.443    0.657    0.045    0.066
#    B10              -0.165    0.128   -1.295    0.195   -0.165   -0.210
#    B15              -0.156    0.124   -1.259    0.208   -0.156   -0.190
#    B20               0.083    0.144    0.575    0.566    0.083    0.105
#    B4                0.192    0.044    4.401    0.000    0.192    0.385
#    B7                0.206    0.096    2.143    0.032    0.206    0.316
#    B14               0.202    0.113    1.787    0.074    0.202    0.306
#    B18               0.208    0.093    2.223    0.026    0.208    0.287
#    B23               0.124    0.086    1.447    0.148    0.124    0.172
#  general =~                                                            
#    B2                0.111    0.042    2.619    0.009    0.111    0.128
#    B3               -0.200    0.094   -2.118    0.034   -0.200   -0.311
#    B4               -0.094    0.059   -1.583    0.113   -0.094   -0.188
#    B5               -0.335    0.035   -9.586    0.000   -0.335   -0.489
#    B7               -0.254    0.071   -3.593    0.000   -0.254   -0.390
#    B8               -0.069    0.041   -1.707    0.088   -0.069   -0.088
#    B9               -0.254    0.100   -2.551    0.011   -0.254   -0.366
#    B10              -0.463    0.054   -8.631    0.000   -0.463   -0.588
#    B12              -0.346    0.105   -3.302    0.001   -0.346   -0.484
#    B13              -0.026    0.040   -0.648    0.517   -0.026   -0.034
#    B14              -0.376    0.065   -5.747    0.000   -0.376   -0.570
#    B15              -0.471    0.059   -7.947    0.000   -0.471   -0.575
#    B16              -0.063    0.037   -1.707    0.088   -0.063   -0.085
#    B18              -0.251    0.071   -3.527    0.000   -0.251   -0.347
#    B19              -0.037    0.040   -0.913    0.361   -0.037   -0.054
#    B20              -0.470    0.043  -11.016    0.000   -0.470   -0.597
#    B21              -0.460    0.036  -12.769    0.000   -0.460   -0.560
#    B22               0.195    0.042    4.703    0.000    0.195    0.233
#    B23              -0.207    0.051   -4.088    0.000   -0.207   -0.287
#    B24              -0.037    0.037   -1.003    0.316   -0.037   -0.047

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.4244887      0.4789474      0.8129870      0.4576771 

#$FactorLevelIndices
#        ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#bis     0.9640639 0.3666898 0.03593615 0.7922513 0.7921140 0.8014569 0.8964207
#bas     0.3370037 0.2088215 0.66299632 0.8241350 0.1649642 0.6435264 0.8016573
#general 0.4244887 0.4244887 0.42448871 0.8129870 0.4576771 0.7973353 0.8856143






################################################################
### Boutilier https://psyarxiv.com/cmtz3/download?format=pdf
### https://osf.io/28b9q
###
### Not included because already included Hope

library(haven)
Boutilier_RRB_Hope_Optimism_Data <- read_sav("Boutilier_RRB_Hope_Optimism_Data.sav")
colnames(Boutilier_RRB_Hope_Optimism_Data)
mydata <- as.data.frame(Boutilier_RRB_Hope_Optimism_Data[,c(30,32,33,35,37,39,41,43,45,47,49,51,53,55,
                                                            57,59,61,63,65,67,69,71,72,74)])
# take out the control questions
mydata <- mydata[-c(1,6,11,17)]
colnames(mydata)
mydata <- na.omit(mydata)
summary(mydata)
min(mydata)
max(mydata)


library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 components
# Eigenvalue 1 = 4.79; eigenvalue 2 = 2.84

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 5.57; eigenvalue 2 = 3.12

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.24, RMSEA=.147, RMSR=.17, TLI=.441

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.26, RMSEA=.17, RMSR=.19, TLI=.401

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities and response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.091, RMSR=.07, TLI=.784
#     MR1  MR2
#MR1 1.00 0.02
#MR2 0.02 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.114, RMSR=.07, TLI=.731
#     MR1  MR2
#MR1 1.00 0.04
#MR2 0.04 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ BISBAS2+BISBAS3+BISBAS4+BISBAS5+BISBAS7+BISBAS8+BISBAS9+BISBAS10+
           BISBAS12+BISBAS13+BISBAS14+BISBAS15+BISBAS16+BISBAS18+BISBAS19+
           BISBAS20+BISBAS21+BISBAS22+BISBAS23+BISBAS24
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.503       0.477
#Tucker-Lewis Index (TLI)                       0.456       0.428
#Robust Comparative Fit Index (CFI)                         0.503
#Robust Tucker-Lewis Index (TLI)                            0.456
#RMSEA                                          0.125       0.119
#Robust RMSEA                                               0.123
#SRMR                                           0.140       0.140

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .208

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 bis =~ BISBAS2+BISBAS8+BISBAS13+BISBAS16+BISBAS19+BISBAS22+BISBAS24
 bas =~ BISBAS3+BISBAS4+BISBAS5+BISBAS7+BISBAS9+BISBAS10+BISBAS12+BISBAS14+BISBAS15+
        BISBAS18+BISBAS20+BISBAS21+BISBAS23
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.788       0.779
#Tucker-Lewis Index (TLI)                       0.761       0.751
#Robust Comparative Fit Index (CFI)                         0.792
#Robust Tucker-Lewis Index (TLI)                            0.766
#RMSEA                                          0.099       0.093
#Robust RMSEA                                               0.097
#SRMR                                           0.094       0.094

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#bis ~~                                                                
#        bas                -0.032    0.127   -0.251    0.802   -0.032   -0.032


fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .388

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.788       0.782
#Tucker-Lewis Index (TLI)                       0.763       0.757
#Robust Comparative Fit Index (CFI)                         0.793
#Robust Tucker-Lewis Index (TLI)                            0.769
#RMSEA                                          0.098       0.092
#Robust RMSEA                                               0.096
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .389

# Bifactor model
BIFmodel= '
 bis =~ BISBAS2+BISBAS8+BISBAS13+BISBAS16+BISBAS19+BISBAS22+BISBAS24
 bas =~ BISBAS3+BISBAS4+BISBAS5+BISBAS7+BISBAS9+BISBAS10+BISBAS12+BISBAS14+BISBAS15+
        BISBAS18+BISBAS20+BISBAS21+BISBAS23
 general=~ BISBAS2+BISBAS3+BISBAS4+BISBAS5+BISBAS7+BISBAS8+BISBAS9+BISBAS10+
           BISBAS12+BISBAS13+BISBAS14+BISBAS15+BISBAS16+BISBAS18+BISBAS19+
           BISBAS20+BISBAS21+BISBAS22+BISBAS23+BISBAS24
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .447

# loadings based on Std.all because Std.lv is larger than 1
#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  bis =~                                                                
#    BISBAS2           1.161    0.106   10.990    0.000    1.161    0.711
#    BISBAS8          -1.165    0.106  -10.973    0.000   -1.165   -0.709
#    BISBAS13         -1.220    0.108  -11.307    0.000   -1.220   -0.737
#    BISBAS16         -1.188    0.095  -12.507    0.000   -1.188   -0.700
#    BISBAS19         -0.750    0.101   -7.390    0.000   -0.750   -0.575
#    BISBAS22          1.007    0.131    7.690    0.000    1.007    0.573
#    BISBAS24         -0.973    0.103   -9.438    0.000   -0.973   -0.652
#  bas =~                                                                
#    BISBAS3           1.159    0.130    8.901    0.000    1.159    0.732
#    BISBAS4           0.340    0.144    2.361    0.018    0.340    0.292
#    BISBAS5           0.321    0.220    1.462    0.144    0.321    0.226
#    BISBAS7           0.289    0.181    1.599    0.110    0.289    0.249
#    BISBAS9           1.366    0.136   10.076    0.000    1.366    0.860
#    BISBAS10          0.218    0.293    0.744    0.457    0.218    0.139
#    BISBAS12          1.019    0.217    4.707    0.000    1.019    0.647
#    BISBAS14          0.538    0.215    2.499    0.012    0.538    0.384
#    BISBAS15          0.290    0.284    1.020    0.308    0.290    0.170
#    BISBAS18          0.433    0.197    2.195    0.028    0.433    0.308
#    BISBAS20          0.181    0.313    0.579    0.563    0.181    0.109
#    BISBAS21          0.734    0.203    3.613    0.000    0.734    0.475
#    BISBAS23          0.282    0.216    1.307    0.191    0.282    0.174
#  general =~                                                            
#    BISBAS2           0.281    0.169    1.661    0.097    0.281    0.172
#    BISBAS3           0.319    0.314    1.016    0.310    0.319    0.202
#    BISBAS4           0.442    0.204    2.172    0.030    0.442    0.381
#    BISBAS5           0.801    0.161    4.962    0.000    0.801    0.562
#    BISBAS7           0.692    0.166    4.162    0.000    0.692    0.597
#    BISBAS8           0.236    0.149    1.583    0.113    0.236    0.143
#    BISBAS9           0.428    0.342    1.250    0.211    0.428    0.269
#    BISBAS10          0.818    0.189    4.316    0.000    0.818    0.519
#    BISBAS12          0.698    0.290    2.405    0.016    0.698    0.443
#    BISBAS13          0.429    0.164    2.616    0.009    0.429    0.259
#    BISBAS14          0.859    0.175    4.901    0.000    0.859    0.612
#    BISBAS15          0.868    0.177    4.898    0.000    0.868    0.509
#    BISBAS16          0.262    0.179    1.466    0.143    0.262    0.155
#    BISBAS18          0.594    0.192    3.088    0.002    0.594    0.423
#    BISBAS19          0.283    0.215    1.317    0.188    0.283    0.217
#    BISBAS20          1.133    0.149    7.627    0.000    1.133    0.681
#    BISBAS21          0.717    0.218    3.296    0.001    0.717    0.464
#    BISBAS22          0.363    0.213    1.707    0.088    0.363    0.207
#    BISBAS23          0.717    0.199    3.610    0.000    0.717    0.442
#    BISBAS24          0.253    0.223    1.136    0.256    0.253    0.170

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.3749281      0.4789474      0.8812578      0.5912325 

#$FactorLevelIndices
#        ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#bis     0.9231826 0.3502621 0.07681736 0.6283522 0.4484878 0.8550563 0.9280437
#bas     0.4428182 0.2748098 0.55718179 0.8893056 0.3366290 0.8483497 0.9085899
#general 0.3749281 0.3749281 0.37492807 0.8812578 0.5912325 0.8258105 0.8869893



