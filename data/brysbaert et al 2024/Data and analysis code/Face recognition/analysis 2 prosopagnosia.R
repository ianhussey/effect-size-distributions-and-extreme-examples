################################################################
### Prosopagnosia
### CMFQ
### Grabman & Dodson (2024) https://osf.io/9kjxp

mydata <- read.csv("Grabman & Dodson 01_CFMT and CFMQ e2 wide.csv")
colnames(mydata)
mydata <- as.data.frame(mydata[,13:29])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .88, omega T = .90

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 components
# Eigenvalue 1 = 5.44; eigenvalue 2 = .78

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 6.73; eigenvalue 2 = .94

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.113, RMSR=.08, TLI=.708

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.153, RMSR=.09, TLI=.646

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias

# Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.1, RMSR=.06, TLI=.769
#     MR1  MR2
#MR1 1.00 0.54
#MR2 0.54 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.139, RMSR=.07, TLI=.706
#     MR1  MR2
#MR1 1.00 0.46
#MR2 0.46 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  CFMQ1+CFMQ2+CFMQ3+CFMQ4+CFMQ5+CFMQ6+CFMQ7+CFMQ8+CFMQ9+CFMQ10+CFMQ11+CFMQ12+CFMQ13+
                CFMQ14+CFMQ15+CFMQ16+CFMQ17
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.955       0.877
#Tucker-Lewis Index (TLI)                       0.948       0.860
#Robust Comparative Fit Index (CFI)                         0.684
#Robust Tucker-Lewis Index (TLI)                            0.639
#RMSEA                                          0.119       0.131
#Robust RMSEA                                               0.156
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .395

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.745       0.742
#Tucker-Lewis Index (TLI)                       0.709       0.705
#Robust Comparative Fit Index (CFI)                         0.750
#Robust Tucker-Lewis Index (TLI)                            0.715
#RMSEA                                          0.113       0.100
#Robust RMSEA                                               0.112
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .333

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors based on EFA Pearson
EGAmodel= '
 factor1 =~  CFMQ15+CFMQ13+CFMQ3+CFMQ9+CFMQ14+CFMQ5+CFMQ17+CFMQ2
 factor2 =~  CFMQ8+CFMQ12+CFMQ1+CFMQ4+CFMQ11+CFMQ10+CFMQ6+CFMQ7+CFMQ16
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.961       0.891
#Tucker-Lewis Index (TLI)                       0.954       0.874
#Robust Comparative Fit Index (CFI)                         0.723
#Robust Tucker-Lewis Index (TLI)                            0.681
#RMSEA                                          0.112       0.125
#Robust RMSEA                                               0.146
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .414

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#       factor2           0.840    0.017   48.883    0.000    0.840    0.840

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.769       0.766
#Tucker-Lewis Index (TLI)                       0.734       0.731
#Robust Comparative Fit Index (CFI)                         0.775
#Robust Tucker-Lewis Index (TLI)                            0.740
#RMSEA                                          0.108       0.095
#Robust RMSEA                                               0.107
#SRMR                                           0.073       0.073

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .339

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.836    0.033   25.508    0.000    0.836    0.836

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.671       0.664
#Tucker-Lewis Index (TLI)                       0.624       0.615
#Robust Comparative Fit Index (CFI)                         0.676
#Robust Tucker-Lewis Index (TLI)                            0.630
#RMSEA                                          0.129       0.114
#Robust RMSEA                                               0.127
#SRMR                                           0.210       0.210

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .314

BIFmodel= '
 factor1 =~  CFMQ15+CFMQ13+CFMQ3+CFMQ9+CFMQ14+CFMQ5+CFMQ17+CFMQ2
 factor2 =~  CFMQ8+CFMQ12+CFMQ1+CFMQ4+CFMQ11+CFMQ10+CFMQ6+CFMQ7+CFMQ16
 general =~  CFMQ1+CFMQ2+CFMQ3+CFMQ4+CFMQ5+CFMQ6+CFMQ7+CFMQ8+CFMQ9+CFMQ10+CFMQ11+CFMQ12+CFMQ13+
                CFMQ14+CFMQ15+CFMQ16+CFMQ17
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.878
#Tucker-Lewis Index (TLI)                       0.847       0.838
#Robust Comparative Fit Index (CFI)                         0.889
#Robust Tucker-Lewis Index (TLI)                            0.852
#RMSEA                                          0.082       0.074
#Robust RMSEA                                               0.080
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .354

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7044833      0.5294118      0.9014899      0.8047772 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3205441 0.1374271 0.6794559 0.7997699 0.1902302 0.6150489 0.7902383
#factor2 0.2767340 0.1580896 0.7232660 0.8628851 0.1389703 0.6811543 0.8760557
#general 0.7044833 0.7044833 0.7044833 0.9014899 0.8047772 0.8884208 0.9304792





################################################################
### Prosopagnosia
### Data: PI20
### Norkaer et al, 2023, https://osf.io/rdfmw, https://www.academia.edu/98757012/The_Danish_Version_of_the_20_Item_Prosopagnosia_Index_PI20_Translation_Validation_and_a_Link_to_Face_Perception

library(readxl)
Danish_PI20_validation_study <- read_excel("Danish PI20 validation study.xlsx")
colnames(Danish_PI20_validation_study)
mydata <- as.data.frame(Danish_PI20_validation_study[,15:34])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
psych::alpha(mydata)
omega(mydata) # alpha = .9, omega = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 component
# Eigenvalue 1 = 6.98; eigenvalue 2 = 1.03

rho <- polychoric(mydata)$rho
# error message
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Eigenvalue 1 = 8.64; eigenvalue 2 = 1.21

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.35, RMSEA=.08, RMSR=.08, TLI=.831

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
# warning
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.372, RMSR=.1, TLI=.186

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.058, RMSR=.06, TLI=.909
#     MR1  MR2
#MR1 1.00 0.51
#MR2 0.51 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
# warning
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.383, RMSR=.08, TLI=.129
#     MR1  MR2
#MR1 1.00 0.56
#MR2 0.56 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.953
#Tucker-Lewis Index (TLI)                       0.990       0.947
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.055       0.076
#Robust RMSEA                                                  NA
#SRMR                                           0.099       0.099

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .469

CFA_model1 <- cfa(UNImodel, data = mydata, estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.837       0.838
#Tucker-Lewis Index (TLI)                       0.818       0.819
#Robust Comparative Fit Index (CFI)                         0.841
#Robust Tucker-Lewis Index (TLI)                            0.823
#RMSEA                                          0.087       0.084
#Robust RMSEA                                               0.086
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .336

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 factor1 =~  PI20_2+PI20_12+PI20_1+PI20_20+PI20_4+PI20_18+PI20_15+PI20_8+PI20_16+PI20_9+PI20_14+PI20_7+
                PI20_17+PI20_19+PI20_3
 factor2 =~  PI20_5+PI20_11+PI20_10+PI20_6
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.999       0.972
#Tucker-Lewis Index (TLI)                       0.999       0.969
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.020       0.061
#Robust RMSEA                                                  NA
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .522

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.741    0.069   10.738    0.000    0.741    0.741


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.879       0.889
#Tucker-Lewis Index (TLI)                       0.863       0.874
#Robust Comparative Fit Index (CFI)                         0.890
#Robust Tucker-Lewis Index (TLI)                            0.875
#RMSEA                                          0.080       0.072
#Robust RMSEA                                               0.075
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .395

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.700    0.090    7.770    0.000    0.700    0.700

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.830       0.841
#Tucker-Lewis Index (TLI)                       0.808       0.821
#Robust Comparative Fit Index (CFI)                         0.841
#Robust Tucker-Lewis Index (TLI)                            0.822
#RMSEA                                          0.094       0.086
#Robust RMSEA                                               0.090
#SRMR                                           0.172       0.172

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .439


# Bifactor model with 2 factors 
BIFmodel= '
 factor1 =~  PI20_2+PI20_12+PI20_1+PI20_20+PI20_4+PI20_18+PI20_15+PI20_8+PI20_16+PI20_9+PI20_14+PI20_7+
                PI20_17+PI20_19+PI20_3
 factor2 =~  PI20_5+PI20_11+PI20_10+PI20_6
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.922       0.901
#Tucker-Lewis Index (TLI)                       0.902       0.875
#Robust Comparative Fit Index (CFI)                         0.912
#Robust Tucker-Lewis Index (TLI)                            0.889
#RMSEA                                          0.064       0.070
#Robust RMSEA                                               0.068
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .446

#Latent Variables:
#                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 =~                                                            
#       PI20_2            0.350    0.316    1.107    0.268    0.350    0.440
#       PI20_12           0.251    0.344    0.732    0.464    0.251    0.295
#       PI20_1            0.415    0.270    1.534    0.125    0.415    0.462
#       PI20_20           0.421    0.173    2.431    0.015    0.421    0.443
#       PI20_4            0.366    0.174    2.101    0.036    0.366    0.379
#       PI20_18           0.064    0.186    0.342    0.733    0.064    0.107
#       PI20_15           0.020    0.271    0.075    0.940    0.020    0.027
#       PI20_8            0.263    0.247    1.065    0.287    0.263    0.243
#       PI20_16           0.008    0.378    0.022    0.982    0.008    0.008
#       PI20_9            0.501    0.153    3.267    0.001    0.501    0.463
#       PI20_14           0.252    0.180    1.396    0.163    0.252    0.278
#       PI20_7           -0.095    0.280   -0.341    0.733   -0.095   -0.100
#       PI20_17           0.438    0.143    3.059    0.002    0.438    0.383
#       PI20_19           0.156    0.140    1.116    0.265    0.156    0.159
#       PI20_3            0.259    0.197    1.312    0.190    0.259    0.248
#factor2 =~                                                            
#       PI20_5            0.269    0.122    2.203    0.028    0.269    0.423
#       PI20_11           0.383    0.099    3.853    0.000    0.383    0.568
#       PI20_10           0.305    0.110    2.782    0.005    0.305    0.397
#       PI20_6            0.305    0.102    2.986    0.003    0.305    0.353
#general =~                                                            
#       PI20_1            0.507    0.219    2.309    0.021    0.507    0.564
#       PI20_2            0.549    0.201    2.726    0.006    0.549    0.690
#       PI20_3            0.070    0.191    0.367    0.714    0.070    0.067
#       PI20_4            0.656    0.139    4.718    0.000    0.656    0.678
#       PI20_5            0.319    0.182    1.753    0.080    0.319    0.503
#       PI20_6            0.500    0.077    6.500    0.000    0.500    0.579
#       PI20_7            0.795    0.136    5.828    0.000    0.795    0.836
#       PI20_8            0.590    0.122    4.836    0.000    0.590    0.545
#       PI20_9            0.379    0.218    1.740    0.082    0.379    0.350
#       PI20_10           0.305    0.076    4.009    0.000    0.305    0.396
#       PI20_11           0.286    0.092    3.117    0.002    0.286    0.424
#       PI20_12           0.645    0.150    4.315    0.000    0.645    0.756
#       PI20_13           0.118    0.107    1.099    0.272    0.118    0.148
#       PI20_14           0.455    0.147    3.094    0.002    0.455    0.503
#       PI20_15           0.605    0.128    4.729    0.000    0.605    0.797
#       PI20_16           0.839    0.089    9.440    0.000    0.839    0.816
#       PI20_17           0.387    0.182    2.125    0.034    0.387    0.338
#       PI20_18           0.334    0.151    2.215    0.027    0.334    0.563
#       PI20_19           0.316    0.103    3.076    0.002    0.316    0.323
#       PI20_20           0.524    0.196    2.672    0.008    0.524    0.552

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#    0.7387271      0.4157895      0.9167435      0.7881892 
#$FactorLevelIndices
#ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2123287 0.16915904 0.7876713 0.9117408 0.1579168 0.6328263 0.8057520
#factor2 0.4589035 0.09211384 0.5410965 0.7437565 0.3392260 0.5059769 0.7446578
#general 0.7387271 0.73872712 0.7387271 0.9167435 0.7881892 0.9282482 0.9537931







################################################################
### Prosopagnosia
### Data: PI20
### Lowes et al, 2024, https://www.sciencedirect.com/science/article/pii/S0010945224000054?via%3Dihub
### Sent by the author

### Not included because dataset too small

library(readxl)
PI20_questionnaire_Jlowes <- read_excel("PI20_questionnaire_Jlowes.xlsx")
colnames(PI20_questionnaire_Jlowes)
mydata <- as.data.frame(PI20_questionnaire_Jlowes[,4:23])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
psych::alpha(mydata)
omega(mydata) # alpha = .97, omega = .98

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 1 component
# Eigenvalue 1 = 12.51; eigenvalue 2 = .6 (estimated from figure)

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Eigenvalue 1 = 13.98; eigenvalue 2 = .6

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.63, RMSEA=.063, RMSR=.05, TLI=.958

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.70, RMSEA=.355, RMSR=.05, TLI=.463

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 1 community

#Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.66, RMSEA=.044, RMSR=.04, TLI=.979
#     MR1  MR2
#MR1 1.00 0.55
#MR2 0.55 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.74, RMSEA=.367, RMSR=.04, TLI=.42
#
print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
CFA_model1 <- cfa(UNImodel, data = mydata, estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.960
#Tucker-Lewis Index (TLI)                       0.951       0.955
#Robust Comparative Fit Index (CFI)                         0.963
#Robust Tucker-Lewis Index (TLI)                            0.959
#RMSEA                                          0.073       0.063
#Robust RMSEA                                               0.066
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .395

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
#EGAmodel= '
# factor1 =~  PI20_2+PI20_12+PI20_1+PI20_20+PI20_4+PI20_18+PI20_15+PI20_8+PI20_16+PI20_9+PI20_14+PI20_7+
#                PI20_17+PI20_19+PI20_3
# factor2 =~  PI20_5+PI20_11+PI20_10+PI20_6
#'

#library(lavaan)
#CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
#summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)

#fitRsquares <- lavInspect(CFA_model2, what='rsquare')
#summary(fitRsquares)[-4]

#library(semPlot)
#semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
#         label.cex=2)


# Analysis with orthogonal factors
#CFA_model3 <- cfa(EGAmodel,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
#summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)

#fitRsquares <- lavInspect(CFA_model3, what='rsquare')
#summary(fitRsquares)[-4]


# Bifactor model with 2 factors 
#BIFmodel= '
# factor1 =~  PI20_2+PI20_12+PI20_1+PI20_20+PI20_4+PI20_18+PI20_15+PI20_8+PI20_16+PI20_9+PI20_14+PI20_7+
#                PI20_17+PI20_19+PI20_3
# factor2 =~  PI20_5+PI20_11+PI20_10+PI20_6
# general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
#                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
#'

#options(max.print=999999)
#CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
#summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#fitRsquares <- lavInspect(CFA_model4, what='rsquare')
#summary(fitRsquares)[-4]

#library(semPlot)
#semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
#         label.cex=2,bifactor="general")


#library(BifactorIndicesCalculator)
#bifactorIndices(CFA_model4)





################################################################
### Prosopagnosia
### Data: PI20
### Kramer, 2023, https://peerj.com/articles/14821/ (sent by the author)

library(readxl)
Kramer_2023_PI20 <- read_excel("Kramer_2023_PI20.xlsx")
colnames(Kramer_2023_PI20)
mydata <- as.data.frame(Kramer_2023_PI20[,1:20])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
psych::alpha(mydata)
omega(mydata) # alpha = .93, omega = .95
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 components
# Eigenvalue 1 = 9.68; eigenvalue 2 = 1.68

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 10.37; eigenvalue 2 = 1.97

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.48, RMSEA=.10, RMSR=.10, TLI=.864

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.13, RMSR=.12, TLI=.814

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.59, RMSEA=.06, RMSR=.03, TLI=.951
#     MR1  MR2
#MR1  1.00 -0.21
#MR2 -0.21  1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.64, RMSEA=.091, RMSR=.03, TLI=.909
#      MR1   MR2
#MR1  1.00 -0.22
#MR2 -0.22  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.983       0.915
#Tucker-Lewis Index (TLI)                       0.981       0.905
#Robust Comparative Fit Index (CFI)                         0.838
#Robust Tucker-Lewis Index (TLI)                            0.819
#RMSEA                                          0.152       0.176
#Robust RMSEA                                               0.131
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .676

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.878       0.887
#Tucker-Lewis Index (TLI)                       0.864       0.874
#Robust Comparative Fit Index (CFI)                         0.890
#Robust Tucker-Lewis Index (TLI)                            0.877
#RMSEA                                          0.102       0.082
#Robust RMSEA                                               0.096
#SRMR                                           0.095       0.095

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .628

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 factor1 =~  PI20_7+PI20_5+PI20_4+PI20_20+PI20_10+PI20_16+PI20_18+PI20_14+PI20_15+PI20_2+PI20_6+PI20_1+
                PI20_12+PI20_11
 factor2 =~  PI20_8+PI20_19+PI20_17+PI20_13+PI20_9+PI20_3
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.983
#Tucker-Lewis Index (TLI)                       0.996       0.981
#Robust Comparative Fit Index (CFI)                         0.914
#Robust Tucker-Lewis Index (TLI)                            0.903
#RMSEA                                          0.070       0.078
#Robust RMSEA                                               0.096
#SRMR                                           0.064       0.064

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2          -0.298    0.043   -6.855    0.000   -0.298   -0.298

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .677

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.952       0.962
#Tucker-Lewis Index (TLI)                       0.946       0.958
#Robust Comparative Fit Index (CFI)                         0.964
#Robust Tucker-Lewis Index (TLI)                            0.959
#RMSEA                                          0.064       0.047
#Robust RMSEA                                               0.055
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .627

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2          -0.245    0.071   -3.460    0.001   -0.245   -0.245

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.949       0.959
#Tucker-Lewis Index (TLI)                       0.943       0.955
#Robust Comparative Fit Index (CFI)                         0.961
#Robust Tucker-Lewis Index (TLI)                            0.956
#RMSEA                                          0.066       0.049
#Robust RMSEA                                               0.057
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .626


# Bifactor model with 2 factors 
BIFmodel= '
 factor1 =~  PI20_7+PI20_5+PI20_4+PI20_20+PI20_10+PI20_16+PI20_18+PI20_14+PI20_15+PI20_2+PI20_6+PI20_1+
                PI20_12+PI20_11
 factor2 =~  PI20_8+PI20_19+PI20_17+PI20_13+PI20_9+PI20_3
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.971       0.974
#Tucker-Lewis Index (TLI)                       0.963       0.968
#Robust Comparative Fit Index (CFI)                         0.978
#Robust Tucker-Lewis Index (TLI)                            0.972
#RMSEA                                          0.053       0.041
#Robust RMSEA                                               0.046
#SRMR                                           0.043       0.043

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
#     0.4235063      0.4421053      0.9395847      0.3718267 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.5129817 0.4177424 0.4870183 0.9693615 0.5019170 0.9001591 0.8779922
#factor2 0.8550727 0.1587512 0.1449273 0.6529732 0.5606422 0.7525059 0.8674283
#general 0.4235063 0.4235063 0.4235063 0.9395847 0.3718267 0.8930739 0.8746272







################################################################
### Prosopagnosia
### Data: PI20
### Estudillo & Wong, 2021, https://peerj.com/articles/10629/ (sent by the author)

Estudillo_2021_PI20_PEERJ <- read.csv("Estudillo_2021_PI20_PEERJ.csv")
library(stringr)
Estudillo_2021_PI20_PEERJ$head2<-str_extract(Estudillo_2021_PI20_PEERJ$head, "[^.]+")
Estudillo_2021_PI20_PEERJ <- Estudillo_2021_PI20_PEERJ[,c(2,3,5)]
library(tidyr)
mydata <- spread(Estudillo_2021_PI20_PEERJ, head2, response)
colnames(mydata) <- c("id","PI20_1","PI20_2","PI20_3","PI20_4","PI20_5","PI20_6","PI20_7","PI20_8","PI20_9","PI20_10",
                      "PI20_11","PI20_12","PI20_13","PI20_14","PI20_15","PI20_16","PI20_17","PI20_18","PI20_19","PI20_20")
mydata <- mydata[,2:21]
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
psych::alpha(mydata)
omega(mydata) # alpha = .88, omega = .90
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 components
# Eigenvalue 1 = 6.30; eigenvalue 2 = .62

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 7.32; eigenvalue 2 = .76

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.071, RMSR=.06, TLI=.851

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.103, RMSR=.07, TLI=.776

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 3 communities

# Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.35, RMSEA=.059, RMSR=.05, TLI=.896
#     MR1  MR2
#MR1 1.00 0.67
#MR2 0.67 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.089, RMSR=.05, TLI=.830
#      MR1   MR2
#MR1 1.00 0.68
#MR2 0.68 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.990       0.958
#Tucker-Lewis Index (TLI)                       0.989       0.953
#Robust Comparative Fit Index (CFI)                         0.825
#Robust Tucker-Lewis Index (TLI)                            0.804
#RMSEA                                          0.050       0.068
#Robust RMSEA                                               0.098
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .308

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.865       0.873
#Tucker-Lewis Index (TLI)                       0.849       0.858
#Robust Comparative Fit Index (CFI)                         0.877
#Robust Tucker-Lewis Index (TLI)                            0.862
#RMSEA                                          0.073       0.065
#Robust RMSEA                                               0.069
#SRMR                                           0.058       0.058

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .258

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 factor1 =~  PI20_12+PI20_1+PI20_15+PI20_9+PI20_4+PI20_19+PI20_11+PI20_20+PI20_17
 factor2 =~  PI20_3+PI20_8+PI20_10+PI20_13+PI20_2+PI20_7+PI20_16+PI20_18+PI20_6+
                PI20_14+PI20_5
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.995       0.970
#Tucker-Lewis Index (TLI)                       0.994       0.967
#Robust Comparative Fit Index (CFI)                         0.867
#Robust Tucker-Lewis Index (TLI)                            0.850
#RMSEA                                          0.036       0.058
#Robust RMSEA                                               0.086
#SRMR                                           0.060       0.060

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2           0.865    0.022   39.080    0.000    0.865    0.865

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .322

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.898       0.908
#Tucker-Lewis Index (TLI)                       0.886       0.896
#Robust Comparative Fit Index (CFI)                         0.911
#Robust Tucker-Lewis Index (TLI)                            0.899
#RMSEA                                          0.064       0.055
#Robust RMSEA                                               0.059
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .268

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2           0.858    0.029   29.715    0.000    0.858    0.858

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.780       0.784
#Tucker-Lewis Index (TLI)                       0.754       0.758
#Robust Comparative Fit Index (CFI)                         0.790
#Robust Tucker-Lewis Index (TLI)                            0.765
#RMSEA                                          0.093       0.084
#Robust RMSEA                                               0.090
#SRMR                                           0.207       0.207

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .267


# Bifactor model with 2 factors 
BIFmodel= '
 factor1 =~  PI20_12+PI20_1+PI20_15+PI20_9+PI20_4+PI20_19+PI20_11+PI20_20+PI20_17
 factor2 =~  PI20_3+PI20_8+PI20_10+PI20_13+PI20_2+PI20_7+PI20_16+PI20_18+PI20_6+
                PI20_14+PI20_5
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.928       0.928
#Tucker-Lewis Index (TLI)                       0.909       0.909
#Robust Comparative Fit Index (CFI)                         0.934
#Robust Tucker-Lewis Index (TLI)                            0.916
#RMSEA                                          0.057       0.052
#Robust RMSEA                                               0.054
#SRMR                                           0.048       0.048

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .285

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7844652      0.5210526      0.7752699      0.7534942 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.2964580 0.15165303 0.7035420 0.3699382 0.07320276 0.5892335 0.8083001
#factor2 0.1307845 0.06388172 0.8692155 0.7715468 0.01902736 0.3447445 0.7406909
#general 0.7844652 0.78446525 0.7844652 0.7752699 0.75349420 0.9148531 0.9530336







################################################################
### Prosopagnosia
### Data: PI20
### Estudillo, 2021, https://www.emerald.com/insight/content/doi/10.1108/JCP-06-2020-0025/full/html (sent by the author)

Estudillo_2021_PI20_JCP <- read.csv("Estudillo_2021_PI20_JCP.csv")
Estudillo_2021_PI20_JCP$head2<-(as.numeric(rownames(Estudillo_2021_PI20_JCP))-1) %% 20
Estudillo_2021_PI20_JCP <- Estudillo_2021_PI20_JCP[,c(2,4,5)]
library(tidyr)
mydata <- spread(Estudillo_2021_PI20_JCP, head2, response)
colnames(mydata) <- c("id","PI20_1","PI20_2","PI20_3","PI20_4","PI20_5","PI20_6","PI20_7","PI20_8","PI20_9","PI20_10",
                      "PI20_11","PI20_12","PI20_13","PI20_14","PI20_15","PI20_16","PI20_17","PI20_18","PI20_19","PI20_20")
mydata <- mydata[,2:21]
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
psych::alpha(mydata)
omega(mydata) # alpha = .90, omega = .91
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 1 component
# Eigenvalue 1 = 6.54; eigenvalue 2 = .85 (estimated on the basis of the figure)

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 7.43; eigenvalue 2 = 1.00

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.33, RMSEA=.058, RMSR=.08, TLI=.885

fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.125, RMSR=.09, TLI=.676

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.044, RMSR=.07, TLI=.934
#     MR1  MR2
#MR1 1.00 0.55
#MR2 0.55 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.12, RMSR=.07, TLI=.696
#      MR1   MR2
#MR1 1.00 0.41
#MR2 0.41 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.995       0.950
#Tucker-Lewis Index (TLI)                       0.994       0.944
#Robust Comparative Fit Index (CFI)                         0.782
#Robust Tucker-Lewis Index (TLI)                            0.756
#RMSEA                                          0.039       0.073
#Robust RMSEA                                               0.120
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.869       0.877
#Tucker-Lewis Index (TLI)                       0.853       0.862
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.866
#RMSEA                                          0.073       0.067
#Robust RMSEA                                               0.069
#SRMR                                           0.077       0.077

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .237

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 factor1 =~  PI20_2+PI20_1+PI20_15+PI20_7+PI20_12+PI20_4+PI20_5+PI20_16+PI20_19+
                PI20_14+PI20_9+PI20_8+PI20_17
 factor2 =~  PI20_20+PI20_3+PI20_6+PI20_13+PI20_11+PI20_18+PI20_10
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    1.000       0.967
#Tucker-Lewis Index (TLI)                       1.003       0.963
#Robust Comparative Fit Index (CFI)                         0.826
#Robust Tucker-Lewis Index (TLI)                            0.804
#RMSEA                                          0.000       0.059
#Robust RMSEA                                               0.107
#SRMR                                           0.083       0.083

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .407

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2           0.756    0.052   14.507    0.000    0.756    0.756

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.909       0.917
#Tucker-Lewis Index (TLI)                       0.897       0.906
#Robust Comparative Fit Index (CFI)                         0.919
#Robust Tucker-Lewis Index (TLI)                            0.909
#RMSEA                                          0.061       0.055
#Robust RMSEA                                               0.057
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#    factor2           0.734    0.084    8.748    0.000    0.734    0.734

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.845       0.850
#Tucker-Lewis Index (TLI)                       0.827       0.832
#Robust Comparative Fit Index (CFI)                         0.854
#Robust Tucker-Lewis Index (TLI)                            0.837
#RMSEA                                          0.079       0.074
#Robust RMSEA                                               0.076
#SRMR                                           0.178       0.178

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .299


# Bifactor model with 2 factors 
BIFmodel= '
 factor1 =~  PI20_2+PI20_1+PI20_15+PI20_7+PI20_12+PI20_4+PI20_5+PI20_16+PI20_19+
                PI20_14+PI20_9+PI20_8+PI20_17
 factor2 =~  PI20_20+PI20_3+PI20_6+PI20_13+PI20_11+PI20_18+PI20_10
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.951       0.946
#Tucker-Lewis Index (TLI)                       0.938       0.932
#Robust Comparative Fit Index (CFI)                         0.950
#Robust Tucker-Lewis Index (TLI)                            0.936
#RMSEA                                          0.047       0.047
#Robust RMSEA                                               0.047
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .433

#Latent Variables:
#                       Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 =~                                                            
#    PI20_2            0.815    0.136    5.988    0.000    0.815    0.663
#    PI20_1            0.836    0.135    6.193    0.000    0.836    0.780
#    PI20_15           0.590    0.140    4.207    0.000    0.590    0.544
#    PI20_7            0.445    0.245    1.815    0.069    0.445    0.347
#    PI20_12           0.449    0.140    3.207    0.001    0.449    0.389
#    PI20_4            0.361    0.218    1.659    0.097    0.361    0.311
#    PI20_5            0.331    0.138    2.406    0.016    0.331    0.306
#    PI20_16           0.128    0.288    0.443    0.658    0.128    0.114
#    PI20_19          -0.243    0.196   -1.238    0.216   -0.243   -0.221
#    PI20_14           0.218    0.271    0.803    0.422    0.218    0.163
#    PI20_9           -0.359    0.137   -2.614    0.009   -0.359   -0.315
#    PI20_8           -0.194    0.136   -1.428    0.153   -0.194   -0.168
#    PI20_17          -0.277    0.177   -1.569    0.117   -0.277   -0.226
#factor2 =~                                                            
#    PI20_20           0.263    0.187    1.409    0.159    0.263    0.244
#    PI20_3           -0.629    0.262   -2.402    0.016   -0.629   -0.721
#    PI20_6            0.138    0.215    0.642    0.521    0.138    0.120
#    PI20_13          -0.443    0.266   -1.663    0.096   -0.443   -0.378
#    PI20_11           0.380    0.216    1.756    0.079    0.380    0.358
#    PI20_18           0.182    0.171    1.064    0.287    0.182    0.187
#    PI20_10           0.274    0.130    2.111    0.035    0.274    0.291
#general =~                                                            
#    PI20_1            0.453    0.183    2.472    0.013    0.453    0.422
#    PI20_2            0.629    0.174    3.618    0.000    0.629    0.512
#    PI20_3           -0.102    0.160   -0.639    0.523   -0.102   -0.117
#    PI20_4            0.713    0.172    4.144    0.000    0.713    0.614
#    PI20_5            0.662    0.130    5.083    0.000    0.662    0.611
#    PI20_6            0.813    0.193    4.224    0.000    0.813    0.711
#    PI20_7            0.841    0.150    5.607    0.000    0.841    0.656
#    PI20_8           -0.535    0.148   -3.621    0.000   -0.535   -0.462
#    PI20_9           -0.345    0.140   -2.459    0.014   -0.345   -0.302
#    PI20_10           0.366    0.127    2.881    0.004    0.366    0.388
#    PI20_11           0.422    0.161    2.619    0.009    0.422    0.398
#    PI20_12           0.759    0.149    5.108    0.000    0.759    0.658
#    PI20_13          -0.297    0.215   -1.379    0.168   -0.297   -0.254
#    PI20_14           0.663    0.180    3.691    0.000    0.663    0.497
#    PI20_15           0.693    0.150    4.622    0.000    0.693    0.639
#    PI20_16           0.746    0.193    3.873    0.000    0.746    0.664
#    PI20_17          -0.357    0.181   -1.971    0.049   -0.357   -0.291
#    PI20_18           0.415    0.167    2.484    0.013    0.415    0.428
#    PI20_19          -0.452    0.135   -3.336    0.001   -0.452   -0.411
#    PI20_20           0.641    0.163    3.924    0.000    0.641    0.595

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.6268695      0.4789474      0.7833458      0.6505494 
#$FactorLevelIndices
#          ECV_SS    ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#factor1 0.3586379 0.2530532 0.6413621 0.7504804 0.249779261 0.7813406 0.8690917
#factor2 0.4078647 0.1200773 0.5921353 0.5020218 0.001146625 0.6154497 0.7792451
#general 0.6268695 0.6268695 0.6268695 0.7833458 0.650549353 0.8889532 0.9130226






################################################################
### Prosopagnosia
### Data: PI20
### Dal Lago et al, 2023, https://www.tandfonline.com/doi/pdf/10.1080/10826084.2023.2247059 (data sent by authors)

library(readxl)
Dal_Lago_2023_PI20_raw_data <- read_excel("Dal_Lago_2023_PI20_raw_data.xlsx")
colnames(Dal_Lago_2023_PI20_raw_data)
mydata <- as.data.frame(Dal_Lago_2023_PI20_raw_data[,4:23])
mydata <- na.omit(mydata)
colnames(mydata) <- c("PI20_1","PI20_2","PI20_3","PI20_4","PI20_5","PI20_6","PI20_7","PI20_8","PI20_9","PI20_10",
                      "PI20_11","PI20_12","PI20_13","PI20_14","PI20_15","PI20_16","PI20_17","PI20_18","PI20_19","PI20_20")
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
psych::alpha(mydata)
omega(mydata) # alpha = .89, omega = .91
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 6.34; eigenvalue 2 = .85

rho <- polychoric(mydata)$rho
# does not run; so, alternative
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 8.05; eigenvalue 2 = 1.07

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.32, RMSEA=.086, RMSR=.07, TLI=.791

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.143, RMSR=.09, TLI=.668

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 4 communities

# Give solution with 2 factors (EGA-based)
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.069, RMSR=.06, TLI=.867
#     MR1  MR2
#MR1 1.00 0.58
#MR2 0.58 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.124, RMSR=.06, TLI=.747
#     MR1  MR2
#MR1 1.00 0.54
#MR2 0.54 1.00
print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.916
#Tucker-Lewis Index (TLI)                       0.977       0.906
#Robust Comparative Fit Index (CFI)                         0.744
#Robust Tucker-Lewis Index (TLI)                            0.714
#RMSEA                                          0.072       0.089
#Robust RMSEA                                               0.135
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .439

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.810       0.809
#Tucker-Lewis Index (TLI)                       0.788       0.786
#Robust Comparative Fit Index (CFI)                         0.819
#Robust Tucker-Lewis Index (TLI)                            0.797
#RMSEA                                          0.089       0.081
#Robust RMSEA                                               0.086
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .317

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 2 factors 
EGAmodel= '
 factor1 =~  PI20_2+PI20_1+PI20_9+PI20_4+PI20_6+PI20_8+PI20_16+PI20_14+PI20_19+PI20_3+PI20_20+PI20_17
 factor2 =~  PI20_11+PI20_10+PI20_18+PI20_5+PI20_7+PI20_15+PI20_12+PI20_13
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.986       0.937
#Tucker-Lewis Index (TLI)                       0.985       0.929
#Robust Comparative Fit Index (CFI)                         0.792
#Robust Tucker-Lewis Index (TLI)                            0.767
#RMSEA                                          0.059       0.077
#Robust RMSEA                                               0.122
#SRMR                                           0.078       0.078

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.840    0.027   31.689    0.000    0.840    0.840

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .470

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.852       0.852
#Tucker-Lewis Index (TLI)                       0.833       0.834
#Robust Comparative Fit Index (CFI)                         0.860
#Robust Tucker-Lewis Index (TLI)                            0.843
#RMSEA                                          0.079       0.071
#Robust RMSEA                                               0.075
#SRMR                                           0.065       0.065

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .346

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.830    0.040   20.493    0.000    0.830    0.830

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.751       0.745
#Tucker-Lewis Index (TLI)                       0.722       0.715
#Robust Comparative Fit Index (CFI)                         0.759
#Robust Tucker-Lewis Index (TLI)                            0.730
#RMSEA                                          0.101       0.093
#Robust RMSEA                                               0.099
#SRMR                                           0.207       0.207

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .335


# Bifactor model with 2 factors 
BIFmodel= '
 factor1 =~  PI20_2+PI20_1+PI20_9+PI20_4+PI20_6+PI20_8+PI20_16+PI20_14+PI20_19+PI20_3+PI20_20+PI20_17
 factor2 =~  PI20_11+PI20_10+PI20_18+PI20_5+PI20_7+PI20_15+PI20_12+PI20_13
 general =~  PI20_1+PI20_2+PI20_3+PI20_4+PI20_5+PI20_6+PI20_7+PI20_8+PI20_9+PI20_10+PI20_11+PI20_12+
                PI20_13+PI20_14+PI20_15+PI20_16+PI20_17+PI20_18+PI20_19+PI20_20
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.904       0.884
#Tucker-Lewis Index (TLI)                       0.879       0.854
#Robust Comparative Fit Index (CFI)                         0.902
#Robust Tucker-Lewis Index (TLI)                            0.876
#RMSEA                                          0.067       0.067
#Robust RMSEA                                               0.067
#SRMR                                           0.055       0.055

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .363

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.7709716      0.5052632      0.8104346      0.7293679 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.2300568 0.1223517 0.7699432 0.5803781 0.1065512 0.5764170 0.8329725
#factor2 0.2278602 0.1066767 0.7721398 0.8070112 0.1426633 0.5012175 0.7212661
#general 0.7709716 0.7709716 0.7709716 0.8104346 0.7293679 0.9143930 0.9425224






