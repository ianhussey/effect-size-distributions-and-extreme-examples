################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 1. Australia

library(haven)
mydata <- read_sav("australia.sav")
colnames(mydata)
mydata <- as.data.frame(mydata[,3:41])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .84, omega T = .88

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 components
# Eigenvalue 1 = 6.51; eigenvalue 2 = 3.91

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 components
# Eigenvalue 1 = 7.47; eigenvalue 2 = 4.40

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.17, RMSEA=.111, RMSR=.15, TLI=.297

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.19, RMSEA=.148, RMSR=.16, TLI=.216

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities if response bias


# Give solution with 5 factors (EGA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.059, RMSR=.05, TLI=.799
#     MR3   MR2   MR4  MR1   MR5
#MR3 1.00  0.10  0.12 0.29  0.09
#MR2 0.10  1.00  0.28 0.19 -0.09
#MR4 0.12  0.28  1.00 0.19 -0.11
#MR1 0.29  0.19  0.19 1.00  0.27
#MR5 0.09 -0.09 -0.11 0.27  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.108, RMSR=.06, TLI=.57
#     MR1   MR2   MR4  MR5   MR3
#MR1 1.00  0.10  0.14 0.33  0.09
#MR2 0.10  1.00  0.29 0.20 -0.08
#MR4 0.14  0.29  1.00 0.20 -0.10
#MR5 0.33  0.20  0.20 1.00  0.29
#MR3 0.09 -0.08 -0.10 0.29  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.666       0.506
#Tucker-Lewis Index (TLI)                       0.647       0.478
#Robust Comparative Fit Index (CFI)                         0.279
#Robust Tucker-Lewis Index (TLI)                            0.239
#RMSEA                                          0.197       0.122
#Robust RMSEA                                               0.158
#SRMR                                           0.166       0.166

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .246

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.326       0.293
#Tucker-Lewis Index (TLI)                       0.289       0.254
#Robust Comparative Fit Index (CFI)                         0.322
#Robust Tucker-Lewis Index (TLI)                            0.284
#RMSEA                                          0.120       0.114
#Robust RMSEA                                               0.118
#SRMR                                           0.142       0.142

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .141

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.941       0.875
#Tucker-Lewis Index (TLI)                       0.937       0.866
#Robust Comparative Fit Index (CFI)                         0.640
#Robust Tucker-Lewis Index (TLI)                            0.614
#RMSEA                                          0.083       0.062
#Robust RMSEA                                               0.112
#SRMR                                           0.096       0.096

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.364    0.073    4.988    0.000    0.364    0.364
#    nonjudging       -0.055    0.075   -0.724    0.469   -0.055   -0.055
#    nonreactivity     0.517    0.078    6.649    0.000    0.517    0.517
#    acting            0.014    0.079    0.182    0.856    0.014    0.014
#  describing ~~                                                         
#    nonjudging        0.286    0.066    4.342    0.000    0.286    0.286
#    nonreactivity     0.482    0.060    8.066    0.000    0.482    0.482
#    acting            0.207    0.064    3.215    0.001    0.207    0.207
#  nonjudging ~~                                                         
#    nonreactivity     0.293    0.077    3.826    0.000    0.293    0.293
#    acting            0.407    0.056    7.295    0.000    0.407    0.407
#  nonreactivity ~~                                                      
#    acting            0.315    0.070    4.522    0.000    0.315    0.315

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .479

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.758       0.766
#Tucker-Lewis Index (TLI)                       0.741       0.750
#Robust Comparative Fit Index (CFI)                         0.773
#Robust Tucker-Lewis Index (TLI)                            0.757
#RMSEA                                          0.072       0.066
#Robust RMSEA                                               0.069
#SRMR                                           0.087       0.087

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.282    0.124    2.285    0.022    0.282    0.282
#    nonjudging       -0.050    0.113   -0.437    0.662   -0.050   -0.050
#    nonreactivity     0.488    0.122    4.010    0.000    0.488    0.488
#    acting           -0.050    0.118   -0.423    0.673   -0.050   -0.050
#  describing ~~                                                         
#    nonjudging        0.263    0.094    2.784    0.005    0.263    0.263
#    nonreactivity     0.460    0.079    5.850    0.000    0.460    0.460
#    acting            0.198    0.087    2.273    0.023    0.198    0.198
#  nonjudging ~~                                                         
#    nonreactivity     0.291    0.119    2.447    0.014    0.291    0.291
#    acting            0.391    0.089    4.392    0.000    0.391    0.391
#  nonreactivity ~~                                                      
#    acting            0.271    0.108    2.501    0.012    0.271    0.271

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .396

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.721       0.731
#Tucker-Lewis Index (TLI)                       0.706       0.716
#Robust Comparative Fit Index (CFI)                         0.737
#Robust Tucker-Lewis Index (TLI)                            0.722
#RMSEA                                          0.077       0.071
#Robust RMSEA                                               0.074
#SRMR                                           0.140       0.140


fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .391


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.796       0.801
#Tucker-Lewis Index (TLI)                       0.772       0.777
#Robust Comparative Fit Index (CFI)                         0.808
#Robust Tucker-Lewis Index (TLI)                            0.785
#RMSEA                                          0.068       0.062
#Robust RMSEA                                               0.065
#SRMR                                           0.109       0.109

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .431

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing =~                                                          
#    EM_1              0.533    0.092    5.769    0.000    0.533    0.509
#    EM_6              0.581    0.113    5.133    0.000    0.581    0.511
#    EM_11             0.533    0.113    4.700    0.000    0.533    0.445
#    EM_15             0.680    0.072    9.403    0.000    0.680    0.770
#    EM_20             0.565    0.101    5.621    0.000    0.565    0.558
#    EM_26             0.514    0.108    4.775    0.000    0.514    0.531
#    EM_31             0.552    0.090    6.160    0.000    0.552    0.540
#    EM_36             0.445    0.107    4.173    0.000    0.445    0.408
#  describing =~                                                         
#    EM_2              0.644    0.066    9.732    0.000    0.644    0.728
#    EM_7              0.621    0.069    9.068    0.000    0.621    0.692
#    EM_12             0.575    0.065    8.917    0.000    0.575    0.672
#    EM_16             0.602    0.074    8.124    0.000    0.602    0.692
#    EM_22             0.521    0.075    6.929    0.000    0.521    0.619
#    EM_27             0.520    0.081    6.450    0.000    0.520    0.542
#    EM_32             0.660    0.081    8.171    0.000    0.660    0.620
#    EM_37             0.697    0.068   10.294    0.000    0.697    0.715
#  nonjudging =~                                                         
#    EM_3              0.542    0.074    7.349    0.000    0.542    0.536
#    EM_10             0.548    0.081    6.781    0.000    0.548    0.547
#    EM_14             0.521    0.077    6.746    0.000    0.521    0.553
#    EM_17             0.543    0.071    7.644    0.000    0.543    0.581
#    EM_25             0.578    0.069    8.370    0.000    0.578    0.625
#    EM_30             0.714    0.072    9.976    0.000    0.714    0.678
#    EM_35             0.605    0.078    7.784    0.000    0.605    0.594
#    EM_39             0.630    0.075    8.363    0.000    0.630    0.590
#  nonreactivity =~                                                      
#    EM_4              0.549    0.070    7.805    0.000    0.549    0.604
#    EM_9              0.519    0.089    5.808    0.000    0.519    0.545
#    EM_19             0.529    0.098    5.403    0.000    0.529    0.541
#    EM_21             0.464    0.085    5.486    0.000    0.464    0.545
#    EM_24             0.597    0.081    7.366    0.000    0.597    0.593
#    EM_29             0.615    0.069    8.866    0.000    0.615    0.687
#    EM_33             0.639    0.079    8.070    0.000    0.639    0.651
#  acting =~                                                             
#    EM_5              0.494    0.299    1.655    0.098    0.494    0.532
#    EM_8              0.377    0.295    1.275    0.202    0.377    0.396
#    EM_13             0.577    0.329    1.752    0.080    0.577    0.593
#    EM_18             0.297    0.241    1.231    0.218    0.297    0.335
#    EM_23            -0.163    0.290   -0.561    0.575   -0.163   -0.176
#    EM_28            -0.060    0.248   -0.244    0.807   -0.060   -0.069
#    EM_34            -0.341    0.279   -1.222    0.222   -0.341   -0.379
#    EM_38            -0.152    0.370   -0.411    0.681   -0.152   -0.187
#  general =~                                                            
#    EM_1             -0.030    0.178   -0.168    0.866   -0.030   -0.029
#    EM_2             -0.026    0.084   -0.314    0.754   -0.026   -0.030
#    EM_3              0.247    0.104    2.369    0.018    0.247    0.245
#    EM_4             -0.058    0.120   -0.489    0.625   -0.058   -0.064
#    EM_5              0.611    0.258    2.369    0.018    0.611    0.657
#    EM_6              0.126    0.216    0.584    0.559    0.126    0.111
#    EM_7              0.082    0.088    0.934    0.350    0.082    0.092
#    EM_8              0.603    0.215    2.810    0.005    0.603    0.635
#    EM_9              0.242    0.104    2.325    0.020    0.242    0.254
#    EM_10             0.376    0.107    3.506    0.000    0.376    0.375
#    EM_11            -0.130    0.112   -1.159    0.246   -0.130   -0.108
#    EM_12             0.155    0.080    1.922    0.055    0.155    0.181
#    EM_13             0.609    0.301    2.023    0.043    0.609    0.626
#    EM_14             0.333    0.120    2.781    0.005    0.333    0.353
#    EM_15            -0.072    0.088   -0.817    0.414   -0.072   -0.082
#    EM_16             0.221    0.093    2.377    0.017    0.221    0.254
#    EM_17             0.215    0.090    2.393    0.017    0.215    0.230
#    EM_18             0.511    0.170    3.010    0.003    0.511    0.576
#    EM_19             0.163    0.124    1.323    0.186    0.163    0.167
#    EM_20             0.045    0.104    0.435    0.663    0.045    0.045
#    EM_21             0.220    0.095    2.315    0.021    0.220    0.259
#    EM_22             0.142    0.091    1.562    0.118    0.142    0.168
#    EM_23             0.641    0.113    5.654    0.000    0.641    0.692
#    EM_24             0.201    0.091    2.201    0.028    0.201    0.200
#    EM_25             0.226    0.099    2.281    0.023    0.226    0.244
#    EM_26             0.066    0.114    0.581    0.561    0.066    0.069
#    EM_27             0.184    0.087    2.122    0.034    0.184    0.191
#    EM_28             0.476    0.075    6.365    0.000    0.476    0.546
#    EM_29             0.233    0.086    2.720    0.007    0.233    0.260
#    EM_30             0.106    0.131    0.807    0.420    0.106    0.100
#    EM_31             0.150    0.138    1.088    0.277    0.150    0.147
#    EM_32             0.045    0.098    0.456    0.648    0.045    0.042
#    EM_33             0.229    0.095    2.410    0.016    0.229    0.233
#    EM_34             0.565    0.219    2.582    0.010    0.565    0.628
#    EM_35             0.348    0.102    3.420    0.001    0.348    0.342
#    EM_36            -0.082    0.106   -0.781    0.435   -0.082   -0.076
#    EM_37             0.070    0.081    0.860    0.390    0.070    0.071
#    EM_38             0.647    0.080    8.111    0.000    0.647    0.795
#    EM_39             0.220    0.103    2.124    0.034    0.220    0.206

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.2700791      0.8205128      0.8897208      0.4588120 
#$FactorLevelIndices
#               ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.9730195 0.14052479 0.02698053 0.7661772 0.76592852 0.7943014 0.8924660
#describing    0.9517310 0.20877698 0.04826896 0.8697878 0.84144101 0.8664915 0.9319621
#nonjudging    0.8206954 0.16542100 0.17930463 0.8518810 0.71076397 0.8131239 0.9032699
#nonreactivity 0.8846684 0.14867584 0.11533161 0.8204455 0.74668261 0.8002662 0.8978773
#acting        0.2495276 0.06652228 0.75047236 0.8872452 0.03503782 0.5985774 0.8621141
#general       0.2700791 0.27007911 0.27007911 0.8897208 0.45881198 0.8832829 0.9444687






################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 2. Austria

library(haven)
mydata <- read_sav("austria.sav")
colnames(mydata)
mydata <- as.data.frame(mydata[,5:43])
mydata <- na.omit(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .84, omega T = .87

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 components
# Eigenvalue 1 = 6.00; eigenvalue 2 = 3.80

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 components
# Eigenvalue 1 = 6.77; eigenvalue 2 = 4.28

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.15, RMSEA=.102, RMSR=.12, TLI=.353

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.17, RMSEA=.117, RMSR=.14, TLI=.327

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities if response bias


# Give solution with 5 factors (EGA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.38, RMSEA=.043, RMSR=.03, TLI=.886
#      MR1   MR2   MR4  MR3   MR5
#MR1  1.00 -0.21 -0.28 0.31  0.10
#MR2 -0.21  1.00  0.41 0.06  0.00
#MR4 -0.28  0.41  1.00 0.04 -0.05
#MR3  0.31  0.06  0.04 1.00  0.16
#MR5  0.10  0.00 -0.05 0.16  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.056, RMSR=.03, TLI=.845
#      MR1   MR2  MR3   MR4   MR5
#MR1  1.00 -0.21 0.32 -0.28  0.11
#MR2 -0.21  1.00 0.05  0.41  0.01
#MR3  0.32  0.05 1.00  0.04  0.16
#MR4 -0.28  0.41 0.04  1.00 -0.05
#MR5  0.11  0.01 0.16 -0.05  1.00

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.690       0.538
#Tucker-Lewis Index (TLI)                       0.673       0.512
#Robust Comparative Fit Index (CFI)                         0.361
#Robust Tucker-Lewis Index (TLI)                            0.325
#RMSEA                                          0.167       0.127
#Robust RMSEA                                               0.119
#SRMR                                           0.142       0.142

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .191

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.395       0.387
#Tucker-Lewis Index (TLI)                       0.361       0.353
#Robust Comparative Fit Index (CFI)                         0.398
#Robust Tucker-Lewis Index (TLI)                            0.365
#RMSEA                                          0.102       0.093
#Robust RMSEA                                               0.101
#SRMR                                           0.125       0.125

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .086

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.949       0.904
#Tucker-Lewis Index (TLI)                       0.946       0.898
#Robust Comparative Fit Index (CFI)                         0.823
#Robust Tucker-Lewis Index (TLI)                            0.810
#RMSEA                                          0.068       0.058
#Robust RMSEA                                               0.063
#SRMR                                           0.066       0.066

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.362    0.031   11.649    0.000    0.362    0.362
#    nonjudging        0.141    0.035    4.065    0.000    0.141    0.141
#    nonreactivity     0.395    0.035   11.172    0.000    0.395    0.395
#    acting            0.037    0.034    1.090    0.276    0.037    0.037
#  describing ~~                                                         
#    nonjudging       -0.231    0.031   -7.355    0.000   -0.231   -0.231
#    nonreactivity     0.316    0.032    9.977    0.000    0.316    0.316
#    acting           -0.354    0.031  -11.410    0.000   -0.354   -0.354
#  nonjudging ~~                                                         
#    nonreactivity    -0.176    0.038   -4.686    0.000   -0.176   -0.176
#    acting            0.504    0.026   19.611    0.000    0.504    0.504
#  nonreactivity ~~                                                      
#    acting           -0.232    0.034   -6.744    0.000   -0.232   -0.232

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .429

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.852       0.859
#Tucker-Lewis Index (TLI)                       0.841       0.849
#Robust Comparative Fit Index (CFI)                         0.861
#Robust Tucker-Lewis Index (TLI)                            0.851
#RMSEA                                          0.051       0.045
#Robust RMSEA                                               0.049
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .363

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.337    0.038    8.810    0.000    0.337    0.337
#    nonjudging        0.115    0.042    2.716    0.007    0.115    0.115
#    nonreactivity     0.328    0.054    6.110    0.000    0.328    0.328
#    acting            0.064    0.042    1.505    0.132    0.064    0.064
#  describing ~~                                                         
#    nonjudging       -0.236    0.038   -6.188    0.000   -0.236   -0.236
#    nonreactivity     0.265    0.047    5.685    0.000    0.265    0.265
#    acting           -0.333    0.041   -8.076    0.000   -0.333   -0.333
#  nonjudging ~~                                                         
#    nonreactivity    -0.170    0.049   -3.487    0.000   -0.170   -0.170
#    acting            0.482    0.036   13.247    0.000    0.482    0.482
#  nonreactivity ~~                                                      
#    acting           -0.221    0.049   -4.542    0.000   -0.221   -0.221

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.813       0.819
#Tucker-Lewis Index (TLI)                       0.802       0.809
#Robust Comparative Fit Index (CFI)                         0.821
#Robust Tucker-Lewis Index (TLI)                            0.812
#RMSEA                                          0.057       0.051
#Robust RMSEA                                               0.055
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .360


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.890       0.897
#Tucker-Lewis Index (TLI)                       0.877       0.885
#Robust Comparative Fit Index (CFI)                         0.899
#Robust Tucker-Lewis Index (TLI)                            0.887
#RMSEA                                          0.045       0.039
#Robust RMSEA                                               0.043
#SRMR                                           0.069       0.069

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .364

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.2851314      0.8205128      0.7989964      0.3472705 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.9806285 0.15813603 0.01937145 0.7717591 0.77103716 0.7888087 0.8889505
#describing    0.8187504 0.21957225 0.18124955 0.3540645 0.35382606 0.8556421 0.9197876
#nonjudging    0.7091581 0.16068869 0.29084185 0.8569398 0.60718245 0.7877175 0.8796644
#nonreactivity 0.9168749 0.10224949 0.08312515 0.6770050 0.65390787 0.6730454 0.8217035
#acting        0.3193034 0.07422211 0.68069657 0.8292623 0.09575562 0.6435436 0.8375043
#general       0.2851314 0.28513143 0.28513143 0.7989964 0.34727051 0.8513233 0.9001514










################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 3. China

library(readxl)
China <- read_excel("China.xlsx")
colnames(China)
mydata <- as.data.frame(China)
mydata <- na.omit(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .84, omega T = .88

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 components
# Eigenvalue 1 = 7.11; eigenvalue 2 = 4.55

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5components
# Eigenvalue 1 = 7.91; eigenvalue 2 = 5.05

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.111, RMSR=.15, TLI=.353

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.14, RMSR=.16, TLI=.281

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EGA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.048, RMSR=.04, TLI=.873
#      MR1   MR3   MR4  MR2  MR5
#MR1  1.00 -0.27 -0.21 0.26 0.19
#MR3 -0.27  1.00  0.39 0.06 0.08
#MR4 -0.21  0.39  1.00 0.12 0.10
#MR2  0.26  0.06  0.12 1.00 0.20
#MR5  0.19  0.08  0.10 0.20 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.088, RMSR=.05, TLI=.712
#      MR1   MR3   MR4  MR2  MR5
#MR1  1.00 -0.27 -0.19 0.27 0.17
#MR3 -0.27  1.00  0.39 0.06 0.09
#MR4 -0.19  0.39  1.00 0.11 0.10
#MR2  0.27  0.06  0.11 1.00 0.18
#MR5  0.17  0.09  0.10 0.18 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.710       0.530
#Tucker-Lewis Index (TLI)                       0.693       0.504
#Robust Comparative Fit Index (CFI)                         0.322
#Robust Tucker-Lewis Index (TLI)                            0.285
#RMSEA                                          0.197       0.130
#Robust RMSEA                                               0.149
#SRMR                                           0.167       0.167

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .187

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.387       0.358
#Tucker-Lewis Index (TLI)                       0.353       0.322
#Robust Comparative Fit Index (CFI)                         0.385
#Robust Tucker-Lewis Index (TLI)                            0.351
#RMSEA                                          0.116       0.112
#Robust RMSEA                                               0.115
#SRMR                                           0.146       0.146

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .117

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.940       0.877
#Tucker-Lewis Index (TLI)                       0.936       0.868
#Robust Comparative Fit Index (CFI)                         0.708
#Robust Tucker-Lewis Index (TLI)                            0.687
#RMSEA                                          0.090       0.067
#Robust RMSEA                                               0.099
#SRMR                                           0.095       0.095

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.418    0.053    7.886    0.000    0.418    0.418
#    nonjudging        0.153    0.071    2.154    0.031    0.153    0.153
#    nonreactivity     0.536    0.072    7.431    0.000    0.536    0.536
#    acting            0.067    0.067    1.007    0.314    0.067    0.067
#  describing ~~                                                         
#    nonjudging       -0.272    0.057   -4.793    0.000   -0.272   -0.272
#    nonreactivity     0.564    0.059    9.526    0.000    0.564    0.564
#    acting           -0.332    0.056   -5.918    0.000   -0.332   -0.332
#  nonjudging ~~                                                         
#    nonreactivity    -0.243    0.078   -3.136    0.002   -0.243   -0.243
#    acting            0.523    0.051   10.242    0.000    0.523    0.523
#  nonreactivity ~~                                                      
#    acting           -0.344    0.065   -5.331    0.000   -0.344   -0.344

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .428

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.805       0.816
#Tucker-Lewis Index (TLI)                       0.791       0.803
#Robust Comparative Fit Index (CFI)                         0.819
#Robust Tucker-Lewis Index (TLI)                            0.807
#RMSEA                                          0.066       0.060
#Robust RMSEA                                               0.063
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .382

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.397    0.086    4.629    0.000    0.397    0.397
#    nonjudging        0.179    0.098    1.827    0.068    0.179    0.179
#    nonreactivity     0.487    0.098    4.989    0.000    0.487    0.487
#    acting            0.081    0.088    0.916    0.360    0.081    0.081
#  describing ~~                                                         
#    nonjudging       -0.277    0.092   -2.999    0.003   -0.277   -0.277
#    nonreactivity     0.497    0.105    4.747    0.000    0.497    0.497
#    acting           -0.308    0.086   -3.567    0.000   -0.308   -0.308
#  nonjudging ~~                                                         
#    nonreactivity    -0.195    0.108   -1.814    0.070   -0.195   -0.195
#    acting            0.496    0.072    6.849    0.000    0.496    0.496
#  nonreactivity ~~                                                      
#    acting           -0.227    0.112   -2.032    0.042   -0.227   -0.227

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.760       0.770
#Tucker-Lewis Index (TLI)                       0.746       0.757
#Robust Comparative Fit Index (CFI)                         0.773
#Robust Tucker-Lewis Index (TLI)                            0.761
#RMSEA                                          0.073       0.067
#Robust RMSEA                                               0.070
#SRMR                                           0.148       0.148

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .353


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.857       0.862
#Tucker-Lewis Index (TLI)                       0.840       0.846
#Robust Comparative Fit Index (CFI)                         0.867
#Robust Tucker-Lewis Index (TLI)                            0.851
#RMSEA                                          0.058       0.053
#Robust RMSEA                                               0.055
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .424

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing =~                                                          
#    EM_1              0.405    0.068    5.930    0.000    0.405    0.465
#    EM_6              0.567    0.079    7.152    0.000    0.567    0.566
#    EM_11             0.509    0.077    6.595    0.000    0.509    0.479
#    EM_15             0.810    0.062   13.150    0.000    0.810    0.809
#    EM_20             0.701    0.069   10.149    0.000    0.701    0.699
#    EM_26             0.542    0.078    6.949    0.000    0.542    0.506
#    EM_31             0.540    0.072    7.541    0.000    0.540    0.540
#    EM_36             0.468    0.078    6.009    0.000    0.468    0.485
#  describing =~                                                         
#    EM_2             -0.749    0.055  -13.511    0.000   -0.749   -0.788
#    EM_7             -0.658    0.061  -10.724    0.000   -0.658   -0.692
#    EM_12             0.506    0.082    6.165    0.000    0.506    0.598
#    EM_16             0.529    0.082    6.431    0.000    0.529    0.606
#    EM_22             0.440    0.083    5.312    0.000    0.440    0.505
#    EM_27            -0.654    0.058  -11.274    0.000   -0.654   -0.728
#    EM_32            -0.767    0.051  -15.123    0.000   -0.767   -0.807
#    EM_37            -0.714    0.052  -13.743    0.000   -0.714   -0.792
#  nonjudging =~                                                         
#    EM_3              0.549    0.084    6.564    0.000    0.549    0.538
#    EM_10             0.497    0.098    5.086    0.000    0.497    0.528
#    EM_14             0.612    0.112    5.456    0.000    0.612    0.535
#    EM_17             0.634    0.090    7.059    0.000    0.634    0.612
#    EM_25             0.391    0.099    3.961    0.000    0.391    0.393
#    EM_30             0.640    0.084    7.609    0.000    0.640    0.595
#    EM_35             0.183    0.085    2.149    0.032    0.183    0.216
#    EM_39             0.239    0.091    2.629    0.009    0.239    0.243
#  nonreactivity =~                                                      
#    EM_4              0.221    0.085    2.607    0.009    0.221    0.264
#    EM_9              0.395    0.075    5.238    0.000    0.395    0.448
#    EM_19             0.436    0.090    4.842    0.000    0.436    0.501
#    EM_21             0.454    0.071    6.365    0.000    0.454    0.518
#    EM_24             0.349    0.081    4.296    0.000    0.349    0.406
#    EM_29             0.341    0.088    3.890    0.000    0.341    0.410
#    EM_33             0.532    0.072    7.387    0.000    0.532    0.606
#  acting =~                                                             
#    EM_5              0.834    0.077   10.816    0.000    0.834    0.774
#    EM_8              0.286    0.099    2.894    0.004    0.286    0.263
#    EM_13             0.757    0.068   11.061    0.000    0.757    0.749
#    EM_18             0.480    0.064    7.536    0.000    0.480    0.506
#    EM_23             0.109    0.110    0.986    0.324    0.109    0.108
#    EM_28             0.363    0.122    2.974    0.003    0.363    0.347
#    EM_34             0.099    0.123    0.807    0.420    0.099    0.107
#    EM_38             0.558    0.118    4.723    0.000    0.558    0.543
#  general =~                                                            
#    EM_1             -0.031    0.086   -0.358    0.720   -0.031   -0.036
#    EM_2              0.188    0.120    1.567    0.117    0.188    0.198
#    EM_3             -0.460    0.078   -5.928    0.000   -0.460   -0.451
#    EM_4             -0.273    0.072   -3.796    0.000   -0.273   -0.327
#    EM_5             -0.411    0.116   -3.530    0.000   -0.411   -0.381
#    EM_6             -0.159    0.098   -1.634    0.102   -0.159   -0.159
#    EM_7              0.174    0.108    1.608    0.108    0.174    0.183
#    EM_8             -0.571    0.098   -5.828    0.000   -0.571   -0.525
#    EM_9              0.270    0.094    2.867    0.004    0.270    0.306
#    EM_10            -0.305    0.100   -3.062    0.002   -0.305   -0.324
#    EM_11             0.013    0.109    0.123    0.902    0.013    0.013
#    EM_12            -0.399    0.121   -3.291    0.001   -0.399   -0.472
#    EM_13            -0.449    0.103   -4.368    0.000   -0.449   -0.444
#    EM_14            -0.636    0.124   -5.137    0.000   -0.636   -0.556
#    EM_15            -0.070    0.100   -0.695    0.487   -0.070   -0.070
#    EM_16            -0.443    0.122   -3.636    0.000   -0.443   -0.507
#    EM_17            -0.222    0.089   -2.501    0.012   -0.222   -0.215
#    EM_18            -0.517    0.071   -7.271    0.000   -0.517   -0.545
#    EM_19             0.162    0.101    1.604    0.109    0.162    0.186
#    EM_20            -0.122    0.108   -1.129    0.259   -0.122   -0.122
#    EM_21             0.052    0.103    0.508    0.611    0.052    0.060
#    EM_22            -0.470    0.113   -4.175    0.000   -0.470   -0.539
#    EM_23            -0.674    0.099   -6.819    0.000   -0.674   -0.671
#    EM_24             0.223    0.074    2.996    0.003    0.223    0.260
#    EM_25            -0.391    0.089   -4.382    0.000   -0.391   -0.393
#    EM_26             0.007    0.112    0.063    0.950    0.007    0.007
#    EM_27             0.227    0.103    2.209    0.027    0.227    0.252
#    EM_28            -0.658    0.109   -6.027    0.000   -0.658   -0.629
#    EM_29            -0.207    0.084   -2.471    0.013   -0.207   -0.248
#    EM_30            -0.469    0.103   -4.557    0.000   -0.469   -0.436
#    EM_31             0.017    0.105    0.163    0.871    0.017    0.017
#    EM_32             0.013    0.117    0.113    0.910    0.013    0.014
#    EM_33            -0.061    0.088   -0.690    0.490   -0.061   -0.069
#    EM_34            -0.596    0.108   -5.505    0.000   -0.596   -0.642
#    EM_35            -0.343    0.084   -4.066    0.000   -0.343   -0.404
#    EM_36             0.010    0.109    0.091    0.927    0.010    0.010
#    EM_37             0.153    0.116    1.328    0.184    0.153    0.170
#    EM_38            -0.495    0.122   -4.055    0.000   -0.495   -0.481
#    EM_39            -0.455    0.089   -5.109    0.000   -0.455   -0.463

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.3029755      0.8205128      0.8562354      0.4645339 
#$FactorLevelIndices
#               ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#observing     0.9828700 0.15852188 0.01713004 0.7982159 0.7937886 0.8322639 0.9128668
#describing    0.8060077 0.22889195 0.19399232 0.6063832 0.5455050 0.8970249 0.9445208
#nonjudging    0.5707653 0.10850769 0.42923474 0.8336044 0.4670619 0.7249210 0.8449204
#nonreactivity 0.8002821 0.08781324 0.19971787 0.6601888 0.6583428 0.6678342 0.8232258
#acting        0.4444547 0.11328971 0.55554530 0.8916069 0.3408628 0.7903981 0.8905046
#general       0.3029755 0.30297552 0.30297552 0.8562354 0.4645339 0.8767996 0.9103451






################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 4. Chile

library(haven)
Chile <- read_sav("Chile.sav")
colnames(Chile)
mydata <- as.data.frame(Chile[,22:60])
mydata[mydata == 0] <- NA
mydata <- na.omit(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .81, omega T = .85

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 5.44; eigenvalue 2 = 3.44

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 6.03; eigenvalue 2 = 3.90

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.096, RMSR=.12, TLI=.343

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.15, RMSEA=.114, RMSR=.14, TLI=.297

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EGA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.047, RMSR=.04, TLI=.844
#      MR1   MR3   MR4  MR2   MR5
#MR1  1.00  0.29 -0.14 0.24  0.14
#MR3  0.29  1.00 -0.33 0.13 -0.04
#MR4 -0.14 -0.33  1.00 0.13  0.13
#MR2  0.24  0.13  0.13 1.00  0.12
#MR5  0.14 -0.04  0.13 0.12  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.066, RMSR=.04, TLI=.759
#      MR1   MR3   MR4  MR2   MR5
#MR1  1.00  0.29 -0.14 0.24  0.14
#MR3  0.29  1.00 -0.34 0.13 -0.04
#MR4 -0.14 -0.34  1.00 0.14  0.13
#MR2  0.24  0.13  0.14 1.00  0.12
#MR5  0.14 -0.04  0.13 0.12  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.651       0.484
#Tucker-Lewis Index (TLI)                       0.631       0.456
#Robust Comparative Fit Index (CFI)                         0.343
#Robust Tucker-Lewis Index (TLI)                            0.306
#RMSEA                                          0.152       0.111
#Robust RMSEA                                               0.116
#SRMR                                           0.135       0.135

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .098

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.376       0.353
#Tucker-Lewis Index (TLI)                       0.342       0.317
#Robust Comparative Fit Index (CFI)                         0.376
#Robust Tucker-Lewis Index (TLI)                            0.341
#RMSEA                                          0.098       0.094
#Robust RMSEA                                               0.097
#SRMR                                           0.118       0.118

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .071

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.907       0.840
#Tucker-Lewis Index (TLI)                       0.900       0.829
#Robust Comparative Fit Index (CFI)                         0.732
#Robust Tucker-Lewis Index (TLI)                            0.713
#RMSEA                                          0.079       0.062
#Robust RMSEA                                               0.075
#SRMR                                           0.083       0.083

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.174    0.060    2.922    0.003    0.174    0.174
#    nonjudging        0.372    0.050    7.410    0.000    0.372    0.372
#    nonreactivity     0.391    0.060    6.564    0.000    0.391    0.391
#    acting            0.115    0.054    2.121    0.034    0.115    0.115
#  describing ~~                                                         
#    nonjudging       -0.258    0.052   -4.946    0.000   -0.258   -0.258
#    nonreactivity     0.200    0.062    3.242    0.001    0.200    0.200
#    acting           -0.456    0.045  -10.232    0.000   -0.456   -0.456
#  nonjudging ~~                                                         
#    nonreactivity     0.149    0.060    2.488    0.013    0.149    0.149
#    acting            0.350    0.047    7.412    0.000    0.350    0.350
#  nonreactivity ~~                                                      
#    acting           -0.109    0.056   -1.949    0.051   -0.109   -0.109

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .362

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.782       0.786
#Tucker-Lewis Index (TLI)                       0.767       0.770
#Robust Comparative Fit Index (CFI)                         0.791
#Robust Tucker-Lewis Index (TLI)                            0.776
#RMSEA                                          0.058       0.055
#Robust RMSEA                                               0.057
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .301

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.155    0.102    1.522    0.128    0.155    0.155
#    nonjudging        0.325    0.065    5.015    0.000    0.325    0.325
#    nonreactivity     0.356    0.091    3.914    0.000    0.356    0.356
#    acting            0.129    0.071    1.826    0.068    0.129    0.129
#  describing ~~                                                         
#    nonjudging       -0.259    0.078   -3.311    0.001   -0.259   -0.259
#    nonreactivity     0.188    0.110    1.709    0.087    0.188    0.188
#    acting           -0.444    0.066   -6.694    0.000   -0.444   -0.444
#  nonjudging ~~                                                         
#    nonreactivity     0.125    0.090    1.385    0.166    0.125    0.125
#    acting            0.350    0.060    5.837    0.000    0.350    0.350
#  nonreactivity ~~                                                      
#    acting           -0.105    0.077   -1.355    0.176   -0.105   -0.105

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.744       0.748
#Tucker-Lewis Index (TLI)                       0.730       0.734
#Robust Comparative Fit Index (CFI)                         0.753
#Robust Tucker-Lewis Index (TLI)                            0.740
#RMSEA                                          0.063       0.059
#Robust RMSEA                                               0.061
#SRMR                                           0.112       0.112

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .300


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.810       0.806
#Tucker-Lewis Index (TLI)                       0.788       0.784
#Robust Comparative Fit Index (CFI)                         0.815
#Robust Tucker-Lewis Index (TLI)                            0.793
#RMSEA                                          0.056       0.053
#Robust RMSEA                                               0.055
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .331

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.2280439      0.8205128      0.8620383      0.5375660 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.5025650 0.08933932 0.49743500 0.7526281 0.34909998 0.6196681 0.7656550
#describing    0.8471809 0.18369239 0.15281907 0.4213966 0.07489212 0.7948309 0.8981669
#nonjudging    0.7302116 0.18877321 0.26978841 0.8591825 0.62495632 0.8149443 0.8861824
#nonreactivity 0.7667791 0.08480160 0.23322091 0.6238643 0.47795259 0.6185183 0.7819476
#acting        0.9536896 0.22534964 0.04631037 0.8412628 0.80802532 0.8437870 0.9166792
#general       0.2280439 0.22804385 0.22804385 0.8620383 0.53756602 0.7860163 0.8470513






################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 5. Croatia

library(haven)
croatia <- read_sav("croatia.sav")
mydata <- as.data.frame(croatia[,4:42])
mydata <- na.omit(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .82, omega T = .87

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 6.26; eigenvalue 2 = 3.53

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 6 components
# Eigenvalue 1 = 7.03; eigenvalue 2 = 4.00

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.16, RMSEA=.107, RMSR=.13, TLI=.328

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.18, RMSEA=.134, RMSR=.15, TLI=.264

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.052, RMSR=.05, TLI=.838
#      MR2   MR1  MR4   MR3   MR5
#MR2  1.00 -0.18 0.03  0.23 -0.18
#MR1 -0.18  1.00 0.20 -0.20  0.34
#MR4  0.03  0.20 1.00  0.02  0.11
#MR3  0.23 -0.20 0.02  1.00 -0.28
#MR5 -0.18  0.34 0.11 -0.28  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.085, RMSR=.05, TLI=.699
#      MR2   MR1  MR4   MR5   MR3
#MR2  1.00 -0.18 0.04 -0.19  0.23
#MR1 -0.18  1.00 0.19  0.34 -0.18
#MR4  0.04  0.19 1.00  0.10  0.04
#MR5 -0.19  0.34 0.10  1.00 -0.28
#MR3  0.23 -0.18 0.04 -0.28  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.686       0.548
#Tucker-Lewis Index (TLI)                       0.668       0.522
#Robust Comparative Fit Index (CFI)                         0.318
#Robust Tucker-Lewis Index (TLI)                            0.280
#RMSEA                                          0.182       0.119
#Robust RMSEA                                               0.139
#SRMR                                           0.153       0.153

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .229

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.360       0.349
#Tucker-Lewis Index (TLI)                       0.325       0.313
#Robust Comparative Fit Index (CFI)                         0.364
#Robust Tucker-Lewis Index (TLI)                            0.329
#RMSEA                                          0.112       0.104
#Robust RMSEA                                               0.110
#SRMR                                           0.131       0.131

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .146

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.924       0.873
#Tucker-Lewis Index (TLI)                       0.919       0.864
#Robust Comparative Fit Index (CFI)                         0.715
#Robust Tucker-Lewis Index (TLI)                            0.695
#RMSEA                                          0.090       0.063
#Robust RMSEA                                               0.090
#SRMR                                           0.094       0.094

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.272    0.060    4.578    0.000    0.272    0.272
#    nonjudging        0.070    0.062    1.124    0.261    0.070    0.070
#    nonreactivity     0.361    0.063    5.740    0.000    0.361    0.361
#    acting            0.027    0.065    0.413    0.679    0.027    0.027
#  describing ~~                                                         
#    nonjudging       -0.270    0.060   -4.501    0.000   -0.270   -0.270
#    nonreactivity     0.209    0.058    3.591    0.000    0.209    0.209
#    acting           -0.317    0.053   -5.962    0.000   -0.317   -0.317
#  nonjudging ~~                                                         
#    nonreactivity    -0.150    0.066   -2.267    0.023   -0.150   -0.150
#    acting            0.474    0.049    9.721    0.000    0.474    0.474
#  nonreactivity ~~                                                      
#    acting           -0.310    0.054   -5.753    0.000   -0.310   -0.310

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .466

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.787       0.801
#Tucker-Lewis Index (TLI)                       0.772       0.787
#Robust Comparative Fit Index (CFI)                         0.805
#Robust Tucker-Lewis Index (TLI)                            0.791
#RMSEA                                          0.065       0.058
#Robust RMSEA                                               0.061
#SRMR                                           0.085       0.085

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .399

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.251    0.097    2.582    0.010    0.251    0.251
#    nonjudging        0.033    0.084    0.393    0.694    0.033    0.033
#    nonreactivity     0.320    0.100    3.207    0.001    0.320    0.320
#    acting            0.058    0.092    0.628    0.530    0.058    0.058
#  describing ~~                                                         
#    nonjudging       -0.250    0.088   -2.828    0.005   -0.250   -0.250
#    nonreactivity     0.204    0.086    2.371    0.018    0.204    0.204
#    acting           -0.329    0.087   -3.796    0.000   -0.329   -0.329
#  nonjudging ~~                                                         
#    nonreactivity    -0.176    0.101   -1.752    0.080   -0.176   -0.176
#    acting            0.473    0.073    6.472    0.000    0.473    0.473
#  nonreactivity ~~                                                      
#    acting           -0.353    0.089   -3.960    0.000   -0.353   -0.353

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.754       0.767
#Tucker-Lewis Index (TLI)                       0.740       0.755
#Robust Comparative Fit Index (CFI)                         0.772
#Robust Tucker-Lewis Index (TLI)                            0.759
#RMSEA                                          0.069       0.062
#Robust RMSEA                                               0.066
#SRMR                                           0.127       0.127

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .390


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.858       0.869
#Tucker-Lewis Index (TLI)                       0.842       0.853
#Robust Comparative Fit Index (CFI)                         0.874
#Robust Tucker-Lewis Index (TLI)                            0.859
#RMSEA                                          0.054       0.048
#Robust RMSEA                                               0.050
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .416

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3102273      0.8205128      0.8213029      0.4138322 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.9716772 0.13114593 0.02832278 0.7366191 0.73400506 0.7592038 0.8724737
#describing    0.8653646 0.21073168 0.13463536 0.3754832 0.36181144 0.8634823 0.9318566
#nonjudging    0.7258954 0.16526910 0.27410464 0.8691360 0.63818203 0.8125993 0.9078169
#nonreactivity 0.8681792 0.10734634 0.13182076 0.7177289 0.66473053 0.7110183 0.8492521
#acting        0.2786155 0.07527969 0.72138448 0.8736416 0.02220557 0.6159620 0.8529891
#general       0.3102273 0.31022726 0.31022726 0.8213029 0.41383218 0.8894819 0.9457382





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 6. Germany

library(haven)
Germany <- read_sav("Germany.sav")
colnames(Germany)
mydata <- as.data.frame(Germany[,16:54])
mydata[mydata==6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .88, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 6 components
# Eigenvalue 1 = 7.65; eigenvalue 2 = 2.67

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 6 components
# Eigenvalue 1 = 8.64; eigenvalue 2 = 3.06

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.20, RMSEA=.106, RMSR=.11, TLI=.417

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.22, RMSEA=.128, RMSR=.13, TLI=.365

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# error message


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.049, RMSR=.03, TLI=.875
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.36  0.32 -0.04 -0.29
#MR2 -0.36  1.00 -0.30  0.17  0.18
#MR4  0.32 -0.30  1.00 -0.02 -0.27
#MR3 -0.04  0.17 -0.02  1.00  0.10
#MR5 -0.29  0.18 -0.27  0.10  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.072, RMSR=.04, TLI=.801
#      MR1   MR2   MR4   MR3   MR5
#MR1  1.00 -0.36  0.32 -0.04 -0.29
#MR2 -0.36  1.00 -0.31  0.19  0.19
#MR4  0.32 -0.31  1.00 -0.02 -0.26
#MR3 -0.04  0.19 -0.02  1.00  0.10
#MR5 -0.29  0.19 -0.26  0.10  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.786       0.594
#Tucker-Lewis Index (TLI)                       0.774       0.572
#Robust Comparative Fit Index (CFI)                         0.385
#Robust Tucker-Lewis Index (TLI)                            0.351
#RMSEA                                          0.165       0.129
#Robust RMSEA                                               0.132
#SRMR                                           0.137       0.137

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .265

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.450       0.447
#Tucker-Lewis Index (TLI)                       0.419       0.416
#Robust Comparative Fit Index (CFI)                         0.454
#Robust Tucker-Lewis Index (TLI)                            0.423
#RMSEA                                          0.107       0.100
#Robust RMSEA                                               0.106
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .175

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.960       0.910
#Tucker-Lewis Index (TLI)                       0.957       0.904
#Robust Comparative Fit Index (CFI)                         0.802
#Robust Tucker-Lewis Index (TLI)                            0.788
#RMSEA                                          0.072       0.061
#Robust RMSEA                                               0.075
#SRMR                                           0.075       0.075

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .468

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.199    0.045    4.435    0.000    0.199    0.199
#    nonjudging       -0.033    0.046   -0.718    0.473   -0.033   -0.033
#    nonreactivity     0.163    0.046    3.546    0.000    0.163    0.163
#    acting           -0.034    0.046   -0.743    0.458   -0.034   -0.034
#  describing ~~                                                         
#    nonjudging       -0.404    0.037  -10.833    0.000   -0.404   -0.404
#    nonreactivity     0.329    0.041    7.983    0.000    0.329    0.329
#    acting           -0.408    0.037  -11.083    0.000   -0.408   -0.408
#  nonjudging ~~                                                         
#    nonreactivity    -0.433    0.039  -11.163    0.000   -0.433   -0.433
#    acting            0.472    0.037   12.790    0.000    0.472    0.472
#  nonreactivity ~~                                                      
#    acting           -0.365    0.040   -9.155    0.000   -0.365   -0.365

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.849       0.854
#Tucker-Lewis Index (TLI)                       0.838       0.843
#Robust Comparative Fit Index (CFI)                         0.857
#Robust Tucker-Lewis Index (TLI)                            0.847
#RMSEA                                          0.057       0.052
#Robust RMSEA                                               0.055
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .412

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.167    0.055    3.062    0.002    0.167    0.167
#    nonjudging       -0.056    0.059   -0.953    0.340   -0.056   -0.056
#    nonreactivity     0.135    0.058    2.324    0.020    0.135    0.135
#    acting           -0.005    0.058   -0.082    0.935   -0.005   -0.005
#  describing ~~                                                         
#    nonjudging       -0.405    0.046   -8.882    0.000   -0.405   -0.405
#    nonreactivity     0.298    0.054    5.538    0.000    0.298    0.298
#    acting           -0.366    0.054   -6.745    0.000   -0.366   -0.366
#  nonjudging ~~                                                         
#    nonreactivity    -0.433    0.048   -9.030    0.000   -0.433   -0.433
#    acting            0.408    0.065    6.248    0.000    0.408    0.408
#  nonreactivity ~~                                                      
#    acting           -0.384    0.054   -7.055    0.000   -0.384   -0.384

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.813       0.817
#Tucker-Lewis Index (TLI)                       0.802       0.807
#Robust Comparative Fit Index (CFI)                         0.821
#Robust Tucker-Lewis Index (TLI)                            0.811
#RMSEA                                          0.063       0.058
#Robust RMSEA                                               0.061
#SRMR                                           0.141       0.141


fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .406


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.906       0.912
#Tucker-Lewis Index (TLI)                       0.895       0.901
#Robust Comparative Fit Index (CFI)                         0.914
#Robust Tucker-Lewis Index (TLI)                            0.904
#RMSEA                                          0.046       0.041
#Robust RMSEA                                               0.043
#SRMR                                           0.059       0.059

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .418

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing =~                                                          
#    EM_1              0.437    0.051    8.627    0.000    0.437    0.439
#    EM_6              0.461    0.058    7.896    0.000    0.461    0.460
#    EM_11             0.513    0.055    9.260    0.000    0.513    0.412
#    EM_15             0.681    0.046   14.802    0.000    0.681    0.696
#    EM_20             0.541    0.052   10.335    0.000    0.541    0.536
#    EM_26             0.586    0.049   11.970    0.000    0.586    0.616
#    EM_31             0.580    0.050   11.534    0.000    0.580    0.594
#    EM_36             0.349    0.045    7.761    0.000    0.349    0.388
#  describing =~                                                         
#    EM_2              0.655    0.036   18.125    0.000    0.655    0.716
#    EM_7              0.464    0.046   10.087    0.000    0.464    0.523
#    EM_12            -0.601    0.037  -16.369    0.000   -0.601   -0.658
#    EM_16            -0.675    0.041  -16.617    0.000   -0.675   -0.723
#    EM_22            -0.412    0.042   -9.771    0.000   -0.412   -0.470
#    EM_27             0.507    0.042   11.980    0.000    0.507    0.551
#    EM_32             0.500    0.042   11.813    0.000    0.500    0.535
#    EM_37             0.543    0.043   12.682    0.000    0.543    0.596
#  nonjudging =~                                                         
#    EM_3              0.547    0.049   11.230    0.000    0.547    0.538
#    EM_10             0.620    0.052   12.014    0.000    0.620    0.588
#    EM_14             0.547    0.049   11.160    0.000    0.547    0.542
#    EM_17             0.621    0.051   12.216    0.000    0.621    0.562
#    EM_25             0.746    0.042   17.768    0.000    0.746    0.731
#    EM_30             0.688    0.047   14.674    0.000    0.688    0.666
#    EM_35             0.537    0.051   10.530    0.000    0.537    0.498
#    EM_39             0.475    0.052    9.195    0.000    0.475    0.470
#  nonreactivity =~                                                      
#    EM_4              0.391    0.042    9.199    0.000    0.391    0.457
#    EM_9              0.338    0.046    7.338    0.000    0.338    0.378
#    EM_19             0.592    0.043   13.837    0.000    0.592    0.627
#    EM_21             0.313    0.044    7.139    0.000    0.313    0.371
#    EM_24             0.458    0.046   10.062    0.000    0.458    0.506
#    EM_29             0.653    0.043   15.339    0.000    0.653    0.730
#    EM_33             0.504    0.040   12.476    0.000    0.504    0.611
#  acting =~                                                             
#    EM_5              0.667    0.053   12.636    0.000    0.667    0.701
#    EM_8              0.225    0.055    4.107    0.000    0.225    0.232
#    EM_13             0.697    0.046   15.211    0.000    0.697    0.749
#    EM_18             0.376    0.055    6.810    0.000    0.376    0.401
#    EM_23            -0.225    0.053   -4.231    0.000   -0.225   -0.240
#    EM_28            -0.042    0.055   -0.752    0.452   -0.042   -0.051
#    EM_34            -0.165    0.056   -2.921    0.003   -0.165   -0.203
#    EM_38            -0.093    0.066   -1.402    0.161   -0.093   -0.104
#  general =~                                                            
#    EM_1              0.093    0.058    1.608    0.108    0.093    0.094
#    EM_2              0.282    0.049    5.738    0.000    0.282    0.308
#    EM_3             -0.382    0.053   -7.204    0.000   -0.382   -0.376
#    EM_4              0.115    0.047    2.453    0.014    0.115    0.134
#    EM_5             -0.462    0.068   -6.773    0.000   -0.462   -0.485
#    EM_6              0.134    0.053    2.530    0.011    0.134    0.134
#    EM_7              0.336    0.045    7.457    0.000    0.336    0.379
#    EM_8             -0.584    0.045  -12.927    0.000   -0.584   -0.603
#    EM_9              0.314    0.049    6.391    0.000    0.314    0.350
#    EM_10            -0.455    0.055   -8.273    0.000   -0.455   -0.431
#    EM_11            -0.185    0.069   -2.668    0.008   -0.185   -0.149
#    EM_12            -0.399    0.048   -8.219    0.000   -0.399   -0.437
#    EM_13            -0.485    0.066   -7.396    0.000   -0.485   -0.522
#    EM_14            -0.485    0.054   -9.046    0.000   -0.485   -0.480
#    EM_15             0.116    0.057    2.030    0.042    0.116    0.118
#    EM_16            -0.344    0.049   -7.060    0.000   -0.344   -0.368
#    EM_17            -0.387    0.057   -6.746    0.000   -0.387   -0.350
#    EM_18            -0.578    0.050  -11.521    0.000   -0.578   -0.616
#    EM_19             0.372    0.049    7.572    0.000    0.372    0.394
#    EM_20            -0.033    0.058   -0.582    0.561   -0.033   -0.033
#    EM_21             0.220    0.043    5.075    0.000    0.220    0.261
#    EM_22            -0.407    0.043   -9.493    0.000   -0.407   -0.464
#    EM_23            -0.457    0.046   -9.828    0.000   -0.457   -0.487
#    EM_24             0.246    0.048    5.137    0.000    0.246    0.271
#    EM_25            -0.479    0.052   -9.128    0.000   -0.479   -0.469
#    EM_26             0.097    0.055    1.771    0.077    0.097    0.102
#    EM_27             0.301    0.045    6.645    0.000    0.301    0.328
#    EM_28            -0.462    0.041  -11.387    0.000   -0.462   -0.572
#    EM_29             0.162    0.049    3.316    0.001    0.162    0.181
#    EM_30            -0.546    0.055   -9.997    0.000   -0.546   -0.528
#    EM_31             0.068    0.057    1.195    0.232    0.068    0.070
#    EM_32             0.213    0.052    4.132    0.000    0.213    0.228
#    EM_33             0.118    0.044    2.679    0.007    0.118    0.144
#    EM_34            -0.508    0.039  -12.931    0.000   -0.508   -0.625
#    EM_35            -0.490    0.052   -9.366    0.000   -0.490   -0.454
#    EM_36             0.055    0.055    1.011    0.312    0.055    0.062
#    EM_37             0.323    0.049    6.596    0.000    0.323    0.354
#    EM_38            -0.463    0.049   -9.427    0.000   -0.463   -0.517
#    EM_39            -0.407    0.053   -7.631    0.000   -0.407   -0.403

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.3346571      0.8205128      0.7936647      0.2783811 
#
#$FactorLevelIndices
#               ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.9640979 0.13182984 0.03590208 0.7527335 0.74586835 0.7721424 0.8796350
#describing    0.7319419 0.17188947 0.26805806 0.2369951 0.21661381 0.8321006 0.9081301
#nonjudging    0.6348973 0.15918008 0.36510272 0.8985676 0.56966517 0.8138817 0.8974834
#nonreactivity 0.8064737 0.12088907 0.19352627 0.7875127 0.64432303 0.7684200 0.8775044
#acting        0.3580020 0.08155445 0.64199801 0.8402040 0.08499351 0.7230605 0.8935834
#general       0.3346571 0.33465710 0.33465710 0.7936647 0.27838109 0.8825456 0.9126358





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 7. Hong Kong

library(haven)
Hong_Kong <- read_sav("Hong Kong.sav")
colnames(Hong_Kong)
mydata <- as.data.frame(Hong_Kong[,17:55])
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .80, omega T = .85

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 5 components
# Eigenvalue 1 = 5.60; eigenvalue 2 = 3.76

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 5 components
# Eigenvalue 1 = 6.34; eigenvalue 2 = 4.27

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.095, RMSR=.13, TLI=.355

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.16, RMSEA=.115, RMSR=.14, TLI=.308

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.041, RMSR=.04, TLI=.881
#      MR1   MR2   MR4  MR3   MR5
#MR1  1.00 -0.20  0.42 0.14 -0.03
#MR2 -0.20  1.00 -0.06 0.16  0.07
#MR4  0.42 -0.06  1.00 0.27  0.17
#MR3  0.14  0.16  0.27 1.00  0.30
#MR5 -0.03  0.07  0.17 0.30  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.064, RMSR=.04, TLI=.782
#      MR1   MR2   MR4  MR3   MR5
#MR1  1.00 -0.20  0.42 0.15 -0.03
#MR2 -0.20  1.00 -0.07 0.16  0.08
#MR4  0.42 -0.07  1.00 0.27  0.17
#MR3  0.15  0.16  0.27 1.00  0.31
#MR5 -0.03  0.08  0.17 0.31  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.612       0.424
#Tucker-Lewis Index (TLI)                       0.591       0.392
#Robust Comparative Fit Index (CFI)                         0.353
#Robust Tucker-Lewis Index (TLI)                            0.317
#RMSEA                                          0.166       0.118
#Robust RMSEA                                               0.118
#SRMR                                           0.144       0.144

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .126

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.389       0.384
#Tucker-Lewis Index (TLI)                       0.355       0.349
#Robust Comparative Fit Index (CFI)                         0.394
#Robust Tucker-Lewis Index (TLI)                            0.361
#RMSEA                                          0.098       0.091
#Robust RMSEA                                               0.096
#SRMR                                           0.124       0.124

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .093

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.821       0.713
#Tucker-Lewis Index (TLI)                       0.808       0.693
#Robust Comparative Fit Index (CFI)                         0.660
#Robust Tucker-Lewis Index (TLI)                            0.636
#RMSEA                                          0.114       0.084
#Robust RMSEA                                               0.086
#SRMR                                           0.112       0.112

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .336

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.214    0.052    4.086    0.000    0.214    0.214
#    nonjudging        0.341    0.054    6.306    0.000    0.341    0.341
#    nonreactivity     0.215    0.093    2.302    0.021    0.215    0.215
#    acting            0.192    0.056    3.418    0.001    0.192    0.192
#  describing ~~                                                         
#    nonjudging       -0.152    0.053   -2.854    0.004   -0.152   -0.152
#    nonreactivity    -0.511    0.087   -5.871    0.000   -0.511   -0.511
#    acting           -0.287    0.048   -5.967    0.000   -0.287   -0.287
#  nonjudging ~~                                                         
#    nonreactivity     0.690    0.089    7.774    0.000    0.690    0.690
#    acting            0.546    0.043   12.638    0.000    0.546    0.546
#  nonreactivity ~~                                                      
#    acting            0.829    0.071   11.738    0.000    0.829    0.829

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.738       0.747
#Tucker-Lewis Index (TLI)                       0.719       0.730
#Robust Comparative Fit Index (CFI)                         0.752
#Robust Tucker-Lewis Index (TLI)                            0.734
#RMSEA                                          0.065       0.059
#Robust RMSEA                                               0.062
#SRMR                                           0.102       0.102

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .271

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.284    0.077    3.687    0.000    0.284    0.284
#    nonjudging        0.283    0.078    3.635    0.000    0.283    0.283
#    nonreactivity     0.456    0.090    5.072    0.000    0.456    0.456
#    acting            0.161    0.074    2.167    0.030    0.161    0.161
#  describing ~~                                                         
#    nonjudging       -0.063    0.082   -0.764    0.445   -0.063   -0.063
#    nonreactivity     0.238    0.130    1.835    0.067    0.238    0.238
#    acting           -0.189    0.075   -2.519    0.012   -0.189   -0.189
#  nonjudging ~~                                                         
#    nonreactivity     0.099    0.181    0.543    0.587    0.099    0.099
#    acting            0.537    0.054    9.902    0.000    0.537    0.537
#  nonreactivity ~~                                                      
#    acting           -0.115    0.203   -0.567    0.571   -0.115   -0.115

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.696       0.708
#Tucker-Lewis Index (TLI)                       0.679       0.692
#Robust Comparative Fit Index (CFI)                         0.711
#Robust Tucker-Lewis Index (TLI)                            0.694
#RMSEA                                          0.069       0.063
#Robust RMSEA                                               0.066
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .283


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.835       0.845
#Tucker-Lewis Index (TLI)                       0.816       0.827
#Robust Comparative Fit Index (CFI)                         0.849
#Robust Tucker-Lewis Index (TLI)                            0.831
#RMSEA                                          0.052       0.047
#Robust RMSEA                                               0.049
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .343

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3318778      0.8205128      0.8689634      0.6181534 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#observing     0.8054129 0.13463494 0.1945871 0.7535626 0.6160899 0.7206286 0.8459624
#describing    0.7712100 0.21641779 0.2287900 0.6503876 0.4585870 0.8535272 0.9254777
#nonjudging    0.5588460 0.10045638 0.4411540 0.7753667 0.4421179 0.6325659 0.7807410
#nonreactivity 0.7390945 0.09060666 0.2609055 0.5999394 0.5133589 0.6329321 0.8041208
#acting        0.5042922 0.12600647 0.4957078 0.8432321 0.3758051 0.7391478 0.8537258
#general       0.3318778 0.33187776 0.3318778 0.8689634 0.6181534 0.8515257 0.8956363





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 8. India

library(haven)
India <- read_sav("India.sav")
colnames(India)
mydata <- as.data.frame(India[,8:46])
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .78, omega T = .83

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 5.71; eigenvalue 2 = 3.35

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 factors and 5 components
# Eigenvalue 1 = 6.47; eigenvalue 2 = 3.9

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.15, RMSEA=.083, RMSR=.11, TLI=.434

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.17, RMSEA=.105, RMSR=.13, TLI=.364

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.34, RMSEA=.04, RMSR=.04, TLI=.87
#       MR5   MR4   MR2  MR1   MR3
#MR5 1.00  0.17  0.13 0.43  0.25
#MR4 0.17  1.00 -0.25 0.28  0.24
#MR2 0.13 -0.25  1.00 0.05 -0.07
#MR1 0.43  0.28  0.05 1.00  0.23
#MR3 0.25  0.24 -0.07 0.23  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.065, RMSR=.05, TLI=.751
#     MR5   MR4   MR2   MR3  MR1
#MR5 1.00  0.15  0.16  0.26 0.43
#MR4 0.15  1.00 -0.26  0.24 0.27
#MR2 0.16 -0.26  1.00 -0.05 0.05
#MR3 0.26  0.24 -0.05  1.00 0.24
#MR1 0.43  0.27  0.05  0.24 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.696       0.528
#Tucker-Lewis Index (TLI)                       0.679       0.501
#Robust Comparative Fit Index (CFI)                         0.414
#Robust Tucker-Lewis Index (TLI)                            0.382
#RMSEA                                          0.134       0.097
#Robust RMSEA                                               0.107
#SRMR                                           0.127       0.127

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .162

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.459       0.459
#Tucker-Lewis Index (TLI)                       0.429       0.429
#Robust Comparative Fit Index (CFI)                         0.466
#Robust Tucker-Lewis Index (TLI)                            0.436
#RMSEA                                          0.087       0.082
#Robust RMSEA                                               0.085
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .103

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.802
#Tucker-Lewis Index (TLI)                       0.875       0.788
#Robust Comparative Fit Index (CFI)                         0.679
#Robust Tucker-Lewis Index (TLI)                            0.656
#RMSEA                                          0.083       0.063
#Robust RMSEA                                               0.080
#SRMR                                           0.090       0.090

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .305

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.310    0.059    5.232    0.000    0.310    0.310
#    nonjudging       -0.317    0.072   -4.428    0.000   -0.317   -0.317
#    nonreactivity     0.511    0.068    7.490    0.000    0.511    0.511
#    acting            0.057    0.068    0.838    0.402    0.057    0.057
#  describing ~~                                                         
#    nonjudging        0.340    0.058    5.880    0.000    0.340    0.340
#    nonreactivity     0.399    0.060    6.692    0.000    0.399    0.399
#    acting            0.636    0.040   16.044    0.000    0.636    0.636
#  nonjudging ~~                                                         
#    nonreactivity    -0.101    0.074   -1.358    0.174   -0.101   -0.101
#    acting            0.622    0.043   14.591    0.000    0.622    0.622
#  nonreactivity ~~                                                      
#    acting            0.129    0.066    1.971    0.049    0.129    0.129

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.746       0.757
#Tucker-Lewis Index (TLI)                       0.728       0.740
#Robust Comparative Fit Index (CFI)                         0.760
#Robust Tucker-Lewis Index (TLI)                            0.743
#RMSEA                                          0.060       0.055
#Robust RMSEA                                               0.057
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .273

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.216    0.090    2.402    0.016    0.216    0.216
#    nonjudging       -0.295    0.093   -3.161    0.002   -0.295   -0.295
#    nonreactivity     0.481    0.093    5.157    0.000    0.481    0.481
#    acting            0.004    0.091    0.048    0.962    0.004    0.004
#  describing ~~                                                         
#    nonjudging        0.379    0.078    4.867    0.000    0.379    0.379
#    nonreactivity     0.262    0.100    2.625    0.009    0.262    0.262
#    acting            0.678    0.055   12.256    0.000    0.678    0.678
#  nonjudging ~~                                                         
#    nonreactivity    -0.097    0.106   -0.913    0.361   -0.097   -0.097
#    acting            0.624    0.064    9.773    0.000    0.624    0.624
#  nonreactivity ~~                                                      
#    acting            0.125    0.111    1.121    0.262    0.125    0.125

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.656       0.666
#Tucker-Lewis Index (TLI)                       0.637       0.648
#Robust Comparative Fit Index (CFI)                         0.669
#Robust Tucker-Lewis Index (TLI)                            0.651
#RMSEA                                          0.069       0.064
#Robust RMSEA                                               0.067
#SRMR                                           0.131       0.131

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .267


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.804       0.816
#Tucker-Lewis Index (TLI)                       0.781       0.794
#Robust Comparative Fit Index (CFI)                         0.818
#Robust Tucker-Lewis Index (TLI)                            0.797
#RMSEA                                          0.054       0.049
#Robust RMSEA                                               0.051
#SRMR                                           0.081       0.081

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .301

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.4080620      0.8205128      0.8579911      0.5842421 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#observing     0.9835213 0.14957191 0.01647874 0.7119648 0.711794256 0.7190487 0.8484420
#describing    0.6024639 0.14067767 0.39753610 0.8121118 0.519878390 0.7082815 0.8442849
#nonjudging    0.6361095 0.12978499 0.36389055 0.7805460 0.511684541 0.6863589 0.8337152
#nonreactivity 0.9287170 0.10654347 0.07128296 0.6215369 0.619433445 0.6398355 0.8018017
#acting        0.2210588 0.06535994 0.77894118 0.8446757 0.002980213 0.4971241 0.8040992
#general       0.4080620 0.40806201 0.40806201 0.8579911 0.584242118 0.8870062 0.9426206





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 9. New Zealand

library(haven)
New.Zealand <- read.csv("New Zealand")
colnames(New.Zealand)
mydata <- as.data.frame(New.Zealand[,1:39])
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .90, omega T = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 8.39; eigenvalue 2 = 3.79

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 9.28; eigenvalue 2 = 4.17

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.127, RMSR=.15, TLI=.351

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.24, RMSEA=.149, RMSR=.16, TLI=.308

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.041, RMSR=.03, TLI=.931
#      MR2  MR1  MR4  MR3   MR5
#MR2  1.00 0.27 0.28 0.27 -0.07
#MR1  0.27 1.00 0.40 0.25  0.20
#MR4  0.28 0.40 1.00 0.18  0.06
#MR3  0.27 0.25 0.18 1.00  0.10
#MR5 -0.07 0.20 0.06 0.10  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.56, RMSEA=.066, RMSR=.03, TLI=.865
#      MR2  MR1  MR4  MR3   MR5
#MR2  1.00 0.27 0.28 0.27 -0.07
#MR1  0.27 1.00 0.41 0.25  0.20
#MR4  0.28 0.41 1.00 0.18  0.06
#MR3  0.27 0.25 0.18 1.00  0.09
#MR5 -0.07 0.20 0.06 0.09  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.442       0.353
#Tucker-Lewis Index (TLI)                       0.411       0.317
#Robust Comparative Fit Index (CFI)                         0.065
#Robust Tucker-Lewis Index (TLI)                            0.013
#RMSEA                                          0.360       0.200
#Robust RMSEA                                               0.183
#SRMR                                           0.283       0.283

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .182

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.387       0.351
#Tucker-Lewis Index (TLI)                       0.353       0.315
#Robust Comparative Fit Index (CFI)                         0.387
#Robust Tucker-Lewis Index (TLI)                            0.353
#RMSEA                                          0.129       0.122
#Robust RMSEA                                               0.128
#SRMR                                           0.151       0.151

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .156

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.983       0.952
#Tucker-Lewis Index (TLI)                       0.982       0.949
#Robust Comparative Fit Index (CFI)                         0.864
#Robust Tucker-Lewis Index (TLI)                            0.855
#RMSEA                                          0.063       0.055
#Robust RMSEA                                               0.070
#SRMR                                           0.066       0.066

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .566

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.293    0.045    6.458    0.000    0.293    0.293
#    nonjudging       -0.038    0.049   -0.785    0.433   -0.038   -0.038
#    nonreactivity     0.134    0.050    2.709    0.007    0.134    0.134
#    acting            0.135    0.051    2.666    0.008    0.135    0.135
#  describing ~~                                                         
#    nonjudging        0.310    0.044    6.965    0.000    0.310    0.310
#    nonreactivity     0.316    0.045    7.071    0.000    0.316    0.316
#    acting            0.447    0.042   10.637    0.000    0.447    0.447
#  nonjudging ~~                                                         
#    nonreactivity     0.349    0.043    8.092    0.000    0.349    0.349
#    acting            0.339    0.042    8.017    0.000    0.339    0.339
#  nonreactivity ~~                                                      
#    acting            0.253    0.045    5.598    0.000    0.253    0.25

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.909       0.919
#Tucker-Lewis Index (TLI)                       0.903       0.913
#Robust Comparative Fit Index (CFI)                         0.921
#Robust Tucker-Lewis Index (TLI)                            0.915
#RMSEA                                          0.050       0.043
#Robust RMSEA                                               0.046
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .489

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.266    0.063    4.230    0.000    0.266    0.266
#    nonjudging       -0.031    0.064   -0.482    0.630   -0.031   -0.031
#    nonreactivity     0.124    0.067    1.846    0.065    0.124    0.124
#    acting            0.111    0.065    1.702    0.089    0.111    0.111
#  describing ~~                                                         
#    nonjudging        0.298    0.054    5.521    0.000    0.298    0.298
#    nonreactivity     0.310    0.059    5.232    0.000    0.310    0.310
#    acting            0.440    0.053    8.283    0.000    0.440    0.440
#  nonjudging ~~                                                         
#    nonreactivity     0.349    0.056    6.205    0.000    0.349    0.349
#    acting            0.317    0.054    5.907    0.000    0.317    0.317
#  nonreactivity ~~                                                      
#    acting            0.244    0.061    4.007    0.000    0.244    0.244

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.884       0.893
#Tucker-Lewis Index (TLI)                       0.878       0.887
#Robust Comparative Fit Index (CFI)                         0.896
#Robust Tucker-Lewis Index (TLI)                            0.890
#RMSEA                                          0.056       0.049
#Robust RMSEA                                               0.053
#SRMR                                           0.143       0.143

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .484


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.927       0.933
#Tucker-Lewis Index (TLI)                       0.918       0.925
#Robust Comparative Fit Index (CFI)                         0.936
#Robust Tucker-Lewis Index (TLI)                            0.928
#RMSEA                                          0.046       0.040
#Robust RMSEA                                               0.043
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .490


library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3005459      0.8205128      0.9370352      0.5988137 
#
#$FactorLevelIndices
#                 ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#observing     0.9691818 0.1327240 0.03081818 0.7949121 0.7815851 0.8097444 0.8997999
#describing    0.6464924 0.1547840 0.35350755 0.9175838 0.5968037 0.8445475 0.8974021
#nonjudging    0.6784203 0.1752521 0.32157967 0.9325567 0.6374475 0.8610137 0.9028829
#nonreactivity 0.8307508 0.1327625 0.16924925 0.8478881 0.7095609 0.8161426 0.8972418
#acting        0.5057488 0.1039316 0.49425124 0.8828095 0.4212263 0.7940665 0.8731721
#general       0.3005459 0.3005459 0.30054585 0.9370352 0.5988137 0.8867508 0.8724313





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 10. Norway

library(haven)
Norway <- read_sav("Norway.sav")
colnames(Norway)
mydata <- as.data.frame(Norway[,5:43])
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .89, omega T = .92

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 7.95; eigenvalue 2 = 5.22

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 9.03; eigenvalue 2 = 5.82

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.20, RMSEA=.121, RMSR=.16, TLI=.360

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.151, RMSR=.18, TLI=.305

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.051, RMSR=.04, TLI=.882
#      MR1   MR3   MR4  MR2   MR5
#MR1  1.00 -0.19  0.52 0.17 -0.13
#MR3 -0.19  1.00 -0.21 0.22  0.28
#MR4  0.52 -0.21  1.00 0.09 -0.11
#MR2  0.17  0.22  0.09 1.00  0.30
#MR5 -0.13  0.28 -0.11 0.30  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.55, RMSEA=.091, RMSR=.05, TLI=.742
#      MR1   MR3   MR4  MR2   MR5
#MR1  1.00 -0.20  0.52 0.17 -0.12
#MR3 -0.20  1.00 -0.22 0.22  0.26
#MR4  0.52 -0.22  1.00 0.10 -0.10
#MR2  0.17  0.22  0.10 1.00  0.31
#MR5 -0.12  0.26 -0.10 0.31  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.709       0.543
#Tucker-Lewis Index (TLI)                       0.693       0.518
#Robust Comparative Fit Index (CFI)                         0.329
#Robust Tucker-Lewis Index (TLI)                            0.292
#RMSEA                                          0.234       0.149
#Robust RMSEA                                               0.161
#SRMR                                           0.193       0.193

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .289

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.399       0.389
#Tucker-Lewis Index (TLI)                       0.366       0.355
#Robust Comparative Fit Index (CFI)                         0.403
#Robust Tucker-Lewis Index (TLI)                            0.370
#RMSEA                                          0.126       0.118
#Robust RMSEA                                               0.124
#SRMR                                           0.163       0.163

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .129

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.929       0.847
#Tucker-Lewis Index (TLI)                       0.924       0.836
#Robust Comparative Fit Index (CFI)                         0.728
#Robust Tucker-Lewis Index (TLI)                            0.709
#RMSEA                                          0.117       0.087
#Robust RMSEA                                               0.103
#SRMR                                           0.115       0.115

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .515

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.418    0.051    8.229    0.000    0.418    0.418
#    nonjudging        0.140    0.065    2.159    0.031    0.140    0.140
#    nonreactivity     0.438    0.063    6.973    0.000    0.438    0.438
#    acting            0.060    0.067    0.894    0.371    0.060    0.060
#  describing ~~                                                         
#    nonjudging       -0.278    0.057   -4.899    0.000   -0.278   -0.278
#    nonreactivity     0.421    0.057    7.420    0.000    0.421    0.421
#    acting           -0.307    0.055   -5.609    0.000   -0.307   -0.307
#  nonjudging ~~                                                         
#    nonreactivity    -0.243    0.060   -4.060    0.000   -0.243   -0.243
#    acting            0.636    0.039   16.395    0.000    0.636    0.636
#  nonreactivity ~~                                                      
#    acting           -0.199    0.066   -3.034    0.002   -0.199   -0.199

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.809       0.814
#Tucker-Lewis Index (TLI)                       0.796       0.801
#Robust Comparative Fit Index (CFI)                         0.820
#Robust Tucker-Lewis Index (TLI)                            0.807
#RMSEA                                          0.071       0.066
#Robust RMSEA                                               0.069
#SRMR                                           0.102       0.102

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .467

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.452    0.087    5.186    0.000    0.452    0.452
#    nonjudging        0.094    0.087    1.089    0.276    0.094    0.094
#    nonreactivity     0.407    0.092    4.419    0.000    0.407    0.407
#    acting            0.038    0.102    0.368    0.713    0.038    0.038
#  describing ~~                                                         
#    nonjudging       -0.192    0.084   -2.284    0.022   -0.192   -0.192
#    nonreactivity     0.437    0.083    5.236    0.000    0.437    0.437
#    acting           -0.230    0.088   -2.614    0.009   -0.230   -0.230
#  nonjudging ~~                                                         
#    nonreactivity    -0.251    0.084   -2.984    0.003   -0.251   -0.251
#    acting            0.625    0.049   12.757    0.000    0.625    0.625
#  nonreactivity ~~                                                      
#    acting           -0.219    0.097   -2.262    0.024   -0.219   -0.219

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.765       0.769
#Tucker-Lewis Index (TLI)                       0.752       0.756
#Robust Comparative Fit Index (CFI)                         0.776
#Robust Tucker-Lewis Index (TLI)                            0.764
#RMSEA                                          0.079       0.073
#Robust RMSEA                                               0.076
#SRMR                                           0.170       0.170

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .475


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.846       0.850
#Tucker-Lewis Index (TLI)                       0.828       0.832
#Robust Comparative Fit Index (CFI)                         0.856
#Robust Tucker-Lewis Index (TLI)                            0.839
#RMSEA                                          0.066       0.060
#Robust RMSEA                                               0.063
#SRMR                                           0.146       0.146

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .496

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing =~                                                          
#    EM_1              0.358    0.117    3.050    0.002    0.358    0.310
#    EM_6              0.485    0.147    3.286    0.001    0.485    0.383
#    EM_11            -0.073    0.145   -0.503    0.615   -0.073   -0.056
#    EM_15             0.451    0.130    3.472    0.001    0.451    0.414
#    EM_20             0.509    0.105    4.829    0.000    0.509    0.466
#    EM_26             0.501    0.127    3.958    0.000    0.501    0.443
#    EM_31             0.693    0.105    6.613    0.000    0.693    0.542
#    EM_36            -0.182    0.140   -1.300    0.194   -0.182   -0.175
#  describing =~                                                         
#    EM_2              0.834    0.075   11.148    0.000    0.834    0.723
#    EM_7              0.749    0.079    9.419    0.000    0.749    0.714
#    EM_12            -0.817    0.070  -11.732    0.000   -0.817   -0.764
#    EM_16            -0.766    0.076  -10.107    0.000   -0.766   -0.677
#    EM_22            -0.616    0.071   -8.621    0.000   -0.616   -0.647
#    EM_27             0.680    0.065   10.523    0.000    0.680    0.647
#    EM_32             0.694    0.092    7.581    0.000    0.694    0.575
#    EM_37             0.574    0.093    6.175    0.000    0.574    0.485
#  nonjudging =~                                                         
#    EM_3              0.899    0.062   14.520    0.000    0.899    0.777
#    EM_10             0.961    0.055   17.452    0.000    0.961    0.800
#    EM_14             0.906    0.068   13.270    0.000    0.906    0.790
#    EM_17             0.803    0.070   11.493    0.000    0.803    0.668
#    EM_25             1.015    0.058   17.564    0.000    1.015    0.819
#    EM_30             1.009    0.059   17.222    0.000    1.009    0.850
#    EM_35             0.730    0.071   10.339    0.000    0.730    0.652
#    EM_39             0.726    0.075    9.660    0.000    0.726    0.589
#  nonreactivity =~                                                      
#    EM_4              0.439    0.083    5.262    0.000    0.439    0.441
#    EM_9              0.489    0.101    4.867    0.000    0.489    0.460
#    EM_19             0.509    0.087    5.834    0.000    0.509    0.450
#    EM_21             0.258    0.100    2.589    0.010    0.258    0.240
#    EM_24             0.793    0.089    8.937    0.000    0.793    0.670
#    EM_29             0.783    0.074   10.560    0.000    0.783    0.711
#    EM_33             0.633    0.093    6.779    0.000    0.633    0.589
#  acting =~                                                             
#    EM_5              0.739    0.084    8.767    0.000    0.739    0.684
#    EM_8              0.734    0.070   10.485    0.000    0.734    0.711
#    EM_13             0.844    0.075   11.195    0.000    0.844    0.750
#    EM_18             0.782    0.076   10.230    0.000    0.782    0.719
#    EM_23             0.680    0.081    8.436    0.000    0.680    0.626
#    EM_28             0.584    0.071    8.245    0.000    0.584    0.623
#    EM_34             0.489    0.086    5.657    0.000    0.489    0.517
#    EM_38             0.560    0.085    6.595    0.000    0.560    0.566
#  general =~                                                            
#    EM_1              0.423    0.091    4.628    0.000    0.423    0.367
#    EM_2              0.425    0.110    3.853    0.000    0.425    0.369
#    EM_3              0.223    0.120    1.859    0.063    0.223    0.192
#    EM_4              0.324    0.093    3.471    0.001    0.324    0.326
#    EM_5              0.117    0.117    1.004    0.315    0.117    0.109
#    EM_6              0.553    0.106    5.233    0.000    0.553    0.437
#    EM_7              0.408    0.103    3.975    0.000    0.408    0.389
#    EM_8              0.120    0.107    1.124    0.261    0.120    0.116
#    EM_9              0.394    0.107    3.681    0.000    0.394    0.370
#    EM_10             0.194    0.130    1.498    0.134    0.194    0.162
#    EM_11             0.617    0.083    7.436    0.000    0.617    0.474
#    EM_12             0.204    0.118    1.726    0.084    0.204    0.191
#    EM_13             0.094    0.118    0.790    0.429    0.094    0.083
#    EM_14             0.068    0.138    0.496    0.620    0.068    0.060
#    EM_15             0.610    0.096    6.355    0.000    0.610    0.561
#    EM_16             0.132    0.120    1.103    0.270    0.132    0.117
#    EM_17             0.321    0.112    2.875    0.004    0.321    0.267
#    EM_18             0.171    0.129    1.330    0.184    0.171    0.157
#    EM_19             0.391    0.095    4.140    0.000    0.391    0.346
#    EM_20             0.554    0.094    5.891    0.000    0.554    0.506
#    EM_21             0.432    0.093    4.666    0.000    0.432    0.401
#    EM_22             0.075    0.103    0.734    0.463    0.075    0.079
#    EM_23             0.056    0.123    0.454    0.650    0.056    0.051
#    EM_24             0.254    0.113    2.250    0.024    0.254    0.215
#    EM_25             0.162    0.125    1.298    0.194    0.162    0.131
#    EM_26             0.624    0.098    6.348    0.000    0.624    0.552
#    EM_27             0.525    0.107    4.896    0.000    0.525    0.499
#    EM_28             0.012    0.102    0.121    0.904    0.012    0.013
#    EM_29             0.308    0.102    3.020    0.003    0.308    0.279
#    EM_30             0.056    0.130    0.431    0.666    0.056    0.047
#    EM_31             0.577    0.099    5.802    0.000    0.577    0.450
#    EM_32             0.629    0.110    5.739    0.000    0.629    0.521
#    EM_33             0.275    0.109    2.512    0.012    0.275    0.256
#    EM_34             0.027    0.098    0.270    0.787    0.027    0.028
#    EM_35             0.290    0.123    2.360    0.018    0.290    0.259
#    EM_36             0.754    0.069   10.876    0.000    0.754    0.726
#    EM_37             0.730    0.083    8.754    0.000    0.730    0.617
#    EM_38             0.107    0.086    1.234    0.217    0.107    0.108
#    EM_39             0.369    0.114    3.251    0.001    0.369    0.300

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.2369179      0.8205128      0.9112988      0.5496472 
#$FactorLevelIndices
#               ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.3487998 0.06078319 0.65120023 0.8243881 0.20280907 0.5929649 0.7900302
#describing    0.7362411 0.18305682 0.26375894 0.7300999 0.09205168 0.8688900 0.9381925
#nonjudging    0.9344260 0.23564481 0.06557403 0.9208704 0.87132454 0.9210124 0.9594290
#nonreactivity 0.7341191 0.10358870 0.26588089 0.8019827 0.58160882 0.7615661 0.8733808
#acting        0.9793104 0.18000863 0.02068957 0.8589027 0.84502052 0.8643901 0.9299654
#general       0.2369179 0.23691786 0.23691786 0.9112988 0.54964720 0.8615937 0.9226033





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 11. Poland

library(readxl)
Poland <- read_excel("Poland.xlsx")
colnames(Poland)
mydata <- as.data.frame(Poland[,1:39])
#mydata[mydata<=0] <- NA
mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .76, omega T = .80

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 7 components
# Eigenvalue 1 = 4.63; eigenvalue 2 = 2.30

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 9 factors and 7 components
# Eigenvalue 1 = 5.11; eigenvalue 2 = 2.56

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.12, RMSEA=.08, RMSR=.09, TLI=.385

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.13, RMSEA=.093, RMSR=.10, TLI=.349

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 7 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 6 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.29, RMSEA=.051, RMSR=.04, TLI=.751
#     MR1   MR4   MR2  MR3  MR5
#MR1 1.00  0.42  0.02 0.06 0.06
#MR4 0.42  1.00 -0.05 0.14 0.02
#MR2 0.02 -0.05  1.00 0.18 0.10
#MR3 0.06  0.14  0.18 1.00 0.14
#MR5 0.06  0.02  0.10 0.14 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.32, RMSEA=.064, RMSR=.05, TLI=.69
#     MR1   MR5   MR2  MR3  MR4
#MR1 1.00  0.42  0.03 0.06 0.06
#MR5 0.42  1.00 -0.05 0.15 0.02
#MR2 0.03 -0.05  1.00 0.18 0.11
#MR3 0.06  0.15  0.18 1.00 0.13
#MR4 0.06  0.02  0.11 0.13 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.650       0.459
#Tucker-Lewis Index (TLI)                       0.630       0.428
#Robust Comparative Fit Index (CFI)                         0.386
#Robust Tucker-Lewis Index (TLI)                            0.352
#RMSEA                                          0.115       0.093
#Robust RMSEA                                               0.094
#SRMR                                           0.102       0.102

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .039

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.416       0.419
#Tucker-Lewis Index (TLI)                       0.384       0.386
#Robust Comparative Fit Index (CFI)                         0.423
#Robust Tucker-Lewis Index (TLI)                            0.391
#RMSEA                                          0.081       0.075
#Robust RMSEA                                               0.080
#SRMR                                           0.091       0.091

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .027

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.685       0.529
#Tucker-Lewis Index (TLI)                       0.663       0.496
#Robust Comparative Fit Index (CFI)                         0.466
#Robust Tucker-Lewis Index (TLI)                            0.429
#RMSEA                                          0.110       0.087
#Robust RMSEA                                               0.089
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .203

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.098    0.066    1.487    0.137    0.098    0.098
#    nonjudging        0.943    0.046   20.392    0.000    0.943    0.943
#    nonreactivity     0.648    0.063   10.265    0.000    0.648    0.648
#    acting            0.531    0.052   10.240    0.000    0.531    0.531
#  describing ~~                                                         
#    nonjudging        0.079    0.049    1.621    0.105    0.079    0.079
#    nonreactivity     0.048    0.056    0.870    0.384    0.048    0.048
#    acting            0.116    0.047    2.476    0.013    0.116    0.116
#  nonjudging ~~                                                         
#    nonreactivity     0.895    0.031   29.297    0.000    0.895    0.895
#    acting            0.790    0.024   32.276    0.000    0.790    0.790
#  nonreactivity ~~                                                      
#    acting            0.820    0.032   25.384    0.000    0.820    0.820

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.485       0.474
#Tucker-Lewis Index (TLI)                       0.449       0.437
#Robust Comparative Fit Index (CFI)                         0.490
#Robust Tucker-Lewis Index (TLI)                            0.454
#RMSEA                                          0.077       0.071
#Robust RMSEA                                               0.076
#SRMR                                           0.097       0.097

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .152

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.476    0.125    3.798    0.000    0.476    0.476
#    nonjudging        0.111    0.324    0.341    0.733    0.111    0.111
#    nonreactivity     0.201    0.226    0.890    0.373    0.201    0.201
#    acting            0.096    0.152    0.630    0.529    0.096    0.096
#  describing ~~                                                         
#    nonjudging        0.009    0.077    0.113    0.910    0.009    0.009
#    nonreactivity    -0.029    0.078   -0.373    0.709   -0.029   -0.029
#    acting            0.048    0.066    0.722    0.470    0.048    0.048
#  nonjudging ~~                                                         
#    nonreactivity     0.871    0.062   14.018    0.000    0.871    0.871
#    acting            0.724    0.050   14.472    0.000    0.724    0.724
#  nonreactivity ~~                                                      
#    acting            0.793    0.066   12.001    0.000    0.793    0.793


library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.398       0.394
#Tucker-Lewis Index (TLI)                       0.365       0.361
#Robust Comparative Fit Index (CFI)                         0.404
#Robust Tucker-Lewis Index (TLI)                            0.370
#RMSEA                                          0.083       0.076
#Robust RMSEA                                               0.081
#SRMR                                           0.124       0.124

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .125


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.515       0.501
#Tucker-Lewis Index (TLI)                       0.458       0.442
#Robust Comparative Fit Index (CFI)                         0.519
#Robust Tucker-Lewis Index (TLI)                            0.462
#RMSEA                                          0.076       0.071
#Robust RMSEA                                               0.075
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .209

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.2638195      0.8205128      0.6780412      0.3323439 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#observing     0.7737994 0.09758040 0.22620064 0.3450817 0.262892775 0.6085237 0.7849765
#describing    0.3615743 0.09004829 0.63842569 0.6112301 0.002389308 0.4875785 0.7535540
#nonjudging    0.8487613 0.23176280 0.15123869 0.7233389 0.704254546 0.7681445 0.8765835
#nonreactivity 0.8748334 0.10964351 0.12516662 0.4875705 0.462981658 0.5763404 0.7638007
#acting        0.9147169 0.20714555 0.08528308 0.6322696 0.617758721 0.7705437 0.8781194
#general       0.2638195 0.26381947 0.26381947 0.6780412 0.332343859 0.7565748 0.8808934





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 12. Portugal

library(haven)
Portugal <- read_sav("Portugal.sav")
colnames(Portugal)
mydata <- as.data.frame(Portugal[,2:40])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .91, omega T = .93

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 8.18; eigenvalue 2 = 3.61

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 9.35; eigenvalue 2 = 3.97

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.21, RMSEA=.112, RMSR=.13, TLI=.392

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.24, RMSEA=.137, RMSR=.15, TLI=.336

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.46, RMSEA=.044, RMSR=.04, TLI=.904
#     MR1  MR4  MR3  MR2  MR5
#MR1 1.00 0.32 0.35 0.14 0.17
#MR4 0.32 1.00 0.19 0.29 0.29
#MR3 0.35 0.19 1.00 0.09 0.03
#MR2 0.14 0.29 0.09 1.00 0.26
#MR5 0.17 0.29 0.03 0.26 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.077, RMSR=.04, TLI=.787
#     MR1  MR4  MR2  MR3  MR5
#MR1 1.00 0.33 0.36 0.16 0.17
#MR4 0.33 1.00 0.20 0.30 0.29
#MR2 0.36 0.20 1.00 0.10 0.03
#MR3 0.16 0.30 0.10 1.00 0.26
#MR5 0.17 0.29 0.03 0.26 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.797       0.583
#Tucker-Lewis Index (TLI)                       0.786       0.559
#Robust Comparative Fit Index (CFI)                         0.376
#Robust Tucker-Lewis Index (TLI)                            0.341
#RMSEA                                          0.197       0.138
#Robust RMSEA                                               0.144
#SRMR                                           0.154       0.154

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .240

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.421       0.403
#Tucker-Lewis Index (TLI)                       0.389       0.370
#Robust Comparative Fit Index (CFI)                         0.425
#Robust Tucker-Lewis Index (TLI)                            0.393
#RMSEA                                          0.116       0.108
#Robust RMSEA                                               0.114
#SRMR                                           0.130       0.130

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .171

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.979       0.939
#Tucker-Lewis Index (TLI)                       0.978       0.934
#Robust Comparative Fit Index (CFI)                         0.797
#Robust Tucker-Lewis Index (TLI)                            0.783
#RMSEA                                          0.063       0.053
#Robust RMSEA                                               0.083
#SRMR                                           0.071       0.071

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .511

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.475    0.049    9.651    0.000    0.475    0.475
#    nonjudging        0.195    0.061    3.174    0.002    0.195    0.195
#    nonreactivity     0.603    0.046   13.169    0.000    0.603    0.603
#    acting            0.276    0.055    5.038    0.000    0.276    0.276
#  describing ~~                                                         
#    nonjudging        0.269    0.054    5.003    0.000    0.269    0.269
#    nonreactivity     0.515    0.047   10.989    0.000    0.515    0.515
#    acting            0.388    0.049    7.945    0.000    0.388    0.388
#  nonjudging ~~                                                         
#    nonreactivity     0.147    0.059    2.477    0.013    0.147    0.147
#    acting            0.445    0.048    9.237    0.000    0.445    0.445
#  nonreactivity ~~                                                      
#    acting            0.346    0.055    6.255    0.000    0.346    0.346

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.874       0.888
#Tucker-Lewis Index (TLI)                       0.865       0.880
#Robust Comparative Fit Index (CFI)                         0.892
#Robust Tucker-Lewis Index (TLI)                            0.884
#RMSEA                                          0.055       0.047
#Robust RMSEA                                               0.050
#SRMR                                           0.065       0.065

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.430    0.071    6.078    0.000    0.430    0.430
#    nonjudging        0.175    0.080    2.171    0.030    0.175    0.175
#    nonreactivity     0.575    0.069    8.347    0.000    0.575    0.575
#    acting            0.240    0.077    3.119    0.002    0.240    0.240
#  describing ~~                                                         
#    nonjudging        0.238    0.080    2.992    0.003    0.238    0.238
#    nonreactivity     0.520    0.064    8.118    0.000    0.520    0.520
#    acting            0.372    0.073    5.106    0.000    0.372    0.372
#  nonjudging ~~                                                         
#    nonreactivity     0.138    0.098    1.405    0.160    0.138    0.138
#    acting            0.417    0.066    6.345    0.000    0.417    0.417
#  nonreactivity ~~                                                      
#    acting            0.308    0.086    3.574    0.000    0.308    0.308

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.830       0.841
#Tucker-Lewis Index (TLI)                       0.820       0.832
#Robust Comparative Fit Index (CFI)                         0.847
#Robust Tucker-Lewis Index (TLI)                            0.838
#RMSEA                                          0.063       0.056
#Robust RMSEA                                               0.059
#SRMR                                           0.153       0.153

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .434


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.902       0.914
#Tucker-Lewis Index (TLI)                       0.891       0.904
#Robust Comparative Fit Index (CFI)                         0.918
#Robust Tucker-Lewis Index (TLI)                            0.908
#RMSEA                                          0.049       0.042
#Robust RMSEA                                               0.044
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .453

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3424645      0.8205128      0.9266706      0.6254258 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega     OmegaH         H        FD
#observing     0.7273020 0.12142322 0.27269798 0.8231284 0.60739709 0.7582799 0.8746823
#describing    0.2158658 0.05913116 0.78413423 0.9158507 0.08892475 0.6063834 0.8412444
#nonjudging    0.9451910 0.19513371 0.05480899 0.8703078 0.83331917 0.8733967 0.9359742
#nonreactivity 0.6645139 0.07127991 0.33548608 0.7196454 0.47688688 0.6206096 0.7961392
#acting        0.8580269 0.21056750 0.14197307 0.9063300 0.77718405 0.8904472 0.9467919
#general       0.3424645 0.34246449 0.34246449 0.9266706 0.62542577 0.9200413 0.9513992





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 13. Romania

### not included

library(readxl)
Romania <- read_excel("Romania.xlsx")
colnames(Romania)
mydata <- as.data.frame(Romania[,7:45])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .87, omega T = .91

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 7.63; eigenvalue 2 = 4.2

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# error

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.20, RMSEA=.102, RMSR=.14, TLI=.389

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# error

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# warning
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.45, RMSEA=.051, RMSR=.05, TLI=.839
#      MR2   MR1   MR4  MR5   MR3
#MR2  1.00 -0.31 -0.11 0.05  0.11
#MR1 -0.31  1.00  0.34 0.09 -0.13
#MR4 -0.11  0.34  1.00 0.34 -0.01
#MR5  0.05  0.09  0.34 1.00  0.18
#MR3  0.11 -0.13 -0.01 0.18  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# NA

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.397       0.368
#Tucker-Lewis Index (TLI)                       0.363       0.333
#Robust Comparative Fit Index (CFI)                         0.391
#Robust Tucker-Lewis Index (TLI)                            0.357
#RMSEA                                          0.117       0.115
#Robust RMSEA                                               0.116
#SRMR                                           0.141       0.141

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .146

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.719       0.710
#Tucker-Lewis Index (TLI)                       0.699       0.690
#Robust Comparative Fit Index (CFI)                         0.721
#Robust Tucker-Lewis Index (TLI)                            0.701
#RMSEA                                          0.081       0.078
#Robust RMSEA                                               0.079
#SRMR                                           0.120       0.120

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .391

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.171    0.155    1.104    0.270    0.171    0.171
#    nonjudging        0.406    0.104    3.892    0.000    0.406    0.406
#    nonreactivity     0.245    0.166    1.478    0.140    0.245    0.245
#    acting            0.092    0.163    0.565    0.572    0.092    0.092
#  describing ~~                                                         
#    nonjudging       -0.238    0.102   -2.340    0.019   -0.238   -0.238
#    nonreactivity     0.305    0.175    1.745    0.081    0.305    0.305
#    acting           -0.365    0.123   -2.977    0.003   -0.365   -0.365
#  nonjudging ~~                                                         
#    nonreactivity    -0.126    0.187   -0.677    0.498   -0.126   -0.126
#    acting            0.566    0.128    4.416    0.000    0.566    0.566
#  nonreactivity ~~                                                      
#    acting           -0.356    0.179   -1.986    0.047   -0.356   -0.356

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE,orthogonal=TRUE)
# warning
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.681       0.676
#Tucker-Lewis Index (TLI)                       0.663       0.658
#Robust Comparative Fit Index (CFI)                         0.685
#Robust Tucker-Lewis Index (TLI)                            0.668
#RMSEA                                          0.085       0.082
#Robust RMSEA                                               0.084
#SRMR                                           0.169       0.169

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .372


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
# warning
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 14. Spain

library(haven)
Spain <- read_sav("Spain.sav")
colnames(Spain)
mydata <- as.data.frame(Spain[,3:41])
#mydata[mydata<=0] <- NA
mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .91, omega T = .94

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 8.70; eigenvalue 2 = 4.15

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 9.72; eigenvalue 2 = 4.65

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.129, RMSR=.15, TLI=.358

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.15, RMSR=.17, TLI=.323

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.049, RMSR=.03, TLI=.906
#      MR1   MR3   MR4   MR2   MR5
#MR1  1.00 -0.26  0.43 -0.02 -0.33
#MR3 -0.26  1.00 -0.31  0.19  0.30
#MR4  0.43 -0.31  1.00  0.04 -0.21
#MR2 -0.02  0.19  0.04  1.00  0.45
#MR5 -0.33  0.30 -0.21  0.45  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.066, RMSR=.03, TLI=.868
#      MR1   MR3   MR4   MR2   MR5
#MR1  1.00 -0.26  0.44 -0.02 -0.33
#MR3 -0.26  1.00 -0.31  0.20  0.30
#MR4  0.44 -0.31  1.00  0.04 -0.21
#MR2 -0.02  0.20  0.04  1.00  0.46
#MR5 -0.33  0.30 -0.21  0.46  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.756       0.537
#Tucker-Lewis Index (TLI)                       0.743       0.511
#Robust Comparative Fit Index (CFI)                         0.347
#Robust Tucker-Lewis Index (TLI)                            0.310
#RMSEA                                          0.227       0.165
#Robust RMSEA                                               0.153
#SRMR                                           0.176       0.176

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .366

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.395       0.385
#Tucker-Lewis Index (TLI)                       0.362       0.351
#Robust Comparative Fit Index (CFI)                         0.397
#Robust Tucker-Lewis Index (TLI)                            0.364
#RMSEA                                          0.130       0.116
#Robust RMSEA                                               0.129
#SRMR                                           0.148       0.148

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .210

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.974       0.929
#Tucker-Lewis Index (TLI)                       0.973       0.924
#Robust Comparative Fit Index (CFI)                         0.856
#Robust Tucker-Lewis Index (TLI)                            0.845
#RMSEA                                          0.074       0.065
#Robust RMSEA                                               0.072
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .598

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.293    0.026   11.441    0.000    0.293    0.293
#    nonjudging       -0.091    0.030   -3.080    0.002   -0.091   -0.091
#    nonreactivity     0.570    0.022   26.120    0.000    0.570    0.570
#    acting           -0.039    0.029   -1.343    0.179   -0.039   -0.039
#  describing ~~                                                         
#    nonjudging       -0.326    0.026  -12.513    0.000   -0.326   -0.326
#    nonreactivity     0.370    0.025   14.726    0.000    0.370    0.370
#    acting           -0.360    0.025  -14.266    0.000   -0.360   -0.360
#  nonjudging ~~                                                         
#    nonreactivity    -0.414    0.025  -16.413    0.000   -0.414   -0.414
#    acting            0.455    0.023   19.923    0.000    0.455    0.455
#  nonreactivity ~~                                                      
#    acting           -0.278    0.026  -10.523    0.000   -0.278   -0.278

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.886       0.889
#Tucker-Lewis Index (TLI)                       0.878       0.881
#Robust Comparative Fit Index (CFI)                         0.892
#Robust Tucker-Lewis Index (TLI)                            0.884
#RMSEA                                          0.057       0.050
#Robust RMSEA                                               0.055
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .517

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.261    0.033    7.921    0.000    0.261    0.261
#    nonjudging       -0.095    0.034   -2.782    0.005   -0.095   -0.095
#    nonreactivity     0.556    0.029   19.267    0.000    0.556    0.556
#    acting           -0.034    0.037   -0.925    0.355   -0.034   -0.034
#  describing ~~                                                         
#    nonjudging       -0.308    0.035   -8.832    0.000   -0.308   -0.308
#    nonreactivity     0.364    0.032   11.216    0.000    0.364    0.364
#    acting           -0.358    0.033  -10.735    0.000   -0.358   -0.358
#  nonjudging ~~                                                         
#    nonreactivity    -0.403    0.035  -11.404    0.000   -0.403   -0.403
#    acting            0.459    0.030   15.129    0.000    0.459    0.459
#  nonreactivity ~~                                                      
#    acting           -0.279    0.038   -7.378    0.000   -0.279   -0.279

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.845       0.847
#Tucker-Lewis Index (TLI)                       0.836       0.839
#Robust Comparative Fit Index (CFI)                         0.851
#Robust Tucker-Lewis Index (TLI)                            0.842
#RMSEA                                          0.066       0.058
#Robust RMSEA                                               0.064
#SRMR                                           0.161       0.161

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .518


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.888       0.889
#Tucker-Lewis Index (TLI)                       0.875       0.876
#Robust Comparative Fit Index (CFI)                         0.893
#Robust Tucker-Lewis Index (TLI)                            0.881
#RMSEA                                          0.057       0.051
#Robust RMSEA                                               0.056
#SRMR                                           0.096       0.096

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .521

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.2846515      0.8205128      0.8431283      0.1796541 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#observing     0.6111316 0.10215615 0.38886837 0.8443794 0.5080431 0.7440930 0.8527504
#describing    0.7977039 0.17767516 0.20229610 0.4341918 0.1939365 0.8702489 0.9273512
#nonjudging    0.8609433 0.20611823 0.13905670 0.9196293 0.7920225 0.8974403 0.9435858
#nonreactivity 0.2273190 0.03643373 0.77268102 0.8353787 0.1041090 0.4945081 0.7296037
#acting        0.9170385 0.19296528 0.08296149 0.8946345 0.8267232 0.8827952 0.9377817
#general       0.2846515 0.28465145 0.28465145 0.8431283 0.1796541 0.8814359 0.9006682





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 15. Sweden

library(haven)
sweden <- read_sav("sweden.sav")
colnames(sweden)
mydata <- as.data.frame(sweden[,6:44])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .86, omega T = .89

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 6.87; eigenvalue 2 = 4.13

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 7.71; eigenvalue 2 = 4.59

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.18, RMSEA=.113, RMSR=.14, TLI=.346

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.20, RMSEA=.111, RMSR=.15, TLI=.31

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.049, RMSR=.03, TLI=.876
#      MR1  MR2   MR4   MR5  MR3
#MR1  1.00 0.18  0.49 -0.13 0.13
#MR2  0.18 1.00  0.27  0.21 0.21
#MR4  0.49 0.27  1.00 -0.05 0.02
#MR5 -0.13 0.21 -0.05  1.00 0.26
#MR3  0.13 0.21  0.02  0.26 1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.071, RMSR=.04, TLI=.800
#      MR1  MR2   MR4   MR5  MR3
#MR1  1.00 0.18  0.49 -0.14 0.13
#MR2  0.18 1.00  0.27  0.21 0.21
#MR4  0.49 0.27  1.00 -0.05 0.02
#MR5 -0.14 0.21 -0.05  1.00 0.25
#MR3  0.13 0.21  0.02  0.25 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.711       0.572
#Tucker-Lewis Index (TLI)                       0.695       0.549
#Robust Comparative Fit Index (CFI)                         0.344
#Robust Tucker-Lewis Index (TLI)                            0.307
#RMSEA                                          0.190       0.130
#Robust RMSEA                                               0.137
#SRMR                                           0.159       0.159

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .291

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.385       0.375
#Tucker-Lewis Index (TLI)                       0.351       0.340
#Robust Comparative Fit Index (CFI)                         0.389
#Robust Tucker-Lewis Index (TLI)                            0.355
#RMSEA                                          0.115       0.105
#Robust RMSEA                                               0.113
#SRMR                                           0.139       0.139

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .151

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.948       0.893
#Tucker-Lewis Index (TLI)                       0.945       0.885
#Robust Comparative Fit Index (CFI)                         0.786
#Robust Tucker-Lewis Index (TLI)                            0.771
#RMSEA                                          0.081       0.066
#Robust RMSEA                                               0.079
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .508

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.357    0.040    8.984    0.000    0.357    0.357
#    nonjudging       -0.130    0.047   -2.771    0.006   -0.130   -0.130
#    nonreactivity     0.417    0.043    9.752    0.000    0.417    0.417
#    acting            0.016    0.046    0.345    0.730    0.016    0.016
#  describing ~~                                                         
#    nonjudging        0.248    0.043    5.811    0.000    0.248    0.248
#    nonreactivity     0.300    0.044    6.778    0.000    0.300    0.300
#    acting            0.359    0.040    9.035    0.000    0.359    0.359
#  nonjudging ~~                                                         
#    nonreactivity     0.223    0.050    4.451    0.000    0.223    0.223
#    acting            0.547    0.035   15.800    0.000    0.547    0.547
#  nonreactivity ~~                                                      
#    acting            0.091    0.046    2.002    0.045    0.091    0.091

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.837       0.846
#Tucker-Lewis Index (TLI)                       0.826       0.835
#Robust Comparative Fit Index (CFI)                         0.850
#Robust Tucker-Lewis Index (TLI)                            0.839
#RMSEA                                          0.059       0.052
#Robust RMSEA                                               0.057
#SRMR                                           0.072       0.072

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.302    0.061    4.969    0.000    0.302    0.302
#    nonjudging       -0.094    0.057   -1.634    0.102   -0.094   -0.094
#    nonreactivity     0.406    0.060    6.771    0.000    0.406    0.406
#    acting            0.011    0.059    0.186    0.853    0.011    0.011
#  describing ~~                                                         
#    nonjudging        0.228    0.055    4.113    0.000    0.228    0.228
#    nonreactivity     0.274    0.063    4.335    0.000    0.274    0.274
#    acting            0.334    0.059    5.669    0.000    0.334    0.334
#  nonjudging ~~                                                         
#    nonreactivity     0.190    0.072    2.651    0.008    0.190    0.190
#    acting            0.550    0.043   12.827    0.000    0.550    0.550
#  nonreactivity ~~                                                      
#    acting            0.065    0.073    0.895    0.371    0.065    0.065

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.799       0.807
#Tucker-Lewis Index (TLI)                       0.787       0.796
#Robust Comparative Fit Index (CFI)                         0.811
#Robust Tucker-Lewis Index (TLI)                            0.800
#RMSEA                                          0.066       0.058
#Robust RMSEA                                               0.063
#SRMR                                           0.131       0.131

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .426


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, ordered=TRUE,orthogonal=TRUE)
#warning
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.859       0.866
#Tucker-Lewis Index (TLI)                       0.843       0.851
#Robust Comparative Fit Index (CFI)                         0.871
#Robust Tucker-Lewis Index (TLI)                            0.855
#RMSEA                                          0.056       0.050
#Robust RMSEA                                               0.054
#SRMR                                           0.113       0.113

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .438

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing =~                                                          
#    EM_1              0.449    0.076    5.893    0.000    0.449    0.407
#    EM_6              0.731    0.069   10.593    0.000    0.731    0.585
#    EM_11             0.413    0.085    4.850    0.000    0.413    0.348
#    EM_15             0.685    0.055   12.510    0.000    0.685    0.651
#    EM_20             0.351    0.071    4.968    0.000    0.351    0.344
#    EM_26             0.288    0.053    5.435    0.000    0.288    0.347
#    EM_31             0.375    0.067    5.600    0.000    0.375    0.335
#    EM_36             0.036    0.074    0.491    0.623    0.036    0.038
#  describing =~                                                         
#    EM_2              0.516    0.080    6.439    0.000    0.516    0.541
#    EM_7              0.526    0.057    9.224    0.000    0.526    0.584
#    EM_12             0.779    0.041   19.087    0.000    0.779    0.807
#    EM_16             0.718    0.045   15.962    0.000    0.718    0.770
#    EM_22             0.551    0.048   11.526    0.000    0.551    0.635
#    EM_27             0.560    0.066    8.453    0.000    0.560    0.547
#    EM_32             0.391    0.110    3.571    0.000    0.391    0.382
#    EM_37             0.496    0.076    6.549    0.000    0.496    0.515
#  nonjudging =~                                                         
#    EM_3              0.699    0.045   15.563    0.000    0.699    0.671
#    EM_10             0.742    0.041   18.199    0.000    0.742    0.709
#    EM_14             0.825    0.040   20.595    0.000    0.825    0.798
#    EM_17             0.651    0.048   13.660    0.000    0.651    0.594
#    EM_25             0.835    0.036   23.231    0.000    0.835    0.790
#    EM_30             0.827    0.038   21.488    0.000    0.827    0.807
#    EM_35             0.671    0.041   16.254    0.000    0.671    0.648
#    EM_39             0.754    0.048   15.752    0.000    0.754    0.691
#  nonreactivity =~                                                      
#    EM_4              0.346    0.053    6.570    0.000    0.346    0.372
#    EM_9              0.441    0.055    7.956    0.000    0.441    0.487
#    EM_19             0.399    0.063    6.376    0.000    0.399    0.387
#    EM_21             0.279    0.056    5.004    0.000    0.279    0.307
#    EM_24             0.523    0.049   10.588    0.000    0.523    0.556
#    EM_29             0.604    0.046   13.032    0.000    0.604    0.636
#    EM_33             0.567    0.049   11.546    0.000    0.567    0.629
#  acting =~                                                             
#    EM_5              0.618    0.048   12.877    0.000    0.618    0.654
#    EM_8              0.661    0.042   15.896    0.000    0.661    0.676
#    EM_13             0.597    0.053   11.229    0.000    0.597    0.655
#    EM_18             0.640    0.040   16.183    0.000    0.640    0.700
#    EM_23             0.515    0.051   10.092    0.000    0.515    0.576
#    EM_28             0.522    0.041   12.589    0.000    0.522    0.630
#    EM_34             0.554    0.039   14.202    0.000    0.554    0.671
#    EM_38             0.543    0.044   12.394    0.000    0.543    0.611
#  general =~                                                            
#    EM_1              0.288    0.075    3.833    0.000    0.288    0.261
#    EM_2              0.538    0.084    6.435    0.000    0.538    0.564
#    EM_3             -0.007    0.086   -0.078    0.938   -0.007   -0.006
#    EM_4              0.280    0.054    5.236    0.000    0.280    0.301
#    EM_5              0.093    0.089    1.040    0.298    0.093    0.098
#    EM_6              0.310    0.088    3.519    0.000    0.310    0.248
#    EM_7              0.383    0.081    4.751    0.000    0.383    0.425
#    EM_8             -0.014    0.086   -0.162    0.871   -0.014   -0.014
#    EM_9              0.214    0.065    3.317    0.001    0.214    0.237
#    EM_10            -0.040    0.097   -0.412    0.680   -0.040   -0.038
#    EM_11             0.272    0.090    3.019    0.003    0.272    0.229
#    EM_12             0.208    0.109    1.910    0.056    0.208    0.215
#    EM_13             0.083    0.091    0.908    0.364    0.083    0.091
#    EM_14             0.035    0.098    0.354    0.723    0.035    0.034
#    EM_15             0.439    0.066    6.674    0.000    0.439    0.417
#    EM_16             0.120    0.108    1.108    0.268    0.120    0.129
#    EM_17            -0.238    0.094   -2.517    0.012   -0.238   -0.217
#    EM_18            -0.077    0.087   -0.892    0.373   -0.077   -0.085
#    EM_19             0.454    0.067    6.818    0.000    0.454    0.440
#    EM_20             0.271    0.077    3.533    0.000    0.271    0.266
#    EM_21             0.279    0.064    4.385    0.000    0.279    0.307
#    EM_22             0.180    0.096    1.883    0.060    0.180    0.207
#    EM_23             0.053    0.074    0.709    0.479    0.053    0.059
#    EM_24             0.132    0.061    2.153    0.031    0.132    0.140
#    EM_25            -0.078    0.095   -0.826    0.409   -0.078   -0.074
#    EM_26             0.290    0.054    5.369    0.000    0.290    0.349
#    EM_27             0.359    0.085    4.224    0.000    0.359    0.352
#    EM_28            -0.004    0.069   -0.059    0.953   -0.004   -0.005
#    EM_29             0.299    0.063    4.712    0.000    0.299    0.314
#    EM_30            -0.028    0.101   -0.274    0.784   -0.028   -0.027
#    EM_31             0.648    0.054   11.976    0.000    0.648    0.579
#    EM_32             0.726    0.080    9.071    0.000    0.726    0.710
#    EM_33             0.068    0.061    1.110    0.267    0.068    0.076
#    EM_34             0.002    0.065    0.033    0.974    0.002    0.003
#    EM_35            -0.042    0.095   -0.439    0.661   -0.042   -0.040
#    EM_36             0.484    0.068    7.150    0.000    0.484    0.512
#    EM_37             0.540    0.080    6.715    0.000    0.540    0.560
#    EM_38            -0.057    0.070   -0.813    0.416   -0.057   -0.064
#    EM_39             0.025    0.090    0.283    0.777    0.025    0.023

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#     0.1972237      0.8205128      0.8790290      0.3145632 
#$FactorLevelIndices
#               ECV_SS     ECV_SG      ECV_GS     Omega    OmegaH         H        FD
#observing     0.5509852 0.08296954 0.449014817 0.7627414 0.4064982 0.6658130 0.8082421
#describing    0.6594729 0.17670895 0.340527120 0.9047441 0.6296028 0.8550172 0.9085027
#nonjudging    0.9861034 0.24294216 0.013896643 0.8952432 0.8919657 0.9028638 0.9503826
#nonreactivity 0.7556751 0.10207105 0.244324877 0.7570264 0.5870916 0.7170273 0.8427749
#acting        0.9903103 0.19808456 0.009689655 0.8530181 0.8528017 0.8544171 0.9248387
#general       0.1972237 0.19722375 0.197223747 0.8790290 0.3145632 0.8220918 0.8771877





################################################################
### Five-Facet Mindfulness Questionnaire
### FFMQ
### cross cultural analysis reported Karl et al. (2020) https://link.springer.com/article/10.1007/s12671-020-01333-6
### data at https://osf.io/3yhp5 
###
### 16. US

library(readxl)
us1 <- read_excel("us1.xlsx")
colnames(us1)
mydata <- as.data.frame(us1[,3:41])
#mydata[mydata<=0] <- NA
#mydata[mydata>=6] <- NA
mydata <- na.omit(mydata)
colnames(mydata)
colnames(mydata) <- c("EM_1","EM_2","EM_3","EM_4","EM_5","EM_6","EM_7","EM_8","EM_9","EM_10","EM_11","EM_12","EM_13","EM_14","EM_15",
                      "EM_16","EM_17","EM_18","EM_19","EM_20","EM_21","EM_22","EM_23","EM_24","EM_25","EM_26","EM_27","EM_28","EM_29","EM_30",
                      "EM_31","EM_32","EM_33","EM_34","EM_35","EM_36","EM_37","EM_38","EM_39")						
colnames(mydata)
mydata[mydata=="Never or rarely true"] <- 1
mydata[mydata=="Rarely true"] <- 2
mydata[mydata=="Sometimes true"] <- 3
mydata[mydata=="Often true"] <- 4
mydata[mydata=="Very often or always true"] <- 5
mydata <- as.data.frame(sapply(mydata, as.numeric))
min(mydata)
max(mydata) # response alternatives = 5

library(psych)
omega(mydata) # alpha = .93, omega T = .95

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 10.65; eigenvalue 2 = 5.59

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 6 factors and 5 components
# Eigenvalue 1 = 12.07; eigenvalue 2 = 6.18

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.27, RMSEA=.145, RMSR=.18, TLI=.368

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.31, RMSEA=.177, RMSR=.19, TLI=.318

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias
EGAL <- EGA(mydata, algorithm="louvain")
EGAL
# suggests 5 communities


# Give solution with 5 factors (EFA-based)
fit4 <- fa(mydata,5)
fit4
diagram(fit4)
# %variance explained=.58, RMSEA=.068, RMSR=.04, TLI=.857
#      MR1   MR3   MR5   MR2   MR4
#MR1  1.00 -0.34  0.61  0.17 -0.17
#MR3 -0.34  1.00 -0.39  0.24  0.28
#MR5  0.61 -0.39  1.00 -0.05 -0.19
#MR2  0.17  0.24 -0.05  1.00  0.40
#MR4 -0.17  0.28 -0.19  0.40  1.00

fit4 <- fa(rho,5,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.65, RMSEA=.106, RMSR=.04, TLI=.752
#      MR1   MR3   MR5   MR2   MR4
#MR1  1.00 -0.35  0.60  0.18 -0.18
#MR3 -0.35  1.00 -0.39  0.24  0.29
#MR5  0.60 -0.39  1.00 -0.06 -0.20
#MR2  0.18  0.24 -0.06  1.00  0.39
#MR4 -0.18  0.29 -0.20  0.39  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.804       0.621
#Tucker-Lewis Index (TLI)                       0.794       0.600
#Robust Comparative Fit Index (CFI)                         0.358
#Robust Tucker-Lewis Index (TLI)                            0.322
#RMSEA                                          0.288       0.175
#Robust RMSEA                                               0.184
#SRMR                                           0.204       0.204

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .439

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.408       0.398
#Tucker-Lewis Index (TLI)                       0.376       0.364
#Robust Comparative Fit Index (CFI)                         0.415
#Robust Tucker-Lewis Index (TLI)                            0.382
#RMSEA                                          0.149       0.131
#Robust RMSEA                                               0.146
#SRMR                                           0.177       0.177

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .246

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 5 factors (theory based) 
EGAmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.977       0.922
#Tucker-Lewis Index (TLI)                       0.976       0.916
#Robust Comparative Fit Index (CFI)                         0.771
#Robust Tucker-Lewis Index (TLI)                            0.755
#RMSEA                                          0.099       0.080
#Robust RMSEA                                               0.111
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .641

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.327    0.053    6.217    0.000    0.327    0.327
#    nonjudging        0.120    0.052    2.306    0.021    0.120    0.120
#    nonreactivity     0.490    0.048   10.307    0.000    0.490    0.490
#    acting           -0.114    0.053   -2.138    0.033   -0.114   -0.114
#  describing ~~                                                         
#    nonjudging       -0.472    0.042  -11.280    0.000   -0.472   -0.472
#    nonreactivity     0.392    0.046    8.449    0.000    0.392    0.392
#    acting           -0.520    0.039  -13.488    0.000   -0.520   -0.520
#  nonjudging ~~                                                         
#    nonreactivity    -0.283    0.051   -5.504    0.000   -0.283   -0.283
#    acting            0.666    0.031   21.160    0.000    0.666    0.666
#  nonreactivity ~~                                                      
#    acting           -0.287    0.055   -5.249    0.000   -0.287   -0.287

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.840       0.852
#Tucker-Lewis Index (TLI)                       0.829       0.842
#Robust Comparative Fit Index (CFI)                         0.858
#Robust Tucker-Lewis Index (TLI)                            0.847
#RMSEA                                          0.078       0.065
#Robust RMSEA                                               0.072
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .546

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  observing ~~                                                          
#    describing        0.299    0.087    3.431    0.001    0.299    0.299
#    nonjudging        0.110    0.070    1.579    0.114    0.110    0.110
#    nonreactivity     0.475    0.065    7.290    0.000    0.475    0.475
#    acting           -0.113    0.077   -1.461    0.144   -0.113   -0.113
#  describing ~~                                                         
#    nonjudging       -0.435    0.075   -5.761    0.000   -0.435   -0.435
#    nonreactivity     0.367    0.070    5.227    0.000    0.367    0.367
#    acting           -0.509    0.069   -7.383    0.000   -0.509   -0.509
#  nonjudging ~~                                                         
#    nonreactivity    -0.252    0.083   -3.050    0.002   -0.252   -0.252
#    acting            0.673    0.044   15.479    0.000    0.673    0.673
#  nonreactivity ~~                                                      
#    acting           -0.291    0.082   -3.556    0.000   -0.291   -0.291

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# Analysis with orthogonal factors (5 factors)
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.795       0.804
#Tucker-Lewis Index (TLI)                       0.784       0.793
#Robust Comparative Fit Index (CFI)                         0.811
#Robust Tucker-Lewis Index (TLI)                            0.801
#RMSEA                                          0.087       0.075
#Robust RMSEA                                               0.083
#SRMR                                           0.213       0.213

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median = .537


# Bifactor model with 5 factors
BIFmodel= '
 observing =~  EM_1+EM_6+EM_11+EM_15+EM_20+EM_26+EM_31+EM_36
 describing =~  EM_2+EM_7+EM_12+EM_16+EM_22+EM_27+EM_32+EM_37
 nonjudging =~  EM_3+EM_10+EM_14+EM_17+EM_25+EM_30+EM_35+EM_39
 nonreactivity =~  EM_4+EM_9+EM_19+EM_21+EM_24+EM_29+EM_33
 acting =~  EM_5+EM_8+EM_13+EM_18+EM_23+EM_28+EM_34+EM_38
 general =~  EM_1+EM_2+EM_3+EM_4+EM_5+EM_6+EM_7+EM_8+EM_9+EM_10+EM_11+EM_12+EM_13+EM_14+EM_15+
                EM_16+EM_17+EM_18+EM_19+EM_20+EM_21+EM_22+EM_23+EM_24+EM_25+EM_26+EM_27+EM_28+EM_29+EM_30+
                EM_31+EM_32+EM_33+EM_34+EM_35+EM_36+EM_37+EM_38+EM_39						
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.876       0.888
#Tucker-Lewis Index (TLI)                       0.862       0.875
#Robust Comparative Fit Index (CFI)                         0.893
#Robust Tucker-Lewis Index (TLI)                            0.880
#RMSEA                                          0.070       0.058
#Robust RMSEA                                               0.064
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .595

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3723357      0.8205128      0.9056296      0.4307381 
#
#$FactorLevelIndices
#                 ECV_SS     ECV_SG     ECV_GS     Omega      OmegaH         H        FD
#observing     0.9863983 0.16879923 0.01360172 0.8847190 0.877950608 0.8913593 0.9447450
#describing    0.7144353 0.15459263 0.28556473 0.5605526 0.558088107 0.8813424 0.9474771
#nonjudging    0.5342474 0.12123674 0.46575263 0.9375688 0.502022424 0.8156512 0.9279971
#nonreactivity 0.9175775 0.14394550 0.08242253 0.8816409 0.810737102 0.8666441 0.9340304
#acting        0.1709358 0.03909016 0.82906416 0.9279843 0.002302805 0.5175222 0.8643479
#general       0.3723357 0.37233574 0.37233574 0.9056296 0.430738118 0.9398080 0.9681489




