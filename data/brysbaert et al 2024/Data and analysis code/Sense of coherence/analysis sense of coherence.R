################################################################
### Analysis Questionnaires Sense of coherence
### Zajenkowska et al. Polish study with 29 items
### Research on COVID situation in Poland https://osf.io/kney3

library(haven)
Coherence_Zajenkowska <- read_sav("Coherence_Zajenkowska.sav")
colnames(Coherence_Zajenkowska)
mydata  <- as.data.frame(Coherence_Zajenkowska[,86:114])
print(colnames(mydata))
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factor and 3 components
# Eigenvalue 1 = 7.31; eigenvalue 2 = 2.4

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 7.88; eigenvalue 2 = 2.56

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.25, RMSEA=.104, RMSR=.1, TLI=.577

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.27, RMSEA=.115, RMSR=.11, TLI=.553

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities 

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.059, RMSR=.04, TLI=.863
#     MR1  MR2  MR3
#MR1 1.00 0.37 0.15
#MR2 0.37 1.00 0.20
#MR3 0.15 0.20 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.068, RMSR=.04, TLI=.842
#     MR1  MR2  MR3
#MR1 1.00 0.37 0.14
#MR2 0.37 1.00 0.20
#MR3 0.14 0.20 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ soc_29_11_r1+soc_29_14_r1+soc_29_13_r1+soc_29_23_r1+soc_29_7_r1+soc_29_20_r1+soc_29_27_r1+
            soc_29_16_r1+soc_29_22_r1+soc_29_8_r1+soc_29_4_r1+soc_29_1_r1+soc_29_19_r1+soc_29_24_r1+
            soc_29_21_r1+soc_29_15_r1+soc_29_18_r1+soc_29_12_r1+soc_29_10_r1+soc_29_29_r1+soc_29_28_r1+
            soc_29_26_r1+soc_29_17_r1+soc_29_9_r1+soc_29_2_r1+soc_29_3_r1+soc_29_6_r1+soc_29_5_r1+
            soc_29_25_r1
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.866       0.722
#Tucker-Lewis Index (TLI)                       0.855       0.701
#Robust Comparative Fit Index (CFI)                         0.579
#Robust Tucker-Lewis Index (TLI)                            0.547
#RMSEA                                          0.158       0.118
#Robust RMSEA                                               0.118
#SRMR                                           0.111       0.111

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .297

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.609       0.616
#Tucker-Lewis Index (TLI)                       0.579       0.586
#Robust Comparative Fit Index (CFI)                         0.620
#Robust Tucker-Lewis Index (TLI)                            0.591
#RMSEA                                          0.104       0.086
#Robust RMSEA                                               0.102
#SRMR                                           0.102       0.102

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .237

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on EGA analysis
EGAmodel= '
 factor1 =~ soc_29_11_r1+soc_29_14_r1+soc_29_13_r1+soc_29_23_r1+soc_29_7_r1+soc_29_20_r1+soc_29_27_r1+
            soc_29_16_r1+soc_29_22_r1+soc_29_8_r1+soc_29_4_r1+soc_29_1_r1
 factor2 =~ soc_29_19_r1+soc_29_24_r1+
            soc_29_21_r1+soc_29_15_r1+soc_29_18_r1+soc_29_12_r1+soc_29_10_r1+soc_29_29_r1+soc_29_28_r1+
            soc_29_26_r1+soc_29_17_r1+soc_29_9_r1+soc_29_2_r1+soc_29_3_r1
 factor3 =~ soc_29_6_r1+soc_29_5_r1+
            soc_29_25_r1
'
library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.935       0.812
#Tucker-Lewis Index (TLI)                       0.930       0.795
#Robust Comparative Fit Index (CFI)                         0.719
#Robust Tucker-Lewis Index (TLI)                            0.695
#RMSEA                                          0.110       0.098
#Robust RMSEA                                               0.097
#SRMR                                           0.082       0.082

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .427

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.552    0.028   19.482    0.000    0.552    0.552
#    factor3           0.465    0.027   17.146    0.000    0.465    0.465
#  factor2 ~~                                                            
#    factor3           0.509    0.035   14.385    0.000    0.509    0.509

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.798       0.807
#Tucker-Lewis Index (TLI)                       0.781       0.790
#Robust Comparative Fit Index (CFI)                         0.812
#Robust Tucker-Lewis Index (TLI)                            0.796
#RMSEA                                          0.075       0.061
#Robust RMSEA                                               0.072
#SRMR                                           0.088       0.088

#Covariances:
#  Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
#  factor2           0.552    0.028   19.482    0.000    0.552    0.552
#  factor3           0.465    0.027   17.146    0.000    0.465    0.465
#factor2 ~~                                                            
#  factor3           0.509    0.035   14.385    0.000    0.509    0.509

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .395


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.761       0.770
#Tucker-Lewis Index (TLI)                       0.743       0.752
#Robust Comparative Fit Index (CFI)                         0.775
#Robust Tucker-Lewis Index (TLI)                            0.758
#RMSEA                                          0.082       0.066
#Robust RMSEA                                               0.078
#SRMR                                           0.160       0.160

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .40

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1 =~ soc_29_11_r1+soc_29_14_r1+soc_29_13_r1+soc_29_23_r1+soc_29_7_r1+soc_29_20_r1+soc_29_27_r1+
            soc_29_16_r1+soc_29_22_r1+soc_29_8_r1+soc_29_4_r1+soc_29_1_r1
 factor2 =~ soc_29_19_r1+soc_29_24_r1+
            soc_29_21_r1+soc_29_15_r1+soc_29_18_r1+soc_29_12_r1+soc_29_10_r1+soc_29_29_r1+soc_29_28_r1+
            soc_29_26_r1+soc_29_17_r1+soc_29_9_r1+soc_29_2_r1+soc_29_3_r1
 factor3 =~ soc_29_6_r1+soc_29_5_r1+
            soc_29_25_r1
 general=~ soc_29_11_r1+soc_29_14_r1+soc_29_13_r1+soc_29_23_r1+soc_29_7_r1+soc_29_20_r1+soc_29_27_r1+
            soc_29_16_r1+soc_29_22_r1+soc_29_8_r1+soc_29_4_r1+soc_29_1_r1+soc_29_19_r1+soc_29_24_r1+
            soc_29_21_r1+soc_29_15_r1+soc_29_18_r1+soc_29_12_r1+soc_29_10_r1+soc_29_29_r1+soc_29_28_r1+
            soc_29_26_r1+soc_29_17_r1+soc_29_9_r1+soc_29_2_r1+soc_29_3_r1+soc_29_6_r1+soc_29_5_r1+
            soc_29_25_r1
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.854       0.851
#Tucker-Lewis Index (TLI)                       0.829       0.826
#Robust Comparative Fit Index (CFI)                         0.864
#Robust Tucker-Lewis Index (TLI)                            0.841
#RMSEA                                          0.066       0.056
#Robust RMSEA                                               0.063
#SRMR                                           0.063       0.063

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .406

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 =~                                                            
#    soc_29_11_r1      0.781    0.274    2.855    0.004    0.781    0.581
#    soc_29_14_r1      0.767    0.281    2.732    0.006    0.767    0.465
#    soc_29_13_r1      0.864    0.113    7.661    0.000    0.864    0.596
#    soc_29_23_r1      0.690    0.279    2.472    0.013    0.690    0.406
#    soc_29_7_r1       0.633    0.265    2.386    0.017    0.633    0.427
#    soc_29_20_r1      0.699    0.107    6.525    0.000    0.699    0.462
#    soc_29_27_r1      0.562    0.124    4.517    0.000    0.562    0.408
#    soc_29_16_r1      0.570    0.107    5.317    0.000    0.570    0.434
#    soc_29_22_r1      0.160    0.532    0.300    0.764    0.160    0.109
#    soc_29_8_r1       0.115    0.487    0.236    0.814    0.115    0.083
#    soc_29_4_r1       0.233    0.231    1.009    0.313    0.233    0.149
#    soc_29_1_r1       0.249    0.140    1.780    0.075    0.249    0.163
#  factor2 =~                                                            
#    soc_29_19_r1      0.821    0.198    4.144    0.000    0.821    0.560
#    soc_29_24_r1      0.777    0.131    5.939    0.000    0.777    0.541
#    soc_29_21_r1      0.833    0.234    3.555    0.000    0.833    0.524
#    soc_29_15_r1      0.518    0.158    3.287    0.001    0.518    0.421
#    soc_29_18_r1      0.676    0.200    3.373    0.001    0.676    0.392
#    soc_29_12_r1      0.491    0.252    1.950    0.051    0.491    0.338
#    soc_29_10_r1      0.681    0.144    4.716    0.000    0.681    0.438
#    soc_29_29_r1      0.527    0.246    2.141    0.032    0.527    0.360
#    soc_29_28_r1      0.357    0.285    1.250    0.211    0.357    0.235
#    soc_29_26_r1      0.355    0.160    2.222    0.026    0.355    0.295
#    soc_29_17_r1      0.438    0.132    3.318    0.001    0.438    0.314
#    soc_29_9_r1       0.385    0.283    1.362    0.173    0.385    0.243
#    soc_29_2_r1       0.048    0.244    0.196    0.844    0.048    0.034
#    soc_29_3_r1      -0.007    0.204   -0.037    0.971   -0.007   -0.006
#  factor3 =~                                                            
#    soc_29_6_r1       1.134    0.107   10.608    0.000    1.134    0.896
#    soc_29_5_r1       0.777    0.075   10.302    0.000    0.777    0.676
#    soc_29_25_r1      0.364    0.136    2.676    0.007    0.364    0.269
#  general =~                                                            
#    soc_29_11_r1      0.622    0.277    2.243    0.025    0.622    0.463
#    soc_29_14_r1      1.149    0.199    5.776    0.000    1.149    0.696
#    soc_29_13_r1      0.628    0.204    3.082    0.002    0.628    0.433
#    soc_29_23_r1      0.919    0.226    4.067    0.000    0.919    0.541
#    soc_29_7_r1       0.724    0.219    3.313    0.001    0.724    0.488
#    soc_29_20_r1      0.689    0.150    4.599    0.000    0.689    0.456
#    soc_29_27_r1      0.763    0.102    7.444    0.000    0.763    0.554
#    soc_29_16_r1      0.657    0.116    5.646    0.000    0.657    0.499
#    soc_29_22_r1      1.017    0.232    4.380    0.000    1.017    0.697
#    soc_29_8_r1       0.739    0.239    3.095    0.002    0.739    0.535
#    soc_29_4_r1       0.618    0.133    4.640    0.000    0.618    0.395
#    soc_29_1_r1       0.692    0.081    8.487    0.000    0.692    0.452
#    soc_29_19_r1      0.631    0.228    2.764    0.006    0.631    0.431
#    soc_29_24_r1      0.414    0.202    2.048    0.041    0.414    0.289
#    soc_29_21_r1      0.745    0.244    3.057    0.002    0.745    0.468
#    soc_29_15_r1      0.573    0.130    4.401    0.000    0.573    0.466
#    soc_29_18_r1      0.786    0.160    4.922    0.000    0.786    0.456
#    soc_29_12_r1      0.691    0.182    3.802    0.000    0.691    0.475
#    soc_29_10_r1      0.159    0.209    0.760    0.447    0.159    0.102
#    soc_29_29_r1      0.771    0.180    4.270    0.000    0.771    0.526
#    soc_29_28_r1      1.008    0.119    8.446    0.000    1.008    0.663
#    soc_29_26_r1      0.408    0.123    3.325    0.001    0.408    0.339
#    soc_29_17_r1      0.194    0.120    1.618    0.106    0.194    0.140
#    soc_29_9_r1       0.663    0.188    3.520    0.000    0.663    0.419
#    soc_29_2_r1       0.442    0.114    3.861    0.000    0.442    0.317
#    soc_29_3_r1       0.491    0.078    6.302    0.000    0.491    0.390
#    soc_29_6_r1       0.387    0.147    2.628    0.009    0.387    0.306
#    soc_29_5_r1       0.227    0.100    2.269    0.023    0.227    0.198
#    soc_29_25_r1      0.701    0.154    4.548    0.000    0.701    0.519

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.5428460      0.6059113      0.9206976      0.7245832 

#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3634505 0.1671113 0.6365495 0.8932392 0.2880082 0.7164642 0.8063855
#factor2 0.4453254 0.1723411 0.5546746 0.8440081 0.3568228 0.7135631 0.8178435
#factor3 0.7682435 0.1177016 0.2317565 0.7778750 0.5945272 0.8329812 0.9377808
#general 0.5428460 0.5428460 0.5428460 0.9206976 0.7245832 0.8976039 0.9070113






################################################################
### Analysis Questionnaires Sense of coherence
### Marzena Lelek-Kratiuk https://onlinelibrary.wiley.com/doi/abs/10.1111/sjop.12813
### Dataset Szczygieł https://osf.io/rh2ja

library(readxl)
Sense_of_coherence_Lelek <- read_excel("Sense_of_coherence_Lelek.xlsx")
colnames(Sense_of_coherence_Lelek)
mydata  <- as.data.frame(Sense_of_coherence_Lelek[,2:30])
print(colnames(mydata))
mydata <- na.omit(mydata)
min(mydata)
max(mydata) #7 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 7.72; eigenvalue 2 = 1.91

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# # Eigenvalue 1 = 8.47; eigenvalue 2 = 2.01

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.27, RMSEA=.089, RMSR=.09, TLI=.667

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.29, RMSEA=.104, RMSR=.09, TLI=.625

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities but one with only two items

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.37, RMSEA=.057, RMSR=.05, TLI=.859
#     MR1  MR2  MR3
#MR1 1.00 0.37 0.25
#MR2 0.37 1.00 0.37
#MR3 0.25 0.37 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.074, RMSR=.05, TLI=.808
#     MR1  MR2  MR3
#MR1 1.00 0.37 0.28
#MR2 0.37 1.00 0.36
#MR3 0.28 0.36 1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ SOC1+SOC2+SOC3+SOC4+SOC5+SOC6+SOC7+SOC8+SOC9+SOC10+SOC11+SOC12+SOC13+SOC14+SOC15+
                SOC16+SOC17+SOC18+SOC19+SOC20+SOC21+SOC22+SOC23+SOC24+SOC25+SOC26+SOC27+SOC28+SOC29
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.924       0.790
#Tucker-Lewis Index (TLI)                       0.918       0.774
#Robust Comparative Fit Index (CFI)                         0.657
#Robust Tucker-Lewis Index (TLI)                            0.631
#RMSEA                                          0.110       0.104
#Robust RMSEA                                               0.107
#SRMR                                           0.094       0.094

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .33

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.688       0.693
#Tucker-Lewis Index (TLI)                       0.665       0.669
#Robust Comparative Fit Index (CFI)                         0.697
#Robust Tucker-Lewis Index (TLI)                            0.673
#RMSEA                                          0.091       0.084
#Robust RMSEA                                               0.089
#SRMR                                           0.086       0.086

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .307

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on EGA analysis
EGAmodel= '
 factor1=~ SOC7+SOC11+SOC14+SOC22+SOC16+SOC23+SOC8+SOC20+SOC13+SOC27+SOC28+SOC4+SOC2+SOC3
 factor2=~ SOC12+SOC10+SOC24+SOC21+SOC19+SOC15+SOC26+SOC29+SOC17+SOC25+SOC18+SOC1
 factor3=~ SOC6+SOC5+SOC9
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.961       0.881
#Tucker-Lewis Index (TLI)                       0.958       0.871
#Robust Comparative Fit Index (CFI)                         0.778
#Robust Tucker-Lewis Index (TLI)                            0.759
#RMSEA                                          0.079       0.079
#Robust RMSEA                                               0.086
#SRMR                                           0.075       0.075

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#factor1 ~~                                                            
# factor2           0.698    0.031   22.709    0.000    0.698    0.698
# factor3           0.439    0.050    8.739    0.000    0.439    0.439
#factor2 ~~                                                            
# factor3           0.572    0.048   11.844    0.000    0.572    0.572

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .42

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.823       0.830
#Tucker-Lewis Index (TLI)                       0.808       0.816
#Robust Comparative Fit Index (CFI)                         0.833
#Robust Tucker-Lewis Index (TLI)                            0.819
#RMSEA                                          0.069       0.063
#Robust RMSEA                                               0.066
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .370

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#    factor2           0.655    0.045   14.440    0.000    0.655    0.655
#    factor3           0.393    0.060    6.553    0.000    0.393    0.393
#  factor2 ~~                                                            
#    factor3           0.543    0.071    7.649    0.000    0.543    0.543


CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.762       0.767
#Tucker-Lewis Index (TLI)                       0.744       0.749
#Robust Comparative Fit Index (CFI)                         0.771
#Robust Tucker-Lewis Index (TLI)                            0.753
#RMSEA                                          0.079       0.073
#Robust RMSEA                                               0.077
#SRMR                                           0.175       0.175

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .354

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 factor1=~ SOC7+SOC11+SOC14+SOC22+SOC16+SOC23+SOC8+SOC20+SOC13+SOC27+SOC28+SOC4+SOC2+SOC3
 factor2=~ SOC12+SOC10+SOC24+SOC21+SOC19+SOC15+SOC26+SOC29+SOC17+SOC25+SOC18+SOC1
 factor3=~ SOC6+SOC5+SOC9
 general=~ SOC1+SOC2+SOC3+SOC4+SOC5+SOC6+SOC7+SOC8+SOC9+SOC10+SOC11+SOC12+SOC13+SOC14+SOC15+
                SOC16+SOC17+SOC18+SOC19+SOC20+SOC21+SOC22+SOC23+SOC24+SOC25+SOC26+SOC27+SOC28+SOC29
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.882       0.887
#Tucker-Lewis Index (TLI)                       0.862       0.869
#Robust Comparative Fit Index (CFI)                         0.890
#Robust Tucker-Lewis Index (TLI)                            0.872
#RMSEA                                          0.058       0.053
#Robust RMSEA                                               0.056
#SRMR                                           0.057       0.057

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .377

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5873328      0.6059113      0.9232034      0.7523368 
#
#$FactorLevelIndices
#           ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#factor1 0.5154480 0.23872479 0.4845520 0.8819048 0.43928356 0.7976875 0.8858626
#factor2 0.2231088 0.09311221 0.7768912 0.8584698 0.08775506 0.6532887 0.8138611
#factor3 0.6762917 0.08083017 0.3237083 0.6827897 0.43378675 0.6189485 0.8214527
#general 0.5873328 0.58733282 0.5873328 0.9232034 0.75233682 0.9091270 0.9397234




