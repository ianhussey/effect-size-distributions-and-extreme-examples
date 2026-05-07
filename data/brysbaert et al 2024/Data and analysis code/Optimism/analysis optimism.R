################################################################
### Analysis Optimism
### Coelho et al. 2018 https://www.sciencedirect.com/science/article/abs/pii/S0191886918303623
### For 2-factor model, take factorial structure from both UK and Brazilian group

library(haven)
Optimism_Coelho_UK <- read_sav("Optimism_Coelho_UK.sav")
colnames(Optimism_Coelho_UK)
mydata  <- as.matrix(Optimism_Coelho_UK[,11:19])
names <- colnames(mydata)
mydata <- na.omit(mydata)

Optimism_Coelho_Bra <- read_sav("Optimism_Coelho_Bra.sav")
colnames(Optimism_Coelho_Bra)
mydata2  <- as.matrix(Optimism_Coelho_Bra[,2:10])
colnames(mydata2) <- names
mydata2 <- na.omit(mydata2)

mydata_all <- rbind(mydata,mydata2)

library(psych)
fit <- fa(mydata_all,2)
fit
diagram(fit)
# 2,9,1,3,5,4 -> factor 1
# 7,6,8 -> factor 2

rm()



################################################################

### Coelho et al. 2018 https://pdf.sciencedirectassets.com/271782/1-s2.0-S0191886918X00099/1-s2.0-S0191886918303623/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEH8aCXVzLWVhc3QtMSJHMEUCIQCDS4QpaofSiVW%2BFHG8qHNzWK%2FU3zqoXiqNltIwA8OVNAIgPB6oxQLnptnMSC%2BA3Y8gBiK5cFLxCOMdV4Cp5LY2vEkquwUI2P%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDNxdNdxA9YyotQ5yKCqPBdHdZM%2F%2F%2Fj%2BotQx%2FbPfJs0bKKaRf2yXTOGdMR2u5ZI5NNR7kW%2BrNPGD5xH2JiXTWLQEd%2FCpUy43nVU7JIAl6tD7NcRtj2geMorr06%2FckH%2F7X3YEj08bllzF9%2BQFCG9M%2FMQTpHGlg076fbLZDdujjtJzwURl%2B99a9HVu7%2BaC%2Fdg4OFVbR9Aq8DZi4RRx8t2gySu7eXqtcCloXLZPUE%2FWBcb7jL5S6TT994h%2BQ5Jn5YVRzedD6v%2Fm%2B%2F%2Fipzkr21H%2FcNRFJv%2BkjMK5bTHDVkQXY1uhiOH5YThhTaV4N9wjiqCtyLa%2BhskhuFPj1xrIwYQvoLLny0JFS%2B3H%2FfwimXkRmpIoqa0Hiwh%2Bi%2BwHqgJudy7wkarIB1BKFdSkkBBSa410hD4TX3beUc2uJDh1bqx3Aj5G5G6MtFcpRct3bnkUEy%2FyI9xQhwEwdZGFnrHSDKqalCSPGRVYTGhiFdCSMZ1tdIJwLti6GMNEkqgCOQrb94qxhH%2Fc5wUCI7etYJ7NlQtXT4l0xbQiNML6C1a%2BClFbM6PX2E6rTe4xtZMVXnah4m0YAeNw10LWDUfazZJMmUNJdXsFPT10SZp4o6UjihTaMwGBlYV8V%2BgoDGPqb1GR5YB%2FsyNQJkYdrM%2BQU%2F5cHZPv2PsZeak%2FrAI4WeuEex15ov14609aKGlYh2avVIhsOm7%2FUkaJt9DJJmcg3byukIe99OcyN5Bs0lwzpRTx4jzDK38xpErsRGnpIeG20PIBQSBfBHTc3dnYbYxxUFvuI3BbGoE3qkCPQzKU1EQQne9sIOYYbwSOziWEmbcNUdSSQA3YCT%2BKDSReltHR3dQu54N3%2BsLVo4RcEfQmgEaM26RlEY6dkdcIJ7NKwNQGLLZsfr3Iwl8LRpAY6sQFr73KCLasp8o76%2FM6CmVmeuG2XeMT%2FK3olT61QtTMpOyPAhmwrWxkPLBU%2BDB11ajfpd4oDQE86k2US1xloXkKTYLnLvsNKjNwRa3WBnEsSwLbutuJJ8moqOjlOLse1eELnwvfLWVMMgn8Jze5jxdeJj6TswnGnDgJ%2B3IvE6Vt8MoWejGJTSAORHUzp71tCAzplO0hJ9gt472NOFQ1lrvgCddHmZJ9lBRFS7rfTxvklugo%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230622T152805Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQDYMSG3B%2F20230622%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=be147c12ca6fd5e1aa076585613a9e28e63f3cfb0ddacfa8e481a47d0fb7d5df&hash=92ae8e150cf24a3a7da7bbedac910e6dc8039a3cd4ea11a8a639f8e0c6defaeb&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0191886918303623&tid=spdf-6c91882a-2fe4-4725-83ce-42f6973f495a&sid=6551f5137aa8174ba29ac3f7d2fe0b003b38gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0703520b065759510504&rr=7db584dfcb512e56&cc=be
### UK sample https://osf.io/9kjxs/?view_only=d532eaf8b5234c8fb1b5eeeb27f0190f

min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 1 component
# Eigenvalue 1 = 5.54; eigenvalue 2 = 0.57

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 6.18; eigenvalue 2 = 0.37

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.62, RMSEA=.184, RMSR=.06, TLI=.821

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.69, RMSEA=.242, RMSR=.06, TLI=.77

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.68, RMSEA=.126, RMSR=.04, TLI=.917
#     MR1  MR2
#MR1 1.00 0.73
#MR2 0.73 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.75, RMSEA=.182, RMSR=.04, TLI=.87
#     MR1  MR2
#MR1 1.00 0.76
#MR2 0.76 1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi06+Otimi07+Otimi08+Otimi09
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.994       0.972
#Tucker-Lewis Index (TLI)                       0.992       0.962
#Robust Comparative Fit Index (CFI)                         0.796
#Robust Tucker-Lewis Index (TLI)                            0.728
#RMSEA                                          0.144       0.197
#Robust RMSEA                                               0.272
#SRMR                                           0.067       0.067

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .707

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.868       0.868
#Tucker-Lewis Index (TLI)                       0.823       0.824
#Robust Comparative Fit Index (CFI)                         0.872
#Robust Tucker-Lewis Index (TLI)                            0.829
#RMSEA                                          0.185       0.151
#Robust RMSEA                                               0.181
#SRMR                                           0.058       0.058

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .571

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on joint analysis UK and Brasil
EGAmodel= '
 factor1=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi09
 factor2=~  Otimi06+Otimi07+Otimi08
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.978
#Tucker-Lewis Index (TLI)                       0.994       0.970
#Robust Comparative Fit Index (CFI)                         0.858
#Robust Tucker-Lewis Index (TLI)                            0.803
#RMSEA                                          0.122       0.177
#Robust RMSEA                                               0.231
#SRMR                                           0.057       0.057

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#   factor2           0.894    0.014   62.579    0.000    0.894    0.894

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .755

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.907       0.908
#Tucker-Lewis Index (TLI)                       0.871       0.873
#Robust Comparative Fit Index (CFI)                         0.911
#Robust Tucker-Lewis Index (TLI)                            0.877
#RMSEA                                          0.158       0.129
#Robust RMSEA                                               0.154
#SRMR                                           0.058       0.058

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#   factor2           0.854    0.032   26.961    0.000    0.854    0.854

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .671

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.785       0.780
#Tucker-Lewis Index (TLI)                       0.713       0.706
#Robust Comparative Fit Index (CFI)                         0.788
#Robust Tucker-Lewis Index (TLI)                            0.718
#RMSEA                                          0.235       0.196
#Robust RMSEA                                               0.233
#SRMR                                           0.377       0.377

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# min %variance given = .611


# Bifactor model
BIFmodel= '
 factor1=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi09
 factor2=~  Otimi07+Otimi06+Otimi08
 general=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi06+Otimi07+Otimi08+Otimi09
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.972       0.974
#Tucker-Lewis Index (TLI)                       0.943       0.948
#Robust Comparative Fit Index (CFI)                         0.975
#Robust Tucker-Lewis Index (TLI)                            0.950
#RMSEA                                          0.105       0.083
#Robust RMSEA                                               0.098
#SRMR                                           0.032       0.032

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .612

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#    0.03764379     0.50000000     4.72713567     0.83513412 
#
#$FactorLevelIndices
#            ECV_SS      ECV_SG     ECV_GS      Omega     OmegaH           H         FD
#factor1 0.22004822 0.006564727 0.77995178  0.8527633  0.1880442   0.6444823  1.0787648
#factor2 0.98518256 0.955791483 0.01481744 19.0270877 18.2114516 134.9966665 18.6598970
#general 0.03764379 0.037643790 0.03764379  4.7271357  0.8351341   0.9327244  0.9572119







################################################################

### Coelho et al. 2018 https://pdf.sciencedirectassets.com/271782/1-s2.0-S0191886918X00099/1-s2.0-S0191886918303623/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEH8aCXVzLWVhc3QtMSJHMEUCIQCDS4QpaofSiVW%2BFHG8qHNzWK%2FU3zqoXiqNltIwA8OVNAIgPB6oxQLnptnMSC%2BA3Y8gBiK5cFLxCOMdV4Cp5LY2vEkquwUI2P%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDNxdNdxA9YyotQ5yKCqPBdHdZM%2F%2F%2Fj%2BotQx%2FbPfJs0bKKaRf2yXTOGdMR2u5ZI5NNR7kW%2BrNPGD5xH2JiXTWLQEd%2FCpUy43nVU7JIAl6tD7NcRtj2geMorr06%2FckH%2F7X3YEj08bllzF9%2BQFCG9M%2FMQTpHGlg076fbLZDdujjtJzwURl%2B99a9HVu7%2BaC%2Fdg4OFVbR9Aq8DZi4RRx8t2gySu7eXqtcCloXLZPUE%2FWBcb7jL5S6TT994h%2BQ5Jn5YVRzedD6v%2Fm%2B%2F%2Fipzkr21H%2FcNRFJv%2BkjMK5bTHDVkQXY1uhiOH5YThhTaV4N9wjiqCtyLa%2BhskhuFPj1xrIwYQvoLLny0JFS%2B3H%2FfwimXkRmpIoqa0Hiwh%2Bi%2BwHqgJudy7wkarIB1BKFdSkkBBSa410hD4TX3beUc2uJDh1bqx3Aj5G5G6MtFcpRct3bnkUEy%2FyI9xQhwEwdZGFnrHSDKqalCSPGRVYTGhiFdCSMZ1tdIJwLti6GMNEkqgCOQrb94qxhH%2Fc5wUCI7etYJ7NlQtXT4l0xbQiNML6C1a%2BClFbM6PX2E6rTe4xtZMVXnah4m0YAeNw10LWDUfazZJMmUNJdXsFPT10SZp4o6UjihTaMwGBlYV8V%2BgoDGPqb1GR5YB%2FsyNQJkYdrM%2BQU%2F5cHZPv2PsZeak%2FrAI4WeuEex15ov14609aKGlYh2avVIhsOm7%2FUkaJt9DJJmcg3byukIe99OcyN5Bs0lwzpRTx4jzDK38xpErsRGnpIeG20PIBQSBfBHTc3dnYbYxxUFvuI3BbGoE3qkCPQzKU1EQQne9sIOYYbwSOziWEmbcNUdSSQA3YCT%2BKDSReltHR3dQu54N3%2BsLVo4RcEfQmgEaM26RlEY6dkdcIJ7NKwNQGLLZsfr3Iwl8LRpAY6sQFr73KCLasp8o76%2FM6CmVmeuG2XeMT%2FK3olT61QtTMpOyPAhmwrWxkPLBU%2BDB11ajfpd4oDQE86k2US1xloXkKTYLnLvsNKjNwRa3WBnEsSwLbutuJJ8moqOjlOLse1eELnwvfLWVMMgn8Jze5jxdeJj6TswnGnDgJ%2B3IvE6Vt8MoWejGJTSAORHUzp71tCAzplO0hJ9gt472NOFQ1lrvgCddHmZJ9lBRFS7rfTxvklugo%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20230622T152805Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYQDYMSG3B%2F20230622%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=be147c12ca6fd5e1aa076585613a9e28e63f3cfb0ddacfa8e481a47d0fb7d5df&hash=92ae8e150cf24a3a7da7bbedac910e6dc8039a3cd4ea11a8a639f8e0c6defaeb&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0191886918303623&tid=spdf-6c91882a-2fe4-4725-83ce-42f6973f495a&sid=6551f5137aa8174ba29ac3f7d2fe0b003b38gxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0703520b065759510504&rr=7db584dfcb512e56&cc=be
### BR sample https://osf.io/9kjxs/?view_only=d532eaf8b5234c8fb1b5eeeb27f0190f

mydata <- mydata2
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 2 factors and 1 component
# Eigenvalue 1 = 3.75; eigenvalue 2 = 0.53

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 1 component
# Eigenvalue 1 = 4.34; eigenvalue 2 = 0.43

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.42, RMSEA=.142, RMSR=.07, TLI=.791

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained= .48, RMSEA=.185, RMSR=.07, TLI=.742

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 1 community

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.091, RMSR=.04, TLI=.914
#     MR1  MR2
#MR1 1.00 0.66
#MR2 0.66 1.00

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.128, RMSR=.04, TLI=.876
#     MR1  MR2
#MR1 1.00 0.66
#MR2 0.66 1.00

colnames(mydata)

# Single factor model lavaan
UNImodel= '
 general=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi06+Otimi07+Otimi08+Otimi09
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.929
#Tucker-Lewis Index (TLI)                       0.974       0.905
#Robust Comparative Fit Index (CFI)                         0.789
#Robust Tucker-Lewis Index (TLI)                            0.718
#RMSEA                                          0.130       0.178
#Robust RMSEA                                               0.196
#SRMR                                           0.076       0.076

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .504

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.844       0.831
#Tucker-Lewis Index (TLI)                       0.791       0.774
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.795
#RMSEA                                          0.143       0.127
#Robust RMSEA                                               0.141
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .464

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on joint analysis UK and Brasil
EGAmodel= '
 factor1=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi09
 factor2=~  Otimi06+Otimi07+Otimi08
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.956
#Tucker-Lewis Index (TLI)                       0.985       0.939
#Robust Comparative Fit Index (CFI)                         0.885
#Robust Tucker-Lewis Index (TLI)                            0.841
#RMSEA                                          0.099       0.143
#Robust RMSEA                                               0.147
#SRMR                                           0.060       0.060

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#   factor2           0.820    0.020   40.325    0.000    0.820    0.820

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .58

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.910       0.907
#Tucker-Lewis Index (TLI)                       0.876       0.871
#Robust Comparative Fit Index (CFI)                         0.914
#Robust Tucker-Lewis Index (TLI)                            0.881
#RMSEA                                          0.111       0.096
#Robust RMSEA                                               0.108
#SRMR                                           0.064       0.064

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  factor1 ~~                                                            
#   factor2           0.761    0.035   21.554    0.000    0.761    0.761

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .369

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.762       0.752
#Tucker-Lewis Index (TLI)                       0.683       0.670
#Robust Comparative Fit Index (CFI)                         0.765
#Robust Tucker-Lewis Index (TLI)                            0.687
#RMSEA                                          0.176       0.154
#Robust RMSEA                                               0.174
#SRMR                                           0.260       0.260

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# min %variance given = .315


# Bifactor model
BIFmodel= '
 factor1=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi09
 factor2=~  Otimi07+Otimi06+Otimi08
 general=~  Otimi01+Otimi02+Otimi03+Otimi04+Otimi05+Otimi06+Otimi07+Otimi08+Otimi09
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.981       0.983
#Tucker-Lewis Index (TLI)                       0.963       0.965
#Robust Comparative Fit Index (CFI)                         0.984
#Robust Tucker-Lewis Index (TLI)                            0.968
#RMSEA                                          0.061       0.050
#Robust RMSEA                                               0.056
#SRMR                                           0.031       0.031

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .452

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6894471      0.5000000      0.8500378      0.7243635 
#
#$FactorLevelIndices
#           ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#factor1 0.3348073 0.1767693 0.6651927 0.6840454 0.2397763 0.5871059 0.8289946
#factor2 0.2834238 0.1337837 0.7165762 0.8982959 0.1264839 0.5582114 1.0124414
#general 0.6894471 0.6894471 0.6894471 0.8500378 0.7243635 0.8713634 0.9293554






