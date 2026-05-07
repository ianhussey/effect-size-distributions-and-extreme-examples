################################################################
### Analysis Perfectionism 



################################################################
### Mackinnon et al https://osf.io/gduy4/
### Big Three Perfectionism Scale (BTPS) https://cruxpsychology.ca/wp-content/uploads/2017/07/Smithetal.2016.670-687.pdf 

MacKinnon_Perfectionism1 <- read.csv("MacKinnon_Perfectionism1.csv")
colnames(MacKinnon_Perfectionism1)
mydata <- as.data.frame(MacKinnon_Perfectionism1[,259:304])
colnames(mydata)
mydata <- na.omit(mydata)
mydata <- subset(mydata, select = -c(rigid.mean,self.critical.mean))
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 9 factors and 7 components
# Eigenvalue 1 = 14.6; eigenvalue 2 = 5.64

rho <- polychoric(mydata)$rho
# warning
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 7 components
# Eigenvalue 1 = 17.39; eigenvalue 2 = 6.68

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# does not work properly
# %variance explained=.33, RMSEA=.156, RMSR=.16, TLI=.387

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.40, RMSEA=.257, RMSR=.18, TLI=.212

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities with response bias

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.121, RMSR=.07, TLI=.632
#     MR2  MR1  MR3
#MR2 1.00 0.30 0.23
#MR1 0.30 1.00 0.56
#MR3 0.23 0.56 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.61, RMSEA=.234, RMSR=.07, TLI=.346
#     MR2  MR1  MR3
#MR2 1.00 0.31 0.25
#MR1 0.31 1.00 0.56
#MR3 0.25 0.56 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  btps.SOP1+btps.SOP2+btps.SOP3+btps.SOP4+btps.SOP5+btps.SWC1+btps.SWC2+btps.SWC3+ 
            btps.SWC4+btps.SWC5+btps.COM1+btps.COM2+btps.COM3+btps.COM4+btps.COM5+btps.DAA1+ 
            btps.DAA2+btps.DAA3+btps.DAA4+btps.SC1+btps.SC2+btps.SC3+btps.SC4+btps.SPP1+ 
            btps.SPP2+btps.SPP3+btps.SPP4+btps.OOP1+btps.OOP2+btps.OOP3+btps.OOP4+btps.OOP5+ 
            btps.HC1+btps.HC2+btps.HC3+btps.HC4+btps.ENT1+btps.ENT2+btps.ENT3+btps.ENT4+ 
            btps.GRAN1+btps.GRAN2+btps.GRAN3+btps.GRAN4
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.872       0.682
#Tucker-Lewis Index (TLI)                       0.866       0.667
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.256       0.163
#Robust RMSEA                                                  NA
#SRMR                                           0.208       0.208

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .531

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.422       0.421
#Tucker-Lewis Index (TLI)                       0.394       0.393
#Robust Comparative Fit Index (CFI)                         0.422
#Robust Tucker-Lewis Index (TLI)                            0.394
#RMSEA                                          0.155       0.145
#Robust RMSEA                                               0.155
#SRMR                                           0.158       0.158

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .36

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on theoretical analysis
EGAmodel= '
 rigid  =~  btps.SOP1+btps.SOP2+btps.SOP3+btps.SOP4+btps.SOP5+btps.SWC1+btps.SWC2+btps.SWC3+ 
            btps.SWC4+btps.SWC5
 selfcrit=~ btps.COM1+btps.COM2+btps.COM3+btps.COM4+btps.COM5+btps.DAA1+ 
            btps.DAA2+btps.DAA3+btps.DAA4+btps.SC1+btps.SC2+btps.SC3+btps.SC4+btps.SPP1+ 
            btps.SPP2+btps.SPP3+btps.SPP4
 narcist =~ btps.OOP1+btps.OOP2+btps.OOP3+btps.OOP4+btps.OOP5+ 
            btps.HC1+btps.HC2+btps.HC3+btps.HC4+btps.ENT1+btps.ENT2+btps.ENT3+btps.ENT4+ 
            btps.GRAN1+btps.GRAN2+btps.GRAN3+btps.GRAN4
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.959       0.878
#Tucker-Lewis Index (TLI)                       0.957       0.872
#Robust Comparative Fit Index (CFI)                            NA
#Robust Tucker-Lewis Index (TLI)                               NA
#RMSEA                                          0.146       0.101
#Robust RMSEA                                                  NA
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .656

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#rigid ~~                                                              
#  selfcrit          0.732    0.006  121.578    0.000    0.732    0.732
#  narcist           0.401    0.012   34.503    0.000    0.401    0.401
#selfcrit ~~                                                           
#  narcist           0.382    0.011   34.229    0.000    0.382    0.382

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.642       0.636
#Tucker-Lewis Index (TLI)                       0.623       0.617
#Robust Comparative Fit Index (CFI)                         0.642
#Robust Tucker-Lewis Index (TLI)                            0.624
#RMSEA                                          0.123       0.115
#Robust RMSEA                                               0.122
#SRMR                                           0.089       0.089

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#rigid ~~                                                              
#  selfcrit          0.732    0.006  121.578    0.000    0.732    0.732
#  narcist           0.401    0.012   34.503    0.000    0.401    0.401
#selfcrit ~~                                                           
#  narcist           0.382    0.011   34.229    0.000    0.382    0.382

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .532

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.619       0.612
#Tucker-Lewis Index (TLI)                       0.601       0.593
#Robust Comparative Fit Index (CFI)                         0.620
#Robust Tucker-Lewis Index (TLI)                            0.601
#RMSEA                                          0.126       0.118
#Robust RMSEA                                               0.126
#SRMR                                           0.225       0.225

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .536

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 rigid  =~  btps.SOP1+btps.SOP2+btps.SOP3+btps.SOP4+btps.SOP5+btps.SWC1+btps.SWC2+btps.SWC3+ 
            btps.SWC4+btps.SWC5
 selfcrit=~ btps.COM1+btps.COM2+btps.COM3+btps.COM4+btps.COM5+btps.DAA1+ 
            btps.DAA2+btps.DAA3+btps.DAA4+btps.SC1+btps.SC2+btps.SC3+btps.SC4+btps.SPP1+ 
            btps.SPP2+btps.SPP3+btps.SPP4
 narcist =~ btps.OOP1+btps.OOP2+btps.OOP3+btps.OOP4+btps.OOP5+ 
            btps.HC1+btps.HC2+btps.HC3+btps.HC4+btps.ENT1+btps.ENT2+btps.ENT3+btps.ENT4+ 
            btps.GRAN1+btps.GRAN2+btps.GRAN3+btps.GRAN4
 general=~  btps.SOP1+btps.SOP2+btps.SOP3+btps.SOP4+btps.SOP5+btps.SWC1+btps.SWC2+btps.SWC3+ 
            btps.SWC4+btps.SWC5+btps.COM1+btps.COM2+btps.COM3+btps.COM4+btps.COM5+btps.DAA1+ 
            btps.DAA2+btps.DAA3+btps.DAA4+btps.SC1+btps.SC2+btps.SC3+btps.SC4+btps.SPP1+ 
            btps.SPP2+btps.SPP3+btps.SPP4+btps.OOP1+btps.OOP2+btps.OOP3+btps.OOP4+btps.OOP5+ 
            btps.HC1+btps.HC2+btps.HC3+btps.HC4+btps.ENT1+btps.ENT2+btps.ENT3+btps.ENT4+ 
            btps.GRAN1+btps.GRAN2+btps.GRAN3+btps.GRAN4
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.689       0.678
#Tucker-Lewis Index (TLI)                       0.657       0.645
#Robust Comparative Fit Index (CFI)                         0.689
#Robust Tucker-Lewis Index (TLI)                            0.657
#RMSEA                                          0.117       0.111
#Robust RMSEA                                               0.117
#SRMR                                           0.078       0.078

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .543

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  rigid =~                                                              
#    btps.SOP1         0.637    0.019   33.336    0.000    0.637    0.523
#    btps.SOP2         0.586    0.016   36.449    0.000    0.586    0.526
#    btps.SOP3         0.508    0.017   29.771    0.000    0.508    0.472
#    btps.SOP4         0.614    0.018   34.278    0.000    0.614    0.511
#    btps.SOP5         0.469    0.021   22.202    0.000    0.469    0.417
#    btps.SWC1         0.639    0.019   33.338    0.000    0.639    0.522
#    btps.SWC2         0.430    0.021   20.375    0.000    0.430    0.376
#    btps.SWC3         0.336    0.026   12.799    0.000    0.336    0.307
#    btps.SWC4         0.543    0.029   18.493    0.000    0.543    0.455
#    btps.SWC5         0.367    0.026   14.294    0.000    0.367    0.311
#  selfcrit =~                                                           
#    btps.COM1         0.294    0.038    7.805    0.000    0.294    0.230
#    btps.COM2         0.183    0.039    4.726    0.000    0.183    0.158
#    btps.COM3         0.179    0.041    4.320    0.000    0.179    0.150
#    btps.COM4         0.263    0.039    6.725    0.000    0.263    0.218
#    btps.COM5         0.262    0.038    6.906    0.000    0.262    0.211
#    btps.DAA1         0.876    0.027   31.895    0.000    0.876    0.668
#    btps.DAA2         0.944    0.022   43.879    0.000    0.944    0.730
#    btps.DAA3         0.866    0.023   37.841    0.000    0.866    0.692
#    btps.DAA4         0.803    0.025   32.602    0.000    0.803    0.679
#    btps.SC1          0.211    0.027    7.766    0.000    0.211    0.172
#    btps.SC2          0.159    0.026    6.021    0.000    0.159    0.126
#    btps.SC3          0.173    0.027    6.317    0.000    0.173    0.140
#    btps.SC4          0.130    0.026    5.052    0.000    0.130    0.105
#    btps.SPP1         0.086    0.047    1.819    0.069    0.086    0.069
#    btps.SPP2         0.052    0.036    1.467    0.142    0.052    0.047
#    btps.SPP3         0.004    0.040    0.100    0.920    0.004    0.003
#    btps.SPP4        -0.041    0.042   -0.982    0.326   -0.041   -0.034
#  narcist =~                                                            
#    btps.OOP1         0.530    0.014   36.746    0.000    0.530    0.619
#    btps.OOP2         0.539    0.014   37.678    0.000    0.539    0.703
#    btps.OOP3         0.589    0.015   38.957    0.000    0.589    0.720
#    btps.OOP4         0.650    0.018   35.640    0.000    0.650    0.674
#    btps.OOP5         0.601    0.015   40.797    0.000    0.601    0.683
#    btps.HC1          0.632    0.017   37.613    0.000    0.632    0.560
#    btps.HC2          0.425    0.019   22.443    0.000    0.425    0.373
#    btps.HC3          0.573    0.015   37.938    0.000    0.573    0.548
#    btps.HC4          0.608    0.016   37.168    0.000    0.608    0.547
#    btps.ENT1         0.529    0.013   41.903    0.000    0.529    0.695
#    btps.ENT2         0.575    0.016   35.831    0.000    0.575    0.621
#    btps.ENT3         0.575    0.017   34.072    0.000    0.575    0.609
#    btps.ENT4         0.599    0.017   35.966    0.000    0.599    0.632
#    btps.GRAN1        0.521    0.014   37.058    0.000    0.521    0.652
#    btps.GRAN2        0.505    0.015   34.768    0.000    0.505    0.705
#    btps.GRAN3        0.532    0.014   37.636    0.000    0.532    0.623
#    btps.GRAN4        0.621    0.014   44.277    0.000    0.621    0.585
#  general =~                                                            
#    btps.SOP1         0.789    0.019   42.110    0.000    0.789    0.647
#    btps.SOP2         0.607    0.019   31.892    0.000    0.607    0.544
#    btps.SOP3         0.657    0.016   40.305    0.000    0.657    0.610
#    btps.SOP4         0.773    0.017   44.899    0.000    0.773    0.644
#    btps.SOP5         0.603    0.018   33.085    0.000    0.603    0.536
#    btps.SWC1         0.773    0.020   37.862    0.000    0.773    0.631
#    btps.SWC2         0.630    0.020   31.361    0.000    0.630    0.551
#    btps.SWC3         0.768    0.019   40.147    0.000    0.768    0.703
#    btps.SWC4         0.699    0.023   30.061    0.000    0.699    0.586
#    btps.SWC5         0.888    0.018   50.200    0.000    0.888    0.752
#    btps.COM1         0.953    0.021   45.650    0.000    0.953    0.746
#    btps.COM2         0.830    0.022   38.157    0.000    0.830    0.716
#    btps.COM3         0.829    0.023   36.003    0.000    0.829    0.694
#    btps.COM4         0.840    0.023   36.810    0.000    0.840    0.698
#    btps.COM5         0.944    0.019   49.232    0.000    0.944    0.759
#    btps.DAA1         0.716    0.032   22.554    0.000    0.716    0.545
#    btps.DAA2         0.705    0.031   22.718    0.000    0.705    0.545
#    btps.DAA3         0.626    0.030   20.952    0.000    0.626    0.500
#    btps.DAA4         0.580    0.030   19.136    0.000    0.580    0.491
#    btps.SC1          0.998    0.015   65.984    0.000    0.998    0.815
#    btps.SC2          0.996    0.016   61.095    0.000    0.996    0.786
#    btps.SC3          1.002    0.018   56.174    0.000    1.002    0.808
#    btps.SC4          1.017    0.015   67.056    0.000    1.017    0.823
#    btps.SPP1         0.636    0.029   21.895    0.000    0.636    0.511
#    btps.SPP2         0.604    0.022   27.605    0.000    0.604    0.546
#    btps.SPP3         0.612    0.025   24.581    0.000    0.612    0.507
#    btps.SPP4         0.611    0.025   23.998    0.000    0.611    0.501
#    btps.OOP1         0.325    0.015   21.432    0.000    0.325    0.379
#    btps.OOP2         0.272    0.015   18.466    0.000    0.272    0.354
#    btps.OOP3         0.275    0.015   18.372    0.000    0.275    0.336
#    btps.OOP4         0.251    0.015   17.203    0.000    0.251    0.260
#    btps.OOP5         0.349    0.017   20.699    0.000    0.349    0.397
#    btps.HC1          0.274    0.019   14.098    0.000    0.274    0.243
#    btps.HC2          0.298    0.021   13.966    0.000    0.298    0.261
#    btps.HC3          0.242    0.018   13.721    0.000    0.242    0.232
#    btps.HC4          0.151    0.021    7.094    0.000    0.151    0.136
#    btps.ENT1         0.193    0.014   14.129    0.000    0.193    0.254
#    btps.ENT2         0.203    0.014   14.176    0.000    0.203    0.219
#    btps.ENT3         0.303    0.015   19.921    0.000    0.303    0.321
#    btps.ENT4         0.186    0.015   12.439    0.000    0.186    0.196
#    btps.GRAN1        0.091    0.013    6.826    0.000    0.091    0.113
#    btps.GRAN2        0.126    0.012   10.668    0.000    0.126    0.176
#    btps.GRAN3        0.146    0.015   10.043    0.000    0.146    0.172
#    btps.GRAN4        0.142    0.017    8.280    0.000    0.142    0.134

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.5331950      0.6649049      0.9670070      0.7281569 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#rigid    0.3414372 0.08663706 0.6585628 0.9341543 0.3144219 0.7236296 0.8694145
#selfcrit 0.2293085 0.09422409 0.7706915 0.9495669 0.1293364 0.8002284 0.9172330
#narcist  0.8526672 0.28594388 0.1473328 0.9334164 0.8066846 0.9203706 0.9599148
#general  0.5331950 0.53319497 0.5331950 0.9670070 0.7281569 0.9601325 0.9664178







################################################################
### Workye et al https://osf.io/dxspb/
### Big Three Perfectionism Scale (BTPS) https://cruxpsychology.ca/wp-content/uploads/2017/07/Smithetal.2016.670-687.pdf 

Workye_Perfectionism <- read.csv("Workye_Perfectionism.csv")
colnames(Workye_Perfectionism)
mydata <- as.data.frame(Workye_Perfectionism[,37:74])
colnames(mydata)
mydata <- na.omit(mydata)

min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 4 components
# Eigenvalue 1 = 13.17; eigenvalue 2 = 2.97

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 components
# Eigenvalue 1 = 14.97; eigenvalue 2 = 3.24

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.35, RMSEA=.123, RMSR=.12, TLI=.54

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.39, RMSEA=.149, RMSR=.13, TLI=.493

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 5 communities with response bias

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.51, RMSEA=.082, RMSR=.05, TLI=.795
#     MR1  MR3  MR2
#MR1 1.00 0.41 0.33
#MR3 0.41 1.00 0.33
#MR2 0.33 0.33 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.108, RMSR=.05, TLI=.73
#     MR1  MR3  MR2
#MR1 1.00 0.42 0.35
#MR3 0.42 1.00 0.34
#MR2 0.35 0.34 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  trtper1+trtper2+trtper3+trtper4+trtper5+trtper6+trtper7+trtper8+trtper9+trtper10+
            trtper11+trtper12+trtper13+trtper14+trtper15+trtper16+trtper17+trtper18+trtper19+
            trtper20+trtper21+trtper22+trtper23+trtper24+trtper25+trtper26+trtper27+trtper28+
            trtper29+trtper30+trtper31+trtper32+trtper33+trtper34+trtper35+trtper36+trtper37+
            trtper38
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.898       0.690
#Tucker-Lewis Index (TLI)                       0.892       0.673
#Robust Comparative Fit Index (CFI)                         0.525
#Robust Tucker-Lewis Index (TLI)                            0.498
#RMSEA                                          0.195       0.150
#Robust RMSEA                                               0.152
#SRMR                                           0.143       0.143

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .449

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.565       0.566
#Tucker-Lewis Index (TLI)                       0.540       0.541
#Robust Comparative Fit Index (CFI)                         0.570
#Robust Tucker-Lewis Index (TLI)                            0.545
#RMSEA                                          0.126       0.117
#Robust RMSEA                                               0.124
#SRMR                                           0.118       0.118

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .334

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on theoretical and EGA analysis
EGAmodel= '
 rigid  =~  trtper1+trtper2+trtper3+trtper4+trtper5+trtper6+trtper7+trtper8+trtper9+trtper10
 selfcrit=~ trtper11+trtper12+trtper13+trtper14+trtper15+trtper16+trtper17+trtper18+trtper19+
            trtper20+trtper21+trtper22+trtper23+trtper24+trtper34+trtper35+trtper36+trtper37+
            trtper38
 expect  =~ trtper25+trtper26+trtper27+trtper28+trtper29+trtper30+trtper31+trtper32+trtper33
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.963       0.866
#Tucker-Lewis Index (TLI)                       0.960       0.858
#Robust Comparative Fit Index (CFI)                         0.725
#Robust Tucker-Lewis Index (TLI)                            0.708
#RMSEA                                          0.118       0.099
#Robust RMSEA                                               0.116
#SRMR                                           0.093       0.093

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#rigid ~~                                                              
#        selfcrit          0.634    0.027   23.287    0.000    0.634    0.634
#        expect            0.471    0.038   12.529    0.000    0.471    0.471
#selfcrit ~~                                                           
#        expect            0.483    0.036   13.275    0.000    0.483    0.483

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .594

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.771       0.776
#Tucker-Lewis Index (TLI)                       0.757       0.762
#Robust Comparative Fit Index (CFI)                         0.778
#Robust Tucker-Lewis Index (TLI)                            0.764
#RMSEA                                          0.091       0.084
#Robust RMSEA                                               0.090
#SRMR                                           0.087       0.08

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .490

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  rigid ~~                                                              
#    selfcrit          0.637    0.038   16.949    0.000    0.637    0.637
#    expect            0.449    0.047    9.534    0.000    0.449    0.449
#  selfcrit ~~                                                           
#    expect            0.427    0.053    8.060    0.000    0.427    0.427

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.745       0.749
#Tucker-Lewis Index (TLI)                       0.730       0.734
#Robust Comparative Fit Index (CFI)                         0.751
#Robust Tucker-Lewis Index (TLI)                            0.737
#RMSEA                                          0.096       0.089
#Robust RMSEA                                               0.095
#SRMR                                           0.222       0.222

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .485

# Bifactor model
BIFmodel= '
 rigid  =~  trtper1+trtper2+trtper3+trtper4+trtper5+trtper6+trtper7+trtper8+trtper9+trtper10
 selfcrit=~ trtper11+trtper12+trtper13+trtper14+trtper15+trtper16+trtper17+trtper18+trtper19+
            trtper20+trtper21+trtper22+trtper23+trtper24+trtper34+trtper35+trtper36+trtper37+
            trtper38
 expect  =~ trtper25+trtper26+trtper27+trtper28+trtper29+trtper30+trtper31+trtper32+trtper33
 general=~  trtper1+trtper2+trtper3+trtper4+trtper5+trtper6+trtper7+trtper8+trtper9+trtper10+
            trtper11+trtper12+trtper13+trtper14+trtper15+trtper16+trtper17+trtper18+trtper19+
            trtper20+trtper21+trtper22+trtper23+trtper24+trtper25+trtper26+trtper27+trtper28+
            trtper29+trtper30+trtper31+trtper32+trtper33+trtper34+trtper35+trtper36+trtper37+
            trtper38
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.856       0.859
#Tucker-Lewis Index (TLI)                       0.839       0.842
#Robust Comparative Fit Index (CFI)                         0.862
#Robust Tucker-Lewis Index (TLI)                            0.846
#RMSEA                                          0.075       0.068
#Robust RMSEA                                               0.072
#SRMR                                           0.060       0.060

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .536

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.5802227      0.6415363      0.9648176      0.8129851 
#
#$FactorLevelIndices
#            ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#rigid    0.4853747 0.1228789 0.5146253 0.9083104 0.4246160 0.7926519 0.9038076
#selfcrit 0.2407987 0.1251566 0.7592013 0.9526777 0.1404820 0.8199884 0.9157720
#expect   0.7563023 0.1717419 0.2436977 0.8976869 0.6712794 0.8746369 0.9423077
#general  0.5802227 0.5802227 0.5802227 0.9648176 0.8129851 0.9544550 0.9671214







################################################################
### Linnett & Kibowski https://osf.io/n2fhu/
### self developed questionnaire from existing questionnaires; measures 3 components
### (concerns over mistakes, striving for excellence, discrepancy)

library(haven)
Perfectionism_Linnett <- read_sav("Perfectionism_Linnett.sav")
colnames(Perfectionism_Linnett)
mydata <- as.data.frame(Perfectionism_Linnett[,7:36])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 3 components
# Eigenvalue 1 = 12.96; eigenvalue 2 = 2.31

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 14.55; eigenvalue 2 = 2.51

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.43, RMSEA=.116, RMSR=.1, TLI=.691

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.142, RMSR=.11, TLI=.642

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities (does not work with response bias)

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.57, RMSEA=.056, RMSR=.03, TLI=.929
#     MR1  MR2  MR3
#MR1 1.00 0.57 0.68
#MR2 0.57 1.00 0.50
#MR3 0.68 0.50 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.082, RMSR=.03, TLI=.88
#     MR1  MR2  MR3
#MR1 1.00 0.58 0.70
#MR2 0.58 1.00 0.52
#MR3 0.70 0.52 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~ RAPS03C+PI04C+Fr09C+HF02S+PI08C+RAPS07C+RAPS05C+HF08S+PI05C+PI06S+
           HF05S+RAPS04C+RAPS06C+Fr02C+HF10S+HF12S+RAPS10C+Fr01C+RAPS11C+Fr08C+
           RAPS08C+PI06C+HF01S+PI04S+Fr04C+RAPS12C+PI01S+PI02S+PI01C+PI07C  
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.960       0.862
#Tucker-Lewis Index (TLI)                       0.957       0.851
#Robust Comparative Fit Index (CFI)                         0.675
#Robust Tucker-Lewis Index (TLI)                            0.651
#RMSEA                                          0.154       0.130
#Robust RMSEA                                               0.144
#SRMR                                           0.110       0.110

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .545

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.712       0.711
#Tucker-Lewis Index (TLI)                       0.691       0.690
#Robust Comparative Fit Index (CFI)                         0.716
#Robust Tucker-Lewis Index (TLI)                            0.695
#RMSEA                                          0.119       0.110
#Robust RMSEA                                               0.117
#SRMR                                           0.095       0.095

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .439

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 3 factors based on a priori scales 
EGAmodel= '
 mistakes  =~  PI04C+Fr09C+PI08C+PI05C+Fr02C+Fr01C+Fr08C+PI06C+Fr04C+PI01C+PI07C 
 excell    =~  HF02S+HF08S+HF05S+PI06S+HF10S+HF12S+HF01S+PI04S+PI01S+PI02S            
 discrep   =~  RAPS03C+RAPS07C+RAPS05C+RAPS04C+RAPS06C+RAPS10C+RAPS11C+RAPS08C+RAPS12C            
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.996       0.968
#Tucker-Lewis Index (TLI)                       0.995       0.965
#Robust Comparative Fit Index (CFI)                         0.892
#Robust Tucker-Lewis Index (TLI)                            0.883
#RMSEA                                          0.050       0.063
#Robust RMSEA                                               0.083
#SRMR                                           0.049       0.049

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .660

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#mistakes ~~                                                           
#        excell            0.667    0.031   21.527    0.000    0.667    0.667
#        discrep           0.769    0.021   36.879    0.000    0.769    0.769
#excell ~~                                                             
#        discrep           0.609    0.034   17.930    0.000    0.609    0.609

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.926       0.929
#Tucker-Lewis Index (TLI)                       0.920       0.924
#Robust Comparative Fit Index (CFI)                         0.931
#Robust Tucker-Lewis Index (TLI)                            0.926
#RMSEA                                          0.060       0.055
#Robust RMSEA                                               0.058
#SRMR                                           0.049       0.049

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#mistakes ~~                                                           
#        excell            0.667    0.031   21.527    0.000    0.667    0.667
#        discrep           0.769    0.021   36.879    0.000    0.769    0.769
#excell ~~                                                             
#        discrep           0.609    0.034   17.930    0.000    0.609    0.609

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .579

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.867       0.868
#Tucker-Lewis Index (TLI)                       0.857       0.858
#Robust Comparative Fit Index (CFI)                         0.872
#Robust Tucker-Lewis Index (TLI)                            0.862
#RMSEA                                          0.081       0.074
#Robust RMSEA                                               0.079
#SRMR                                           0.305       0.305

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .579

library(semPlot)
semPaths(CFA_model3,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# Bifactor model
BIFmodel= '
 mistakes  =~  PI04C+Fr09C+PI08C+PI05C+Fr02C+Fr01C+Fr08C+PI06C+Fr04C+PI01C+PI07C 
 excell    =~  HF02S+HF08S+HF05S+PI06S+HF10S+HF12S+HF01S+PI04S+PI01S+PI02S            
 discrep   =~  RAPS03C+RAPS07C+RAPS05C+RAPS04C+RAPS06C+RAPS10C+RAPS11C+RAPS08C+RAPS12C            
 general=~ RAPS03C+PI04C+Fr09C+HF02S+PI08C+RAPS07C+RAPS05C+HF08S+PI05C+PI06S+
           HF05S+RAPS04C+RAPS06C+Fr02C+HF10S+HF12S+RAPS10C+Fr01C+RAPS11C+Fr08C+
           RAPS08C+PI06C+HF01S+PI04S+Fr04C+RAPS12C+PI01S+PI02S+PI01C+PI07C  
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.940       0.940
#Tucker-Lewis Index (TLI)                       0.930       0.930
#Robust Comparative Fit Index (CFI)                         0.943
#Robust Tucker-Lewis Index (TLI)                            0.934
#RMSEA                                          0.056       0.052
#Robust RMSEA                                               0.054
#SRMR                                           0.043       0.043

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .592

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.6625351      0.6873563      0.9675948      0.8357996 
#
#$FactorLevelIndices
#            ECV_SS     ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#mistakes 0.1243448 0.04667721 0.8756552 0.9360634 0.07493199 0.4814188 0.7381559
#excell   0.5153370 0.16808853 0.4846630 0.9261514 0.47887641 0.8061398 0.9051094
#discrep  0.4111312 0.12269912 0.5888688 0.9212158 0.37224868 0.7426339 0.8787210
#general  0.6625351 0.66253515 0.6625351 0.9675948 0.83579963 0.9559302 0.9599780








################################################################
### Kaçar-Başaran et al https://link.springer.com/article/10.1007/s12144-020-01131-2
### Big Three Perfectionism Scale (BTPS) 

library(haven)
Perfectionism_Kaçar_Başaran <- read_sav("Perfectionism_Kaçar-Başaran.sav")
colnames(Perfectionism_Kaçar_Başaran)
mydata <- as.data.frame(Perfectionism_Kaçar_Başaran[,2:46])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 7 factors and 4 components
# Eigenvalue 1 = 16.45; eigenvalue 2 = 3.4

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components
# Eigenvalue 1 = 18.57; eigenvalue 2 = 3.69

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.37, RMSEA=.102, RMSR=.1, TLI=.624

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.128, RMSR=.11, TLI=.559

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- EGA(mydata, algorithm="louvain")
EGALV
# suggests 6 communities no response bias as there are no reverse items

# Give solution with 3 factors
fit4 <- fa(mydata,3)
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.074, RMSR=.05, TLI=.8
#     MR1  MR3  MR2
#MR1 1.00 0.58 0.38
#MR3 0.58 1.00 0.08
#MR2 0.38 0.08 1.00

fit4 <- fa(rho,3,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.102, RMSR=.05, TLI=.719
#     MR1  MR3  MR2
#MR1 1.00 0.60 0.36
#MR3 0.60 1.00 0.08
#MR2 0.36 0.08 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  btps1+btps2+btps3+btps4+btps5+btps6+btps7+btps8+btps9+btps10+
            btps11+btps12+btps13+btps14+btps15+btps16+btps17+btps18+btps19+btps20+
            btps21+btps22+btps23+btps24+btps25+btps26+btps27+btps28+btps29+btps30+
            btps31+btps32+btps33+btps34+btps35+btps36+btps37+btps38+btps39+btps40+
            btps41+btps42+btps43+btps44+btps45
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.942       0.761
#Tucker-Lewis Index (TLI)                       0.939       0.749
#Robust Comparative Fit Index (CFI)                         0.588
#Robust Tucker-Lewis Index (TLI)                            0.569
#RMSEA                                          0.152       0.124
#Robust RMSEA                                               0.131
#SRMR                                           0.114       0.114

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .465

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.640       0.643
#Tucker-Lewis Index (TLI)                       0.622       0.626
#Robust Comparative Fit Index (CFI)                         0.646
#Robust Tucker-Lewis Index (TLI)                            0.629
#RMSEA                                          0.105       0.097
#Robust RMSEA                                               0.103
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .357

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 3 factors based on theoretical analysis
EGAmodel= '
 selfcrit=~ btps1+btps2+btps3+btps4+btps5+btps6+btps7+btps8+btps9+btps10+
            btps11+btps12+btps13+btps14+btps15+btps16+btps17+btps18
 rigid  =~  btps19+btps20+btps21+btps22+btps23+btps24+btps25+btps26+btps27+btps28
 narcist =~ btps29+btps30+btps31+btps32+btps33+btps34+btps35+btps36+btps37+btps38+
            btps39+btps40+btps41+btps42+btps43+btps44+btps45
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.942       0.761
#Tucker-Lewis Index (TLI)                       0.939       0.749
#Robust Comparative Fit Index (CFI)                         0.588
#Robust Tucker-Lewis Index (TLI)                            0.569
#RMSEA                                          0.152       0.124
#Robust RMSEA                                               0.131
#SRMR                                           0.114       0.114

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  selfcrit ~~                                                           
#    rigid             1.001    0.007  146.507    0.000    1.001    1.001
#    narcist           0.907    0.011   85.426    0.000    0.907    0.907
#  rigid ~~                                                              
#    narcist           0.914    0.010   96.095    0.000    0.914    0.914

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .397

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.654       0.657
#Tucker-Lewis Index (TLI)                       0.636       0.640
#Robust Comparative Fit Index (CFI)                         0.661
#Robust Tucker-Lewis Index (TLI)                            0.643
#RMSEA                                          0.103       0.095
#Robust RMSEA                                               0.101
#SRMR                                           0.098       0.098

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#selfcrit ~~                                                           
#        rigid             0.982    0.010   99.560    0.000    0.982    0.982
#        narcist           0.895    0.018   50.039    0.000    0.895    0.895
#rigid ~~                                                              
#        narcist           0.907    0.017   53.841    0.000    0.907    0.907

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .397

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.554       0.556
#Tucker-Lewis Index (TLI)                       0.533       0.535
#Robust Comparative Fit Index (CFI)                         0.559
#Robust Tucker-Lewis Index (TLI)                            0.538
#RMSEA                                          0.117       0.108
#Robust RMSEA                                               0.115
#SRMR                                           0.302       0.302

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .396

# Bifactor model
BIFmodel= '
 selfcrit=~ btps1+btps2+btps3+btps4+btps5+btps6+btps7+btps8+btps9+btps10+
            btps11+btps12+btps13+btps14+btps15+btps16+btps17+btps18
 rigid  =~  btps19+btps20+btps21+btps22+btps23+btps24+btps25+btps26+btps27+btps28
 narcist =~ btps29+btps30+btps31+btps32+btps33+btps34+btps35+btps36+btps37+btps38+
            btps39+btps40+btps41+btps42+btps43+btps44+btps45
 general=~  btps1+btps2+btps3+btps4+btps5+btps6+btps7+btps8+btps9+btps10+
            btps11+btps12+btps13+btps14+btps15+btps16+btps17+btps18+btps19+btps20+
            btps21+btps22+btps23+btps24+btps25+btps26+btps27+btps28+btps29+btps30+
            btps31+btps32+btps33+btps34+btps35+btps36+btps37+btps38+btps39+btps40+
            btps41+btps42+btps43+btps44+btps45
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.741       0.741
#Tucker-Lewis Index (TLI)                       0.715       0.715
#Robust Comparative Fit Index (CFI)                         0.747
#Robust Tucker-Lewis Index (TLI)                            0.722
#RMSEA                                          0.091       0.085
#Robust RMSEA                                               0.089
#SRMR                                           0.091       0.091

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .472

#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  selfcrit =~                                                           
#    btps1             0.110    0.042    2.639    0.008    0.110    0.110
#    btps2            -0.404    0.035  -11.556    0.000   -0.404   -0.404
#    btps3             0.580    0.038   15.441    0.000    0.580    0.580
#    btps4             0.325    0.039    8.417    0.000    0.325    0.325
#    btps5             0.057    0.040    1.423    0.155    0.057    0.057
#    btps6            -0.059    0.040   -1.485    0.137   -0.059   -0.059
#    btps7             0.616    0.033   18.558    0.000    0.616    0.616
#    btps8            -0.150    0.040   -3.726    0.000   -0.150   -0.150
#    btps9             0.093    0.039    2.369    0.018    0.093    0.093
#    btps10            0.025    0.034    0.728    0.467    0.025    0.025
#    btps11           -0.409    0.036  -11.298    0.000   -0.409   -0.409
#    btps12            0.190    0.034    5.663    0.000    0.190    0.190
#    btps13            0.168    0.040    4.198    0.000    0.168    0.168
#    btps14            0.394    0.037   10.624    0.000    0.394    0.394
#    btps15            0.468    0.037   12.667    0.000    0.468    0.468
#    btps16            0.093    0.037    2.511    0.012    0.093    0.093
#    btps17           -0.156    0.035   -4.405    0.000   -0.156   -0.156
#    btps18           -0.094    0.037   -2.516    0.012   -0.094   -0.094
#  rigid =~                                                              
#    btps19           -0.034    0.030   -1.142    0.253   -0.034   -0.034
#    btps20           -0.191    0.030   -6.312    0.000   -0.191   -0.191
#    btps21           -0.505    0.051   -9.984    0.000   -0.505   -0.505
#    btps22            0.513    0.033   15.519    0.000    0.513    0.513
#    btps23           -0.072    0.042   -1.722    0.085   -0.072   -0.072
#    btps24            0.106    0.037    2.852    0.004    0.106    0.106
#    btps25            0.658    0.039   16.859    0.000    0.658    0.658
#    btps26           -0.034    0.030   -1.120    0.263   -0.034   -0.034
#    btps27           -0.229    0.037   -6.175    0.000   -0.229   -0.229
#    btps28           -0.342    0.048   -7.147    0.000   -0.342   -0.342
#  narcist =~                                                            
#    btps29            0.137    0.034    4.078    0.000    0.137    0.137
#    btps30            0.398    0.040    9.921    0.000    0.398    0.398
#    btps31            0.335    0.039    8.621    0.000    0.335    0.335
#    btps32           -0.321    0.041   -7.764    0.000   -0.321   -0.321
#    btps33            0.263    0.043    6.130    0.000    0.263    0.263
#    btps34            0.443    0.038   11.656    0.000    0.443    0.443
#    btps35           -0.059    0.040   -1.483    0.138   -0.059   -0.059
#    btps36            0.488    0.037   13.021    0.000    0.488    0.488
#    btps37            0.570    0.030   19.219    0.000    0.570    0.570
#    btps38            0.047    0.032    1.446    0.148    0.047    0.047
#    btps39            0.520    0.032   16.465    0.000    0.520    0.520
#    btps40            0.216    0.033    6.478    0.000    0.216    0.216
#    btps41            0.479    0.035   13.640    0.000    0.479    0.479
#    btps42            0.277    0.039    7.090    0.000    0.277    0.277
#    btps43            0.525    0.035   14.798    0.000    0.525    0.525
#    btps44            0.033    0.032    1.031    0.302    0.033    0.033
#    btps45            0.477    0.035   13.767    0.000    0.477    0.477
#  general =~                                                            
#    btps1             0.645    0.029   21.899    0.000    0.645    0.645
#    btps2             0.665    0.025   26.740    0.000    0.665    0.665
#    btps3             0.402    0.041    9.743    0.000    0.402    0.402
#    btps4             0.325    0.043    7.630    0.000    0.325    0.325
#    btps5             0.677    0.027   25.392    0.000    0.677    0.677
#    btps6             0.760    0.021   35.651    0.000    0.760    0.760
#    btps7             0.430    0.039   11.103    0.000    0.430    0.430
#    btps8             0.769    0.021   36.873    0.000    0.769    0.769
#    btps9             0.689    0.025   27.477    0.000    0.689    0.689
#    btps10            0.834    0.016   52.492    0.000    0.834    0.834
#    btps11            0.695    0.024   29.415    0.000    0.695    0.695
#    btps12            0.680    0.027   24.940    0.000    0.680    0.680
#    btps13            0.598    0.032   18.440    0.000    0.598    0.598
#    btps14            0.561    0.033   16.741    0.000    0.561    0.561
#    btps15            0.433    0.038   11.494    0.000    0.433    0.433
#    btps16            0.740    0.024   30.505    0.000    0.740    0.740
#    btps17            0.780    0.019   41.470    0.000    0.780    0.780
#    btps18            0.741    0.022   33.800    0.000    0.741    0.741
#    btps19            0.839    0.016   53.118    0.000    0.839    0.839
#    btps20            0.841    0.017   50.579    0.000    0.841    0.841
#    btps21            0.434    0.039   11.174    0.000    0.434    0.434
#    btps22            0.734    0.022   33.301    0.000    0.734    0.734
#    btps23            0.760    0.022   33.864    0.000    0.760    0.760
#    btps24            0.777    0.020   39.537    0.000    0.777    0.777
#    btps25            0.657    0.027   24.442    0.000    0.657    0.657
#    btps26            0.863    0.014   61.868    0.000    0.863    0.863
#    btps27            0.751    0.024   31.557    0.000    0.751    0.751
#    btps28            0.484    0.037   12.932    0.000    0.484    0.484
#    btps29            0.776    0.020   38.930    0.000    0.776    0.776
#    btps30            0.495    0.039   12.773    0.000    0.495    0.495
#    btps31            0.478    0.036   13.310    0.000    0.478    0.478
#    btps32            0.755    0.021   36.817    0.000    0.755    0.755
#    btps33            0.331    0.041    8.117    0.000    0.331    0.331
#    btps34            0.459    0.036   12.646    0.000    0.459    0.459
#    btps35            0.610    0.029   20.977    0.000    0.610    0.610
#    btps36            0.425    0.040   10.759    0.000    0.425    0.425
#    btps37            0.631    0.030   21.292    0.000    0.631    0.631
#    btps38            0.733    0.023   31.712    0.000    0.733    0.733
#    btps39            0.647    0.029   22.399    0.000    0.647    0.647
#    btps40            0.750    0.023   32.924    0.000    0.750    0.750
#    btps41            0.525    0.034   15.256    0.000    0.525    0.525
#    btps42            0.599    0.031   19.208    0.000    0.599    0.599
#    btps43            0.573    0.033   17.438    0.000    0.573    0.573
#    btps44            0.771    0.020   38.061    0.000    0.771    0.771
#    btps45            0.514    0.035   14.575    0.000    0.514    0.514

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.7861986      0.6626263      0.9715186      0.9450354 

#$FactorLevelIndices
#         ECV_SS     ECV_SG    ECV_GS     Omega       OmegaH         H        FD
#selfcrit 0.1809374 0.06903644 0.8190626 0.9282449 0.0207872181 0.6930672 0.8835765
#rigid    0.1817352 0.04822408 0.8182648 0.9168740 0.0003876667 0.6307287 0.9169595
#narcist  0.2734114 0.09654091 0.7265886 0.9232472 0.1732125454 0.7527126 0.9109435
#general  0.7861986 0.78619857 0.7861986 0.9715186 0.9450353650 0.9767692 0.9890114







