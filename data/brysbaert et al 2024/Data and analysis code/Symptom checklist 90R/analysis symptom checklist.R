######################################################################
### Analysis symptom checklist to test the usefulness of the criteria
###

################################################################
### Chinese version Chen et al. (2020) https://www.sciencedirect.com/science/article/pii/S2352340920300962
### data https://www.sciencedirect.com/science/article/pii/S2352340920300962

SCL90R_Chen <- read.csv("SCL90R_Chen.csv")
colnames(SCL90R_Chen)
mydata <- SCL90R_Chen[,55:144]
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 19 factors and 7 components
# Eigenvalue 1 = 29.59; eigenvalue 2 = 1.87

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 19 factors and 7 components
# Eigenvalue 1 = 38.19; eigenvalue 2 = 2.18

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.33, RMSEA=.049, RMSR=.04, TLI=.769

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.42, RMSEA=.066, RMSR=.05, TLI=.724

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 7 communities

# Give solution with 9 factors (scale-based)
fit4 <- fa(mydata,9)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.033, RMSR=.02, TLI=.899
#     MR2  MR3  MR1  MR7  MR4  MR6  MR8  MR5  MR9
#MR2 1.00 0.52 0.43 0.40 0.48 0.38 0.30 0.34 0.27
#MR3 0.52 1.00 0.46 0.50 0.48 0.38 0.38 0.29 0.32
#MR1 0.43 0.46 1.00 0.53 0.50 0.34 0.44 0.36 0.32
#MR7 0.40 0.50 0.53 1.00 0.45 0.40 0.43 0.34 0.31
#MR4 0.48 0.48 0.50 0.45 1.00 0.41 0.38 0.33 0.24
#MR6 0.38 0.38 0.34 0.40 0.41 1.00 0.32 0.27 0.29
#MR8 0.30 0.38 0.44 0.43 0.38 0.32 1.00 0.19 0.30
#MR5 0.34 0.29 0.36 0.34 0.33 0.27 0.19 1.00 0.23
#MR9 0.27 0.32 0.32 0.31 0.24 0.29 0.30 0.23 1.00

fit4 <- fa(rho,9,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.54, RMSEA=.048, RMSR=.02, TLI=.855
#     MR2  MR1  MR3  MR7  MR5  MR8  MR6  MR4  MR9
#MR2 1.00 0.44 0.54 0.45 0.48 0.43 0.35 0.32 0.27
#MR1 0.44 1.00 0.47 0.58 0.48 0.51 0.39 0.28 0.32
#MR3 0.54 0.47 1.00 0.52 0.44 0.42 0.37 0.26 0.27
#MR7 0.45 0.58 0.52 1.00 0.48 0.46 0.44 0.31 0.35
#MR5 0.48 0.48 0.44 0.48 1.00 0.38 0.31 0.29 0.33
#MR8 0.43 0.51 0.42 0.46 0.38 1.00 0.37 0.27 0.24
#MR6 0.35 0.39 0.37 0.44 0.31 0.37 1.00 0.19 0.21
#MR4 0.32 0.28 0.26 0.31 0.29 0.27 0.19 1.00 0.08
#MR9 0.27 0.32 0.27 0.35 0.33 0.24 0.21 0.08 1.00


print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  SCL_b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11+b12+b13+b14+b15+b16+b17+b18+b19+b20+
            b21+b22+b23+b24+b25+b26+b27+b28+b29+b30+b31+b32+b33+b34+b35+b36+b37+b38+b39+b40+
            b41+b42+b43+b44+b45+b46+b47+b48+b49+b50+b51+b52+b53+b54+b55+b56+b57+b58+b59+b60+
            b61+b62+b63+b64+b65+b66+b67+b68+b69+b70+b71+b72+b73+b74+b75+b76+b77+b78+b79+b80+
            b81+b82+b83+b84+b85+b86+b87+b88+b89+b90
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.989       0.917
#Tucker-Lewis Index (TLI)                       0.989       0.915
#Robust Comparative Fit Index (CFI)                         0.751
#Robust Tucker-Lewis Index (TLI)                            0.745
#RMSEA                                          0.043       0.041
#Robust RMSEA                                               0.064
#SRMR                                           0.046       0.046

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .446

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.774       0.790
#Tucker-Lewis Index (TLI)                       0.769       0.785
#Robust Comparative Fit Index (CFI)                         0.790
#Robust Tucker-Lewis Index (TLI)                            0.785
#RMSEA                                          0.050       0.035
#Robust RMSEA                                               0.047
#SRMR                                           0.040       0.040

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .321

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

colnames(mydata)

# CFA analysis model with 9 factors based on the classification of the authors
EGAmodel= '
 somatization   =~ SCL_b1+b4+b12+b27+b40+b42+b48+b49+b52+b53+b56+b58
 obs_compuls    =~ b3+b9+b10+b28+b38+b45+b46+b51+b55+b65
 interpers_sens =~ b6+b21+b34+b36+b37+b41+b61+b69+b73
 depression     =~ b5+b14+b15+b20+b22+b26+b29+b30+b31+b32+b54+b71+b79
 anxiety        =~ b2+b17+b23+b33+b39+b57+b72+b78+b80+b86
 hostility      =~ b11+b24+b63+b67+b74+b81
 phobic_anx     =~ b13+b25+b47+b50+b70+b75+b82
 paranoid       =~ b8+b18+b43+b68+b76+b83
 psychoticism   =~ b7+b16+b35+b62+b77+b84+b85+b87+b88+b90
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.992       0.939
#Tucker-Lewis Index (TLI)                       0.992       0.937
#Robust Comparative Fit Index (CFI)                         0.811
#Robust Tucker-Lewis Index (TLI)                            0.804
#RMSEA                                          0.037       0.037
#Robust RMSEA                                               0.058
#SRMR                                           0.041       0.041

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .394

#Covariances:
#                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  somatization ~~                                                        
#    obs_compuls        0.820    0.008   98.887    0.000    0.820    0.820
#    interpers_sens     0.769    0.010   78.489    0.000    0.769    0.769
#    depression         0.788    0.009   89.916    0.000    0.788    0.788
#    anxiety            0.882    0.006  136.361    0.000    0.882    0.882
#    hostility          0.748    0.011   67.483    0.000    0.748    0.748
#    phobic_anx         0.777    0.012   66.022    0.000    0.777    0.777
#    paranoid           0.817    0.010   78.839    0.000    0.817    0.817
#    psychoticism       0.849    0.009   95.903    0.000    0.849    0.849
#  obs_compuls ~~                                                         
#    interpers_sens     0.917    0.006  162.638    0.000    0.917    0.917
#    depression         0.922    0.004  209.956    0.000    0.922    0.922
#    anxiety            0.921    0.005  172.683    0.000    0.921    0.921
#    hostility          0.823    0.008   99.432    0.000    0.823    0.823
#    phobic_anx         0.805    0.010   76.848    0.000    0.805    0.805
#    paranoid           0.916    0.007  131.412    0.000    0.916    0.916
#    psychoticism       0.930    0.006  161.420    0.000    0.930    0.930
#  interpers_sens ~~                                                      
#    depression         0.930    0.004  223.691    0.000    0.930    0.930
#    anxiety            0.909    0.005  178.135    0.000    0.909    0.909
#    hostility          0.828    0.008  100.257    0.000    0.828    0.828
#    phobic_anx         0.830    0.009   89.021    0.000    0.830    0.830
#    paranoid           0.983    0.005  188.508    0.000    0.983    0.983
#    psychoticism       0.950    0.005  180.883    0.000    0.950    0.950
#  depression ~~                                                          
#    anxiety            0.927    0.004  234.302    0.000    0.927    0.927
#    hostility          0.816    0.008   99.597    0.000    0.816    0.816
#    phobic_anx         0.802    0.010   82.726    0.000    0.802    0.802
#    paranoid           0.903    0.006  143.073    0.000    0.903    0.903
#    psychoticism       0.940    0.005  203.615    0.000    0.940    0.940
#  anxiety ~~                                                             
#    hostility          0.845    0.008  104.866    0.000    0.845    0.845
#    phobic_anx         0.882    0.008  107.725    0.000    0.882    0.882
#    paranoid           0.913    0.007  134.241    0.000    0.913    0.913
#    psychoticism       0.948    0.005  188.168    0.000    0.948    0.948
#  hostility ~~                                                           
#    phobic_anx         0.722    0.013   55.701    0.000    0.722    0.722
#    paranoid           0.891    0.008  108.042    0.000    0.891    0.891
#    psychoticism       0.858    0.008  105.120    0.000    0.858    0.858
#  phobic_anx ~~                                                          
#    paranoid           0.819    0.012   69.368    0.000    0.819    0.819
#    psychoticism       0.825    0.011   76.069    0.000    0.825    0.825
#  paranoid ~~                                                            
#    psychoticism       0.988    0.006  174.642    0.000    0.988    0.988

# average correlation = .867


CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.828       0.842
#Tucker-Lewis Index (TLI)                       0.821       0.836
#Robust Comparative Fit Index (CFI)                         0.843
#Robust Tucker-Lewis Index (TLI)                            0.837
#RMSEA                                          0.046       0.032
#Robust RMSEA                                               0.043
#SRMR                                           0.038       0.038

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .396

#Covariances:
#                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  somatization ~~                                                        
#    obs_compuls        0.797    0.011   69.365    0.000    0.797    0.797
#    interpers_sens     0.755    0.014   55.191    0.000    0.755    0.755
#    depression         0.766    0.013   58.460    0.000    0.766    0.766
#    anxiety            0.860    0.010   82.160    0.000    0.860    0.860
#    hostility          0.734    0.015   47.915    0.000    0.734    0.734
#    phobic_anx         0.761    0.017   44.027    0.000    0.761    0.761
#    paranoid           0.801    0.015   55.173    0.000    0.801    0.801
#    psychoticism       0.832    0.014   59.825    0.000    0.832    0.832
#  obs_compuls ~~                                                         
#    interpers_sens     0.920    0.009  105.800    0.000    0.920    0.920
#    depression         0.929    0.007  137.402    0.000    0.929    0.929
#    anxiety            0.906    0.008  110.599    0.000    0.906    0.906
#    hostility          0.809    0.012   66.458    0.000    0.809    0.809
#    phobic_anx         0.792    0.015   54.324    0.000    0.792    0.792
#    paranoid           0.903    0.011   84.766    0.000    0.903    0.903
#    psychoticism       0.918    0.008  111.425    0.000    0.918    0.918
#  interpers_sens ~~                                                      
#    depression         0.951    0.006  149.017    0.000    0.951    0.951
#    anxiety            0.917    0.008  115.798    0.000    0.917    0.917
#    hostility          0.829    0.013   63.995    0.000    0.829    0.829
#    phobic_anx         0.829    0.014   58.425    0.000    0.829    0.829
#    paranoid           0.993    0.008  124.493    0.000    0.993    0.993
#    psychoticism       0.956    0.008  112.578    0.000    0.956    0.956
#  depression ~~                                                          
#    anxiety            0.940    0.006  148.797    0.000    0.940    0.940
#    hostility          0.819    0.013   64.876    0.000    0.819    0.819
#    phobic_anx         0.803    0.014   56.347    0.000    0.803    0.803
#    paranoid           0.909    0.009   96.867    0.000    0.909    0.909
#    psychoticism       0.947    0.007  130.094    0.000    0.947    0.947
#  anxiety ~~                                                             
#    hostility          0.838    0.013   66.628    0.000    0.838    0.838
#    phobic_anx         0.880    0.012   72.569    0.000    0.880    0.880
#    paranoid           0.907    0.010   89.294    0.000    0.907    0.907
#    psychoticism       0.939    0.008  115.045    0.000    0.939    0.939
#  hostility ~~                                                           
#    phobic_anx         0.719    0.018   38.901    0.000    0.719    0.719
#    paranoid           0.890    0.012   71.215    0.000    0.890    0.890
#    psychoticism       0.857    0.012   71.225    0.000    0.857    0.857
#  phobic_anx ~~                                                          
#    paranoid           0.810    0.017   47.489    0.000    0.810    0.810
#    psychoticism       0.809    0.016   49.674    0.000    0.809    0.809
#  paranoid ~~                                                            
#    psychoticism       0.979    0.009  113.851    0.000    0.979    0.979

# average correlation = .861

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.596       0.606
#Tucker-Lewis Index (TLI)                       0.586       0.596
#Robust Comparative Fit Index (CFI)                         0.607
#Robust Tucker-Lewis Index (TLI)                            0.597
#RMSEA                                          0.069       0.050
#Robust RMSEA                                               0.068
#SRMR                                           0.308       0.308

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .388


# Bifactor model
BIFmodel= '
 somatization   =~ SCL_b1+b4+b12+b27+b40+b42+b48+b49+b52+b53+b56+b58
 obs_compuls    =~ b3+b9+b10+b28+b38+b45+b46+b51+b55+b65
 interpers_sens =~ b6+b21+b34+b36+b37+b41+b61+b69+b73
 depression     =~ b5+b14+b15+b20+b22+b26+b29+b30+b31+b32+b54+b71+b79
 anxiety        =~ b2+b17+b23+b33+b39+b57+b72+b78+b80+b86
 hostility      =~ b11+b24+b63+b67+b74+b81
 phobic_anx     =~ b13+b25+b47+b50+b70+b75+b82
 paranoid       =~ b8+b18+b43+b68+b76+b83
 psychoticism   =~ b7+b16+b35+b62+b77+b84+b85+b87+b88+b90
 general=~  SCL_b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11+b12+b13+b14+b15+b16+b17+b18+b19+b20+
            b21+b22+b23+b24+b25+b26+b27+b28+b29+b30+b31+b32+b33+b34+b35+b36+b37+b38+b39+b40+
            b41+b42+b43+b44+b45+b46+b47+b48+b49+b50+b51+b52+b53+b54+b55+b56+b57+b58+b59+b60+
            b61+b62+b63+b64+b65+b66+b67+b68+b69+b70+b71+b72+b73+b74+b75+b76+b77+b78+b79+b80+
            b81+b82+b83+b84+b85+b86+b87+b88+b89+b90
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.841       0.856
#Tucker-Lewis Index (TLI)                       0.834       0.849
#Robust Comparative Fit Index (CFI)                         0.857
#Robust Tucker-Lewis Index (TLI)                            0.851
#RMSEA                                          0.042       0.029
#Robust RMSEA                                               0.039
#SRMR                                           0.035       0.035

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .402

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.8282928      0.9116105      0.9793524      0.9659481 
#
#$FactorLevelIndices
#                   ECV_SS      ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#somatization   0.30831718 0.038795786 0.6916828 0.8738274 0.26572964 0.6126491 0.8164745
#obs_compuls    0.15974939 0.017345822 0.8402506 0.8481831 0.07786503 0.4419949 0.7125567
#interpers_sens 0.16853218 0.018134013 0.8314678 0.8516213 0.05787477 0.4761648 0.8114609
#depression     0.12342155 0.019505159 0.8765785 0.9019914 0.03018794 0.5021772 0.9236971
#anxiety        0.14012114 0.018398434 0.8598789 0.8866425 0.05460216 0.4500753 0.7817370
#hostility      0.27424085 0.020940038 0.7257592 0.8211666 0.19088109 0.4871148 0.7848366
#phobic_anx     0.34143908 0.023069781 0.6585609 0.7685568 0.23842234 0.5024151 0.7489161
#paranoid       0.10388674 0.006563520 0.8961133 0.7736272 0.06328581 0.2025254 0.5091696
#psychoticism   0.09290522 0.008954653 0.9070948 0.8238376 0.01344378 0.2694238 0.6058607
#general        0.82829280 0.828292795 0.8282928 0.9793524 0.96594812 0.9790764 0.9868140





################################################################
### Urban et al. (2014) https://www.sciencedirect.com/science/article/pii/S0165178114000687
### data sent by the author

d = read.fwf("SCL_90_Urban.dat", widths= c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                                           1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
mydata  <- as.data.frame(d[,2:91])
mydata <- na.omit(mydata)
colnames(mydata) <- paste0('b', 1:ncol(mydata))
colnames(mydata)
#mydata[mydata>4] <- NA
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
packageVersion("psych")

fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 16 factors and 7 components
# Eigenvalue 1 = 31.26; eigenvalue 2 = 3.14

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 16 factors and 7 components
# Eigenvalue 1 = 43.89; eigenvalue 2 = 3.54

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.35, RMSEA=.06, RMSR=.05, TLI=.707

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.49, RMSEA=.09, RMSR=.06, TLI=.634

library(EGAnet); library(foreign); library(ggplot2)
EGAL <- EGA(mydata)
EGAL
# suggests 7 communities

# Give solution with 9 factors (scale-based)
fit4 <- fa(mydata,9)
fit4
diagram(fit4)
# %variance explained=.48, RMSEA=.037, RMSR=.02, TLI=.885
#     MR2  MR3  MR7  MR8  MR4  MR5  MR1  MR6  MR9
#MR2 1.00 0.49 0.32 0.47 0.37 0.46 0.40 0.16 0.13
#MR3 0.49 1.00 0.34 0.42 0.36 0.35 0.30 0.40 0.19
#MR7 0.32 0.34 1.00 0.41 0.44 0.42 0.34 0.30 0.32
#MR8 0.47 0.42 0.41 1.00 0.34 0.38 0.39 0.26 0.27
#MR4 0.37 0.36 0.44 0.34 1.00 0.35 0.27 0.38 0.22
#MR5 0.46 0.35 0.42 0.38 0.35 1.00 0.43 0.21 0.13
#MR1 0.40 0.30 0.34 0.39 0.27 0.43 1.00 0.11 0.17
#MR6 0.16 0.40 0.30 0.26 0.38 0.21 0.11 1.00 0.23
#MR9 0.13 0.19 0.32 0.27 0.22 0.13 0.17 0.23 1.00
#
# Average absolute correlation = .320
fit4 <- fa(rho,9,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.63, RMSEA=.068, RMSR=.02, TLI=.787
#     MR2  MR3  MR4  MR5  MR7  MR1  MR8  MR9  MR6
#MR2 1.00 0.49 0.38 0.38 0.31 0.36 0.30 0.55 0.38
#MR3 0.49 1.00 0.38 0.43 0.35 0.36 0.47 0.46 0.26
#MR4 0.38 0.38 1.00 0.47 0.46 0.37 0.47 0.32 0.24
#MR5 0.38 0.43 0.47 1.00 0.38 0.35 0.39 0.37 0.32
#MR7 0.31 0.35 0.46 0.38 1.00 0.47 0.40 0.26 0.31
#MR1 0.36 0.36 0.37 0.35 0.47 1.00 0.40 0.34 0.34
#MR8 0.30 0.47 0.47 0.39 0.40 0.40 1.00 0.31 0.14
#MR9 0.55 0.46 0.32 0.37 0.26 0.34 0.31 1.00 0.35
#MR6 0.38 0.26 0.24 0.32 0.31 0.34 0.14 0.35 1.00
#
# Average absolute correlation = .370

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11+b12+b13+b14+b15+b16+b17+b18+b19+b20+
            b21+b22+b23+b24+b25+b26+b27+b28+b29+b30+b31+b32+b33+b34+b35+b36+b37+b38+b39+b40+
            b41+b42+b43+b44+b45+b46+b47+b48+b49+b50+b51+b52+b53+b54+b55+b56+b57+b58+b59+b60+
            b61+b62+b63+b64+b65+b66+b67+b68+b69+b70+b71+b72+b73+b74+b75+b76+b77+b78+b79+b80+
            b81+b82+b83+b84+b85+b86+b87+b88+b89+b90
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.986       0.902
#Tucker-Lewis Index (TLI)                       0.985       0.900
#Robust Comparative Fit Index (CFI)                         0.683
#Robust Tucker-Lewis Index (TLI)                            0.676
#RMSEA                                          0.051       0.043
#Robust RMSEA                                               0.085
#SRMR                                           0.061       0.061

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .507

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.713       0.736
#Tucker-Lewis Index (TLI)                       0.706       0.730
#Robust Comparative Fit Index (CFI)                         0.737
#Robust Tucker-Lewis Index (TLI)                            0.731
#RMSEA                                          0.060       0.040
#Robust RMSEA                                               0.057
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .352

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


colnames(mydata)

# CFA analysis model with 9 factors 
EGAmodel= '
 somatization   =~ b1+b4+b12+b27+b40+b42+b48+b49+b52+b53+b56+b58
 obs_compuls    =~ b3+b9+b10+b28+b38+b45+b46+b51+b55+b65
 interpers_sens =~ b6+b21+b34+b36+b37+b41+b61+b69+b73
 depression     =~ b5+b14+b15+b20+b22+b26+b29+b30+b31+b32+b54+b71+b79
 anxiety        =~ b2+b17+b23+b33+b39+b57+b72+b78+b80+b86
 hostility      =~ b11+b24+b63+b67+b74+b81
 phobic_anx     =~ b13+b25+b47+b50+b70+b75+b82
 paranoid       =~ b8+b18+b43+b68+b76+b83
 psychoticism   =~ b7+b16+b35+b62+b77+b84+b85+b87+b88+b90
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.991       0.932
#Tucker-Lewis Index (TLI)                       0.991       0.929
#Robust Comparative Fit Index (CFI)                         0.770
#Robust Tucker-Lewis Index (TLI)                            0.762
#RMSEA                                          0.042       0.038
#Robust RMSEA                                               0.075
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .579

#Covariances:
#                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  somatization ~~                                                        
#    obs_compuls        0.829    0.010   82.194    0.000    0.829    0.829
#    interpers_sens     0.732    0.014   52.369    0.000    0.732    0.732
#    depression         0.836    0.009   91.487    0.000    0.836    0.836
#    anxiety            0.923    0.006  145.451    0.000    0.923    0.923
#    hostility          0.672    0.017   40.400    0.000    0.672    0.672
#    phobic_anx         0.809    0.013   62.044    0.000    0.809    0.809
#    paranoid           0.691    0.017   41.740    0.000    0.691    0.691
#    psychoticism       0.761    0.013   58.688    0.000    0.761    0.761
#  obs_compuls ~~                                                         
#    interpers_sens     0.929    0.007  141.553    0.000    0.929    0.929
#    depression         0.953    0.005  209.134    0.000    0.953    0.953
#    anxiety            0.956    0.005  174.066    0.000    0.956    0.956
#    hostility          0.838    0.012   71.003    0.000    0.838    0.838
#    phobic_anx         0.859    0.011   75.619    0.000    0.859    0.859
#    paranoid           0.892    0.010   92.725    0.000    0.892    0.892
#    psychoticism       0.920    0.008  121.833    0.000    0.920    0.920
#  interpers_sens ~~                                                      
#    depression         0.920    0.007  139.365    0.000    0.920    0.920
#    anxiety            0.921    0.007  129.856    0.000    0.921    0.921
#    hostility          0.840    0.012   70.684    0.000    0.840    0.840
#    phobic_anx         0.832    0.012   68.429    0.000    0.832    0.832
#    paranoid           0.954    0.007  140.305    0.000    0.954    0.954
#    psychoticism       0.927    0.008  123.366    0.000    0.927    0.927
#  depression ~~                                                          
#    anxiety            0.969    0.004  238.087    0.000    0.969    0.969
#    hostility          0.802    0.012   65.214    0.000    0.802    0.802
#    phobic_anx         0.835    0.011   75.095    0.000    0.835    0.835
#    paranoid           0.868    0.009   91.427    0.000    0.868    0.868
#    psychoticism       0.873    0.009   99.989    0.000    0.873    0.873
#  anxiety ~~                                                             
#    hostility          0.869    0.011   80.847    0.000    0.869    0.869
#    phobic_anx         0.927    0.008  122.416    0.000    0.927    0.927
#    paranoid           0.892    0.010   85.540    0.000    0.892    0.892
#    psychoticism       0.924    0.007  129.846    0.000    0.924    0.924
#  hostility ~~                                                           
#    phobic_anx         0.762    0.017   45.099    0.000    0.762    0.762
#    paranoid           0.863    0.012   71.761    0.000    0.863    0.863
#    psychoticism       0.851    0.012   71.877    0.000    0.851    0.851
#  phobic_anx ~~                                                          
#    paranoid           0.770    0.016   47.372    0.000    0.770    0.770
#    psychoticism       0.881    0.011   80.871    0.000    0.881    0.881
#  paranoid ~~                                                            
#    psychoticism       0.902    0.010   90.307    0.000    0.902    0.902
#
# Average absolute correlation = .861

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.796       0.819
#Tucker-Lewis Index (TLI)                       0.789       0.813
#Robust Comparative Fit Index (CFI)                         0.821
#Robust Tucker-Lewis Index (TLI)                            0.815
#RMSEA                                          0.053       0.034
#Robust RMSEA                                               0.049
#SRMR                                           0.054       0.054

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .415

#Covariances:
#                    Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  somatization ~~                                                        
#    obs_compuls        0.810    0.015   52.423    0.000    0.810    0.810
#    interpers_sens     0.705    0.020   35.471    0.000    0.705    0.705
#    depression         0.810    0.014   56.373    0.000    0.810    0.810
#    anxiety            0.884    0.012   72.737    0.000    0.884    0.884
#    hostility          0.612    0.026   23.547    0.000    0.612    0.612
#    phobic_anx         0.744    0.024   30.974    0.000    0.744    0.744
#    paranoid           0.648    0.024   27.406    0.000    0.648    0.648
#    psychoticism       0.681    0.020   33.461    0.000    0.681    0.681
#  obs_compuls ~~                                                         
#    interpers_sens     0.928    0.010   95.324    0.000    0.928    0.928
#    depression         0.957    0.007  134.729    0.000    0.957    0.957
#    anxiety            0.952    0.009  102.232    0.000    0.952    0.952
#    hostility          0.797    0.022   36.893    0.000    0.797    0.797
#    phobic_anx         0.788    0.023   34.386    0.000    0.788    0.788
#    paranoid           0.873    0.015   56.929    0.000    0.873    0.873
#    psychoticism       0.883    0.013   68.875    0.000    0.883    0.883
#  interpers_sens ~~                                                      
#    depression         0.913    0.010   87.486    0.000    0.913    0.913
#    anxiety            0.901    0.012   75.179    0.000    0.901    0.901
#    hostility          0.804    0.021   38.623    0.000    0.804    0.804
#    phobic_anx         0.746    0.024   31.167    0.000    0.746    0.746
#    paranoid           0.942    0.011   85.791    0.000    0.942    0.942
#    psychoticism       0.896    0.013   67.415    0.000    0.896    0.896
#  depression ~~                                                          
#    anxiety            0.961    0.007  144.187    0.000    0.961    0.961
#    hostility          0.751    0.023   32.320    0.000    0.751    0.751
#    phobic_anx         0.752    0.021   35.118    0.000    0.752    0.752
#    paranoid           0.847    0.015   56.244    0.000    0.847    0.847
#    psychoticism       0.821    0.016   52.618    0.000    0.821    0.821
#  anxiety ~~                                                             
#    hostility          0.815    0.022   37.432    0.000    0.815    0.815
#    phobic_anx         0.888    0.017   53.378    0.000    0.888    0.888
#    paranoid           0.844    0.019   44.943    0.000    0.844    0.844
#    psychoticism       0.875    0.015   58.430    0.000    0.875    0.875
#  hostility ~~                                                           
#    phobic_anx         0.664    0.032   20.536    0.000    0.664    0.664
#    paranoid           0.829    0.020   40.765    0.000    0.829    0.829
#    psychoticism       0.825    0.020   41.390    0.000    0.825    0.825
#  phobic_anx ~~                                                          
#    paranoid           0.651    0.027   23.861    0.000    0.651    0.651
#    psychoticism       0.805    0.027   29.695    0.000    0.805    0.805
#  paranoid ~~                                                            
#    psychoticism       0.856    0.017   50.712    0.000    0.856    0.856
#
# Average absolute correlation = .818

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# MLR model with orthogonal factors
CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.599       0.614
#Tucker-Lewis Index (TLI)                       0.588       0.604
#Robust Comparative Fit Index (CFI)                         0.617
#Robust Tucker-Lewis Index (TLI)                            0.608
#RMSEA                                          0.074       0.050
#Robust RMSEA                                               0.071
#SRMR                                           0.327       0.327

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .415


# Bifactor model
BIFmodel= '
 somatization   =~ b1+b4+b12+b27+b40+b42+b48+b49+b52+b53+b56+b58
 obs_compuls    =~ b3+b9+b10+b28+b38+b45+b46+b51+b55+b65
 interpers_sens =~ b6+b21+b34+b36+b37+b41+b61+b69+b73
 depression     =~ b5+b14+b15+b20+b22+b26+b29+b30+b31+b32+b54+b71+b79
 anxiety        =~ b2+b17+b23+b33+b39+b57+b72+b78+b80+b86
 hostility      =~ b11+b24+b63+b67+b74+b81
 phobic_anx     =~ b13+b25+b47+b50+b70+b75+b82
 paranoid       =~ b8+b18+b43+b68+b76+b83
 psychoticism   =~ b7+b16+b35+b62+b77+b84+b85+b87+b88+b90
 general=~  b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11+b12+b13+b14+b15+b16+b17+b18+b19+b20+
            b21+b22+b23+b24+b25+b26+b27+b28+b29+b30+b31+b32+b33+b34+b35+b36+b37+b38+b39+b40+
            b41+b42+b43+b44+b45+b46+b47+b48+b49+b50+b51+b52+b53+b54+b55+b56+b57+b58+b59+b60+
            b61+b62+b63+b64+b65+b66+b67+b68+b69+b70+b71+b72+b73+b74+b75+b76+b77+b78+b79+b80+
            b81+b82+b83+b84+b85+b86+b87+b88+b89+b90
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.797       0.819
#Tucker-Lewis Index (TLI)                       0.787       0.811
#Robust Comparative Fit Index (CFI)                         0.822
#Robust Tucker-Lewis Index (TLI)                            0.814
#RMSEA                                          0.051       0.033
#Robust RMSEA                                               0.047
#SRMR                                           0.047       0.047

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .421

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.8036089      0.9086142      0.9817527      0.9642132 
#
#$FactorLevelIndices
#                   ECV_SS     ECV_SG    ECV_GS     Omega      OmegaH         H        FD
#somatization   0.32484084 0.04711496 0.6751592 0.9086284 0.283524322 0.6928541 0.8740455
#obs_compuls    0.09272652 0.01011732 0.9072735 0.8724902 0.043587672 0.3121182 0.6346119
#interpers_sens 0.15054435 0.01459025 0.8494557 0.8536262 0.072621742 0.4263220 0.7751628
#depression     0.13273938 0.02058846 0.8672606 0.9094776 0.055394242 0.5680344 0.9388797
#anxiety        0.09175332 0.01105237 0.9082467 0.8867066 0.007754996 0.3120059 0.6598366
#hostility      0.37667833 0.02742056 0.6233217 0.8237648 0.260417320 0.6109936 0.8376919
#phobic_anx     0.33339743 0.02963054 0.6666026 0.8596001 0.241308192 0.6095213 0.8432700
#paranoid       0.21915624 0.01375244 0.7808438 0.7950532 0.155187988 0.3778031 0.6834587
#psychoticism   0.22368745 0.02212416 0.7763126 0.8518470 0.161004935 0.4989634 0.7565797
#general        0.80360895 0.80360895 0.8036089 0.9817527 0.964213198 0.9801315 0.9867234






