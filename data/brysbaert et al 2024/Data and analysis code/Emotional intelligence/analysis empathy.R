################################################################
### Emotion perception
### Data: Interpersonal Reactivity Index Briganti
### https://donaldrwilliams.github.io/BGGM/reference/iri.html#examples

library(BGGM)
data(iri)
head(iri)
mydata <- iri[,1:28]
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 6 factors 5 components
# Eigenvalue 1 = 3.79; eigenvalue 2 = 1.93

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components
# Eigenvalue 1 = 4.44; eigenvalue 2 = 2.28

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.14, RMSEA=.095, RMSR=.1, TLI=.63

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.16, RMSEA=.113, RMSR=.12, TLI=.365

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.31, RMSEA=.048, RMSR=.03, TLI=.845
#     MR3  MR1  MR2  MR4
#MR3 1.00 0.23 0.21 0.17
#MR1 0.23 1.00 0.19 0.23
#MR2 0.21 0.19 1.00 0.00
#MR4 0.17 0.23 0.00 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.36, RMSEA=.061, RMSR=.04, TLI=.816
#     MR3  MR1   MR2   MR4
#MR3 1.00 0.24  0.20  0.18
#MR1 0.24 1.00  0.16  0.23
#MR2 0.20 0.16  1.00 -0.02
#MR4 0.18 0.23 -0.02  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  davi1+davi2+davi3+davi4+davi5+davi6+davi7+davi8+davi9+davi10+davi11+davi12+davi13+davi14+
             davi15+davi16+davi17+davi18+davi19+davi20+davi21+davi22+davi23+davi24+davi25+davi26+
             davi27+davi28
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.655       0.456
#Tucker-Lewis Index (TLI)                       0.628       0.412
#Robust Comparative Fit Index (CFI)                         0.407
#Robust Tucker-Lewis Index (TLI)                            0.360
#RMSEA                                          0.141       0.131
#Robust RMSEA                                               0.114
#SRMR                                           0.117       0.117

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .179

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.438       0.431
#Tucker-Lewis Index (TLI)                       0.393       0.385
#Robust Comparative Fit Index (CFI)                         0.441
#Robust Tucker-Lewis Index (TLI)                            0.396
#RMSEA                                          0.095       0.085
#Robust RMSEA                                               0.095
#SRMR                                           0.098       0.098

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .127

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors 
EGAmodel= '
 fantasy =~  davi1+davi5+davi7+davi12+davi16+davi23+davi26
 perspect=~  davi3+davi8+davi11+davi15+davi21+davi25+davi28
 empathy =~  davi2+davi4+davi9+davi14+davi18+davi20+davi22
 distress=~  davi6+davi10+davi13+davi17+davi19+davi24+davi27
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.888       0.817
#Tucker-Lewis Index (TLI)                       0.877       0.799
#Robust Comparative Fit Index (CFI)                         0.758
#Robust Tucker-Lewis Index (TLI)                            0.734
#RMSEA                                          0.081       0.076
#Robust RMSEA                                               0.073
#SRMR                                           0.073       0.073

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  fantasy ~~                                                            
#    perspect          0.198    0.025    7.818    0.000    0.198    0.198
#    empathy           0.401    0.022   18.124    0.000    0.401    0.401
#    distress          0.307    0.022   13.910    0.000    0.307    0.307
#  perspect ~~                                                           
#    empathy           0.494    0.023   21.927    0.000    0.494    0.494
#    distress         -0.074    0.026   -2.883    0.004   -0.074   -0.074
#  empathy ~~                                                            
#    distress          0.348    0.023   15.174    0.000    0.348    0.348

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .352

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.789       0.790
#Tucker-Lewis Index (TLI)                       0.768       0.769
#Robust Comparative Fit Index (CFI)                         0.794
#Robust Tucker-Lewis Index (TLI)                            0.774
#RMSEA                                          0.059       0.052
#Robust RMSEA                                               0.058
#SRMR                                           0.062       0.062

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .293

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  fantasy ~~                                                            
#    perspect          0.208    0.034    6.136    0.000    0.208    0.208
#    empathy           0.371    0.031   11.965    0.000    0.371    0.371
#    distress          0.304    0.034    9.033    0.000    0.304    0.304
#  perspect ~~                                                           
#    empathy           0.476    0.032   14.969    0.000    0.476    0.476
#    distress         -0.059    0.043   -1.360    0.174   -0.059   -0.059
#  empathy ~~                                                            
#    distress          0.361    0.039    9.296    0.000    0.361    0.361

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.730       0.730
#Tucker-Lewis Index (TLI)                       0.709       0.709
#Robust Comparative Fit Index (CFI)                         0.735
#Robust Tucker-Lewis Index (TLI)                            0.714
#RMSEA                                          0.066       0.059
#Robust RMSEA                                               0.065
#SRMR                                           0.101       0.101

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .279


# Bifactor model
BIFmodel= '
 fantasy =~  davi1+davi5+davi7+davi12+davi16+davi23+davi26
 perspect=~  davi3+davi8+davi11+davi15+davi21+davi25+davi28
 empathy =~  davi2+davi4+davi9+davi14+davi18+davi20+davi22
 distress=~  davi6+davi10+davi13+davi17+davi19+davi24+davi27
 general =~  davi1+davi2+davi3+davi4+davi5+davi6+davi7+davi8+davi9+davi10+davi11+davi12+davi13+davi14+
             davi15+davi16+davi17+davi18+davi19+davi20+davi21+davi22+davi23+davi24+davi25+davi26+
             davi27+davi28
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.862       0.863
#Tucker-Lewis Index (TLI)                       0.838       0.839
#Robust Comparative Fit Index (CFI)                         0.867
#Robust Tucker-Lewis Index (TLI)                            0.844
#RMSEA                                          0.049       0.044
#Robust RMSEA                                               0.048
#SRMR                                           0.058       0.058

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .324

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.2829155      0.7777778      0.8213444      0.4278784 
#
#$FactorLevelIndices
#            ECV_SS    ECV_SG    ECV_GS     Omega    OmegaH         H        FD
#fantasy  0.7436516 0.1904492 0.2563484 0.7499708 0.5648689 0.7265509 0.8438162
#perspect 0.8148380 0.1549603 0.1851620 0.6703606 0.5928836 0.6394185 0.7973384
#empathy  0.6370328 0.1744549 0.3629672 0.7692859 0.5052744 0.7053925 0.8257448
#distress 0.7046810 0.1972201 0.2953190 0.7576033 0.6378327 0.7174746 0.8461552
#general  0.2829155 0.2829155 0.2829155 0.8213444 0.4278784 0.7567983 0.8395357








################################################################
### Emotion perception
### Variant of Interpersonal Reactivity Index
### Bainbridge https://osf.io/f9hmg/?view_only=  dataset s2

Bainbridge_s2_IRI <- read.csv("Bainbridge_s2_IRI.csv")
colnames(Bainbridge_s2_IRI)
mydata <- as.data.frame(Bainbridge_s2_IRI[,224:239])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 3 components
# Eigenvalue 1 = 3.54; eigenvalue 2 = 1.12

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 3 components
# Eigenvalue 1 = 4.3; eigenvalue 2 = 1.32

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.22, RMSEA=.125, RMSR=.12, TLI=.52

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.27, RMSEA=.159, RMSR=.13, TLI=.479

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 3 communities if response bias

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.042, RMSR=.03, TLI=.946
#     MR3  MR1  MR4  MR2
#MR3 1.00 0.30 0.39 0.13
#MR1 0.30 1.00 0.49 0.10
#MR4 0.39 0.49 1.00 0.10
#MR2 0.13 0.10 0.10 1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.52, RMSEA=.076, RMSR=.03, TLI=.879
#     MR3  MR1  MR4  MR2
#MR3 1.00 0.33 0.43 0.12
#MR1 0.33 1.00 0.52 0.10
#MR4 0.43 0.52 1.00 0.09
#MR2 0.12 0.10 0.09 1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  empec_1+empec_2+empec_3+empec_4+
             empf_1+empf_2+empf_3+empf_4+
             emppd_1+emppd_2+emppd_3+emppd_4+
             emppt_1+emppt_2+emppt_3+emppt_4
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.824       0.662
#Tucker-Lewis Index (TLI)                       0.797       0.610
#Robust Comparative Fit Index (CFI)                         0.551
#Robust Tucker-Lewis Index (TLI)                            0.482
#RMSEA                                          0.169       0.169
#Robust RMSEA                                               0.160
#SRMR                                           0.129       0.129

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .357

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.584       0.558
#Tucker-Lewis Index (TLI)                       0.521       0.491
#Robust Comparative Fit Index (CFI)                         0.587
#Robust Tucker-Lewis Index (TLI)                            0.523
#RMSEA                                          0.127       0.118
#Robust RMSEA                                               0.125
#SRMR                                           0.108       0.108

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .263

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors 
EGAmodel= '
 concern =~  empec_1+empec_2+empec_3+empec_4
 fantasy =~  empf_1+empf_2+empf_3+empf_4
 distress=~  emppd_1+emppd_2+emppd_3+emppd_4
 perspect=~  emppt_1+emppt_2+emppt_3+emppt_4
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.983       0.956
#Tucker-Lewis Index (TLI)                       0.979       0.946
#Robust Comparative Fit Index (CFI)                         0.900
#Robust Tucker-Lewis Index (TLI)                            0.878
#RMSEA                                          0.054       0.063
#Robust RMSEA                                               0.078
#SRMR                                           0.059       0.059

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  concern ~~                                                            
#    fantasy           0.525    0.043   12.092    0.000    0.525    0.525
#    distress          0.124    0.056    2.208    0.027    0.124    0.124
#    perspect          0.679    0.036   19.071    0.000    0.679    0.679
#  fantasy ~~                                                            
#    distress          0.194    0.054    3.602    0.000    0.194    0.194
#    perspect          0.416    0.049    8.544    0.000    0.416    0.416
#  distress ~~                                                           
#    perspect          0.112    0.062    1.812    0.070    0.112    0.112

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .498

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.934       0.938
#Tucker-Lewis Index (TLI)                       0.919       0.925
#Robust Comparative Fit Index (CFI)                         0.942
#Robust Tucker-Lewis Index (TLI)                            0.929
#RMSEA                                          0.052       0.045
#Robust RMSEA                                               0.048
#SRMR                                           0.052       0.052

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .441

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  concern ~~                                                            
#    fantasy           0.492    0.063    7.762    0.000    0.492    0.492
#    distress          0.122    0.077    1.580    0.114    0.122    0.122
#    perspect          0.650    0.061   10.654    0.000    0.650    0.650
#  fantasy ~~                                                            
#    distress          0.177    0.069    2.553    0.011    0.177    0.177
#    perspect          0.377    0.071    5.304    0.000    0.377    0.377
#  distress ~~                                                           
#    perspect          0.090    0.079    1.147    0.251    0.090    0.090

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.826       0.826
#Tucker-Lewis Index (TLI)                       0.799       0.799
#Robust Comparative Fit Index (CFI)                         0.833
#Robust Tucker-Lewis Index (TLI)                            0.808
#RMSEA                                          0.082       0.074
#Robust RMSEA                                               0.080
#SRMR                                           0.152       0.152

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .431


# Bifactor model
BIFmodel= '
 concern =~  empec_1+empec_2+empec_3+empec_4
 fantasy =~  empf_1+empf_2+empf_3+empf_4
 distress=~  emppd_1+emppd_2+emppd_3+emppd_4
 perspect=~  emppt_1+emppt_2+emppt_3+emppt_4
 general =~  empec_1+empec_2+empec_3+empec_4+
             empf_1+empf_2+empf_3+empf_4+
             emppd_1+emppd_2+emppd_3+emppd_4+
             emppt_1+emppt_2+emppt_3+emppt_4
'
options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.970       0.973
#Tucker-Lewis Index (TLI)                       0.959       0.963
#Robust Comparative Fit Index (CFI)                         0.975
#Robust Tucker-Lewis Index (TLI)                            0.966
#RMSEA                                          0.037       0.032
#Robust RMSEA                                               0.033
#SRMR                                           0.040       0.040

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
#     0.3551897      0.8000000      0.8804781      0.6214389 
#
#$FactorLevelIndices
#            ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#concern  0.5655910 0.2162105 0.43440896 0.9129827 0.2992132 1.8005646 1.6570995
#fantasy  0.7081861 0.1656396 0.29181394 0.7881396 0.5560420 0.7031356 0.8375937
#distress 0.9482166 0.1620291 0.05178343 0.6875462 0.6558181 0.6838369 0.8270832
#perspect 0.4739528 0.1009310 0.52604724 0.7531779 0.3390067 0.5683952 0.7591687
#general  0.3551897 0.3551897 0.35518972 0.8804781 0.6214389 0.8082967 0.8649965








################################################################
### Emotion perception
### Data: Interpersonal Reactivity Index
### Martingano & Konrath, 2021 https://osf.io/mc8sj?view_only=c6b8c4647da042efab016435b29d2125

Martingano_IRI <- read_sav("Martingano_IRI.sav")
colnames(Martingano_IRI)
mydata <- as.data.frame(Martingano_IRI[,187:214])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 5 factors and 5 components
# Eigenvalue 1 = 5.93; eigenvalue 2 = 2.68

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 5 components
# Eigenvalue 1 = 7.04; eigenvalue 2 = 3.02

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.21, RMSEA=.119, RMSR=.13, TLI=.441

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.25, RMSEA=.144, RMSR=.15, TLI=.411

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities if response bias

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.43, RMSEA=.057, RMSR=.04, TLI=.87
#      MR1   MR3   MR4   MR2
#MR1  1.00 -0.39 -0.35  0.05
#MR3 -0.39  1.00  0.25  0.06
#MR4 -0.35  0.25  1.00 -0.11
#MR2  0.05  0.06 -0.11  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.082, RMSR=.05, TLI=.81
#      MR4   MR3   MR1   MR2
#MR4  1.00 -0.41 -0.36  0.05
#MR3 -0.41  1.00  0.27  0.05
#MR1 -0.36  0.27  1.00 -0.12
#MR2  0.05  0.05 -0.12  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  iri1+iri2+iri3+iri4+iri5+iri6+iri7+iri8+iri9+iri10+iri11+iri12+iri13+iri14+
             iri15+iri16+iri17+iri18+iri19+iri20+iri21+iri22+iri23+iri24+iri25+iri26+
             iri27+iri28
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.758       0.537
#Tucker-Lewis Index (TLI)                       0.739       0.500
#Robust Comparative Fit Index (CFI)                         0.457
#Robust Tucker-Lewis Index (TLI)                            0.414
#RMSEA                                          0.195       0.153
#Robust RMSEA                                               0.147
#SRMR                                           0.146       0.146

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .312

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.482       0.469
#Tucker-Lewis Index (TLI)                       0.440       0.426
#Robust Comparative Fit Index (CFI)                         0.486
#Robust Tucker-Lewis Index (TLI)                            0.445
#RMSEA                                          0.121       0.111
#Robust RMSEA                                               0.120
#SRMR                                           0.128       0.128

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .224

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors 
EGAmodel= '
 fantasy =~  iri1+iri5+iri7+iri12+iri16+iri23+iri26
 perspect=~  iri3+iri8+iri11+iri15+iri21+iri25+iri28
 empathy =~  iri2+iri4+iri9+iri14+iri18+iri20+iri22
 distress=~  iri6+iri10+iri13+iri17+iri19+iri24+iri27
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.951       0.883
#Tucker-Lewis Index (TLI)                       0.946       0.872
#Robust Comparative Fit Index (CFI)                         0.794
#Robust Tucker-Lewis Index (TLI)                            0.774
#RMSEA                                          0.088       0.078
#Robust RMSEA                                               0.091
#SRMR                                           0.082       0.082

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  fantasy ~~                                                            
#    perspect         -0.393    0.042   -9.349    0.000   -0.393   -0.393
#    empathy           0.606    0.036   17.027    0.000    0.606    0.606
#    distress          0.045    0.050    0.898    0.369    0.045    0.045
#  perspect ~~                                                           
#    empathy          -0.595    0.037  -16.007    0.000   -0.595   -0.595
#    distress          0.217    0.046    4.691    0.000    0.217    0.217
#  empathy ~~                                                            
#    distress          0.011    0.050    0.225    0.822    0.011    0.011

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .495

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.835       0.840
#Tucker-Lewis Index (TLI)                       0.819       0.824
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.830
#RMSEA                                          0.069       0.061
#Robust RMSEA                                               0.066
#SRMR                                           0.074       0.074

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .426

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  fantasy ~~                                                            
#    perspect         -0.371    0.064   -5.819    0.000   -0.371   -0.371
#    empathy           0.567    0.060    9.442    0.000    0.567    0.567
#    distress          0.049    0.076    0.645    0.519    0.049    0.049
#  perspect ~~                                                           
#    empathy          -0.586    0.060   -9.819    0.000   -0.586   -0.586
#    distress          0.212    0.066    3.183    0.001    0.212    0.212
#  empathy ~~                                                            
#    distress          0.003    0.080    0.041    0.967    0.003    0.003

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.777       0.781
#Tucker-Lewis Index (TLI)                       0.759       0.763
#Robust Comparative Fit Index (CFI)                         0.787
#Robust Tucker-Lewis Index (TLI)                            0.770
#RMSEA                                          0.079       0.071
#Robust RMSEA                                               0.077
#SRMR                                           0.153       0.153

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .425


# Bifactor model
BIFmodel= '
 fantasy =~  iri1+iri5+iri7+iri12+iri16+iri23+iri26
 perspect=~  iri3+iri8+iri11+iri15+iri21+iri25+iri28
 empathy =~  iri2+iri4+iri9+iri14+iri18+iri20+iri22
 distress=~  iri6+iri10+iri13+iri17+iri19+iri24+iri27
 general =~  iri1+iri2+iri3+iri4+iri5+iri6+iri7+iri8+iri9+iri10+iri11+iri12+iri13+iri14+
             iri15+iri16+iri17+iri18+iri19+iri20+iri21+iri22+iri23+iri24+iri25+iri26+
             iri27+iri28
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.870       0.872
#Tucker-Lewis Index (TLI)                       0.848       0.850
#Robust Comparative Fit Index (CFI)                         0.879
#Robust Tucker-Lewis Index (TLI)                            0.858
#RMSEA                                          0.063       0.057
#Robust RMSEA                                               0.061
#SRMR                                           0.066       0.066

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
#     0.4279538      0.7777778      0.6485029      0.3517940 
#
#$FactorLevelIndices
#            ECV_SS     ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#fantasy  0.6318399 0.15727289 0.36816012 0.5569506 0.4155351 0.7511963 0.8726826
#perspect 0.6600096 0.14996979 0.33999040 0.5678536 0.3492145 0.7171105 0.8505978
#empathy  0.1637972 0.04909555 0.83620278 0.5026393 0.3226788 0.4444891 0.7278688
#distress 0.9624235 0.21570800 0.03757651 0.4913085 0.4813506 0.8194151 0.9072503
#general  0.4279538 0.42795376 0.42795376 0.6485029 0.3517940 0.8925426 0.9340746









################################################################
### Emotion perception
### Interpersonal Reactivity Index
### Stosic et al. (2022) https://www.tandfonline.com/doi/pdf/10.1080/00224545.2021.1985417

library(readxl)
IRI_Stosic <- read_excel("IRI_Stosic.xlsx")
colnames(IRI_Stosic)
mydata <- as.data.frame(IRI_Stosic[,2:29])
mydata <- na.omit(mydata)
colnames(mydata)
min(mydata)
max(mydata) # 5 response alternatives

library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 5.68; eigenvalue 2 = 2.79

rho <- polychoric(mydata)$rho
# Does not run because no full range on all questions
library(EGAnet)
rho <- polychoric.matrix(mydata)
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 4 factors and 4 components
# Eigenvalue 1 = 6.31; eigenvalue 2 = 3.21

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.20, RMSEA=.122, RMSR=.15, TLI=.398

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.23, RMSEA=.155, RMSR=.17, TLI=.318

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 4 communities with response bias

# Give solution with 4 factors (theory-based)
fit4 <- fa(mydata,4)
fit4
diagram(fit4)
# %variance explained=.44, RMSEA=.057, RMSR=.05, TLI=.866
#      MR1   MR3   MR4   MR2
#MR1  1.00 -0.33 -0.36 -0.18
#MR3 -0.33  1.00  0.22 -0.08
#MR4 -0.36  0.22  1.00 -0.14
#MR2 -0.18 -0.08 -0.14  1.00

fit4 <- fa(rho,4,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.095, RMSR=.05, TLI=.737
#      MR3   MR1   MR4   MR2
#MR3  1.00 -0.36  0.20 -0.09
#MR1 -0.36  1.00 -0.35 -0.17
#MR4  0.20 -0.35  1.00 -0.15
#MR2 -0.09 -0.17 -0.15  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general =~  dav1+dav2+dav3+dav4+dav5+dav6+dav7+dav8+dav9+dav10+dav11+dav12+dav13+dav14+
             dav15+dav16+dav17+dav18+dav19+dav20+dav21+dav22+dav23+dav24+dav25+dav26+
             dav27+dav28
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.688       0.443
#Tucker-Lewis Index (TLI)                       0.663       0.399
#Robust Comparative Fit Index (CFI)                         0.387
#Robust Tucker-Lewis Index (TLI)                            0.338
#RMSEA                                          0.208       0.161
#Robust RMSEA                                               0.160
#SRMR                                           0.165       0.165

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .244

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.438       0.409
#Tucker-Lewis Index (TLI)                       0.393       0.361
#Robust Comparative Fit Index (CFI)                         0.436
#Robust Tucker-Lewis Index (TLI)                            0.391
#RMSEA                                          0.128       0.123
#Robust RMSEA                                               0.127
#SRMR                                           0.141       0.141

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .187

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)


# CFA analysis model with 4 factors 
EGAmodel= '
 fantasy =~  dav1+dav5+dav7+dav12+dav16+dav23+dav26
 perspect=~  dav3+dav8+dav11+dav15+dav21+dav25+dav28
 empathy =~  dav2+dav4+dav9+dav14+dav18+dav20+dav22
 distress=~  dav6+dav10+dav13+dav17+dav19+dav24+dav27
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.956       0.893
#Tucker-Lewis Index (TLI)                       0.952       0.882
#Robust Comparative Fit Index (CFI)                         0.788
#Robust Tucker-Lewis Index (TLI)                            0.767
#RMSEA                                          0.079       0.071
#Robust RMSEA                                               0.095
#SRMR                                           0.090       0.090

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  fantasy ~~                                                            
#    perspect         -0.322    0.065   -4.947    0.000   -0.322   -0.322
#    empathy           0.488    0.055    8.855    0.000    0.488    0.488
#    distress         -0.097    0.078   -1.241    0.214   -0.097   -0.097
#  perspect ~~                                                           
#    empathy          -0.495    0.063   -7.835    0.000   -0.495   -0.495
#    distress          0.178    0.070    2.535    0.011    0.178    0.178
#  empathy ~~                                                            
#    distress          0.245    0.070    3.510    0.000    0.245    0.245

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .512

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.842       0.846
#Tucker-Lewis Index (TLI)                       0.827       0.831
#Robust Comparative Fit Index (CFI)                         0.852
#Robust Tucker-Lewis Index (TLI)                            0.837
#RMSEA                                          0.068       0.063
#Robust RMSEA                                               0.065
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .424

#Covariances:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  fantasy ~~                                                            
#    perspect         -0.354    0.091   -3.872    0.000   -0.354   -0.354
#    empathy           0.503    0.076    6.635    0.000    0.503    0.503
#    distress         -0.146    0.109   -1.336    0.182   -0.146   -0.146
#  perspect ~~                                                           
#    empathy          -0.491    0.095   -5.187    0.000   -0.491   -0.491
#    distress          0.176    0.111    1.591    0.112    0.176    0.176
#  empathy ~~                                                            
#    distress          0.223    0.101    2.197    0.028    0.223    0.223

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.791       0.795
#Tucker-Lewis Index (TLI)                       0.774       0.779
#Robust Comparative Fit Index (CFI)                         0.801
#Robust Tucker-Lewis Index (TLI)                            0.785
#RMSEA                                          0.078       0.072
#Robust RMSEA                                               0.075
#SRMR                                           0.148       0.148

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .431


# Bifactor model
BIFmodel= '
 fantasy =~  dav1+dav5+dav7+dav12+dav16+dav23+dav26
 perspect=~  dav3+dav8+dav11+dav15+dav21+dav25+dav28
 empathy =~  dav2+dav4+dav9+dav14+dav18+dav20+dav22
 distress=~  dav6+dav10+dav13+dav17+dav19+dav24+dav27
 general =~  dav1+dav2+dav3+dav4+dav5+dav6+dav7+dav8+dav9+dav10+dav11+dav12+dav13+dav14+
             dav15+dav16+dav17+dav18+dav19+dav20+dav21+dav22+dav23+dav24+dav25+dav26+
             dav27+dav28
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
# warning
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.873       0.874
#Tucker-Lewis Index (TLI)                       0.851       0.852
#Robust Comparative Fit Index (CFI)                         0.880
#Robust Tucker-Lewis Index (TLI)                            0.860
#RMSEA                                          0.063       0.059
#Robust RMSEA                                               0.061
#SRMR                                           0.089       0.089

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .432

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")


library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#$ModelLevelIndices
#   ECV.general            PUC  Omega.general OmegaH.general 
#     0.3255055      0.7777778      0.5490203      0.2504243 
#
#$FactorLevelIndices
#            ECV_SS    ECV_SG     ECV_GS     Omega    OmegaH         H        FD
#fantasy  0.3970365 0.1280818 0.60296346 0.6483650 0.2822089 1.4543927 1.4178621
#perspect 0.8046790 0.1889195 0.19532103 0.5668477 0.4947822 0.8024255 0.8985895
#empathy  0.6868944 0.1659318 0.31310558 0.1578005 0.1122736 0.7617401 0.8769116
#distress 0.9527519 0.1915614 0.04724808 0.3778226 0.3766694 0.8014906 0.8972630
#general  0.3255055 0.3255055 0.32550550 0.5490203 0.2504243 0.8842926 0.9319294






