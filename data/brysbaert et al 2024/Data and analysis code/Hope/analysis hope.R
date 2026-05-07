################################################################
### Hope scale Snyder, Harris, Anderson, Holleran, Irvin, Sigmon, et al., 1991
### https://fetzer.org/sites/default/files/images/stories/pdf/selfmeasures/PURPOSE_MEANING-AdultHopeScale.pdf

### Boutilier https://psyarxiv.com/cmtz3/download?format=pdf
### https://osf.io/28b9q

library(haven)
Boutilier_RRB_Hope_Optimism_Data <- read_sav("Boutilier_RRB_Hope_Optimism_Data.sav")
colnames(Boutilier_RRB_Hope_Optimism_Data)
mydata <- as.data.frame(Boutilier_RRB_Hope_Optimism_Data[,89:100])
colnames(mydata)
mydata <- na.omit(mydata)
min(mydata)
max(mydata) # 8 response alternatives


library(psych)
fit1 <- fa.parallel(mydata)
plot(fit1,show.legend=TRUE,fa="both")
fit1
# Suggests 3 factors and 2 components
# Eigenvalue 1 = 4.71; eigenvalue 2 = .69

rho <- polychoric(mydata)$rho
fit2 <- fa.parallel(rho,n.obs=nrow(mydata))
plot(fit2,show.legend=TRUE,fa="both")
fit2
# Suggests 2 components
# Eigenvalue 1 = 4.94; eigenvalue 2 = .76

# Give solution with 1 factor
fit3 <- fa(mydata,1)
fit3
diagram(fit3)
# %variance explained=.39, RMSEA=.121, RMSR=.09, TLI=.806

# Give solution with 1 factor
fit4 <- fa(rho,1,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.41, RMSEA=.14, RMSR=.09, TLI=.77

library(EGAnet); library(foreign); library(ggplot2)
EGALV <- riEGA(mydata, algorithm="louvain")
EGALV
# suggests 2 communities and response bias

# Give solution with 2 factors
fit4 <- fa(mydata,2)
fit4
diagram(fit4)
# %variance explained=.47, RMSEA=.093, RMSR=.06, TLI=.884
#     MR1  MR2
#MR1  1.0 -0.4
#MR2 -0.4  1.0

fit4 <- fa(rho,2,n.obs=nrow(mydata))
fit4
diagram(fit4)
# %variance explained=.50, RMSEA=.118, RMSR=.06, TLI=.838
#      MR1   MR2
#MR1  1.00 -0.35
#MR2 -0.35  1.00

print(colnames(mydata))

# Single factor model lavaan
UNImodel= '
 general=~  Hope1+Hope2+Hope3+Hope4+Hope5+Hope6+Hope7+Hope8+Hope9+Hope10+Hope11+Hope12
'

library(lavaan)
CFA_model1 <- cfa(UNImodel, data = mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.971       0.874
#Tucker-Lewis Index (TLI)                       0.965       0.846
#Robust Comparative Fit Index (CFI)                         0.806
#Robust Tucker-Lewis Index (TLI)                            0.763
#RMSEA                                          0.128       0.165
#Robust RMSEA                                               0.147
#SRMR                                           0.088       0.088

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .565

CFA_model1 <- cfa(UNImodel, data = mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model1, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.841       0.839
#Tucker-Lewis Index (TLI)                       0.806       0.803
#Robust Comparative Fit Index (CFI)                         0.846
#Robust Tucker-Lewis Index (TLI)                            0.812
#RMSEA                                          0.123       0.113
#Robust RMSEA                                               0.120
#SRMR                                           0.080       0.080

fitRsquares <- lavInspect(CFA_model1, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .501

library(semPlot)
semPaths(CFA_model1,"std",layout="tree2",residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

# CFA analysis model with 2 factors based on theoretical analysis
EGAmodel= '
 agency  =~ Hope2+Hope3+Hope7+Hope9+Hope10+Hope12
 pathway =~ Hope1+Hope4+Hope5+Hope6+Hope8+Hope11
'

library(lavaan)
CFA_model2 <- cfa(EGAmodel,mydata,ordered=TRUE,std.lv=TRUE)
# warning
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.980       0.900
#Tucker-Lewis Index (TLI)                       0.974       0.875
#Robust Comparative Fit Index (CFI)                         0.865
#Robust Tucker-Lewis Index (TLI)                            0.831
#RMSEA                                          0.110       0.148
#Robust RMSEA                                               0.124
#SRMR                                           0.080       0.080

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#agency ~~                                                             
# pathway           0.830    0.027   31.142    0.000    0.830    0.830

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .610

CFA_model2 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE)
summary(CFA_model2,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.900       0.900
#Tucker-Lewis Index (TLI)                       0.876       0.876
#Robust Comparative Fit Index (CFI)                         0.905
#Robust Tucker-Lewis Index (TLI)                            0.882
#RMSEA                                          0.098       0.090
#Robust RMSEA                                               0.095
#SRMR                                           0.074       0.074

#Covariances:
#Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#agency ~~                                                             
# pathway           0.813    0.046   17.676    0.000    0.813    0.813

fitRsquares <- lavInspect(CFA_model2, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .573

library(semPlot)
semPaths(CFA_model2,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2)

CFA_model3 <- cfa(EGAmodel,mydata,estimator="MLR",std.lv=TRUE,orthogonal=TRUE)
summary(CFA_model3,standardized=TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.772       0.765
#Tucker-Lewis Index (TLI)                       0.722       0.713
#Robust Comparative Fit Index (CFI)                         0.776
#Robust Tucker-Lewis Index (TLI)                            0.726
#RMSEA                                          0.147       0.136
#Robust RMSEA                                               0.145
#SRMR                                           0.252       0.252

fitRsquares <- lavInspect(CFA_model3, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .561

# Bifactor model
BIFmodel= '
 agency  =~ Hope2+Hope3+Hope7+Hope9+Hope10+Hope12
 pathway =~ Hope1+Hope4+Hope5+Hope6+Hope8+Hope11
 general=~  Hope1+Hope2+Hope3+Hope4+Hope5+Hope6+Hope7+Hope8+Hope9+Hope10+Hope11+Hope12
'

options(max.print=999999)
CFA_model4 <- cfa(BIFmodel, data = mydata,std.lv=TRUE, estimator="MLR",orthogonal=TRUE)
summary(CFA_model4, standardized = TRUE, fit.measures = TRUE)
#Comparative Fit Index (CFI)                    0.933       0.920
#Tucker-Lewis Index (TLI)                       0.894       0.875
#Robust Comparative Fit Index (CFI)                         0.932
#Robust Tucker-Lewis Index (TLI)                            0.894
#RMSEA                                          0.090       0.090
#Robust RMSEA                                               0.090
#SRMR                                           0.064       0.064

fitRsquares <- lavInspect(CFA_model4, what='rsquare')
summary(fitRsquares)[-4]
# median %variance explained = .61

# loadings based on Std.all because Std.lv is larger than 1
#Latent Variables:
#                   Estimate  Std.Err  z-value  P(>|z|)   Std.lv  Std.all
#  agency =~                                                             
#    Hope2             0.485    0.171    2.840    0.005    0.485    0.290
#    Hope3            -0.493    0.237   -2.081    0.037   -0.493   -0.270
#    Hope7            -0.165    0.243   -0.680    0.497   -0.165   -0.093
#    Hope9             0.566    0.198    2.853    0.004    0.566    0.301
#    Hope10            0.983    0.180    5.463    0.000    0.983    0.533
#    Hope12            0.862    0.147    5.882    0.000    0.862    0.552
#  pathway =~                                                            
#    Hope1             0.459    0.202    2.277    0.023    0.459    0.283
#    Hope4             0.687    0.550    1.248    0.212    0.687    0.472
#    Hope5            -0.271    0.406   -0.668    0.504   -0.271   -0.154
#    Hope6            -0.075    0.520   -0.144    0.886   -0.075   -0.048
#    Hope8             0.336    0.184    1.820    0.069    0.336    0.223
#    Hope11            0.120    0.318    0.379    0.705    0.120    0.068
#  general =~                                                            
#    Hope1             1.305    0.111   11.777    0.000    1.305    0.804
#    Hope2             1.281    0.120   10.656    0.000    1.281    0.767
#    Hope3            -0.534    0.210   -2.540    0.011   -0.534   -0.293
#    Hope4             1.027    0.176    5.847    0.000    1.027    0.706
#    Hope5            -0.538    0.172   -3.130    0.002   -0.538   -0.305
#    Hope6             1.251    0.150    8.330    0.000    1.251    0.808
#    Hope7             0.091    0.182    0.500    0.617    0.091    0.051
#    Hope8             1.115    0.116    9.586    0.000    1.115    0.740
#    Hope9             1.265    0.133    9.503    0.000    1.265    0.673
#    Hope10            1.149    0.139    8.273    0.000    1.149    0.622
#    Hope11           -0.581    0.195   -2.979    0.003   -0.581   -0.327
#    Hope12            0.881    0.126    7.003    0.000    0.881    0.564

library(semPlot)
semPaths(CFA_model4,"std", layout="tree2", residuals=FALSE, thresholds=FALSE,sizeMan=5,sizeLat=,edge.label.cex=1.2,
         label.cex=2,bifactor="general")

library(BifactorIndicesCalculator)
bifactorIndices(CFA_model4)
#Ideally there wouldn't be negative loadings
#$ModelLevelIndices
#ECV.general            PUC  Omega.general OmegaH.general 
#0.7811049      0.5454545      0.8000663      0.7239035 

#$FactorLevelIndices
#         ECV_SS    ECV_SG    ECV_GS     Omega     OmegaH         H        FD
#agency  0.3151914 0.1506039 0.6848086 0.6906482 0.16061599 0.5267019 0.7619593
#pathway 0.1307803 0.0682912 0.8692197 0.6822326 0.07355667 0.3136090 0.6506467
#general 0.7811049 0.7811049 0.7811049 0.8000663 0.72390351 0.9055466 0.9371256





