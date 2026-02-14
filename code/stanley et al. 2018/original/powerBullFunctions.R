
#Load the necessary opackages
library(readstata13)
library(metafor)
library(tidyverse)
library(reshape2)
#library(dplyr)
library(broom)

# Function to return a result data frame either in wide or long format
returnRes <- function(res, long=TRUE, reduce=TRUE) {
  if (is.null(res)) return(NULL)
  
  # convert all factor columns to characters
  res %>% mutate_if(is.factor, as.character) -> res	
  
  if (long==FALSE) {
    # return wide format
    return(res)
  } else {
    # transform to long format
    longRes <- melt(res, id.vars=c("method", "term"))
    if (reduce==TRUE & nrow(res) > 1) {longRes <- longRes %>% filter(!is.na(value)) %>% arrange(method, term, variable)}
    return(longRes)
  }
}

#Function for tidy output for metafor's RMA function
tidyRMA <- function(RMA) {
  res <- data.frame(
    term = if(length(RMA$b)==1) {"b0"} else {c("b0", "b1")},
    estimate = as.vector(RMA$b),
    std.error = RMA$se,
    statistic = RMA$zval,
    p.value = RMA$pval,
    conf.low = RMA$ci.lb,
    conf.high = RMA$ci.ub
  )
  rownames(res) <- NULL
  return(res)
}

#Wrapper function for metafor's RMA function and trim-and-fill
RMA.est <- function(d, v, long=TRUE) {
  
  #analyzes MA data set using standard RE model estimators
  #produces estimate of true effect, CI around estimate,
  #and estimate of tau (could get CIs if wanted)
  
  # adjust stepadj (make it smaller by x0.5) and increase maxiter from 100 to 500 to prevent convergence problems
  #reMA <- rma(d, v, method="REML")
  
  reMA <- tryCatch(
    rma(d, v, method="REML", control = list(stepadj = .5, maxiter=500)),
    error = function(e) {
      NULL
    }
  )
  
  # if it fails: return empty
  if (is.null(reMA)) {
    warning("RMA did not converge.")
    # initialize empty results object
    res <- data.frame(
      method = "reMA",
      term = "b0",
      estimate = NA,
      std.error = NA,
      statistic = NA,
      p.value = NA,	# one-tailed p-value of p-uniform's test of null-hypothesis of no effect
      conf.low = NA,
      conf.high = NA
    )
    return(returnRes(res, long))
  }
  
  # assign NULL to tfMA if an error is raised
  tfMA <- NULL
  try({
    tfMA <- trimfill(reMA, side="left", maxiter=500)
  }, silent=TRUE)  
  
  res <- data.frame(method="reMA", tidyRMA(reMA))
  res <- plyr::rbind.fill(res, data.frame(
    method="reMA",
    term=c("tau2", "I2", "Q"),
    estimate=c(reMA$tau2, reMA$I2, reMA$QE),
    std.error=c(reMA$se.tau2, NA, NA)
  ))		
  
  
  if (!is.null(tfMA)) {	  
    res <- rbind(res, data.frame(method="TF", tidyRMA(tfMA)))
    res <- plyr::rbind.fill(res, data.frame(
      method="TF",
      term="tau2",
      estimate=tfMA$tau2,
      std.error=tfMA$se.tau2
    ))
    res <- plyr::rbind.fill(res, data.frame(
      method="TF",
      term="kFilled",
      estimate=tfMA$k0
    ))
  }
  
  returnRes(res, long)
}


#Function to apply a WLS meta-regression
WLS.est <- function(d, v, long=TRUE) {
  se <- sqrt(v)
  d.precision <- 1/se
  d.stand <- d/se
  l1 <- lm(d.stand ~ 0 + d.precision)
  s1 <- summary(l1)
  res <- data.frame(
    method = "WLS",
    term = "b0",
    estimate = 	s1$coefficients[1, 1],
    std.error = s1$coefficients[1, 2],
    statistic = s1$coefficients[1, 2],
    p.value = s1$coefficients[1, 4],
    conf.low = confint(l1)[1],
    conf.high = confint(l1)[2]
  )
  returnRes(res, long)
}

#Function for WAAP
WAAP.est <- function(d, v, n1, n2, est=c("WAAP-WLS"), long=TRUE) {
  
  # 1. determine which studies are in the top-N set
  
  # "Here, we employ FE (or, equivalently, WLS) as the proxy for 'true' effect."
  WLS.all  <- WLS.est(d, v, long=FALSE)
  true.effect <- WLS.all$estimate
  
  # only select studies that are adequatly powered (Stanley uses a two-sided test)
  #added abs()
  powered <- abs(true.effect/2.8) >= sqrt(v)
  
  # 2. compute the unrestricted weighted average (WLS) rma of either all or only adequatly powered studies
  # "Thus, the estimate we employ here is a hybrid between WAAP and WLS; thereby called, 'WAAP-WLS'"
  # "If there is no or only one adequately powered study in a systematic review, WAAP's standard error and confidence interval are undefined. In this case, we compute WLS across the entire research record."
  # "When there are two or more adequately powered studies, WAAP is calculated using WLS's formula from this subset of adequately powered studies."
  # Combined estimator: WAAP-WLS	
  kAdequate <- sum(powered,na.rm=T)
  
  #4/3/2019
  #Modified the following so that WAAP was used even with only 1 adequately powered study
  #if (kAdequate >= 2) { changed to if (kAdequate >= 1) {
  
  if (kAdequate >= 1) {
    res <- WLS.est(d[powered], v[powered], long=FALSE)
    res$method <- "WAAP-WLS"
    res <- plyr::rbind.fill(res, data.frame(method="WAAP-WLS", term="estimator", type="WAAP", kAdequate=kAdequate))
  } else {
    res <- WLS.all
    res$method <- "WAAP-WLS"
    res <- plyr::rbind.fill(res, data.frame(method="WAAP-WLS", term="estimator", type="WLS", kAdequate=kAdequate))
  }
  returnRes(res, long)
}

#Function for tidy output of LM() function
tidyLM <- function(...) {
  res <- tidy(..., conf.int=TRUE)
  res$term[res$term == "(Intercept)"] <- "b0"
  res$term[res$term == "sqrt(v)"] <- "b1"
  res$term[res$term == "v"] <- "b1"
  return(res)
}


#Function to get estimates using PET-PEESE
PETPEESE.est <- function(d, v, PP.test = "one-sided", runRMA=FALSE, long=TRUE) {
  
  PP.test <- match.arg(PP.test, PP.test)
  
  PET.lm <- lm(d~sqrt(v), weights=1/v)
  PEESE.lm <- lm(d~v, weights=1/v)
  
  res <- rbind(
    data.frame(method="PET.lm", tidyLM(PET.lm)),
    data.frame(method="PEESE.lm", tidyLM(PEESE.lm))
  )  
  
  # conditional PET/PEESE estimator	
  lm.p.value  <- res %>% filter(method == "PET.lm", term == "b0") %>% .[["p.value"]]
  lm.est  <- res %>% filter(method == "PET.lm", term == "b0") %>% .[["estimate"]]
  
  
  # "For the purpose of deciding which meta-regression 
  #  accommodation for selective reporting bias to employ, 
  #  we recommend testing H0:b0 < 0 at the 10% significance 
  #  level." From: Stanley, T. D. (2016). Limitations of 
  #  PET-PEESE and other meta-analysis methods. Abgerufen 
  #  von https://www.hendrix.edu/uploadedFiles/Departments_and_Programs/Business_and_Economics/AMAES/Limitations%20of%20PET-PEESE.pdf
  
  
  #4/3/2019
  #Modified the conditional determining whether PET or PEESE is used.
  # lm.p.value < p.crit & lm.est > 0 became lm.p.value < p.crit & abs(lm.est) > 0
  
  if (PP.test == "one-sided") {
    p.crit <- .10
  } else if (PP.test == "two-sided") {
    p.crit <- .05
  }
  
  # the condition always uses a directional hypothesis (a reversed slope should not happen), but different critical levels (.05 vs .10)
  usePEESE.lm <- ifelse(lm.p.value < p.crit & abs(lm.est) > 0, TRUE, FALSE)
  
  res <- rbind(res,
               data.frame(method="PETPEESE.lm", if (usePEESE.lm == TRUE) {tidyLM(PEESE.lm)} else {tidyLM(PET.lm)})
  )
  
  
  if (runRMA == TRUE) {
    PET.rma <- tryCatch(
      rma(yi = d, vi = v, mods=sqrt(v), method="REML", control = list(stepadj = .5, maxiter=500)),
      error = function(e) rma(yi = d, vi = v, mods=sqrt(v), method="DL")
    )
    
    PEESE.rma <- tryCatch(
      rma(yi = d, vi = v, mods=v, method="REML", control = list(stepadj = .5, maxiter=500)),
      error = function(e) rma(yi = d, vi = v, mods=v, method="DL")
    )  
    
    res <- rbind(
      res,
      data.frame(method="PET.rma", tidyRMA(PET.rma)),
      data.frame(method="PEESE.rma", tidyRMA(PEESE.rma))
    )
    
    rma.p.value <- res %>% filter(method == "PET.rma", term == "b0") %>% .[["p.value"]]	
    rma.est <- res %>% filter(method == "PET.rma", term == "b0") %>% .[["estimate"]]
    usePEESE.rma <- ifelse(rma.p.value < p.crit & rma.est > 0, TRUE, FALSE)
    res <- rbind(res,
                 data.frame(method="PETPEESE.rma", if (usePEESE.rma == TRUE) {tidyRMA(PEESE.rma)} else {tidyRMA(PET.rma)})
    )
  }
  
  
  res <- plyr::rbind.fill(res, data.frame(method="PETPEESE.lm", term="estimator", type=ifelse(usePEESE.lm == TRUE, 4, 3)))
  
  if (runRMA == TRUE) {
    res <- plyr::rbind.fill(res, data.frame(method="PETPEESE.rma", term="estimator", type=ifelse(usePEESE.rma == TRUE, 4, 3)))	
  } else {
    res$method <- gsub(".lm", "", res$method, fixed=TRUE)
  }
  
  returnRes(res, long)
}
