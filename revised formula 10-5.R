
#clean environment
rm(list=ls())

#Set working directory
setwd("~/Summer25")

#Include necessary libraries
library(blme)
library(lme4)
library(dplyr)
library(tidyr)
library(emmeans)

#function declaration
power_mm_8_30 <- function(Formula, Varcomp, Data, Expanded_Data = NULL, Family, Fixed_Effects = NULL, Random_Effects = NULL, Alpha,
                          Resid_var = NULL, tol = 1e-6, Optimizer = "nloptwrap", OptCtrl = list(maxeval = 2e5), NAGQ0initStep = 0) {
  
  set.seed(123)
 
  
  Resid_var_TEMP_notlikely <<- Resid_var
  
 
  new_random_effects <- Random_Effects
  for (k in 1:length(Random_Effects)) {
    if(is.element(":", unlist(strsplit(Random_Effects[k], "")))) {
      new_random_effects[k] <- paste("`", Random_Effects[k], "`", sep = "")
    }
  }
  
  
  form_char <- vector("list", length(Varcomp - 1))
  
  if(is.null(Resid_var_TEMP_notlikely) == TRUE) {
    fn <- function(x, use.var) {
      dnorm(x, sqrt(use.var), tol, log = TRUE)
    }
  } else {
    
    fn <- function(x, use.var) {
      dnorm(x, sqrt(use.var / Resid_var_TEMP_notlikely), tol, log = TRUE)
    }
  }
  
  fn.list <- lapply(1:length(Varcomp),
                    function(y) functional::Curry(fn, use.var = Varcomp[y]))
  
  for (i in 1:length(Varcomp)) {
    
    x <- 1:i
    
    custom_list <- rep("custom_prior", times = i)
    var_list <- paste(custom_list, x, sep = "")
    
    assign(var_list[i], fn.list[[i]])
    
    form_char[[i]] <- paste(new_random_effects[i], " ~ custom(", var_list[i],
                            ", chol = TRUE", ", scale = 'log')", sep = "")
    
  }
  
  form_list <- lapply(form_char, noquote)
  
  if (Family == "gaussian") {
    if (length(new_random_effects > 1)) {
      b2 <- blme::blmer(Formula, data = Data,
                        resid.prior = point(sqrt(Resid_var_TEMP_notlikely)),
                        cov.prior = form_list)
      
    } else {
      b2 <- blme::blmer(Formula, data = Data,
                        resid.prior = point(sqrt(Resid_var_TEMP_notlikely)),
                        cov.prior = new_random_effects ~ custom(custom_prior1,
                                                                chol = TRUE,
                                                                scale = "log"))
      
    }
  } 
  else {
    if (length(new_random_effects > 1)) {
      
      if (is.null(Resid_var_TEMP_notlikely)== TRUE){
        
        b2 <- blme::bglmer(Formula, data = Data, family = Family,
                           control = glmerControl(optimizer = Optimizer, optCtrl = OptCtrl, nAGQ0initStep = NAGQ0initStep),
                           cov.prior = form_list)
        
        
        
        
      }
      else{
        
        b2 <- blme::bglmer(Formula, data = Data, family = Family, resid.prior = point(sqrt(Resid_var_TEMP_notlikely)), 
                           control = glmerControl(optimizer = Optimizer, optCtrl = OptCtrl, nAGQ0initStep = NAGQ0initStep),
                           cov.prior = form_list)
        
      }
      
    }
    
    
    else {
      
      if (length(new_random_effects > 1)) {
        
        if (is.null(Resid_var_TEMP_notlikely) == TRUE){
          
          b2 <- blme::bglmer(Formula, data = Data, family = Family,
                             control = glmerControl(optimizer = Optimizer, optCtrl = OptCtrl, nAGQ0initStep = NAGQ0initStep),
                             cov.prior = new_random_effects ~ custom(custom_prior1,
                                                                     chol = TRUE,
                                                                     scale = "log"))
          
          
          
        }
        else{
          
          
          b2 <- blme::bglmer(Formula, data = Data, family = Family, 
                             control = glmerControl(optimizer = Optimizer, optCtrl = OptCtrl, nAGQ0initStep = NAGQ0initStep),
                             cov.prior = new_random_effects ~ custom(custom_prior1,
                                                                     chol = TRUE,
                                                                     scale = "log"))
          
          
        }
        
      }
      
      
    }
  }
  
  power_mat <- matrix(NA, 1, length(Fixed_Effects)) 
  
  if (Family == "gaussian"){
    
    
    
    b3 <- lmerTest::lmer(Formula, data = Data)
    aov_b3 <- anova(b3, ddf = "Kenward-Roger")
    
    Fixed_Effects <- stringr::str_replace(Fixed_Effects, "\\*", ":")
    aov <- anova(b2)
    reff_index <- vector("numeric", length(Fixed_Effects))
    f <- vector("numeric", length(Fixed_Effects))
    noncent_param <- vector("numeric", length(Fixed_Effects))
    reff_index_b3 <- vector("numeric", length(Fixed_Effects))
    ndf <- vector("numeric", length(Fixed_Effects))
    ddf <- vector("numeric", length(Fixed_Effects))
    FCrit <- vector("numeric", length(Fixed_Effects))

    
    for (i in 1:length(Fixed_Effects)) {
      
      reff_index[i] <- which(rownames(aov) == Fixed_Effects[i])
      f[i] <- aov$`F value`[reff_index[i]]
      print(f[i])
      noncent_param[i] <- aov$npar[reff_index[i]] * f[i]
      reff_index_b3[i] <- which(rownames(aov_b3) == Fixed_Effects[i])
      ndf[i] <- aov$npar[reff_index_b3[i]]
      ddf[i] <- aov_b3$DenDF[reff_index_b3[i]]
      FCrit[i] <- qf(1 - Alpha, ndf[i], ddf[i], 0)
      power_mat[1, i] <- 1 - pf(FCrit[i], ndf[i], ddf[i], noncent_param[i])
      
    }
  }

  else if (Family == "binomial") {
      print(Expanded_Data)
        Data1 <- Expanded_Data
        formula_short <- update(Formula, rnorm(length(Expanded_Data[,1])) ~ .)
        environment(formula_short) <- environment()

    if (is.null(Data1)){
      Data1 = Data
      formula_short <- update(Formula, rnorm(length(Data[,1])) ~ .)
      environment(formula_short) <- environment()
      
    }
    
    
    try(b3 <- lmerTest::lmer(formula_short, Data1))
    aov_b3 <- anova(b3, ddf = "Kenward-Roger")
    
    reff_index <- vector("numeric", length(Fixed_Effects))
    f <- vector("numeric", length(Fixed_Effects))
    noncent_param <- vector("numeric", length(Fixed_Effects))
    reff_index_b3 <- vector("numeric", length(Fixed_Effects))
    ndf <- vector("numeric", length(Fixed_Effects))
    ddf <- vector("numeric", length(Fixed_Effects))
    FCrit <- vector("numeric", length(Fixed_Effects))
    
    for (i in 1:length(Fixed_Effects)) {
      effect_formula <- as.formula(paste("~", Fixed_Effects[i]))
      
      emmeans_result <- emmeans(b2, effect_formula)
      joint_test <- as.data.frame(joint_tests(emmeans_result))
      effects <- joint_test$`model term`
      f[i] <- joint_test[joint_test$`model term` == Fixed_Effects[i], "F.ratio"]
      ndf[i] <- aov_b3[Fixed_Effects[i],"NumDF"]
      
      noncent_param[i] <- ndf[i] * f[i]
      
      ddf[i] <- aov_b3[Fixed_Effects[i],"DenDF"]
      FCrit[i] <- qf(1 - Alpha, ndf[i], ddf[i], 0)
      power_mat[1, i] <- 1 - pf(FCrit[i], ndf[i], ddf[i], noncent_param[i])
    }
  }
  else{
    
    try(b3 <- lmerTest::lmer(formula_short, Data1))
    aov_b3 <- anova(b3, ddf = "Kenward-Roger")
    
    reff_index <- vector("numeric", length(Fixed_Effects))
    f <- vector("numeric", length(Fixed_Effects))
    noncent_param <- vector("numeric", length(Fixed_Effects))
    reff_index_b3 <- vector("numeric", length(Fixed_Effects))
    ndf <- vector("numeric", length(Fixed_Effects))
    ddf <- vector("numeric", length(Fixed_Effects))
    FCrit <- vector("numeric", length(Fixed_Effects))
    
    
    for (i in 1:length(Fixed_Effects)) {
      effect_formula <- as.formula(paste("~", Fixed_Effects[i]))
      emmeans_result <- emmeans(b2, effect_formula)
      joint_test <- as.data.frame(joint_tests(emmeans_result))
      effects <- joint_test$`model term`
      f[i] <- joint_test[joint_test$`model term` == Fixed_Effects[i], "F.ratio"]
      ndf[i] <- aov_b3[Fixed_Effects[i],"NumDF"]
      noncent_param[i] <- ndf[i] * f[i]
      ddf[i] <- aov_b3[Fixed_Effects[i],"DenDF"]
      FCrit[i] <- qf(1 - Alpha, ndf[i], ddf[i], 0)
      power_mat[1, i] <- 1 - pf(FCrit[i], ndf[i], ddf[i], noncent_param[i])
    }
  }

  
power_mat <- as.matrix(power_mat)
  
  colnames(power_mat) <- Fixed_Effects
  rownames(power_mat) <- "Power Approximation"
  print(power_mat)
  
}

poisson_example <- power_mm_8_30(Formula = exp_count ~ rate + trt*rate + (1 | block) + (1 | block:trt), Varcomp = c(0.25,0.15), 
              Resid_var = NULL, Data = ex_poisson, Family = "poisson", 
              Random_Effects = c("block", "block:trt"), Fixed_Effects = c("trt", "rate", "rate:trt"), Alpha = 0.05,
              tol = 1e-6, Optimizer = "nloptwrap", OptCtrl = list(maxeval = 2e5), NAGQ0initStep = 0)


gaussian_example <-  power_mm_8_30(Formula = estY ~ NSource*Thatch + (1 | Field) + (1|Field:NSource), 
                                  Varcomp = c(0.008, 0.07), Resid_var = 0.2, Data = ex_gaussian,
                                  Family = "gaussian", Fixed_Effects =c("NSource", "NSource:Thatch", "Thatch"),
                                  Alpha = 0.05)


binomial_example <- power_mm_8_30(Formula = cbind((expected_y), n - (expected_y)) ~ 1 + trt + (1 | location) + 
                                    (1 | location:trt), Varcomp = c(0.02,0.05), Resid_var = NULL, Data = ex_binomial, Expanded_Data = expanded, 
                                  Family = "binomial", Random_Effects = c("location", "location:trt"), Fixed_Effects = c("trt")  ,Alpha = 0.05, tol = 1e-6, 
                                  Optimizer = "nloptwrap", OptCtrl = list(maxeval = 2e5), NAGQ0initStep = 0)


gaussian_example2 <- power_mm_8_30(Formula = Reaction ~  Days + (1 | Subject), 
                                   Varcomp = 1378.2, Resid_var = 960.5, Data = sleepstudy,
                                   Family = "gaussian", Random_Effects = c("Subject"), Fixed_Effects =c("Days"),
                                   Alpha = 0.05)

