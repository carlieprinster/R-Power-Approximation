#Introduce Datasets

#poisson dataset
trt <- as.factor(c(rep(1, 12), rep(2, 12)))
rate <- as.factor(rep(c(rep(1, 4), rep(2, 4), rep(3, 4)), 2))
exp_count <- c(rep(10, 4), rep(9, 4), rep(8, 4), rep(9, 4), rep(6, 4), rep(3, 4))
block <- as.factor(rep(1:4, 6))
ex_poisson <- cbind.data.frame(trt, rate, exp_count, block)

#binomial dataset
trt <- as.factor(c(rep(0, 4), rep(1, 4)))
n <- rep(65, 8)
pi <- c(rep(0.15, 4), rep(0.25, 4))
location <- as.factor(rep(1:4, 2))
expected_y <- n * pi
ex_binomial <- cbind.data.frame(trt, n, pi, location, expected_y)

#gaussian dataset
Thatch <- as.factor(c(rep(2, 16), rep(5, 16), rep(8, 16)))
NSource <- rep(c(rep("AmmSulph", 4), rep("IBDU", 4), rep("SCUrea", 4), rep("Urea", 4)), 3)
estY <- c(rep(6, 12), rep(4, 4), rep(c(rep(6, 4), rep(7, 4), rep(8, 4), rep(5, 4)), 2))
Field <- as.factor(rep(1:4, 12))
ex_gaussian <- cbind.data.frame(Thatch, NSource, estY, Field)



###############glmmTMB function#####################################

#Import necessary libraries
library(car)
library("glmmTMB")
library(emmeans)

power_F_glmmTMB <- function(model, fixed_terms, alpha = 0.05, Data) {
  print(class(fixed_terms))
  effect_formula <- as.formula(paste("~as.factor(",fixed_terms,")"))
  emm_obj <- emmeans(model, effect_formula)
  jt <- joint_tests(emm_obj)
  jt_df <- as.data.frame(jt)
  f<- as.numeric(as.character(jt_df[jt_df$`model term` == (paste("as.factor(",fixed_terms,")")), "F.ratio"]))
  print(f)
  formula_short <- model$modelInfo$allForm$formula
  environment(formula_short) <- environment()
  
  
  print(formula_short)
  b3 <- lmerTest::lmer(formula_short, Data)
  aov_b3 <- anova(b3, ddf = "Kenward-Roger")
  print(aov_b3)
  ndf <- aov_b3["NSource:Thatch", "NumDF"]
  print(ndf)
  noncent_param <- ndf * f
  print(noncent_param)
  ddf <- aov_b3["NSource:Thatch", "DenDF"]
  print(ddf)
  FCrit <- qf(1 - alpha, ndf, ddf, 0)
  print(FCrit)
  power <- 1 - pf(FCrit, ndf, ddf, noncent_param)
  print(power)
  
  
  
  data_frame = data.frame(power, FCrit, ndf, ddf , noncent_param)
  
  
  
  return(data_frame)
  
}




