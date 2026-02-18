###############glmmTMB function#####################################

#Import necessary libraries
library(car)
library("glmmTMB")
library(emmeans)
power_F_glmmTMB <- function(model, terms, alpha = 0.05, Data) {
  effect_formula <- as.formula(paste("~", terms))
  emm_obj <- emmeans(model, effect_formula)
  jt <- joint_tests(emm_obj)
  jt_df <- as.data.frame(jt)
  f<- as.numeric(as.character(jt_df[jt_df$`model term` == terms, "F.ratio"]))
  print(f)
  formula_short <- model$modelInfo$allForm$formula
  environment(formula_short) <- environment()
  
  
  print(formula_short)
  b3 <- lmerTest::lmer(formula_short, Data)
  print(b3)
  aov_b3 <- anova(b3, ddf = "Kenward-Roger")
  print(aov_b3)
  ndf <- aov_b3[terms, "NumDF"]
  print(ndf)
  noncent_param <- ndf * f
  print(noncent_param)
  ddf <- aov_b3[terms, "DenDF"]
  print(ddf)
  FCrit <- qf(1 - alpha, ndf, ddf, 0)
  print(FCrit)
  power <- 1 - pf(FCrit, ndf, ddf, noncent_param)
  print(power)
  
  
  
  data_frame = data.frame(power, FCrit, ndf, ddf , noncent_param)
  
  
  
  return(data_frame)
  
}

###gaussian power###

#build model
start = list(theta = log(c(sqrt(0.008), sqrt(0.07))), betadisp = log(sqrt(0.2)))
map = list(theta = factor(c(NA, NA)), betadisp = factor(NA))

ex_gaussianmodel <- glmmTMB(estY ~ NSource * Thatch + (1 | Field) + (1 | Field:NSource),
                            data = ex_gaussian,
                            family = gaussian(),
                            start = start,
                            map = map)

#call function
power_F_glmmTMB(ex_gaussianmodel, "NSource:Thatch", 0.05, ex_gaussian)

##gaussian power - sleepstudy##
#build model
start = list(theta = c(log(sqrt(1378.2))), betadisp = log(sqrt(960.5)))
map = list(theta = factor(c(NA)), betadisp = factor(NA))

sleep_model <- glmmTMB(Reaction ~  Days + (1 | Subject), 
                       data = sleepstudy, 
                       family = gaussian(),
                       start = start, 
                       map = map)

power_F_glmmTMB(sleep_model, "Days", 0.05, sleepstudy)

###binomial power###

#build model
start = list(theta = log(c(sqrt(0.02), sqrt(0.05))))
map = list(theta = factor(c(NA, NA)))

ex_binomialmodel <- glmmTMB(y ~ 1 + trt + (1 | location) + (1 | location:trt), 
                            data = expanded_data, 
                            family = binomial(), 
                            REML = FALSE,
                            start = start, 
                            map = map)
#call function
power_F_glmmTMB(ex_binomialmodel, "trt", 0.05, expanded_data)


###poisson power###

#build model
start <- list(theta = log(c(sqrt(0.25), sqrt(0.15))))
map <- list(theta = factor(c(NA, NA)))

ex_poissonmodel <- glmmTMB(exp_count ~ rate + trt*rate  + (1|block) + (1|block:trt), 
                           data = ex_poisson, 
                           family = poisson(), 
                           start = start, 
                           map = map)

#call function
power_F_glmmTMB(ex_poissonmodel, "rate", 0.05, ex_poisson)


