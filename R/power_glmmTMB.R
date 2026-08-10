#' Approximate power for a mixed model using glmmTMB
#'
#' @param model A model formula.
#' @param terms Character vector of fixed effect names.
#' @param Alpha Significance level, default 0.05.
#' @param Data A data frame.
#'
#' @return A data frame of power approximations per fixed effect.
#' @examples
#' # see vignette for full example
#' @export


power_glmmTMB <- function(model, terms, alpha = 0.05, Data) {
  effect_formula <- as.formula(paste("~", terms))
  emm_obj <- emmeans(model, effect_formula)
  jt <- joint_tests(emm_obj)
  jt_df <- as.data.frame(jt)
  f<- as.numeric(as.character(jt_df[jt_df$`model term` == terms, "F.ratio"]))
  formula_short <- model$modelInfo$allForm$formula
  environment(formula_short) <- environment()


  b3 <- lmerTest::lmer(formula_short, Data)
  aov_b3 <- anova(b3, ddf = "Kenward-Roger")
  ndf <- aov_b3[terms, "NumDF"]
  noncent_param <- ndf * f
  ddf <- aov_b3[terms, "DenDF"]
  FCrit <- qf(1 - alpha, ndf, ddf, 0)
  power <- 1 - pf(FCrit, ndf, ddf, noncent_param)



  power_mat <- matrix(power)
  power_mat_t <- t(power_mat)

  colnames(power_mat_t) <-terms
  rownames(power_mat_t) <- "Power Approximation"

  print(power_mat_t)

}

