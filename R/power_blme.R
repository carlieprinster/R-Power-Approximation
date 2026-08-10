#' Approximate power for a mixed model using blme
#'
#' @param Formula A model formula.
#' @param Varcomp Numeric vector of random variance components.
#' @param Resid_var Residual variance.
#' @param Data A data frame.
#' @param Family Model family, e.g. "gaussian".
#' @param Random_Effects Character vector of random effect names.
#' @param Fixed_Effects Character vector of fixed effect names.
#' @param Alpha Significance level, default 0.05.
#'
#' @return A data frame of power approximations per fixed effect.
#' @examples
#' # see vignette for full example
#' @export


power_blme <- function(Formula, Varcomp, Data, Expanded_Data = NULL, Family, Fixed_Effects = NULL, Random_Effects = NULL, Alpha,
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

  fn.list <- lapply(1:length(Varcomp), function(y) {
    uv <- Varcomp[y]              # force/capture the value NOW
    function(x) fn(x, use.var = uv)
  })

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
      noncent_param[i] <- aov$npar[reff_index[i]] * f[i]
      reff_index_b3[i] <- which(rownames(aov_b3) == Fixed_Effects[i])
      ndf[i] <- aov$npar[reff_index_b3[i]]
      ddf[i] <- aov_b3$DenDF[reff_index_b3[i]]
      FCrit[i] <- qf(1 - Alpha, ndf[i], ddf[i], 0)
      power_mat[1, i] <- 1 - pf(FCrit[i], ndf[i], ddf[i], noncent_param[i])

    }
  }

  else if (Family == "binomial") {
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

    try(b3 <- lmerTest::lmer(Formula, Data))
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



