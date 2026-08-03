
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glmpa

<!-- badges: start -->

<!-- badges: end -->

The goal of glmpa is to approximate power of generalized linear mixed
models using the probability distribution method (PDM).

This goal can be achieved using this package, as it includes two
functions that do this.

These two functions both extract all the components needed to implement
the PDM for a power approximation of a given fixed effect.  
The needed components include a numerator and denominator degrees of
freedom, a non-centrality parameter, and an F-statistic for the given
fixed effect.

One function uses the blme (Citation) package to set (steep????) priors
to fix the random variance components in order to implement the PDM.

The other function uses the glmmTMB (Citation) package through its start
and map arguments to fix random variance components in order to
implement the PDM.

Additional considerations for the decision between these two functions…

Implementation issues to consider…

The functions using the blme package is a development on an existing
function created to approximate power of linear mixed models using the
probability distribution method.

## Installation

You can install the development version of glmpa like so:

\##Example

Using these functions first involves creating a data frame for the
dataset that “informs” the model.

The example is to do with an experiment that seeks to assess the effects
of two treatment factors, nitrogen source and years of thatch
accumulation, on chlorophyll content in grass clippings of golf greens
(Kuehl 1994, example 14.1). This split plot designed experiment has an
estimated response distribution from the Gaussian family.

``` r
#Build Dataset
Thatch <- as.factor(c(rep(2, 16), rep(5, 16), rep(8, 16)))
NSource <- rep(c(rep("AmmSulph", 4), rep("IBDU", 4), rep("SCUrea", 4),
                 rep("Urea", 4)), 3)
estY <- c(rep(6, 12), rep(4, 4), rep(c(rep(6, 4), rep(7, 4), rep(8, 4),
                                       rep(5, 4)), 2))
Field <- as.factor(rep(1:4, 12))
ex_gaussian <- cbind.data.frame(Thatch, NSource, estY, Field)
```

# Example for blme function

One can implement the function by defining the following as arguments:
model formula, the values of the random variance components, the
previously created data frame, the model distribution family, the names
of the fixed effects of power approximation interest, and the names of
the random variance components. All other arguments have default values
that can be customized as users need.

``` r

gaussian_example <-  power_mm(Formula = estY ~ NSource*Thatch + (1 | Field) + (1|Field:NSource), 
                                  Varcomp = c(0.008, 0.07), Resid_var = 0.2, Data = ex_gaussian,
                                  Family = "gaussian", Random_Effects = c("Field", "Field:NSource"), Fixed_Effects =c("NSource", "NSource:Thatch", "Thatch"), Alpha = 0.05)
#>                     NSource NSource:Thatch    Thatch
#> Power Approximation       1       0.947531 0.9999965
```

# Example for glmmTMB function

One can implement the other function by first defining a model within
the glmmTMB package framework and defining the start and map arguments
with the fixed random variance components. Then, the user must define
the following: a glmmTMB model, the name of the fixed effect of power
approximation interest, and the model dataset. The last argument, the
alpha value, has a default value of 0.05 that can be customized as users
need.

``` r

#build model
start = list(theta = log(c(sqrt(0.008), sqrt(0.07))), betadisp = log(sqrt(0.2)))
map = list(theta = factor(c(NA, NA)), betadisp = factor(NA))

ex_gaussianmodel <- glmmTMB(estY ~ NSource * Thatch + (1 | Field) + (1 | Field:NSource),
                            data = ex_gaussian,
                            family = gaussian(),
                             control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")),
                            start = start,
                            map = map)

#call function
gaussian_glmmTMB_example <- power_F_glmmTMB(ex_gaussianmodel, "NSource:Thatch", 0.05, ex_gaussian)
#>                     NSource:Thatch
#> Power Approximation      0.9475087
```
