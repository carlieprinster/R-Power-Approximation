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






