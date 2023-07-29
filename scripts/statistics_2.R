setwd("scripts/")
load("../data/e100.RData")
e99 = e100[-which.max(e100)]

library("vcd")
gf1 = goodfit( e99, "poisson")
rootogram(gf1, xlab = "", rect_gp = gpar(fill = "chartreuse4"))

table(e100)
loglikelihood  =  function(lambda, data = e100) {
        sum(log(dpois(data, lambda)))
}

lambdas = seq(0.05, 0.95, length = 100)
loglik = vapply(lambdas, loglikelihood, numeric(1))
plot(lambdas, loglik, type = "l", col = "red", ylab = "", lwd = 2,
     xlab = expression(lambda))
m0 = mean(e100)
abline(v = m0, col = "blue", lwd = 2)
abline(h = loglikelihood(m0), col = "purple", lwd = 2)
m0

gf  =  goodfit(e100, "poisson")
names(gf)
