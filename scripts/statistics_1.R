setwd("scripts/")

#ex 1.3
func1 <- function(lambda = 3.5, n = 30, m = 8) {
        epsilon = 1 - ppois(m - 1, lambda)
        1 - (1 - epsilon)^n
}

#ex 1.5  
func1(lambda = 0.5, m = 9, n = 100)

##Monte carlo
## granularity is inverse of probability -> 1/prob --> 1e6 simulations

maxes <- replicate(1/0.000001, {
        max(rpois(100, 0.5))
})
table(maxes)
mean(maxes >= 9)

#ex 1.7
mean(rpois(100, 3))
var(rpois(100, 3))

#ex 1.8

renv::install("bioc::BSgenome.Celegans.UCSC.ce2")
library(Biostrings)
library(BSgenome.Celegans.UCSC.ce2)

samplesize <- sum(alphabetFrequency(BSgenome.Celegans.UCSC.ce2$chrM))

#H0 -- EVery Letter has the same probability
prob0 <- rep(1/4, 4)
obs0 <- rmultinom(1000, samplesize, prob = prob0)

exp0 = prob0 * samplesize

stats <- function(exp, obs) {
        sum(((exp - obs)^2) / exp)
}

S0 <- apply(obs0, 2, stats, exp = exp0)

q95 <- quantile(S0, probs = 0.95)

obs1 <- alphabetFrequency(BSgenome.Celegans.UCSC.ce2$chrM)[1:4]

