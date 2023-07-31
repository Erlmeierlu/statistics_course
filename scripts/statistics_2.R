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

###2.5

library("Biostrings")
staph = readDNAStringSet("../data/staphsequence.ffn.txt", "fasta")
staph[1]
letterFrequency(staph[[1]], letters = "ACGT", OR = 0)
Biostrings::letterFrequency(staph[1], letters = "ACGT", OR = 0)

letterFrq = vapply(staph, letterFrequency, FUN.VALUE = numeric(4),
                   letters = "ACGT", OR = 0)
colnames(letterFrq) = paste0("gene", seq(along = staph))
tab10 = letterFrq[, 1:10]
computeProportions = function(x) { x/sum(x) }
prop10 = apply(tab10, 2, computeProportions)
round(prop10, digits = 2)


p0 = rowMeans(prop10)
p0

cs = colSums(tab10)
cs

expectedtab10 = outer(p0, cs, FUN = "*")

randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) } )

stat = function(obsvd, exptd) {
        sum((obsvd - exptd)^2 / exptd)
}

B = 1000

simulstat = replicate(B, {
        randomtab10 = sapply(cs, function(s) { rmultinom(1, s, p0) })
        stat(randomtab10, expectedtab10)
})
S1 = stat(tab10, expectedtab10)
sum(simulstat >= S1)



hist(simulstat, col = "lavender", breaks = seq(0, 75, length.out=50))
abline(v = S1, col = "red")
abline(v = quantile(simulstat, probs = c(0.95, 0.99)),
       col = c("darkgreen", "blue"), lty = 2)

ranchi2 <- rchisq(1000, 30)
hist(simulstat, col = "lavender", breaks = 50)
hist(ranchi2, col = "lavender", breaks = 50)
qs = ppoints(100)
quantile(simulstat, qs)
quantile(qchisq(qs, df = 30), qs)

hist(quantile(simulstat, qs), breaks = 50)
hist(quantile(ranchi2, qs), breaks = 50)

qqplot(ranchi2, simulstat, main = "",
       xlab = expression(chi[nu==30]^2), asp = 1, cex = 0.5, pch = 16)
abline(a = 0, b = 1, col = "red")
qqplot(qchisq(ppoints(B), df = 30), simulstat, main = "",
       xlab = expression(chi[nu==30]^2), asp = 1, cex = 0.5, pch = 16)
abline(a = 0, b = 1, col = "red")

1 - pchisq(S1, df = 30)

####2.15
chisq.test(Deuteranopia)

