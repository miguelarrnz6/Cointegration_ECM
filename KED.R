# Kremer, Ericsson, Dolado
# Critical Values ECM known Cointegrating vector
library(zoo)
library(dynlm)

T <- 200
B <- 50
N <- T+B
M <- 10000
aa <- c(0,0.5,1)

cvalues <- matrix(nrow=3, ncol=3)
ectest <- array(M)
q <- c(0.01, 0.05, 0.1)

for(j in 1:3){
  a <- aa[j]
  set.seed(123456)
  for(i in 1:M){
    u <- rnorm(N)
    dx <- u
    epsn <- rnorm(N)
    dy <- a*dx + epsn
    x <- cumsum(dx)
    y <- cumsum(dy)
    x <- x[-(1:B)]
    y <- y[-(1:B)]
    x <- as.zoo(x)
    y <- as.zoo(y)
    z <- y - x
    m <- dynlm(d(y) ~ L(z) + d(x))
    s <- summary(m)
    ectest[i] <- s$coeff[2,3]
  }
  cvalues[j,] <- quantile(ectest, q)
}

colnames(cvalues) <- c("1%", "5%", "10%")
rownames(cvalues) <- c("a=0", "a=0.5", "a=1")
print(cvalues)