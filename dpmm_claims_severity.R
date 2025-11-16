## ----------------------------------------------------------------------------
## This is a R code for the DPMM model applied to log(claims severity) 
## data consists of its MCMC algorithm, its posterior predictive
## distribution and evaluating its clustering performance. In this
## code, we assume that we have two predictor variables and thus
## three BNP regression parameters. The definitions of variables
## used in this code are listed here:
##
## n: size of the training data set
## y: vector of response variables which is log of amount of claims
## X: design matrix of predictor variables
## vn: number of columns in design matrix which is 3 here
## g0mu: mean of the multivariate normal base measure
## shape: shape parameter of the inverse gamma base measure
## rate: rate parameter of the inverse gamma base measure
## n0: coefficient factor in the covariance matrix of the
##     multivariate normal base measure
## ----------------------------------------------------------------------------

# Install and load required packages

library(zoo)
library(xts)

# Install the 'CASdatasets' package, which contains the French motor insurance claims dataset
install.packages("CASdatasets", 
                 repos = "https://dutangc.perso.math.cnrs.fr/RRepository/pub/", 
                 type="source")

install.packages("factoextra", repos = "https://cloud.r-project.org")

library(CASdatasets) # French motor third-party liability insurance datasets
library(mvtnorm) # multivariate normal
library(coda) # diagnostic tests of convergence of the Markov chain
library(tidyverse) # data manipulation and creating graphics
library(cluster) # clustering algorithms
library(factoextra) # clustering visualization
library(parallel) # parallel computation
library(foreach) # parallel computation
library(doMC) # provides a parallel backend

# register the multicore parallel backend with the foreach package
registerDoMC(2)

# check the number of execution workers after registered
getDoParWorkers()
getDoParRegistered()
getDoParName()

# Load data
data("freMTPLfreq")
data("freMTPLsev")

freq <- freMTPLfreq
sev <- freMTPLsev

data <- freq[sev$PolicyID, ]
data <- cbind(data, sev$ClaimAmount)
colnames(data)[11] <- "ClaimAmount"

nrow(data)
head(data)

# Train Test split
ntrain <- round((0.1) * nrow(data))
set.seed(1099)
train <- sample(nrow(data), ntrain)
dataTrain <- data[train, ]
dataTest <- data[-train, ]

# y: vector of response variables which is log of amount of claims
y <- log(dataTrain$ClaimAmount)

# n: size of the training data set
n <- length(y)

# Data preprocessing
X0 <- rep(1, n)

X1 <- dataTrain$CarAge
X1st <- (X1 - mean(X1)) / sd(X1)

X2 <- dataTrain$DriverAge
X2st <- (X2 - mean(X2)) / sd(X2)

xVar <- cbind(X0, X1, X2)

# X: design matrix of predictor variables
X <- cbind(X0, X1st, X2st)

# vn: number of columns in design matrix which is 3 here
vn <- ncol(X)

# Function 'configure' will be used in Neal's Algorithm 8
# to create configuration vector
configure <- function(m) {
  n <- unique(m)
  conf <- rep(0, nrow(m))
  for (j in 1:nrow(n)) {
    for (i in 1:nrow(m)) {
      if (identical(m[i, ], n[j, ])) {
        conf[i] <- j
      }
    }
  }
  return(conf)
}

# Function 'reorder' will be used in Neal's Algorithm 8
reorder <- function(v, k) {
  for (i in 1:length(v)) {
    if (v[i] > v[k]) {
      v[i] <- v[i] - 1
    }
  }
  return(v)
}

# Function 'n_i' will be used in Neal's Algorithm 8
n_i <- function(v, k1, k2) {
  ifelse(k1 == v[k2], length(v[v == k1]) - 1, length(v[v == k1]))
}

# shape: shape parameter of the inverse gamma base measure
shape <- 3

# rate: rate parameter of the inverse gamma base measure
rate <- 5

# g0mu: mean of the multivariate normal base measure
g0mu <- rep(0, vn)

# n0: coefficient factor in the covariance matrix of the
#     multivariate normal base measure
n0 <- 0.5

# Base probability measure, G_0, in DPMM model
g0 <- function(w) {
  s2 <- 1 / rgamma(w, shape = shape, rate = rate)
  bm <- NULL
  for (i in 1:w) {
    bm <- rbind(bm, rmvnorm(1, mean = g0mu, sigma = n0 * s2[i] * diag(vn)))
  }
  return(cbind(bm, s2))
}

#### -------------------------------------------------------------------
#### Neal's Algorithm8 to sample from the posterior distribution of
#### the DPMM
#### -------------------------------------------------------------------

neal_algorithm8 <- function(betasigma0, m, alpha0, n_iter, a1, a2) {

  # The function parameters are:
  # betasigma0: initial value of model parameters matrix
  # m: auxiliary parameter
  # alpha0: initial value of alpha in DP prior
  # n_iter: number of MCMC algorithm iterations
  # a1: shape parameter of gamma prior for updating alpha
  # a2: rate parameter of gamma prior for updating alpha

  betasigma <- betasigma0 # posterior samples of model parameters

  # BNP regression coefficients, 'beta's
  for (i in 1:vn) {
    assign(paste("beta", i - 1, sep = "_"), betasigma0[, i])
  }

  sigma <- betasigma0[, vn + 1] # vector of variances in normal likelihood

  phi <- unique(betasigma0) # distinct components

  ndistinct <- nrow(phi) # number of distinct components

  c <- configure(betasigma0) # configuration vector

  temp <- betasigma0 # temporary matrix

  a <- shape # shape parameter
  b <- rate # rate parameter

  alpha <- alpha0 # precision parameter
  precision <- alpha # posterior samples of alpha

  for (k in 2:n_iter) {
    for (i in 1:n) {
      k_minus <- length(unique(c[-i]))
      h <- k_minus + m

      if (n_i(c, c[i], i) != 0) {
        phi <- rbind(phi, g0(m))
      } else {
        phi <- rbind(phi[-c[i], ], phi[c[i], ])
        c <- reorder(c, i)
        c[i] <- k_minus + 1
        phi <- rbind(phi, g0(m - 1))
      }

      prob <- rep(0, h)

      for (l in 1:k_minus) {
        prob[l] <- (log(n_i(c, l, i)) - (1 / 2) * log(phi[l, vn + 1]) - (1 / 2) *
          (1 / (phi[l, vn + 1])) * (y[i] - X[i, ] %*% phi[l, 1:vn])^2)
      }

      for (l in (k_minus + 1):h) {
        prob[l] <- (log((alpha) / m) - (1 / 2) * log(phi[l, vn + 1]) - (1 / 2) *
          (1 / (phi[l, vn + 1])) * (y[i] - X[i, ] %*% phi[l, 1:vn])^2)
      }

      prob <- prob - max(prob)
      prob <- exp(prob)

      c[i] <- sample(1:h, size = 1, replace = T, prob = prob)

      if (c[i] > k_minus) {
        phi <- rbind(phi[1:k_minus, ], phi[c[i], ])
        c[i] <- k_minus + 1
      } else {
        phi <- phi[1:k_minus, ]
      }
    }

    phi <- as.matrix(phi)
    if (ncol(phi) == 1) {
      phi <- t(phi)
    }

    # This step of the Neal’s algorithm is using the
    # conjugacy of the base probability measure G0 to
    # the likelihood in order to implement the Gibbs
    # sampling update of distinct components
    for (j in 1:nrow(phi)) {
      abc <- (c == j)
      index <- which(abc == TRUE)

      dt <- t(rep(0, vn))
      M <- (1 / n0) * diag(vn)

      for (w in index) {
        dt <- dt + y[w] * t(X[w, ])
        M <- M + X[w, ] %*% t(X[w, ])
      }

      shapePost <- length(index) / 2 + a + vn / 2
      scalePost <- (b + (1 / 2) * sum(y[index]^2) - (1 / 2) * dt %*% solve(M) %*%
        t(dt))

      phi[j, vn + 1] <- 1 / rgamma(1, shape = shapePost, rate = scalePost)

      meanPost <- solve(M) %*% t(dt)
      varPost <- phi[j, vn + 1] * solve(M)

      phi[j, 1:vn] <- rmvnorm(1, mean = meanPost, sigma = varPost)
    }

    temp <- phi[c, ]
    betasigma <- cbind(betasigma, temp)

    for (i in 1:vn) {
      assign(paste("beta", i - 1, sep = "_"), rbind(get(paste("beta", i - 1,
        sep = "_"
      )), temp[, i]))
    }

    sigma <- rbind(sigma, temp[, vn + 1])

    # Number of distinct components
    nd <- nrow(phi)
    ndistinct <- c(ndistinct, nd)

    # Update alpha by Gibbs sampler in Escobar and West (1995)
    eta <- rbeta(1, alpha + 1, n)
    t <- sample(1:2, size = 1, replace = T, prob = c(a1 + nd - 1, n * (a2 - log(eta))))
    alpha <- ifelse(t == 1, rgamma(1, shape = a1 + nd, rate = a2 - log(eta)),
      rgamma(1, shape = a1 + nd - 1, rate = a2 - log(eta))
    )

    precision <- c(precision, alpha)
  }

  b <- list(
    "betasigma" = betasigma,
    "beta_0" = beta_0, "beta_1" = beta_1, "beta_2" = beta_2,
    "sigma" = sigma,
    "ndistinct" = ndistinct,
    "precision" = precision
  )
  return(b)
}

n_iter <- 50000 # number of MCMC algorithm iterations
burnin <- (0.5) * n_iter # burn-in
alpha0 <- 1
m <- 2
a1 <- 1
a2 <- 1

betasigma01 <- g0(n) # initial values

stime <- system.time(b1 <- neal_algorithm8(
  betasigma0 = betasigma01, m = m, alpha0 = alpha0,
  n_iter = n_iter, a1 = a1, a2 = a2
))

betasigma1 <- b1$betasigma

for (i in 1:vn) {
  assign(paste("beta1", i - 1, sep = "_"), b1[[paste("beta", i - 1, sep = "_")]])
}

sigma1 <- b1$sigma

ndistinct1 <- b1$ndistinct
hist(ndistinct1[-(1:burnin)], xlab = "Clusters", main = "Histogram", breaks = 100)
k1 <- round(mean(ndistinct1[-(1:burnin)]))

precision1 <- b1$precision
plot(precision1[-(1:burnin)], type = "l", ylab = "Precision parameter")
boxplot(precision1[-(1:burnin)])

# Diagnostic tests of convergence
mc1 <- as.mcmc(cbind(beta1_0[-(1:burnin), ], beta1_1[-(1:burnin), ],
                     beta1_2[-(1:burnin), ], sigma1[-(1:burnin), ]))

par(mfrow = c(1, 1))

# Trace plots and density function estimate plots from MCMC output
for (i in 1:vn) {
  plot(mc1[, (i - 1) * n + 1], main = substitute(beta[s[t]], list(s = i - 1, t = 1)))
  cat("\\")
  cat("\\linebreak")
  plot(mc1[, (i - 1) * n + n / 2], main = substitute(beta[s[t]], list(
    s = i - 1,
    t = n / 2
  )))
  cat("\\")
  cat("\\linebreak")
  plot(mc1[, (i) * n], main = substitute(beta[s[t]], list(s = i - 1, t = n)))
  cat("\\")
  cat("\\linebreak")
}

plot(mc1[, (vn) * n + 1], main = substitute(sigma[s]^t, list(s = 1, t = 2)))
plot(mc1[, (vn) * n + n / 2], main = substitute(sigma[s]^t, list(s = n / 2, t = 2)))
plot(mc1[, (vn + 1) * n], main = substitute(sigma[s]^t, list(s = n, t = 2)))

par(mfrow = c(1, 1))

# Trace plots of MCMC output
for (i in 1:vn) {
  traceplot(mc1[, (i - 1) * n + 1], main = substitute(beta[s[t]], list(s = i -
    1, t = 1)))
  cat("\\")
  cat("\\linebreak")
  traceplot(mc1[, (i - 1) * n + n / 2], main = substitute(beta[s[t]], list(s = i -
    1, t = n / 2)))
  cat("\\")
  cat("\\linebreak")
  traceplot(mc1[, (i) * n], main = substitute(beta[s[t]], list(s = i - 1, t = n)))
  cat("\\")
  cat("\\linebreak")
}

traceplot(mc1[, (vn) * n + 1], main = substitute(sigma[s]^t, list(s = 1, t = 2)))
traceplot(mc1[, (vn) * n + n / 2], main = substitute(sigma[s]^t, list(s = n / 2, t = 2)))
traceplot(mc1[, (vn + 1) * n], main = substitute(sigma[s]^t, list(s = n, t = 2)))

par(mfrow = c(1, 1))

# Autocorrelation function plots for Markov Chains
for (i in 1:vn) {
  autocorr.plot(mc1[, (i - 1) * n + 1], main = substitute(beta[s[t]], list(s = i -
    1, t = 1)), lag.max = 100)
  cat("\\")
  cat("\\linebreak")
  autocorr.plot(mc1[, (i - 1) * n + n / 2], main = substitute(beta[s[t]], list(s = i -
    1, t = n / 2)), lag.max = 100)
  cat("\\")
  cat("\\linebreak")
  autocorr.plot(mc1[, (i) * n],
    main = substitute(beta[s[t]], list(s = i - 1, t = n)),
    lag.max = 100
  )
  cat("\\")
  cat("\\linebreak")
}

autocorr.plot(mc1[, (vn) * n + 1],
  main = substitute(sigma[s]^t, list(s = 1, t = 2)),
  lag.max = 100
)
autocorr.plot(mc1[, (vn) * n + n / 2], main = substitute(sigma[s]^t, list(
  s = n / 2,
  t = 2
)), lag.max = 100)
autocorr.plot(mc1[, (vn + 1) * n],
  main = substitute(sigma[s]^t, list(s = n, t = 2)),
  lag.max = 100
)

par(mfrow = c(1, 1))

# Boxplots
for (i in 1:vn) {
  boxplot(get(paste("beta1", i - 1, sep = "_"))[-(1:burnin), 1], main = paste("beta",
    i - 1, 1,
    sep = "_"
  ))
  cat("\\")
  cat("\\linebreak")
  boxplot(get(paste("beta1", i - 1, sep = "_"))[-(1:burnin), n / 2], main = paste("beta",
    i - 1, n / 2,
    sep = "_"
  ))
  cat("\\")
  cat("\\linebreak")
  boxplot(get(paste("beta1", i - 1, sep = "_"))[-(1:burnin), n], main = paste("beta",
    i - 1, n,
    sep = "_"
  ))
  cat("\\")
  cat("\\linebreak")
}

boxplot(sigma1[-(1:burnin), 1], main = paste("sigma2", 1, sep = "_"))
boxplot(sigma1[-(1:burnin), n / 2], main = paste("sigma2", n / 2, sep = "_"))
boxplot(sigma1[-(1:burnin), n], main = paste("sigma2", n, sep = "_"))


#### -------------------------------------------------------------------
#### Posterior predictive distribution
#### -------------------------------------------------------------------

pred <- function(ynew, xnew) {

  # ynew: new data point
  # xnew: vector of risk factors

  dtn1 <- ynew * t(xnew)
  Mn1 <- (1 / n0) * diag(vn) + xnew %*% t(xnew)

  a <- shape
  b <- rate

  coef1 <- mean((precision1[-(1:burnin)]) / (precision1[-(1:burnin)] + n))

  part1 <- ((coef1) * (1 / sqrt(2 * pi)) * ((n0)^(-vn / 2)) * (b^a) * (1 / gamma(a)) *
    (1 / sqrt(det(Mn1))) * gamma(1 / 2 + a) * (b + (1 / 2) * ynew^2 - (1 / 2) * dtn1 %*%
      solve(Mn1) %*% t(dtn1))^(-(1 / 2 + a)))

  f <- function(p) {
    ((1 / sqrt(2 * pi)) * (1 / sqrt(p[vn + 1])) * exp(-(1 / (2 * p[vn + 1])) * (ynew -
      t(xnew) %*% p[1:vn])^2))
  }

  part2 <- 0

  for (t in (burnin + 1):n_iter) {
    temp_theta <- betasigma1[, ((vn + 1) * t - (vn + 1) + 1):((vn + 1) * t)]
    for (v1 in 1:n) {
      part2 <- (part2 + (f(temp_theta[v1, ]) / (precision1[t] + n)))
    }
  }

  part2 <- (part2 / (n_iter - burnin))

  result <- (part1 + part2)

  return(result = result)
}

xtest <- c(1, 9, 38)

sdx <- apply(xVar, 2, sd)
meanx <- apply(xVar, 2, mean)
xnew <- (xtest - meanx) / sdx
xnew[1] <- 1

lTest <- nrow(dataTest[dataTest$CarAge == xtest[2] & dataTest$DriverAge == xtest[3], ])
yTest <- log(dataTest[dataTest$CarAge == xtest[2] & dataTest$DriverAge == xtest[3], ]$ClaimAmount)
my <- round(max(yTest))

pr <- foreach(ynew = 0:my, .combine = "c") %do% {pred(ynew = ynew, xnew = xnew)}

yHist = hist(yTest, breaks=((0:(my+1))-0.5),plot = F)

# MSE and SE
MSE <- mean((yHist$counts - pr*lTest)^2)
SE_MSE <- sd((yHist$counts - pr * lTest)^2)/sqrt(length((yHist$counts - pr * lTest)))


# Posterior predictive distribution plot

yHist <- hist(
  yTest,
  breaks = ((0:(my+1)) - 0.5),
  plot = FALSE
)

col_obs_fill   <- rgb(33, 158, 188, 90,  maxColorValue = 255)  
col_obs_fill2 <- rgb(33, 158, 188, 200, maxColorValue = 255)
col_obs_border <- rgb(70, 70, 70, 220,    maxColorValue = 255)   
col_pred       <- "#D00000"  

ylim_max <- 0.7

hist(
  yTest,
  breaks = ((0:(my+1)) - 0.5),
  freq = FALSE,
  xlim = c(0, my + 0.5),
  ylim = c(0, ylim_max),
  main = "DPMM",
  xlab = "log(Claim Amount)",
  ylab   = "Density",
  col    = col_obs_fill,
  border = col_obs_border
)

points(0:my, pr, pch = 16, col = col_pred)
lines(0:my, pr, lwd = 2, col = col_pred)

legend(
  x = my * 0.1, y = ylim_max * 0.9,
  legend = c("Observed", "Predicted"),
  pch = c(15, 16),
  pt.cex = c(1.2, 1),
  lwd = c(NA, 2),
  col = c(col_obs_fill2, col_pred),
  bty = "o",
  bg = "white",
  box.lwd= 1
)


#### -------------------------------------------------------------------
#### Hierarchical Clustering
#### -------------------------------------------------------------------

burnin_cluster <- burnin # burn-in

mat1 <- function(v) {
  m <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in (i):n) {
      if (identical(v[i, ], v[j, ])) {
        m[i, j] <- 1
      }
    }
  }
  return(m)
}

# Function to create the Dissimilarity Matrix from the posterior samples of parameters
dissimilarity_mat <- function(b) {
  s <- matrix(0, nrow = n, ncol = n)
  for (i in (burnin_cluster + 1):n_iter) {
    s <- s + mat1(b[, ((vn + 1) * i - (vn + 1) + 1):((vn + 1) * i)])
  }
  s <- s + t(s) - (n_iter - burnin_cluster) * diag(n)
  d <- 1 - (1 / (n_iter - burnin_cluster)) * s
  return(d)
}

stime <- system.time(d1 <- dissimilarity_mat(betasigma1))

dist1 <- as.dist(d1)

# Hierarchical cluster analysis on the set of dissimilarities
hc1 <- hclust(dist1, method = "ward.D2")
k <- k1
sub_grp <- cutree(hc1, k = k)
table(sub_grp)

# Cluster Dendrogram
plot(hc1, cex = 0.6, hang = -1)
if (k > 1) {
  rect.hclust(hc1, k = k, border = 2:(2 + k))
}

# Visualize the heatmap of the distance matrix
fviz_dist(dist1, show_labels = FALSE, gradient = list(
  low = "#00AFBB", mid = "white",
  high = "#FC4E07"
))

# Check plots of response variable vs others based on clusters
ggplot(dataTrain, aes(1:length(y), y, colour = factor(sub_grp))) +
  labs(
    x = "Index",
    y = "log(ClaimAmount)", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(Density, y, colour = factor(sub_grp))) +
  labs(
    x = "Density",
    y = "log(ClaimAmount)", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(log(Density), y, colour = factor(sub_grp))) +
  labs(
    x = "log(Density)",
    y = "log(ClaimAmount)", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(DriverAge, y, colour = factor(sub_grp))) +
  labs(
    x = "DriverAge",
    y = "log(ClaimAmount)", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(CarAge, y, colour = factor(sub_grp))) +
  labs(
    x = "CarAge",
    y = "log(ClaimAmount)", colour = "Cluster"
  ) +
  geom_point()

