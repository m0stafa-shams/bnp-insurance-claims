## ----------------------------------------------------------------------------
## This is a R code for the DPMM model applied to claims frequency
## data consists of its MCMC algorithm, its posterior predictive
## distribution and evaluating its clustering performance. In this
## code, we assume that we have two predictor variables and thus
## three BNP regression parameters. The definitions of variables
## used in this code are listed here:
##
## n: size of the training data set
## ncl: vector of response variables which is number of claims
## X: design matrix of predictor variables
## vn: number of columns in design matrix which is 3 here
## g0mu: mean of the multivariate normal base measure
## ----------------------------------------------------------------------------

# Load libraries
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
freq <- freMTPLfreq

nrow(freq)
head(freq)

length(which(freq$ClaimNb == 0)) / nrow(freq)

# Train Test split
ntrain <- 2000
set.seed(1099)
train <- sample(nrow(freq), ntrain)
dataTrain <- freq[train, ]
dataTest <- freq[-train, ]

# ncl: vector of response variables which is number of claims
ncl <- dataTrain$ClaimNb

# n: size of the training data set
n <- length(ncl)

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

# Laplace approximation to obtain the mean and variance of
# the multivariate normal proposal distribution for the
# Metropolis–Hastings algorithm update to draw new values
# for distinct components in the Neal's Algorithm 8
laplace_approx <- function(model, inits) {
  fit <- optim(
    par = inits, fn = model, method = "BFGS", control = list(fnscale = -1),
    hessian = TRUE
  )
  param_mean <- fit$par
  param_cov_mat <- solve(-fit$hessian)
  list(param_mean = param_mean, param_cov_mat = param_cov_mat)
}

# covM: covariance matrix of the multivariate normal base measure
covM <- diag(vn)

# g0mu: mean of the multivariate normal base measure
g0mu <- rep(0, vn)

# Base probability measure, G_0, in DPMM model
g0 <- function(s) {
  rmvnorm(s, g0mu, covM)
}

#### -------------------------------------------------------------------
#### Neal's Algorithm8 to sample from the posterior distribution of
#### the DPMM
#### -------------------------------------------------------------------

neal_algorithm8 <- function(beta0, m, alpha0, n_iter, a1, a2) {

  # The function parameters are:
  # beta0: initial value of model parameters matrix
  # m: auxiliary parameter
  # alpha0: initial value of alpha in DP prior
  # n_iter: number of MCMC algorithm iterations
  # a1: shape parameter of gamma prior for updating alpha
  # a2: rate parameter of gamma prior for updating alpha

  beta <- beta0 # posterior samples of model parameters

  # BNP regression coefficients, 'beta's
  for (i in 1:(vn)) {
    assign(paste("beta", i - 1, sep = "_"), beta0[, i])
  }

  phi <- unique(beta0) # distinct components

  ndistinct <- nrow(phi) # number of distinct components

  c <- configure(beta0) # configuration vector

  temp <- beta0 # temporary matrix

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
        prob[l] <- (log(n_i(c, l, i)) - exp(X[i, ] %*% phi[l, ]) + ncl[i] *
          (X[i, ] %*% phi[l, ]))
      }

      for (l in (k_minus + 1):h) {
        prob[l] <- (log((alpha) / m) - exp(X[i, ] %*% phi[l, ]) + ncl[i] *
          (X[i, ] %*% phi[l, ]))
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

    # Metropolis–Hastings algorithm update to draw new values
    # for distinct components
    for (j in 1:nrow(phi)) {
      abc <- (c == j)
      index <- which(abc == TRUE)

      model <- function(p) {
        log_lik <- (sum(-exp(X[index, ] %*% p) + ncl[index] * (X[index, ] %*%
          p)))
        log_post <- (log_lik - (1 / 2) * t(p - g0mu) %*% solve(covM) %*% (p -
          g0mu))
        return(log_post)
      }

      inits <- phi[j, ]

      # Laplace approximation to obtain the mean and variance of
      # the multivariate normal proposal distribution for the
      # Metropolis–Hastings algorithm
      report <- laplace_approx(inits = inits, model = model)

      mu <- report$param_mean
      sigma <- (1) * report$param_cov_mat

      proposed <- as.matrix(rmvnorm(1, mu, sigma))

      # Proposal is then accepted or rejected according to 
      # this Metropolis ratio r:
      r <- (sum(-exp(X[index, ] %*% t(proposed)) + ncl[index] * (X[index, ] %*%
        t(proposed))) - sum(-exp(X[index, ] %*% phi[j, ]) + ncl[index] *
        (X[index, ] %*% phi[j, ])) + (-1 / 2) * (proposed - g0mu) %*% solve(covM) %*%
        t(proposed - g0mu) + (-1 / 2) * t(phi[j, ] - mu) %*% solve(sigma) %*%
        (phi[j, ] - mu) - (-1 / 2) * t(phi[j, ] - g0mu) %*% solve(covM) %*%
        (phi[j, ] - g0mu) - (-1 / 2) * (proposed - mu) %*% solve(sigma) %*%
        t(proposed - mu))

      if (log(runif(1)) < r) {
        phi[j, ] <- proposed
      }
    }

    temp <- phi[c, ]
    beta <- cbind(beta, temp)

    for (i in 1:vn) {
      assign(paste("beta", i - 1, sep = "_"), rbind(get(paste("beta", i - 1,
        sep = "_"
      )), temp[, i]))
    }

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
    "beta" = beta,
    "beta_0" = beta_0, "beta_1" = beta_1, "beta_2" = beta_2,
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

beta01 <- g0(n) # initial values

stime <- system.time(b1 <- neal_algorithm8(
  beta0 = beta01, m = m, alpha0 = alpha0,
  n_iter = n_iter, a1 = a1, a2 = a2
))

beta1 <- b1$beta

for (i in 1:vn) {
  assign(paste("beta1", i - 1, sep = "_"), b1[[paste("beta", i - 1, sep = "_")]])
}

ndistinct1 <- b1$ndistinct
hist(ndistinct1[-(1:burnin)], xlab = "Clusters", main = "Histogram", breaks = 100)
k1 <- round(mean(ndistinct1[-(1:burnin)]))

precision1 <- b1$precision
plot(precision1[-(1:burnin)], type = "l", ylab = "Precision parameter")
boxplot(precision1[-(1:burnin)])

#  Diagnostic tests of convergence
mc1 <- as.mcmc(cbind(beta1_0[-(1:burnin), ], beta1_1[-(1:burnin), ], beta1_2[-(1:burnin), ]))

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

#### -------------------------------------------------------------------
#### Posterior predictive distribution
#### -------------------------------------------------------------------

pred <- function(nclnew, xnew) {

  # nclnew: new data point
  # xnew: vector of risk factors

  model <- function(p) {
    log_lik <- (-exp(xnew %*% p) + nclnew * (xnew %*% p))
    log_post <- (log_lik - (1 / 2) * t(p - g0mu) %*% solve(covM) %*% (p - g0mu))
    return(log_post)
  }

  inits <- rep(0, vn)

  report <- laplace_approx(model = model, inits = inits)

  mu <- report$param_mean
  sigma <- report$param_cov_mat

  coef1 <- mean((precision1[-(1:burnin)]) / (precision1[-(1:burnin)] + n))

  part1 <- coef1 * exp(model(mu)) * sqrt(det(sigma))

  f <- function(p) {
    (-exp(xnew %*% p) + nclnew * (xnew %*% p))
  }

  part2 <- 0

  for (t in (burnin + 1):n_iter) {
    temp_theta <- beta1[, ((vn) * t - (vn) + 1):((vn) * t)]
    for (v1 in 1:n) {
      part2 <- (part2 + (exp(f(temp_theta[v1, ])) / (precision1[t] + n)))
    }
  }

  part2 <- (part2 / (n_iter - burnin))

  result <- ((1 / factorial(nclnew)) * (part1 + part2))

  return(result = result)
}

xtest <- c(1, 13, 47)

sdx <- apply(xVar, 2, sd)
meanx <- apply(xVar, 2, mean)
xnew <- (xtest - meanx) / sdx
xnew[1] <- 1
xnew

lTest <- nrow(dataTest[dataTest$CarAge == xtest[2] & dataTest$DriverAge == xtest[3], ])
yTest <- dataTest[dataTest$CarAge == xtest[2] & dataTest$DriverAge == xtest[3], ]$ClaimNb
my <- max(yTest)

stime_par <- system.time(pr <- foreach(nclnew = 0:my, .combine = "c") %dopar% {
  pred(nclnew = nclnew, xnew = xnew)
})

sum(pr)

yTable <- table(yTest)
dyt <- as.data.frame(yTable)
prl <- (pr * lTest)[as.numeric(as.character(dyt[, 1])) + 1]
freqtable <- rbind(prl, yTable)

yHist <- hist(yTest, breaks = ((0:(max(yTest) + 1)) - 0.5), plot = F)

# Posterior predictive density distribution plots
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))

plot(0:my, pr * lTest,
  type = "h", main = "Predictive density distribution", xlab = expression(y[i]),
  ylab = "Frequency", col = "darkblue"
)

hist(yTest, main = "Testing data", xlab = expression(y[i]), xlim = c(0, my), breaks = ((0:(max(yTest) +
  1)) - 0.5), freq = T, col = "red")

barplot(freqtable,
  beside = T, xlab = expression(y[i]), ylab = "Frequency", main = "Bar Plot",
  col = c("darkblue", "red")
)
legend("topright", c("Predictive distribution", "Testing data"), col = c(
  "darkblue",
  "red"
), lwd = 5, cex = 0.5)
par(mfrow = c(1, 1))

barplot(freqtable,
  beside = T, xlab = expression(y[i]), ylab = "Frequency", main = "Bar Plot",
  col = c("darkblue", "red")
)
legend("topright", c("Predictive distribution", "Testing data"), col = c(
  "darkblue",
  "red"
), lwd = 5, cex = 0.5)

# Chi-Square of Goodness of Fit test
ind <- which(prl >= 5)
expected <- prl[ind]
if (length(expected) < length(prl)) {
  expected[length(ind) + 1] <- sum(prl[-ind])
}
obs <- as.numeric(dyt[, 2])
observed <- obs[ind]
if (length(observed) < length(expected)) {
  observed[length(ind) + 1] <- sum(obs[-ind])
}

rbind(expected, observed)

df <- length(observed) - 1
df

chi.test.stat <- sum((observed - expected)^2 / expected)
chi.test.stat

pv <- 1 - pchisq(chi.test.stat, length(observed) - 1)
pv

qchisq(p = 0.05, df = length(observed) - 1, lower.tail = FALSE)

# Mean absolute percentage error (MAPE)
mape <- sum(abs((observed - expected) / observed)) / length(observed)
mape

# Mean absolute error (MAE)
mae <- sum(abs((observed - expected))) / length(observed)
mae

# Root-mean-square deviation (RMSE)
rmse <- sqrt(sum((observed - expected)^2) / length(observed))
rmse
rmse / mean(observed)

# MSE
mean((observed - expected)^2)

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
    s <- s + mat1(b[, ((vn) * i - (vn) + 1):((vn) * i)])
  }
  # Dissimilarity Matrix
  s <- s + t(s) - (n_iter - burnin_cluster) * diag(n)
  d <- 1 - (1 / (n_iter - burnin_cluster)) * s
  return(d)
}

stime <- system.time(d1 <- dissimilarity_mat(beta1))

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
ggplot(dataTrain, aes(1:length(ncl), ncl, colour = factor(sub_grp))) +
  labs(
    x = "Index",
    y = "ClaimNumber", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(Density, ncl, colour = factor(sub_grp))) +
  labs(
    x = "Density",
    y = "ClaimNumber", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(log(Density), ncl, colour = factor(sub_grp))) +
  labs(
    x = "log(Density)",
    y = "ClaimNumber", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(DriverAge, ncl, colour = factor(sub_grp))) +
  labs(
    x = "DriverAge",
    y = "ClaimNumber", colour = "Cluster"
  ) +
  geom_point()

ggplot(dataTrain, aes(CarAge, ncl, colour = factor(sub_grp))) +
  labs(
    x = "CarAge",
    y = "ClaimNumber", colour = "Cluster"
  ) +
  geom_point()

