## ----------------------------------------------------------------------------
## This is a R code for the PYMM model applied to log(claims severity) 
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

# Base probability measure, G_0, in PYM model
g0 <- function(w) {
  s2 <- 1 / rgamma(w, shape = shape, rate = rate)
  bm <- NULL
  for (i in 1:w) {
    bm <- rbind(bm, rmvnorm(1, mean = g0mu, sigma = n0 * s2[i] * diag(vn)))
  }
  return(cbind(bm, s2))
}

# Function 'log_a_s' will be used in the process 
# of updating d and alpha
log_a_s <- function(a, s) {
  if (s == 0) {
    return(0)
  } else {
    return(sum(log(a + (0:(s - 1)))))
  }
}

# The conditional log joint distribution of the
# number of distinct components r.v. and the
# frequencies r.v. given the discount parameter d
# which is calculated in Lijoi et al. (2007) in order
# to update the discount parameter d.
d_logjoint <- function(nd, freq_temp, d, alpha) {
  if (nd == 1) {
    expr1 <- 0
  } else {
    expr1 <- sum(log(alpha + (1:(nd - 1)) * d))
  }
  expr2 <- sum(sapply(freq_temp - 1, log_a_s, a = 1 - d))
  return(expr1 + expr2)
}

# The conditional log joint distribution of the
# number of distinct components r.v. and the
# frequencies r.v. given the precision parameter alpha
# which is calculated in Lijoi et al. (2007) in order
# to update the precision parameter alpha
alpha_logjoint <- function(nd, alpha, d) {
  if (nd == 1) {
    expr1 <- 0
  } else {
    expr1 <- sum(log(alpha + (1:(nd - 1)) * d))
  }
  expr2 <- log_a_s(a = alpha + 1, s = n - 1)
  return(expr1 - expr2)
}

#### -------------------------------------------------------------------
#### Neal's Algorithm8 to sample from the posterior distribution of
#### the PYM
#### -------------------------------------------------------------------

neal_algorithm8 <- function(betasigma0, m, alpha0, d0, n_iter) {

  # The function parameters are:
  # betasigma0: initial value of model parameters matrix
  # m: auxiliary parameter
  # alpha0: initial value of alpha in PY prior
  # d0: initial value of d in PY prior
  # n_iter: number of MCMC algorithm iterations

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

  d <- d0 # discount parameter
  discount <- d0 # posterior samples of d

  # number of accepted proposals in each batch
  # of 50 iterations when updating d
  d_accept <- 0

  # logarithm of the standard deviation of the
  # proposal when updating d
  d_ls <- log(0.1)

  # number of accepted proposals in each batch
  # of 50 iterations when updating alpha
  alpha_accept <- 0

  # logarithm of the standard deviation of the
  # proposal when updating alpha
  alpha_ls <- log(0.1)

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
        prob[l] <- (log(n_i(c, l, i) - d) - (1 / 2) * log(phi[l, vn + 1]) -
          (1 / 2) * (1 / (phi[l, vn + 1])) * (y[i] - X[i, ] %*% phi[l, 1:vn])^2)
      }

      for (l in (k_minus + 1):h) {
        prob[l] <- (log((alpha + d * k_minus) / m) - (1 / 2) * log(phi[l, vn +
          1]) - (1 / 2) * (1 / (phi[l, vn + 1])) * (y[i] - X[i, ] %*% phi[
          l,
          1:vn
        ])^2)
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
      scalePost <- b + (1 / 2) * sum(y[index]^2) - (1 / 2) * dt %*% solve(M) %*%
        t(dt)


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

    # Vector of the frequencies of each distinct component
    freq_temp <- NULL
    for (var1 in 1:nd) {
      freq_temp <- c(freq_temp, sum(c == var1))
    }

    # Random-walk Metropolis–Hastings algorithm to update d
    # and also tuning of associated proposal variances using the
    # Adaptive Metropolis-within-Gibbs algorithm used in
    # Roberts and Rosenthal (2009)
    if (k %% 50 == 0) {
      d_ar <- d_accept / 50 # acceptance rate
      n_batch <- k %/% 50
      if (d_ar < 0.43) {
        d_ls <- d_ls - min(0.01, 1 / sqrt(n_batch))
      } else if (d_ar > 0.45) {
        d_ls <- d_ls + min(0.01, 1 / sqrt(n_batch))
      }
      # put 0 back for the number of accepted proposals
      d_accept <- 0
    }

    proposed <- rnorm(1, mean = d, sd = exp(d_ls))

    # Metropolis ratio r:
    if ((proposed >= 0) & (proposed < 1) & (alpha + proposed) > 0) {
      r <- (d_logjoint(nd = nd, freq_temp = freq_temp, d = proposed, alpha = alpha) -
        d_logjoint(nd = nd, freq_temp = freq_temp, d = d, alpha = alpha))

      if (log(runif(1)) < r) {
        d <- proposed
        d_accept <- d_accept + 1
      }
    }

    discount <- c(discount, d)

    # Random-walk Metropolis–Hastings algorithm to update alpha
    # and also tuning of associated proposal variances using the
    # Adaptive Metropolis-within-Gibbs algorithm used in
    # Roberts and Rosenthal (2009)
    if (k %% 50 == 0) {
      alpha_ar <- alpha_accept / 50 # acceptance rate
      n_batch <- k %/% 50
      if (alpha_ar < 0.43) {
        alpha_ls <- alpha_ls - min(0.01, 1 / sqrt(n_batch))
      } else if (alpha_ar > 0.45) {
        alpha_ls <- alpha_ls + min(0.01, 1 / sqrt(n_batch))
      }
      # put 0 back for the number of accepted proposals
      alpha_accept <- 0
    }

    proposed <- rnorm(1, mean = alpha, sd = exp(alpha_ls))

    # Metropolis ratio
    if ((proposed + d) > 0) {
      r <- (alpha_logjoint(nd = nd, alpha = proposed, d = d) + dnorm(x = log(proposed +
        d), mean = 0, sd = 1, log = TRUE) - log(proposed + d) - alpha_logjoint(
        nd = nd,
        alpha = alpha, d = d
      ) - dnorm(
        x = log(alpha + d), mean = 0, sd = 1,
        log = TRUE
      ) + log(alpha + d))

      if (log(runif(1)) < r) {
        alpha <- proposed
        alpha_accept <- alpha_accept + 1
      }
    }

    precision <- c(precision, alpha)
  }

  b <- list(
    "betasigma" = betasigma,
    "beta_0" = beta_0, "beta_1" = beta_1, "beta_2" = beta_2,
    "sigma" = sigma,
    "ndistinct" = ndistinct,
    "discount" = discount,
    "precision" = precision
  )
  return(b)
}

n_iter <- 50000 # number of MCMC algorithm iterations
burnin <- (0.5) * n_iter # burn-in
alpha0 <- 1
d0 <- 0.5
m <- 2

betasigma01 <- g0(n) # initial values

stime <- system.time(b1 <- neal_algorithm8(
  betasigma0 = betasigma01, m = m, alpha0 = alpha0,
  d0 = d0, n_iter = n_iter
))

betasigma1 <- b1$betasigma

for (i in 1:vn) {
  assign(paste("beta1", i - 1, sep = "_"), b1[[paste("beta", i - 1, sep = "_")]])
}

sigma1 <- b1$sigma

ndistinct1 <- b1$ndistinct
hist(ndistinct1[-(1:burnin)], xlab = "Clusters", main = "Histogram", breaks = 100)
k1 <- round(mean(ndistinct1[-(1:burnin)]))

discount1 <- b1$discount
plot(discount1[-(1:burnin)], type = "l", ylab = "Discount parameter")
boxplot(discount1[-(1:burnin)])

precision1 <- b1$precision
plot(precision1[-(1:burnin)], type = "l", ylab = "Precision parameter")
boxplot(precision1[-(1:burnin)])

# Diagnostic tests of convergence
mc1 <- as.mcmc(cbind(beta1_0[-(1:burnin), ], beta1_1[-(1:burnin), ], beta1_2[-(1:burnin), ], sigma1[-(1:burnin), ]))

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

  coef1 <- (mean((precision1[-(1:burnin)] + discount1[-(1:burnin)] * ndistinct1[-(1:burnin)]) *
    (1 / (precision1[-(1:burnin)] + n))))


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
    part2_a <- 0
    for (v1 in 1:n) {
      part2_a <- (part2_a + f(temp_theta[v1, ]))
    }

    temp_theta_star <- unique(betasigma1[, ((vn + 1) * t - (vn + 1) + 1):((vn +
      1) * t)])
    part2_b <- 0
    for (v2 in 1:ndistinct1[t]) {
      part2_b <- (part2_b + f(temp_theta_star[v2, ]))
    }

    part2 <- (part2 + (part2_a / (precision1[t] + n)) - (discount1[t] * part2_b / (precision1[t] +
      n)))
  }

  part2 <- (part2 / (n_iter - burnin))

  result <- (part1 + part2)

  return(result = result)
}

xtest <- c(1, 13, 47)

sdx <- apply(xVar, 2, sd)
meanx <- apply(xVar, 2, mean)
xnew <- (xtest - meanx) / sdx
xnew[1] <- 1
xnew

lTest <- nrow(dataTest[dataTest$CarAge == xtest[2] & dataTest$DriverAge == xtest[3], ])
yTest <- log(dataTest[dataTest$CarAge == xtest[2] & dataTest$DriverAge == xtest[3], ]$ClaimAmount)
my <- round(max(yTest))

stime_par <- system.time(pr <- foreach(ynew = 0:my, .combine = "c") %dopar% {
  pred(ynew = ynew, xnew = xnew)
})

yHist <- hist(yTest, breaks = ((0:(my + 1)) - 0.5), plot = F)

# Posterior predictive density distribution
hist(yTest,
  main = "Posterior predictive density distribution", xlab = "log(ClaimAmount)",
  xlim = c(0, my + 0.5), freq = FALSE, ylim = c(0, 1), breaks = ((0:(my + 1)) -
    0.5)
)
points(0:my, pr, type = "p", col = "red")
lines(0:my, pr, type = "l", col = "red")

# MSE
mean((yHist$counts - pr * lTest)^2)

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

