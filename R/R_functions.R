#By M.D. Edge, 9/19/2018.
#Functions used in the main text or exercises of
#Statistical Thinking from Scratch.

#Install (if necessary) and load all packages used in the book's code.
#if(!("MASS" %in% installed.packages())){install.packages("MASS")}
#library(MASS)

#' Simulate and plot the distribution of the sample mean
#' 
#' Simulates a specified number of samples of a specified size from the beta distribution and plots the distribution of sample means with the density of an approximating normal distribution overlaid.
#' @param samp.size The size of each simulated sample.
#' @param n.samps The number of samples to simulate.
#' @param shape1 The first parameter of the beta distribution from which the samples will be drawn.
#' @param shape2 The second parameter of the beta distribution.
#' @keywords sample mean, law of large numbers, central limit theorem, distribution of the sample mean
#' @export
#' @examples 
#' dosm.beta.hist(25, 1000, 2, 2) 
dosm.beta.hist <- function(samp.size, n.samps, shape1 = 1, shape2 = 1, ...){
  samps <- rbeta(samp.size*n.samps, shape1, shape2)
  sim.mat <- matrix(samps, nrow = n.samps)
  dosm <- rowMeans(sim.mat)
  hist(dosm, freq=FALSE, ...)#freq=FALSE re-scales the height
  x <- seq(0,1,length.out = 1000)
  lines(x, dnorm(x, mean = mean(dosm), sd = sd(dosm)))
}


#' Sample from a Pareto distribution
#' 
#' Simulates a sample of Pareto-distributed data.
#' @param n The number of independent observations to draw (i.e. sample size).
#' @param a The shape parameter of the Pareto distribution.
#' @param b The scale parameter of the Pareto distribution.
#' @return A vector of Pareto-distributed random numbers of length n.
#' @keywords simulated data, Pareto
#' @export
#' @examples 
#' rpareto(100, 0.5, 1)
rpareto <- function(n, a=0.5, b=1){
  qpareto <- function(u, a=0.5, b=1){b/(1-u)^(1/a)}
  qpareto(runif(n),a,b)
}



#' Compare tails of a sample to normal expectation.
#' 
#' Compares the proportion of observations in a vector beyond k standard deviations to the expectation under a normal distribution. The output is a ratio of the proportion of data beyond k standard deviations from the expectation divided by the Normal-theory probability of observations beyond k standard deviations from the expectation.
#' @param x A numeric vector
#' @param k The number of standard deviations beyond which is considered the "tail" (for the purpose of the computation)
#' @param expec The expectation of a normal distribution to which x will be compared
#' @param st.d The standard deviation of a normal distribution to which x will be compared
#' @return A ratio of the proportion of data beyond k standard deviations from the expectation divided by the Normal-theory probability of observations beyond k standard deviations from the expectation.
#' @keywords extreme observations, central limit theorem
#' @export
#' @examples 
#' compare.tail.to.normal(rnorm(10000, 0, 1), 3, 0, 1)
compare.tail.to.normal <- function(x, k, expec, st.d){
  mean(x < (expec - k* st.d) | x > (expec + k* st.d))/(1 - (pnorm(k) - pnorm(-k)))
}

#' Simulate data from a simple linear model.
#' 
#' Generates one set of n independent pairs of observations using a simple linear model (y = a + b*x + disturbance) with intercept a, slope b, and (by default) normally distributed disturbances. The output is a matrix with x in the first column and y in the second column. The output matrix is sorted by the x values, with the lowest x value in the first row.
#' @param a The intercept parameter of the linear regression model to be simulated.
#' @param b The slope parameter of the linear regression model to be simulated.
#' @param var.eps The variance of the disturbances in the model y = a + b*x + disturbance.
#' @param n The number of pairs of observations to simulate.
#' @param mu.x The expectation of the x values.
#' @param var.x The variance of the x values.
#' @param rx A function for drawing the random x values. rnorm() is the default, which makes x normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @param rdist A function for drawing the random disturbances. rnorm() is the default, which makes the disturbances normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @return A matrix with n rows, x in the first column and y in the second column. The output matrix is sorted by the x values, with the lowest x value in the first row.
#' @keywords Simulation, simple linear regression
#' @export
#' @examples 
#' sim.lm(a = 3, b = 1/2, n = 11)
sim.lm <- function(a, b, var.eps = 1, n = 50, mu.x = 8, var.x = 4, rx = rnorm, rdist = rnorm){
  x <- sort(rx(n, mu.x, sqrt(var.x)))
  disturbs <- rdist(n, 0, sqrt(var.eps))
  y <- a + b*x + disturbs
  cbind(x,y)
}


#' Generate a matrix of samples from a normal distribution.
#' 
#' Draws normal samples and formats them into a matrix, where each row contains a sample.
#' @param mu The expectation of the normal distribution from which to draw samples.
#' @param sigma The standard deviation of the normal distribution from which to draw samples.
#' @param n The number of independent observations to include in each sample.
#' @param nsamps The number of samples to generate.
#' @return A matrix of independent normally distributed random numbers with nsamps rows and n columns.
#' @keywords simulation, normal distribution
#' @export
#' @examples 
#' norm.samps(10, 1, 5, 8)
norm.samps <- function(mu = 0, sigma = 1, n = 25, nsamps = 10000){
  samps <- rnorm(n*nsamps, mu, sigma)
  matrix(samps, nrow = nsamps, ncol = n)
}


#' Sample from a Laplace (also called double-exponential) distribution
#' 
#' Simulates a sample of Laplace-distributed data.
#' @param n The number of independent observations to draw (i.e. sample size).
#' @param mean The expectation of the Laplace distribution.
#' @param sd The standard deviation of the Laplace distribution.
#' @return A vector of Laplace-distributed random numbers of length n.
#' @keywords simulated data, Laplace, double-exponential
#' @export
#' @examples 
#' rlaplace(10, 0, 1)
rlaplace <- function(n, mean = 0, sd = 1){
  exp1 <- rexp(2*n, 1)
  x <- exp1[1:n] - exp1[(n+1):(2*n)]
  x * sd/sqrt(2) + mean
}

#' Generate a matrix of samples from a Laplace distribution.
#' 
#' Draws Laplace samples and formats them into a matrix, where each row contains a sample.
#' @param mu The expectation of the Laplace distribution from which to draw samples.
#' @param sigma The standard deviation of the Laplace distribution from which to draw samples.
#' @param n The number of independent observations to include in each sample.
#' @param nsamps The number of samples to generate.
#' @return A matrix of independent Laplace-distributed random numbers with nsamps rows and n columns.
#' @keywords simulation, Laplace distribution
#' @export
#' @examples 
#' laplace.samps(10, 1, 5, 8)
laplace.samps <- function(mu = 0, sigma = 1, n = 25, nsamps = 10000){
  samps <- rlaplace(n*nsamps, mu, sigma)
  matrix(samps, nrow = nsamps, ncol = n)
}

#' A matrix of samples from a Normal distribution with outliers
#' 
#' Simulates samples of normally distributed data, each "contaminated" with outliers from a different normal distribution according to a user-specified probability. 
#' @param nsamps The number of samples to simulate.
#' @param n The number of independent observations per sample
#' @param gamma The probability that each observation is from the non-target distribution.
#' @param theta The expectation of the target distribution.
#' @param lambda The expectation of the non-target / "contamination" / outlier distribution.
#' @param targ.sd The standard deviation of the target distribution.
#' @param out.sd The standard deviation of the outlier distribution.
#' @return A matrix of independent random numbers with nsamps rows and n columns.
#' @keywords simulated data, normal distribution, outliers
#' @export
#' @examples 
#' rnorm.out(5, 20, .2, 0, 10)
rnorm.out <- function(nsamps, n, gamma, theta = 0, lambda = 10, targ.sd = 1, out.sd = 1){
  n.target <- rbinom(1, n*nsamps, 1-gamma)
  n.nontarget <- n*nsamps - n.target
  samp.target <- rnorm(n.target, theta, targ.sd)
  samp.nontarget <- rnorm(n.nontarget, lambda, out.sd)
  allsamps <- sample(c(samp.target, samp.nontarget))
  matrix(allsamps, nrow = n, ncol = nsamps)
}


#' Simulate data from a linear model and estimate the parameters.
#' 
#' Simulates data from a linear model using the sim.lm() function, each time estimating the parameters using a user-supplied function.
#' @param nsims The number of samples to simulate.
#' @param n The number of pairs of observations to simulate per sample.
#' @param a The intercept parameter of the linear regression model to be simulated.
#' @param b The slope parameter of the linear regression model to be simulated.
#' @param var.eps The variance of the disturbances in the model y = a + b*x + disturbance.
#' @param mu.x The expectation of the x values.
#' @param var.x The variance of the x values.
#' @param rx A function for drawing the random x values. rnorm() is the default, which makes x normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @param rdist A function for drawing the random disturbances. rnorm() is the default, which makes the disturbances normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @param estfun A function for estimating the parameters. lm() is the default, which gives least-squares estimates. Other functions can be supplied so long as they accept the same formula-style input in the first argument and return a list that has the coefficients stored in a vector of length 2 accessible by $coef. rq() in the quantreg package is another function that works.
#' @return A matrix with 2 columns and nsamps rows. Each row contains estimation results from a different sample, with the intercept estimate in the first column and the slope estimate in the second column.
#' @keywords simulated data, simple linear regression
#' @export
#' @examples 
#' sim.lm.ests(10)
#' sim.lm.ests(10, estfun = rq) #if package quantreg is loaded.
sim.lm.ests <- function(nsims, n = 50, a = 0, b = 0, var.eps = 1, mu.x = 0, var.x = 1, rx = rnorm, rdist = rnorm, estfun = lm){
  ests <- matrix(nrow = nsims, ncol = 2)
  for(i in 1:nsims){
    dat <- sim.lm(n = n, a = a, b = b, var.eps = var.eps, mu.x = 		mu.x, var.x = var.x, rx = rx, rdist = rdist)
    ests[i,] <- estfun(dat[,2] ~ dat[,1])$coef
  }
  ests
}


#' Sample from a Normal distribution with outliers
#' 
#' Simulates a samples of normally distributed data, each "contaminated" with outliers from a different normal distribution according to a user-specified probability. Similar to rnorm.out, except that the order of arguments is different and the output is formatted as a vector.
#' @param n The number of independent observations per sample
#' @param mean The expectation of the target distribution.
#' @param sd The standard deviation of the target distribution.
#' @param contam.p The probability that each observation is from the non-target distribution.
#' @param contam.mean The expectation of the non-target / "contamination" / outlier distribution.
#' @param contam.sd The standard deviation of the outlier distribution.
#' @return A vector of independent random numbers n entries.
#' @keywords simulated data, normal distribution, outliers
#' @export
#' @examples 
#' rnorm.mix(20, 0, 1, .2, 10, 1)
rnorm.mix <- function(n, mean = 0, sd = 1, contam.p = 0.01, contam.mean = -5, contam.sd = 0){
  ncontam <- rbinom(1, n, contam.p)
  c(rnorm(n - ncontam, mean, sd), rnorm(ncontam, contam.mean, contam.sd))
}


#' Simulate a multiple-testing scenario
#' 
#' Simulates type I error inflation from multiple testing for a test comparing two group means. The groups are compared on n.meas outcome measurements, each mutually correlated at a user-specified level. A matrix of p values is returned with n.meas columns and n.sims rows.
#' @param n The number of people measured in each group per simulation.
#' @param n.meas The number of outcome measurements taken per person.
#' @param correl The correlation of the (normally distributed) outcome measurements---all pairs of measurements have the same correlation.
#' @param n.sims The number of simulated samples to run.
#' @return A matrix of p values with n.meas columns and n.sims rows. Each row corresponds to a different sample, and the p values from the respective outcome measures are in distinct columns.
#' @keywords simulated data, multiple testing, p hacking
#' @export
#' @examples 
#' many.outcome.sim(100, 3, .5, 10)
many.outcome.sim <- function(n, n.meas, correl, n.sims){
  #A function for computing a p-value comparing the means
  #of two groups whose scores are entered in a vector.
  #Assumes an input vector with an even number of entries,
  #the first half of which are in one group, with the second
  #half in a different group..
  t.p <- function(x){
    n.tot <- length(x)
    n.1 <- n.tot/2
    x1 <- x[1:n.1]
    x2 <- x[(n.1+1):n.tot]
    t.st <- (mean(x1) - mean(x2))/sqrt((sd(x1)^2 + sd(x2)^2)/n.1)
    2*pt(-abs(t.st), n.tot - 2)
  }
  require(MASS)
  #Specify a covariance matrix for mvrnorm().
  sigma <- matrix(correl, nrow = n.meas, ncol = n.meas)
  diag(sigma) <- 1
  #Initialize a matrix for holding p values.
  sim.ps <- matrix(-9, nrow = n.sims, ncol = n.meas)
  #simulate all the necessary observations.
  sim.dat <- MASS::mvrnorm(n.sims*n*2, mu = rep(0, n.meas), Sigma = sigma)
  ps.mat <- matrix(NA, nrow = n.sims, ncol = n.meas)
  #For each measurement, conduct hypothesis tests for each sample.
  for(i in 1:n.meas){
    test.dat <- matrix(sim.dat[,i], nrow = n.sims)
    ps.mat[,i] <- apply(test.dat, 1, t.p)
  }
  ps.mat
}

#' Simulate an optional stopping p-value scenario
#' 
#' Simulates type I error inflation from testing a difference in group means many times as the trial is ongoing, stopping if the test is significant. A matrix of p values is returned with length(ns) columns and n.sims rows.
#' @param ns A vector of numbers, each of which is a number of people per group at which the group difference should be tested.
#' @param n.sims The number of simulated samples to run.
#' @return A matrix of p values with length(ns) columns and n.sims rows. Each row corresponds to a different sample, and the p values resulting from tests of the nested subsamples of size specified by ns are in different columns.
#' @keywords simulated data, p hacking
#' @export
#' @examples 
#' sims <- serial.testing.sim(n.sims = 100)
#' mean(apply(sims, 1, min) < .05) #How often would optional stopping lead to a rejection of the null at the .05 level?
serial.testing.sim <- function(ns = c(20, 30, 40, 50), n.sims = 10000){
  #A function for computing a p-value comparing the means
  #of two groups whose scores are entered in a vector.
  #Assumes an input vector with an even number of entries,
  #the first half of which are in one group, with the second
  #half in a different group..
  t.p <- function(x){
    n.tot <- length(x)
    n.1 <- n.tot/2
    x1 <- x[1:n.1]
    x2 <- x[(n.1+1):n.tot]
    t.st <- (mean(x1) - mean(x2))/sqrt((sd(x1)^2 + sd(x2)^2)/n.1)
    2*pt(-abs(t.st), n.tot - 2)
  }
  #Initialize a matrix for holding p values.
  sim.ps <- matrix(-9, nrow = n.sims, ncol = length(ns))
  #simulate all the necessary observations.
  sim.dat <- matrix(rnorm(n.sims*max(ns)*2, 0, 1), nrow = n.sims)
  ps.mat <- matrix(NA, nrow = n.sims, ncol = length(ns))
  max.n <- max(ns)
  #For each measurement, conduct hypothesis tests for each sample.
  for(i in 1:length(ns)){
    test.dat <- sim.dat[,c(1:ns[i], (max.n + 1):(max.n + ns[i]))]
    ps.mat[,i] <- apply(test.dat, 1, t.p)
  }
  ps.mat
}


#' Estimate the power of a one-sample z test by simulation
#' 
#' Simulates nsims normal samples with expectations differing from a value specified by the null hypothesis by d (known) standard deviations. For each sample, the null hypothesis is tested by a one-sample z test, and the p value is compared with the specified significance level.
#' @param d An effect size: the number of (known) standard deviations by which the expectation of the sampled data differs from a hypothesized expectation.
#' @param n The size of each simulated sample.
#' @param level The significance level of the one-sample z-test.
#' @param nsim The number of simulated samples to run.
#' @return The proportion of p values less than the specified level (i.e. the power).
#' @keywords simulation, power, z test
#' @export
#' @examples 
#' ps.1sz(d = 0.5, n = 20)
ps.1sz <- function(d, n, level = 0.05, nsim = 10000){
  simmat <- matrix(rnorm(n*nsim, d, 1), nrow = nsim)
  samp.means <- rowMeans(simmat)
  neg.devs <- -abs(samp.means)
  ps <- 2*pnorm(neg.devs, 0, 1/sqrt(n))
  mean(ps < level)
}


#' Simulate the winner's curse with a one-sample z-test
#' 
#' Simulates the "winner's curse"---the effect that especially in low-power situations, estimates that result in signficant tests for the null hypothesis of theta = 0 tend also to be overestimates. In each simulation, a single normal sample is drawn and the hypothesis that the true effect size is zero is tested(by one-sample t-test). The output is a named vector with a a true effect size (measured in number of standard deviations from the value under the null hypothesis), the "estimated" effect size, which is the mean effect size from simulations that produced significant results, and the power, which is the proportion of simulations that achieved significance. A histogram is also produced, with all the estimated effect sizes shown, and the estimates associated with significant tests colored grey.
#' @param d The true effect size in the simulations.
#' @param n The size of each simulated sample.
#' @param lev The significance level of the one-sample z-test. Effect size estimates are averaged in the "estimated d" part of the output only if the one-sample z test produces a p value less than the level.
#' @param nsim The number of simulated samples to run.
#' @param abs.vals controls whether we pay attention to the sign of the estimate (if FALSE, the default, we do).
#' @param br the breaks parameter for the histogram.
#' @return a named vector with a a true effect size (measured in number of standard deviations from the value under the null hypothesis), the "estimated" effect size, which is the mean effect size from simulations that produced significant results, and the power, which is the proportion of simulations that achieved significance. A histogram is also produced, with all the estimated effect sizes shown, and the estimates associated with significant tests colored grey.
#' @keywords simulation, power, z test, winner's curse
#' @export
#' @examples 
#' wc.1sz(d = 0.5, n = 20, nsim = 1000)
wc.1sz <- function(d, n, lev = 0.05, nsim = 10000, abs.vals = FALSE, br = 50){
  samps <- rnorm(n*nsim, d, 1)
  simmat <- matrix(samps, nrow = nsim)
  samp.means <- rowMeans(simmat)
  neg.devs <- -abs(samp.means)
  ps <- 2*pnorm(neg.devs, 0, 1/sqrt(n))
  power <- mean(ps < lev)
  if(abs.vals == TRUE){
    ests.out <- abs(samp.means)
  }
  if(abs.vals == FALSE){
    ests.out <- samp.means
  }
  est.d <- mean(ests.out[ps < lev])
  cut.xbar.n <- qnorm(lev/2, 0, 1/sqrt(n))
  h <- hist(ests.out, breaks = br)
  cuts <- cut(h$breaks, c(-Inf, cut.xbar.n, -cut.xbar.n, Inf))
  plot(h, col = c("grey", "white", "grey")[cuts], xlab = "Estimated effect sizes", main = "")
  return(c("true d" = d, "estimated d" = est.d, "power" = power))
}




#' Draw a bootstrap sample from a matrix or vector
#' 
#' Given a matrix with n rows (or a vector of length n), samples n rows with replacement from the original matrix.
#' @param x A matrix or vector.
#' @return A matrix containing a bootstrap sample of the original matrix or vector.
#' @keywords Bootstrapping
#' @export
#' @examples 
#' boot.samp(rnorm(10))
boot.samp <- function(x){
  #If x is a vector, convert it to a matrix with one column.
  if(is.null(dim(x))){
    x <- matrix(x, ncol = 1)
  }
  n <- nrow(x)
  boot.inds <- sample(1:n, replace = TRUE)
  x[boot.inds,]
}

#' Compute the least-squares slope estimate for simple linear regression.
#' 
#' Given two numeric vectors, computes a least-squares slope estimate for the model y = a + b*x + disturbance. lm() also returns this number, but beta.mm() is faster, which is advantageous for simulating bootstrapping or permutation testing.
#' @param x A numeric vector containing entries for the independent variable.
#' @param y A numeric vector containing entries for the dependent variable, must be the same length as x.
#' @return The least-squares slope estimate (i.e. the estimate of b) for the model y = a + b*x + disturbance.
#' @keywords least-squares, simple linear regression, method of moments
#' @export
#' @examples 
#' beta.mm(rnorm(10), rnorm(10))
beta.mm <- function(x, y){
  n <- length(x)
  (sum(x*y) - (1/n)*sum(x)*sum(y)) / 
    (sum(x^2) - (1/n)*sum(x)^2)
}


#' Randomly permute the columns of a matrix independently
#' 
#' Given a matrix, return a version with its columns randomly and independently permuted.
#' @param x A matrix.
#' @return A matrix resulting from randomly permuting the columns of x independently.
#' @keywords permutation
#' @export
#' @examples 
#' perm.samp(matrix(1:25, nrow = 5))
perm.samp <- function(x){
  apply(x, 2, sample)
}




#' Simulate data from a a linear model with an unobserved confounder.
#' 
#' Generates nsims samples of n independent pairs of observations using a linear model y = alpha + beta*x + gamma*z + disturbance where z is unobserved. x and z are jointly normally distributed, and the disturbances are also normal.
#' @param n The number of independent pairs of observations to generate per simulation.
#' @param nsims The number of simulations to run.
#' @param alpha The intercept.
#' @param beta The slope for the observed independent variable.
#' @param gamma The slope for the unobserved confounder.
#' @param d.sd The standard deviation of the disturbances.
#' @param rho The correlation of (observed) x and (unobserved) z.
#' @return A list of two matrices, each with nsims rows and n columns. The first matrix contains the x values from each simulation, with one simuluation in each row. The second matrix contains the y values in the same configuration.
#' @keywords Simulation, simple linear regression, confounding, omitted variable bias.
#' @export
#' @examples 
#' sim.2var(10, 50, 3, 1/2, 1/5, 1, .5)
sim.2var <- function(n, nsims, alpha, beta, gamma, d.sd, rho){
  require(MASS)
  sig <- matrix(c(1, rho, rho, 1), nrow = 2)
  ivs <- MASS::mvrnorm(n*nsims, mu = c(0,0), sig)
  x <- ivs[,1]
  z <- ivs[,2]
  disturb <- rnorm(n*nsims, 0, d.sd)
  y <- alpha + beta * x + gamma * z + disturb
  xmat <- matrix(x, nrow = nsims)
  ymat <- matrix(y, nrow = nsims)
  cbind(xmat, ymat)
}

#' Generate a bootstrap distribution for a statistic computed on a vector.
#' 
#' Given data x, generates a bootstrap distribution with B bootstrap samples, computing a statistic speficied by FUN each time.
#' @param x A numeric vector of data.
#' @param B The number of bootstrap samples to draw.
#' @param FUN A function that computes the statistic of interest. FUN must have a vector of data as its first argument and must return a numeric vector of length 1 (i.e. a scalar).
#' @param ... Additional arguments passed to FUN.
#' @return A numeric vector containing the bootstrap distribution of the statistic specified by FUN.
#' @keywords Bootstrapping, Bootstrap distribution
#' @export
#' @examples 
#' boot.dist.1d(rnorm(50, 0, 1), 1000, median)
boot.dist.1d <- function(x, B, FUN = mean, ...){
  boot.samps <- sample(x, length(x)*B, replace = TRUE)
  boot.mat <- matrix(boot.samps, ncol = B)
  apply(boot.mat, 2, FUN, ...)
}

#' Simulate a normal sample and compute a bootstrap distribution
#' 
#' Wraps boot.dist.1d, generating a normal sample and computing a bootstrap distribution for a statistic of interest. Plots a histogram of the bootstrap distribution and returns a list containing the mean and standard deviation of the bootstrap distribution.
#' @param n The size of the normal sample to draw.
#' @param B The number of bootstrap samples to draw from the normal sample.
#' @param mu The expectation of the normal distribution from which the data are drawn.
#' @param sigma The standard deviation of the normal distribution from which the data are drawn.
#' @param FUN A function that computes the statistic of interest. FUN must have a vector of data as its first argument and must return a numeric vector of length 1 (i.e. a scalar).
#' @param ... Additional arguments passed to FUN.
#' @return A list containing the mean and standard deviation of the bootstrap distribution.
#' @keywords Bootstrapping, Bootstrap distribution
#' @export
#' @examples 
#' wrap.bm(50, 1000, FUN = median)
wrap.bm <- function(n, B, mu = 0, sigma = 1, FUN = mean, ...){
  sim <- rnorm(n, mu, sigma)
  boots <- boot.dist.1d(sim, B, FUN = FUN, ...)
  hist(boots, main = "", xlab = "Bootstrap distribution")
  list("boot m"= mean(boots), "boot se"= sd(boots))
}

#' Simulate simple linear regression data and apply a permutation test
#' 
#' Simulates data from a linear model using the sim.lm() function, each time applying a permutation test to test the null hypothesis that the slope is zero.
#' @param a The intercept parameter of the linear regression model to be simulated.
#' @param b The slope parameter of the linear regression model to be simulated.
#' @param n.perm The number of permutations to use per simulation to conduct the permutation test.
#' @param n.sim The number of simulations to conduct.
#' @param var.eps The variance of the disturbances in the model y = a + b*x + disturbance.
#' @param n The numer of independent pairs of observations to simulate per sample.
#' @param mu.x The expectation of the x values.
#' @param var.x The variance of the x values.
#' @param rdist A function for drawing the random disturbances. rnorm() is the default, which makes the disturbances normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @param rx A function for drawing the random x values. rnorm() is the default, which makes x normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @return A vector of p values of length nsims, one from each permutation test of the hypothesis that the slope is zero.
#' @keywords simulated data, simple linear regression, permutation test
#' @export
#' @examples 
#' sim.perm.B(3, .1, n.sim = 20)
sim.perm.B <- function(a, b, n.perm = 500, n.sim = 500, var.eps = 1, 
                           n = 50, mu.x = 8, var.x = 4, rdist = rnorm, rx = rnorm){
  #This version differs from the one in the text and is about twice as fast.
  #simulate all data
  dat <- sim.lm(a, b, var.eps, n*n.sim, mu.x, var.x, rdist = rdist, rx = rx)[sample(1:(n*n.sim)),]
  #reorganize for beta.mm.vec
  xs <- matrix(dat[,1], nrow = n)
  ys <- matrix(dat[,2], nrow = n)
  dat.all <- rbind(xs, ys)
  #compute all the beta estimates
  beta.mm.vec <- function(z){
    n <- length(z)/2
    x <- z[1:n]
    y <- z[(n+1):(2*n)]
    (sum(x*y) - (1/n)*sum(x)*sum(y)) / 
      (sum(x^2) - (1/n)*sum(x)^2)
  }
  betas <- apply(dat.all, 2, beta.mm.vec)
  #for each permutation, permute the y values and recompute
  #beta estimates
  perm.second.half <- function(x){
    n <- length(x)/2
    c(x[1:n], sample(x[(n+1):(2*n)]))
  }
  perm.dists <- matrix(nrow = n.perm, ncol = n.sim)
  for(i in 1:n.perm){
    perm.all <- apply(dat.all, 2, perm.second.half)
    perm.dists[i,] <- apply(perm.all, 2, beta.mm.vec)
  }
  #compare the betas to their permutation distributions to
  #arrive at p values.
  getp <- function(x){
    mean(abs(x[2:length(x)]) > abs(x[1]))
  }
  apply(rbind(matrix(betas, nrow = 1), perm.dists), 2, getp)
}



#' Simulate simple linear regression data and apply a Wald test
#' 
#' Simulates data from a linear model using the sim.lm() function, each time applying a Wald test test the null hypothesis that the slope = B0.
#' @param a The intercept parameter of the linear regression model to be simulated.
#' @param b The slope parameter of the linear regression model to be simulated.
#' @param B0 The value of the slope under the null hypothesis to be tested.
#' @param n.sim The number of simulations to conduct.
#' @param var.eps The variance of the disturbances in the model y = a + b*x + disturbance.
#' @param n The numer of independent pairs of observations to simulate per sample.
#' @param mu.x The expectation of the x values.
#' @param var.x The variance of the x values.
#' @param rdist A function for drawing the random disturbances. rnorm() is the default, which makes the disturbances normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @param rx A function for drawing the random x values. rnorm() is the default, which makes x normally distributed, but you can use any function for random number generation with first argument the sample size, second argument the expectation, and third argument the standard deviation.
#' @param pfun A cumulative distribution function to which to compare the Wald statistic.
#' @param ... Additional arguments to pfun
#' @return A vector of p values of length nsims, one from each Wald test.
#' @keywords simulated data, simple linear regression, Wald
#' @export
#' @examples 
#' sim.Wald.B(3, .1, n.sim = 20)
sim.Wald.B <- function(a, b, B0 = 0, n.sim = 1000, var.eps = 1, n = 50, 
                       mu.x = 8, var.x = 4, rdist = rnorm, rx = rnorm, pfun = pnorm, ...){
  #Initialize variables.
  ps <- numeric(n.sim)
  for(i in 1:n.sim){
    #Simulate data and compute p value.
    dat <- sim.lm(a, b, var.eps, n, mu.x, var.x, rdist = rdist, rx = rx)    
    x <- dat[,1]
    y <- dat[,2]
    #compute MLEs of beta and alpha
    B.hat <- (sum(x*y)-sum(x)*sum(y)/n)/( sum(x^2) - sum(x)^2/n)
    A.hat <- (sum(y) - B.hat*sum(x))/n
    #Compute estimated variance of MLE of beta
    vhat.dists <- sum((y - A.hat - B.hat*x)^2)/(n-2)
    vhat.Bhat <- vhat.dists/sum((x - mean(x))^2)
    #Wald statistic
    wald <- (B.hat - B0)/sqrt(vhat.Bhat)
    ps[i] <- 2*pfun(-abs(wald), ...)
  }
  #Return the p values
  return(ps)
}


#' Compute normal conjugate posterior
#' 
#' Given a sample of data from a normal distribution with known variance / standard deviation, and a normal prior for the expectation of the data, compute the parameters of the posterior for the expectation.
#' @param z A numeric vector (assumed to be drawn from a normal distribution with known variance / standard deviation)
#' @param known.sd The known standard deviation of the data distribution.
#' @param prior.mean The expectation of the normal prior.
#' @param prior.sd The standard devation of the normal prior.
#' @return A list containing the mean and standard deviation of the conjugate posterior distribution.
#' @keywords Conjugate prior, Normal, posterior.
#' @export
#' @examples 
#' post.conj.norm.norm(rnorm(100), 1, 2, 4)
post.conj.norm.norm <- function(z, known.sd, prior.mean, prior.sd){
  xbar <- mean(z)
  post.expec <- (prior.mean / prior.sd^2 + xbar*length(z) / 
                   known.sd^2)/(1 /   prior.sd^2 + length(z) / known.sd^2)
  post.var <- 1 / (1 /   prior.sd^2 + length(z) / known.sd^2)
  list("posterior.expectation" = post.expec, "posterior.variance" = post.var)
}


#' Perform rejection sampling for normally distributed data with a normal prior on the expectation.
#' 
#' Given a sample of data from a normal distribution with known variance / standard deviation, and a normal prior for the expectation of the data, draw samples from the posterior by rejection sampling.
#' @param z A numeric vector (assumed to be drawn from a normal distribution with known variance / standard deviation)
#' @param known.sd The known standard deviation of the data distribution.
#' @param prior.mn The expectation of the normal prior.
#' @param prior.sd The standard devation of the normal prior.
#' @param nsamps The number of observations to draw from the posterior by rejection sampling.
#' @return A numeric vector of length nsamps containing simulated data from the posterior distribution.
#' @keywords Rejection sampling, normal, posterior.
#' @export
#' @examples 
#' reject.samp.norm(rnorm(100), 1, 2, 4)
reject.samp.norm <- function(z, known.sd = 1, prior.mn = 0, prior.sd = 1, nsamps = 10000){
  #Get 1 sample under rejection sampling from normal with known sd.
  #Prior is a normal.
  #z is a vector of data.
  get.1.samp.norm <- function(z, known.sd = 1, prior.mn = 0, prior.sd = 1){
    accepted <- FALSE
    max.like <- exp(sum(log(dnorm(z, mean = mean(z), sd = known.sd))))
    while(accepted == FALSE){
      cand <- rnorm(1, prior.mn, prior.sd)
      like <- exp(sum(log(dnorm(z, mean = cand, sd = known.sd))))
      crit <- like / max.like
      xunif <- runif(1,0,1)
      if(xunif <= crit){accepted <- TRUE}
    }
    cand
  }
  samps <- numeric(nsamps)
  for(i in seq_along(samps)){
    samps[i] <- get.1.samp.norm(z, known.sd, prior.mn, prior.sd)
  }
  samps
}

