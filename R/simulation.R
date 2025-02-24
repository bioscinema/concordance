#' Simulate Microbial Sequence Data with Estimated Hyper-parameters
#'
#' This function simulates microbial sequence count data based on a set of estimated Dirichlet parameters.
#' It generates simulated OTU counts for a specified number of samples (\code{nSam}) and OTUs (\code{nOTU}),
#' incorporating differential abundance according to the specified proportion and mode. Effects of a covariate
#' (and optionally a confounder) are modeled, sequencing depth is simulated using a negative binomial distribution,
#' and the final OTU counts are obtained via multinomial sampling.
#'
#' @param para A list of estimated hyper-parameters. It is expected that \code{para$ref.otu.tab} contains a reference OTU table.
#' @param nSam Integer. Number of samples to simulate. Default is 100.
#' @param nOTU Integer. Number of OTUs to simulate. Default is 500.
#' @param diff.otu.pct Numeric. Proportion of OTUs that are differentially abundant. Default is 0.1.
#' @param diff.otu.direct Character. Direction of differential OTU effects. Must be either \code{"balanced"} or \code{"unbalanced"}.
#' @param diff.otu.mode Character. Mode for differential OTU effects. Options include \code{"abundant"}, \code{"rare"}, \code{"mix"}, or \code{"user_specified"}.
#' @param user_specified_otu Optional. A vector of OTU names or indices to be set as differentially abundant when \code{diff.otu.mode} is \code{"user_specified"}.
#' @param covariate.type Character. Type of covariate to simulate; either \code{"binary"} or \code{"continuous"}.
#' @param grp.ratio Numeric. For binary covariates, the ratio between groups. Default is 1.
#' @param covariate.eff.mean Numeric. Mean effect size for the covariate. Default is 1.
#' @param covariate.eff.sd Numeric. Standard deviation of the covariate effect. Default is 0.
#' @param confounder.type Character. Type of confounder to simulate. Options are \code{"none"}, \code{"binary"}, \code{"continuous"}, or \code{"both"}.
#' @param conf.cov.cor Numeric. Correlation between the covariate and the confounder. Default is 0.6.
#' @param conf.diff.otu.pct Numeric. Proportion of differential OTUs affected by the confounder. Default is 0.
#' @param conf.nondiff.otu.pct Numeric. Proportion of non-differential OTUs affected by the confounder. Default is 0.1.
#' @param confounder.eff.mean Numeric. Mean effect size for the confounder. Default is 0.
#' @param confounder.eff.sd Numeric. Standard deviation of the confounder effect. Default is 0.
#' @param error.sd Numeric. Standard deviation for the error term. Default is 0.
#' @param depth.mu Numeric. Mean sequencing depth. Default is 10000.
#' @param depth.theta Numeric. Theta parameter for the negative binomial distribution modeling sequencing depth. Default is 5.
#' @param depth.conf.factor Numeric. Factor for the covariate effect on sequencing depth. Default is 0.
#' @param cont.conf Numeric vector. Continuous confounder values. Required if \code{confounder.type} is \code{"continuous"} or \code{"both"}.
#' @param epsilon Numeric vector. Error values for the covariate simulation; typically generated with \code{rnorm(nSam)}.
#'
#' @return A list containing:
#'   \describe{
#'     \item{otu.tab.sim}{A matrix of simulated OTU counts (rows correspond to OTUs, columns to samples).}
#'     \item{covariate}{A matrix or vector of simulated covariate values for each sample.}
#'     \item{confounder}{A matrix of simulated confounder values (if applicable), or \code{NULL} if not used.}
#'     \item{diff.otu.ind}{A logical vector indicating which OTUs are truly differential.}
#'     \item{otu.names}{A character vector of OTU names used in the simulation.}
#'     \item{conf.otu.ind}{A logical vector indicating which OTUs are affected by the confounder.}
#'   }
#'
#' @details
#' The simulation process involves the following steps:
#' \enumerate{
#'   \item Subsetting the reference OTU table to the specified numbers of OTUs (\code{nOTU}) and samples (\code{nSam}).
#'   \item Simulating a confounding variable based on the specified \code{confounder.type}.
#'   \item Generating a covariate with a given correlation (\code{conf.cov.cor}) with the confounder and incorporating error (\code{epsilon}).
#'   \item Introducing differential abundance among OTUs according to the chosen direction and mode.
#'   \item Simulating sequencing depth using a negative binomial distribution, followed by generation of
#'         absolute OTU counts via multinomial sampling.
#' }
#'
#' @export
SimulateMSeqU <- function (para = EstPara, nSam = 100, nOTU = 500, diff.otu.pct = 0.1,
                           diff.otu.direct = c("balanced", "unbalanced"),
                           diff.otu.mode = c("abundant", "rare", "mix", "user_specified"),
                           user_specified_otu = NULL,
                           covariate.type = c("binary", "continuous"),
                           grp.ratio = 1, covariate.eff.mean = 1, covariate.eff.sd = 0,
                           confounder.type = c("none", "binary", "continuous", "both"),
                           conf.cov.cor = 0.6, conf.diff.otu.pct = 0, conf.nondiff.otu.pct = 0.1,
                           confounder.eff.mean = 0, confounder.eff.sd = 0, error.sd = 0,
                           depth.mu = 10000, depth.theta = 5, depth.conf.factor = 0, cont.conf, epsilon) {
  ## Estimated Dirichlet Parameters
  model.paras <- para

  ## Select number of OTUs
  ref.otu.tab <- model.paras$ref.otu.tab[(1:(nOTU)), ]

  ## Select number of Samples
  idx.otu <- rownames(ref.otu.tab)
  idx.sample <- colnames(model.paras$ref.otu.tab)[1:nSam]
  ref.otu.tab <- ref.otu.tab[, idx.sample]

  ## Confounding variable
  if (confounder.type == "none") {
    confounder.type <- "continuous"
    confounder.eff.mean <- 0
    confounder.eff.sd <- 0
    Z <- NULL
  }
  if (confounder.type == "continuous") {
    Z <- cbind(cont.conf)
  }
  if (confounder.type == "binary") {
    Z <- cbind(c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)))
  }
  if (confounder.type == "both") {
    Z <- cbind(rnorm(nSam), c(rep(0, nSam %/% 2), rep(1, nSam - nSam %/% 2)))
  }

  ## Covariate of Interest
  rho <- sqrt(conf.cov.cor^2/(1 - conf.cov.cor^2))
  if (covariate.type == "continuous") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + epsilon
  }
  if (covariate.type == "binary") {
    X <- rho * scale(scale(Z) %*% rep(1, ncol(Z))) + epsilon
    X <- cbind(ifelse(X <= quantile(X, grp.ratio/(1 + grp.ratio)), 0, 1))
  }
  rownames(X) <- colnames(ref.otu.tab)
  covariate.eff.mean1 <- covariate.eff.mean
  covariate.eff.mean2 <- covariate.eff.mean

  ## Simulate OTU Absolute Abundance
  ## Balanced OTU Counts
  if (diff.otu.direct == "balanced") {
    if (diff.otu.mode == "user_specified") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*%
        t(scale(X))
    }
    if (diff.otu.mode == "abundant") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*%
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- sample(c(rnorm(floor(nOTU/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
                           rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd))) %*%
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/4), mean = -covariate.eff.mean1, sd = covariate.eff.sd),
                             rnorm(floor(nOTU/2) - floor(nOTU/4), mean = covariate.eff.mean1, sd = covariate.eff.sd))),
                    sample(c(rnorm(floor((nOTU - floor(nOTU/2))/2), mean = -covariate.eff.mean2, sd = covariate.eff.sd),
                             rnorm(nOTU - floor(nOTU/2) - floor((nOTU - floor(nOTU/2))/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*%
        t(scale(X))
    }
  }

  ## Unbalanced OTU Counts
  if (diff.otu.direct == "unbalanced") {
    if (diff.otu.mode == "user_specified") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*%
        t(scale(X))
    }
    if (diff.otu.mode == "abundant") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*%
        t(scale(X))
    }
    else if (diff.otu.mode == "rare") {
      eta.diff <- rnorm(nOTU, mean = covariate.eff.mean2, sd = covariate.eff.sd) %*%
        t(scale(X))
    }
    else {
      eta.diff <- c(sample(c(rnorm(floor(nOTU/2), mean = covariate.eff.mean1, sd = covariate.eff.sd))),
                    sample(c(rnorm(nOTU - floor(nOTU/2), mean = covariate.eff.mean2, sd = covariate.eff.sd)))) %*%
        t(scale(X))
    }
  }
  eta.conf <- sample(c(rnorm(floor(nOTU/2), mean = -confounder.eff.mean, sd = confounder.eff.sd),
                       rnorm(nOTU - floor(nOTU/2), mean = confounder.eff.mean, sd = confounder.eff.sd))) %*%
    t(scale(scale(Z) %*% rep(1, ncol(Z))))
  otu.ord <- 1:(nOTU)
  diff.otu.ind <- NULL
  diff.otu.num <- round(diff.otu.pct * nOTU)
  if (diff.otu.mode == "user_specified")
    diff.otu.ind <- c(diff.otu.ind, which(idx.otu %in% user_specified_otu))
  if (diff.otu.mode == "mix")
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord, diff.otu.num))
  if (diff.otu.mode == "abundant")
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[1:round(length(otu.ord)/4)], diff.otu.num))
  if (diff.otu.mode == "rare")
    diff.otu.ind <- c(diff.otu.ind, sample(otu.ord[round(3 * length(otu.ord)/4):length(otu.ord)], diff.otu.num))
  if (length(diff.otu.ind) >= round(nOTU * conf.diff.otu.pct)) {
    conf.otu.ind1 <- sample(diff.otu.ind, round(nOTU * conf.diff.otu.pct))
  } else {
    conf.otu.ind1 <- diff.otu.ind
  }
  conf.otu.ind <- c(conf.otu.ind1,
                    sample(setdiff(1:(nOTU), diff.otu.ind), round(conf.nondiff.otu.pct * nOTU)))

  ## Calculate the new composition
  eta.diff[setdiff(1:(nOTU), diff.otu.ind), ] <- 0
  eta.conf[setdiff(1:(nOTU), conf.otu.ind), ] <- 0
  eta.error <- matrix(rnorm(nOTU * nSam, 0, error.sd), nOTU, nSam)
  eta.exp <- exp(t(eta.diff + eta.conf + eta.error))
  eta.exp <- eta.exp * t(ref.otu.tab)
  ref.otu.tab.prop <- eta.exp/rowSums(eta.exp) ## Normalized abundance
  ref.otu.tab.prop <- t(ref.otu.tab.prop)

  ## Simulate Sequencing Depth using Negative Binomial Distribution
  nSeq <- rnegbin(nSam,
                  mu = depth.mu * exp(scale(X) * depth.conf.factor),
                  theta = depth.theta)

  ## Simulate Absolute Abundance of OTU using Multinomial Sampling
  otu.tab.sim <- sapply(1:ncol(ref.otu.tab.prop),
                        function(i) rmultinom(1, nSeq[i], ref.otu.tab.prop[, i]))
  colnames(otu.tab.sim) <- rownames(eta.exp)
  rownames(otu.tab.sim) <- rownames(ref.otu.tab)
  diff.otu.ind <- (1:nOTU) %in% diff.otu.ind
  conf.otu.ind <- (1:nOTU) %in% conf.otu.ind

  ## Output results as a list
  return(list(otu.tab.sim = otu.tab.sim, covariate = X, confounder = Z,
              diff.otu.ind = diff.otu.ind, otu.names = idx.otu, conf.otu.ind = conf.otu.ind))
}

