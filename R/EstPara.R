#' Sample from the Dirichlet Distribution (Matrix Version)
#'
#' Generates samples from the Dirichlet distribution using gamma distributions.
#' This is a helper function for simulating OTU proportions in a matrix format.
#'
#' @param alpha A matrix of Dirichlet parameters (same dimension as the OTU table).
#' Each column corresponds to one sample's Dirichlet parameters for all OTUs.
#'
#' @return A matrix of proportions with the same dimensions as \code{alpha}, where each column sums to 1.
#'
#' @keywords internal
rdirichlet.m <- function (alpha) {
  Gam <- matrix(rgamma(length(alpha), shape = alpha), nrow(alpha), ncol(alpha))
  t(t(Gam) / colSums(Gam))
}

#' Estimate Dirichlet Parameters and Simulate OTU Table
#'
#' This function estimates the Dirichlet hyperparameters using the \code{dirmult} package
#' and generates a simulated OTU count matrix and proportion matrix. It adds a pseudocount
#' based on the estimated Dirichlet parameters and applies sample-specific size factors.
#'
#' @param ref.otu.tab A matrix of raw OTU counts (rows = OTUs, columns = samples).
#' Row names should be OTU names and column names should be sample names.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{mu}}{Estimated OTU proportions matrix (OTUs × Samples) after Dirichlet smoothing.}
#'   \item{\code{ref.otu.tab}}{Simulated OTU count matrix (OTUs × Samples) after applying size factors.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Estimates Dirichlet parameters using the \code{dirmult} package.
#'   \item Adds a pseudocount to each OTU count based on the estimated Dirichlet parameters.
#'   \item Converts the adjusted counts into OTU proportions using the Dirichlet distribution.
#'   \item Applies randomly sampled size factors (from a log-normal distribution) to simulate sequencing depth variability.
#' }
#'
#' OTUs are ordered by mean abundance across all samples.
#'
#' @seealso \code{\link[dirmult]{dirmult}}
#'
#' @importFrom dirmult dirmult
#' @export
EstPara <- function (ref.otu.tab) {

  if (is.null(rownames(ref.otu.tab))) {
    rownames(ref.otu.tab) <- paste0('OTU', 1 : nrow(ref.otu.tab))
  } # otu * sample
  samplenames = colnames(ref.otu.tab)
  taxnames = rownames(ref.otu.tab)

  dirmult.paras <- dirmult::dirmult(t(ref.otu.tab))

  gamma = dirmult.paras$gamma
  names(gamma) = names(dirmult.paras$pi)

  # Add pseduo count(each OTU add gamma estimated from dirmult)
  ref.otu.tab = sapply(1:ncol(ref.otu.tab), function (i) gamma + ref.otu.tab[,i]) # C_ij otu * sample

  # back to Dirichlet, calculate the true proportion
  ref.otu.tab.p <- rdirichlet.m(ref.otu.tab) # P_ij nOTU*nSam
  colnames(ref.otu.tab.p) = samplenames
  rownames(ref.otu.tab.p) = taxnames

  # order OTUs by mean OTU proportion, for later selection
  ord = order(rowMeans(ref.otu.tab.p), decreasing = TRUE)
  ref.otu.tab.p =  ref.otu.tab.p[ord,]

  # apply size factor
  Si = exp(rnorm(ncol(ref.otu.tab.p)))
  ref.otu.tab0 = t(t(ref.otu.tab.p)*Si)
  colnames(ref.otu.tab0) = colnames(ref.otu.tab.p)
  return(list(mu = ref.otu.tab.p, ref.otu.tab = ref.otu.tab0))
}
