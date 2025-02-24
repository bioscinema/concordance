#' Calculate Required Sample Size and Total Cost for 16S and WGS to Achieve Target Power
#'
#' This function runs power calculations using the simulation function \code{sampleCalculationRange} for both
#' 16S (amplicon) and whole-genome shotgun (WGS) sequencing data. It estimates the minimum sample size required
#' to achieve a specified target power and calculates the total cost based on per-sample cost estimates for each method.
#'
#' @param target_power Numeric. The desired statistical power (e.g., 0.8). Default is 0.8.
#' @param lower_bound Numeric. The minimum sample size to consider. Default is 50.
#' @param upper_bound Numeric. The maximum sample size to consider. Default is 500.
#' @param step Numeric. The increment between sample sizes. Default is 5.
#' @param cost16s Numeric. The estimated cost per sample for 16S sequencing. Default is 100.
#' @param costWGS Numeric. The estimated cost per sample for WGS sequencing. Default is 500.
#' @param data.type Character. The type of data to simulate (e.g., "gut", "oral", "infant"). Default is "gut".
#' @param nSim Numeric. Number of simulations to run per sample size. Default is 100.
#' @param alpha Numeric. Significance level for power calculation. Default is 0.05.
#' @param nOTU Numeric. Number of OTUs to simulate. Default is 500.
#' @param diff.otu.pct Numeric. Proportion of OTUs that are differentially abundant. Default is 0.1.
#' @param diff.otu.direct Character. Direction of differential abundance ("balanced" by default).
#' @param diff.otu.mode Character. Mode of differential abundance ("abundant" by default).
#' @param user_specified_otu Optional. User-specified OTU data; default is NULL.
#' @param covariate.type Character. Type of covariate ("binary" by default).
#' @param grp.ratio Numeric. Group ratio for the covariate. Default is 1.
#' @param covariate.eff.mean Numeric. Mean effect of the covariate. Default is 1.
#' @param covariate.eff.sd Numeric. Standard deviation of the covariate effect. Default is 0.
#' @param confounder.type Character. Type of confounder ("none" by default).
#' @param conf.cov.cor Numeric. Correlation between confounders and covariates. Default is 0.6.
#' @param conf.diff.otu.pct Numeric. Proportion of OTUs differentially affected by confounders. Default is 0.
#' @param conf.nondiff.otu.pct Numeric. Proportion of non-differential OTUs affected by confounders. Default is 0.1.
#' @param confounder.eff.mean Numeric. Mean effect of the confounder. Default is 0.
#' @param confounder.eff.sd Numeric. Standard deviation of the confounder effect. Default is 0.
#' @param error.sd Numeric. Standard deviation of error. Default is 0.
#' @param depth.mu Numeric. Mean sequencing depth. Default is 10000.
#' @param depth.theta Numeric. Dispersion parameter for sequencing depth. Default is 5.
#' @param depth.conf.factor Numeric. Depth configuration factor. Default is 0.
#' @param cont.conf Numeric vector. Continuous confounder values; if NULL, generated randomly. Default is NULL.
#' @param epsilon Numeric vector. Error vector; if NULL, generated randomly. Default is NULL.
#'
#' @return A list containing:
#' \describe{
#'   \item{required_sample_size_16s}{Minimum sample size required to achieve the target power for 16S data.}
#'   \item{total_cost_16s}{Total estimated cost for 16S data at the required sample size.}
#'   \item{required_sample_size_wgs}{Minimum sample size required to achieve the target power for WGS data.}
#'   \item{total_cost_wgs}{Total estimated cost for WGS data at the required sample size.}
#' }
#'
#' @examples
#' \dontrun{
#' price_estimates <- priceCalculation(target_power = 0.8,
#'                                     lower_bound = 50,
#'                                     upper_bound = 500,
#'                                     step = 5,
#'                                     cost16s = 100,
#'                                     costWGS = 500,
#'                                     data.type = "gut",
#'                                     nSim = 100,
#'                                     alpha = 0.05,
#'                                     nOTU = 500,
#'                                     diff.otu.pct = 0.1,
#'                                     diff.otu.direct = "balanced",
#'                                     diff.otu.mode = "abundant",
#'                                     user_specified_otu = NULL,
#'                                     covariate.type = "binary",
#'                                     grp.ratio = 1,
#'                                     covariate.eff.mean = 1,
#'                                     covariate.eff.sd = 0,
#'                                     confounder.type = "none",
#'                                     conf.cov.cor = 0.6,
#'                                     conf.diff.otu.pct = 0,
#'                                     conf.nondiff.otu.pct = 0.1,
#'                                     confounder.eff.mean = 0,
#'                                     confounder.eff.sd = 0,
#'                                     error.sd = 0,
#'                                     depth.mu = 10000,
#'                                     depth.theta = 5,
#'                                     depth.conf.factor = 0,
#'                                     cont.conf = NULL,
#'                                     epsilon = NULL)
#' print(price_estimates)
#' }
#'
#' @export
priceCalculation <- function(target_power = 0.8,
                             lower_bound = 50,
                             upper_bound = 500,
                             step = 5,
                             cost16s = 100,
                             costWGS = 500,
                             data.type = "gut",
                             nSim = 100,
                             alpha = 0.05,
                             nOTU = 500,
                             diff.otu.pct = 0.1,
                             diff.otu.direct = "balanced",
                             diff.otu.mode = "abundant",
                             user_specified_otu = NULL,
                             covariate.type = "binary",
                             grp.ratio = 1,
                             covariate.eff.mean = 1,
                             covariate.eff.sd = 0,
                             confounder.type = "none",
                             conf.cov.cor = 0.6,
                             conf.diff.otu.pct = 0,
                             conf.nondiff.otu.pct = 0.1,
                             confounder.eff.mean = 0,
                             confounder.eff.sd = 0,
                             error.sd = 0,
                             depth.mu = 10000,
                             depth.theta = 5,
                             depth.conf.factor = 0,
                             cont.conf = NULL,
                             epsilon = NULL) {

  # If continuous confounders or error vectors are not provided, generate random ones with length equal to lower_bound.
  # if (is.null(cont.conf)) {
  #   cont.conf <- rnorm(lower_bound)
  # }
  # if (is.null(epsilon)) {
  #   epsilon <- rnorm(lower_bound)
  # }

  # Run simulation for 16S (amplicon) method
  result16s <- sampleCalculationRange(lower_bound = lower_bound,
                                      upper_bound = upper_bound,
                                      step = step,
                                      data.type = data.type,
                                      method = "amp",
                                      nSim = nSim,
                                      alpha = alpha,
                                      nOTU = nOTU,
                                      diff.otu.pct = diff.otu.pct,
                                      diff.otu.direct = diff.otu.direct,
                                      diff.otu.mode = diff.otu.mode,
                                      user_specified_otu = user_specified_otu,
                                      covariate.type = covariate.type,
                                      grp.ratio = grp.ratio,
                                      covariate.eff.mean = covariate.eff.mean,
                                      covariate.eff.sd = covariate.eff.sd,
                                      confounder.type = confounder.type,
                                      conf.cov.cor = conf.cov.cor,
                                      conf.diff.otu.pct = conf.diff.otu.pct,
                                      conf.nondiff.otu.pct = conf.nondiff.otu.pct,
                                      confounder.eff.mean = confounder.eff.mean,
                                      confounder.eff.sd = confounder.eff.sd,
                                      error.sd = error.sd,
                                      depth.mu = depth.mu,
                                      depth.theta = depth.theta,
                                      depth.conf.factor = depth.conf.factor,
                                      cont.conf = NULL,
                                      epsilon = NULL)

  df16s <- result16s$data

  # Determine required sample size for 16S to achieve target power
  required_size_16s <- min(df16s$sample_size[df16s$mean_power >= target_power], na.rm = TRUE)
  total_cost_16s <- required_size_16s * cost16s

  # Run simulation for WGS method
  resultWGS <- sampleCalculationRange(lower_bound = lower_bound,
                                      upper_bound = upper_bound,
                                      step = step,
                                      data.type = data.type,
                                      method = "wgs",
                                      nSim = nSim,
                                      alpha = alpha,
                                      nOTU = nOTU,
                                      diff.otu.pct = diff.otu.pct,
                                      diff.otu.direct = diff.otu.direct,
                                      diff.otu.mode = diff.otu.mode,
                                      user_specified_otu = user_specified_otu,
                                      covariate.type = covariate.type,
                                      grp.ratio = grp.ratio,
                                      covariate.eff.mean = covariate.eff.mean,
                                      covariate.eff.sd = covariate.eff.sd,
                                      confounder.type = confounder.type,
                                      conf.cov.cor = conf.cov.cor,
                                      conf.diff.otu.pct = conf.diff.otu.pct,
                                      conf.nondiff.otu.pct = conf.nondiff.otu.pct,
                                      confounder.eff.mean = confounder.eff.mean,
                                      confounder.eff.sd = confounder.eff.sd,
                                      error.sd = error.sd,
                                      depth.mu = depth.mu,
                                      depth.theta = depth.theta,
                                      depth.conf.factor = depth.conf.factor,
                                      cont.conf = NULL,
                                      epsilon = NULL)

  dfWGS <- resultWGS$data

  required_size_wgs <- min(dfWGS$sample_size[dfWGS$mean_power >= target_power], na.rm = TRUE)
  total_cost_wgs <- required_size_wgs * costWGS

  # Create a data frame to return the results
  result_df <- data.frame(
    platform = c("16S", "WGS"),
    sample_size_needed = c(required_size_16s, required_size_wgs),
    price = c(total_cost_16s, total_cost_wgs)
  )

  return(result_df)
}
