#' Sample Calculation and Plot of Power Versus Sample Size (Range)
#'
#' This function calculates simulation-based power for a range of sample sizes. The user supplies a
#' lower bound and an upper bound, and the function computes the mean power at every 5-sample increment
#' by calling the existing \code{powerCalculation} function. A plot is generated showing the relationship
#' between sample size and power, and a data frame with the detailed results is returned.
#'
#' @param lower_bound An integer specifying the minimum sample size to evaluate.
#' @param upper_bound An integer specifying the maximum sample size to evaluate.
#' @param step An integer specifying the increment in sample size. Default is 5.
#' @param data.type A character string specifying the data type. Must be one of \code{"gut"}, \code{"oral"}, or \code{"infant"}.
#' @param method A character string specifying the sequencing/estimation method. Must be one of \code{"amp"} or \code{"wgs"}.
#' @param nSim Number of simulation replicates for each sample size. Default is 100.
#' @param alpha Significance level for the PERMANOVA test. Default is 0.05.
#' @param nOTU Number of OTUs to simulate. Default is 500.
#' @param diff.otu.pct Proportion of OTUs that are differentially abundant. Default is 0.1.
#' @param diff.otu.direct Direction of differential OTUs. Either \code{"balanced"} or \code{"unbalanced"}. Default is \code{"balanced"}.
#' @param diff.otu.mode Mode of differential OTUs. Options include \code{"abundant"}, \code{"rare"}, \code{"mix"}, or \code{"user_specified"}. Default is \code{"abundant"}.
#' @param user_specified_otu Optional vector specifying the OTUs that are differentially abundant when \code{diff.otu.mode} is \code{"user_specified"}. Default is \code{NULL}.
#' @param covariate.type Type of covariate; either \code{"binary"} or \code{"continuous"}. Default is \code{"binary"}.
#' @param grp.ratio Group ratio for binary covariates. Default is 1.
#' @param covariate.eff.mean Mean effect size for the covariate. Default is 1.
#' @param covariate.eff.sd Standard deviation of the covariate effect. Default is 0.
#' @param confounder.type Type of confounding variable. Options are \code{"none"}, \code{"binary"}, \code{"continuous"}, or \code{"both"}. Default is \code{"none"}.
#' @param conf.cov.cor Correlation between the covariate and the confounder. Default is 0.6.
#' @param conf.diff.otu.pct Proportion of differential OTUs affected by the confounder. Default is 0.
#' @param conf.nondiff.otu.pct Proportion of non-differential OTUs affected by the confounder. Default is 0.1.
#' @param confounder.eff.mean Mean effect size for the confounder. Default is 0.
#' @param confounder.eff.sd Standard deviation of the confounder effect. Default is 0.
#' @param error.sd Standard deviation for the error term. Default is 0.
#' @param depth.mu Mean sequencing depth. Default is 10000.
#' @param depth.theta Theta parameter for the negative binomial distribution of sequencing depth. Default is 5.
#' @param depth.conf.factor Factor for the covariate effect on sequencing depth. Default is 0.
#' @param cont.conf A vector of continuous confounder values. If not provided, new values will be generated for each simulation using \code{rnorm(n)}.
#' @param epsilon A vector of random error values. If not provided, new values will be generated for each simulation using \code{rnorm(n)}.
#'
#' @return A data frame with columns \code{sample_size} and \code{mean_power} showing the estimated power at each sample size.
#'
#' @examples
#' \dontrun{
#'   # Evaluate power for sample sizes from 50 to 200 (in increments of 5) for gut data using amp parameters,
#'   # with 50 simulation replicates and 200 OTUs.
#'   res <- sampleCalculationRange(lower_bound = 50, upper_bound = 200,
#'                                 data.type = "gut", method = "amp", nSim = 50, nOTU = 200)
#'   print(res)
#' }
#'
#' @import ggplot2
#' @export
sampleCalculationRange <- function(lower_bound, upper_bound, step = 5,
                                   data.type = c("gut", "oral", "infant"),
                                   method = c("amp", "wgs"),
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
  data.type <- match.arg(data.type)
  method <- match.arg(method)

  # Generate a sequence of sample sizes from lower_bound to upper_bound.
  sample_sizes <- seq(lower_bound, upper_bound, by = step)

  # For each sample size, call powerCalculation and extract the overall power.
  mean_powers <- sapply(sample_sizes, function(n) {
    # print(n)
    # Generate new confounder and error vectors if not provided.
    if (is.null(cont.conf)) {
      cont_conf_vec <- rnorm(n)
    } else {
      cont_conf_vec <- cont.conf
    }
    if (is.null(epsilon)) {
      epsilon_vec <- rnorm(n)
    } else {
      epsilon_vec <- epsilon
    }

    result <- powerCalculation(data.type = data.type,
                               method = method,
                               nSim = nSim,
                               alpha = alpha,
                               nSam = n,  # current sample size
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
                               cont.conf = cont_conf_vec,
                               epsilon = epsilon_vec)
    result$overall_power
  })

  # Create a data frame with the sample sizes and their corresponding mean power.
  df <- data.frame(sample_size = sample_sizes,
                   mean_power = mean_powers)

  # Use ggplot2 to create the plot.
  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample_size, y = mean_power)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Sample Size", y = "Mean Power",
                  title = paste("Power vs. Sample Size for", data.type, method, "data")) +
    ggplot2::theme_minimal()

  return(list(data = df, plot = p))
}
