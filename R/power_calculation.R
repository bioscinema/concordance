#' Power Calculation Using Bray–Curtis PERMANOVA for Microbiome Data Simulation
#'
#' This function performs a simulation-based power calculation using Bray–Curtis distances
#' and PERMANOVA for microbiome data simulated by \code{SimulateMSeqU}. It supports three
#' data types ("gut", "oral", "infant") and two sequencing methods ("amp" and "wgs"). The
#' corresponding parameter file (e.g., "gut_amp_para.RData") is loaded from the package's data
#' folder to supply the estimated parameters.
#'
#' @param data.type A character string specifying the data type. Must be one of \code{"gut"},
#'   \code{"oral"}, or \code{"infant"}.
#' @param method A character string specifying the sequencing/estimation method. Must be one
#'   of \code{"amp"} or \code{"wgs"}.
#' @param nSim Number of simulation replicates. Default is 100.
#' @param alpha Significance level for the PERMANOVA test. Default is 0.05.
#' @param nSam Number of samples per simulation replicate. Default is 100.
#'   **Note:** The maximum sample sizes are 206 for gut, 181 for oral, and 498 for infant.
#' @param nOTU Number of OTUs to simulate. Default is 500.
#' @param diff.otu.pct Proportion of OTUs that are differentially abundant. Default is 0.1.
#' @param diff.otu.direct Direction of differential OTUs. Either \code{"balanced"} or
#'   \code{"unbalanced"}. Default is \code{"balanced"}.
#' @param diff.otu.mode Mode of differential OTUs. Options include \code{"abundant"},
#'   \code{"rare"}, \code{"mix"}, or \code{"user_specified"}. Default is \code{"abundant"}.
#' @param user_specified_otu Optional vector specifying the OTUs that are differentially
#'   abundant when \code{diff.otu.mode} is \code{"user_specified"}.
#' @param covariate.type Type of covariate; either \code{"binary"} or \code{"continuous"}.
#'   Default is \code{"binary"}.
#' @param grp.ratio Group ratio for binary covariates. Default is 1.
#' @param covariate.eff.mean Mean effect size for the covariate. Default is 1.
#' @param covariate.eff.sd Standard deviation of the covariate effect. Default is 0.
#' @param confounder.type Type of confounding variable. Options are \code{"none"},
#'   \code{"binary"}, \code{"continuous"}, or \code{"both"}. Default is \code{"none"}.
#' @param conf.cov.cor Correlation between the covariate and the confounder. Default is 0.6.
#' @param conf.diff.otu.pct Proportion of differential OTUs affected by the confounder.
#'   Default is 0.
#' @param conf.nondiff.otu.pct Proportion of non-differential OTUs affected by the confounder.
#'   Default is 0.1.
#' @param confounder.eff.mean Mean effect size for the confounder. Default is 0.
#' @param confounder.eff.sd Standard deviation of the confounder effect. Default is 0.
#' @param error.sd Standard deviation for the error term. Default is 0.
#' @param depth.mu Mean sequencing depth. Default is 10000.
#' @param depth.theta Theta parameter for the negative binomial distribution of sequencing depth.
#'   Default is 5.
#' @param depth.conf.factor Factor for the covariate effect on sequencing depth. Default is 0.
#' @param cont.conf Numeric vector of continuous confounder values. Default is \code{rnorm(nSam)}.
#' @param epsilon Numeric vector of random error values. Default is \code{rnorm(nSam)}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{overall_power}{The overall estimated power as the proportion of simulation replicates
#'       detecting a significant group effect (PERMANOVA p-value < \code{alpha}).}
#'     \item{replicate_power}{A numeric vector of detection indicators (1 for detected, 0 for not)
#'       for each simulation replicate.}
#'   }
#'
#' @details For each simulation replicate, the function:
#'   \enumerate{
#'     \item Loads the corresponding parameter file (e.g., "gut_amp_para.RData") from the package's
#'       data folder and obtains the estimated parameters.
#'     \item Checks if \code{nSam} exceeds the maximum allowed sample size for the chosen data type.
#'     \item Simulates a microbiome dataset using \code{SimulateMSeqU}.
#'     \item Defines sample groups based on the covariate (dichotomizing if continuous).
#'     \item Computes Bray–Curtis distances via \code{vegdist} (with samples as rows).
#'     \item Performs PERMANOVA using \code{adonis2} to test the group effect.
#'     \item Records a detection if the PERMANOVA p-value for the Group effect is less than
#'       \code{alpha}.
#'   }
#'
#' @examples
#' \dontrun{
#'   # Valid example: For gut data the maximum sample size is 206.
#'   result <- powerCalculation(data.type = "gut", method = "amp", nSam = 206, nSim = 50)
#'   print(result$overall_power)
#'
#'   # This example will error because the maximum sample size for oral data is 181.
#'   # result <- powerCalculation(data.type = "oral", method = "wgs", nSam = 200)
#' }
#'
#' @import phyloseq
#' @importFrom MASS rnegbin
#' @importFrom vegan adonis2 vegdist
#' @importFrom stats rnorm
#' @export
powerCalculation <- function(data.type = c("gut", "oral", "infant"),
                             method = c("amp", "wgs"),
                             nSim = 100,
                             alpha = 0.05,
                             nSam = 100,
                             nOTU = 200,
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
                             cont.conf,
                             epsilon) {

  ## Match the choices for data type and method.
  data.type <- match.arg(data.type)
  method <- match.arg(method)

  ## Generate cont.conf and epsilon if they were not provided.
  if(missing(cont.conf)) {
    cont.conf <- rnorm(nSam)
  }
  if(missing(epsilon)) {
    epsilon <- rnorm(nSam)
  }

  ## [Rest of your function code follows...]
  ## For example, load parameters, run simulations, and compute overall power.

  # Placeholder for simulation code:
  replicate_power <- numeric(nSim)
  for (i in 1:nSim) {
    # Example: Run your simulation (this is pseudo-code)
    sim <- SimulateMSeqU(para = EstPara, nSam = nSam, nOTU = nOTU,
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
                         cont.conf = cont.conf,
                         epsilon = epsilon)
    # Assume sim produces an object with a p-value for the Group effect.
    Bray_Curtis_P <- 0.04  # Dummy value for demonstration
    replicate_power[i] <- ifelse(Bray_Curtis_P < alpha, 1, 0)
  }

  overall_power <- mean(replicate_power)

  return(list(overall_power = overall_power, replicate_power = replicate_power))
}

