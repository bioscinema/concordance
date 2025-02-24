#' Calculate Overall Concordance Metrics and Plot Their Distributions
#'
#' This function computes overall concordance metrics (IPAC and IWPAC) for a merged phyloseq object
#' by calling \code{ALL_Concordance}. The resulting metrics are then pivoted to a long format, and group-level
#' summary statistics (mean and median) are computed for each concordance measure. Finally, a faceted ggplot2
#' histogram and density plot is generated to display the distribution of concordance values with vertical
#' lines indicating the mean and median.
#'
#' @param physeq A merged phyloseq object for which overall concordance metrics are to be computed.
#'   The object should contain sample data and an OTU table with subject names, where each subject is
#'   represented by both WGS and 16S data.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{results}{A data frame in long format containing the overall concordance metrics (IPAC and IWPAC)
#'       for each subject, along with the mean and median concordance values for each metric.}
#'     \item{plot}{A ggplot2 object displaying histograms and density curves for the concordance metrics,
#'       with vertical dashed lines for the mean and median values.}
#'   }
#'
#' @details
#' The function first calls \code{ALL_Concordance} to compute the overall inverse concordance metrics
#' for each subject from the merged phyloseq object. The output is then reshaped to long format using
#' \code{tidyr::pivot_longer}, so that the two metrics (IPAC and IWPAC) are contained in a single column.
#' For each concordance measure, the mean and median are computed and added to the data frame. Finally, a
#' ggplot2 object is created with a histogram (binwidth = 0.02) and density curve, and vertical dashed lines
#' are added to indicate the computed mean and median. The plot is faceted by the concordance measure.
#'
#' @export
concordance_overall <- function(physeq) {
  # Calculate overall concordance metrics using ALL_Concordance
  results_combined <- ALL_Concordance(physeq)

  # Pivot the results to long format and compute summary statistics (mean and median) for each measure.
  results_combined <- results_combined %>%
    tidyr::pivot_longer(cols = c("IPAC", "IWPAC"),
                        names_to = "Concordance_Measure",
                        values_to = "Concordance") %>%
    dplyr::group_by(Concordance_Measure) %>%
    dplyr::mutate(mean_concordance = mean(Concordance, na.rm = TRUE),
                  median_concordance = median(Concordance, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # Create a ggplot2 plot of the concordance distributions
  Concordance_plots <- ggplot2::ggplot(results_combined, ggplot2::aes(x = Concordance)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..), binwidth = 0.02, fill = "lightblue", color = "black") +
    ggplot2::geom_density(color = "black", size = 0.5) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean_concordance, color = "Mean", linetype = "Mean"), size = 0.5) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = median_concordance, color = "Median", linetype = "Median"), size = 0.5) +
    ggplot2::facet_grid(Concordance_Measure ~ .) +
    ggplot2::labs(title = NULL, x = NULL, y = "Density") +
    ggplot2::scale_color_manual(values = c("Mean" = "#E69F00", "Median" = "#CC79A7"), name = "Statistics") +
    ggplot2::scale_linetype_manual(values = c("Mean" = "dashed", "Median" = "dashed"), name = "Statistics") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = ggplot2::element_rect(color = "black", fill = "gray90")
    )

  return(list(results = results_combined, plot = Concordance_plots))
}
