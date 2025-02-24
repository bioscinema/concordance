#' Calculate and Plot Taxa Concordance for Merged Phyloseq Data
#'
#' This function merges two phyloseq objects (one from 16S data and one from WGS data) using
#' the \code{merge_phy} function, and then calculates taxa concordance metrics (PAC and WPAC)
#' for a set of specified taxonomic levels. The results are combined into a long-format data frame,
#' and a ggplot2 object is generated displaying histograms and density curves for each concordance
#' measure, with vertical lines indicating the mean and median.
#'
#' @param phyloseq1 A phyloseq object containing 16S sequencing data.
#' @param phyloseq2 A phyloseq object containing WGS sequencing data.
#' @param levels A character vector specifying the taxonomic levels at which to calculate concordance.
#'   Default is \code{c("Phylum", "Class", "Order", "Family", "Genus")}.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{results_combined}{A data frame in long format containing the calculated concordance
#'       metrics (PAC and WPAC) for each taxonomic level.}
#'     \item{Concordance_plots}{A ggplot2 object showing the distributions of concordance values with
#'       histograms, density curves, and vertical lines for the mean and median.}
#'   }
#'
#' @details
#' The function first calls \code{merge_phy} (which must be defined elsewhere in your package)
#' to merge the two input phyloseq objects. Then, for each taxonomic level in \code{levels}, it computes
#' the concordance using the \code{Taxa_Concordance} function. The results are combined, pivoted to a long format,
#' and summarized by taxonomic level and concordance measure. Finally, a faceted ggplot is generated to visualize
#' the distributions of concordance metrics.
#'
#' @importFrom stats median quantile rmultinom
#' @export
concordance_split <- function(phyloseq1, phyloseq2, levels = c("Phylum", "Class", "Order", "Family", "Genus")) {

  # Merge the two phyloseq objects using merge_phy (which merges at Genus level by default)
  merged_phy <- merge_phy(phyloseq1, phyloseq2, level = "Genus")

  # Calculate concordance for each level using Taxa_Concordance
  results <- lapply(levels, function(lvl) {
    Taxa_Concordance(merged_phy, lvl)
  })
  names(results) <- levels

  # Combine the results into a single data frame, adding a column for taxonomic level
  results_combined <- do.call(rbind, lapply(names(results), function(lvl) {
    df <- results[[lvl]]
    df$Level <- lvl
    return(df)
  }))

  # Pivot the data to long format for the two concordance measures
  results_combined <- results_combined %>%
    dplyr::mutate(Level = factor(Level, levels = levels)) %>%
    tidyr::pivot_longer(cols = c("PAC", "WPAC"),
                        names_to = "Concordance_Measure",
                        values_to = "Concordance") %>%
    dplyr::group_by(Level, Concordance_Measure) %>%
    dplyr::mutate(mean_concordance = mean(Concordance, na.rm = TRUE),
                  median_concordance = median(Concordance, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # Create the plot using ggplot2
  Concordance_plots <- ggplot2::ggplot(results_combined, ggplot2::aes(x = Concordance)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ..density..), binwidth = 0.05, fill = "lightblue", color = "black") +
    ggplot2::geom_density(color = "black", size = 0.5) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean_concordance, color = "Mean", linetype = "Mean"), size = 0.5) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = median_concordance, color = "Median", linetype = "Median"), size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::facet_grid(Level ~ Concordance_Measure) +
    ggplot2::labs(title = NULL, x = NULL, y = "Density") +
    ggplot2::scale_color_manual(values = c("Mean" = "#E69F00", "Median" = "#CC79A7"), name = "Statistics") +
    ggplot2::scale_linetype_manual(values = c("Mean" = "dashed", "Median" = "dashed"), name = "Statistics") +
    ggplot2::theme(
      legend.position = "bottom",
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      strip.background = ggplot2::element_rect(color = "black", fill = "gray90")
    )

  return(list(results_combined = results_combined, Concordance_plots = Concordance_plots))
}
