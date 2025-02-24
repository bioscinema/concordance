#' Calculate Taxa Concordance between WGS and 16S Sequencing Data
#'
#' This function computes the concordance between taxa detection in whole genome shotgun sequencing
#' (WGS) and 16S rRNA gene sequencing. It first aggregates taxonomic counts to a specified level (e.g., "Phylum")
#' and transforms the counts to relative abundances. Then, for each subject, the function calculates two concordance
#' metrics: the Presence/Absence Concordance (PAC) and the Weighted PAC (WPAC).
#'
#' @param physeq A phyloseq object containing sequencing data, taxonomic assignments, and sample metadata.
#' @param level Character string specifying the taxonomic level to which the data should be aggregated.
#'   Default is "Phylum".
#'
#' @return A data frame with one row per subject containing:
#'   \describe{
#'     \item{Subject}{The subject identifier.}
#'     \item{PAC}{The Presence/Absence Concordance between WGS and 16S data, computed as
#'       \eqn{1 - \frac{\sum |(count\_wgs > 0) - (count\_16s > 0)|}{\sum ((count\_wgs > 0) + (count\_16s > 0))}}.}
#'     \item{WPAC}{The Weighted Presence/Absence Concordance, computed as
#'       \eqn{1 - \frac{\sum |count\_wgs - count\_16s|}{2}}.}
#'   }
#'
#' @details
#' The function aggregates taxa using \code{aggregate_taxa} and converts the OTU counts to relative abundances
#' via \code{transform_sample_counts}. It assumes that the OTU table contains columns for both WGS and 16S data,
#' with column names following the format \code{<subject>-wgs} and \code{<subject>-16s} respectively. For each subject,
#' the function calculates the PAC metric based on the presence or absence of taxa and the WPAC metric based on the
#' absolute differences in relative abundances.
#'
#' @export
Taxa_Concordance <- function(physeq, level = "Phylum"){
  ## Relative abundance table
  physeq <- aggregate_taxa(physeq, level = level)
  physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
  subjects <- physeq@sam_data$sample_names
  count_tab <- physeq@otu_table

  PAC <- list()
  WPAC <- list()

  for (subject in subjects) {
    wgs_id <- paste0(subject, "-wgs")
    # print(wgs_id)
    sixs_id <- paste0(subject, "-16s")
    count_wgs <- count_tab[, wgs_id]
    count_sixs <- count_tab[, sixs_id]

    # PAC: 1 - (abs diff in presence/absence) / (total number of presences)
    PAC[[subject]] <- 1 - sum(abs((count_wgs > 0) - (count_sixs > 0))) / sum((count_wgs > 0) + (count_sixs > 0))

    # WPAC: 1 - (absolute differences in counts) / 2
    WPAC[[subject]] <- 1 - sum(abs(count_wgs - count_sixs)) / 2
  }

  results_df <- data.frame(Subject = subjects,
                           PAC = unlist(PAC),
                           WPAC = unlist(WPAC))

  return(results_df)
}


##' @title fix_duplicate_tax
##'
##' @param physeq a phyloseq object

##' @export
##' @description fix the duplicatae taxonomy names of a phyloseq object

fix_duplicate_tax = function(physeq){
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package \"phyloseq\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  taxtab <- phyloseq::tax_table(physeq)
  for(i in 3:ncol(taxtab)){
    uniqs = unique(taxtab[,i])
    for(j in 1:length(uniqs)){
      if(is.na(uniqs[j])) next
      ind = which(taxtab[,i]== as.character(uniqs[j]))
      if(length(unique(taxtab[ind,i-1]))>1){
        taxtab[ind,i] = paste(taxtab[ind,i-1], taxtab[ind,i], sep="_")
      }
    }
  }
  phyloseq::tax_table(physeq) = taxtab
  return(physeq)
}


#' Merge 16S and WGS Phyloseq Objects with Taxonomy and Phylogeny Adjustment
#'
#' This function processes two phyloseq objects—one containing 16S sequencing data and
#' one containing WGS sequencing data. It aggregates the data to a specified taxonomic level,
#' updates sample names and metadata to indicate the sequencing platform, standardizes the taxonomic
#' tables, merges the two objects, fixes duplicate taxa using \code{BAZE::fix_duplicate_tax}, and
#' updates the phylogenetic tree tip labels by removing any leading taxonomic prefix.
#'
#' @param phylo_16s A phyloseq object containing 16S sequencing data.
#' @param phylo_wgs A phyloseq object containing WGS sequencing data.
#' @param level A character string specifying the taxonomic level to which the data should be aggregated.
#'   Default is \code{"Genus"}.
#'
#' @return A merged and processed phyloseq object with updated sample names, taxonomy, and phylogenetic tree.
#'
#' @details The function performs the following steps:
#' \enumerate{
#'   \item Aggregates taxa for both 16S and WGS objects to the specified level using \code{aggregate_taxa}.
#'   \item Modifies sample names by appending \code{"-16s"} or \code{"-wgs"} to indicate the data source.
#'   \item Updates sample metadata to include a \code{Platform} column.
#'   \item Standardizes the taxonomic tables by renaming columns to \code{"Kingdom"}, \code{"Phylum"},
#'         \code{"Class"}, \code{"Order"}, \code{"Family"}, \code{"Genus"}, and \code{"unique"}.
#'   \item Merges the two phyloseq objects and processes the merged taxonomic table by removing the
#'         \code{"unique"} column and converting taxonomic columns to factors.
#'   \item Fixes duplicate taxa using \code{BAZE::fix_duplicate_tax}.
#'   \item Adjusts the phylogenetic tree by removing any leading taxonomic prefix (matching \code{"^[a-z]__"})
#'         from the tip labels.
#' }
#' @importFrom dplyr filter select
#' @importFrom BAZE phy_to_tax
#' @export
merge_phy <- function(phylo_16s, phylo_wgs, level = "Genus") {

  # Aggregate taxa to the specified level for both datasets.
  physeq.16s <- aggregate_taxa(phylo_16s, level = level)
  physeq.wgs <- aggregate_taxa(phylo_wgs, level = level)

  # Update sample names to indicate platform.
  sample_16s <- sample_names(physeq.16s)
  sample_names(physeq.16s) <- paste(sample_16s, "-16s", sep = "")

  sample_wgs <- sample_names(physeq.wgs)
  sample_names(physeq.wgs) <- paste(sample_wgs, "-wgs", sep = "")

  # Update sample metadata to include a Platform column.
  mysample_16s <- data.frame(physeq.16s@sam_data, check.names = FALSE)
  mysample_16s$Platform <- rep("16s", length(sample_16s))
  sample_data(physeq.16s) <- sample_data(mysample_16s)

  mysample_wgs <- data.frame(physeq.wgs@sam_data, check.names = FALSE)
  mysample_wgs$Platform <- rep("wgs", length(sample_wgs))
  sample_data(physeq.wgs) <- sample_data(mysample_wgs)

  # Standardize taxonomic tables by renaming columns.
  tax_16s <- data.frame(physeq.16s@tax_table)
  colnames(tax_16s) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "unique")
  tax_table(physeq.16s) <- tax_table(as.matrix(tax_16s))

  tax_wgs <- data.frame(physeq.wgs@tax_table)
  colnames(tax_wgs) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "unique")
  tax_table(physeq.wgs) <- tax_table(as.matrix(tax_wgs))

  # Merge the two phyloseq objects.
  ps_merge <- merge_phyloseq(physeq.16s, physeq.wgs)

  # Process the merged taxonomic table: remove the 'unique' column and convert columns to factors.
  tax_merge <- data.frame(ps_merge@tax_table)
  tax_merge <- tax_merge %>% dplyr::select(-unique)
  tax_merge$Kingdom <- as.factor(tax_merge$Kingdom)
  tax_merge$Phylum  <- as.factor(tax_merge$Phylum)
  tax_merge$Class   <- as.factor(tax_merge$Class)
  tax_merge$Order   <- as.factor(tax_merge$Order)
  tax_merge$Family  <- as.factor(tax_merge$Family)
  tax_merge$Genus   <- as.factor(tax_merge$Genus)
  tax_table(ps_merge) <- tax_table(as.matrix(tax_merge))

  # Fix duplicate taxa using the BAZE package function.
  ps_merge_fixed <- BAZE::fix_duplicate_tax(ps_merge)

  # Generate a phylogenetic tree and update tip labels.
  tree_m <- BAZE::phy_to_tax(ps_merge_fixed)
  phylo <- tree_m@phylo
  tip <- phylo$tip.label
  phylo$tip.label <- gsub("^[a-z]__", "", tip)
  phy_tree(ps_merge_fixed) <- phylo

  return(ps_merge_fixed)
}


#' Calculate Overall Taxa Concordance Across Multiple Taxonomic Levels
#'
#' This function computes overall concordance metrics between WGS and 16S data across
#' several taxonomic levels (Phylum, Class, Order, Family, and Genus) for a given phyloseq object.
#' For each level, the function aggregates the data to relative abundances and calculates two metrics:
#' the Presence/Absence Concordance (PAC) and the Weighted Presence/Absence Concordance (WPAC). Each metric
#' is weighted by a branch length factor (BL), which is set to 1 for Genus, 2 for Family, 3 for Order,
#' 4 for Class, and 5 for Phylum. The overall metrics, IPAC and IWPAC, are computed by summing the numerators
#' and denominators across all levels and then calculating 1 minus the ratio.
#'
#' @param physeq A phyloseq object containing sample data and an OTU table. It is assumed that sample names
#'   in the phyloseq object correspond to subjects, and that each subject has associated samples with suffixes
#'   "-wgs" and "-16s" for WGS and 16S data, respectively.
#'
#' @return A data frame with one row per subject containing:
#'   \describe{
#'     \item{Subject}{The subject identifier.}
#'     \item{IPAC}{Overall inverse Presence/Absence Concordance, defined as
#'         \eqn{1 - \frac{\sum (BL \times |(wgs > 0) - (16s > 0)|)}{\sum (BL \times ((wgs > 0) + (16s > 0)))}}, summed over levels.}
#'     \item{IWPAC}{Overall inverse Weighted Presence/Absence Concordance, defined as
#'         \eqn{1 - \frac{\sum (BL \times |wgs - 16s|)}{\sum (BL \times 2)}}, summed over levels.}
#'   }
#'
#' @details
#' The function iterates over the taxonomic levels "Phylum", "Class", "Order", "Family", and "Genus".
#' For each level, the phyloseq object is aggregated and transformed to relative abundances. For each subject,
#' the function extracts count data for the corresponding WGS and 16S samples (using sample names with suffixes
#' "-wgs" and "-16s") and computes:
#' \itemize{
#'   \item \strong{PAC:} The numerator is \code{BL * sum(abs((wgs > 0) - (16s > 0)))}, and the denominator is
#'         \code{BL * sum((wgs > 0) + (16s > 0))}.
#'   \item \strong{WPAC:} The numerator is \code{BL * sum(abs(wgs - 16s))} and the denominator is \code{BL * 2}.
#' }
#' The branch length factor BL is assigned according to the taxonomic level:
#' \itemize{
#'   \item Genus: BL = 1
#'   \item Family: BL = 2
#'   \item Order: BL = 3
#'   \item Class: BL = 4
#'   \item Phylum: BL = 5
#' }
#' The overall inverse concordance metrics (IPAC and IWPAC) for each subject are computed as 1 minus the ratio
#' of summed numerators to summed denominators across all levels.
#'
#' @export
ALL_Concordance <- function(physeq) {
  levels <- c("Phylum", "Class", "Order", "Family", "Genus")
  subjects <- physeq@sam_data$sample_names

  all_level_weighted <- lapply(levels, function(level) {

    if (level == "Genus") {
      BL <- 1
    } else if (level == "Family") {
      BL <- 2
    } else if (level == "Order") {
      BL <- 3
    } else if (level == "Class") {
      BL <- 4
    } else if (level == "Phylum") {
      BL <- 5
    }

    ## Relative abundance table
    physeq_level <- aggregate_taxa(physeq, level = level)
    physeq_level <- transform_sample_counts(physeq_level, function(x) x / sum(x))
    count_tab <- physeq_level@otu_table

    PAC_Num <- list()
    PAC_Den <- list()
    WPAC_Num <- list()
    WPAC_Den <- list()

    for (subject in subjects) {
      wgs_id <- paste0(subject, "-wgs")
      sixs_id <- paste0(subject, "-16s")
      count_wgs <- count_tab[, wgs_id]
      count_sixs <- count_tab[, sixs_id]

      # PAC Calculation
      PAC_Num[[subject]] <- BL * sum(abs((count_wgs > 0) - (count_sixs > 0)))
      PAC_Den[[subject]] <- BL * sum((count_wgs > 0) + (count_sixs > 0))

      # Weighted PAC Calculation
      WPAC_Num[[subject]] <- BL * sum(abs(count_wgs - count_sixs))
      WPAC_Den[[subject]] <- BL * 2
    }

    results_df <- data.frame(Subject = subjects,
                             PAC_Num = unlist(PAC_Num), PAC_Den = unlist(PAC_Den),
                             WPAC_Num = unlist(WPAC_Num), WPAC_Den = unlist(WPAC_Den))
    return(results_df)
  })

  names(all_level_weighted) <- c("Phylum", "Class", "Order", "Family", "Genus")

  sum_column <- function(data_list, column_name) {
    column_data <- lapply(data_list, function(df) df[[column_name]])
    Reduce(`+`, column_data)
  }

  Sum_PAC_Num <- sum_column(all_level_weighted, "PAC_Num")
  Sum_PAC_Den <- sum_column(all_level_weighted, "PAC_Den")
  Sum_WPAC_Num <- sum_column(all_level_weighted, "WPAC_Num")
  Sum_WPAC_Den <- sum_column(all_level_weighted, "WPAC_Den")

  out_df <- data.frame(Subject = subjects,
                       IPAC = 1 - Sum_PAC_Num / Sum_PAC_Den,
                       IWPAC = 1 - Sum_WPAC_Num / Sum_WPAC_Den)

  return(out_df)
}



