# Using concordance

*concordance* is a comprehensive package including study design for your microbiome data analysis study, and a concordance measurement method to evaluate the concordance between the same datasets from different sequencing platforms (16S rRNA and shotgun data).

If you want to design a new study to explore the relationship between microbiome and human disease, *concordance* package provide power calculation methods for specific sample size range and it can help researcher calculate potential cost for different sequencing paltforms with expected power.

## Table of Contents

-   [Citation](#citation)
-   [Installation](#installation)
-   [Dependent packages](#dependencies)
-   [Function Details](#function-details)
    -   [sampleCalculationRange()](#sample-calculation)
    -   [powerCalculation()](#power-calculation)
    -   [priceCalculation()](#price-calculation)
    -   [concordance_overall()](#concordance-overall)
    -   [concordance_split()](#concordance_split)

## Citation

## Installation

You can install **concordance** R package by:

```{r}
# # Install devtools if you haven't already
# install.packages("devtools")
# library(devtools)
# 
# # Install the package from GitHub
# install_github("bioscinema/concordance")
# 
# # Load your package
# library(concordance)
```

## Dependent packages 

```{r}
library(knitr)

# Create a data frame with dependencies and functions
dependencies <- data.frame(
  Package = c("ggplot2", "phyloseq", "BAZE", "MASS", "dplyr", "stats", "vegan"),
  Functions = c("All",
                "All",
                "phy_to_tax",
                "rnegbin",
                "filter, select",
                "median, quantile, rmultinom, rnorm",
                "adonis2, vegdist")
)

# Display the table using knitr::kable
kable(dependencies, caption = "Dependencies for the Concordance Package")

```

## Function details

### sampleCalculationRange() 

```{r sampleCalculationRange sample,echo=TRUE, eval=FALSE}
library(concordance)

# Evaluate power for sample sizes from 50 to 200 (in increments of 5)
# for gut data using the 16S (amplicon) method.
result <- sampleCalculationRange(
  lower_bound = 50,
  upper_bound = 200,
  data.type = "gut",
  method = "amp",
  nSim = 100,       # Number of simulation replicates
  nOTU = 300,       # Number of OTUs to simulate
  diff.otu.pct = 0.1  # Proportion of OTUs that are differentially abundant
)

# Display the results:
print(result$data)
print(result$plot)

```

### powerCalculation

```{r}
# Load the package containing powerCalculation (adjust if needed)
library(concordance)

# Example: Evaluate power for gut data using amplicon sequencing ("amp")
# with 50 simulation replicates, a significance level of 0.05,
# 100 samples per simulation, and 200 OTUs simulated with 10% differential abundance.
result <- powerCalculation(
  data.type = "gut",     # Choose data type: "gut", "oral", or "infant"
  method = "amp",        # Choose sequencing method: "amp" (16S) or "wgs" (WGS)
  nSim = 50,             # Number of simulation replicates
  alpha = 0.05,          # Significance level
  nSam = 100,            # Number of samples per simulation replicate
  nOTU = 200,            # Number of OTUs to simulate
  diff.otu.pct = 0.1      # Proportion of OTUs that are differentially abundant
  # All other parameters will use their default values.
)

# Print the overall estimated power (proportion of replicates with significant results)
print(result$overall_power)

# Optionally, print the detection indicator for each replicate
print(result$replicate_power)
```

### priceCalculation

```{r}
library(concordance)

# Minimal example call to priceCalculation:
price_estimates <- priceCalculation(
  target_power = 0.8,       # Desired power (80%)
  lower_bound = 50,         # Minimum sample size to consider
  upper_bound = 200,        # Maximum sample size to consider
  data.type = "gut",        # Data type ("gut", "oral", or "infant")
  nOTU=300
)

# Print the resulting data frame
print(price_estimates)

```


### concordance_overall()
```{r}
result <- concordance_overall(phy)
```

### concordance_split()
```{r}
result <- concordance_split(phy1,phy2)
```
