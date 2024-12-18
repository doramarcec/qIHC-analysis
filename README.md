# High-throughput image analysis of quantitative immunohistochemistry data reveals age-dependent changes in mitochondrial and extracellular matrix proteins in multiple human tissues

The process of ageing induces physiological changes at cellular and molecular levels, often leading to complex chronic diseases. Whilst a large body of literature describes age-induced changes at a tissue level, not as much is known about the impact of ageing on the tissue microenvironment, defined by its extracellular matrix. In this project, I investigated the changes in matrisome and mitochondrial proteome expression with age across 44 human tissue types. 

This project was my Bachelor's thesis, and its purpose was to identify a subset of mitochondrial and extracellular matrix proteins whose expression levels significantly change during ageing and curate a list of proteins that are differentially expressed across multiple age groups, in a single or multiple tissue types. The objectives to achieve this include i) development of a differential expression workflow which will perform a thorough statistical testing for significance between staining intensity levels, representative of protein expression, across age groups, and ii) the implementation of the developed workflow on the quantitative data from mitochondrial proteome and the matrisome to identify the proteins of interest.

**Data structure**
Data mining techniques were used by the host lab to obtain i) 153,919 pathologist-annotated immunohistochemistry images spanning 44 tissue types and ii) patient metadata from the Human Protein Atlas database. The images were batch processed using a custom MATLAB script, generating a dataset of the following structure:

Gene | Tissue_1 | Antibody | Image | Sex | Age | Tissue_2 | Subject | Intensity | SD | Max | 75 | Median | 25 | Min | Quantity | Staining
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
BGN | skin | CAB003678 | http://images.proteinatlas.org/3678/10249_B_6_8.jpg | F | 15 | Normal tissue | 2027 | 62.38 | 22.72 | 225 | 73 | 57 | 46 | 10 | 19.49 | 12.16

## Data analysis in R

Before importing the data, we will load the libraries required in my analysis and set my custom theme set for visualisations.

```
library(tidyverse)
library(here)
library(openxlsx)
library(rstatix)
library(ggpubr)
theme_set(
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13),
    title = element_text(size = 14),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5))
  )
library(ggsci)
library(scales)
library(correlation)
```
Next, we import the data using the 'here' package, followed by filtering and tidying the data.

```
full_ECM <- read.xlsx(here("data", "ECM_data.xlsx"))
full_ECM

# Filter the data to only contain healthy tissue samples
ECM <- full_ECM %>%
  remove_rownames() %>%
  dplyr::filter(Tissue_type_2 == "Normal tissue, NOS (M-00100)")

# Remove unneccessary columns
ECM$Unknown_code_1 <- NULL
ECM$Unknown_code_2 <- NULL

# Tidy the antibody names
tidy_antibody_names <- function(ECM) {
  ECM$Antibody_name <- gsub("^Antibody ", "", ECM$Antibody_name)
  return(ECM)
}
ECM <- tidy_antibody_names(ECM)
```

This filtering and tidying generates a female-dominated data set containing the data from 1,268 patients and 153,919 immunohistochemistry image samples, spanning 44 tissue types and including 734 genes and 1,206 antibodies. 

Now that the data is tidy, we are going to group it according to age and remove any missing values. 

```
ECM$Age2 <- as.numeric(as.character(ECM$Age))

ECM <- ECM %>%
  mutate(Age2 = case_when(
    Age2 < 40 ~ "< 40", 
    Age2 >= 40 & Age2 < 60 ~ "40-60", 
    Age2 >= 60 ~ "≥ 60" 
  )) %>%
  drop_na()

ECM %>%
  count(Age2)

level_order <- c('< 40', '40-60', '≥ 60')
```
The remainder of the analysis will center around inferential statistics, aiming to identify whether the protein expression across the three age categories statistically differs. 

We start off with normality testing to know which statistical tests will be the most appropriate in the upcoming steps. As the sample size was to big for formal normality test (i.e. Shapiro-Wilk test, with a sample size limit of 5000), density and QQ plots were used. 

```
dp <- ECM %>%
  select(Age2, Intensity) %>%
  ggplot() + aes(x = Intensity, colour = Age2) + geom_density() + 
  scale_colour_npg() + scale_fill_npg() + ggtitle("Intensity distribution") + ylab("Density")

qqp <- ECM %>%
  select(Age2, Intensity) %>%
  ggplot() + aes(sample = Intensity, colour = Age2) + geom_qq() + geom_qq_line() + 
  scale_colour_npg() + scale_fill_npg() + ggtitle("Q-Q plot of intensities") + 
  ylab("Sample quantiles") + xlab("Theoretical quantiles")

ggarrange(
  dp1, qqp1,
  common.legend = FALSE, legend = "top"
  )

ggsave(here("output/density-and-qq2.tiff"), height = 4, width = 9)
```

This produced the following visual:
![density-and-qq](https://github.com/user-attachments/assets/14758f9d-0940-413e-9b52-62b352b3db71)

