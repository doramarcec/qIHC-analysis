# High-throughput image analysis of quantitative immunohistochemistry data reveals age-dependent changes in mitochondrial and extracellular matrix proteins in multiple human tissues

The process of ageing induces physiological changes at cellular and molecular levels, often leading to complex chronic diseases. Whilst a large body of literature describes age-induced changes at a tissue level, not as much is known about the impact of ageing on the tissue microenvironment, defined by its extracellular matrix. In this project, I investigated the changes in matrisome and mitochondrial proteome expression with age across 44 human tissue types. 

This project was the basis of my Bachelor's thesis. The host lab used data mining to obtain i) 153,919 pathologist-annotated immunohistochemistry images and ii) corresponding patient metadata from the Human Protein Atlas (HPA) database. To my knowledge, the patient metadata from the HPA database, including age and sex, has never been explored in the context of ageing, making this project the first exploration of the HPA data in an ageing study. 

The purpose of the project was to identify a subset of mitochondrial and extracellular matrix proteins whose expression levels significantly change during ageing and curate a list of proteins that are differentially expressed across multiple age groups, in single or multiple tissue types. The objectives to achieve this include i) the development of a differential expression workflow which will perform thorough statistical testing for significance between staining intensity levels, representative of protein expression, across age groups, and ii) the implementation of the developed workflow on the quantitative data from mitochondrial proteome and the matrisome to identify the proteins of interest.

**Data structure**
The images were batch-processed using a custom MATLAB script, generating a dataset of the following structure:

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
The remainder of the analysis will centre around inferential statistics, aiming to identify whether the protein expression across the three age categories statistically differs. 

We begin with normality testing to understand which statistical tests are the most appropriate for the upcoming steps. As the sample size was too big for a formal normality test (i.e. Shapiro-Wilk test, with a sample size limit of 5000), density and QQ plots were used. 

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

This produced the following figure:
![density-and-qq](https://github.com/user-attachments/assets/b8ab2651-ded0-4035-b8d0-6feda1c2ddf4)

Based on the output above, we can deduce that the data is not normally distributed. It is positively skewed, as visible in the density plot, complemented by the upward curve shape in the QQ plot, representative of positively skewed data. 

This means that statistical testing throughout this project will have to be non-parametric. My analysis and the workflow I've developed for this type of data were guided by the following questions:
1. Is there a significant difference in intensity between different age groups across tissue types?
2. If so, which age groups significantly statistically differ across tissue types?
3. Among the tissue types with significant differences, which genes do significantly statistically differ in intensity?
4. Which age groups significantly statistically differ in each of those genes?

## 1. Is there a significant difference in intensity between different age groups across tissue types?
The optimal statistical test to answer this question was the Kruskal-Wallis test. First, we are going to define the function that will perform the test, followed by creating a nested data frame grouped by our column of interest. Then, we perform the test on the data and save the output into the nested data frame. Lastly, we unnest the data, move the columns of interest into a results data frame, and filter the data based on our desired significance threshold.

```
doKruskalTest <- function(x) {
  x %>%
    kruskal_test(Intensity~Age2) %>%
  print(kruskalTest)
}

a <- ECM %>%
  group_by(Tissue_type_1) %>%
  nest()

a$data[[1]]

a <- ECM %>%
  group_by(Tissue_type_1) %>%
  nest() %>%
  mutate(Kruskal = map(data, doKruskalTest))

a$Kruskal[[1]]

# Tidy the data frames
a <- a %>%
  unnest()

kruskalResults <- a %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, statistic, df, p) %>%
  unique() %>%
  dplyr::filter(p < 0.001) %>%
  group_by(Tissue_type_1) %>%
  arrange(.by_group = TRUE)
kruskalResults

kruskalResults %>%
  count(Tissue_type_1) 
```
This testing identified that 34/44 tissue types have a statistically significant difference in intensity levels, representative of protein expression, between the age groups. Knowing this, we can move on to the next part of the workflow i.e. the next question!

## 2. Which age groups significantly statistically differ across tissue types?
To answer this question, we will perform a post-hoc analysis of the Kruskall-Wallis test results, for which Dunn's test is used, using a similar approach as above. 

Note how we are using *a* data frame, generated in the previous code chunk, to perform Dunn's test only on the data stored in that data frame, before storing it into a new results data frame, *b*. We will keep using that same approach throughout this workflow. 

```
doDunnsTest <- function(x) {
    x %>%
    dunn_test(Intensity~Age2) %>%
  print(dunnsTest)
  }

b <- a %>%
  group_by(Tissue_type_1) %>%
  dplyr::filter(p < 0.001) %>%
  nest()

b$data[[1]]

b <- a %>%
  group_by(Tissue_type_1) %>%
  dplyr::filter(p < 0.001) %>%
  nest() %>%
  mutate(Dunn = map(data, doDunnsTest))

b$Dunn[[1]]

# Tidy the data frames
b <- b %>%
  unnest(data) %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, Dunn) %>%
  unnest(Dunn)

dunnResults <- b %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, group1, group2, n1, n2, statistic, p, p.adj.signif) %>%
  dplyr::filter(p < 0.001) %>%
  group_by(Tissue_type_1) %>%
  arrange(.by_group = TRUE)
dunnResults 
```

By doing so, we get the 219,703-row tibble with the following output structure:
Tissue | Gene | Antibody | Sex | Age | Age2 | Subject | Intensity | SD | group1 | group2 | n1 | n2 | statistic | p | p.adj.sifnif 
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- 
adipose+tissue | BGN | CAB003678 | M | 26 | < 40 | 1995 | 62.74 | 25.57 | < 40 | 40-60 | 1297 | 1189 | 5.17 | 2.33e-07 | ****

Now that we know which tissues have significant differences in protein expression and between which age groups these differences occur, we want to know which genes and proteins contribute to it.

## 3. Among the tissue types with significant differences, which genes do significantly statistically differ in intensity?
Here, we will use the same Kruskall-Wallis function we defined in the first step, but this time we will apply it to the Gene column instead of the Tissue column, and include the filtering step based on our set significance threshold. 

```
c <- b %>%
  group_by(Gene_name) %>%
  dplyr::filter(p < 0.001) %>%
  nest()

c$data[[1]]

c <- b %>%
  group_by(Gene_name) %>%
  dplyr::filter(p < 0.001) %>%
  nest() %>%
  mutate(Kruskal = map(data, doKruskalTest))

c$Kruskal[[1]]

# Tidy the data frames
c <- c %>%
  unnest(data) %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, Kruskal) %>%
  unnest(Kruskal)

kruskalGenes <- c %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, n, statistic, df, p) %>%
  dplyr::filter(p < 0.001)
kruskalGenes
```
Running this chunk of code identified 73 genes (out of 734) which significantly differ in intensity across two or more age groups, in the 34 tissues identified above, which brings us to the last part of our workflow. 

## 4. Which age groups significantly statistically differ in each of those genes?
Once again, as we used the Kruskall-Wallis test in the previous step, we need to use Dunn's test in the post-hoc analysis to identify the specific genes which significantly differ.

```
d <- c %>%
  group_by(Gene_name) %>%
  dplyr::filter(p < 0.001) %>%
  nest()

d$data[[1]]

d <- c %>%
  group_by(Gene_name) %>%
  dplyr::filter(p < 0.001) %>%
  nest() %>%
  mutate(Dunn = map(data, doDunnsTest))

d$Dunn[[1]]

# Tidy the data frames
d <- d %>%
  unnest(data) %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, Dunn) %>%
  unnest(Dunn)

dunnResults <- d %>%
  select(Tissue_type_1, Gene_name, Antibody_name, Sex, Age, Age2, Subject_ID, Intensity, SD, group1, group2, n1, n2, statistic, p, p.adj.signif) %>%
  dplyr::filter(p < 0.001) %>%
  group_by(Gene_name) %>%
  arrange(.by_group = TRUE)
dunnResults 
```
Finally, this generates a 39,265-row tibble of the following format:
Tissue | Gene | Antibody | Sex | Age | Age2 | Subject | Intensity | SD | group1 | group2 | n1 | n2 | statistic | p | p.adj.sifnif 
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- 
adipose+tissue | ABAT | HPA041528 | F | 26 | < 40 | 2744 | 68.91 | 28.46 | < 40 | 40-60 | 104 | 100 | 4.73 | 2.24e-07 | ****

From this tibble, we start the exploration and visualisation of differentially expressed proteins. By counting the number of genes in the descending order as below, we discover that matrix metalloproteinase 9 (MMP9) is the protein that is differentially expressed in the greatest number of samples. 

```
dunnResults %>%
  ungroup() %>%
  count(Gene_name) %>%
  arrange(desc(n))
```
To explore the data further, we will separate the dunnResults tibble into new data frames where we group the tissues with similar functions or that make a single organ (e.g. basal ganglia, cerebellum, cerebral cortex and hippocampal formation into the *brain* data frame) as such:

```
metabolicRes <- dunnResults %>%
  dplyr::filter(Tissue_type_1 == c("spleen", "liver"))
```

Next, we will count the number of genes in the newly created data frames to identify the genes differentially expressed in the greatest number of image samples. 
```
metabolicRes %>%
  count(Gene_name) %>%
  arrange(desc(n))
```

Visualisation of the relationship between intensity and age in each of the 73 genes would be a long and laborious task, so I defined a function that automatically plots the linear regression. 
```
plotLinReg <- function(x) {
  x %>%
    ggplot() + aes(x = Age, y = Intensity, colour = Antibody_name) + 
    geom_point() + facet_wrap(~Tissue_type_1, scales = "free_x") + geom_smooth(method = "lm")
}
```

The output is going to be stored in a newly created, nested data frame, and we are going to write a loop to display all linear regression plots simultaneously. Note that we are colouring the plot by the antibody names. That is because the purpose of this step is to identify whether or not the different antibodies that were used follow the same pattern or deviate from each other. 
```
metabolicLinReg <- metabolicRes %>%
  group_by(Gene_name) %>%
  nest()

metabolicLinReg$data[[1]]

metabolicLinReg <- metabolicRes %>%
  group_by(Gene_name) %>%
  nest() %>%
  mutate(LinReg = map(data, plotLinReg))

metabolicLinReg$LinReg[[1]]

# Display all the linear regression plots at once
for (a in seq_along(metabolicLinReg$LinReg)) {
  print(metabolicLinReg$LinReg[[a]])
}
```

During this step, we are also going to filter out any antibodies that produced outliers in the data. 



