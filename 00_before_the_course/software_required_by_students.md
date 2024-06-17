# Software required to be installed by students on their local computer <!-- omit from toc -->

## Table of contents <!-- omit from toc -->
- [FileZilla client](#filezilla-client)
- [R](#r)
- [RStudio](#rstudio)
- [R packages](#r-packages)

## FileZilla client
Follow the instructions given [here](https://filezilla-project.org) for the operating system on your machine.

## R
Follow the instructions given [here](https://cran.r-project.org/mirrors.html) for the operating system on your machine.

## RStudio
Follow the instructions given [here](https://posit.co/downloads/) for the operating system on your machine.

## R packages
Open RStudio, and execute these commands:
```R
# for packages in CRAN
install.packages(c("ade4", "adegenet", "corrplot", "data.table", "dplyr", "geosphere", "ggVennDiagram", "gplots", "ggplot2", "magrittr", "reshape2", "RColorBrewer", "splitstackshape", "tibble", "tidyverse", "varhandle", "vcfR", "vegan"))

# for packages in Bioconductor
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!("devtools" %in% installed.packages())){install.packages(devtools)}

library(devtools)

BiocManager::install("edgeR")

BiocManager::install("LEA")

BiocManager::install("goseq")

devtools::install_github("petrelharp/local_pca/lostruct")

BiocManager::install("qvalue")

devtools::install_github("whitlock/OutFLANK")

BiocManager::install("StAMPP") # install via BioConductor

```
>Answer "yes" if asked to install additional libraries.

Check that the packages were installed correctly:
```R
library(ade4)
library(adegenet)
library(corrplot)
library(data.table)
library(dplyr)
library(edgeR)
library(LEA)
library(lostruct)
library(geosphere)
library(ggplot2)
library(ggVennDiagram)
library(goseq)
library(gplots)
library(magrittr)
library(tibble)
library(tidyverse)
library(OutFLANK)
library(RColorBrewer)
library(reshape2)
library(splitstackshape)
library(StAMPP)
library(varhandle)
library(vcfR)
library(vegan)

```
If no error messages appear, you are all set!

If you get error messages examine them carefully, to identify which additional packages might be missing/outdated, and need to be installed (answer YES when the prompt asks to install additional packages). 

Follow equivalent commands to install missing packages as shown before (directly from CRAN, Bioconductor, or a dedicated GitHub repository).

If you encounter further challenges, perhaps you are not alone and someone else already found a solution in the internet. 

You can also direct your questions to the dedicated Slack channel of the course for this purpose. For this, please make sure to document the problem as much as you can. For example, describe what the problem is, document it (a screenshot would be great), and provide context (e.g., which operating system you are using, which R version you have, which commands were used, etc.). We are happy to help!

Best of luck!
