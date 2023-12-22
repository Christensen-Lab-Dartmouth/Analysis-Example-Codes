# title: "Running REMP For LINE-1 Prediction Example"
# author: "Curtis Petersen"
# date: "10/23/2019"


# This script assumes that you have matrix of normalized beta values and that the array type is 450k. If you have EPIC make sure to change the arrayType to 'EPIC' where it says '450k'. 
# It takes those beta values and estimates what the methylation status is for each of the defined LINE-1 sites.
# Finally, we wrangle this data, filtering out those elements that might be low quality, summarize accross 2 CpGs per element, join it with the LINE-1 element annotation. 
# All of these are saved as data frames in an RData file to be used in analysis and visualization

# Load Packages & Initialize
library(REMP)
library(tidyverse)
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
N_Cores <- detectCores() - 1

# Load Data
load("../Data/Example_Normalized_Betas.RData")


# LINE-1 Prediction

## Initalize REMP For LINE-1: Specify parcel with array and compenents to return
rempL1 <- initREMP(arrayType = '450k', REtype = 'L1')

## Run REMP for LINE-1 
L1Results <- remp(methyDat = Beta_Matrix, REtype = 'L1', parcel = rempL1, seed = 1234, ncore = N_Cores)
# THIS TAKES A LONG TIME! ~ 10- secionds per sample running locally...
# On EPIC this can take ~ 300 seconds per sample....

# Wrangle LINE-1 Data

# Filter results
L1Aggregate <- rempTrim(L1Results, threshold = 1.7, missingRate = 0.2)

# Aggregate results
L1Aggregate <- rempAggregate(L1Aggregate, NCpG = 2)

# Add genomic region annotation
L1Aggregate <- decodeAnnot(L1Aggregate)

# Extract annotation
L1Annot <- as.data.frame(rempAnnot(L1Aggregate))

# Extract beta 
L1Betas <- as.data.frame(rempB(L1Aggregate))


# Save results
save(rempL1, L1Results, L1Aggregate, L1Annot, L1Betas, file = '../Data/Repeat_Elements_LINE1.RData')