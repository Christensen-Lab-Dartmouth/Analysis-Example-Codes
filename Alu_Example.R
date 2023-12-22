# title: "Running REMP For Alu Prediction Example"
# author: "Curtis Petersen"
# date: "10/23/2019"


# This script assumes that you have matrix of normalized beta values and that the array type is 450k. If you have EPIC make sure to change the arrayType to 'EPIC' where it says '450k'. 
# It takes those beta values and estimates what the methylation status is for each of the defined Alu sites.
# Finally, we wrangle this data, filtering out those elements that might be low quality, summarize accross 2 CpGs per element, join it with the Alu element annotation. 
# All of these are saved as data frames in an RData file to be used in analysis and visualization.

# Load Packages & Initialize
library(REMP)
library(tidyverse)
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
N_Cores <- detectCores() - 1

# Load Data
load("Example_Normalized_Betas.RData")


# Alu Prediction

## Initalize REMP For Alu: Specify parcel with array and compenents to return
rempAlu <- initREMP(arrayType = '450k', REtype = 'Alu')

## Run REMP for Alu 
AluResults <- remp(methyDat = Beta_Matrix, REtype = 'Alu', parcel = rempAlu, seed = 1234, ncore = N_Cores)
# THIS TAKES A LONG TIME! ~ 80 secionds per sample running locally...
# On EPIC this can take ~ 300 seconds per sample....

# Wrangle Alu Data

# Filter results
AluAggregate <- rempTrim(AluResults, threshold = 1.7, missingRate = 0.2)

# Aggregate results
AluAggregate <- rempAggregate(AluAggregate, NCpG = 2)

# Add genomic region annotation
AluAggregate <- decodeAnnot(AluAggregate)

# Extract annotation
AluAnnot <- as.data.frame(rempAnnot(AluAggregate))

# Extract beta values
AluBetas <- as.data.frame(rempB(AluAggregate))


# Save results
save(rempAlu, AluResults, AluAggregate, AluAnnot, AluBetas, file = '../Data/Repeat_Elements_Alu.RData')