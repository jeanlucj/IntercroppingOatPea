# Randomization for the winter oat expt, fall planting 2023
# Four sets of lines:
# 1. Lines replicated in NY and KY
# 2. Lines unreplicated in KY but replicated in NY
# 3. Lines replicated in NY and not present in KY
# 4. Lines unreplicated in NY and not present in KY

here::i_am("code/WinterOatPea2023.R")
library(AlgDesign)
library(tidyverse)
set.seed(1234)

# Approach #1
# Set this all up in a data.frame with line and state columns and
# apply AlgDesign optBlock.

# Get the oat line names from 2024_winter_oat_planning_02092023.csv
oatPlanning <- read_csv(
  here::here("data", "2024_winter_oat_planning_02092023.csv"))
oatPlanning <- dplyr::select(oatPlanning, source, Genotype,
                             oat_reps, "JLJ--30Aug23_Keep4KY")
colnames(oatPlanning) <- c("source", "accession", "nRepsNY",
                           "nRepsKY")

# For some reason hangs on to a lot of NA rows
oatPlanning <- dplyr::filter(oatPlanning, !is.na(source))
oatPlanning <- dplyr::mutate(oatPlanning,
                             nRepsKY=if_else(is.na(nRepsKY), 0, nRepsKY))
# There are some slaches in the Nobel line names. Replacing "/" with "_"
oatPlanning <- dplyr::mutate(oatPlanning,
  accession=gsub("/", "_", accession, fixed=T))

# NCSU has replicates named like "Horizon 578-1" and "Horizon 578-2"
# Collapse those
ncsuAcc <- filter(oatPlanning, source == "NCSU-ALT")
acc <- strsplit(ncsuAcc$accession, "-")
noLast1or2 <- function(svec){
  ncat <- length(svec)
  if (dplyr::last(svec) %in% c("1", "2")) ncat <- ncat - 1
  return(paste(svec[1:ncat], collapse="-"))
}
acc <- sapply(acc, noLast1or2)
ncsuAcc$accession <- acc
ncsuAcc <- dplyr::group_by(ncsuAcc, accession) %>%
  dplyr::summarize(nRepsNY=sum(nRepsNY)) %>%
  dplyr::mutate(source="NCSU-ALT", nRepsKY=0) %>%
  dplyr::relocate(source, accession, nRepsNY, nRepsKY)
oatPlanning <- dplyr::filter(oatPlanning, source != "NCSU-ALT") %>%
  dplyr::bind_rows(ncsuAcc)

# The reps are still given as plots.  I want to pair it down to "split plots"
# So divide by two and also I don't want to have > 2 reps
#oatPlanning <- dplyr::mutate(oatPlanning,
#                             nRepsNY=pmin(2, floor(nRepsNY/2)),
#                             nRepsKY=pmin(2, floor(nRepsKY/2)))
oatPlanning <- dplyr::mutate(oatPlanning,
                             nRepsNY=floor(nRepsNY/2),
                             nRepsKY=floor(nRepsKY/2))

# Six barley lines from Mark's program will be included. Each replicated in NY
# I'm calling one "NY_expt" because Mark/David haven't told us what yet
barleyNames <- c("Lightning", "Calypson", "Scala",
                 "Saturn", "Violetta", "NY_Expt")
barleyNames <- paste0(barleyNames, "Barley")
blyPlanning <- tibble(source="Barley", accession=barleyNames,
                      nRepsNY=2, nRepsKY=0)
oatPlanning <- dplyr::bind_rows(oatPlanning, blyPlanning)

# To make this pretty, I'm going to take out NF13-4173-4_6 which has what look
# like two fullsibs in there and put in 8 plots of monocrop pea
oatPlanning <- dplyr::filter(oatPlanning, accession != "NF13-4173-4_6")
peaPlanning <- tibble(source="Pea", accession="BlazePea",
                      nRepsNY=8, nRepsKY=0)
oatPlanning <- dplyr::bind_rows(oatPlanning, peaPlanning)

# How many plots are we talking about in NY?
nPlots <- sum(oatPlanning$nRepsNY)

# Make the data.frame for optBlock
plotsNY <- dplyr::tibble(
  accession=rep(oatPlanning$accession, oatPlanning$nRepsNY),
  state="NY") %>%
  mutate(crop=dplyr::if_else(grepl("Barley", accession), "Barley", "Oat")) %>%
  mutate(crop=dplyr::if_else(grepl("BlazePea", accession), "Pea", crop))

plotsKY <- dplyr::tibble(
  accession=rep(oatPlanning$accession, oatPlanning$nRepsKY),
  state="KY", crop="Oat")
plotsAll <- dplyr::bind_rows(plotsNY, plotsKY)

# Ideally, I would do AlgDesign with plotsAll, while forcing the KY plots to be
# fixed. But I don't know how to do that. I would probably need to write my own
# optBlock function.
# Instead, I will apply it to the NY plots and just concatenate with the KY
# plots. Hmmm: should I have one NY block that is the same as the KY block? No.
# I don't think I need to have that.
blkOP <- AlgDesign::optBlock(frml=as.formula("~ accession"),
                  withinData=plotsNY,
                  blocksizes=rep(24, 8))

# Assemble the design
winterOatDesign <- plotsKY %>%
  dplyr::arrange(sample(nrow(plotsKY))) %>%
  dplyr::mutate(block="Block1") %>%
  dplyr::relocate(block, .before=accession)
for (blockNum in 1:8){
  thisBlock <- as_tibble(blkOP$Blocks[[blockNum]]) %>%
    dplyr::arrange(sample(nrow(blkOP$Blocks[[blockNum]]))) %>%
    dplyr::mutate(block=paste0("Block", blockNum+1)) %>%
    dplyr::relocate(block, .before=accession)
  winterOatDesign <- dplyr::bind_rows(winterOatDesign, thisBlock)
}

# Some quick checks
monoPea <- dplyr::filter(winterOatDesign, accession=="BlazePea")
table(monoPea$block)
hasBarley <- dplyr::slice(winterOatDesign, grep("Barley", accession))
table(hasBarley$block)

readr::write_csv(winterOatDesign, here::here("output", "WinterOatDesign.csv"))
saveRDS(winterOatDesign,
        here::here("output", "WinterOatDesign.rds"))
