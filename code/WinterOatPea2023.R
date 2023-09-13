# Randomization for the winter oat expt, fall planting 2023
# Four sets of lines:
# 1. Lines replicated in NY and KY
# 2. Lines unreplicated in KY but replicated in NY
# 3. Lines replicated in NY and not present in KY
# 4. Lines unreplicated in NY and not present in KY
#
# The sources are
# 1. Noble population from Lauren
# 2. North Carolina Advanced Trial from Paul Murphy
# 3. Uniform Winter Oat Trial (UWOT) from Tan Tuong
# 4. "Winter Hayden" from Klaas Martens
# The first two are in 2024_winter_oat_planning_02092023.csv
# In that file there are also stand-ins for the UWOT but I will
# replace those from UON24_Entries.csv
# I will manually add in Klaas Martens' line

here::i_am("code/WinterOatPea2023.R")
library(AlgDesign)
library(tidyverse)
ss <- floor(runif(1, 10000, 1000000))
set.seed(521891) # I dink around with the seed until optBlock gives me what I want

# Set this all up in a data.frame with line and state columns and
# apply AlgDesign optBlock.

# Accessions from the UWOT
uwot <- read_csv(
  here::here("data", "UON24_Entries_(11Sep23).csv"), skip=1)
uwot <- dplyr::select(uwot, EXPT, DESIG) %>%
  mutate(nRepsNY=6, nRepsKY=0)
colnames(uwot)[1:2] <- c("source", "accession")

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
# There are some slashes in the Nobel line names. Replacing "/" with "_"
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

# Find out what uwot accessions are in oatPlanning
uwotInOP <- dplyr::filter(uwot, accession %in% oatPlanning$accession)
# This looks right so no manual curation here
# NOTE: one accession in the UWOT is duplicated
rm(uwotInOP)

oatPlanning <- dplyr::filter(oatPlanning, source != "UWOT") %>%
  dplyr::bind_rows(uwot) %>%
  dplyr::group_by(accession) %>%
  dplyr::summarize(source=dplyr::last(source),
                   nRepsNY=sum(nRepsNY),
                   nRepsKY=dplyr::last(nRepsKY)) %>%
  dplyr::relocate(source, accession, nRepsNY, nRepsKY) %>%
  dplyr::bind_rows(tibble(source="Martens", accession="Winter Hayden",
                          nRepsNY=6, nRepsKY=0)) %>%
  dplyr::arrange(source)

#### Done with entry curation work

# The reps are still given as plots.  I want to pair it down to "split plots"
# So divide by two and also I don't want to have > 2 reps
oatPlanning <- dplyr::mutate(oatPlanning,
                             nRepsNY=pmin(2, floor(nRepsNY/2)),
                             nRepsKY=pmin(2, floor(nRepsKY/2)))
# oatPlanning$nRepsNY[which(oatPlanning$source=="Martens")] <- 3

#oatPlanning <- dplyr::mutate(oatPlanning,
#                             nRepsNY=floor(nRepsNY/2),
#                             nRepsKY=floor(nRepsKY/2))

# Six barley lines from Mark's program will be included. Each replicated in NY
# I'm calling one "NY_expt" because Mark/David haven't told us what yet
barleyNames <- c("Lightning", "Calypson", "Scala",
                 "Saturn", "Violetta", "NY_Expt")
barleyNames <- paste0(barleyNames, "Barley")
blyPlanning <- tibble(source="Barley", accession=barleyNames,
                      nRepsNY=2, nRepsKY=0)
oatPlanning <- dplyr::bind_rows(oatPlanning, blyPlanning)

# Add in monoculture peas
peaPlanning <- tibble(source="Pea", accession="BlazePea",
                      nRepsNY=9, nRepsKY=0)
oatPlanning <- dplyr::bind_rows(oatPlanning, peaPlanning)

# How many plots are we talking about in NY?
nPlots <- sum(oatPlanning$nRepsNY)
# This all generates 180 plots.  So the plan will be for 9 blocks of 20 plots
nBlksNY <- 9
nPlotsPerBlk <- 20

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
                  blocksizes=rep(nPlotsPerBlk, nBlksNY))

getOatSource <- function(accession){
  return(oatPlanning$source[which(oatPlanning$accession == accession)])
}
# Assemble the design
winterOatDesign <- plotsKY %>%
  dplyr::arrange(sample(nrow(plotsKY))) %>%
  dplyr::mutate(block="Block1") %>%
  dplyr::mutate(source=sapply(accession, getOatSource)) %>%
  dplyr::relocate(block, .before=accession)
for (blockNum in 1:nBlksNY){
  thisBlock <- as_tibble(blkOP$Blocks[[blockNum]]) %>%
    dplyr::arrange(sample(nrow(blkOP$Blocks[[blockNum]]))) %>%
    dplyr::mutate(block=paste0("Block", blockNum+1)) %>%
    dplyr::mutate(source=sapply(accession, getOatSource)) %>%
    dplyr::relocate(block, .before=accession)
  winterOatDesign <- dplyr::bind_rows(winterOatDesign, thisBlock)
}

# Some quick checks
monoPea <- dplyr::filter(winterOatDesign, accession=="BlazePea")
table(monoPea$block)
hasBarley <- dplyr::slice(winterOatDesign, grep("Barley", accession))
table(hasBarley$block)
print(ss)

readr::write_csv(winterOatDesign, here::here("output", "WinterOatDesign.csv"))
saveRDS(winterOatDesign,
        here::here("output", "WinterOatDesign.rds"))
