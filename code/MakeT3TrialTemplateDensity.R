# Make template to upload the Density design to T3 based on AlgDesignDensity.R
# The key output from there is densDesign2
# To not mess that up, I'm going to copy it to T3Design here
library(tidyverse)

T3Design <- densDesign2
oatAccDens <- c(checkOats$germplasmName,
                oatToTest %>% filter(entryNumber %in% c(41, 11, 18, 32)) %>%
                  pull(germplasmName))
names(oatAccDens) <- c(1:4, 11, 18, 32, 41)
peaAccDens <- c("AUSTRIAN FIELD PEAS", "delta field pea",
                "DS Admiral", "Icicle Winter Peas",
                "ND Victory", "Organic 4010")
names(peaAccDens) <- c(1:4, 11, 5)

locationNames <- c("NY")
breeding_programs <- c("Cornell University")
trialNameBase <- "_OatPeaDensity_2023_"
trialNameBreedProg <- c("Cornell")
trialNameLocation <- c("Ithaca")
trial_names <- paste0(trialNameBreedProg, trialNameBase, trialNameLocation)
locations <- c("Ithaca, NY - Caldwell")
planting_dates <- c("2023-April-13")
plot_widths <- c(1.26)
plot_lengths <- c(3)

trialTemplateHeader <- "trial_name	breeding_program	location	year	design_type	description	trial_type	plot_width	plot_length	field_size	planting_date	harvest_date	plot_name	accession_name	plot_number	block_number	is_a_control	rep_number	range_number	row_number	col_number	seedlot_name	num_seed_per_plot	weight_gram_seed_per_plot	is_private"
cn <- strsplit(trialTemplateHeader, "\t")

# FIELD LAYOUTS
# New York
nyMap <- read_csv(here::here("data", "NY_Oat-Pea23map.csv"), col_names=F,
                  skip=4, n_max=5) %>%
  dplyr::select(where(is.numeric))
rowCol <- sapply(101:200, function(p) which(nyMap==p, arr.ind=T)) %>% t
rowCol[,1] <- nrow(nyMap) - rowCol[,1] + 1
plot_names <- paste0(trial_names, "-PLOT_", T3Design %>% pull(plotNo))

# Set haveLayoutInfo to whatever location you want to load now.
allTbl <- T3Design
allTbl <- allTbl %>%
  mutate(peaName=peaAccDens[as.character(peaEntry)],
         oatDensity=if_else(is.na(oatDens),
                                  "75",
                                  c("0", "50", "100")[oatDens]),
         peaDensity=if_else(is.na(oatDens),
                            "75",
                            c("0", "50", "100")[peaDens])) %>%
  mutate(trial_name=trial_names,
         breeding_program=breeding_programs,
         location=locations,
         year=2023,
         design_type="Augmented",
         description="MultiEnvironment Oat Pea Intercropping",
         trial_type="phenotyping_trial",
         plot_width=plot_widths,
         plot_length=plot_lengths,
         field_size="",
         planting_date=planting_dates,
         harvest_date="",
         plot_name=plot_names,
         accession_name=oatAccDens[as.character(oatEntry)],
         plot_number=plotNo,
         block_number=incompBlk,
         is_a_control=isCheck,
         rep_number=1,
         range_number=rowCol[,1],
         row_number=rowCol[,1],
         col_number=rowCol[,2],
         seedlot_name="",
         num_seed_per_plot="",
         weight_gram_seed_per_plot="")


# Add in the pea genotypes
peaIncMat <- model.matrix(as.formula("~ -1 + peaName"), data=allTbl)
peaIncMat <- matrix(as.character(peaIncMat),
                    nrow=nrow(peaIncMat),
                    ncol=ncol(peaIncMat),
                    dimnames=dimnames(peaIncMat))
peaIncMat[peaIncMat == "0"] <- ""
allTbl <- allTbl %>% bind_cols(peaIncMat)
correctPeaVar <- function(peaName) return(paste0("pea_genotype_", substr(peaName, 8, nchar(peaName))))
allTbl <- allTbl %>% rename_with(.fn=correctPeaVar, .cols=contains("peaName"))

oDensIncMat <- model.matrix(as.formula("~ -1 + oatDensity"), data=allTbl)
oDensIncMat <- matrix(as.character(oDensIncMat),
                      nrow=nrow(oDensIncMat),
                      ncol=ncol(oDensIncMat),
                      dimnames=dimnames(oDensIncMat))
oDensIncMat[oDensIncMat == "0"] <- ""
allTbl <- allTbl %>% bind_cols(oDensIncMat)
correctODens <- function(oDens) return(paste0("oat_density_", substr(oDens, 11, nchar(oDens))))
allTbl <- allTbl %>% rename_with(.fn=correctODens, .cols=contains("oatDensity"))

pDensIncMat <- model.matrix(as.formula("~ -1 + peaDensity"), data=allTbl)
pDensIncMat <- matrix(as.character(pDensIncMat),
                      nrow=nrow(pDensIncMat),
                      ncol=ncol(pDensIncMat),
                      dimnames=dimnames(pDensIncMat))
pDensIncMat[pDensIncMat == "0"] <- ""
allTbl <- allTbl %>% bind_cols(pDensIncMat)
correctPDens <- function(pDens) return(paste0("pea_density_", substr(pDens, 11, nchar(pDens))))
allTbl <- allTbl %>% rename_with(.fn=correctPDens, .cols=contains("peaDensity"))

# Take out the metDesign columns that are no longer needed
allTblNoMD <- allTbl %>% dplyr::select(-(1:10)) %>%
  dplyr::relocate(location, .after=breeding_program)

readr::write_csv(allTblNoMD,
                 file=here::here("output", "oatPeaDensityT3upload.csv"))
