# Make template to upload to T3 based on ExptDesign2023
library(tidyverse)
# You need to enter this having run ExptDesign2023.R
# The key output from there is metDesign
# To not mess that up, I'm going to copy it to T3Design here

T3Design <- metDesign

locationNames <- c("NY", "IL", "SD", "ND", "WI", "AL", "MN")
breeding_programs <- c("Cornell University",
                      "University of Illinois",
                      "South Dakota State University",
                      "North Dakota State University",
                      "University of Wisconsin",
                      "Auburn University",
                      "University of Minnesota")
trialNameBase <- "_OatPeaIntercropPilot_2023_"
trialNameBreedProg <- c("Cornell", "Illinois", "SDSU", "NDSU",
                        "Wisconsin", "Auburn", "Minnesota")
trialNameLocation <- c("Ithaca", "Urbana", "Volga", "Fargo",
                        "Madison", "Auburn", "Saint_Paul")
trial_names <- paste0(trialNameBreedProg, trialNameBase, trialNameLocation)
locations <- c("Ithaca, NY - Caldwell", "Urbana, IL", "Volga, SD",
               "Fargo, ND", "Madison, WI", "Auburn, AL", "Saint Paul, MN")
planting_dates <- c("2023-April-13", "2023-March-28", "2023-April-25",
                    "2022-April-13", "2022-April-13", "2022-April-13",
                    "2022-April-13")
plot_widths <- c(1.26, 1.33, 1.52, -1, -1, -1, -1)
plot_lengths <- c(3, 3.98, 3.66, -1, -1, -1, -1)
plot_names <- list()

trialTemplateHeader <- "trial_name	breeding_program	location	year	design_type	description	trial_type	plot_width	plot_length	field_size	planting_date	harvest_date	plot_name	accession_name	plot_number	block_number	is_a_control	rep_number	range_number	row_number	col_number	seedlot_name	num_seed_per_plot	weight_gram_seed_per_plot	is_private"
cn <- strsplit(trialTemplateHeader, "\t")

# FIELD LAYOUTS
rowCol <- list()

# New York
nyMap <- read_csv(here::here("data", "NY_Oat-Pea23map.csv"), col_names=F,
                  skip=4, n_max=5) %>%
  dplyr::select(where(is.numeric))
rowCol[[1]] <- sapply(1:100, function(p) which(nyMap==p, arr.ind=T)) %>% t
rowCol[[1]][,1] <- nrow(nyMap) - rowCol[[1]][,1] + 1
plot_names[[1]] <- paste0(trial_names[1], "-PLOT_",
                          T3Design %>% dplyr::filter(location=="NY") %>%
                            pull(plotNo))

# Illinois
ilMap <- read_csv(
  here::here("data", "2023_Urbana_M461-iCrop_Map-label-02.csv"),
  col_names=T) %>%
  dplyr::select("PASS", "RANGE", "UNIQUE", "trial", "plot") %>%
  dplyr::filter(trial=="iCrop-UIUC-23")
ilPlotToOrigPlot <- function(p){
  r <- floor(p/100)
  if (r %in% c(1,3)){
    o <- p - r*100 + ifelse(r == 3, 50, 0)
  } else{
    o <- 51 - p + r*100 + ifelse(r == 4, 50, 0)
  }
}
ilPlotToOrigPlot <- function(p){
  r <- floor(p/100)
  o <- dplyr::if_else(r %in% c(1,3),
                      p - r*100 + ifelse(r == 3, 50, 0),
                      51 - p + r*100 + ifelse(r == 4, 50, 0)
  )
  return(o)
}
ilMap <- ilMap %>% dplyr::mutate(origPlot=plot %>% ilPlotToOrigPlot) %>%
  dplyr::arrange(origPlot) %>%
  bind_cols(T3Design %>% dplyr::filter(location=="IL")) %>%
  dplyr::mutate(UNIQUE=gsub(" ", "_", UNIQUE))
ilPlot <- integer(nrow(T3Design))
ilPlot[T3Design$location=="IL"] <- as.integer(ilMap$plot)
T3Design <- T3Design %>%
  dplyr::mutate(plotNo=if_else(location=="IL", ilPlot, plotNo))

getILrowcol <- function(s){
  u <- strsplit(s, "_", fixed=T)
  return(as.numeric(u[[1]][2]))
}
ilCol <- sapply(ilMap$PASS, getILrowcol)
ilRow <- sapply(ilMap$RANGE, getILrowcol)
rowCol[[2]] <- cbind(ilRow, ilCol)
plot_names[[2]] <- ilMap$UNIQUE

# South Dakota
sdMap <- read_csv(
  here::here("data", "SDSU_OatPeaIntercropPilot_2023_Volga Map and trial info.csv"),
  col_names=F,
  skip=3, n_max=15) %>%
  dplyr::select(where(is.numeric))
rowCol[[3]] <- sapply(101:200, function(p) which(sdMap==p, arr.ind=T)) %>% t
rowCol[[3]][,1] <- nrow(sdMap) - rowCol[[3]][,1] + 1
T3Design <- T3Design %>%
  dplyr::mutate(plotNo=if_else(location=="SD", plotNo+100L, plotNo))
plot_names[[3]] <- paste0(trial_names[3], "-PLOT_",
                          T3Design %>% dplyr::filter(location=="SD") %>%
                            pull(plotNo))
# Set haveLayoutInfo to whatever location you want to load now.
haveLayoutInfo <- c("IL")
allTbl <- tibble()
for (loc in locationNames){
  if (loc %in% haveLayoutInfo){
    locTbl <- T3Design %>% filter(location == loc)
    wchLoc <- which(locationNames == loc)
    locTbl <- locTbl %>%
      mutate(trial_name=trial_names[wchLoc],
             breeding_program=breeding_programs[wchLoc],
             location=locations[wchLoc],
             year=2023,
             design_type="Augmented",
             description="MultiEnvironment Oat Pea Intercropping",
             trial_type="phenotyping_trial",
             plot_width=plot_widths[wchLoc],
             plot_length=plot_lengths[wchLoc],
             field_size="",
             planting_date=planting_dates[wchLoc],
             harvest_date="",
             plot_name=plot_names[[wchLoc]],
             accession_name=oatName,
             plot_number=plotNo,
             block_number=incompBlk,
             is_a_control=isCheck,
             rep_number=1,
             range_number=rowCol[[wchLoc]][,1],
             row_number=rowCol[[wchLoc]][,1],
             col_number=rowCol[[wchLoc]][,2],
             seedlot_name="",
             num_seed_per_plot="",
             weight_gram_seed_per_plot="")
    allTbl <- dplyr::bind_rows(allTbl, locTbl)
  }
}

# Add in the pea genotypes
peaIncMat <- model.matrix(as.formula("~ -1 + peaName"), data=allTbl)
peaIncMat <- matrix(as.character(peaIncMat),
                    nrow=nrow(peaIncMat),
                    ncol=ncol(peaIncMat),
                    dimnames=dimnames(peaIncMat))
peaIncMat[peaIncMat == "0"] <- ""
allTbl <- allTbl %>% bind_cols(peaIncMat)
correctPeaVar <- function(peaName) return(paste0("pea_genotype_",
                                                 substr(peaName, 8, nchar(peaName))))
allTbl <- allTbl %>% rename_with(.fn=correctPeaVar, .cols=contains("peaName"))

# Take out the metDesign columns that are no longer needed
allTblNoMD <- allTbl %>% dplyr::select(-(2:8)) %>%
  dplyr::relocate(location, .after=breeding_program)

readr::write_csv(allTblNoMD, file=here::here("output",
  paste0("oatPeaT3upload", paste0(haveLayoutInfo, collapse=""), ".csv")))
