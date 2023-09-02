args <- commandArgs(trailingOnly = TRUE)
iter <- 1
if (length(args) > 0) iter <- args[1]

here::i_am("code/AlgDesign2023Intercrop.R")

library("AlgDesign")
library("tidyverse")

plotsPerLoc <- 96
plotsPerIncBlk <- 24
incBlkPerLoc <- plotsPerLoc / plotsPerIncBlk

# Design just for oat and pea then block by location
fullFactOP <- gen.factorial(c(nOat, nPea), factors="all",
                            varNames=c("oatEntry", "peaEntry"))

# Intercalate rows of this tibble into the incomplete blocks
# checkInter comes from ExptDesign2023.R
checkInter <- checkInter %>%
  mutate(oatEntry=factor(oatEntry, levels=levels(fullFactOP$oatEntry)),
         peaEntry=factor(peaEntry, levels=levels(fullFactOP$peaEntry)),
         isCheck=1)

# Make sure you have enough possibilities to not get this error:
# "nTrials must not be greater than the number of rows in data"
allLoc <- list()
maxnLoc <- nrow(fullFactOP) %/% plotsPerLoc
for (i in 0:(nLoc %/% maxnLoc)){
  designOP <- optFederov(frml=as.formula("~ oatEntry + peaEntry"),
                         data=fullFactOP, nTrials=plotsPerLoc*maxnLoc)

  design <- designOP$design

  # Divide the overall experiment into the number of locations
  # Then subdivide the number of locations into incomplete blocks
  # NOTE: this is imperfect: the incomplete blocks should be coordinated
  # across locations, but they are not
  blkOP <- optBlock(frml=as.formula("~ oatEntry + peaEntry"), withinData=design,
                    blocksizes=rep(plotsPerLoc, maxnLoc))
  for (blk in 1:maxnLoc){
    blkOP$Blocks[[blk]] <- blkOP$Blocks[[blk]] %>% arrange(sample(plotsPerLoc))
    withinLocBlk <- optBlock(frml=as.formula("~ oatEntry + peaEntry"),
                             withinData=blkOP$Blocks[[blk]],
                             blocksizes=rep(plotsPerIncBlk, incBlkPerLoc))
    locDesign <- tibble()
    whichCheckWhere <- sample(incBlkPerLoc)
    for (incBlk in 1:incBlkPerLoc){
      ibDesign <- as_tibble(withinLocBlk$Blocks[[incBlk]]) %>%
        mutate(isCheck=0)
      ibDesign <- bind_rows(ibDesign,
                            checkInter %>% slice(whichCheckWhere[incBlk]))
      ibDesign <- ibDesign %>% dplyr::arrange(sample(plotsPerIncBlk+1))
      ibDesign <- ibDesign %>% dplyr::mutate(incompBlk=incBlk)
      locDesign <- dplyr::bind_rows(locDesign, ibDesign)
    }
    locDesign <- locDesign %>% dplyr::mutate(plotNo=1:nrow(locDesign))
    print(table(locDesign$oatEntry))
    print(table(locDesign$peaEntry))
    blkOP$Blocks[[blk]] <- locDesign
  }
  allLoc <- c(allLoc, blkOP$Blocks)
}

saveRDS(list(designOP, allLoc),
        here::here("output", "OatPeaDesignByLocation.rds"))
