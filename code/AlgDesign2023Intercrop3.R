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
if (nLoc > maxnLoc){
  # For some reason, this doesn't work with just one location
  # partLoc <- optFederov(frml=as.formula("~ oatEntry + peaEntry"),
  #                      data=fullFactOP, nTrials=plotsPerLoc*(nLoc - maxnLoc))
  partLoc <- optFederov(frml=as.formula("~ oatEntry + peaEntry"),
                        data=fullFactOP, nTrials=2*plotsPerLoc)
}
plDesign <- partLoc$design
plBlkOP <- optBlock(frml=as.formula("~ oatEntry + peaEntry"), withinData=plDesign,
                  blocksizes=rep(plotsPerLoc, 2))
for (i in 1:2){
  print(table(plBlkOP$Blocks[[i]]$oatEntry))
  print(table(plBlkOP$Blocks[[i]]$peaEntry))
}

designOP <- optFederov(frml=as.formula("~ oatEntry + peaEntry"),
                       data=fullFactOP, nTrials=plotsPerLoc*maxnLoc)
alDesign <- designOP$design
alBlkOP <- optBlock(frml=as.formula("~ oatEntry + peaEntry"), withinData=alDesign,
                    blocksizes=rep(plotsPerLoc, maxnLoc))
for (i in 1:maxnLoc){
  print(table(alBlkOP$Blocks[[i]]$oatEntry))
  print(table(alBlkOP$Blocks[[i]]$peaEntry))
}

newData <- plBlkOP$Blocks[[1]]
for (i in 1:maxnLoc){
  newData <- rbind(newData, alBlkOP$Blocks[[i]])
}
blkOP <- optBlock(frml=as.formula("~ oatEntry*peaEntry"), withinData=newData,
                  rows=1:672, blocksizes=rep(plotsPerLoc, nLoc), nRepeats=50)
# blkOP <- optBlock(withinData=newData,
#                   blocksizes=rep(plotsPerLoc, nLoc))
for (i in 1:nLoc){
  print(table(blkOP$Blocks[[i]]$oatEntry))
  print(table(blkOP$Blocks[[i]]$peaEntry))
}

  # Divide the overall experiment into the number of locations
  # Then subdivide the number of locations into incomplete blocks
  # NOTE: this is imperfect: the incomplete blocks should be coordinated
  # across locations, but they are not
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
