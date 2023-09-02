# To implement:
# 1. Run ExptDesign2023.R up to row 200, so that
# AlgDesign2023Intercrop2.R
# gets run
# 2. Run AlgDesignDensity.R

# PROBLEM
# Since I do a complete factorial of different oat genotypes, including ones at
# density zero, there are a number of times when the density is zero but of
# different oat (or pea) varieties. DUMB
# SOLUTION
# 1. Create all the meaningful combinations
#   -- Full factorial of oat and pea varieties and densities when above zero
#   -- For oat density at zero, factorial of pea varieties and densities
#   -- For pea density at zero, factorial of oat varieties and densities
# 2. Pin in place the ones that Peter has already packed
# 3. Add ones to optimize
# ~ oatEntry + peaEntry + oatEntry:peaEntry + oatDens + peaDens
# WARNING: not sure how to make sure I get the same randomization each time.
here::i_am("code/AlgDesignDensity.R")

library("AlgDesign")
library("tidyverse")

nOat <- 4
nPea <- 4
oatDens <- 0:2 * 50
peaDens <- 0:2 * 50

# Design just for oat and pea then block by location
fullFactDens <- gen.factorial(c(nOat, nPea, length(oatDens), length(peaDens)),
                factors="all", center=F,
                varNames=c("oatEntry", "peaEntry", "oatDens", "peaDens"))
# Remove the entries where both oat and pea density are zero
fullFactDens <- fullFactDens %>% dplyr::filter(!(oatDens == 1 & peaDens == 1))
designDens32 <- optFederov(
  frml=as.formula("~ oatEntry + peaEntry + oatEntry:peaEntry + oatDens + peaDens"),
  data=fullFactDens, nTrials=32)
# I want to have a total of 96 plots, with 32 of them replicated, so keep the
# 32 entries from designDens32
designDens64 <- optFederov(
  frml=as.formula("~ oatEntry + peaEntry + oatEntry:peaEntry + oatDens + peaDens"),
  data=fullFactDens, nTrials=64,
  augment=T, rows=as.numeric(rownames(designDens32$design)))
# Now put it together and try to block
designDensTot <- rbind(designDens32$design, designDens64$design)
blkDens <- optBlock(
  frml=as.formula("~ oatEntry + peaEntry + oatEntry:peaEntry + oatDens + peaDens"),
  withinData=designDensTot, blocksizes=rep(24, 4))

for (blk in 1:4){
  for (i in 1:4){
    print(table(blkDens$Blocks[[blk]][,i]))
  }
}

clustUMOPN <- T
# NOTE: This will stop at "Just doing density"
source(here::here("code", "ExptDesign2023.R"))

nOat <- 4
nPea <- 4
oatDens <- 0:2 * 50
peaDens <- 0:2 * 50

# Intercalate rows of this tibble into the incomplete blocks
# checkInter comes from ExptDesign2023.R
checkInter <- checkInter %>%
  mutate(oatEntry=factor(oatEntry, levels=levels(fullFactOP$oatEntry)),
         peaEntry=factor(peaEntry, levels=levels(fullFactOP$peaEntry)),
         isCheck=1)

incBlkPerLoc <- 4
plotsPerIncBlk <- 24
densDesign <- tibble()
whichCheckWhere <- sample(incBlkPerLoc)
for (incBlk in 1:incBlkPerLoc){
  ibDesign <- as_tibble(blkDens$Blocks[[incBlk]]) %>% mutate(isCheck=0)
  ibDesign <- bind_rows(ibDesign,
                        checkInter %>% slice(whichCheckWhere[incBlk]))
  ibDesign <- ibDesign %>% dplyr::arrange(sample(plotsPerIncBlk+1))
  ibDesign <- ibDesign %>% dplyr::mutate(incompBlk=incBlk)
  densDesign <- dplyr::bind_rows(densDesign, ibDesign)
}
densDesign <- densDesign %>% dplyr::mutate(plotNo=1:nrow(densDesign))

saveWithMistake <- FALSE
if (saveWithMistake){
  saveRDS(list(densDesign, allLoc),
          here::here("output", "OatPeaDensityDesign.rds"))
  write_csv(densDesign, file=here::here("output", "OatPeaDensityDesign.csv"),
            col_names=T, quote="none")
# Check that densDesign actually has the 32 repeated combinations
checkRep <- densDesign %>%
  dplyr::mutate(combo=paste0(oatEntry, peaEntry, oatDens, peaDens))
print(sum(duplicated(checkRep$combo)))
}

#############################################################################
# Implement fix
# Design just for oat and pea then block by location
fullFactDens <- gen.factorial(c(nOat, nPea, length(oatDens), length(peaDens)),
                              factors=1:2, center=F,
                              varNames=c("oatEntry", "peaEntry", "oatDens", "peaDens"))
# Remove the entries where both oat and pea density are zero
fullFactDens <- fullFactDens %>% dplyr::filter(oatDens > 1 & peaDens > 1)
fullFactDens$oatDens <- as.factor(fullFactDens$oatDens)
fullFactDens$peaDens <- as.factor(fullFactDens$peaDens)

# Find the correct unique combos in the existing design
densDnum <- densDesign
for (i in 1:4) densDnum[,i] <- densDnum[,i] %>% unlist %>% as.numeric
densDic <- densDnum %>% dplyr::filter(oatDens > 1 & peaDens > 1)
densDmp <- densDnum %>% dplyr::filter(oatDens == 1 & peaDens > 1)
densDmo <- densDnum %>% dplyr::filter(oatDens > 1 & peaDens == 1)
# Make code to find duplicates
densDic <- densDic %>% dplyr::mutate(code=paste0(oatEntry, peaEntry, oatDens, peaDens))
densDmp <- densDmp %>% dplyr::mutate(code=paste0(peaEntry, peaDens))
densDmo <- densDmo %>% dplyr::mutate(code=paste0(oatEntry, oatDens))
# remove duplicates
densDic2 <- densDic %>% dplyr::filter(!duplicated(code))
densDmp2 <- densDmp %>% dplyr::filter(!duplicated(code))
densDmp2$oatEntry <- 1
densDmo2 <- densDmo %>% dplyr::filter(!duplicated(code))
densDmo2$peaEntry <- 1

# Bring it all together
densDgood <- dplyr::bind_rows(densDic2, densDmp2, densDmo2)
densDgood$oatEntry <- as.factor(densDgood$oatEntry)
densDgood$peaEntry <- as.factor(densDgood$peaEntry)
densDgood$oatDens <- as.factor(densDgood$oatDens)
densDgood$peaDens <- as.factor(densDgood$peaDens)

designDens32 <- optFederov(
  frml=as.formula("~ oatEntry + peaEntry + oatDens + peaDens"),
  augment=TRUE, rows=33:48,
  data=densDgood, nTrials=32)
View(designDens32$design)

densDesGood <- bind_rows(designDens32$design[,1:4], fullFactDens)

blkDens <- optBlock(
  frml=as.formula("~ oatEntry + peaEntry + oatEntry:peaEntry + oatDens + peaDens"),
  withinData=densDesGood, blocksizes=rep(24, 4))

for (blk in 1:4){
  for (i in 1:4){
    print(table(blkDens$Blocks[[blk]][,i]))
  }
}

incBlkPerLoc <- 4
plotsPerIncBlk <- 24
densDesign2 <- tibble()
whichCheckWhere <- sample(incBlkPerLoc)
for (incBlk in 1:incBlkPerLoc){
  ibDesign <- as_tibble(blkDens$Blocks[[incBlk]]) %>% mutate(isCheck=0)
  ibDesign <- bind_rows(ibDesign,
                        checkInter %>% slice(whichCheckWhere[incBlk]))
  ibDesign <- ibDesign %>% dplyr::arrange(sample(plotsPerIncBlk+1))
  ibDesign <- ibDesign %>% dplyr::mutate(incompBlk=incBlk)
  densDesign2 <- dplyr::bind_rows(densDesign2, ibDesign)
}
densDesign2 <- densDesign2 %>% dplyr::mutate(plotNo=1:nrow(densDesign2))

saveRDS(list(densDesign2, allLoc),
        here::here("output", "OatPeaDensityDesign.rds"))
write_csv(densDesign2, file=here::here("output", "OatPeaDensityDesign.csv"),
          col_names=T, quote="none")
# Check that densDesign2 actually has the 32 repeated combinations
checkRep <- densDesign2 %>%
  dplyr::mutate(combo=paste0(oatEntry, peaEntry, oatDens, peaDens))
print(sum(duplicated(checkRep$combo)))

# The pea accessions are:
# 1 AUSTRIAN FIELD PEAS
# 2 delta field pea
# 3 DS Admiral
# 4 Icicle Winter Peas
