library(tidyverse)
library(lme4)

here::i_am("code/ExptDesign2023.R")

set.seed(12345)

locationNames <- c("NY", "IL", "SD", "ND", "WI", "AL", "MN")
nLoc <- length(locationNames)
nOat <- 48
nPea <- 12

# Get the accessions Juan is working on
juan <- read_csv(here::here("data", "Jean-Luc_2023_Urb_ICROP.csv")) %>%
  rename(germplasmName=accession_name)
juan$germplasmName <- toupper(juan$germplasmName)
juanOat <- juan %>% filter(species_name=="Avena sativa")


calcTraitBLUPs <- function(trait, data){
  frmla <- paste(trait, "~ locationName + (1 | replicate:locationName) + (1 | germplasmName)") %>% as.formula
  fitUO <- lmer(frmla, data)
  traitBLUPs <- fitUO %>% ranef %>% .$germplasmName
  colnames(traitBLUPs) <- trait
  traitBLUPs <- traitBLUPs %>% mutate(germplasmName=rownames(traitBLUPs),
                                      .before=1)
  return(traitBLUPs)
}
renameTraits <- function(v){
  whichTraits <- grep("CO_350", v)
  v[whichTraits] <- v[whichTraits] %>%
    strsplit(" - ", fixed=T) %>%
    sapply(FUN=dplyr::first)
  v[whichTraits] <- gsub(" ", "", v[whichTraits])
  return(v)
}
# Get the data from 2022 UOPN
umopn22 <- read_csv(here::here("data", "UMOPN2022ALL.csv"), skip=3)
traits <- umopn22 %>% select(contains("CO_350")) %>% colnames %>%
  strsplit(" - ", fixed=T) %>% sapply(FUN=dplyr::first)
traits <- gsub(" ", "", traits)
umopn22 <- umopn22 %>% rename_with(renameTraits)
umopn22$germplasmName <- toupper(umopn22$germplasmName)
umopn22B <- lapply(traits, calcTraitBLUPs, data=umopn22)
umopn22BLUPs <- umopn22B[[1]]
for (i in 2:length(umopn22B)){
  umopn22BLUPs <- umopn22BLUPs %>% full_join(umopn22B[[i]], by="germplasmName")
}

ueopn22 <- read_csv(here::here("data", "UEOPN2022ALL.csv"), skip=3)
traits <- ueopn22 %>% select(contains("CO_350")) %>% colnames %>%
  strsplit(" - ", fixed=T) %>% sapply(FUN=dplyr::first)
traits <- gsub(" ", "", traits)
ueopn22$germplasmName <- toupper(ueopn22$germplasmName)
ueopn22 <- ueopn22 %>% rename_with(renameTraits)
ueopn22$germplasmName <- toupper(ueopn22$germplasmName)
ueopn22B <- lapply(traits, calcTraitBLUPs, data=ueopn22)
ueopn22BLUPs <- ueopn22B[[1]]
for (i in 2:length(ueopn22B)){
  ueopn22BLUPs <- ueopn22BLUPs %>% full_join(ueopn22B[[i]], by="germplasmName")
}

uopn22 <- bind_rows(ueopn22BLUPs, umopn22BLUPs)
uopn22 <- uopn22 %>% mutate(Grainyield=scale(Grainyield) %>% as.numeric,
                              Headingdate=scale(Headingdate) %>% as.numeric,
                              Plantheight=scale(Plantheight) %>% as.numeric,
                              Testweight=scale(Testweight) %>% as.numeric)

# Find accessions Juan is working on that were NOT in 2022 UOPN
uopnAcc <- union(umopn22$germplasmName, ueopn22$germplasmName)
juanNotUOPN <- setdiff(juanOat %>% pull(germplasmName), uopnAcc)
writeLines(juanNotUOPN, here::here("output", "juanNotUOPN.txt"))

# Analyze other experiments to get these traits for Juan accessions not in UOPN
urbMon22 <- read_csv(here::here("data", "UrbMon2022OatsAll.csv"), skip=3)
urbMon22 <- urbMon22 %>% select(locationName, replicate, germplasmName,
                                dplyr::contains("CO_350"))
urbMon22 <- urbMon22 %>% dplyr::mutate(blkInLoc=paste0(locationName, replicate))
urbMon22 <- urbMon22 %>% rename_with(renameTraits)

urbMon22B <- tibble(germplasmName=(urbMon22 %>% pull(germplasmName) %>%
                                       unique %>% sort))

# Heading date analyzed separately because only recorded in Urbana
frmla <- "Headingdate ~ (1 | blkInLoc) + (1 | germplasmName)" %>%
  as.formula
# fitUM <- lmer(formula = frmla,
#               data=urbMon22 %>% filter(locationName=="Urbana, IL"))
# trtBLUP <- ranef(fitUM)$germplasmName %>% scale %>% tibble %>%
#   mutate(germplasmName=rownames(trtBLUP))
# Headingdate has no variance so take means
hd <- urbMon22 %>% group_by(germplasmName) %>%
  summarize(Headingdate=mean(Headingdate, na.rm=T)) %>%
  mutate(Headingdate=scale(Headingdate) %>% as.numeric)
urbMon22B <- urbMon22B %>% dplyr::full_join(hd, by="germplasmName")

traits <- c("Grainyield", "Plantheight", "Testweight")
for (trait in traits){
  frmla <- paste0(trait, " ~ locationName + (1 | blkInLoc) + (1 | germplasmName)") %>% as.formula
  fitUM <- lmer(formula = frmla, data=urbMon22)
  trtBLUP <- ranef(fitUM)$germplasmName %>% scale %>% as.data.frame
  trtBLUP$germplasmName=rownames(trtBLUP)
  trtBLUP <- trtBLUP %>% as_tibble
  trtBLUP <- trtBLUP %>% rename_with(.fn=function(n) return(trait), .cols=contains("Intercept"))
  urbMon22B <- urbMon22B %>% dplyr::full_join(trtBLUP, by="germplasmName")
}
urbMon22BLUPs <- urbMon22B %>% filter(germplasmName %in% juanNotUOPN) %>%
  relocate(Grainyield, .before=Headingdate)

allAcc <- bind_rows(uopn22, urbMon22BLUPs)
allAcc <- allAcc %>% dplyr::mutate(
  juanAcc=if_else(germplasmName %in% juanOat$germplasmName, 1, 0),
  uEopn22=if_else(germplasmName %in% ueopn22$germplasmName, 1, 0),
  uMopn22=if_else(germplasmName %in% umopn22$germplasmName, 1, 0)
)

allAcc$germplasmName %>% unique %>% length
allAcc <- allAcc %>% dplyr::group_by(germplasmName)
allAccS <- allAcc %>% dplyr::summarize(Grainyield=mean(Grainyield),
                                Plantheight=mean(Plantheight),
                                Headingdate=mean(Headingdate),
                                Testweight=mean(Testweight),
                                juanAcc=mean(juanAcc),
                                uEopn22=mean(uEopn22),
                                uMopn22=mean(uMopn22)
)

keepRows <- (allAccS$juanAcc == 1) %>% which
nJuanMissing <- nrow(juanOat) - length(keepRows)
frmla <- "~ Grainyield + Plantheight + Headingdate + Testweight" %>% as.formula
divDesign <- AlgDesign::optFederov(frml=frmla, data=as.data.frame(allAccS),
                                   nTrials=nOat-nJuanMissing,
                                   rows=keepRows, center=T, augment=T)
oatToTest <- divDesign$design
# Figure out which from Juan is missing
which(!(juanOat$germplasmName %in% divDesign$design$germplasmName))

missing <- juanOat %>% filter(!(juanOat$germplasmName %in%
                       divDesign$design$germplasmName))
missing <- tibble(germplasmName=missing$germplasmName,
                  GrainYield=as.numeric(NA),
                  PlantHeight=as.numeric(NA),
                  HeadingDate=as.numeric(NA),
                  TestWeight=as.numeric(NA),
                  juanAcc=1,
  uEopn22=if_else(missing$germplasmName %in% ueopn22$germplasmName, 1, 0),
  uMopn22=if_else(missing$germplasmName %in% umopn22$germplasmName, 1, 0)
)
oatToTest <- bind_rows(oatToTest, missing) %>%
  select(germplasmName, juanAcc, uEopn22, uMopn22)
oatToTest <- oatToTest %>% dplyr::arrange(desc(juanAcc),
                                          germplasmName)
oatToTest <- oatToTest %>% dplyr::mutate(entryNumber=1:48)

# Get the pea information
# 2023.2.20_pea_seed_variety_info.csv has been superceded
# peaInfo <- read_csv(here::here("data", "2023.2.20_pea_seed_variety_info.csv"))
peaInfo <- read_csv(here::here("data", "Pea Source Status - Sheet1.csv"))
toUse <- c(grep("arrived", peaInfo$status), grep("shipped", peaInfo$status))
peaToTest <- peaInfo %>% slice(toUse) %>%
  mutate(entryNumber=1:12) %>%
  rename(germplasmName=variety)

# Choose repeated checks
###################################################################
# NOTE: When we know what peas we have, swap out peaToTest %>% slice(sort(sample
toClust <- uopn22 %>% filter(germplasmName %in% oatToTest$germplasmName)
nChecks <- 4

if (!exists("clustUMOPN")) clustUMOPN <- F
if (clustUMOPN){
  toClust <- umopn22BLUPs %>%
    dplyr::select(germplasmName, Headingdate, Plantheight)
}
uopnClust <- kmeans(toClust[,-1], centers=nChecks)
calcDistToCenter <- function(accVec, clust, centerNo){
  return(sqrt(crossprod(as.numeric(accVec - clust$centers[centerNo,]))))
}
findOatAtCenter <- function(centerNo){
  return(which.min(apply(toClust[,-1], 1, calcDistToCenter,
                         clust=uopnClust, centerNo)))
}
checkOats <- sapply(1:nChecks, findOatAtCenter)
checkOats <- toClust %>% slice(checkOats) %>% pull(germplasmName)
if (clustUMOPN){
  checkOats <- umopn22BLUPs %>% dplyr::filter(germplasmName %in% checkOats)
  write_csv(checkOats, file=here::here("output", "oatAccDensity.csv"),
            col_names=T, quote="none")
  stop("Just doing density")
}
checkOats <- oatToTest %>% filter(germplasmName %in% checkOats)

checkPeas <- peaToTest %>% slice(c(3,4,5,11))

checkInter <- tibble(oatEntry=sample(checkOats %>% pull(entryNumber)),
                     peaEntry=sample(checkPeas %>% pull(entryNumber)))

source(here::here("code", "AlgDesign2023Intercrop2.R"))

###################################################################
# NOTE: When we know what peas we have, swap out peaName=paste0("Pea" ...
metDesign <- tibble()
for (loc in 1:nLoc){
#  locDesign <- as_tibble(blkOP$Blocks[[loc]])
  locDesign <- as_tibble(allLoc[[loc]])
  locDesign <- locDesign %>% mutate(location=locationNames[loc])
  locDesign <- locDesign %>% mutate(oatName=oatToTest$germplasmName[oatEntry])
  locDesign <- locDesign %>% mutate(peaName=peaToTest$germplasmName[peaEntry])
  # locDesign <- locDesign %>% mutate(peaName=paste0("Pea", peaEntry))
  metDesign <- bind_rows(metDesign, locDesign)
}
metDesign <- metDesign %>% dplyr::relocate(location, plotNo, incompBlk,
                                           isCheck, oatEntry, peaEntry)
saveRDS(list(metDesign=metDesign, oatToTest=oatToTest, peaToTest=peaToTest),
        here::here("output", "oatPeaMETdesign.rds"))
write_csv(metDesign, file=here::here("output", "oatPeaMETdesign.csv"),
          col_names=T, quote="none")
write_csv(oatToTest, file=here::here("output", "oatAccessions.csv"),
          col_names=T, quote="none")
write_csv(peaToTest, file=here::here("output", "peaAccessions.csv"),
          col_names=T, quote="none")

print(head(metDesign))
print(checkInter)

print(table(metDesign$location))
for (loc in metDesign$location %>% unique){
  metDesign %>% filter(location==loc) %>% .$oatName %>% table %>% print
  metDesign %>% filter(location==loc) %>% .$peaName %>% table %>% print
}
