library(tidyverse)
library(lme4)

here::i_am("code/ExptDesign2023.R")

set.seed(12345)

# Get the accessions Juan is working on
juan <- read_csv(here::here("data", "Jean-Luc_2023_Urb_ICROP.csv")) %>%
  rename(germplasmName=accession_name)
juan$germplasmName <- toupper(juan$germplasmName)
juanOat <- juan %>% filter(species_name=="Avena sativa")

# Get the data from 2022 UOPN
umopn22 <- read_csv(here::here("data", "UMOPN2022.csv")) %>%
  rename(germplasmName=`Accession -- UMOPN`)
umopn22$germplasmName <- toupper(umopn22$germplasmName)
umopn22$germplasmName[umopn22$germplasmName=="SD120419"] <- "WARRIOR"

umopn22 <- umopn22 %>% filter(!(germplasmName %in% c("LSD", "HSD"))) %>%
  dplyr::rename(GrainYield=`Grain yield - g/m2`,
                HeadingDate=`Heading date - Julian day`,
                PlantHeight=`Plant height - cm`,
                TestWeight=`Test weight - g/L`)

ueopn22 <- read_csv(here::here("data", "UEOPN2022.csv")) %>%
  rename(germplasmName=`Accession -- UEOPN`)
ueopn22$germplasmName <- toupper(ueopn22$germplasmName)
ueopn22 <- ueopn22 %>% filter(!(germplasmName %in% c("LSD", "HSD"))) %>%
  dplyr::rename(GrainYield=`Grain yield - g/m2`,
                HeadingDate=`Heading date - Julian day`,
                PlantHeight=`Plant height - cm`,
                TestWeight=`Test weight - g/L`)

uopn22 <- bind_rows(ueopn22, umopn22)
uopn22 <- uopn22 %>% mutate(GrainYield=scale(GrainYield) %>% as.numeric,
                              HeadingDate=scale(HeadingDate) %>% as.numeric,
                              PlantHeight=scale(PlantHeight) %>% as.numeric,
                              TestWeight=scale(TestWeight) %>% as.numeric)

# Find accessions Juan is working on that were NOT in 2022 UOPN
uopnAcc <- union(umopn22$germplasmName, ueopn22$germplasmName)
juanNotUOPN <- setdiff(juanOat %>% pull(germplasmName), uopnAcc)
writeLines(juanNotUOPN, here::here("output", "juanNotUOPN.txt"))

# Analyze other experiments to get these traits for Juan accessions not in UOPN
urbMon22 <- read_csv(here::here("data", "UrbMon2022Oats.csv"), skip=3)
urbMon22 <- urbMon22 %>% select(locationName, replicate, germplasmName,
                                dplyr::contains("CO_350"))
urbMon22 <- urbMon22 %>% dplyr::mutate(blkInLoc=paste0(locationName, replicate))
urbMon22 <- urbMon22 %>%
  dplyr::rename(GrainYield=`Grain yield - g/m2|CO_350:0000260`,
                HeadingDate=`Heading date - Julian day|CO_350:0000270`,
                PlantHeight=`Plant height - cm|CO_350:0000232`,
                TestWeight=`Test weight - g/L|CO_350:0000259`)

notUOPNblup <- tibble(germplasmName=(urbMon22 %>% pull(germplasmName) %>%
                                       unique %>% sort))
frmla <- "HeadingDate ~ (1 | blkInLoc) + (1 | germplasmName)" %>%
  as.formula
fitUM <- lmer(formula = frmla,
              data=urbMon22 %>% filter(locationName=="Urbana, IL"))
trtBLUP <- ranef(fitUM)$germplasmName %>% scale
notUOPNblup <- bind_cols(notUOPNblup, trtBLUP)
colnames(notUOPNblup)[colnames(notUOPNblup)=="(Intercept)"] <- "HeadingDate"

traits <- c("GrainYield", "PlantHeight", "TestWeight")
for (trait in traits){
  frmla <- paste0(trait, " ~ locationName + (1 | blkInLoc) + (1 | germplasmName)") %>% as.formula
  fitUM <- lmer(formula = frmla, data=urbMon22)
  trtBLUP <- ranef(fitUM)$germplasmName %>% scale
  notUOPNblup <- bind_cols(notUOPNblup, trtBLUP)
  colnames(notUOPNblup)[colnames(notUOPNblup)=="(Intercept)"] <- trait
}

allAcc <- bind_rows(uopn22, notUOPNblup)
allAcc <- allAcc %>% dplyr::mutate(juanAcc=if_else(germplasmName %in% juanOat$germplasmName, 1, 0))

allAcc$germplasmName %>% unique %>% length
allAcc <- allAcc %>% dplyr::group_by(germplasmName)
allAccS <- allAcc %>% dplyr::summarize(GrainYield=mean(GrainYield),
                                PlantHeight=mean(PlantHeight),
                                HeadingDate=mean(HeadingDate),
                                TestWeight=mean(TestWeight),
                                juanAcc=mean(juanAcc))

keepRows <- (allAccS$juanAcc == 1) %>% which
frmla <- "~ GrainYield + PlantHeight + HeadingDate + TestWeight" %>% as.formula
divDesign <- AlgDesign::optFederov(frml=frmla, data=allAccS, nTrials=47,
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
                  juanAcc=1)
oatToTest <- bind_rows(oatToTest, missing) %>% select(germplasmName, juanAcc)
oatToTest <- oatToTest %>% dplyr::arrange(desc(juanAcc), germplasmName)
oatToTest <- oatToTest %>% dplyr::mutate(entryNumber=1:48)

# Get the pea information
peaInfo <- read_csv(here::here("data", "2023.2.20_pea_seed_variety_info.csv"))
peaToTest <- peaInfo %>% filter(use=="y") %>%
  dplyr::select(Variety, company, use_rationale, color) %>%
  mutate(entryNumber=1:16) %>%
  rename(germplasmName=Variety)

# Choose repeated checks
###################################################################
# NOTE: When we know what peas we have, swap out peaToTest %>% slice(sort(sample
toClust <- uopn22 %>% filter(germplasmName %in% oatToTest$germplasmName)
nChecks <- 4
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
checkOats <- oatToTest %>% filter(germplasmName %in% checkOats)

checkPeas <- peaToTest %>% slice(sort(sample(nPea, nChecks)))

checkInter <- tibble(oatEntry=sample(checkOats %>% pull(entryNumber)),
                     peaEntry=sample(checkPeas %>% pull(entryNumber)))

source(here::here("code", "AlgDesign2023Intercrop.R"))

###################################################################
# NOTE: When we know what peas we have, swap out peaName=paste0("Pea" ...
locationNames <- c("NY", "IL", "SD", "ND", "WI")
metDesign <- tibble()
for (loc in 1:nLoc){
  locDesign <- as_tibble(blkOP$Blocks[[loc]])
  locDesign <- locDesign %>% mutate(location=locationNames[loc])
  locDesign <- locDesign %>% mutate(oatName=oatToTest$germplasmName[oatEntry])
  #  locDesign <- locDesign %>% mutate(peaName=peaToTest$germplasmName[peaEntry])
  locDesign <- locDesign %>% mutate(peaName=paste0("Pea", peaEntry))
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
