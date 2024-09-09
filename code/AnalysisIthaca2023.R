# Analysis of Ithaca location of the 2023 oat pea intercrop diversity experiment
library(tidyverse)
require(BGLR)

here::i_am("code/AnalysisIthaca2023.R")

# 7 Oct 2023 data from Ithaca downloaded from T3
opiFromT3 <- readr::read_csv(
  here::here("data", "IntercroppingOatPeaDiversityPhenotypes13Sept23.csv"),
  skip=3)
opiFromT3 <- dplyr::filter(
  opiFromT3, studyName == "Cornell_OatPeaIntercropPilot_2023_Ithaca")

# 7 Oct 2023 data from Urbana for spreadsheet sent by Milcah Kigoni Sept. 18th
opiUrbana <- openxlsx::read.xlsx(
  here::here("data", "20230918_pheno_iCropYield_Urb_23.xlsx")
)
opiUrbana <- dplyr::rename(opiUrbana, plotNumber=plotNo, blockNumber=incompBlk,
                           locationName=location, germplasmName=oatName,
                           peaAcc=peaName, replicate=rep_number,
                           rowNumber=row_number, colNumber=col_number,
                           oatYield=oat.yld.gm2, peaYield=pea.yld.gm2) %>%
  dplyr::select(locationName, germplasmName,
                peaAcc, replicate,
                blockNumber, plotNumber,
                rowNumber, colNumber,
                oatYield, peaYield)


# Make a list of the accessions
oatAcc <- opiFromT3$germplasmName %>% unique
# At the moment, these accessions are not well genotyped.  Look to pedigree
source(here::here("code", "calcRelationshipMatrices.R"))
pedFromT3 <- readr::read_delim(
  here::here("data", "IntercropOatPeaDiversity_Pedigree.txt"),
  delim="\t")
pedFromT3 <- select(pedFromT3, Accession, Female_Parent, Male_Parent)
pedFromT3[is.na(pedFromT3)] <- "0"
threeCol <- matrix(unlist(pedFromT3), nrow=nrow(pedFromT3)) %>%
  convertNamesToRows
pedRelMat <- pedigreeToCCmatrix(threeCol)
oatRelMat <- 2 * pedRelMat[threeCol[oatAcc, 1], threeCol[oatAcc, 1]]
rownames(oatRelMat) <- colnames(oatRelMat) <- oatAcc
# Hmmm.  Not much pedigree connectedness either.  Out of ~1200 pairwise
# relationships, 30 are non zero. I will use this and hope to get marker data

# Get the accessions names of the pea
# "Cornell_OatPeaIntercropPilot_2023_Ithaca_pea_genotype_Maxum"
peaAcc <- dplyr::select(opiFromT3, contains("OatPeaIntercropPilot")) %>%
  colnames %>%
  strsplit(split="genotype_")
peaAcc <- sapply(peaAcc, function(v) dplyr::last(v))
peaRelMat <- diag(length(peaAcc))
rownames(peaRelMat) <- colnames(peaRelMat) <- peaAcc

# Make sure oat and pea accessions from Ithaca and Urbana agree
oatAccIL <- unique(opiUrbana$oatName)
sum(oatAccIL %in% oatAcc) # 48 yay!
peaAccIL <- unique(opiUrbana$peaName)
sum(peaAccIL %in% peaAcc) # 12 yay!

# Function to extract the pea accession name from the column name
peaAccFromColName <- function(cn){
  return(strsplit(cn, split="genotype_") %>% sapply(function(v) return(v[2])))
}

# Set up the data on grain weight of oat and pea
grainWgt <- dplyr::filter(opiFromT3, observationLevel == "plot") %>%
  dplyr::select(locationName, germplasmName, replicate, blockNumber,
                plotNumber, rowNumber, colNumber,
                contains("weight"), -contains("223"),
                contains("OatPeaIntercropPilot")) %>%
  dplyr::arrange(germplasmName)

grainWgt[is.na(grainWgt)] <- 0
getPeaName <- function(incVec){
  return(peaAcc[which(incVec == 1)[1]])
}
grainWgt <- dplyr::mutate(grainWgt,
  peaAcc=apply(dplyr::select(grainWgt, contains("OatPeaIntercropPilot")), 1, getPeaName)) %>%
  dplyr::select(-contains("OatPeaIntercropPilot")) %>%
  dplyr::relocate(peaAcc, .after=germplasmName) %>%
  dplyr::arrange(peaAcc) %>%
  rename(oatYield=contains("000029"), peaYield=contains("000038")) %>%
  dplyr::bind_rows(opiUrbana)

grainWgt <- dplyr::mutate(grainWgt,
  blockNumberF=as.factor(paste(locationName, blockNumber)))


yTraits <- as.matrix(dplyr::select(grainWgt, contains("Yield")))
incLocations <- model.matrix(~ -1 + locationName, grainWgt)
incBlocks <- model.matrix(~ -1 + blockNumberF, grainWgt)
incOatAcc <- model.matrix(~ -1 + germplasmName, grainWgt)
incPeaAcc <- model.matrix(~ -1 + peaAcc, grainWgt)
ETA <- list(list(X=incLocations, model="FIXED"),
            list(X=incBlocks, model="BRR"),
            list(X=incOatAcc, model="BRR"),
            list(X=incPeaAcc, model="BRR"))
tst2 <- BGLR::Multitrait(yTraits, ETA, intercept=TRUE,
                  resCov=list(df0=4,S0=NULL,type="UN"),
                  R2=0.5,
                  nIter=100000, burnIn=20000,
                  thin=10, saveAt="",verbose=FALSE)

oatEff <- tst2$ETA[[3]]$beta
oatEffSD <- tst2$ETA[[3]]$SD.beta
fitOE <- lm(oatEff[,1] ~ oatEff[,2])
summary(fitOE)
anova(fitOE) # Pr(>F) = 0.24 (so NS...)
pdf(here::here("output", "OatEffectScatter.pdf"))
plot(oatEff, xlab="Oat on oat yield", ylab="Oat on pea yield",
     cex.lab=1.3, cex.axis=1.3, pch=16)
dev.off()
oatEffCov <- tst2$ETA[[3]]$Cov$Omega
cov2cor(oatEffCov)

peaEff <- tst2$ETA[[4]]$beta
plot(peaEff)
peaEffCov <- tst2$ETA[[4]]$Cov$Omega
cov2cor(peaEffCov)

residYldCov <- tst2$resCov$R
cov2cor(residYldCov)
