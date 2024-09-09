# Return a coefficient of coancestry matrix from a three-column pedigree
# The first column has to be the row number
# Sire and dam columns refer directly to rows
# Unknown parents need to be set to 0 (zero)
pedigreeToCCmatrix <- function(threeColPed){
  nInd <- nrow(threeColPed)
  ccMat <- matrix(0, nInd, nInd)
  # the very first individual in the pedigree has to be a founder
  ccMat[1, 1] <- 0.5
  for (prog in 2:nInd){
    sire <- threeColPed[prog, 2]
    dam <- threeColPed[prog, 3]
    prog_1 <- prog - 1
    if (sire){
      sireRow <- ccMat[sire, 1:prog_1]
    } else{
      sireRow <- rep(0, prog_1)
    }
    if (dam){
      damRow <- ccMat[dam, 1:prog_1]
    } else{
      damRow <- rep(0, prog_1)
    }
    ccMat[prog, 1:prog_1] <- ccMat[1:prog_1, prog] <- (sireRow + damRow) / 2
    ccSelf <- 0.5
    if (sire > 0 & dam > 0) ccSelf <- ccSelf + ccMat[sire, dam] / 2
    ccMat[prog, prog] <- ccSelf
  }
  rownames(ccMat) <- colnames(ccMat) <- 1:nInd
  return(ccMat)
}

# Utility function for pedigreeToAmat. Since pedigreeToCCmatrix needs IDs
# to refer directly to rows, if you have an ID matrix with character names,
# this function converts to integer IDs that point to rows.
convertNamesToRows <- function(nameMat){
  nameToRow <- 1:nrow(nameMat)
  names(nameToRow) <- nameMat[,1]
  parVecToRow <- function(parVec){
    rowID <- integer(length(parVec))
    rowID[parVec != "0"] <- nameToRow[parVec[parVec != "0"]]
    return(rowID)
  }
  return(cbind(nameToRow, parVecToRow(nameMat[,2]), parVecToRow(nameMat[,3])))
}

# Calculate GRM as (W %*% W^T) / sum(2pq)
# To calculate p, the function expects dosage coding to be 0, 1, 2
calcGenomicRelationshipMatrix <- function(locusMat){
  freq <- colMeans(locusMat) / 2
  locusMat <- scale(locusMat, center=T, scale=F)
  return(tcrossprod(locusMat) / sum(2*freq*(1-freq)))
}
