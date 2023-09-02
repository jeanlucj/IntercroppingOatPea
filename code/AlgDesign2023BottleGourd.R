# Simple incomplete block design with two complete blocks
here::i_am("code/AlgDesign2023BottleGourd.R")

library("crossdes")
library("tidyverse")

# find.BIB uses AlgDesign::optBlock
bibNoRepChk <- find.BIB(trt=176, b=16, k=22, iter=30)
checkEntry <- 176
addRepChk <- function(bibBlk, checkEntry){
  nEntries <- length(bibBlk)
  whereChk <- sample(0:nEntries, size=1) + 1
  withChk <- integer(nEntries+1)
  withChk[whereChk] <- checkEntry
  withChk[-whereChk] <- sample(bibBlk)
  return(withChk)
}
test <- apply(bibNoRepChk, MARGIN=1, FUN=addRepChk, checkEntry=176)
