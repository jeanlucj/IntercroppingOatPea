library(BGLR)
data(wheat)
K <- tcrossprod(scale(wheat.X, center=TRUE))
K <- K/mean(diag(K))
Y <- wheat.Y # 4 traits

# Fitting a GBLUP un-structured cov-matrices
LP <- list(mar=list(K=K, model="RKHS"))
set.seed(123)
fmUN <- Multitrait(y=Y, ETA=LP, nIter=1000, burnIn=500,
                   saveAt="UN_",verbose=FALSE)

# For no particluar reason, just work with the first 200 accessions
k2 <- K[1:200, 1:200]
y2 <- Y[1:200,]
y2 <- rbind(y2, y2) # Pretend there are two reps.  Add noise.
y2 <- y2 + matrix(rnorm(length(y2), 0, 0.5), nrow=nrow(y2))
z2 <- rbind(diag(200), diag(200)) # Make an incidence matrix for the two reps
k2e <- z2 %*% k2 %*% t(z2) # Expand the covariance matrix to cover the two reps
LP2 <- list(mar=list(K=k2e, model="RKHS"))
set.seed(123)
fmUN <- Multitrait(y=y2, ETA=LP2, nIter=2000, burnIn=500, saveAt="UN2_", verbose=FALSE)

# Retrieving estimates and posterior SD
fmUN$resCov$R # residual covariance matrix
fmUN$resCov$SD.R

fmUN$ETA$mar$Cov$Omega # genomic covariance matrix
fmUN$ETA$mar$Cov$SD.Omega

fmUN$ETA$mar$u # predicted random effects
