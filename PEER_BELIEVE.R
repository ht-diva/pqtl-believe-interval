require(RNOmni)
install.packages("RNOmni")
library(RNOmni)
data <- read.csv("/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE/cleaned_Believe_restricted.csv")
table(duplicated(colnames(data)))
DATA_COLS <- colnames(data)[68:ncol(data)]
t_data <- apply(data[,DATA_COLS],2, RankNorm)

# scale the data
s_data <- apply(t_data, 2, scale)
sample.id <- data$SampleId
# Set the maximum number of unobserved factors to model
n_peer <- 10
maxItr <- 100

library(peer)

model=PEER()
PEER_setNmax_iterations(model, maxItr)
PEER_setPhenoMean(model, as.matrix(s_data))

PEER_setNk(model, n_peer)
PEER_setPhenoMean(model, as.matrix(s_data))
PEER_setAdd_mean(model, TRUE)
# 7| (Optional) Set expression data uncertainty.
# PEER_setPhenoVar(model, as.matrix(expression_variance))
# 8|(Optional) Set observed covariates.
# PEER_setCovariates(model, as.matrix(covariates))
# 10| Train the model, observing convergence:
PEER_update(model)
#If the model is not converged after the default 1,000 iterations, and the variance of the residuals keeps decreasing, choose a higher value of iterations, e.g., 10,000.
# PEER_setNMax_iterations(model, 10000)
# A total of 100 iterations should be sufficient to reach convergence on most data sets.
PEER_plotModel(model)
#consider only including the more relevant factors, or rerun the model with this number of factors 
# 13| Correcting for learned determinants in eQTL scans using option B: SHOULD USE 15 PEER FACTORS
# This approach can only be used with parametric models. To perform per-marker parametric test for genes 7–12, using inferred PEER factors as covariates, use:
# lod_scores=scanone(cross, method=“hk”, model=“normal”, pheno.col=7:12, addcovar=PEER_getX(model))
factors = PEER_getX(model)
dim(factors)
factors <- as.data.frame(factors)
colnames(factors) <- c('MeanFactor',paste('peer', 1:n_peer, sep=''))
factors$SampleId<- sample.id
write.csv(factors, paste0('/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE_peers', '_', n_peer, '.txt'))

residuals = PEER_getResiduals(model); colnames(residuals) <- colnames(DATA_COLS)
resid <- data.frame(SampleId = sample.id, residuals)
write.csv(resid, paste0('/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE_resid', '_', n_peer, '.txt'))

precision = PEER_getAlpha(model)
write.table(precision, paste0('/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE_precision.', n_peer, '.txt'), col.names = F, row.names = F)
weights = PEER_getW(model)
write.table(weights, paste0('/group/diangelantonio/users/alessia_mapelli/PEER_prot_BELIEVE_weights.', n_peer, '.txt'), col.names = F, row.names = F)


