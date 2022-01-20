# DREAM R Script
library(Rcpp)
library(knitr)
library(tidyverse)
library(robsel)
library(qgraph)
library(igraph)
library(glasso)
library(gridExtra)
library(CVglasso)
library(Matrix) # convenient to deal with sparse matrices, has sparseMatrix() in it
setwd("~/robsel-reproducible/gene-network")
source(file="../simulation/simulation_source_functions.R")
########
### read in, clean, and construct data. Un-comment code to use corresponding files of each network
########


# In Silico Network
data.net1.orig <- read.table("Net1-In-silico/input-data/net1_expression_data.tsv", header=T, sep="\t")
data.net1.orig <- scale(data.net1.orig)
net1.TF <- read.table("Net1-In-silico/input-data/net1_transcription_factors.tsv")$V1
net1.gold.standard <- read.table("Net1-In-silico/gold-standard/DREAM5_Network1_edges_only.tsv", header=F, sep="\t")

# Ecoli Network
# data.net1.orig <- read.table("Net3-Ecoli/input-data/net3_expression_data.tsv", header=T, sep="\t")
# data.net1.orig <- scale(data.net1.orig)
# net1.TF <- read.table("Net3-Ecoli/input-data/net3_transcription_factors.tsv")$V1
# net1.gold.standard <- read.table("Net3-Ecoli/gold-standard/DREAM5_Network3_edges_only.tsv", header=F, sep="\t")


# S. cerevisiae Network
# data.net1.orig <- read.table("Net4-Scerevisiae/input-data/net4_expression_data.tsv", header=T, sep="\t")
# data.net1.orig <- scale(data.net1.orig)
# net1.TF <- read.table("Net4-Scerevisiae/input-data/net4_transcription_factors.tsv")$V1
# net1.gold.standard <- read.table("Net4-Scerevisiae/gold-standard/DREAM5_Network4_edges_only.tsv", header=F, sep="\t")


net1.gene.names <- colnames(data.net1.orig)

true.net1 <- diag(length(net1.gene.names))
colnames(true.net1) <- net1.gene.names
rownames(true.net1) <- net1.gene.names
for (i in 1:length(net1.gold.standard$V1)) {
    tf.row <- net1.gold.standard[i,]
    true.net1[tf.row$V1, tf.row$V2] <- tf.row$V3
    true.net1[tf.row$V2, tf.row$V1] <- tf.row$V3
}
data.net1.subset <- data.net1.orig[,net1.TF]
true.net1.subset <- true.net1[net1.TF,net1.TF]
data.net1.subset.cov <- cov(data.net1.subset)
net1.n <- nrow(data.net1.subset)

ptm <- proc.time()
robsel.lambda <- robsel(as.matrix(data.net1.subset), alpha = 0.9)
net1.robsel <- glasso(data.net1.subset.cov, rho=robsel.lambda, penalize.diagonal=F)
net1.robsel.graph <- net1.robsel$wi
net1.robsel.time <- proc.time() - ptm
net1.robsel.time
net1.robsel.TP <- sum(net1.robsel.graph[upper.tri(net1.robsel.graph)]!=0 & true.net1.subset[upper.tri(true.net1.subset)]!=0)
net1.robsel.TP

net1.robsel.total <- sum(net1.robsel.graph[upper.tri(net1.robsel.graph)]!=0)
net1.robsel.total

ptm <- proc.time()
net1.EBIC.graph.10 <- EBICglasso(S=data.net1.subset.cov, n=net1.n, gamma=0.5, penalize.diagonal = F, nlambda=10, lambda.min.ratio=0.05)
net1.EBIC.10.time <- proc.time() - ptm
net1.EBIC.10.time

net1.EBIC.graph.10.TP <- sum(net1.EBIC.graph.10[upper.tri(net1.robsel.graph)]!=0 & true.net1.subset[upper.tri(true.net1.subset)]!=0)
net1.EBIC.graph.10.TP

net1.EBIC.graph.10.total <- sum(net1.EBIC.graph.10[upper.tri(net1.robsel.graph)]!=0)
net1.EBIC.graph.10.total

cv.numlamb <- 10
ptm <- proc.time()
net1.cv.graph <- CVglasso(X=data.net1.subset, diagonal = F, nlam=cv.numlamb, lam.min.ratio=0.05)$Omega
net1.cv.time <- proc.time() - ptm
net1.cv.time

net1.cv.graph.TP <- sum(net1.cv.graph[upper.tri(net1.cv.graph)]!=0 & true.net1.subset[upper.tri(true.net1.subset)]!=0)
net1.cv.graph.TP

net1.cv.graph.total <- sum(net1.cv.graph[upper.tri(net1.cv.graph)]!=0)
net1.cv.graph.total

ptm <- proc.time()
net1.testing.graph <- testing_fit(as.matrix(data.net1.subset), alpha = 0.9)
net1.testing.time <- proc.time() - ptm
net1.testing.time

net1.testing.graph.TP <- sum(net1.testing.graph[upper.tri(net1.testing.graph)]!=0 & true.net1.subset[upper.tri(true.net1.subset)]!=0)
net1.testing.graph.TP

net1.testing.graph.total <- sum(net1.testing.graph[upper.tri(net1.testing.graph)]!=0)
net1.testing.graph.total