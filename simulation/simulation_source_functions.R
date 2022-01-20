library(qgraph)
#Multiple testing method for graphical model selection
testing_fit <- function(data, alpha){
    n <- dim(data)[1]
    p <- dim(data)[2]
    estimate_graph <- matrix(0, p, p)
    upper_index <- upper.tri(estimate_graph)
    sample_covariance <- (t(data) %*% data) / n
    sample_precision <- solve(sample_covariance)
    sample_partial_cor <- as.matrix(wi2net(sample_precision))
    sample_partial_cor_vec <- sample_partial_cor[upper_index]
    num_test = choose(p,2)
    p_values_vec <- rep(NA, num_test)
    for (i in 1:num_test) {
        z <- atanh(sample_partial_cor_vec[i])
        p_values_vec[i] <- 2 * (1 - pnorm(sqrt(n - p + 2 - 3) * abs(z)))
    }
    pvalues_sorted_index = order(p_values_vec)
    edge_list <- rep(0, num_test)
    m = 0
    for (i in pvalues_sorted_index) {
        if (p_values_vec[i] <= alpha/(num_test - m)){
            edge_list[i] <- 1
            m <-  m+1
        } else {
            break
        }
    }
    #browser()
    estimate_graph[upper_index] <- edge_list
    return(estimate_graph + t(estimate_graph))
}

# FWER
graph_fwer <- function(graph_list, true_omega) {
    num_graph <- length(graph_list)
    fwer <- c()
    for (i in 1:num_graph){
        fp <- (graph_list[[i]][upper.tri(true_omega)]!=0) & (true_omega[upper.tri(true_omega)]==0)
        fp <- sum(fp) != 0
        fwer <- c(fwer, fp)
    }
    fwer <- mean(fwer)
    return(fwer)
}

jaccard_index <- function(graph_list_1, graph_list_2) {
    num_graph <- length(graph_list_1)
    jaccard_index <- c()
    for (i in 1:num_graph){
        intersection <- sum((graph_list_1[[i]][upper.tri(Omega)]!=0) & (graph_list_2[[i]][upper.tri(Omega)]!=0))
        union <- sum((graph_list_1[[i]][upper.tri(Omega)]!=0) | (graph_list_2[[i]][upper.tri(Omega)]!=0))
        if (union == 0) {
            jaccard <- 1
        } else {
            jaccard <- intersection/union
        }
        jaccard_index <- c(jaccard_index, jaccard)
    }
    #browser()
    jaccard_index <- mean(jaccard_index)
    return(jaccard_index)
}

graph_rev_perf <- function(graph_list, true_graph) {
  num_sims <- length(graph_list)
  grp_ind <- upper.tri(graph_list[[1]])
  graph_pred <- lapply(graph_list, function(x) factor((x!=0)[grp_ind], levels=c("FALSE", "TRUE")))
  graph_ref <- factor((true_graph!=0)[grp_ind], levels=c("FALSE", "TRUE"))
  confMat_list <- lapply(1:num_sims, function(x) caret::confusionMatrix(graph_pred[[x]], graph_ref, positive = "TRUE"))
  TP <- c()
  FP <- c()
  True_P <- c()
  True_N <- c()
  mcc <- c()
  for (i in 1:num_sims) {
    TP <- c(TP, confMat_list[[i]]$table[2,2])
    FP <- c(FP, confMat_list[[i]]$table[2,1])
    True_P <- c(True_P, sum(confMat_list[[i]]$table[,2]))
    True_N <- c(True_N, sum(confMat_list[[i]]$table[,1]))
    mcc <- c(mcc, mcc(graph_pred[[i]], graph_ref))
  }
  #browser()
  FDR <- FP/pmax((FP + TP), 1)
  TPR <- TP/(True_P)
  FPR <- FP/(True_N)
  return(list(TP = TP, FP = FP, FDR = FDR, MCC=mean(mcc), TPR=mean(TPR), FPR=mean(FPR)))
}






