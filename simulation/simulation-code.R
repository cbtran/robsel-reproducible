library(knitr)
library(tidyverse)
library(robsel)
library(ggpubr)
library(igraph)
library(qgraph)
library(scales)
library(CVglasso)
library(mvtnorm)
library(caret)
library(latex2exp)
library(mltools)
library(gridExtra)
library(parallel)
library(Matrix) # convenient to deal with sparse matrices, has sparseMatrix() in it
setwd("~/robsel-reproducible/simulation")
source(file="simulation_source_functions.R")

#ER_graph
d=100
prob=0.02
set.seed(1)
ER_graph = sample_gnp(n=d, p=prob, directed = FALSE, loops = F)
adj_matrix_ER = as_adjacency_matrix(ER_graph, sparse=F)
p = d
#Edge weights
edge_weights =  matrix(runif(p*p, min=0.5, max=1), ncol=p)
edge_sign =  matrix(sample(c(-1,1), size=p*p, replace=T), ncol=p)
Omega = edge_weights * edge_sign
Omega = Omega * as_adjacency_matrix(ER_graph, sparse=F)
Omega_0 = Omega
#Positive Definite procedure
off_sum = rowSums(abs(Omega_0)) * 1.5
off_sum[off_sum==0] = 1
off_sum = matrix(rep(off_sum, p), ncol=p)
Omega = Omega_0/off_sum
Omega = (Omega + t(Omega))/2
set.seed(1)
diag(Omega) = runif(d, min=1, max=1.5)
Sigma = solve(Omega)
min(diag(Sigma))
max(diag(Sigma))

alphas <- seq(0.05, 0.95, 0.1)
#Generate data
n50.data <- list()
for (i in 1:200) {
    set.seed(i)
    n50.data[[i]] <- rmvnorm(n=50, mean=rep(0,p), sigma=Sigma)
}
robsel_fwer_50 <- rep(NA, 10)
robsel_tpr_50 <- rep(NA, 10)
robsel_fpr_50 <- rep(NA, 10)
robsel_mcc_50 <- rep(NA, 10)
robsel.fit.n50 <- lapply(n50.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n50[[j]]$Omega[[i]]
    }
    robsel_fwer_50[i] <- graph_fwer(robsel.fit.alpha, Omega)
    fit_res_50 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_50[i] <- fit_res_50$TPR
    robsel_fpr_50[i] <- fit_res_50$FPR
    robsel_mcc_50[i] <- fit_res_50$MCC
}
# Sample size n = 100
n100.data <- list()
for (i in 1:200) {
    set.seed(i)
    n100.data[[i]] <- rmvnorm(n=100, mean=rep(0,p), sigma=Sigma)
}
robsel_fwer_100 <- rep(NA, 10)
robsel_tpr_100 <- rep(NA, 10)
robsel_fpr_100 <- rep(NA, 10)
robsel_mcc_100 <- rep(NA, 10)
robsel.fit.n100 <- lapply(n100.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n100[[j]]$Omega[[i]]
    }
    robsel_fwer_100[i] <- graph_fwer(robsel.fit.alpha, Omega)
    fit_res_100 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_100[i] <- fit_res_100$TPR
    robsel_fpr_100[i] <- fit_res_100$FPR
    robsel_mcc_100[i] <- fit_res_100$MCC
}
# Sample size n = 200
n200.data <- list()
for (i in 1:200) {
    set.seed(i)
    n200.data[[i]] <- rmvnorm(n=200, mean=rep(0,p), sigma=Sigma)
}
testing_fwer_200 <- rep(NA, 10)
robsel_fwer_200 <- rep(NA, 10)
testing_tpr_200 <- rep(NA, 10)
robsel_tpr_200 <- rep(NA, 10)
testing_fpr_200 <- rep(NA, 10)
robsel_fpr_200 <- rep(NA, 10)
testing_mcc_200 <- rep(NA, 10)
robsel_mcc_200 <- rep(NA, 10)
jaccard_index_200 <- rep(NA, 10)
robsel.fit.n200 <- lapply(n200.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    testing <- lapply(n200.data, testing_fit, alpha=alphas[i])
    testing_fwer_200[i] <- graph_fwer(testing, Omega)
    testing_res_200 <- graph_rev_perf(testing, Omega)
    testing_tpr_200[i] <- testing_res_200$TPR
    testing_fpr_200[i] <- testing_res_200$FPR
    testing_mcc_200[i] <- testing_res_200$MCC
    
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n200[[j]]$Omega[[i]]
    }
    robsel_fwer_200[i] <- graph_fwer(robsel.fit.alpha, Omega)
    robsel_res_200 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_200[i] <- robsel_res_200$TPR
    robsel_fpr_200[i] <- robsel_res_200$FPR
    robsel_mcc_200[i] <- robsel_res_200$MCC
    
    
    
    jaccard_index_200[i] <- jaccard_index(robsel.fit.alpha, testing)
}

# Sample size n = 400
n400.data <- list()
for (i in 1:200) {
    set.seed(i)
    n400.data[[i]] <- rmvnorm(n=400, mean=rep(0,p), sigma=Sigma)
}
testing_fwer_400 <- rep(NA, 10)
testing_tpr_400 <- rep(NA, 10)
testing_fpr_400 <- rep(NA, 10)
testing_mcc_400 <- rep(NA, 10)

robsel_fwer_400 <- rep(NA, 10)
robsel_tpr_400 <- rep(NA, 10)
robsel_fpr_400 <- rep(NA, 10)
robsel_mcc_400 <- rep(NA, 10)

jaccard_index_400 <- rep(NA, 10)
robsel.fit.n400 <- lapply(n400.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    testing <- lapply(n400.data, testing_fit, alpha=alphas[i])
    testing_fwer_400[i] <- graph_fwer(testing, Omega)
    testing_res_400 <- graph_rev_perf(testing, Omega)
    testing_tpr_400[i] <- testing_res_400$TPR
    testing_fpr_400[i] <- testing_res_400$FPR
    testing_mcc_400[i] <- testing_res_400$MCC
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n400[[j]]$Omega[[i]]
    }
    robsel_fwer_400[i] <- graph_fwer(robsel.fit.alpha, Omega)
    robsel_res_400 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_400[i] <- robsel_res_400$TPR
    robsel_fpr_400[i] <- robsel_res_400$FPR
    robsel_mcc_400[i] <- robsel_res_400$MCC
    
    jaccard_index_400[i] <- jaccard_index(robsel.fit.alpha, testing)
}


# Sample size n = 800
n800.data <- list()
for (i in 1:200) {
    set.seed(i)
    n800.data[[i]] <- rmvnorm(n=800, mean=rep(0,p), sigma=Sigma)
}
testing_fwer_800 <- rep(NA, 10)
testing_tpr_800 <- rep(NA, 10)
testing_fpr_800 <- rep(NA, 10)
testing_mcc_800 <- rep(NA, 10)

robsel_fwer_800 <- rep(NA, 10)
robsel_tpr_800 <- rep(NA, 10)
robsel_fpr_800 <- rep(NA, 10)
robsel_mcc_800 <- rep(NA, 10)

jaccard_index_800 <- rep(NA, 10)
robsel.fit.n800 <- lapply(n800.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    testing <- lapply(n800.data, testing_fit, alpha=alphas[i])
    testing_fwer_800[i] <- graph_fwer(testing, Omega)
    testing_res_800 <- graph_rev_perf(testing, Omega)
    testing_tpr_800[i] <- testing_res_800$TPR
    testing_fpr_800[i] <- testing_res_800$FPR
    testing_mcc_800[i] <- testing_res_800$MCC
    
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n800[[j]]$Omega[[i]]
    }
    robsel_fwer_800[i] <- graph_fwer(robsel.fit.alpha, Omega)
    robsel_res_800 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_800[i] <- robsel_res_800$TPR
    robsel_fpr_800[i] <- robsel_res_800$FPR
    robsel_mcc_800[i] <- robsel_res_800$MCC
    
    jaccard_index_800[i] <- jaccard_index(robsel.fit.alpha, testing)
}

# Sample size n = 1600
n1600.data <- list()
for (i in 1:200) {
    set.seed(i)
    n1600.data[[i]] <- rmvnorm(n=1600, mean=rep(0,p), sigma=Sigma)
}
testing_fwer_1600 <- rep(NA, 10)
testing_tpr_1600 <- rep(NA, 10)
testing_fpr_1600 <- rep(NA, 10)
testing_mcc_1600 <- rep(NA, 10)

robsel_fwer_1600 <- rep(NA, 10)
robsel_tpr_1600 <- rep(NA, 10)
robsel_fpr_1600 <- rep(NA, 10)
robsel_mcc_1600 <- rep(NA, 10)

jaccard_index_1600 <- rep(NA, 10)
robsel.fit.n1600 <- lapply(n1600.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    testing <- lapply(n1600.data, testing_fit, alpha=alphas[i])
    testing_fwer_1600[i] <- graph_fwer(testing, Omega)
    testing_res_1600 <- graph_rev_perf(testing, Omega)
    testing_tpr_1600[i] <- testing_res_1600$TPR
    testing_fpr_1600[i] <- testing_res_1600$FPR
    testing_mcc_1600[i] <- testing_res_1600$MCC
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n1600[[j]]$Omega[[i]]
    }
    robsel_fwer_1600[i] <- graph_fwer(robsel.fit.alpha, Omega)
    robsel_res_1600 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_1600[i] <- robsel_res_1600$TPR
    robsel_fpr_1600[i] <- robsel_res_1600$FPR
    robsel_mcc_1600[i] <- robsel_res_1600$MCC
    
    jaccard_index_1600[i] <- jaccard_index(robsel.fit.alpha, testing)
}


# Sample size n = 3200
n3200.data <- list()
for (i in 1:200) {
    set.seed(i)
    n3200.data[[i]] <- rmvnorm(n=3200, mean=rep(0,p), sigma=Sigma)
}
testing_fwer_3200 <- rep(NA, 10)
testing_tpr_3200 <- rep(NA, 10)
testing_fpr_3200 <- rep(NA, 10)
testing_mcc_3200 <- rep(NA, 10)

robsel_fwer_3200 <- rep(NA, 10)
robsel_tpr_3200 <- rep(NA, 10)
robsel_fpr_3200 <- rep(NA, 10)
robsel_mcc_3200 <- rep(NA, 10)

jaccard_index_3200 <- rep(NA, 10)
robsel.fit.n3200 <- lapply(n3200.data, robsel.glasso, alpha=alphas, penalize.diagonal=F)
for (i in 1:10) {
    testing <- lapply(n3200.data, testing_fit, alpha=alphas[i])
    testing_fwer_3200[i] <- graph_fwer(testing, Omega)
    testing_res_3200 <- graph_rev_perf(testing, Omega)
    testing_tpr_3200[i] <- testing_res_3200$TPR
    testing_fpr_3200[i] <- testing_res_3200$FPR
    testing_mcc_3200[i] <- testing_res_3200$MCC
    robsel.fit.alpha <- list()
    for (j in 1:200) {
        robsel.fit.alpha[[j]] <- robsel.fit.n3200[[j]]$Omega[[i]]
    }
    robsel_fwer_3200[i] <- graph_fwer(robsel.fit.alpha, Omega)
    robsel_res_3200 <- graph_rev_perf(robsel.fit.alpha, Omega)
    robsel_tpr_3200[i] <- robsel_res_3200$TPR
    robsel_fpr_3200[i] <- robsel_res_3200$FPR
    robsel_mcc_3200[i] <- robsel_res_3200$MCC
    
    jaccard_index_3200[i] <- jaccard_index(robsel.fit.alpha, testing)
}

n50.df <- data.frame(n=50,
                     alpha=alphas,
                     FWER=robsel_fwer_50,
                     TPR=robsel_tpr_50,
                     FPR=robsel_fpr_50,
                     MCC=robsel_mcc_50,
                     method="RobSel")
n100.df <- data.frame(n=100,
                      alpha=alphas,
                      FWER=robsel_fwer_100,
                      TPR=robsel_tpr_100,
                      FPR=robsel_fpr_100,
                      MCC=robsel_mcc_100,
                      method="RobSel")

n200.df <- rbind(data.frame(n=200,
                            alpha=alphas,
                            FWER=testing_fwer_200,
                            TPR=testing_tpr_200,
                            FPR=testing_fpr_200,
                            MCC=testing_mcc_200,
                            method="Holm"),
                 data.frame(n=200,
                            alpha=alphas,
                            FWER=robsel_fwer_200,
                            TPR=robsel_tpr_200,
                            FPR=robsel_fpr_200,
                            MCC=robsel_mcc_200,
                            method="RobSel"))
n400.df <- rbind(data.frame(n=400,
                            alpha=alphas,
                            FWER=testing_fwer_400,
                            TPR=testing_tpr_400,
                            FPR=testing_fpr_400,
                            MCC=testing_mcc_400,
                            method="Holm"),
                 data.frame(n=400,
                            alpha=alphas,
                            FWER=robsel_fwer_400,
                            TPR=robsel_tpr_400,
                            FPR=robsel_fpr_400,
                            MCC=robsel_mcc_400,
                            method="RobSel"))
n800.df <- rbind(data.frame(n=800,
                            alpha=alphas,
                            FWER=testing_fwer_800,
                            TPR=testing_tpr_800,
                            FPR=testing_fpr_800,
                            MCC=testing_mcc_800,
                            method="Holm"),
                 data.frame(n=800,
                            alpha=alphas,
                            FWER=robsel_fwer_800,
                            TPR=robsel_tpr_800,
                            FPR=robsel_fpr_800,
                            MCC=robsel_mcc_800,
                            method="RobSel"))

n1600.df <- rbind(data.frame(n=1600,
                             alpha=alphas,
                             FWER=testing_fwer_1600,
                             TPR=testing_tpr_1600,
                             FPR=testing_fpr_1600,
                             MCC=testing_mcc_1600,
                             
                             method="Holm"),
                  data.frame(n=1600,
                             alpha=alphas,
                             FWER=robsel_fwer_1600,
                             TPR=robsel_tpr_1600,
                             FPR=robsel_fpr_1600,
                             MCC=robsel_mcc_1600,
                             method="RobSel"))

n3200.df <- rbind(data.frame(n=3200,
                             alpha=alphas,
                             FWER=testing_fwer_3200,
                             TPR=testing_tpr_3200,
                             FPR=testing_fpr_3200,
                             MCC=testing_mcc_3200,
                             method="Holm"),
                  data.frame(n=3200,
                             alpha=alphas,
                             FWER=robsel_fwer_3200,
                             TPR=robsel_tpr_3200,
                             FPR=robsel_fpr_3200,
                             MCC=robsel_mcc_3200,
                             method="RobSel"))

plot.df <- rbind(n50.df,n100.df,n200.df, n400.df, n800.df, n1600.df, n3200.df)
plot.df$n <- factor(plot.df$n)
plot.df$method <- factor(plot.df$method, levels=c("RobSel", "Holm")) 
colors <- hue_pal()(7)
sim.plot.TPR <- ggplot(data=plot.df, aes(x=alpha, y=TPR, color=n, shape=n, linetype=method)) + 
    geom_point(size=2.5) +
    geom_line(size=1) +
    theme_classic() +
    labs(color = "n", linetype = "method", shape = "n", x=TeX("$\\alpha$"))+
    scale_x_continuous(breaks=alphas, limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_shape_manual(values=seq(0,6)) +
    scale_color_manual(values=colors) +
    theme(axis.text.x = element_text(color = "grey20", size = 15),
          axis.text.y = element_text(color = "grey20", size = 15),  
          axis.title.x = element_text(color = "grey20", size = 15),
          axis.title.y = element_text(color = "grey20", size = 15),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          aspect.ratio=1,
          legend.position = "none",
          legend.key.size = grid::unit(2, "lines"))
sim.plot.TPR
pdf("new_TPR.pdf")
par(mar = c(0, 0, 0, 0))
sim.plot.TPR
dev.off()

sim.plot.FPR <- ggplot(data=plot.df, aes(x=alpha, y=FPR, color=n, shape=n, linetype=method)) + 
  geom_point(size=2.5) +
  geom_line(size=1) +
  theme_classic() +
  labs(color = "n", linetype = "method", shape = "n", x=TeX("$\\alpha$"))+
  scale_x_continuous(breaks=alphas, limits = c(0,1)) +
  scale_y_continuous(limits = c(0,max(plot.df$FPR))) +
  scale_shape_manual(values=seq(0,6)) +
  scale_color_manual(values=colors) +
  theme(axis.text.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(color = "grey20", size = 15),  
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15),
        legend.text.align = 0,
        legend.title = element_text(size=15),
        legend.text=element_text(size=15),
        aspect.ratio=1,
        legend.position = "none",
        legend.key.size = grid::unit(2, "lines"))
sim.plot.FPR
pdf("new_FPR.pdf")
par(mar = c(0, 0, 0, 0))
sim.plot.FPR
dev.off()

sim.plot.FWER <- ggplot(data=plot.df, aes(x=alpha, y=FWER, color=n, shape=n, linetype=method)) + 
    geom_point(size=2.5) +
    geom_line(size=1) +
    theme_classic() +
    labs(color = "n", linetype = "method", shape = "n", x=TeX("$\\alpha$"))+
    scale_x_continuous(breaks=alphas, limits = c(0,1)) +
    scale_y_continuous(breaks=alphas, limits = c(0,1)) +
    geom_abline(intercept = 0, slope = 1) +
    scale_shape_manual(values=seq(0,6)) +
    scale_color_manual(values=colors) +
    theme(axis.text.x = element_text(color = "grey20", size = 15),
          axis.text.y = element_text(color = "grey20", size = 15),  
          axis.title.x = element_text(color = "grey20", size = 15),
          axis.title.y = element_text(color = "grey20", size = 15),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          aspect.ratio=1,
          legend.position = c(0.13, 0.63),
          legend.key.size = grid::unit(2, "lines"))
sim.plot.FWER
pdf("new_FWER.pdf")
par(mar = c(0, 0, 0, 0))
sim.plot.FWER
dev.off()


sim.plot.MCC <- ggplot(data=plot.df, aes(x=alpha, y=MCC, color=n, shape=n, linetype=method)) + 
  geom_point(size=2.5) +
  geom_line(size=1) +
  theme_classic() +
  labs(color = "n", linetype = "method", shape = "n", x=TeX("$\\alpha$"))+
  scale_x_continuous(breaks=alphas, limits = c(0,1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_shape_manual(values=seq(0,6)) +
  scale_color_manual(values=colors) +
  theme(axis.text.x = element_text(color = "grey20", size = 15),
        axis.text.y = element_text(color = "grey20", size = 15),  
        axis.title.x = element_text(color = "grey20", size = 15),
        axis.title.y = element_text(color = "grey20", size = 15),
        legend.text.align = 0,
        legend.title = element_text(size=15),
        legend.text=element_text(size=15),
        aspect.ratio=1,
        legend.position = "none",
        legend.key.size = grid::unit(2, "lines"))
sim.plot.MCC
pdf("new_MCC.pdf")
par(mar = c(0, 0, 0, 0))
sim.plot.MCC
dev.off()


jaccard.index.200df <- data.frame(n=200, alpha=alphas, jaccard=jaccard_index_200)
jaccard.index.400df <- data.frame(n=400, alpha=alphas, jaccard=jaccard_index_400)
jaccard.index.800df <- data.frame(n=800, alpha=alphas, jaccard=jaccard_index_800)
jaccard.index.1600df <- data.frame(n=1600, alpha=alphas, jaccard=jaccard_index_1600)
jaccard.index.3200df <- data.frame(n=3200, alpha=alphas, jaccard=jaccard_index_3200)
jaccard.index.df <- rbind(jaccard.index.200df,jaccard.index.400df,jaccard.index.800df,
                          jaccard.index.1600df,jaccard.index.3200df)
jaccard.index.df$n <- factor(jaccard.index.df$n)
jaccard.plot <- ggplot(data=jaccard.index.df, aes(x=alpha, y=jaccard, group=n, color=n, shape=n)) +
    geom_point(size=2.5) +
    geom_line(size=1, linetype = "dotdash") +
    theme_classic() +
    labs(y="Jaccard Index", x=TeX("$\\alpha$"))+
    scale_x_continuous(breaks=alphas, limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_shape_manual(values=seq(2,6)) +
    scale_color_manual(values=c("#53B400", "#00C094", "#00B6EB", "#A58AFF", "#FB61D7")) +
    theme(axis.text.x = element_text(color = "grey20", size = 15),
          axis.text.y = element_text(color = "grey20", size = 15),  
          axis.title.x = element_text(color = "grey20", size = 15),
          axis.title.y = element_text(color = "grey20", size = 15),
          legend.text.align = 0,
          legend.title = element_text(size=15),
          legend.text=element_text(size=15),
          aspect.ratio=1,
          legend.position = "none",
          legend.direction="horizontal",
          legend.key.size = grid::unit(2, "lines"))
jaccard.plot

pdf("new_jaccard.pdf")
par(mar = c(0, 0, 0, 0))
jaccard.plot
dev.off()




# Graph plot
n800.data <- rmvnorm(n=3200, mean=rep(0,d), sigma=Sigma)
robsel.fit <- robsel.glasso(n800.data, alpha=0.05, penalize.diagonal=F)
testing.graph.fit <- testing_fit(n800.data, 0.05)
cv.graph.fit <- CVglasso::CVglasso(n800.data)$Omega
diag(cv.graph.fit ) <- 0

#True Graph
pdf("true_graph.pdf", width=12, height=12)
g <- make_empty_graph(n = d, directed = F) %>%
    add_edges(c(t(which(adj_matrix_ER!=0,arr.ind = T)))) %>%
    set_edge_attr("color", value = "black") %>% 
    set_edge_attr("curved", value = 0)
par(mar = c(0, 0, 0, 0))
plot(g, 
     layout=layout.circle, 
     vertex.color = "lightblue",
     vertex.label.cex=0.5,
     vertex.size=6)
dev.off()

robsel.graph.1 <- robsel.fit$Omega[[1]]
diag(robsel.graph.1) <- 0
robsel_true_index <- which((robsel.graph.1!=0 & adj_matrix_ER!=0),arr.ind = T)
robsel_false_index <- which((robsel.graph.1!=0 & adj_matrix_ER==0),arr.ind = T)
g.robsel.1 <- make_empty_graph(n = d, directed = F) %>%
    add_edges(c(t(robsel_false_index)), color="red", curved=0) %>% 
    add_edges(c(t(robsel_true_index)), color="black", curved=0)
pdf("robsel_graph.pdf", width=12, height=12)
par(mar = c(0, 0, 0, 0))
plot(g.robsel.1, 
     layout=layout.circle, 
     vertex.color = "lightblue",
     vertex.label.cex=0.5,
     vertex.size=6
     #,main="Robsel graph"
)
dev.off()

testing_true_index <- which((testing.graph.fit!=0 & adj_matrix_ER!=0),arr.ind = T)
testing_false_index <- which((testing.graph.fit!=0 & adj_matrix_ER==0),arr.ind = T)
g.testing.1 <- make_empty_graph(n = d, directed = F) %>%
    add_edges(c(t(testing_false_index)), color="red", curved=0) %>% 
    add_edges(c(t(testing_true_index)), color="black", curved=0)
pdf("testing_graph.pdf", width=12, height=12)
par(mar = c(0, 0, 0, 0))
plot(g.testing.1, 
     layout=layout.circle, 
     vertex.color = "lightblue",
     vertex.label.cex=0.5,
     vertex.size=6
     #,main="testing graph"
)
dev.off()


cv_true_index <- which((cv.graph.fit!=0 & adj_matrix_ER!=0),arr.ind = T)
cv_false_index <- which((cv.graph.fit!=0 & adj_matrix_ER==0),arr.ind = T)
g.cv.1 <- make_empty_graph(n = d, directed = F) %>%
    add_edges(c(t(cv_false_index)), color="red", curved=0) %>% 
    add_edges(c(t(cv_true_index)), color="black", curved=0)
pdf("cv_graph.pdf", width=12, height=12)
par(mar = c(0, 0, 0, 0))
plot(g.cv.1, 
     layout=layout.circle, 
     vertex.color = "lightblue",
     vertex.label.cex=0.5,
     vertex.size=6
     #, main="CV graph"
)
dev.off()



