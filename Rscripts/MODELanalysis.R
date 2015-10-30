library(reshape2)
library(ggplot2)
library(data.table)
library(animation)

source("./Rscripts/MODELfunctions.R")
#source("./Rscripts/fournodeSUBGRAPHS.R")


# create model network

nm1 <- niche.model(S = 100, C = .25)

rel.rand(nm1, 11, 10)

# dynamic model

dyn2 <- Crmod(Adj = nm1, t = 1:500, G = G.i, method = conres, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = .2, B.o =.5, plot = FALSE)

rel.rand.z(mat = nm1, dyn = dyn2, iter = 10)

# general crmod
results <- list()
results2 <- list()
for(i in 1:200){
  test <- CRmod.gen(100, .2, "erdosrenyi", xpar = 1)
  test2 <- CRmod.gen(100, .2, "erdosrenyi", xpar = .2)
  results[[i]] <- data.frame(q = 1, iter = i, test)
  results2[[i]] <- data.frame(q = .2, iter = i, test2)
  print(i)
}
z.all <- melt(rbindlist(results), id.vars = c(1,2))
z.all2 <- melt(rbindlist(results2), id.vars = c(1,2))


ggplot(rbind(z.all, z.all2), aes(x = variable, y = value, fill = factor(q))) + geom_boxplot()

# visualize results

vis_dyn(dyn = dyn1, t.start = 100, t.end = 500)

plot.igraph(graph.adjacency(nm2), vertex.size = -1/log(dyn1[i, c(which(tail(dyn1, 1) > 0)[-1])])*100, edge.arrow.size = .5, layout = lay)

# analysis

coef.var <- apply(dyn1[,-1][,which(dyn1[200,-1] > 0)], 2, function(x){sd(x)/mean(x)})
