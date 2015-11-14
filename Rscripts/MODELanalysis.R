library(reshape2)
library(ggplot2)
library(data.table)
library(animation)
library(rnetcarto)

source("./Rscripts/MODELfunctions.R")
#source("./Rscripts/fournodeSUBGRAPHS.R")


# create model network

nm1 <- niche.model(S = 100, C = .2)

rel.rand(nm1, 11, 10)

nm1 <- matrix(c(0,0,1,0), nrow = 2)
# dynamic model

dyn2 <- Crmod(Adj = nm1, t = 1:500, G = G.i, method = conres, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = 10, B.o =.5, plot = TRUE)

web_props(nm1)
web_props(nm1[which(tail(dyn2, 1)[-1] > 0),which(tail(dyn2, 1)[-1] > 0)])

netcarto(nm1[which(tail(dyn2, 1)[-1] > 0),which(tail(dyn2, 1)[-1] > 0)])


K = 1
x.i = .5
yij = 6
xpar = 0
B.o = .5
A =  matrix(c(0,0,1,0), nrow = 2)
FR = Fij

prey <- seq(0, 2, .01)
pred <- .5
x <- seq(0, 5, .2)
eaten <- matrix(nrow = length(prey), ncol = length(x))
for(j in 1:length(x)){
  for(i in 1:length(prey)){
    states = c(prey[i], pred)
    eaten[i,j] <- rowSums((x.i * yij * FR(states, A, B.o, xpar = x[j]) * states))[2]
  }
}
matplot(prey, eaten, typ = "l")

FR = Fbd
prey <- seq(0, 2, .01)
pred <- .5
x <- seq(0, 5, .2)
eaten <- matrix(nrow = length(prey), ncol = length(x))
for(j in 1:length(x)){
  for(i in 1:length(prey)){
    states = c(prey[i], pred)
    eaten[i,j] <- rowSums((x.i * yij * FR(states, A, B.o, xpar = x[j]) * states))[2]
  }
}
matplot(prey, eaten, typ = "l")

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
