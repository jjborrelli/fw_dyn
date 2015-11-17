library(reshape2)
library(ggplot2)
library(data.table)
library(animation)
library(rnetcarto)

source("./Rscripts/MODELfunctions.R")
#source("./Rscripts/fournodeSUBGRAPHS.R")


# create model network

nm1 <- niche.model(S = 30, C = .2)

rel.rand(nm1, 11, 10)

nm1 <- matrix(c(0,0,1,0), nrow = 2)
# dynamic model

dyn2 <- Crmod(Adj = nm1, t = 1:500, G = G.i, method = conres, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = 10, B.o =.5, plot = TRUE)

web_props(nm1)
web_props(nm1[which(tail(dyn2, 1)[-1] > 0),which(tail(dyn2, 1)[-1] > 0)])

nc <- netcarto(nm1)
g <- graph.adjacency(nm1)
V(g)$color[nc[[1]]$name] <- nc[[1]]$module

TrophInd(nm1)

plot(g)

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



## Parallel lists
n.mod.webs <- web_maker("niche", 100, 60, .1)
er.mod.webs <- web_maker("erg", 100, 60, .1)

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("Crmod", "get.r", "G.i", "Fij", "eventfun", "conres"))
registerDoSNOW(cl)
dynRESniche <- parLapply(cl, n.mod.webs, Crmod)
dynRESerg <- parLapply(cl, er.mod.webs, Crmod)
stopCluster(cl)

final.webs.niche <- lapply(1:length(n.mod.webs), function(x){n.mod.webs[[x]][which(tail(dynRESniche[[x]], 1)[-1] > 0),which(tail(dynRESniche[[x]], 1)[-1] > 0)]})
final.webs.erg <- lapply(1:length(er.mod.webs), function(x){er.mod.webs[[x]][which(tail(dynRESerg[[x]], 1)[-1] > 0),which(tail(dynRESerg[[x]], 1)[-1] > 0)]})

conNiche <- sapply(lapply(final.webs.niche, graph.adjacency), is.connected)
conERG <- sapply(lapply(final.webs.erg, graph.adjacency), is.connected)

niche.modular <- lapply(final.webs.niche, netcarto)
niche.modularIN <- lapply(n.mod.webs, netcarto)
erg.modular <- lapply(final.webs.erg, netcarto)
erg.modularIN <- lapply(er.mod.webs, netcarto)

df1 <- data.frame(mod = sapply(erg.modular, "[[", 2), typ = rep("erg", 100), time = rep("final", 100))
df2 <- data.frame(mod = sapply(erg.modularIN, "[[", 2), typ = rep("erg", 100), time = rep("initial", 100))
df3 <- data.frame(mod = sapply(niche.modular, "[[", 2), typ = rep("niche", 100), time = rep("final", 100))
df4 <- data.frame(mod = sapply(niche.modularIN, "[[", 2), typ = rep("niche", 100), time = rep("initial", 100))

mods <- rbind(df1, df2, df3, df4)
ggplot(mods, aes(x = typ, y = mod, fill = time)) + geom_boxplot()

n.props.fin <- t(sapply(final.webs.niche, web_props))
n.props.in <- t(sapply(n.mod.webs, web_props))

e.props.fin <- t(sapply(final.webs.erg, web_props))
e.props.in <- t(sapply(er.mod.webs, web_props))

wp1 <- cbind(e.props.fin, df1)
wp2 <- cbind(e.props.in, df2)
wp3 <- cbind(n.props.fin, df3)
wp4 <- cbind(n.props.in, df4)

wprops <- rbind(wp1, wp2, wp3, wp4)
colnames(wprops) <- c("N", "C", "Ltot", "LD", "clust", "apl", "diam", "bas", "top", "mod", "typ", "time")
wp.df <- melt(wprops)
ggplot(wp.df, aes(x = variable, y = value, fill = time)) + geom_boxplot() + facet_grid(~typ)
