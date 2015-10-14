library(reshape2)
library(ggplot2)


source("./Rscripts/MODELfunctions.R")
#source("./Rscripts/fournodeSUBGRAPHS.R")


# create model network

nm1 <- niche.model(S = 100, C = .25)

# dynamic model

dyn1 <- Crmod(Adj = nm1, t = 1:200, G = G.i, method = conres, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = .2, B.o =.5, plot = FALSE)

matplot(dyn1[,-1], typ = "l")
dyn1[200,-1][which(dyn1[200,-1] > 0)]

ggplot(melt(dyn1[,-1]), aes(x = Var1, y = value, col = factor(Var2))) + geom_line()