library(reshape2)
library(ggplot2)
library(data.table)
library(animation)
library(rnetcarto)
library(igraph)
library(NetIndices)
library(rend)

# source("./Rscripts/MODELfunctions.R")

niche_model <- function(S,C){
  cond <- FALSE
  while(!cond){
    new.mat<-matrix(0,nrow=S,ncol=S)
    ci<-vector()
    niche<-runif(S,0,1)
    r<-rbeta(S,1,((1/(2*C))-1))*niche
    
    for(i in 1:S){
      ci[i]<-runif(1,r[i]/2,niche[i])
    }
    
    r[which(niche==min(niche))]<-.00000001
    
    for(i in 1:S){
      
      for(j in 1:S){
        if(niche[j]>(ci[i]-(.5*r[i])) && niche[j]<(ci[i]+.5*r[i])){
          new.mat[j,i]<-1
        }
      }
    }
    cond <- is.connected(graph.adjacency(new.mat))
  }
  
  new.mat<-new.mat[order(apply(new.mat,2,sum)),order(apply(new.mat,2,sum))]
  
  return(new.mat)
}

strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

RESULTS <- foreach(i = 1:10, .packages = c("rend", "igraph")) %dopar% {
  nm05 <- niche_model(60, .05)
  nm15 <- niche_model(60, .15)
  nm25 <- niche_model(60, .25)
  
  sim05 <- CRsimulator(nm.05)
  sim15 <- CRsimulator(nm.15)
  sim25 <- CRsimulator(nm.25)
  
  wi05 <- WEBind(sim05, nm.05)
  wi15 <- WEBind(sim15, nm.15)
  wi25 <- WEBind(sim25, nm.25)
  
  tc05 <- trophicChange(sim05, nm.05)
  tc15 <- trophicChange(sim15, nm.15)
  tc25 <- trophicChange(sim25, nm.25)
  
  mo05 <- motifCounter3(sim05, nm.05)
  mo15 <- motifCounter3(sim15, nm.15)
  mo25 <- motifCounter3(sim25, nm.25)
  
  return(list(models = list(nm05, nm15, nm25), sims = list(sim05, sim15, sim25), wind = list(wi05, wi15, wi25), 
              troph = list(tc05, tc15, tc25), motif = list(mo05, mo15, mo25)))
}

stopCluster(cl)
end <- Sys.time()
strt - end


rand.niche <- function(n, wi){
  nlist <- lapply(1:nrow(wi), function(x){nm <- list();for(i in 1:n){nm[[i]] <- niche_model(wi[x,"N"], wi[x,"C"])};return(nm)})
  return(nlist)
}


mysim <- function(iter, S, C){
  strt <- Sys.time()
  cl <- makeCluster(detectCores() - 1)
  registerDoSNOW(cl)
  
  RESULTS <- foreach(i = 1:iter, .packages = c("rend", "igraph")) %dopar% {
    
    nm <- niche_model(60, .05)
    sim <- CRsimulator(nm)
    wi <- WEBind(sim, nm)
    tc <- trophicChange(sim, nm)
    mo <- motifCounter3(sim, nm)
    
    rnm <- rand.niche(10, wi)

    return(list(nm, sim, wi, tc, mo))
  }
  
  stopCluster(cl)
  end <- Sys.time()
  strt - end
  
  
  
}





nm.05 <- lapply(1:10, function(x){niche_model(60, .05)})
nm.15 <- lapply(1:10, function(x){niche_model(60, .15)})
nm.25 <- lapply(1:10, function(x){niche_model(60, .25)})

library(parallel)
library(doSNOW)

strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

RESULTS <- foreach(i = 1:10, .packages = "rend") %dopar% {
  sim05 <- CRsimulator(nm.05[[i]])
  sim15 <- CRsimulator(nm.15[[i]])
  sim25 <- CRsimulator(nm.25[[i]])
  
  return(list(sim05 = sim05, sim15 = sim15, sim25 = sim25))
}

stopCluster(cl)
end <- Sys.time()
strt - end


s05 <- lapply(RESULTS, "[[", 1)
s15 <- lapply(RESULTS, "[[", 2)
s25 <- lapply(RESULTS, "[[", 3)

wi05 <- lapply(1:10, function(x){WEBind(s05[[x]], nm.05[[x]])})

strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

RESULTS <- foreach(i = 1:10, .packages = c("rend", "igraph")) %dopar% {
  wi05 <- WEBind(s05[[i]], nm.05[[i]])
  wi15 <- WEBind(s15[[i]], nm.15[[i]])
  wi25 <- WEBind(s25[[i]], nm.25[[i]])
  
  return(list(wi05 = wi05, wi15 = wi15, wi25 = wi25))
}

stopCluster(cl)
end <- Sys.time()
strt - end


#############################################################################
## Parallel lists
n.mod.webs <- web_maker("niche", 100, 60, .1)
er.mod.webs <- web_maker("erg", 1000, 60, .1)

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("Crmod", "get.r", "G.i", "Fij", "eventfun", "conres"))
registerDoSNOW(cl)

# niche model webs with type II dynamics
dynRESniche1 <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = 0)
dynRESniche2 <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = .2)
dynRESniche3 <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = 1)
dynRESniche4 <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = 10)

# niche model webs with interference dynamics
dynRESniche1.bd <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = 0, FuncRes = Fbd)
dynRESniche2.bd <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = .2, FuncRes = Fbd)
dynRESniche3.bd <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = 1, FuncRes = Fbd)
dynRESniche4.bd <- parLapply(cl, n.mod.webs, Crmod, t = 1:500, xpar = 10, FuncRes = Fbd)


# random model webs with type II dynamics
dynRESerg1 <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = 0, FuncRes = Fij)
dynRESerg2 <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = .2, FuncRes = Fij)
dynRESerg3 <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = 1, FuncRes = Fij)
dynRESerg4 <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = 10, FuncRes = Fij)

# random model webs with interference dynamics
dynRESerg1.bd <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = 0, FuncRes = Fbd)
dynRESerg2.bd <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = .2, FuncRes = Fbd)
dynRESerg3.bd <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = 1, FuncRes = Fbd)
dynRESerg4.bd <- parLapply(cl, er.mod.webs, Crmod, t = 1:500, xpar = 10, FuncRes = Fbd)

stopCluster(cl)

###################################
###################################
test <- dynRESniche2[[1]]
mc <- list()
for(i in 1:nrow(test)){
  dfin <- test[i, -1]
  nfin <- init[dfin > 0, dfin > 0]
  
  mc[[i]] <- motif_counter(list(graph.adjacency(nfin)))
}



motif_loss <- function(dmat, init){
  nfin <- lapply(1:nrow(dmat), function(x){init[dmat[x, -1] > 0, dmat[x, -1] > 0]})
  mlc <- motif_counter(lapply(nfin, graph.adjacency))
  return(mlc)
}
library(compiler)
motif_loss.cmp <- cmpfun(motif_loss)

rand_motif_loss <- function(d, comp, init, runs){
  
  tmp <- lapply(1:runs, function(x) d[,c(1,sample(2:61))])
  ml <- lapply(tmp, function(x) motif_loss.cmp(x, init))
  
  tmp1 <- t(sapply(ml, function(x){unlist(apply(x, 2, function(x) unlist(which.min(x))))}))
  return((comp - colMeans(tmp1))/apply(tmp1, 2, sd))
}

test <- lapply(1:100, function(x) motif_loss.cmp(dynRESniche2a[[x]], n.mod.webs[[x]]))
test3 <- t(sapply(test, function(x){unlist(apply(x, 2, function(x) unlist(which.min(x))))}))
boxplot(test3)

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("dynRESniche2a", "test3", "n.mod.webs", "rand_motif_loss", "motif_loss.cmp", "motif_loss", "motif_counter", "rbindlist"))
registerDoSNOW(cl)
p.rml <- parLapply(cl, 1:100, function(x){rand_motif_loss(dynRESniche2a[[x]], comp = test3[x,], init = n.mod.webs[[x]], runs = 200)})
stopCluster(cl)




temp <- lapply(dynRESerg2, motif_loss)
temp3 <- t(sapply(temp, function(x){unlist(apply(x, 2, function(x) unlist(which.min(x))))}))
boxplot(temp3)

mc1 <- log(rbindlist(mc))
mc1[mc1 == -Inf] <- 0
matplot(mc1, typ = "l")

dynprops <- function(dyn, init){
  dfin <- dyn[500, -1]
  nfin <- init[dfin > 0, dfin > 0]
  
  mc <- motif_counter(list(graph.adjacency(nfin)))
  
  if(sum(dfin > 0) < 10){
    return(data.frame(N = NA, L = NA, tl = NA, mod = NA, mc))
  }
  
  if(!is.connected(graph.adjacency(nfin))){
    return(data.frame(N = NA, L = NA, tl = NA, mod = NA, mc))
  }
  
  N <- sum(dfin > 0)
  L <- sum(nfin)
  
  ti <- average.path.length(graph.adjacency(nfin))
  
  nc <- netcarto(conversion(nfin))
  
  return(data.frame(N = N, L = L, tl = ti, mod = nc[[2]], mc))
} 

initp <- lapply(n.mod.webs, function(x){
  N <- nrow(x)
  L <- sum(x)
  ti <- average.path.length(graph.adjacency(x))
  nc <- netcarto(conversion(x))
  return(data.frame(N = N, L = L, tl = ti, mod = nc[[2]]))
})

r.initp <- lapply(er.mod.webs, function(x){
  N <- nrow(x)
  L <- sum(x)
  ti <- average.path.length(graph.adjacency(x))
  nc <- netcarto(conversion(x))
  return(data.frame(N = N, L = L, tl = ti, mod = nc[[2]]))
})



initp.r <- cbind(rbindlist(initp), typ = "init", model = "niche")
r.initp.r <- cbind(rbindlist(r.initp), typ = "init", model = "er")

drn1p <- lapply(1:1000,function(x) dynprops(dynRESniche1[[x]], n.mod.webs[[x]]))
drn1p.r <- cbind(rbindlist(drn1p), typ = "drn1", model = "niche", fr = "LV", par = 0)


drn2p <- lapply(1:1000,function(x) dynprops(dynRESniche2[[x]], n.mod.webs[[x]]))
drn2p.r <- cbind(rbindlist(drn2p), typ = "drn2", model = "niche", fr = "LV", par = .2)

drn3p <- lapply(1:1000,function(x) dynprops(dynRESniche3[[x]], n.mod.webs[[x]]))
drn3p.r <- cbind(rbindlist(drn3p), typ = "drn3", model = "niche", fr = "LV", par = 1)

drn4p <- lapply(1:1000,function(x) dynprops(dynRESniche4[[x]], n.mod.webs[[x]]))
drn4p.r <- cbind(rbindlist(drn4p), typ = "drn4", model = "niche", fr = "LV", par = 10)


drnbd1p <- lapply(1:1000,function(x) dynprops(dynRESniche1.bd[[x]], n.mod.webs[[x]]))
drnbd1p.r <- cbind(rbindlist(drnbd1p), typ = "drnbd1", model = "niche", fr = "Interference", par = 0)

drnbd2p <- lapply(1:1000,function(x) dynprops(dynRESniche2.bd[[x]], n.mod.webs[[x]]))
drnbd2p.r <- cbind(rbindlist(drnbd2p), typ = "drnbd2", model = "niche", fr = "Interference", par = 0)

drnbd3p <- lapply(1:1000,function(x) dynprops(dynRESniche3.bd[[x]], n.mod.webs[[x]]))
drnbd3p.r <- cbind(rbindlist(drnbd3p), typ = "drnbd3", model = "niche", fr = "Interference", par = 0)

drnbd4p <- lapply(1:1000,function(x) dynprops(dynRESniche4.bd[[x]], n.mod.webs[[x]]))
drnbd4p.r <- cbind(rbindlist(drnbd4p), typ = "drnbd4", model = "niche", fr = "Interference", par = 0)

#######################

dre1p <- lapply(1:1000,function(x) dynprops(dynRESerg1[[x]], n.mod.webs[[x]]))
dre1p.r <- cbind(rbindlist(dre1p), typ = "dre1", model = "erg", fr = "LV", par = 0)

dre2p <- lapply(1:1000,function(x) dynprops(dynRESerg2[[x]], n.mod.webs[[x]]))
dre2p.r <- cbind(rbindlist(dre2p), typ = "dre2", model = "erg", fr = "LV", par = .2)

dre3p <- lapply(1:1000,function(x) dynprops(dynRESerg3[[x]], n.mod.webs[[x]]))
dre3p.r <- cbind(rbindlist(dre3p), typ = "dre3", model = "erg", fr = "LV", par = 1)

dre4p <- lapply(1:1000,function(x) dynprops(dynRESerg4[[x]], n.mod.webs[[x]]))
dre4p.r <- cbind(rbindlist(dre4p), typ = "dre4", model = "erg", fr = "LV", par = 10)


drebd1p <- lapply(1:1000,function(x) dynprops(dynRESerg1.bd[[x]], n.mod.webs[[x]]))
drebd1p.r <- cbind(rbindlist(drebd1p), typ = "drebd1", model = "erg", fr = "Interference", par = 0)

drebd2p <- lapply(1:1000,function(x) dynprops(dynRESerg2.bd[[x]], n.mod.webs[[x]]))
drebd2p.r <- cbind(rbindlist(drebd2p), typ = "drebd2", model = "erg", fr = "Interference", par = .2)

drebd3p <- lapply(1:1000,function(x) dynprops(dynRESerg3.bd[[x]], n.mod.webs[[x]]))
drebd3p.r <- cbind(rbindlist(drebd3p), typ = "drebd3", model = "erg", fr = "Interference", par = 1)

drebd4p <- lapply(1:1000,function(x) dynprops(dynRESerg4.bd[[x]], n.mod.webs[[x]]))
drebd4p.r <- cbind(rbindlist(drebd4p), typ = "drebd4", model = "erg", fr = "Interference", par = 10)


allDAT <- rbind(drn1p.r, drn2p.r, drn3p.r, drn4p.r,
                drnbd1p.r, drnbd2p.r, drnbd3p.r, drnbd4p.r,
                dre1p.r, dre2p.r, dre3p.r, dre4p.r,
                drebd1p.r, drebd2p.r, drebd3p.r, drebd4p.r)

allDAT <- allDAT[which(!is.na(allDAT$N)),]


ggplot(allDAT) + facet_grid(fr~.) + geom_boxplot(aes(x = factor(par), y = N, fill = model))
###################################
###################################





final.webs.niche <- lapply(1:length(n.mod.webs), function(x){n.mod.webs[[x]][which(tail(dynRESniche[[x]], 1)[-1] > 0),which(tail(dynRESniche[[x]], 1)[-1] > 0)]})
final.webs.erg <- lapply(1:length(er.mod.webs), function(x){er.mod.webs[[x]][which(tail(dynRESerg[[x]], 1)[-1] > 0),which(tail(dynRESerg[[x]], 1)[-1] > 0)]})

conNiche <- sapply(lapply(final.webs.niche, graph.adjacency), is.connected)
conERG <- sapply(lapply(final.webs.erg, graph.adjacency), is.connected)

fwn.con <- final.webs.niche[conNiche]
nmw.vm <- lapply(n.mod.webs[conNiche], vert.mot)

ex <- lapply(which(conNiche), function(x){which(tail(dynRESniche[[x]], 1)[-1] == 0)})
ex2 <- lapply(which(conNiche), function(x){
  l1 <- length(which(tail(dynRESniche[[x]], 1)[-1] == 0))
  sample(1:60, l1)
  })

ex.mot <- lapply(1:63, function(x){nmw.vm[[x]][ex[[x]],]})
ex.mot2 <- lapply(1:63, function(x){nmw.vm[[x]][ex2[[x]],]})

diffs <- list()
for(i in 1:63){
  diffs[[i]] <- ex.mot[[i]] - ex.mot2[[i]]
}

#############################################################
niche.modular <- lapply(final.webs.niche, netcarto)
niche.modularIN <- lapply(n.mod.webs, netcarto)
erg.modular <- lapply(final.webs.erg, netcarto)
erg.modularIN <- lapply(er.mod.webs, netcarto)


#####

nm.mod <- lapply(niche.modularIN, "[[", 1)[conNiche]
nmex.m <- lapply(1:63, function(x){nm.mod[[x]][ex[[x]],]})


c1 <- lapply(fwn.con, function(x){sum(x)/(nrow(x)*(nrow(x)-1))})
n1 <- lapply(fwn.con, nrow)

rgs <- list()
for(i in 1:length(c1)){
  rgs[[i]] <- erdos.renyi.game(n1[[i]], c1[[i]], "gnp", directed = T)
}

mod.rgs <- lapply(lapply(rgs, get.adjacency, sparse = F), netcarto)
#####

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
# need to append final to initial
wp.df <- melt(wprops, id.var)
ggplot(wp.df, aes(x = variable, y = value, fill = time)) + geom_boxplot() + facet_grid(~typ)




##############################################
##############################################
CRcurve <- function(iter){
  nm <- niche.model(100, .1)
  N <- c()
  mc <- matrix(nrow = iter, ncol = 13)
  for(i in 1:iter){
    dyn <- Crmod(nm)
    N[i] <- sum(dyn[200, -1] > 0)
    mc[i,]<- unlist(motif_counter(list(graph.adjacency(nm))))
    nm <- curve_ball(nm)
    print(i)
  }
  return(list(N, mc))
}

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c("Crmod", "get.r", "G.i", "Fij", "eventfun", "conres", "CRcurve", "niche.model", "motif_counter", "curve_ball"))
registerDoSNOW(cl)

crc <- parLapply(cl, 1:500, function(x){CRcurve(1000)})

stopCluster(cl)


assemble <- function(S, C, iter, new.ints = 5){
  n1 <- niche.model(S, C)
  d1 <- Crmod(n1, t = 1:500)

  n2 <- n1[d1[500, -1] > 0, d1[500, -1] > 0]
  n.int <- nrow(n2)*2+1
  
  mc.init <- motif_counter(list(graph.adjacency(n2)))
  
  mc.diff <- matrix(nrow = iter, ncol = 13)
  dyn <- matrix(nrow = iter, ncol = nrow(n2)+1)
  init.b <- matrix(nrow = iter, ncol = nrow(n2)+1)
  new.d <- matrix(nrow = 500, ncol = iter)
  tl <- matrix(nrow = iter, ncol = 4)
  ints <- c(rep(1, new.ints), rep(0, n.int - new.ints))
  s.ints <- replicate(iter, sample(ints))
  for(i in 1:iter){
    n3 <- rbind(cbind(n2, s.ints[1:nrow(n2),i]), s.ints[(nrow(n2)+1):length(ints),i])
    mc.diff[i,] <- unlist(motif_counter(list(graph.adjacency(n3))) - mc.init)
    cm1 <- Crmod(n3, t = 1:500)
    init.b[i,] <- cm1[1, -1] 
    dyn[i,] <- cm1[500,-1] 
    new.d[,i] <- cm1[,nrow(n3)]
    t1 <- TrophInd(n3)
    t2 <- TrophInd(n3[dyn[i,] > 0, dyn[i,] > 0])$TL
    if(cm1[500,nrow(n3)] != 0){t3 <- tail(t2, 1)}else{t3 = 0}
    tl[i,] <- c(mean(t1$TL), t1[nrow(n3),1], mean(t2), t3)
    print(i)
  }
  return(list(mc.diff = mc.diff, dyn = dyn, init.b = init.b, new.d = new.d, tl = tl))
}
extinct <- apply(dyn, 1, function(x) sum(x == 0))  
cor.test(apply(new.d, 2, sd)/colMeans(new.d),rowSums(mc.diff[,c(1,2,4,5)])) 






###############################
###############################

test <- Fij(dynRESniche2[[1]][1,-1], n.mod.webs[[1]], .5, .2)
test[60, 22]
test2 <- melt(test)
test3 <- test2[test2$value != 0,]
un.t <- unique(apply(test2[test2$value != 0,1:2], 1, sort))
qtest <- quantile(test3$value)
el <- list()
for(i in 1:4){
  el[[i]] <- test3[test3$value >= qtest[i], 1:2]
}

qweb <- function(dyn, web, time, xpar = .2){
  test <- Fij(dyn[time,-1], web, .5, xpar)
  test2 <- melt(test)
  test3 <- test2[test2$value != 0,]
  
  qtest <- quantile(test3$value)
  el <- list()
  for(i in 1:4){
    el[[i]] <- test3[test3$value >= qtest[i], 1:2]
  }
  return(motif_counter(lapply(el, function(x) graph.edgelist(as.matrix(x)))))
}


qtest <- lapply(1:1000, function(x){qweb(dynRESniche2[[x]], n.mod.webs[[x]], time = 500, xpar = .2)})
res <- list()
for(i in 1:4){
  res[[i]] <- cbind(melt(do.call(rbind, lapply(qtest, function(x) x[i,]))), i, t="final")
}

qtest2 <- lapply(1:1000, function(x){qweb(dynRESniche2[[x]], n.mod.webs[[x]], time = 1, xpar = .2)})
res2 <- list()
for(i in 1:4){
  res2[[i]] <- cbind(melt(do.call(rbind, lapply(qtest2, function(x) x[i,]))), i, t="init")
}

ggplot(rbindlist(list(rbindlist(res), rbindlist(res2))), aes(x = variable, y = value, fill = t)) + geom_boxplot() + 
  facet_wrap(~i, scales = "free_y")

intstrength <- function(dyn, web, xpar){
  dbls <- c()
  sing <- c()
  for(i in 1:500){
    test <- Fij(dyn[i,-1], web, .5, xpar)
    ints1 <- cbind(test[upper.tri(test)], t(test)[upper.tri(test)])
    dbls[i] <- mean(ints1[which(ints1[,1] > 0 & ints1[,2] > 0),])
    sing[i] <- mean(ints1[!(ints1[,1] > 0 & ints1[,2] > 0),][rowSums(ints1) > 0][ints1[!(ints1[,1] > 0 & ints1[,2] > 0),][rowSums(ints1) > 0] > 0])
  }
  
  return(sing-dbls)
}

is1 <- lapply(1:1000, function(x){intstrength(dynRESniche2[[x]], n.mod.webs[[x]], xpar = .2)})


pathdyn <- function(dyn, web){
  apl <- c()
  for(i in 1:500){
    nw <- web[dyn[i, -1] > 0, dyn[i, -1] > 0]
    apl[i] <- diameter(graph.adjacency(nw)) # average.path.length(graph.adjacency(nw))
  }
  apl2 <- c()
  for(i in 1:200){
    s <- sample(dyn[i, -1] > 0)
    nw <- web[s, s]
    apl2[i] <- diameter(graph.adjacency(nw))
  }
  
  return(list(apl, apl2))
}

pld <- lapply(1:1000, function(x){pathdyn(dynRESniche2[[x]], n.mod.webs[[x]])})
plot(rowMeans(t(do.call(rbind, pld))))

pdyn <- lapply(pld, "[[", 1)
pnull <- lapply(pld, "[[", 2)