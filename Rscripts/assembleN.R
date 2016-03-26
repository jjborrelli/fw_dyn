library(igraph)
library(NetIndices)
library(reshape2)
library(rend)
library(ggplot2)

niche.model <- function(S,C){
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
  
  new.mat<-new.mat[order(apply(new.mat,2,sum)),order(apply(new.mat,2,sum))]
  
  return(new.mat)
}

motif_counter <- function(graph.lists){
  require(igraph)
  
  if(!is.list(graph.lists)){
    stop("The input should be a list of graph objects")
  }
  
  triad.count <- lapply(graph.lists, triad.census)
  triad.matrix <- matrix(unlist(triad.count), nrow = length(graph.lists), ncol = 16, byrow = T)
  colnames(triad.matrix) <- c("empty", "single", "mutual", "s5", "s4", "s1", "d4",
                              "d3", "s2", "s3","d8", "d2", "d1", "d5", "d7", "d6")
  
  triad.df <- as.data.frame(triad.matrix)
  
  motif.data.frame <- data.frame(s1 = triad.df$s1, s2 = triad.df$s2, s3 = triad.df$s3, s4 = triad.df$s4, 
                                 s5 = triad.df$s5, d1 = triad.df$d1, d2 = triad.df$d2, d3 = triad.df$d3, d4 = triad.df$d4,
                                 d5 = triad.df$d5, d6 = triad.df$d6, d7 = triad.df$d7, d8 = triad.df$d8)
  
  return(motif.data.frame)
}

nullN <- function(n, s, bas, iter = 200){
  lapply(1:iter, function(x){
    s1 <- sample((1:1000)[-bas], s)
    s2 <- c(bas, s1)
    n[s2, s2]
  })
}


n1 <- niche.model(1000, .15)
sppNAMES <- as.character(1:1000)
colnames(n1) <- sppNAMES
rownames(n1) <- sppNAMES
cond <- FALSE
#get a connected initial network
while(!cond){
  #initial basal spp
  inBAS <- sample(which(colSums(n1) == 0), 5)
  #initial consumer spp
  inCON <- sample((1:1000)[-which(colSums(n1) == 0)], 40)
  #all spp
  inSPP <- c(inBAS, inCON)
  #intial matrix
  inN <- n1[inSPP, inSPP]
  cond <- is.connected(graph.adjacency(inN))
}

sppINTRO <- sample((1:1000)[-inSPP], 100)

dyn1 <- CRsimulator(inN)
extant <- which(tail(dyn1, 1)[,-1] > 0) 
extantSPP <- list(as.numeric(colnames(inN)[extant]))
eqABUND <- list(tail(dyn1, 1)[,-1])

for(i in 2:10){
  invaded <- c(extantSPP[[i-1]], sppINTRO[i-1])
  newN <- n1[invaded, invaded]
  inABUND <- c(as.vector(eqABUND[[i-1]][as.vector(eqABUND[[i-1]]) > 0]), runif(1, .5, 1))
  dyn <- CRsimulator(newN, states = inABUND)
  ext1 <- which(tail(dyn, 1)[,-1] > 0)
  extantSPP[[i]] <- as.numeric(colnames(newN)[ext1])
  eqABUND[[i]] <- tail(dyn, 1)[,-1]
  
  cat(i-1, "------->", length(extantSPP[[i]]), "species", "\n")
}

nSPP <- sapply(extantSPP, length)
whichbas <- lapply(extantSPP, function(x){which(x %in% which(colSums(n1) == 0))})
allnets <- lapply(extantSPP, function(x){n1[x,x]})
allg <- lapply(allnets, graph.adjacency)



mc <- motif_counter(allg)
nulls <- lapply(1:length(allg), function(x){
  nn <-nullN(n1, nSPP[x], whichbas[[x]])
  m1 <- motif_counter(lapply(nn, graph.adjacency))
  return(list(colMeans(m1), apply(m1, 2, sd)))
})

nullmeans <- lapply(nulls, "[[", 1)
nullsds <- lapply(nulls, "[[", 2)

zmc <- t(sapply(1:nrow(mc), function(x){(unlist(mc[x,]) - unlist(nullmeans[[x]]))/unlist(nullsds[[x]])}))
zmc[is.nan(zmc)] <- 0
boxplot(zmc)
ggplot(melt(zmc), aes(x = Var1, y = value)) + geom_point() + stat_smooth(method = "lm") + facet_wrap(~Var2)


####################################################################
####################################################################
####################################################################

initialNiche <- function(niche, nbasal){
  sppNAMES <- as.character(1:nrow(niche))
  colnames(niche) <- sppNAMES
  rownames(niche) <- sppNAMES
  cond <- FALSE
  #get a connected initial network
  while(!cond){
    #initial basal spp
    inBAS <- sample(which(colSums(niche) == 0), nbasal)
    #initial consumer spp
    inCON <- sample((1:1000)[-which(colSums(niche) == 0)], 40)
    #all spp
    inSPP <- c(inBAS, inCON)
    #intial matrix
    inN <- niche[inSPP, inSPP]
    cond <- sum(colSums(inN)[names(colSums(inN)) %in% as.character(inCON)] == 0) == 0
  }
  
  sppINTRO <- sample((1:1000)[-inSPP], 100)
  
  return(list(inN = inN, sppINTRO = sppINTRO))
}

assembly <- function(niche, initial, intro, n.inv){
  dyn1 <- CRsimulator(initial)
  extant <- which(tail(dyn1, 1)[,-1] > 0) 
  extantSPP <- list(as.numeric(colnames(initial)[extant]))
  eqABUND <- list(tail(dyn1, 1)[,-1])
  
  for(i in 2:n.inv){
    invaded <- c(extantSPP[[i-1]], intro[i-1])
    newN <- niche[invaded, invaded]
    if(sum(newN[,ncol(newN)]) == 0 & sum(niche[, intro[i-1]]) != 0){
      
    }
    
    # invaders initial abundance set to 0.5
    inABUND <- c(as.vector(eqABUND[[i-1]][as.vector(eqABUND[[i-1]]) > 0]), 0.5)
    dyn <- CRsimulator(newN, states = inABUND)
    ext1 <- which(tail(dyn, 1)[,-1] > 0)
    extantSPP[[i]] <- as.numeric(colnames(newN)[ext1])
    eqABUND[[i]] <- tail(dyn, 1)[,-1]
  }
  
  return(list(community = extantSPP, abundances = eqABUND))
}

motifChanges <- function(community, niche){
  whichbas <- lapply(community, function(x){which(x %in% which(colSums(niche) == 0))})
  nSPP <- sapply(community, length)
  allnets <- lapply(community, function(x){niche[x,x]})
  allg <- lapply(allnets, graph.adjacency)
  mc <- motif_counter(allg)
  nulls <- lapply(1:length(allg), function(x){
    nn <-nullN(n1, nSPP[x], whichbas[[x]])
    m1 <- motif_counter(lapply(nn, graph.adjacency))
    return(list(colMeans(m1), apply(m1, 2, sd)))
  })
  
  nullmeans <- lapply(nulls, "[[", 1)
  nullsds <- lapply(nulls, "[[", 2)
  
  zmc <- t(sapply(1:nrow(mc), function(x){(unlist(mc[x,]) - unlist(nullmeans[[x]]))/unlist(nullsds[[x]])}))
  zmc[is.nan(zmc)] <- 0
  
  return(zmc)
}

absMOTIF <- function(community, niche){
  allnets <- lapply(community, function(x){niche[x,x]})
  allg <- lapply(allnets, graph.adjacency)
  mc <- motif_counter(allg)
}


# generates regional species pool food web
system.time(
  regional <- niche.model(1000, .15)
) # 3.69 user

# generates local equilibrium web from subset of regional
system.time(
  local <- initialNiche(regional, nbasal = 5)
) # 0.04 user

# simulate a series of 100 invasions of random spp from regional pool
system.time(
  comm.assem <- assembly(regional, local$inN, local$sppINTRO, n.inv = 100) 
) # 1567.3 user (26 min)

# get motif changes
system.time(
  moch <- motifChanges(comm.assem$community, regional)
) # 46.5 user

ggplot(melt(moch), aes(x = Var1, y = value)) + geom_point() + stat_smooth(method = "lm") + facet_wrap(~Var2)


############################################################
############################################################
############################################################

regional <- niche.model(1000, .15)
colnames(regional) <- as.character(1:nrow(regional))
rownames(regional) <- as.character(1:nrow(regional))

locals <- lapply(1:100, function(x){initialNiche(regional, nbasal = 5)})

initialnets <- lapply(locals, "[[", 1)
sppintro <- lapply(locals, "[[", 2)


library(parallel)
library(doSNOW)

strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("assembly", "initialnets", "sppintro", "regional"))
registerDoSNOW(cl)

clusterCall(cl, function() library(rend))
combuild <- parLapply(cl, 1:100, function(x){assembly(regional, initialnets[[x]], sppintro[[x]], n.inv = 100)})

stopCluster(cl)
ends <- Sys.time() # 5 hrs for 100 webs and 100 invasions

deltaS <- lapply(lapply(1:100, function(x) combuild[[x]][[1]]), function(y){sapply(y, length)})
dS <- t(do.call(rbind, deltaS))
mean(dS[100,])
sd(dS[100,])

moch <- lapply(1:100, function(x){motifChanges(combuild[[x]][[1]], regional)})
mots <- lapply(1:100, function(x){absMOTIF(combuild[[x]][[1]], regional)})
motifdiff <- lapply(mots, function(x){
  newmat <- matrix(0, nrow = nrow(x), ncol = 13)
  for(i in 2:nrow(x)){
    newmat[i,] <- unlist(x[i, ]) - unlist(x[i-1,])
  }
  return(newmat)
})

comms <- lapply(1:100, function(x) combuild[[x]][[1]])
invader <- lapply(1:100, function(x){
  inv <- c()
  delS <- c()
  delL <- c()
  apl <- c()
  mot <- matrix(0, nrow = 99, ncol = 13)
  colnames(mot) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
  for(i in 2:100){
    inv[i-1] <- sppintro[[x]][i-1] %in% comms[[x]][[i]]
    delS[i-1] <- length(comms[[x]][[i]]) - length(comms[[x]][[i-1]])  
    rnet <- regional[c(sppintro[[x]][[i-1]], comms[[x]][[i-1]]), c(sppintro[[x]][[i-1]], comms[[x]][[i-1]])]
    rnet2 <- regional[comms[[x]][[i-1]], comms[[x]][[i-1]]]
    delL[i-1] <- sum(rnet) - sum(rnet2)
    apl[i-1] <- average.path.length(graph.adjacency(rnet)) - average.path.length(graph.adjacency(rnet2))
    mot[i-1,] <- unlist(motif_counter(list(graph.adjacency(rnet)))) - unlist(motif_counter(list(graph.adjacency(rnet2)))) 
  }
  return(data.frame(inv, delS, delL, apl, mot))
})


allinv <- do.call(rbind, invader)

inv.ext <- allinv[allinv$inv & allinv$delS < 0,]
inv.next <- allinv[allinv$inv & allinv$delS > 0,]
ninv.ext <- allinv[!allinv$inv & allinv$delS < 0,]
ninv.next <- allinv[!allinv$inv & allinv$delS == 0,]

outcomes <- list(IE = inv.ext, InE = inv.next, nIE = ninv.ext, nInE = ninv.next)

websub <- lapply(outcomes, function(x){(x[,5:17])})
wsub <- melt(websub)
ymin <- function(y){mean(y) - 1.96*(sd(y)/length(y))}
ymax <- function(y){mean(y) + 1.96*(sd(y)/length(y))}

a1 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), mean)
a2 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymax)
a3 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymin)

ws.a <- data.frame(a1, ymin = a3$x, ymax = a2$x)

ggplot(ws.a, aes(x = Group.2, y = x)) + geom_point() + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + facet_wrap(~Group.1, scales = "free_y")


wp <- melt(webp.out)
ggplot(wp, aes(x = L1, y = value)) + geom_point(position = "jitter", alpha = .1) + geom_errorbar() + facet_wrap(~variable, scales = "free_y")

mu.out <- lapply(outcomes, function(x){sapply(x, function(y){c(m = mean(y), up = mean(y) + sd(y)/sqrt(length(y)), low = mean(y) - sd(y)/sqrt(length(y)))})})
d1 <- melt(mu.out)
ggplot(d1, aes(x = L1, y = value, col = Var1)) + geom_point() + facet_wrap(~Var2, scales = "free_y")

# means and errors of web props for four outcomes
ggplot(d1[d1$Var2 == "delS" | d1$Var2 == "delL" | d1$Var2 == "apl", ], aes(x = L1, y = value, col = Var1)) + geom_point(size= 5) + facet_wrap(~Var2, scales = "free_y")

# means and error of subgraphs for four outcomes
ggplot(d1[d1$Var2 != "delS" & d1$Var2 != "delL" & d1$Var2 != "apl" & d1$Var2 != "inv", ], aes(x = L1, y = value, col = Var1)) + geom_point(size= 5) + facet_wrap(~Var2, scales = "free_y")



invaderMOT <- lapply(1:100, function(x){
  mot <- matrix(0, nrow = 99, ncol = 13)
  colnames(mot) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
  for(i in 2:100){
    inv[i-1] <- sppintro[[x]][i-1] %in% comms[[x]][[i]]
    rnet <- graph.adjacency(regional[c(sppintro[[x]][[i-1]], comms[[x]][[i-1]]), c(sppintro[[x]][[i-1]], comms[[x]][[i-1]])])
    rnet2 <- graph.adjacency(regional[comms[[x]][[i-1]], comms[[x]][[i-1]]])
    rn2mot <- motif_counter(list(rnet2))
    rmot <- matrix(nrow = 200, ncol = 13)
    for(j in 1:200){
      spp <- c(sample((1:1000)[-comms[[x]][[i-1]]]), comms[[x]][[i-1]])
      rand <- graph.adjacency(regional[spp,spp])
      rmot[j,] <- unlist(motif_counter(list(rand))) - unlist(rn2mot)
    }
    diffs <- unlist(motif_counter(list(rnet))) - unlist(rn2mot)
    mot[i-1,] <- (diffs - colMeans(rmot))/apply(rmot, 2, sd)
  }
  return(data.frame(mot))
})

allinv2 <- do.call(rbind, invaderMOT)

inv.ext <- allinv2[allinv$inv & allinv$delS < 0,]
inv.next <- allinv2[allinv$inv & allinv$delS > 0,]
ninv.ext <- allinv2[!allinv$inv & allinv$delS < 0,]
ninv.next <- allinv2[!allinv$inv & allinv$delS == 0,]

outcomes <- list(IE = inv.ext, InE = inv.next, nIE = ninv.ext, nInE = ninv.next)

websub <- lapply(outcomes, function(x){(x[,5:17])})
wsub <- melt(websub)
ymin <- function(y){mean(y) - 1.96*(sd(y)/length(y))}
ymax <- function(y){mean(y) + 1.96*(sd(y)/length(y))}

a1 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), mean)
a2 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymax)
a3 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymin)

ws.a <- data.frame(a1, ymin = a3$x, ymax = a2$x)

ggplot(ws.a, aes(x = Group.2, y = x)) + geom_point() + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + facet_wrap(~Group.1, scales = "free_y")