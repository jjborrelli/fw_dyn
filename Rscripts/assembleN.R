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

for(i in 2:101){
  invaded <- c(extantSPP[[i-1]], sppINTRO[i-1])
  newN <- n1[invaded, invaded]
  
  dyn <- CRsimulator(newN)
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
  colnames(n1) <- sppNAMES
  rownames(n1) <- sppNAMES
  cond <- FALSE
  #get a connected initial network
  while(!cond){
    #initial basal spp
    inBAS <- sample(which(colSums(n1) == 0), nbasal)
    #initial consumer spp
    inCON <- sample((1:1000)[-which(colSums(n1) == 0)], 40)
    #all spp
    inSPP <- c(inBAS, inCON)
    #intial matrix
    inN <- n1[inSPP, inSPP]
    cond <- is.connected(graph.adjacency(inN))
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
    
    dyn <- CRsimulator(newN)
    ext1 <- which(tail(dyn, 1)[,-1] > 0)
    extantSPP[[i]] <- as.numeric(colnames(newN)[ext1])
    eqABUND[[i]] <- tail(dyn, 1)[,-1]
  }
  
  return(list(community = extantSPP, abundances = eqABUND))
}

motifChanges <- function(community, niche){
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

cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("assembly", "initialnets", "sppintro", "regional"))
registerDoSNOW(cl)

clusterCall(cl, function() library(rend))
combuild <- parLapply(cl, 1:100, function(x){assembly(regional, initialnets[[x]], sppintro[[x]], n.inv = 100)})

stopCluster(cl)
