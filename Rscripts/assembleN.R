library(igraph)
library(NetIndices)
library(reshape2)
library(rend)


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
  
  return(inN)
}

wprop <- function(m){
  N <- nrow(m)
  L <- sum(m)
  C <- L/(N*(N-1))
  meanGen <- mean(colSums(m))
  meanVul <- mean(rowSums(m))
  sdGen <- sd(colSums(m))
  sdVul <- sd(rowSums(m))
  ti <- TrophInd(m)
  meanTP <- mean(ti$TL)
  sdTP <- sd(ti$TL)
  meanOI <- mean(ti$OI)
  sdOI <- sd(ti$OI)
  
  g <- graph.adjacency(m)
  APL <- average.path.length(g)
  D <- diameter(g)
  
  mc <- motif_counter(list(g))
  
  return(data.frame(N, L, C, meanGen, meanVul, sdGen, sdVul, meanTP, sdTP, meanOI, sdOI, APL, D, mc))
}

iprop <- function(m, m.post, invader){
  invGen <- sum(m[,ncol(m)])
  invVul <- sum(m[nrow(m),])
  ti <- TrophInd(m)
  invTP <- tail(ti$TL, 1)
  invOI <- tail(ti$OI, 1)
  invEST <- invader %in% colnames(m.post)
  spGain <- nrow(m.post) - nrow(m) 
  
  return(data.frame(invID = invader, invEST, spGain, invGen, invVul, invTP, invOI))
}

assembly <- function(niche, initial, n.inv){
  dyn1 <- CRsimulator(initial, t = 1:300)
  extant <- which(tail(dyn1, 1)[,-1] > 0) 
  extantSPP <- list(as.numeric(colnames(initial)[extant]))
  eqABUND <- tail(dyn1, 1)[,-1][which(tail(dyn1, 1)[,-1] > 0)]
  
  wp.before <- list()
  wp.after <- list()
  invp <- list()
  for(i in 2:n.inv){
    intro <- sample((1:1000)[-extantSPP[[i-1]]], 1)
    invaded <- c(extantSPP[[i-1]], intro)
    newN <- niche[invaded, invaded]
    
    if(sum(niche[,intro]) != 0){
      while(sum(newN[,nrow(newN)]) == 0){
        intro <- sample((1:1000)[-extantSPP[[i-1]]], 1)
        invaded <- c(extantSPP[[i-1]], intro)
        newN <- niche[invaded, invaded]
      }
    }
    
    r.in <- as.numeric(colnames(newN) %in% which(colSums(niche) == 0)) 
    
    # invaders initial abundance set to 0.5
    inABUND <- c(as.vector(eqABUND), 0.5)
    #dynamics
    dyn <- CRsimulator(newN, states = inABUND, r = r.in, t = 1:300)
    # who is extant
    ext1 <- which(tail(dyn, 1)[,-1] > 0)
    extantSPP[[i]] <- as.numeric(colnames(newN)[ext1])
    #abundance of extants
    eqABUND <- tail(dyn, 1)[,-1][which(tail(dyn, 1)[,-1] > 0)]
    
    wp.before[[i-1]] <- wprop(newN)
    wp.after[[i-1]] <- wprop(newN[ext1, ext1])
    invp[[i-1]] <- iprop(newN, newN[ext1, ext1], intro)
    
  }
  
  wp.pre <- melt(wp.before, measure.vars = colnames(wp.before[[1]]))
  wp.post <- melt(wp.after, measure.vars = colnames(wp.after[[1]]))
  invprop <- do.call(rbind, invp)
  
  return(list(community = extantSPP, pre = wp.pre, post = wp.post, inv = invprop))
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

trophINFO <- function(mat){
  ti <- TrophInd(mat)
  
  maxTP <- max(ti$TL)
  meanTP <- mean(ti$TL)
  sdTP <- sd(ti$TL)
  maxOI <- max(ti$OI)
  meanOI <- mean(ti$OI)
  sdOI <- sd(ti$OI)
  
  consumers <- ti$TL > 1
  
  maxTPc <- max(ti$TL[consumers])
  meanTPc <- mean(ti$TL[consumers])
  sdTPc <- sd(ti$TL[consumers])
  maxOIc <- max(ti$OI[consumers])
  meanOIc <- mean(ti$OI[consumers])
  sdOIc <- sd(ti$OI[consumers])
  
  return(c(meanTP, meanTPc, maxTP, maxTPc, sdTP, sdTPc, meanOI, meanOIc, maxOI, maxOIc, sdOI, sdOIc))
}

curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:(5*l_hp)){
    AB=sample(1:l_hp,2)
    a=hp[[AB[1]]]
    b=hp[[AB[2]]]
    ab=intersect(a,b)
    l_ab=length(ab)
    l_a=length(a)
    l_b=length(b)
    if ((l_ab %in% c(l_a,l_b))==F){
      tot=setdiff(c(a,b),ab)
      l_tot=length(tot)
      tot=sample(tot, l_tot, replace = FALSE, prob = NULL)
      L=l_a-l_ab
      hp[[AB[1]]] = c(ab,tot[1:L])
      hp[[AB[2]]] = c(ab,tot[(L+1):l_tot])}
    
  }
  rm=matrix(0,R,C)
  for (row in 1:R){rm[row,hp[[row]]]=1}
  rm
}
