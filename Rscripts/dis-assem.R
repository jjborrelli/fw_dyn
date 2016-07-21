library(igraph)
library(NetIndices)
library(rend)
library(data.table)

filepath.sink <- "D:/jjborrelli/AssemblyDATA/invdel/"

niche.model<-function(S,C){
  connected = FALSE
  while(!connected){  
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
    
    new.mat<-new.mat[order(apply(new.mat,2,sum)), order(apply(new.mat,2,sum))]
    
    connected <- is.connected(graph.adjacency(new.mat))
  }
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


netdiff <- function(mat1, mat2){
  w.in <- wprop(mat1)
  w.fi <- wprop(mat2)
  
  return(w.fi - w.in)
}

initialNiche <- function(niche, S.local, nbasal){
  sppNAMES <- as.character(1:nrow(niche))
  colnames(niche) <- sppNAMES
  rownames(niche) <- sppNAMES
  cond <- FALSE
  #get a connected initial network
  while(!cond){
    #initial basal spp
    inBAS <- sample(which(colSums(niche) == 0), nbasal)
    #initial consumer spp
    inCON <- sample((1:1000)[-which(colSums(niche) == 0)], S.local)
    #all spp
    inSPP <- c(inBAS, inCON)
    #intial matrix
    inN <- niche[inSPP, inSPP]
    cond <- sum(colSums(inN)[names(colSums(inN)) %in% as.character(inCON)] == 0) == 0
  }
  
  return(inN)
}

initialize <- function(Sregion, Cregion, n.initial, Slocal = 45, Sbasal = 5, times = 1000, fr = Fij, frtune = 0.2){
  n.regional <- niche.model(Sregion, Cregion)
  
  in.n <- lapply(1:n.initial, function(x) initialNiche(n.regional, S.local = Slocal, nbasal = Sbasal))
  
  cl <- makeCluster(detectCores() - 1)
  registerDoSNOW(cl)
  in.dyn <- parLapply(cl, in.n, CRsimulator, t = 1:times, FuncRes = fr, xpar = frtune)
  stopCluster(cl)
  
  return(list(REGIONAL = n.regional, INITIAL = in.n, DYNAMICS = in.dyn))
}


iprop <- function(m, invader){
  invGen <- sum(m[,ncol(m)])
  invVul <- sum(m[nrow(m),])
  ti <- TrophInd(m)
  invTP <- tail(ti$TL, 1)
  invOI <- tail(ti$OI, 1)
  
  return(data.frame(invID = invader, invGen, invVul, invTP, invOI))
}

delprop <- function(m, d){
  delGen <- sum(m[,d])
  delVul <- sum(m[d,])
  ti <- TrophInd(m)
  delTP <- ti$TL[d]
  delOI <- ti$OI[d]
  
  return(data.frame(delID = d, delGen, delVul, delTP, delOI))
}


invading <- function(initial, regional){
  dyn1 <- CRsimulator(initial, t = 1:500)
  extant <- which(tail(dyn1, 1)[,-1] > 0) 
  extantSPP <- as.numeric(colnames(initial)[extant])
  eqABUND <- tail(dyn1, 1)[,-1][extant]
  
  eq.n <- initial[extant, extant]
  intro <- sample((1:nrow(regional))[-extantSPP], 1)
  invaded <- c(extantSPP, intro)
  
  newN <- regional[invaded, invaded]

  
  if(sum(regional[,intro]) != 0){
    while(sum(newN[,nrow(newN)]) == 0){
      intro <- sample((1:nrow(regional))[-extantSPP], 1)
      invaded <- c(extantSPP, intro)
      newN <- regional[invaded, invaded]
    }
  }
  colnames(newN) <- as.character(invaded)
  rownames(newN) <- as.character(invaded)
  
  nd1 <- netdiff(eq.n, newN)
  
  r.in <- as.numeric(colnames(newN) %in% which(colSums(regional) == 0))
  inABUND <- c(as.vector(eqABUND), 0.5)
  
  dyn <- CRsimulator(newN, states = inABUND, r = r.in, t = 1:500)
  
  ext1 <- which(tail(dyn, 1)[,-1] > 0)
  eq.n2 <- newN[ext1, ext1]
  
  in.success <- intro %in% colnames(eq.n2)
  sp.change <- nrow(eq.n2) - nrow(newN) 
  inv.tte <- which.min(dyn[,nrow(newN)])

  inv.prop <- iprop(newN, intro)
  inv.prop$dN <- sp.change
  inv.prop$tte <- inv.tte
  inv.prop$I <- in.success
  
  return(list(nd1, inv.prop))
}

invading2 <- function(dyn1, initial, regional, FR = Fij, frtune = 0.2){
  extant <- which(tail(dyn1, 1)[,-1] > 0) 
  extantSPP <- as.numeric(colnames(initial)[extant])
  eqABUND <- tail(dyn1, 1)[,-1][extant]
  
  eq.n <- initial[extant, extant]
  intro <- sample((1:nrow(regional))[-extantSPP], 1)
  invaded <- c(extantSPP, intro)
  
  newN <- regional[invaded, invaded]
  
  
  if(sum(regional[,intro]) != 0){
    while(sum(newN[,nrow(newN)]) == 0){
      intro <- sample((1:nrow(regional))[-extantSPP], 1)
      invaded <- c(extantSPP, intro)
      newN <- regional[invaded, invaded]
    }
  }
  colnames(newN) <- as.character(invaded)
  rownames(newN) <- as.character(invaded)
  
  nd1 <- netdiff(eq.n, newN)
  
  r.in <- as.numeric(colnames(newN) %in% which(colSums(regional) == 0))
  inABUND <- c(as.vector(eqABUND), 0.5)
  
  dyn <- CRsimulator(newN, states = inABUND, r = r.in, t = 1:500, FuncRes = FR, xpar = frtune)
  
  ext1 <- which(tail(dyn, 1)[,-1] > 0)
  eq.n2 <- newN[ext1, ext1]
  
  in.success <- intro %in% colnames(eq.n2)
  sp.change <- nrow(eq.n2) - nrow(newN) 
  
  invdyn <- dyn[,(nrow(newN) + 1)]
  if(invdyn[500] != 0){inv.tte <- 500}else{inv.tte <- min(which(invdyn == 0))}
 
  inv.prop <- iprop(newN, intro)
  inv.prop$dN <- sp.change
  inv.prop$tte <- inv.tte
  inv.prop$I <- in.success
  
  return(list(nd1, inv.prop))
}

deleting <- function(initial, regional){
  dyn1 <- CRsimulator(initial, t = 1:500)
  extant <- which(tail(dyn1, 1)[,-1] > 0) 
  extantSPP <- as.numeric(colnames(initial)[extant])
  eqABUND <- tail(dyn1, 1)[,-1][extant]
  
  eq.n <- initial[extant, extant]
  
  d.ch <- list()
  nd1 <- list()
  persist <- c()
  s.initial <- c() 
  s.final <- c()
  for(i in 1:nrow(eq.n)){
    d.ch[[i]] <- delprop(eq.n, i)
    newN <- eq.n[-i,-i]
    nd1[[i]] <- netdiff(eq.n, newN)
    r.in <- as.numeric(colnames(newN) %in% which(colSums(regional) == 0))
    inABUND <- eqABUND[-i]
    
    dyn <- CRsimulator(newN, states = inABUND, r = r.in, t = 1:500)
    
    ext1 <- which(tail(dyn, 1)[,-1] > 0)
    persist[i] <- length(ext1)/nrow(newN)
    s.initial[i] <- nrow(newN)
    s.final[i] <- length(ext1)
  }
  
  nd2 <- data.frame(persist = persist, do.call(rbind, nd1), do.call(rbind, d.ch))  
  nd2$s.in <- s.initial
  nd2$s.fi <- s.final 
  return(nd2)
}

deleting2 <- function(dyn1, initial, regional, FR = Fij, frtune = 0.2){
  extant <- which(tail(dyn1, 1)[,-1] > 0) 
  extantSPP <- as.numeric(colnames(initial)[extant])
  eqABUND <- tail(dyn1, 1)[,-1][extant]
  
  eq.n <- initial[extant, extant]
  
  d.ch <- list()
  nd1 <- list()
  persist <- c()
  s.initial <- c() 
  s.final <- c()
  for(i in 1:nrow(eq.n)){
    d.ch[[i]] <- delprop(eq.n, i)
    newN <- eq.n[-i,-i]
    nd1[[i]] <- netdiff(eq.n, newN)
    r.in <- as.numeric(colnames(newN) %in% which(colSums(regional) == 0))
    inABUND <- eqABUND[-i]
    
    dyn <- CRsimulator(newN, states = inABUND, r = r.in, t = 1:500, FuncRes = FR, xpar = frtune)
    
    ext1 <- which(tail(dyn, 1)[,-1] > 0)
    persist[i] <- length(ext1)/nrow(newN)
    s.initial[i] <- nrow(newN)
    s.final[i] <- length(ext1)
  }
  
  nd2 <- data.frame(persist = persist, do.call(rbind, nd1), do.call(rbind, d.ch))  
  nd2$s.in <- s.initial
  nd2$s.fi <- s.final 
  return(nd2)
}


