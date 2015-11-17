library(igraph)
library(NetIndices)
library(modMax)
library(deSolve)
library(animation)
library(reshape2)
library(ggplot2)
library(data.table)
library(rnetcarto)
library(parallel)
library(doSNOW)

# C-R-model Functions

# growth (nonzero for basal spp only)
G.i <- function(r, B, K){return(r * B * (1 - (B/K)))}

# functional response 
Fij <- function(B, A, B.0, xpar){
  sum.bk <- rowSums(sapply(1:nrow(A), function(x){B[x] * A[x,]}))^(1+xpar)
  denom <- sum.bk + B.0^(1+xpar)
  
  F1 <- sapply(1:nrow(A), function(x){(B[x] * A[x,])^(1+xpar)})/denom
  
  return(F1)
}

# functional response with interference
Fbd <- function(B, A, B.0, xpar){
  sum.bk <- rowSums(sapply(1:nrow(A), function(x){B[x] * A[x,]}))
  denom <- sum.bk + (1 + (xpar * B)) * B.0
  
  F1 <- sapply(1:nrow(A), function(x){(B[x] * A[x,])})/denom
  
  return(F1)
}

# get r.i from adjacency matrix
get.r <- function(amat){
  r.i <- c()
  r.i[colSums(amat) == 0] <- 1
  r.i[colSums(amat) != 0] <- 0
  
  if(sum(r.i) == 0){r.i[sample(1:length(r.i), 1)] <- 1}
  return(r.i)
}

# dynamic function for ode solver
conres <- function(t,states,par){
  
  with(as.list(c(states, par)), {
    dB <- G.i(r = r.i, B = states, K = K) - x.i*states + rowSums((x.i * yij * FR(states, A, B.o, xpar = xpar) * states)) - rowSums((x.i * yij * t(FR(states, A, B.o, xpar = xpar)* states))/eij)
    
    list(c(dB))
  })
  
}

# wrapper  to input params, run dynamics, and plot
Crmod <- function(Adj, t = 1:200, G = G.i, method = conres, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = .2, B.o =.5, plot = FALSE){
  require(deSolve)
  
  grow <- get.r(Adj)
  
  par <- list(
    K = K,
    x.i = x.i,
    yij = yij,
    eij = 1,
    xpar = xpar,
    B.o = B.o,
    r.i = grow,
    A = Adj,
    G.i = G,
    FR = FuncRes
  )
  
  states <- runif(nrow(Adj), .5, 1)
 
  out <- ode(y=states, times=t, func=method, parms=par, events = list(func = eventfun, time = t))
  
  if(plot) print(matplot(out[,-1], typ = "l", lwd = 2, xlab = "Time", ylab = "Biomass"))
  
  return(out)
}

# event function to make any spp under the threshold equal to 0 (extinction)
eventfun <- function(times, states, parms){
  with(as.list(states), {
    for(i in 1:length(states)){
      if(states[i] < 10^-10){states[i] <- 0}else{states[i]} 
    }
    return(c(states))
  })
}

# Niche model functions
niche.model<-function(S,C){
  require(igraph)
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
    
    new.mat<-new.mat[,order(apply(new.mat,2,sum))]
    
    connected <- is.connected(graph.adjacency(new.mat))
  }
  return(new.mat)
}

erdos_renyi <- function(S, C){
  basal = FALSE
  while(!basal){
    ag <- erdos.renyi.game(n = S, p.or.m = C, type = "gnp", directed = T)
    amat <- get.adjacency(ag, sparse = F)
    
    basal <- sum(colSums(amat) == 0) > 1
  }
  return(amat)
}

web_maker <- function(typ, n, S, C){
  if(typ == "niche"){
    niche.list <- list()
    for (i in 1:n){
      niche.list[[i]]<- niche.model(S, C)
    }
    return(niche.list)
  }else if(typ == "erg"){
    erg.list <- list()
    for (i in 1:n){
      erg.list[[i]]<- erdos_renyi(S, C)
    }
    return(erg.list)
  }else{
    stop("This function only works for niche model ('niche') and erdos-renyi random networks ('erg')")
  }
}


# motif_counter function from "web_functions.R"
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

web_props <- function(mat){
  require(NetIndices)
  require(modMax)
  
  N <- nrow(mat)
  C <- sum(mat)/(nrow(mat)*(nrow(mat)-1))
  
  genind <- GenInd(mat)
  Ltot <- genind$Ltot
  LD <- genind$LD
  
  g <- graph.adjacency(mat)
  
  clust <- transitivity(g)
  apl <- average.path.length(g)
  diam <- diameter(g)
  
  bas <- sum(degree(g, mode = "in") == 0)
  top <- sum(degree(g, mode = "out") == 0)
  
  #mot <- motif_counter(list(graph.adjacency(mat)))
  
  df <- matrix(c(N, C, Ltot, LD, clust, apl, diam, bas, top), nrow = 1)
  colnames(df) <- c("N", "C", "Ltot", "LD", "clust", "apl", "diam", "bas", "top")
  return(df)
}

# curveball function for null model
curve_ball<-function(m){
  RC=dim(m)
  R=RC[1]
  C=RC[2]
  hp=list()
  for (row in 1:dim(m)[1]) {hp[[row]]=(which(m[row,]==1))}
  l_hp=length(hp)
  for (rep in 1:5*l_hp){
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

# wrapper for n curveball null matrices 
curving <- function(adjmat, n){
  mot <- motif_counter(list(graph.adjacency(adjmat)))
  newmat <- adjmat
  
  for(i in 1:n){
    newmat <- curve_ball(newmat)
    m <- motif_counter(list(graph.adjacency(newmat)))
    mot <- rbind(mot, m)
  }
  return(mot[-1,])
}

# analyze the network

get.zscore <- function(mat, iter = 10000){
  mot <- motif_counter(graph.lists = list(graph.adjacency(mat)))
  nullmot <- curving(adjmat = mat, iter)
  
  z.score <- as.matrix((mot-colMeans(nullmot))/apply(nullmot, 2, sd))
  return(z.score)
}

rel.rand <- function(mat, n, iter = 10000){
  require(igraph)
  rands <- list()
  for(i in 1:iter){
    spp <- sample(1:nrow(mat), n)
    rands[[i]] <- mat[spp, spp]
  }
  rands.g <- lapply(rands, graph.adjacency)
  mot <- motif_counter(rands.g)
  return(mot)
}

rel.rand.z <- function(mat, dyn, iter = 10000){
  m1 <- motif_counter(list(graph.adjacency(mat[which(tail(dyn, 1)[-1] > 0),which(tail(dyn, 1)[-1] > 0)])))
  n <- sum(tail(dyn, 1)[-1] > 0)
  m2 <- rel.rand(mat = mat, n = n, iter = iter)
  
  zscore <- (m1 - colMeans(m2))/apply(m2, 2, sd)
  
  return(zscore)
}

# multiple iterations

CRmod.gen <- function(S, con, network, xpar, ...){
  # create model network
  if(network == "niche"){amat <- niche.model(S = S, C = con)}
  if(network == "erdosrenyi"){
    ag <- erdos.renyi.game(n = S, p.or.m = con, type = "gnp", directed = T)
    amat <- get.adjacency(ag, sparse = F)
  }
  
  # dynamic model
  
  dyn1 <- Crmod(Adj = amat, t = 1:500, G = G.i, method = conres, FuncRes = Fij, K = 1, x.i = .5, yij = 6, eij = 1, xpar = xpar, B.o =.5, plot = FALSE)
  
  # get motif scores
  
  amat2 <- amat[which(tail(dyn1[,-1], 1) > 0),which(tail(dyn1[,-1], 1) > 0)]
  
  #z1 <- data.frame(time = factor("init"), get.zscore(amat))
  #z2 <- data.frame(time = factor("final"), get.zscore(amat2))
  #return(rbind(z1, z2))
  
  z <- rel.rand.z(amat, dyn1, iter = 10000)
  z1 <- z/sqrt(sum(z^2, na.rm = T))
  t(apply(zscore.t, 1, function(x){x/sqrt(sum(x^2, na.rm = T))}))
  wp <- cbind(melt(web_props(amat), variable.name = "property", value.name = "initial"), final = melt(web_props(amat2))[,2])
  
  return(list(wp = wp, z = z1))
}


# visualize 

vis_dyn <- function(dyn, t.start = 1, t.end = 200){
  require(ggplot2)
  g <- ggplot(melt(dyn[t.start:t.end,-1]), aes(x = Var1, y = value, col = factor(Var2))) + geom_line() + theme(legend.position="none")
  print(g)
}

netHTML <- function(mat, dyn, path1 = getwd()){
  require(animation)
  
  lay <- matrix(c(layout.sphere(graph.adjacency(nm1))[,1], TrophInd(nm1)$TL), ncol = 2)
  s <- matrix(0, nrow = nrow(dyn), ncol = ncol(mat))
  
  ani.options(interval = .25)
  saveHTML(
    {
      for(i in 1:50){
        fr <- Fij(dyn2[i,-1], nm1, .5, .2)
        strength <- melt(fr)[,3][melt(fr)[,3] > 0]
        fr[fr > 0 ] <- 1
        
        g.new <- graph.adjacency(t(fr))
        E(g.new)$weight <- strength/max(strength)*10
        s[i,c(which(dyn2[i,-1] > 0))] <-log(dyn2[i, c(which(dyn2[i,] > 0)[-1])])+abs(min(log(dyn2[i, c(which(dyn2[i,] > 0)[-1])])))
        
        plot.igraph(g.new, vertex.size = s[i,], edge.width = E(g.new)$weight, layout = lay)
      }
    },
    htmlfile = "fwdyn.html", interval = .25, nmax =500, ani.width = 500, ani.height = 500,
    outdir = path1
  )
}


# QSS functions
conversion <- function(tm){
  for(i in 1:nrow(tm)){
    for(j in 1:ncol(tm)){
      if(tm[i,j] == 1 & tm[j,i] == 0){tm[j,i] <- -1}
    }
  }
  return(tm)
}

ran.unif <- function(motmat){
  newmat <- apply(motmat, c(1,2), function(x){
    if(x==1){runif(1, 0, 10)}else if(x==-1){runif(1, -1, 0)} else{0}
  })
  diag(newmat) <- runif(nrow(newmat), -1, 0)
  return(newmat)
}

maxRE <- function(rmat){
  lam.max <- max(Re(eigen(rmat)$values))
  return(lam.max)
}

eig.analysis <- function(n, matrices){
  cols <- length(matrices)
  rows <- n
  eigenMATRIX <- matrix(0, nrow = rows, ncol = cols)
  for(i in 1:n){
    ranmat <- lapply(matrices, ran.unif)
    
    eigs <- sapply(ranmat, maxRE)
    eigenMATRIX[i,] <- eigs
  }
  return(eigenMATRIX)
}

# For slides
set.seed(12)
nm1 <- niche.model(S = 30, C = .1)

## For parallel
tot = 100
modTYPE <- c(rep("erdosrenyi", tot/2), rep("niche", tot/2))

results <- lapply(1:tot, function(x) data.frame(rep(0, 10), rep(0, 10), rep(0, 10), rep(0, 10), rep(0, 10)))
results2 <- lapply(1:tot, function(x) data.frame(rep(0, 10), rep(0, 10), rep(0, 10), rep(0, 10), rep(0, 10)))
z.results <- lapply(1:tot, function(x) data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
z.results2 <- lapply(1:tot, function(x) data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

z.res <- lapply(1:tot, function(x){rbind(data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))})

res <- lapply(1:1, function(x) data.frame(rep(0, 20), rep(0, 20), rep(0, 20), rep(0, 20), rep(0, 20)))
