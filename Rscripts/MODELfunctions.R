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


niche_model <- function(S, C){
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

rel.rand <- function(mat, n, iter = 100){
  require(igraph)
  rands <- list()
  for(i in 1:iter){
    spp <- sample(1:nrow(mat), n)
    rands[[i]] <- mat[spp, spp]
  }
  rands.g <- lapply(rands, netcarto)
  mod <- sapply(rands.g, "[[", 2)
  return(mod)
}

rel.rand.z <- function(mat, dyn, iter = 10000){
  m1 <- motif_counter(list(graph.adjacency(mat[which(tail(dyn, 1)[-1] > 0),which(tail(dyn, 1)[-1] > 0)])))
  n <- sum(tail(dyn, 1)[-1] > 0)
  m2 <- rel.rand(mat = mat, n = n, iter = iter)
  
  zscore <- (m1 - colMeans(m2))/apply(m2, 2, sd)
  
  return(zscore)
}


# motif loss through time

motif_loss <- function(dmat, init){
  mc <- list()
  for(i in 1:nrow(dmat)){
    dfin <- dmat[i, -1]
    nfin <- init[dfin > 0, dfin > 0]
    
    mc[[i]] <- motif_counter(list(graph.adjacency(nfin)))
  }
  return(rbindlist(mc))
}

