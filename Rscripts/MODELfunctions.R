# C-R-model Functions

# growth (nonzero for basal spp only)
G.i <- function(r, B, K){return(r * B * (1 - (B/K)))}
# death
d.i <- function(x, B){return(x * B)}


# impact of the resource on the consumer
Fij <- function(B, amat, B.o, q){
  denom <- (amat %*% B)^(1+q) + B.o^(1+q)
  numer <- (t(apply(amat, 1, function(x){x * B})))^(1+q)
  return(apply(numer, 2, function(x){x/denom}))
}


# impact of the consumer on the resource
Fji <- function(B, amat, B.o, q){
  denom <- (amat %*% B)^(1+q) + .5^(1+q)
  numer <- (t(apply(amat, 1, function(x){x * B})))^(1+q)
  apply(numer, 2, function(x){x/denom})
}


# Function for change of biomass in one timestep
#amat is input as the adjacency matrix of the network
dBdt <- function(r.i, B.i, K.i, x.i, yij, amat, q, eij){
  growth <- G.i(r.i, B.i, K.i) 
  death <- d.i(x.i, B.i) 
  funct <- x.i * yij * Fij(B.i, t(amat), .5, q) %*% B.i - (x.i * yij * Fji(B.i, amat, .5, q) %*% B.i)/eij
  return(growth - death + funct)
}


# Full model
web_dyn <- function(n.times = 100, amat, B.i, r.i, K.i, x.i = .5, yij = .6, q = .2, eij = 1, stochastic = FALSE, ext.thres = 10^-10){
  r.i2 <- r.i
  # create matrix with initial abundances
  B <- matrix(B.i, nrow = n.times, ncol = ncol(amat), byrow = T)
  for(i in 2:n.times){
    # stochastic growth rate
    if(stochastic){r.i2[r.i == 1] <- sapply(r.i[r.i == 1], function(x){rnorm(1, x, .1)})}
    
    # next time step
    B[i,] <- B[i-1,] + dBdt(r.i2, B.i = B[i-1,], K.i, x.i, yij, amat, q, eij)
    
    # Check extinctions
    B[i,][B[i,] < ext.thres] <- 0
  }
  return(B)
}

# get r.i from adjacency matrix
get.ri <- function(amat){
  r.i <- c()
  r.i[colSums(amat) == 0] <- 1
  r.i[colSums(amat) != 0] <- 0
  
  if(sum(r.i) == 0){r.i[sample(1:length(r.i), 1)] <- 1}
  return(r.i)
}

# Niche model functions
niche.model<-function(S,C){
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
