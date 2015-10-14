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
  
  if(plot) print(matplot(out[,-1], typ = "l", lwd = 2))
  
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



