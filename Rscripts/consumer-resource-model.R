# Predator prey growth model based on Romanuk et al 2009
library(igraph)

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

# Check functional response
res1 <- c()
for(i in 1:length(seq(.1, 10, .025))){res1[i] <- Fij(seq(.1, 10, .025)[i], amat = matrix(1), .5, q = 0)}
res2 <- c()
for(i in 1:length(seq(.1, 10, .025))){res2[i] <- Fij(seq(.1, 10, .025)[i], amat = matrix(1), .5, q = .2)}
res3 <- c()
for(i in 1:length(seq(.1, 10, .025))){res3[i] <- Fij(seq(.1, 10, .025)[i], amat = matrix(1), .5, q = 1)}

plot(res1, typ = "l", col = "orange")
points(res2, typ = "l", col = "blue")
points(res3, typ = "l", col = "purple", xpd = F)


# Function for change of biomass in one timestep
#amat is input as the adjacency matrix of the network
dBdt <- function(r.i, B.i, K.i, x.i, yij, amat, q, eij){
  growth <- G.i(r.i, B.i, K.i) 
  death <- d.i(x.i, B.i) 
  funct <- x.i * yij * Fij(B.i, t(amat), .5, q) %*% B.i - (x.i * yij * Fji(B.i, amat, .5, q) %*% B.i)/eij
  return(growth - death + funct)
}

# Parameters
r.i <- c(1,1,0) # initial abundances
K.i <- c(1,1,1) # carrying capacity
x.i <- c(.5, .5, .5)
yij <- .6
eij <- 1
q <- .2

# adjacency matrix
amat <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 0), nrow = 3)
n.amat <- matrix(c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), nrow = 4)

B.i <- runif(3, .5, 1) # initial abundance
B2 <- matrix(B.i, nrow = 100, ncol = 3, byrow = T)
for(i in 2:100){
  r.i2 <- r.i#sapply(r.i, function(x){rnorm(1, x, .2)})
  B2[i,] <- B2[i-1,] + dBdt(r.i2, B.i = B2[i-1,], K.i, x.i, yij, t(amat), q, eij)
  B2[i,][B2[i,] < .1] <- 0
}
matplot(B2, type = "l")
#b.init <-colMeans(B)
#b.init.v <- apply(B, 2, sd)


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

S = 300
nm <- niche.model(S, .05)
r.i <- c()
r.i[colSums(nm) == 0] <- 1
r.i[colSums(nm) != 0] <- 0
r.i2 <- r.i

K.i <- rep(1, S) # carrying capacity
x.i <- rep(.5, S)
yij <- .6
eij <- 1
q <-  1
amat <- nm

B.i <- runif(S, .5, 1) # initial abundance
B <- matrix(B.i, nrow = 2000, ncol = S, byrow = T)
for(i in 2:2000){
  r.i2[colSums(nm) == 0] <- sapply(r.i[colSums(nm) == 0], function(x){rnorm(1, x, .1)})
  B[i,] <- B[i-1,] + dBdt(r.i2, B.i = B[i-1,], K.i, x.i, yij, amat, q, eij)
  B[i,][B[i,] < 10^-10] <- 0
}
matplot(B, type = "l")

which(B[100,] != 0)
length(which(B[100,] != 0))

apply(B, 2, function(x){sum(x == 0)})
plot(graph.adjacency(nm[which(B[100,] != 0),which(B[100,] != 0)]), layout = layout.kamada.kawai)
