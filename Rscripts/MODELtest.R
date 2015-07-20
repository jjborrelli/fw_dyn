LV.fr <- function(N, A, h){
  aN <- A
  aN[A < 0] <-(A * as.vector(N))[A < 0]
  aN[A > 0] <- (t(A) * as.vector(N))[t(A) > 0]
  return(t(sapply(1:nrow(A), function(x){aN[x,]/(1 + h * rowSums(abs(aN)))[x]})))
}
# something still wrong -- doesnt add up

A <- matrix(c(0,1,0,-1,0,1,0,-1,0), nrow = 3) * matrix(runif(9, 0, 1), nrow = 3)
#A <- matrix(c(0,1,1,-1,0,0,-1,0,0), nrow = 3) * matrix(runif(9, 0, 1), nrow = 3)
h <- 1
K <- 1000
t = 50

N <- c(100, 15, 20)

Nt <- matrix(N, nrow = t+1, ncol = nrow(A), byrow = T)
for(i in 2:(t+1)){
  N <- N + c(1, .1, .1) * LV.fr(N, A, h) %*% N - (.2 * N) + (c(1,0,0) * N * (1 - (N/K))) 
  N[N<1] <- 0
  Nt[i,] <- N
}
matplot(Nt, typ = "l")


library(deSolve)

parameters<-c(r=1,alpha=abs(A[4]), alpha.1=A[2],alpha2=abs(A[8]), alpha2.1=A[6],e=.1,e2=.1,q=.2,K=K)

state<-c(N=100,N1=15,N2=20)

lvmodel<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    #rate of change
    dN<-r*N*(1-N/K)-((alpha.1*N)/(1+alpha.1*N))*N1 - q*N
    dN1<-e*((alpha*N)/(1+alpha*N))*N1-((alpha2.1*N1)/(1+alpha2.1*N1))*N2-q*N1
    dN2<-e2*((alpha2*N1)/(1+alpha2*N1))*N2-q*N2
    
    #return rate of change
    list(c(dN,dN1,dN2))
  })
}


times<-seq(0,50,by=1)

out<-ode(y=state,times=times,func=lvmodel,parms=parameters, method = "adams")

matplot(out[,2:4], typ = "l")


parameters<-c(r=1,K=K)

state<-c(N=100)

ss.test <- function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    #rate of change
    dN<-r*N*(1-N/K)
    #return rate of change
    list(c(dN))
  })
}

times<-seq(1,20,by=1)

out<-ode(y=state,times=times,func=ss.test,parms=parameters)
