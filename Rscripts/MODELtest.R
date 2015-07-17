


A <- matrix(c(0,1,0,-1,0,1,0,-1,0), nrow = 3) * matrix(runif(9, 0, 1), nrow = 3)
#A <- matrix(c(0,1,1,-1,0,0,-1,0,0), nrow = 3) * matrix(runif(9, 0, 1), nrow = 3)
h <- 1
K <- 1000
t = 50

N <- c(100, 15, 20)

Nt <- matrix(nrow = t, ncol = nrow(A))
for(i in 1:t){
  N <- N + c(1, 1, 1) * N * LV.fr(N, A, h) - (.2 * N) + (c(1,0,0) * N * (1 - (N/K))) 
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

out<-ode(y=state,times=times,func=lvmodel,parms=parameters)

matplot(out[,2:4])
