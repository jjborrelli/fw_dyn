source("./Rscripts/MODELfunctions.R")
source("./Rscripts/fournodeSUBGRAPHS.R")


###
### Dynamics on all 199 four species systems
###

S = 4

K.i <- rep(1, S) # carrying capacity
x.i <- rep(.5, S)
yij <- 6
eij <- 1
q <-  .2

B.i <- runif(S, .5, 1)

meandens <- matrix(0, nrow = length(fournode.am), ncol = 4)
enddens <- matrix(0, nrow = length(fournode.am), ncol = 4)
for(i in 1:length(fournode.am)){
  amat <- fournode.am[[i]]
  r.i <- get.ri(amat)
  r.i2 <- r.i
  
  res <- web_dyn(n.times = 50, amat = amat, B.i = B.i, 
                 r.i = r.i2, K.i = K.i, x.i = x.i, yij = yij, eij = eij, q = q, stochastic = T, ext.thres = 10^-10)
  
  matplot(res, type = "l")
  
  meandens[i, ] <- colMeans(res)
  enddens[i, ] <- res[2000, ]
  
  print(i)
}


fract.p <- apply(enddens, 1, function(x){4 - sum(x == 0)})



### 
### Get qss of fournode
### 

mot4 <- motif_counter(fournode.gr)
rownames(mot4) <- names(fournode)


# run the stability analysis on four node subgraphs
fourN.co <- lapply(fournode.am, conversion)
system.time(
  emat <- eig.analysis(10000, fourN.co)
)
#2.33 min

qss4 <- apply(emat, 2, function(x){sum(x < 0)/length(x)})
names(qss4) <- names(fournode.am)



###
###
###
S = 40

K.i <- rep(1, S) # carrying capacity
x.i <- rep(.5, S)
yij <- 6
eij <- 1
q <-  .2

B.i <- runif(S, .5, 1)

amat <- niche.model(S, .1)
r.i <- get.ri(amat)
r.i2 <- r.i

res <- web_dyn(n.times = 50, amat = amat, B.i = B.i, 
               r.i = r.i2, K.i = K.i, x.i = x.i, yij = yij, eij = eij, q = q, stochastic = T, ext.thres = 10^-10)

matplot(res[1:20,], type = "l")
