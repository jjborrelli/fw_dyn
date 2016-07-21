library(rend)
library(ggplot2)
library(igraph)
library(NetIndices)

load("C:/Users/jjborrelli/Desktop/testR.Rdata")
#save.image("C:/Users/jjborrelli/Desktop/testR.Rdata")

niche_model <- function(S, C){
  conn <- FALSE
  while(!conn){
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
    
    conn <- is.connected(graph.adjacency(new.mat))
  }
  df <- data.frame(niche, r, ci)
  return(list(new.mat, df))
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



nm1 <- list()
for(i in 1:10000){nm1[[i]] <- niche_model(S = 50, C = 0.15)}

nm2 <- lapply(nm1, "[[", 1)
nmdf <- lapply(nm1, "[[", 2)

itro <- sapply(nm2, function(x){max(TrophInd(x)$TL)})

cr1 <- list()
for(i in 1:10000){cr1[[i]] <- CRsimulator(nm2[[i]], t = 1:400);print(i)}

per <- sapply(cr1, function(x){sum(tail(x, 1)[-1] > 0)})

mc1 <- motif_counter(lapply(nm2, graph.adjacency))
hist(sapply(nm2, sum))
ibio <- lapply(cr1, function(x){head(x, 1)[-1]})
fbio <- lapply(cr1, function(x){tail(x, 1)[-1][tail(x, 1)[-1] > 0]})
fext <- lapply(1:10000, function(x){w1 <- which(tail(cr1[[x]], 1)[-1] > 0); return(nm2[[x]][w1,w1])})
fg1 <- lapply(fext, graph.adjacency)
fint <- lapply(1:10000, function(x){Fij(fbio[[x]], fext[[x]], .5, .2)})
ftro <- lapply(fext, function(x){TrophInd(x)$TL})
fmo <- motif_counter(fg1)


mc2 <- t(apply(mc1, 1, function(x){x/sum(x)}))
fmo2 <- t(apply(fmo, 1, function(x){x/sum(x)}))
conn <- sapply(nm2, function(x){sum(x)/(50*49)})
s1 <- sample(1:10000, 100)
pred1 <- rowSums(mc2[, c(1,4,5)])
pred2 <- rowSums(fmo2[, c(1,4,5)])
fit1 <- lm(per~conn + pred1)
summary(fit1)
abline(fit1)

ggplot(data.frame(per, itro), aes(x = itro, y = per)) + geom_point() + geom_smooth(method = "lm")
ggplot(data.frame(ftro1 = sapply(ftro, max), itro), aes(x = itro, y = ftro1)) + geom_point() + geom_abline(slope = 1)



mopro <- function(m, bio, net){
  test <- Fij(bio, net, .5, .2)
  
  mo <- m
  for(i in seq(.05, .95, .05)){
    test2 <- test
    test2[test2 < quantile(test[test > 0], probs = i)] <- 0
    mo <- rbind(motif_counter(list(graph.adjacency(ceiling(test2)))), mo)
  }
  
  return(cbind(mo, c(seq(.05, .95, .05),1)))
}


mopro2 <- function(bio, net){
  test <- Fij(bio, net, .5, .2)
  sq1 <- seq(.1, 1, .1)
  sq2 <- seq(0, .9, .1)
  
  
  test2 <- matrix(0, nrow = nrow(test), ncol = ncol(test))
  mo <- matrix(nrow = 10, ncol = 13)
  for(i in 1:10){
    
    co1 <- test > quantile(test[test > 0], probs = sq2[i])
    co2 <- test <= quantile(test[test > 0], probs = sq1[i])
    test2[co1 & co2] <- 1
    mo[i,] <- unlist(motif_counter(list(graph.adjacency(test2))))
  }
  return(mo)
}


qmo <- lapply(1:2, function(x){mopro2(ibio[[x]], nm2[[x]])})
qmo2 <- lapply(1:10000, function(x){mopro2(fbio[[x]], fext[[x]])})

wdub <- sapply(qmo, function(x){sum(rowSums(x[,6:13]) == 0)})
wdub <- sapply(qmo2, function(x){sum(rowSums(x[,6:13]) != 0)})


plot(itro~sapply(momean, "[[",1))

##
mots <- list()
for(i in 1:10000){
  nlist <- lapply(1:1000, function(x){niche_model(per[i], conn[i])[[1]]})
  glist <- lapply(nlist, graph.adjacency)
  mots[[i]] <- motif_counter(glist)
}


momean <- lapply(mots, function(x){colMeans(x)})
mosd <- lapply(mots, function(x){apply(x, 2, sd)})

z1 <- t(sapply(1:10000, function(x){(unlist(fmo[x,]) - momean[[x]])/mosd[[x]]}))
boxplot(z1)
