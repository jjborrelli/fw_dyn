library(igraph)
library(NetIndices)
library(reshape2)
library(rend)
library(ggplot2)
library(parallel)
library(doSNOW)
filepath.sink <- "D:/jjborrelli/AssemblyDATA/Sim2/"
# sim1 Nregion = 1000, Cregion = .15
# sim1 Nregion = 1000, Cregion = .3; regional2

regional2 <- niche.model(1000, .3)
colnames(regional2) <- as.character(1:nrow(regional2))
rownames(regional2) <- as.character(1:nrow(regional2))

locals <- lapply(1:100, function(x){initialNiche(regional2, nbasal = 5)})


strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)

RESULT2 <- foreach(i = 1:100) %dopar% {
  source("./fw_dyn/Rscripts/assembleN.R")
  inv.events <- 200
  sink(file = paste(filepath.sink, "simdata", i,".txt", collapse = ""))
  combuild <- assembly(regional2, locals[[i]], n.inv = inv.events)
  
  pa <- matrix(nrow = length(combuild[[1]]), ncol = length(unique(unlist(combuild[[1]]))))
  colnames(pa) <- as.character(unique(unlist(combuild[[1]])))
  for(n in 1:nrow(pa)){pa[n,] <- as.numeric(colnames(pa) %in% combuild[[1]][[n]])}
  combuild$pa <- pa
  
  print(combuild)

  sink()
  
  return(combuild)
}

stopCluster(cl)
ends <- Sys.time()
ends - strt
run2time <- ends - strt

#run1time <- ends - strt

#tryout <- RESULT # run of simulation for 7 webs
comms <- lapply(RESULT, "[[", 1)
invwebp <- lapply(RESULT, "[[", 2)
eqwebp <- lapply(RESULT, "[[", 3)
invaders <- lapply(RESULT, "[[", 4)
occupancy <- lapply(RESULT, "[[", 5)


test <- invaders[[1]]

low <- function(y){mean(y) - 1.96 * (sd(y)/sqrt(length(y)))}
upp <- function(y){mean(y) + 1.96 * (sd(y)/sqrt(length(y)))}
ggplot(melt(do.call(rbind, invaders)[-c(1,3)], id.vars = "invEST"), aes(x = variable, y = value, fill = invEST)) + geom_bar(stat = "summary", position = "dodge") + stat_summary(fun.y = "mean", fun.ymin = "low", fun.ymax = "upp", geom = "errorbar", position = "dodge") 
ggplot(melt(do.call(rbind, invaders)[-c(1,3)], id.vars = "invEST"), aes(x = invEST, y = value, fill = invEST)) + geom_bar(stat = "summary", position = "dodge") + stat_summary(fun.y = "mean", fun.ymin = "low", fun.ymax = "upp", geom = "errorbar", position = "dodge") + facet_wrap(~variable, scales = "free_y") 

test2 <- invwebp[[1]]
test3 <- cbind(test2, SCENARIO = paste(rep(test$invEST, each = 26),rep(test$spGain < 0, each = 26))) 
ggplot(test3, aes(x = SCENARIO, y = value, fill = SCENARIO)) + geom_bar(stat = "summary", position = "dodge") + stat_summary(fun.y = "mean", fun.ymin = "low", fun.ymax = "upp", geom = "errorbar", position = "dodge") + facet_wrap(~variable, scales = "free_y") 


InE <- lapply(wps, function(x){which(x$inv & x$delS > 0)})
IE <- lapply(wps, function(x){which(x$inv & x$delS < 0)})
nIE <- lapply(wps, function(x){which(!x$inv & x$delS > 0)})
nInE <- lapply(wps, function(x){which(!x$inv & x$delS < 0)})

s.ine <- do.call(rbind, lapply(1:100, function(x){subg[[x]][InE[[x]],]}))
s.ine[is.nan(s.ine)] <- 0


deltaS <- lapply(lapply(1:100, function(x) combuild[[x]][[1]]), function(y){sapply(y, length)})
dS <- t(do.call(rbind, deltaS))
mean(dS[100,])
sd(dS[100,])

moch <- lapply(1:100, function(x){motifChanges(combuild[[x]][[1]], regional)})
mots <- lapply(1:100, function(x){absMOTIF(combuild[[x]][[1]], regional)})
motifdiff <- lapply(mots, function(x){
  newmat <- matrix(0, nrow = nrow(x), ncol = 13)
  for(i in 2:nrow(x)){
    newmat[i,] <- unlist(x[i, ]) - unlist(x[i-1,])
  }
  return(newmat)
})

comms <- lapply(1:100, function(x) combuild[[x]][[1]])
invader <- lapply(1:100, function(x){
  inv <- c()
  delS <- c()
  delL <- c()
  apl <- c()
  mot <- matrix(0, nrow = 99, ncol = 13)
  colnames(mot) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
  for(i in 2:100){
    inv[i-1] <- sppintro[[x]][i-1] %in% comms[[x]][[i]]
    delS[i-1] <- length(comms[[x]][[i]]) - length(comms[[x]][[i-1]])  
    rnet <- regional[c(sppintro[[x]][[i-1]], comms[[x]][[i-1]]), c(sppintro[[x]][[i-1]], comms[[x]][[i-1]])]
    rnet2 <- regional[comms[[x]][[i-1]], comms[[x]][[i-1]]]
    delL[i-1] <- sum(rnet) - sum(rnet2)
    apl[i-1] <- average.path.length(graph.adjacency(rnet)) - average.path.length(graph.adjacency(rnet2))
    mot[i-1,] <- unlist(motif_counter(list(graph.adjacency(rnet)))) - unlist(motif_counter(list(graph.adjacency(rnet2)))) 
  }
  return(data.frame(inv, delS, delL, apl, mot))
})


allinv <- do.call(rbind, invader)

inv.ext <- allinv[allinv$inv & allinv$delS < 0,]
inv.next <- allinv[allinv$inv & allinv$delS > 0,]
ninv.ext <- allinv[!allinv$inv & allinv$delS < 0,]
ninv.next <- allinv[!allinv$inv & allinv$delS == 0,]

outcomes <- list(IE = inv.ext, InE = inv.next, nIE = ninv.ext, nInE = ninv.next)

websub <- lapply(outcomes, function(x){(x[,5:17])})
wsub <- melt(websub)
ymin <- function(y){mean(y) - 1.96*(sd(y)/length(y))}
ymax <- function(y){mean(y) + 1.96*(sd(y)/length(y))}

a1 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), mean)
a2 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymax)
a3 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymin)

ws.a <- data.frame(a1, ymin = a3$x, ymax = a2$x)

ggplot(ws.a, aes(x = Group.2, y = x)) + geom_point() + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + facet_wrap(~Group.1, scales = "free_y")


wp <- melt(webp.out)
ggplot(wp, aes(x = L1, y = value)) + geom_point(position = "jitter", alpha = .1) + geom_errorbar() + facet_wrap(~variable, scales = "free_y")

mu.out <- lapply(outcomes, function(x){sapply(x, function(y){c(m = mean(y), up = mean(y) + sd(y)/sqrt(length(y)), low = mean(y) - sd(y)/sqrt(length(y)))})})
d1 <- melt(mu.out)
ggplot(d1, aes(x = L1, y = value, col = Var1)) + geom_point() + facet_wrap(~Var2, scales = "free_y")

# means and errors of web props for four outcomes
ggplot(d1[d1$Var2 == "delS" | d1$Var2 == "delL" | d1$Var2 == "apl", ], aes(x = L1, y = value, col = Var1)) + geom_point(size= 5) + facet_wrap(~Var2, scales = "free_y")

# means and error of subgraphs for four outcomes
ggplot(d1[d1$Var2 != "delS" & d1$Var2 != "delL" & d1$Var2 != "apl" & d1$Var2 != "inv", ], aes(x = L1, y = value, col = Var1)) + geom_point(size= 5) + facet_wrap(~Var2, scales = "free_y")



invaderMOT <- lapply(1:100, function(x){
  mot <- matrix(0, nrow = 99, ncol = 13)
  colnames(mot) <- c("s1", "s2", "s3", "s4", "s5", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8")
  for(i in 2:100){
    inv[i-1] <- sppintro[[x]][i-1] %in% comms[[x]][[i]]
    rnet <- graph.adjacency(regional[c(sppintro[[x]][[i-1]], comms[[x]][[i-1]]), c(sppintro[[x]][[i-1]], comms[[x]][[i-1]])])
    rnet2 <- graph.adjacency(regional[comms[[x]][[i-1]], comms[[x]][[i-1]]])
    rn2mot <- motif_counter(list(rnet2))
    rmot <- matrix(nrow = 200, ncol = 13)
    for(j in 1:200){
      spp <- c(sample((1:1000)[-comms[[x]][[i-1]]]), comms[[x]][[i-1]])
      rand <- graph.adjacency(regional[spp,spp])
      rmot[j,] <- unlist(motif_counter(list(rand))) - unlist(rn2mot)
    }
    diffs <- unlist(motif_counter(list(rnet))) - unlist(rn2mot)
    mot[i-1,] <- (diffs - colMeans(rmot))/apply(rmot, 2, sd)
  }
  return(data.frame(mot))
})

allinv2 <- do.call(rbind, invaderMOT)

inv.ext <- allinv2[allinv$inv & allinv$delS < 0,]
inv.next <- allinv2[allinv$inv & allinv$delS > 0,]
ninv.ext <- allinv2[!allinv$inv & allinv$delS < 0,]
ninv.next <- allinv2[!allinv$inv & allinv$delS == 0,]

outcomes <- list(IE = inv.ext, InE = inv.next, nIE = ninv.ext, nInE = ninv.next)

websub <- lapply(outcomes, function(x){(x[,5:17])})
wsub <- melt(websub)
ymin <- function(y){mean(y) - 1.96*(sd(y)/length(y))}
ymax <- function(y){mean(y) + 1.96*(sd(y)/length(y))}

a1 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), mean)
a2 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymax)
a3 <- aggregate(wsub$value, list(wsub$variable, wsub$L1), ymin)

ws.a <- data.frame(a1, ymin = a3$x, ymax = a2$x)

ggplot(ws.a, aes(x = Group.2, y = x)) + geom_point() + geom_errorbar(aes(ymin = ymin, ymax = ymax)) + facet_wrap(~Group.1, scales = "free_y")