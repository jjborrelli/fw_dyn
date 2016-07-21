library(dplyr)
library(lme4)
library(igraph)
library(ggplot2)
library(reshape2)
library(data.table)

#load("D:/jjborrelli/AssemblyDATA/invdel/simulation3.Rdata")

# How many webs in each 

nweb1 <- sum(in1.eqS > 5)
nweb2 <- sum(in2.eqS > 5)
nweb3 <- sum(in3.eqS > 5)
nweb4 <- sum(in4.eqS > 5)
nweb5 <- sum(in5.eqS > 5)
nweb6 <- sum(in6.eqS > 5)

n.eqwebs <- c(nweb1, nweb2, nweb3, nweb4, nweb5, nweb6)

# How many species were deleted from each web

ndel1 <- sum(in1.eqS[in1.eqS > 5])
ndel2 <- sum(in2.eqS[in2.eqS > 5])
ndel3 <- sum(in3.eqS[in3.eqS > 5])
ndel4 <- sum(in4.eqS[in4.eqS > 5])
ndel5 <- sum(in5.eqS[in5.eqS > 5])
ndel6 <- sum(in6.eqS[in6.eqS > 5])

n.deletions <- c(ndel1, ndel2, ndel3, ndel4, ndel5, ndel6)

# How many invasions events for each set

n.invasions <- c(nweb1, nweb2, nweb3, nweb4, nweb5, nweb6)*300
n.invasions


# Is the web connected

iscon1 <- sapply(1:200, function(x){is.connected(graph.adjacency(init1[["INITIAL"]][[x]][init1[["DYNAMICS"]][[x]][1000,-1] > 0,init1[["DYNAMICS"]][[x]][1000,-1] > 0]))})
iscon2 <- sapply(1:200, function(x){is.connected(graph.adjacency(init2[["INITIAL"]][[x]][init2[["DYNAMICS"]][[x]][1000,-1] > 0,init2[["DYNAMICS"]][[x]][1000,-1] > 0]))})
iscon3 <- sapply(1:200, function(x){is.connected(graph.adjacency(init3[["INITIAL"]][[x]][init3[["DYNAMICS"]][[x]][1000,-1] > 0,init3[["DYNAMICS"]][[x]][1000,-1] > 0]))})
iscon4 <- sapply(1:200, function(x){is.connected(graph.adjacency(init4[["INITIAL"]][[x]][init4[["DYNAMICS"]][[x]][1000,-1] > 0,init4[["DYNAMICS"]][[x]][1000,-1] > 0]))})
iscon5 <- sapply(1:200, function(x){is.connected(graph.adjacency(init5[["INITIAL"]][[x]][init5[["DYNAMICS"]][[x]][1000,-1] > 0,init5[["DYNAMICS"]][[x]][1000,-1] > 0]))})
iscon6 <- sapply(1:200, function(x){is.connected(graph.adjacency(init6[["INITIAL"]][[x]][init6[["DYNAMICS"]][[x]][1000,-1] > 0,init6[["DYNAMICS"]][[x]][1000,-1] > 0]))})


# How many species in each web set

l1 <- list(in1.eqS[in1.eqS > 5 & iscon1], in2.eqS[in2.eqS > 5 & iscon2], in3.eqS[in3.eqS > 5 & iscon3], in4.eqS[in4.eqS > 5 & iscon4], in5.eqS[in5.eqS > 5 & iscon5], in6.eqS[in6.eqS > 5 & iscon6])
ggplot(melt(l1), aes(x = factor(L1), y = value)) + geom_boxplot()


########################################################
#### DATA SELECTION

# get deletion data from each connected web

spd1 <- spdel1[spdel1$locweb %in% which(iscon1[in1.eqS > 5]),]
spd2 <- spdel2[spdel2$locweb %in% which(iscon2[in2.eqS > 5]),]
spd3 <- spdel3[spdel3$locweb %in% which(iscon3[in3.eqS > 5]),]
spd4 <- spdel4[spdel4$locweb %in% which(iscon4[in4.eqS > 5]),]
spd5 <- spdel5[spdel5$locweb %in% which(iscon5[in5.eqS > 5]),]
spd6 <- spdel6[spdel6$locweb %in% which(iscon6[in6.eqS > 5]),]

# get invasion data from each connectd web

spa1 <- lapply(spadd1, function(x) x[x$locweb %in% which(iscon1[in1.eqS > 5]),])
spa2 <- lapply(spadd2, function(x) x[x$locweb %in% which(iscon2[in1.eqS > 5]),])
spa3 <- lapply(spadd3, function(x) x[x$locweb %in% which(iscon3[in1.eqS > 5]),])
spa4 <- lapply(spadd4, function(x) x[x$locweb %in% which(iscon4[in1.eqS > 5]),])
spa5 <- lapply(spadd5, function(x) x[x$locweb %in% which(iscon5[in1.eqS > 5]),])
spa6 <- lapply(spadd6, function(x) x[x$locweb %in% which(iscon6[in1.eqS > 5]),])


########################################################
#### DELETION
triplot.del <- function(spd){
  spdP <- spd$persist == 1
  spdNP <- spd$persist != 1
  wp.spd1 <- dplyr::select(spd, L:D, -C)
  mo.spd1 <- dplyr::select(spd, s1:d8)
  dp.spd1 <- dplyr::select(spd, delGen:delOI)
  
  print(barplot(colMeans(wp.spd1[spdNP,], na.rm = T) - colMeans(wp.spd1[spdP,], na.rm = T), las = 2))
  print(barplot(colMeans(mo.spd1[spdNP,], na.rm = T) - colMeans(mo.spd1[spdP,], na.rm = T)))
  print(barplot(colMeans(dp.spd1[spdNP,], na.rm = T) - colMeans(dp.spd1[spdP,], na.rm = T)))
}

par(mfcol = c(3, 6), mar = c(5,4,1.5,1.5))
# deletion
triplot.del(spd1)
triplot.del(spd2)
triplot.del(spd3)
triplot.del(spd4)
triplot.del(spd5)
triplot.del(spd6)

wplot.del <- function(spdlist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spdlist)){
    spd <- spdlist[[i]]
    spdP <- spd$persist == 1
    spdNP <- spd$persist != 1
    wp.spd1 <- dplyr::select(spd, L:D, -C)
    t1 <- melt(wp.spd1[spdNP,], measure.vars = names(wp.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = factor(t1$variable, levels = names(wp.spd1)),
                          L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not persistent")
    t2 <- melt(wp.spd1[spdP,], measure.vars = names(wp.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = factor(t2$variable, levels = names(wp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "persistent")
    df[[i]] <- rbind(df1, df2)
  }
 
  return(rbindlist(df))
}

mplot.del <- function(spdlist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spdlist)){
    spd <- spdlist[[i]]
    spdP <- spd$persist == 1
    spdNP <- spd$persist != 1
    mo.spd1 <- dplyr::select(spd, s1:d8)
    t1 <- melt(mo.spd1[spdNP,], measure.vars = names(mo.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = factor(t1$variable, levels = names(mo.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not persistent")
    t2 <- melt(mo.spd1[spdP,], measure.vars = names(mo.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = factor(t2$variable, levels = names(mo.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "persistent")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}

dplot.del <- function(spdlist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spdlist)){
    spd <- spdlist[[i]]
    spdP <- spd$persist == 1
    spdNP <- spd$persist != 1
    dp.spd1 <- dplyr::select(spd, delGen:delOI)
    t1 <- melt(dp.spd1[spdNP,], measure.vars = names(dp.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = factor(t1$variable, levels = names(dp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not persistent")
    t2 <- melt(dp.spd1[spdP,], measure.vars = names(dp.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = factor(t2$variable, levels = names(dp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "persistent")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}

lplot.del <- function(spalist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spalist)){
    spd <- spalist[[i]]
    spdP <- spd$persist == 1
    spdNP <- spd$persist != 1
    wp.spd1 <- dplyr::select(spd, L)
    t1 <- melt(wp.spd1[spdNP,], measure.vars = names(wp.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = "Not Persistent",
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not invaded")
    t2 <- melt(wp.spd1[spdP,], measure.vars = names(wp.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = "Persistent",
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "invaded")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}


max.se <- function(x){mean(x) + 1.96*(sd(x)/sqrt(length(x)))}
min.se <- function(x){mean(x) - 1.96*(sd(x)/sqrt(length(x)))}

spdl <- list(spd1, spd2, spd3, spd4, spd5, spd6)

ggplot(lplot.del(spdl), aes(x = metric, y = values, fill = stability)) + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(par~FR, scales = "free_y") +
  theme_bw() + scale_fill_discrete(guide = FALSE) +
  xlab(NULL) + ylab("Change in Links Following Species Introduction")
#ggsave("./Figs/lplotDEL.png")

wp1 <- wplot.del(spdl)
wp2 <- wp1[wp1$metric == "meanTP" & FR == "Fij" & par == 0.2 | wp1$metric == "D" & FR == "Fij" & par == 0.2 | wp1$metric == "APL" & FR == "Fij" & par == 0.2,]
wp3 <- wp1[wp1$metric == "meanGen" | wp1$metric == "meanVul" | wp1$metric == "meanTP"| wp1$metric == "D" | wp1$metric == "meanOI",]
ggplot(wp2, aes(x = par, y = values, fill = stability)) + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(metric~FR, scales = "free_y") +
  theme_bw() + #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0)) +
  xlab("Model Parameter") + ylab("Change in Value Following Species Removal")
#ggsave("./Figs/wplotDEL.png")

mp1 <- mplot.del(spdl)
mp2 <- mp1[mp1$metric == "s1" & FR == "Fij" & par == 0.2 | mp1$metric == "s2" & FR == "Fij" & par == 0.2 | mp1$metric == "s4" & FR == "Fij" & par == 0.2 | mp1$metric == "s5" & FR == "Fij" & par == 0.2,]
mp3 <- mp1[mp1$metric == "d1" & FR == "Fij" & par == 0.2 | mp1$metric == "d2" & FR == "Fij" & par == 0.2 | mp1$metric == "d3" & FR == "Fij" & par == 0.2 | mp1$metric == "d4" & FR == "Fij" & par == 0.2 | mp1$metric == "d5" & FR == "Fij" & par == 0.2 | mp1$metric == "d6" & FR == "Fij" & par == 0.2 |  mp1$metric == "d7" & FR == "Fij" & par == 0.2 |  mp1$metric == "d8" & FR == "Fij" & par == 0.2 ,]
ggplot(mp3, aes(x = stability, y = values, fill = stability))  + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(~metric, scales = "free_y") +
  theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0)) +
  xlab("Model Parameter") + ylab("Change in Frequency Following Species Removal")
#ggsave("./Figs/mplotDEL.png")

dpd <- dplot.del(spdl)
dpd[dpd$FR == "Fij" & dpd$par == .2]
ggplot(dpd[dpd$FR == "Fij" & dpd$par == .2], aes(x = par, y = values, fill = stability))  + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(metric~FR, scales = "free_y") +
  theme_bw() +
  xlab("Model Parameter") + ylab("Value For Removed Species")
#ggsave("./Figs/dplotDEL.png")

### Species deletion
## linear mixed effects models with principal components of (1) subgraph composition and (2) web properties

get_pc <- function(spdeldata){
  subg1 <- select(spdeldata, s1:d8)
  webp1 <- select(spdeldata, L:D, -C)
  webp1$APL[is.na(webp1$APL)] <- 0
  delp1 <- select(spdeldata, delGen:delOI)
  
  pcsub <- princomp(subg1)
  pcweb <- princomp(webp1)
  pcdel <- princomp(delp1)
  
  return(list(pcsub, pcweb, pcdel))
}


get_df_glms <- function(spdeldata, pcomplist){
  pcsub <- pcomplist[[1]]
  pcweb <- pcomplist[[2]]
  pcdel <- pcomplist[[3]]

  persisted <- select(spdeldata, persist)
  persisted.bin <- matrix(c(spdeldata$s.fi, (spdeldata$s.in - spdeldata$s.fi)), ncol = 2)

  
  df1 <- data.frame(pc1 = pcsub$scores[,1], pc2 = pcsub$scores[,2], pc3 = pcsub$scores[,3], dID = spdeldata$delID, lweb = spdeldata$locweb)
  df2 <- data.frame(pc1 = pcweb$scores[,1], pc2 = pcweb$scores[,2], pc3 = pcweb$scores[,3], dID = spdeldata$delID, lweb = spdeldata$locweb)
  df3 <- data.frame(pc1 = pcdel$scores[,1], pc2 = pcdel$scores[,2], pc3 = pcdel$scores[,3], dID = spdeldata$delID, lweb = spdeldata$locweb)
  
  return(list(df1, df2, df3))
}

SS <- function(x, y){
  ssr <- sum((fitted(x) - y)^2)
  sst <- sum((y- mean(y))^2)
  return(c(r2 = ssr/sst, SSR = ssr))
}


spd1pc <- get_pc(spd1)
spd2pc <- get_pc(spd2)
spd3pc <- get_pc(spd3)
spd4pc <- get_pc(spd4)
spd5pc <- get_pc(spd5)
spd6pc <- get_pc(spd6)

load1sD <- lapply(list(spd1pc, spd2pc, spd3pc, spd4pc, spd5pc, spd6pc), function(q) lapply(q, function(x) loadings(x)[,1]))
subloadingD <- round(sapply(load1sD, "[[",1), 3)
webloadingD <- round(sapply(load1sD, "[[",2), 3)
delloadingD <- round(sapply(load1sD, "[[",3), 3)
#write.csv(subloadingD, "D:/jjborrelli/AssemblyDATA/invdel/subloadD.csv")
#write.csv(webloadingD, "D:/jjborrelli/AssemblyDATA/invdel/webloadD.csv")
#write.csv(delloadingD, "D:/jjborrelli/AssemblyDATA/invdel/delloadD.csv")


spd1gm <- get_df_glms(spd1, spd1pc)
spd2gm <- get_df_glms(spd2, spd2pc)
spd3gm <- get_df_glms(spd3, spd3pc)
spd4gm <- get_df_glms(spd4, spd4pc)
spd5gm <- get_df_glms(spd5, spd5pc)
spd6gm <- get_df_glms(spd6, spd6pc)

allspD <- list(spd1, spd2, spd3, spd4, spd5, spd6)
allgmdf <- list(spd1gm, spd2gm, spd3gm, spd4gm, spd5gm, spd6gm)

allgmtab <- list()
for(i in 1:6){
  # glm web prop
  persisted.bin <- cbind(allspD[[i]]$s.fi, (allspD[[i]]$s.in - allspD[[i]]$s.fi))
  gm1 <- glmer(persisted.bin ~ pc1 + (1 | lweb), data = allgmdf[[i]][[1]], family = "binomial")
  ss <- getME(gm1,c("theta","fixef"))
  gm1.1 <- update(gm1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  gm.sum <- summary(gm1.1)
  gm.ss <- SS(gm1, allspD[[i]]$persist)
  subg.gm <- c(gm.sum$coefficients[2,], gm.sum$AICtab, gm.ss)
  
  gm2 <- glmer(persisted.bin ~ pc1 + (1 | lweb), data = allgmdf[[i]][[2]], family = "binomial")
  gm.sum <- summary(gm2)
  gm.ss <- SS(gm2, allspD[[i]]$persist)
  webp.gm <- c(gm.sum$coefficients[2,], gm.sum$AICtab, gm.ss)
  
  gm3 <- glmer(persisted.bin ~ pc1 + (1 | lweb), data = allgmdf[[i]][[3]], family = "binomial")
  gm.sum <- summary(gm3)
  gm.ss <- SS(gm3, allspD[[i]]$persist)
  delp.gm <- c(gm.sum$coefficients[2,], gm.sum$AICtab, gm.ss)
  
  allgmtab[[i]] <- rbind(subg.gm, webp.gm, delp.gm)
  print(i)
}

allgmtab

########################################################
#### INVASION
triplot.inv <- function(spa){
  
  #par(mfcol = c(2,1))
  barplot(colMeans(spa[[1]][spa[[2]]$I == TRUE & spa[[2]]$dN == 0])[2:26]-colMeans(spa[[1]][spa[[2]]$I == TRUE & spa[[2]]$dN < 0])[2:26], las = 2)
  
  barplot(colMeans(spa[[1]][spa[[2]]$I == TRUE & spa[[2]]$dN < 0])[2:26]-colMeans(spa[[1]][spa[[2]]$I == FALSE & spa[[2]]$dN < 0])[2:26], las = 2)
}

par(mfcol = c(2, 3), mar = c(5,4,1.5,1.5))
# invasion
triplot.inv(spa1)
triplot.inv(spa2)
triplot.inv(spa3)
triplot.inv(spa4)
triplot.inv(spa5)
triplot.inv(spa6)

barplot(colMeans(spa1[[1]][spa1[[2]]$I == TRUE])[2:26] - colMeans(spa1[[1]][spa1[[2]]$I == FALSE])[2:26])
barplot(colMeans(spa2[[1]][spa2[[2]]$I == TRUE])[2:26] - colMeans(spa2[[1]][spa2[[2]]$I == FALSE])[2:26])
barplot(colMeans(spa3[[1]][spa3[[2]]$I == TRUE])[2:26] - colMeans(spa3[[1]][spa3[[2]]$I == FALSE])[2:26])
barplot(colMeans(spa4[[1]][spa4[[2]]$I == TRUE])[2:26] - colMeans(spa4[[1]][spa4[[2]]$I == FALSE])[2:26])
barplot(colMeans(spa5[[1]][spa5[[2]]$I == TRUE])[2:26] - colMeans(spa5[[1]][spa5[[2]]$I == FALSE])[2:26])
barplot(colMeans(spa6[[1]][spa6[[2]]$I == TRUE])[2:26] - colMeans(spa6[[1]][spa6[[2]]$I == FALSE])[2:26])

wplot.inv <- function(spalist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spalist)){
    spd <- spalist[[i]]
    spdP <- spd[[2]]$I == TRUE
    spdNP <- spd[[2]]$I == FALSE
    wp.spd1 <- dplyr::select(spd[[1]], L:D, -C)
    t1 <- melt(wp.spd1[spdNP,], measure.vars = names(wp.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = factor(t1$variable, levels = names(wp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not invaded")
    t2 <- melt(wp.spd1[spdP,], measure.vars = names(wp.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = factor(t2$variable, levels = names(wp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "invaded")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}

mplot.inv <- function(spalist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spalist)){
    spd <- spalist[[i]]
    spdP <- spd[[2]]$I == TRUE
    spdNP <- spd[[2]]$I == FALSE
    mo.spd1 <- dplyr::select(spd[[1]], s1:d8)
    t1 <- melt(mo.spd1[spdNP,], measure.vars = names(mo.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = factor(t1$variable, levels = names(mo.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not invaded")
    t2 <- melt(mo.spd1[spdP,], measure.vars = names(mo.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = factor(t2$variable, levels = names(mo.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "invaded")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}

dplot.inv <- function(spdlist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spdlist)){
    spd <- spdlist[[i]]
    spdP <- spd[[2]]$I == TRUE
    spdNP <- spd[[2]]$I == FALSE
    dp.spd1 <- dplyr::select(spd[[2]], invGen:dN)
    dp.spd1$dN <- dp.spd1$dN * -1
    t1 <- melt(dp.spd1[spdNP,], measure.vars = names(dp.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = factor(t1$variable, levels = names(dp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not invaded")
    t2 <- melt(dp.spd1[spdP,], measure.vars = names(dp.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = factor(t2$variable, levels = names(dp.spd1)),
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "invaded")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}

lplot.inv <- function(spalist){
  dyn <- expand.grid(c("Fij", "Fbd"), c("0", "0.2", "1"))
  df <- list()
  for(i in 1:length(spalist)){
    spd <- spalist[[i]]
    spdP <- spd[[2]]$I == TRUE
    spdNP <- spd[[2]]$I == FALSE
    wp.spd1 <- dplyr::select(spd[[1]], L)
    t1 <- melt(wp.spd1[spdNP,], measure.vars = names(wp.spd1[spdNP,])) 
    df1 <- data.frame(values = t1$value, metric = "Not Invaded",
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "not invaded")
    t2 <- melt(wp.spd1[spdP,], measure.vars = names(wp.spd1[spdP,]))
    df2 <- data.frame(values = t2$value, metric = "Invaded",
                      L1 = i, FR = dyn[i,1], par = dyn[i,2], stability = "invaded")
    df[[i]] <- rbind(df1, df2)
  }
  
  return(rbindlist(df))
}

spal <- list(spa1, spa2, spa3, spa4, spa5, spa6)

ggplot(lplot.inv(spal), aes(x = metric, y = values, fill = stability)) + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(par~FR) +
  theme_bw() + scale_fill_discrete(guide = FALSE) +
  xlab(NULL) + ylab("Change in Links Following Species Introduction")

wp1 <- wplot.inv(spal)
wp2 <- wp1[wp1$metric == "meanTP" & wp1$FR == "Fij" & wp1$par == 0.2 | wp1$metric == "D" & wp1$FR == "Fij" & wp1$par == 0.2,]
ggplot(wp2, aes(x = stability, y = values, fill = stability)) + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(~metric, scales = "free_y") +
  theme_bw() + #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0)) +
  xlab("Model Parameter") + ylab("Change in Value Following Species Introduction")
#ggsave("./Figs/wplotINV.png")

mp1 <- mplot.inv(spal)
mp2 <- mp1[mp1$metric == "s1" & mp1$FR == "Fij" & mp1$par == 0.2 | mp1$metric == "s2" & mp1$FR == "Fij" & mp1$par == 0.2 | mp1$metric == "s4" & mp1$FR == "Fij" & mp1$par == 0.2 | mp1$metric == "s5"& mp1$FR == "Fij" & mp1$par == 0.2,]
mp3 <- mp1[mp1$metric == "d1" | mp1$metric == "d2" | mp1$metric == "d3"  | mp1$metric == "d4" | mp1$metric == "d5" | mp1$metric == "d6" |  mp1$metric == "d7" |  mp1$metric == "d8",]
ggplot(mp3, aes(x = par, y = values, fill = stability)) + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(metric~FR, scales = "free_y") +
  theme_bw() + #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0)) +
  xlab("Model Parameter") + ylab("Change in Frequency Following Species Introduction")
#ggsave("./Figs/mplotINV.png")

dp1 <- dplot.inv(spal)
ggplot(dp1[dp1$FR == "Fij" & dp1$par == 0.2,], aes(x = par, y = values, fill = stability)) + 
  geom_bar(stat = "summary", position = "dodge") + 
  stat_summary(fun.y = "mean", fun.ymin = min.se, fun.ymax = max.se, geom = "errorbar", position = "dodge") + 
  facet_grid(metric~FR, scales = "free_y") +
  theme_bw() +
  xlab("Model Parameter") + ylab("Value for Introduced Species")
#ggsave("./Figs/dplotINV.png")





##################################################################################
##################################################################################
##################################################################################


get_pc_inv <- function(spinvdata){
  subg1 <- select(spinvdata[[1]], s1:d8)
  webp1 <- select(spinvdata[[1]], L:D, -C)
  webp1$APL[is.na(webp1$APL)] <- 0
  invp1 <- select(spinvdata[[2]], invGen:invOI)
  
  pcsub <- princomp(subg1)
  pcweb <- princomp(webp1)
  pcinv <- princomp(invp1)
  
  return(list(pcsub, pcweb, pcinv))
}


get_df_glmsI <- function(spdeldata, pcomplist){
  pcsub <- pcomplist[[1]]
  pcweb <- pcomplist[[2]]
  pcdel <- pcomplist[[3]]
  
  
  df1 <- data.frame(pc1 = pcsub$scores[,1], pc2 = pcsub$scores[,2], pc3 = pcsub$scores[,3], lweb = spdeldata[[1]]$locweb)
  df2 <- data.frame(pc1 = pcweb$scores[,1], pc2 = pcweb$scores[,2], pc3 = pcweb$scores[,3], lweb = spdeldata[[1]]$locweb)
  df3 <- data.frame(pc1 = pcdel$scores[,1], pc2 = pcdel$scores[,2], pc3 = pcdel$scores[,3], lweb = spdeldata[[1]]$locweb)
  
  return(list(df1, df2, df3))
}



spa1pc <- get_pc_inv(spa1)
spa2pc <- get_pc_inv(spa2)
spa3pc <- get_pc_inv(spa3)
spa4pc <- get_pc_inv(spa4)
spa5pc <- get_pc_inv(spa5)
spa6pc <- get_pc_inv(spa6)

load1s <- lapply(list(spa1pc, spa2pc, spa3pc, spa4pc, spa5pc, spa6pc), function(q) lapply(q, function(x) loadings(x)[,1]))
subloadingI <- round(sapply(load1s, "[[",1), 3)
webloadingI <- round(sapply(load1s, "[[",2), 3)
invloadingI <- round(sapply(load1s, "[[",3), 3)
#write.csv(subloadingI, "D:/jjborrelli/AssemblyDATA/invdel/subloadI.csv")
#write.csv(webloadingI, "D:/jjborrelli/AssemblyDATA/invdel/webloadI.csv")
#write.csv(invloadingI, "D:/jjborrelli/AssemblyDATA/invdel/invloadI.csv")


spa1gm <- get_df_glmsI(spa1, spa1pc)
spa2gm <- get_df_glmsI(spa2, spa2pc)
spa3gm <- get_df_glmsI(spa3, spa3pc)
spa4gm <- get_df_glmsI(spa4, spa4pc)
spa5gm <- get_df_glmsI(spa5, spa5pc)
spa6gm <- get_df_glmsI(spa6, spa6pc)

allspI <- list(spa1[[2]], spa2[[2]], spa3[[2]], spa4[[2]], spa5[[2]], spa6[[2]])
allgmdfI <- list(spa1gm, spa2gm, spa3gm, spa4gm, spa5gm, spa6gm)

allgmtabI <- list()
for(i in 1:6){
  # glm web prop
  invading <- allspI[[i]]$I
  gm1 <- glmer(invading ~ pc1 + (1 | lweb), data = allgmdfI[[i]][[1]], family = "binomial")
  gm1.1 <- update(gm1,start=ss,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
  gm.sum <- summary(gm1.1)
  gm.ss <- SS(gm1, invading)
  subg.gm <- c(gm.sum$coefficients[2,], gm.sum$AICtab, gm.ss)
  
  gm2 <- glmer(invading ~ pc1 + (1 | lweb), data = allgmdfI[[i]][[2]], family = "binomial")
  gm.sum <- summary(gm2)
  gm.ss <- SS(gm2, invading)
  webp.gm <- c(gm.sum$coefficients[2,], gm.sum$AICtab, gm.ss)
  
  gm3 <- glmer(invading ~ pc1 + (1 | lweb), data = allgmdfI[[i]][[3]], family = "binomial")
  gm.sum <- summary(gm3)
  gm.ss <- SS(gm3, invading)
  invp.gm <- c(gm.sum$coefficients[2,], gm.sum$AICtab, gm.ss)
  
  allgmtabI[[i]] <- rbind(subg.gm, webp.gm, invp.gm)
  print(i)
}

allgmtabI


frs <- c("Fij", "Fij", "Fij", "Fbd", "Fbd", "Fbd")
xpars <- c(0, .2, 1, 0, .2, 1)
allgmtabI.1 <- allgmtabI
for(i in 1:6){
  allgmtabI.1[[i]] <- data.frame(allgmtabI[[i]], FR = frs[i], Par = xpars[i], mod = rownames(allgmtabI[[i]]))
  allgmtab[[i]] <- data.frame(allgmtab[[i]], FR = frs[i], Par = xpars[i], mod = rownames(allgmtab[[i]]))
}

#write.csv(rbindlist(allgmtabI.1), "D:/jjborrelli/AssemblyDATA/invdel/invSLOPES.csv")
#write.csv(rbindlist(allgmtab), "D:/jjborrelli/AssemblyDATA/invdel/delSLOPES.csv")


la <- lapply(list(spd1pc, spd2pc, spd3pc, spd4pc, spd5pc, spd6pc), function(q) lapply(q, function(x) (x$sdev^2 / sum(x$sdev^2))[1]))
lapply(la, "[[", 1)
lapply(la, "[[", 2)
lapply(la, "[[", 3)

lb <- lapply(list(spa1pc, spa2pc, spa3pc, spa4pc, spa5pc, spa6pc), function(q) lapply(q, function(x) (x$sdev^2 / sum(x$sdev^2))[1]))
range(unlist(lapply(lb, "[[", 1)))
range(unlist(lapply(lb, "[[", 2)))
range(unlist(lapply(lb, "[[", 3)))
#####################################################################
#####################################################################
##### Time to First Extinction
library(rend)
library(NetIndices)

firstEXT <- c()
propDIFF <- matrix(ncol = 26, nrow = 100)
for(i in 1:100){
  cond <- FALSE
  while(!cond){
    nm <- niche.model(50, .15)
    cond <- sum(colSums(nm) == 0) == 5 & is.connected(graph.adjacency(nm))
  }
  
  dyn <- CRsimulator(nm, t = 1:100)
  firstEXT[i] <- min(which(apply(dyn[,-1], 1, function(x) sum(x == 0)) > 1))
  extinctSPP <- which(dyn[firstEXT[i],-1] == 0)
  if(length(extinctSPP) > 1){extinctSPP <- sample(extinctSPP, 1)}
  propDIFF[i,] <- unlist(wprop(nm) - wprop(nm[-extinctSPP, -extinctSPP]))
  print(i)
}