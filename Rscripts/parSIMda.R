library(igraph)
library(parallel)
library(doSNOW)

source("../../../Dropbox/dis-assem.R")
source("./Rscripts/sim-dis-assem.R")

########################
## Get starter webs

init1 <- initialize(Sregion = 1000, Cregion = 0.15, n.initial = 200, Slocal = 45, Sbasal = 5, times = 1000, fr = Fij, frtune = 0)
in1.eqS <- sapply(init1[["DYNAMICS"]], function(x){sum(x[1000, -1] > 0)})
init2 <- initialize(Sregion = 1000, Cregion = 0.15, n.initial = 200, Slocal = 45, Sbasal = 5, times = 1000, fr = Fij, frtune = 0.2)
in2.eqS <- sapply(init2[["DYNAMICS"]], function(x){sum(x[1000, -1] > 0)})
init3 <- initialize(Sregion = 1000, Cregion = 0.15, n.initial = 200, Slocal = 45, Sbasal = 5, times = 1000, fr = Fij, frtune = 1)
in3.eqS <- sapply(init3[["DYNAMICS"]], function(x){sum(x[1000, -1] > 0)})

init4 <- initialize(Sregion = 1000, Cregion = 0.15, n.initial = 200, Slocal = 45, Sbasal = 5, times = 1000, fr = Fbd, frtune = 0)
in4.eqS <- sapply(init4[["DYNAMICS"]], function(x){sum(x[1000, -1] > 0)})
init5 <- initialize(Sregion = 1000, Cregion = 0.15, n.initial = 200, Slocal = 45, Sbasal = 5, times = 1000, fr = Fbd, frtune = 0.2)
in5.eqS <- sapply(init5[["DYNAMICS"]], function(x){sum(x[1000, -1] > 0)})
init6 <- initialize(Sregion = 1000, Cregion = 0.15, n.initial = 200, Slocal = 45, Sbasal = 5, times = 1000, fr = Fbd, frtune = 1)
in6.eqS <- sapply(init6[["DYNAMICS"]], function(x){sum(x[1000, -1] > 0)})

########################
## Deletion simulations
allstrt <- Sys.time()

spdel1 <- sp.deletion(in.dyn = init1[["DYNAMICS"]][in1.eqS> 5], in.n = init1[["INITIAL"]][in1.eqS > 5], n.regional = init1[["REGIONAL"]], funcres = Fij, x = 0, init = "ini1")
spdel2 <- sp.deletion(in.dyn = init2[["DYNAMICS"]][in2.eqS> 5], in.n = init2[["INITIAL"]][in2.eqS> 5], n.regional = init2[["REGIONAL"]], funcres = Fij, x = 0.2, init = "ini2")
spdel3 <- sp.deletion(in.dyn = init3[["DYNAMICS"]][in3.eqS> 5], in.n = init3[["INITIAL"]][in3.eqS> 5], n.regional = init3[["REGIONAL"]], funcres = Fij, x = 1, init = "ini3")
spdel4 <- sp.deletion(in.dyn = init4[["DYNAMICS"]][in4.eqS> 5], in.n = init4[["INITIAL"]][in4.eqS> 5], n.regional = init4[["REGIONAL"]], funcres = Fbd, x = 0, init = "ini4")
spdel5 <- sp.deletion(in.dyn = init5[["DYNAMICS"]][in5.eqS> 5], in.n = init5[["INITIAL"]][in5.eqS> 5], n.regional = init5[["REGIONAL"]], funcres = Fbd, x = 0.2, init = "ini5")
spdel6 <- sp.deletion(in.dyn = init6[["DYNAMICS"]][in6.eqS> 5], in.n = init6[["INITIAL"]][in6.eqS> 5], n.regional = init6[["REGIONAL"]], funcres = Fbd, x = 1, init = "ini6")

allend <- Sys.time()

allend - allstrt

########################
## Invasion simulations
allstrt1 <- Sys.time()

spadd1 <- sp.addition(in.dyn = init1[["DYNAMICS"]][in1.eqS> 5], in.n = init1[["INITIAL"]][in1.eqS> 5], n.regional = init1[["REGIONAL"]], funcres = Fij, x = 0, init = "ini1")
end1 <- Sys.time()
end1 - allstrt1

spadd2 <- sp.addition(in.dyn = init2[["DYNAMICS"]][in2.eqS> 5], in.n = init2[["INITIAL"]][in2.eqS> 5], n.regional = init2[["REGIONAL"]], funcres = Fij, x = 0.2, init = "ini2")
end2 <- Sys.time()
end2 - end1
spadd3 <- sp.addition(in.dyn = init3[["DYNAMICS"]][in3.eqS> 5], in.n = init3[["INITIAL"]][in3.eqS > 5], n.regional = init3[["REGIONAL"]], funcres = Fij, x = 1, init = "ini3")
end3 <- Sys.time()
end3 - end2
spadd4 <- sp.addition(in.dyn = init4[["DYNAMICS"]][in4.eqS> 5], in.n = init4[["INITIAL"]][in4.eqS> 5], n.regional = init4[["REGIONAL"]], funcres = Fbd, x = 0, init = "ini4")
end4 <- Sys.time()
end4 - end3
spadd5 <- sp.addition(in.dyn = init5[["DYNAMICS"]][in5.eqS> 5], in.n = init5[["INITIAL"]][in5.eqS> 5], n.regional = init5[["REGIONAL"]], funcres = Fbd, x = 0.2, init = "ini5")
end5 <- Sys.time()
end5 - end4
spadd6 <- sp.addition(in.dyn = init6[["DYNAMICS"]][in6.eqS> 5], in.n = init6[["INITIAL"]][in6.eqS> 5], n.regional = init6[["REGIONAL"]], funcres = Fbd, x = 1, init = "ini6")

allend1 <- Sys.time()

allend1 - allstrt1

