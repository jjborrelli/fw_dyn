library(ggplot2)
library(parallel)
library(doSNOW)

source("../../../Dropbox/dis-assem.R")
filepath.sink <- "D:/jjborrelli/AssemblyDATA/invdel/"


sp.addition <- function(in.dyn, in.n, n.regional, funcres = Fij, x = 0.2, init){
  cl <- makeCluster(detectCores() - 1)
  registerDoSNOW(cl)
  
  INVASION <- foreach(i = 1:length(in.n)) %dopar% {
    source("../../../Dropbox/dis-assem.R")
    sink(file = paste0(filepath.sink, "invasion/", init, "/", "invading-web", i,".txt", collapse = ""))
    
    # SIMULATION
    numinvs = 300
    
    invaded <- lapply(1:numinvs, function(x) invading2(in.dyn[[i]], in.n[[i]], n.regional, FR = funcres, frtune = x))
    
    webi.props <- do.call(rbind, lapply(invaded, "[[", 1))
    inv.props <- do.call(rbind, lapply(invaded, "[[", 2))
    
    webi.props <- data.frame(webi.props, inum = 1:numinvs, locweb = i)
    inv.props <- data.frame(inv.props, inum = 1:numinvs, locweb = i)
    
    
    print(list(webi.props, inv.props))
    
    sink()
    
    return(list(webi.props, inv.props))
  }
  
  stopCluster(cl)
  
  wpdf <- rbindlist(lapply(INVASION, "[[", 1))
  ipdf <- rbindlist(lapply(INVASION, "[[", 2))
  
  return(list(wpdf = wpdf, ipdf = ipdf))
}



##
##
## Disassembly
##




sp.deletion <- function(in.dyn, in.n, n.regional, funcres = Fij, x = 0.2, init){
  cl <- makeCluster(detectCores() - 1)
  registerDoSNOW(cl)
  
  DELETION <- foreach(i = 1:length(in.n)) %dopar% {
    source("../../../Dropbox/dis-assem.R")
    sink(file = paste0(filepath.sink, "deletion/", init, "/", "deleting-web", i,".txt", collapse = ""))
    
    # SIMULATION
    
    deleted <- deleting2(in.dyn[[i]], in.n[[i]], n.regional, FR = funcres, frtune = x)
    del.props <- data.frame(deleted, locweb = i)
    
    print(del.props)
    
    sink()
    
    return(del.props)
  }
  
  stopCluster(cl)
  
  
  dels <- rbindlist(DELETION)
  
  return(dels)
}


