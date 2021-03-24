
library(rjags)
load.module("dic")

args = commandArgs(trailingOnly = TRUE)
snum = args[1]

#living
###########################################################################################
# load data
load(paste0("data/site", snum, "_alive"))
load("rboone/odcm_alive_inits1.R")

# run without inits
# alive1 <- jags.model("Ring_width_model_FOR_DISTRIBUTION.R",
#                             data = SAM_DL_live,
#                             n.adapt = 1000,
#                             n.chains = 3)
# 
# update(alive1, n.iter = 10000)
# coda1 <- coda.samples(alive1, variable.names = c("alpha", "weight", "weightX", "cum.weight", "yr.w",
#                                                                "a", "mu.a", "sig", "sig.a", "deviance",
#                                                                "weightOrdered"),
#                              n.iter = 1000)


for (i in 1:3){
  alive1 <- jags.model("Ring_width_model_FOR_DISTRIBUTION.R",
                             data = SAM_DL_live,
                             # inits = inits[[i]],
                             n.adapt = 1000,
                             n.chains = 1)
  
  update(alive.odcm1, n.iter = 10000)
  chain <- i
  for(z in 1:9){
    coda<-coda.samples(model = alive1,variable.names = c("a", "mu.a", "sig", "sig.a")
                       ,n.iter = 5000)
    save(coda,file=paste0("codas/", snum, "alive_out", z, "_c", chain,".R"))
    rm(coda)
    gc()
  }
}
library(coda)
res<-mcmc.list()
z=1
load(paste0("codas/", snum, "alive_out",z,"_c1.R")) 
res[[1]]<-as.mcmc(coda[[1]])
load(paste0("codas/", snum, "alive_out",z,"_c2.R"))
res[[2]]<-as.mcmc(coda[[1]])
load(paste0("codas/", snum, "alive_out",z,"_c3.R"))
res[[3]]<-as.mcmc(coda[[1]])
for(z in 2:9){
  load(paste0("codas/", snum, "alive_out",z,"_c1.R"))
  res[[1]]<-as.mcmc(rbind(res[[1]],coda[[1]][seq(1,5000,35),]))
  load(paste0("codas/", snum, "alive_out",z,"_c2.R"))
  res[[2]]<-as.mcmc(rbind(res[[2]],coda[[1]][seq(1,5000,35),]))
  load(paste0("codas/", snum, "alive_out",z,"_c3.R"))
  res[[3]]<-as.mcmc(rbind(res[[3]],coda[[1]][seq(1,5000,35),]))
  rm(coda)
  gc()
}
# source("Z:/common/Methods/Modeling tools/R/Mike/Rjags/MCMCrestart_v2.r")
source("/Rtools/MCMCrestart_v2.r")
newinits <-  initfind(res)
newinits$variables # list variables
# newinits <- removevars(initsin = newinits, variables=
# c(2,3,6,9,10,11)) # remove non-variable nodes #and dump mu_a/sig_a
newinits<-newinits[[2]]
inits<-newinits
save(inits,file="rboone/odcm_alive_inits1.R")



#dead
###########################################################################################
load("rboone/dead")
load("rboone/dead_inits1.R")
dead.model <- jags.model("rboone/Ring_width_model_FOR_DISTRIBUTION.R",
                           data = SAM_DL_dead,
                           #inits = inits,
                           n.adapt = 1000,
                           n.chains = 3)

update(dead.model, n.iter = 10000)
dead.coda <- coda.samples(dead.model, variable.names = c("alpha", "weight", "weightX", "cum.weight", "yr.w",
                                                             "a", "mu.a", "sig", "sig.a", "deviance",
                                                             "weightOrdered"),
                            n.iter = 1000)
save(dead.coda, file="dead_coda.R")
load(dead.coda)
mcmcplots::mcmcplot(dead.coda)


for (i in 1:3){
  odcm.dead.model <- jags.model("rboone/Ring_width_model_FOR_DISTRIBUTION.R",
                             data = SAM_DL_dead,
                             #inits = inits[[i]],
                             n.adapt = 1000,
                             n.chains = 1)
  
  update(odcm.dead.model, n.iter = 10000)
  chain <- i
  for(z in 1:9){
    coda<-coda.samples(model = odcm.dead.model,variable.names = c("a", "mu.a", "sig", "sig.a")
                       ,n.iter = 5000)
    save(coda,file=paste0("rboone/odcm_dead_out", z, "_c", chain,".R"))
    rm(coda)
    gc()
  }
}
library(coda)
res<-mcmc.list()
z=1
load(paste0("rboone/dead_out",z,"_c1.R")) 
res[[1]]<-as.mcmc(coda[[1]])
load(paste0("rboone/dead_out",z,"_c2.R"))
res[[2]]<-as.mcmc(coda[[1]])
load(paste0("rboone/dead_out",z,"_c3.R"))
res[[3]]<-as.mcmc(coda[[1]])
for(z in 2:9){
  load(paste0("rboone/dead_out",z,"_c1.R"))
  res[[1]]<-as.mcmc(rbind(res[[1]],coda[[1]][seq(1,5000,35),]))
  load(paste0("rboone/dead_out",z,"_c2.R"))
  res[[2]]<-as.mcmc(rbind(res[[2]],coda[[1]][seq(1,5000,35),]))
  load(paste0("rboone/dead_out",z,"_c3.R"))
  res[[3]]<-as.mcmc(rbind(res[[3]],coda[[1]][seq(1,5000,35),]))
  rm(coda)
  gc()
}
source("Z:/common/Methods/Modeling tools/R/Mike/Rjags/MCMCrestart_v2.r")
# source("/Volumes/lab-ogle/common/Methods/Modeling tools/R/Mike/Rjags/MCMCrestart_v2.r")
newinits <-  initfind(res)
newinits$variables # list variables
# newinits <- removevars(initsin = newinits, variables=
# c(2,3,6,9,10,11)) # remove non-variable nodes #and dump mu_a/sig_a
newinits<-newinits[[2]]
inits<-newinits
save(inits,file="rboone/dead_inits1.R")

