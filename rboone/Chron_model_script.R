# Script for running the SAM DEAD/LIVE chronology model

# starting with libraries
library(rjags)
library(coda)
library(mcmcplots)

load.module("dic")

# set working directory
# this is my personal computer directory:
# setwd("/Volumes/GoogleDrive/My Drive/SAM_DeadLive_trees")
setwd("G:/My Drive/SAM_DeadLive_trees")

# First we're going to run the model without inits and let JAGS generate them
# then we'll save inits from the posterior and run again
# this is my lazy way of making inints

# load Chronology data for trees that lived
# this dataset has values for deltaX to make the model run fast
load("./rboone/data/SAM_dat_living_chron_with_deltaX.R")

# model test run
alive1 <- jags.model("rboone/Chronology_species_model.R",
                     data = SAM_dat_living_chron_with_deltaX,
                     n.adapt = 1000,
                     n.chains = 3)
update(alive1, n.iter = 10000)
coda<-coda.samples(model = alive1,variable.names = c("a", "mu.a", "sig", "sig.a")
                   ,n.iter = 5000)

for (i in 1:3){
  alive1 <- jags.model("rboone/Chronology_species_model.R",
                       data = SAM_dat_living_chron_with_deltaX,
                       n.adapt = 1000,
                       n.chains = 1)
  
  update(alive1, n.iter = 10000)
  chain <- i
  for(z in 1:9){
    coda<-coda.samples(model = alive1,variable.names = c("a", "mu.a", "sig", "sig.a")
                       ,n.iter = 5000)
    save(coda,file=paste0("rboone/Chronology Model Output/alive_out", z, "_c", chain,".R"))
    rm(coda)
    gc()
  }
}
library(coda)
res<-mcmc.list()
z=1
load(paste0("./rboone/Chronology Model Output/alive_out",z,"_c1.R")) 
res[[1]]<-as.mcmc(coda[[1]])
load(paste0("./rboone/Chronology Model Output/alive_out",z,"_c2.R"))
res[[2]]<-as.mcmc(coda[[1]])
load(paste0("./rboone/Chronology Model Output/alive_out",z,"_c3.R"))
res[[3]]<-as.mcmc(coda[[1]])
for(z in 2:9){
  load(paste0("./rboone/Chronology Model Output/alive_out",z,"_c1.R"))
  res[[1]]<-as.mcmc(rbind(res[[1]],coda[[1]][seq(1,5000,35),]))
  load(paste0("./rboone/Chronology Model Output/alive_out",z,"_c2.R"))
  res[[2]]<-as.mcmc(rbind(res[[2]],coda[[1]][seq(1,5000,35),]))
  load(paste0("./rboone/Chronology Model Output/alive_out",z,"_c3.R"))
  res[[3]]<-as.mcmc(rbind(res[[3]],coda[[1]][seq(1,5000,35),]))
  rm(coda)
  gc()
}

# Rohan Home PC
source("Z:/common/Methods/Modeling tools/R/Mike/Rjags/MCMCrestart_v2.r")
# source("/Volumes/lab-ogle/common/Methods/Modeling tools/R/Mike/Rjags/MCMCrestart_v2.r")
# source("/Rtools/MCMCrestart_v2.r")
newinits <-  initfind(res)
# source("/Rtools/MCMCrestart_v2.r")
newinits <-  initfind(res)
newinits$variables # list variables
# newinits <- removevars(initsin = newinits, variables=
# c(2,3,6,9,10,11)) # remove non-variable nodes #and dump mu_a/sig_a
newinits<-newinits[[2]]
inits<-newinits
save(inits,file="./rboone/Chronology Model Inits/alive_inits1.R")

load("./rboone/data/SAM_dat_living_chron.R")
load("./rboone/Chronology Model Inits/alive_inits1.R")

alive2 <- jags.model("rboone/Chronology_Model2.R",
                     data = SAM_dat_living_chron,
                     n.adapt = 1000,
                     inits = inits,
                     n.chains = 3)

update(alive2, n.iter = 5000)

coda2 <- coda.samples(alive2, variable.names = c("alpha", "weight", "weightX", "cum.weight", "yr.w",
                                                             "a", "mu.a", "sig", "sig.a", "deviance",
                                                             "weightOrdered"),
                      n.iter = 1000)
save(coda2, file = "./rboone/alive_out.R")

mcmcplot(coda2)
