# setwd("/Volumes/cirrus/scratch/rdb273/SAM_DeadLive_trees")
# install.packages("tidyverse", repos = "https://cran.r-project.org/web/packages/dplyr/index.html")
library(tidyverse)
library(dplR)
wd <- getwd()
# args = commandArgs(trailingOnly = TRUE)

#set the site number from command line
# snum <- args[1]
snum <- 1

# for(i in 1:7){ # there are 117 sites
# snum <- i

#read in ring widths for individuals from that site
rw <- read.csv("RingWidths.csv")%>%
  filter(SiteNum == snum) %>%
  select(Year, TreeNum, Width, YearNum, Age)
new_colnames = c("Year", "CoreID", "Width", "yrID", "Age")
colnames(rw)<- new_colnames

Year <- rw %>%
  select(Year, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width) %>%
  select(Year)

rwi <- rw %>%
  select(Year, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width) %>%
  select(-Year) %>%
  detrend(method = "ModNegExp") %>%
  cbind(Year) %>%
  pivot_longer(1:41, names_to = "CoreID", values_to = "RWI") %>%
  group_by(CoreID)

rwi$CoreID <- as.numeric(rwi$CoreID)

rw <- rw %>%
  left_join(rwi, by = c("Year", "CoreID"))
summary(rw)

# trees <- read.csv("RingWidths.csv")

#get the site name
mort_sites <- as.data.frame(read.csv("MortalitySites.csv"))
current_site <- mort_sites$site[mort_sites$SiteNum == snum]

# sites <- read.csv("TreeTable.csv")%>%
#   select(site) %>%
#   count(site)
# 
# sites

#get the status (dead or alive) for each tree
sta <- read.csv("TreeTable.csv") %>%
  filter(site == current_site) %>%
  select(TreeNum, status)
new_colnames <- c("CoreID", "status")
colnames(sta) <- new_colnames
# head(sta)

#join the status data with the ring widths
rw.sta = left_join(rw, sta, by = "CoreID")
# head(rw.sta)

#separate living and dead trees into their own data structures
current_living <- rw.sta %>%
  filter(status == "LIVING")
# head(current_living)

current_dead <- rw.sta %>%
  filter(status == "DEAD")
# head(current_dead)

#create a lagged ring RWI structure
current_living.lag <- current_living %>%
  transmute(lagRWI = lag(RWI))
current_dead.lag <- current_dead %>%
  transmute(lagRWI = lag(RWI))

#convert data frame to matrix for easier handling
current_living <- as.matrix(current_living)
current_dead <- as.matrix(current_dead)

#get first and last yearIDs
lmin = as.numeric(min(current_living[,4]))
lmax = as.numeric(max(current_living[,4]))

dmin = as.numeric(min(current_dead[,4]))
dmax = as.numeric(max(current_dead[,4]))

#set number of years and cores
l.Nyears = lmax - lmin + 1
l.Ncores = length(unique(current_living[,2])) 
l.l = l.Nyears*l.Ncores

d.Nyears = dmax - dmin + 1
d.Ncores = length(unique(current_dead[,2]))
d.l = d.Nyears*d.Ncores

#make matrcies for RWI, lagRWI, and age
l.RWI.matrix <- matrix(nrow = l.Nyears, ncol = l.Ncores)
l.lagRWI.matrix <- matrix(nrow = l.Nyears, ncol = l.Ncores)
l.age <- matrix(nrow = l.Nyears, ncol = l.Ncores)
l.yrID <- as.numeric(current_living[,4])
l.CoreID <- rep(1:l.Ncores, each = l.Nyears)
l.RWI <- as.numeric(current_living[,3])
l.lagRWI <- as.numeric(current_living.lag[,1])
l.Age <- as.numeric(current_living[,5])

d.RWI.matrix <- matrix(nrow = d.Nyears, ncol = d.Ncores)
d.lagRWI.matrix <- matrix(nrow = d.Nyears, ncol = d.Ncores)
d.age <- matrix(nrow = d.Nyears, ncol = d.Ncores)
d.yrID <- as.numeric(current_dead[,4])
d.CoreID <- rep(1:d.Ncores, each = d.Nyears)
d.RWI <- as.numeric(current_dead[,3])
d.lagRWI <- as.numeric(current_dead.lag[,1])
d.Age <- as.numeric(current_dead[,5])

#filling RWI, age, and lagRWI matricies
for(k in 1:l.l){
  l.RWI.matrix[l.yrID[k],l.CoreID[k]] <- l.RWI[k]
  l.age[l.yrID[k],l.CoreID[k]] <- l.Age[k]
}
for(k in 1:l.l){
  l.lagRWI.matrix[l.yrID[k],l.CoreID[k]] <- l.lagRWI[k]
}

for(k in 1:d.l){
  d.RWI.matrix[d.yrID[k],d.CoreID[k]] <- d.RWI[k]
  d.age[d.yrID[k],d.CoreID[k]] <- d.Age[k]
}
for(k in 1:d.l){
  d.lagRWI.matrix[d.yrID[k],d.CoreID[k]] <- d.lagRWI[k]
}

#making logRWI and AR1matrices
l.LogRWI <- matrix(nrow = l.Nyears, ncol = l.Ncores)
for(y in 1:l.Nyears){
  for(c in 1:l.Ncores){
    l.LogRWI[y,c] <- log(l.RWI.matrix[y,c]+1)
  }
}
l.AR1 <- matrix(nrow = l.Nyears, ncol = l.Ncores)
for(y in 1:l.Nyears){
  for(c in 1:l.Ncores){
    l.AR1[y,c] <- log(l.lagRWI.matrix[y,c]+1)
  }
}
d.LogRWI <- matrix(nrow = d.Nyears, ncol = d.Ncores)
for(y in 1:d.Nyears){
  for(c in 1:d.Ncores){
    d.LogRWI[y,c] <- log(d.RWI.matrix[y,c]+1)
  }
}
d.AR1 <- matrix(nrow = d.Nyears, ncol = d.Ncores)
for(y in 1:d.Nyears){
  for(c in 1:d.Ncores){
    d.AR1[y,c] <- log(d.lagRWI.matrix[y,c]+1)
  }
}

#getting climate data 
l.ppt <- read.csv("pre_2017_11_28.csv") %>%
  filter(site_num == snum, year < max(as.numeric(current_living[,1]) + 1)) %>% 
  spread(month, pre) %>%
  select(-year, -site_num, -lon, -lat)
l.ppt <- as.matrix(l.ppt)
# dim(l.ppt)

l.temp <- read.csv("tmp_2017_11_28.csv") %>%
  filter(site_num == 76, year < max(as.numeric(current_living[,1]) + 1)) %>% 
  spread(month, tmp) %>%
  select(-year, -site_num, -lon, -lat)
l.temp <- as.matrix(l.temp)
# dim(l.temp)

d.ppt <- read.csv("pre_2017_11_28.csv") %>%
  filter(site_num == snum, year < max(as.numeric(current_dead[,1]) + 1)) %>% 
  spread(month, pre) %>%
  select(-year, -site_num, -lon, -lat)
d.ppt <- as.matrix(d.ppt)
# dim(d.ppt)

d.temp <- read.csv("tmp_2017_11_28.csv") %>%
  filter(site_num == 76, year < max(as.numeric(current_dead[,1]) + 1)) %>% 
  spread(month, tmp) %>%
  select(-year, -site_num, -lon, -lat)
d.temp <- as.matrix(d.temp)
# dim(d.temp)

SAM_DL_live <- list(Nv=2, Nlag=5, Nyears=l.Nyears, Ncores=l.Ncores,
                    xave=c(as.integer(mean(l.ppt)),as.integer(mean(l.temp))), site = snum,
                    block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                              13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                              25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                              31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                              35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                      .Dim=c(5,12)),
                    # deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                    BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                    LogRWI = l.LogRWI,
                    AR1 = l.AR1,
                    # age = l.age,
                    ppt = l.ppt,
                    tave = l.temp)
save(SAM_DL_live, file = paste0(wd, "/rboone/data/site", snum, "_living_rwi"))

SAM_DL_dead <- list(Nv=2, Nlag=5, Nyears=d.Nyears, Ncores=d.Ncores,
                    xave=c(as.integer(mean(d.ppt)),as.integer(mean(d.temp))), site = snum,
                    block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                              13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                              25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                              31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                              35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                      .Dim=c(5,12)),
                    deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                    BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                    LogRWI = d.LogRWI,
                    AR1 = d.AR1,
                    # age = d.age,
                    ppt = d.ppt,
                    tave = d.temp)
save(SAM_DL_dead, file = paste0(wd,"/rboone/data/site", snum, "_dead_rwi"))

# }