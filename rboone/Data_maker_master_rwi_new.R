# setwd("/Volumes/cirrus/scratch/rdb273/SAM_DeadLive_trees")
# install.packages("tidyverse", repos = "https://cran.r-project.org/web/packages/dplyr/index.html")
library(tidyverse)
library(dplR)
wd <- getwd()
# args = commandArgs(trailingOnly = TRUE)

#set the site number from command line
# snum <- args[1]
# snum <- 8

for(i in 1:116){ # there are 117 sites
snum <- i

#read in ring widths for individuals from that site
rw <- read.csv("RingWidths.csv")%>%
  filter(SiteNum == snum) %>%
  select(SiteNum, Year, TreeNum, Width, YearNum, Age)
new_colnames = c("SiteNum", "Year", "CoreID", "Width", "yrID", "Age")
colnames(rw)<- new_colnames

# rw$Width[rw$Width == 0] <- 0.01

yrID <- sort(unique(rw$yrID))

rwi <- rw %>%
  select(yrID, CoreID, Width) %>%
  pivot_wider(names_from = CoreID, values_from = Width)

rwi <-rwi[order(rwi$yrID),] %>%
  select(-yrID) %>%
  detrend(method = "ModNegExp") %>%
  cbind(yrID) %>%
  pivot_longer(1:length(unique(rw$CoreID)),
               names_to = "CoreID",
               values_to = "RWI") %>%
  group_by(CoreID)

rwi$CoreID <- as.numeric(rwi$CoreID)

rw <- rw %>%
  left_join(rwi, by = c("yrID", "CoreID"))
summary(rw)

# trees <- read.csv("RingWidths.csv")

# get the site name
mort_sites <- as.data.frame(read.csv("MortalitySites.csv"))
current_site <- mort_sites$site[mort_sites$SiteNum == snum]

# sites <- read.csv("TreeTable.csv")%>%
#   select(site) %>%
#   count(site)
#
# sites

#get the status (dead or alive) for each tree
sta <- read.csv("TreeTable2.csv") %>%
  filter(site == current_site) %>%
  select(TreeNum, status, statusCode)
new_colnames <- c("CoreID", "status", "statusCode")
colnames(sta) <- new_colnames
# head(sta)

statusCode <- sta$statusCode
# 
# #join the status data with the ring widths
# rw.sta = left_join(rw, sta, by = "CoreID")
# # head(rw.sta)
# 
rw.lag <- rw %>%
  transmute(lagRWI = lag(RWI))
# 
# 
# #convert data frame to matrix for easier handling
# rw.sta <- as.matrix(rw.sta)

#get first and last yearIDs
min = as.numeric(min(rw[,5]))
max = as.numeric(max(rw[,5]))

#set number of years and cores
Nyears = length(unique(rw$yrID))
Ncores = length(unique(rw$CoreID)) 
l = Nyears*Ncores

#make matrcies for RWI, lagRWI, and age
RWI.matrix <- matrix(nrow = Nyears, ncol = Ncores)
lagRWI.matrix <- matrix(nrow = Nyears, ncol = Ncores)
yrID <- as.numeric(rw[,5])
# CoreID <- rep(1:Ncores, each = table(rw$CoreID))
CoreID <- as.vector(rw$CoreID - rw$CoreID[1] +1)
CoreID <- 1:length(Ncores)
for(i in 1:Ncores){
  
}
RWI <- as.numeric(rw[,7])
lagRWI <- as.numeric(rw.lag[,1])

#filling RWI, age, and lagRWI matricies
n <- 1
for(c in 1:Ncores){
  for(y in 1:Nyears){
    RWI.matrix[y,c] <- RWI[n]
    lagRWI.matrix[y,c] <- lagRWI[n] 
    n <- n+1
  }
}

#making logRWI and AR1matrices
LogRWI <- matrix(nrow = Nyears, ncol = Ncores)
AR1 <- matrix(nrow = Nyears, ncol = Ncores)
for(c in 1:Ncores){
  for(y in 1:Nyears){
    LogRWI[y,c] <- log(RWI.matrix[y,c]+1)
    AR1[y,c] <- log(lagRWI.matrix[y,c]+1)
  }
}


#getting climate data 
ppt <- read.csv("pre_2017_11_28.csv") %>%
  filter(site_num == snum, year < max(as.numeric(rw[,2]) + 1)) %>% 
  spread(month, pre) %>%
  select(-year, -site_num, -lon, -lat)
ppt <- as.matrix(ppt)
# dim(ppt)

temp <- read.csv("tmp_2017_11_28.csv") %>%
  filter(site_num == snum, year < max(as.numeric(rw[,2]) + 1)) %>% 
  spread(month, tmp) %>%
  select(-year, -site_num, -lon, -lat)
temp <- as.matrix(temp)
# dim(temp)


SAM_DL_dat <- list(Nv=2, Nlag=5, Nyears=Nyears, Ncores=Ncores,
                    xave=c(as.integer(mean(ppt)),as.integer(mean(temp))), site = 1,
                    block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                              13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                              25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                              31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                              35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                      .Dim=c(5,12)),
                    # deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                    BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                    LogRWI = LogRWI,
                    status = statusCode,
                    AR1 = AR1,
                    ppt = ppt,
                    tave = temp)
save(SAM_DL_dat, file = paste0(wd, "/rboone/data/site", snum, "dat"))
}
