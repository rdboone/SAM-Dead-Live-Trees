# Script for compiling raw ring width data into chronology-level list compatible
# with RJAGS
# Author: Rohan David Boonce

# Rohan home PC working directory
setwd("G:/My Drive/SAM_DeadLive_trees")

# Load libraries
library(dplR)
library(tidyverse)

# load data file with status (living/dead) data for each tree
status <- read.csv("TreeTable2.csv")%>%
  select(TreeNum, statusCode)

# load data file with raw ring width data, joining status data by coreID
raw <- read.csv("RingWidths.csv")%>%
  select(SiteNum, Year, TreeNum, Width, YearNum, Age) %>%
  left_join(status, by = "TreeNum")

# changing colnames
new_colnames = c("SiteNum", "Year", "CoreID", "Width",
                 "yrID", "Age", "statusCode")
colnames(raw)<- new_colnames

# get some values
(nsites <- length(unique(raw$SiteNum)))
(nyears <- length(unique(raw$Year)))
(yearmin <- min(raw$Year))
(yearmax <- max(raw$Year))

# filter out only living trees
raw.living <- filter(raw, statusCode == 1)
# filter out only dead trees
raw.dead <- filter(raw, statusCode == 2)

# looping through all sites
for(i in 1:nsites){
  site.crn <- raw.living %>%
    filter(SiteNum == i) %>%
    select(Year, CoreID, Width) %>%
    # pivot data into dplR-compatible format
    # CoreID as columns, years as rows
    pivot_wider(names_from = CoreID, values_from = Width)
  
  # we'll have to pull out the year column, but we want to keep track of the 
  # years, so we make vector of the years represented in this site
  Year <- sort(unique(site.crn$Year))
  
  # sometimes the years get funky, so we're going to sort them to be safe
  site.crn <-site.crn[order(site.crn$Year),] %>%
    select(-Year) # remove year
  
  site.crn = as.data.frame(site.crn)%>% # need to convert from tibble to data frame
    detrend(method = "ModNegExp")%>% # detrend the age related growth signal
    chron() %>% # turn ring widths into a chronolgy
    cbind(Year = Year - min(Year) + 1) %>% # put year back on
    filter(samp.depth >= 5) %>% # filter for sample depth
    mutate(logRWI = log(xxxstd)) %>% # take the log of the standardized ring widths
    mutate(lagRWI = lag(logRWI)) %>% # make a lagged logRWI column for autoregressive term
    select(-samp.depth, -xxxstd) # remove sample depth and original rwi columns
  
  # add a column with site ID
  site.crn <- cbind(site.crn, SiteID = rep(i, times = length(site.crn$logRWI))) 
  
  # make an object to stack all of the chronologies
  if(i == 1){
    chronMat.living <- site.crn
  }else{
    chronMat.living <- rbind(chronMat.living, site.crn)
  }
}

# looping through all sites
for(i in 1:nsites){
  site.crn <- raw.dead %>%
    filter(SiteNum == i) %>%
    select(Year, CoreID, Width) %>%
    # pivot data into dplR-compatible format
    # CoreID as columns, years as rows
    pivot_wider(names_from = CoreID, values_from = Width)
  
  # we'll have to pull out the year column, but we want to keep track of the 
  # years, so we make vector of the years represented in this site
  Year <- sort(unique(site.crn$Year))
  
  # sometimes the years get funky, so we're going to sort them to be safe
  site.crn <-site.crn[order(site.crn$Year),] %>%
    select(-Year) # remove year
  
  site.crn = as.data.frame(site.crn)%>% # need to convert from tibble to data frame
    detrend(method = "ModNegExp")%>% # detrend the age related growth signal
    chron() %>% # turn ring widths into a chronolgy
    cbind(Year = Year - min(Year) + 1) %>% # put year back on
    filter(samp.depth >= 5) %>% # filter for sample depth
    mutate(logRWI = log(xxxstd)) %>% # take the log of the standardized ring widths
    mutate(lagRWI = lag(logRWI)) %>% # make a lagged logRWI column for autoregressive term
    select(-samp.depth, -xxxstd) # remove sample depth and original rwi columns
  
  # add a column with site ID
  site.crn <- cbind(site.crn, SiteID = rep(i, times = length(site.crn$logRWI)))
  
  # make an object to stack all of the chronologies
  if(i == 1){
    chronMat.dead <- site.crn
  }else{
    chronMat.dead <- rbind(chronMat.dead, site.crn)
  } 
}

# let's get the climate data

# precipitation
pre <- read.csv("pre_2017_11_28.csv")
# making an empty nyears*12*nsites array to store the precip data
# an array of matrices, with each matrix containing the precip data for a site
# matrix columns as months, rows as years
ppt.array <- array(dim = c(length(unique(pre$year)), 12, nsites))

# fill in the array
# loop through sites
for(i in 1:nsites){
  ppt <- pre %>% # grab precip data
    filter(site_num == i) %>% # filter by site
    spread(month, pre) %>% # spread into correct format
    select(-year, -site_num, -lon, -lat) # remove unwanted columns
  ppt <- as.matrix(ppt) # convert to matrix
  for(j in 1:116){ # year (row) loop
    for(k in 1:12){ # month (column) loop
      ppt.array[j,k,i] <- ppt[j,k] # fill in the ppt array
    }
  }
}

# same procedure for temp
tmp <- read.csv("tmp_2017_11_28.csv")
temp.array <- array(dim = c(length(unique(tmp$year)), 12, nsites))

for(i in 1:nsites){
  temp <- tmp %>%
    filter(site_num == i) %>% 
    spread(month, tmp) %>%
    select(-year, -site_num, -lon, -lat)
  temp <- as.matrix(temp)
  for(j in 1:116){
    for(k in 1:12){
      temp.array[j,k,i] <- temp[j,k]
    }
  }
}

# generating vector of average climate by site for covariate centering
ppt.ave <- c()
temp.ave <- c()
for(i in 1:nsites){
  ppt.ave[i] <- as.integer(mean(ppt.array[,,i], na.rm = T))
  temp.ave[i] <- as.integer(mean(temp.array[,,i], na.rm = T))
}

# lets compile our data lists
# living chron data
SAM_dat_living_chron <- list(Nv=2, Nlag=5,
                             Nyears= nyears,
                             Nring = length(chronMat.living$Year),
                             Nsites = nsites,
                             xave=cbind(ppt.ave, temp.ave),
                             block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                       25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                       31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                       35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                               .Dim=c(12,5)),
                             BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                             LogWidth = chronMat.living$logRWI,
                             AR1 = chronMat.living$lagRWI,
                             Year = chronMat.living$Year,
                             SiteID = chronMat.living$SiteID,
                             ppt = ppt.array,
                             tave = temp.array)
save(SAM_dat_living_chron, file = paste0("./rboone/data/SAM_dat_living_chron.R"))

SAM_dat_dead_chron <- list(Nv=2, Nlag=5,
                           Nyears = nyears,
                           Nring = length(chronMat.dead$Year),
                           Nsites = nsites,
                           xave=cbind(ppt.ave, temp.ave),
                           block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                       25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                       31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                       35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                               .Dim=c(12,5)),
                             BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                             LogWidth = chronMat.dead$logRWI,
                             AR1 = chronMat.dead$lagRWI,
                             Year = chronMat.dead$Year,
                             SiteID = chronMat.dead$SiteID,
                             ppt = ppt.array,
                             tave = temp.array)
save(SAM_dat_dead_chron, file = paste0("./rboone/data/SAM_dat_dead_chron.R"))

# I'm also going to make data lists with values saved for deltaX
# this is useful for generating inits if you are stupid and lazy like me

# living chron data with deltaX
SAM_dat_living_chron_with_deltaX <- list(Nv=2, Nlag=5,
                                         Nyears= nyears,
                                         Nring = length(chronMat.living$Year),
                                         Nsites = nsites,
                                         xave=cbind(ppt.ave, temp.ave),
                                         block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                                   13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                                   25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                                   31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                                   35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                                           .Dim=c(12,5)),
                                         deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                                         BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                                         LogWidth = chronMat.living$logRWI,
                                         AR1 = chronMat.living$lagRWI,
                                         Year = chronMat.living$Year,
                                         SiteID = chronMat.living$SiteID,
                                         ppt = ppt.array,
                                         tave = temp.array)
save(SAM_dat_living_chron_with_deltaX, file = paste0("./rboone/data/SAM_dat_living_chron_with_deltaX.R"))

# dead chron data with deltax
SAM_dat_dead_chron_with_deltaX <- list(Nv=2, Nlag=5,
                                       Nyears= nyears,
                                       Nring = length(chronMat.dead$Year),
                                       Nsites = nsites,
                                       xave=cbind(ppt.ave, temp.ave),
                                       block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                                 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                                 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                                 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                                 35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                                         .Dim=c(12,5)),
                                       deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                                       BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                                       LogWidth = chronMat.dead$logRWI,
                                       AR1 = chronMat.dead$lagRWI,
                                       Year = chronMat.dead$Year,
                                       SiteID = chronMat.dead$SiteID,
                                       ppt = ppt.array,
                                       tave = temp.array)
save(SAM_dat_dead_chron_with_deltaX, file = paste0("./rboone/data/SAM_dat_dead_chron_with_deltaX.R"))

