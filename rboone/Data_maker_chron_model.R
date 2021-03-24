# Script for compiling raw ring width data into chronology-level list compatible
# with RJAGS
# Author: Rohan David Boone

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

# make an empty matrix for the living chrons, with site-chron numbers as 
# columns and years as rows
chronMat.living <- matrix(nrow = nyears, ncol = nsites)
row.names(chronMat.living) <- c(yearmin:yearmax)
colnames(chronMat.living) <- c(1:nsites)

# same for dead chrons
chronMat.dead <- matrix(nrow = nyears, ncol = nsites)
row.names(chronMat.dead) <- c(yearmin:yearmax)
colnames(chronMat.dead) <- c(1:nsites)

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
    cbind(Year) %>% # put year back on
    filter(samp.depth >= 5) %>% # filter for sample depth
    select(-samp.depth) # remove sample depth0
  
  colnames(site.crn) <- c("rwi", "Year") # change column names
  rownames(site.crn) <- site.crn$Year # make the years into rownames for indexing
  site.crn <- select(site.crn, -Year) # remove year column
  
  # loop through all of the years represented in this site-chron
  for(j in rownames(site.crn)){
    # fill in chron mat with site-level standardised ring width inidicies
    # calling specific rows by their names (the year)
    chronMat.living[j,i] <- site.crn[j,] 
  }
}

# now we want to remove sites where there are more than 60 NAs (less than 50
# years of data)

# make an empy vecto in which we will store the amount of NAs in each chron
missing.living <- c()

# loop through columns in our chronology matrix
for(i in 1:nsites){
  # append our vector with the total number of NAs in each column
  missing.living[i] <- sum(is.na(chronMat.living[,i]))
}

# Now we're going to do all of that again on the dead chronologies

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
  
  site.crn <- as.data.frame(site.crn)%>% 
    detrend(method = "ModNegExp") %>% # detrend the age related growth signal
    chron() %>% # turn ring widths into a chronolgy
    cbind(Year) %>% # put year back on
    filter(samp.depth >= 5) %>% # filter for sample depth
    select(-samp.depth) # remove sample depth0
  
  colnames(site.crn) <- c("rwi", "Year") # change column names
  rownames(site.crn) <- site.crn$Year # make the years into rownames for indexing
  site.crn <- select(site.crn, -Year) # remove year column
  
  # loop through all of the years represented in this site-chron
  for(j in rownames(site.crn)){
    # fill in chron mat with site-level standardised ring width inidicies
    # calling specific rows by their names (the year)
    chronMat.dead[j,i] <- site.crn[j,] 
  }
}

# now we want to remove sites where there are more than 60 NAs (less than 50
# years of data)

# make an empy vecto in which we will store the amount of NAs in each chron
missing.dead <- c()

# loop through columns in our chronology matrix
for(i in 1:nsites){
  # append our vector with the total number of NAs in each column
  missing.dead[i] <- sum(is.na(chronMat.dead[,i]))
}

# now we are going to compare the number of missing years between missing and 
# dead chronologies. We want to remove the same sites from both living and dead 
# matrices, so we'll look at the missing.living and missing.dead vectors and 
# take the larger number from each

# start by making a new vector that will be the 'master'
missing.years <- c()

# loop through all sites again
for(i in 1:nyears){
  # check to see which one is bigger, and append our 'master' with the larger #
  if(missing.living[i] >= missing.dead[i]){
    missing.years[i] <- missing.living[i]
  }else{
    missing.years[i] <- missing.dead[i]
  }
}

# pull out a the positions of values in our master whose values are greater than
# 60
bad.chrons <- which(missing.years>60)

# filtering out the columns in our chron mats that correspond to our 'bad.chrons'
chronMat.living <- chronMat.living[,-bad.chrons]
chronMat.dead <- chronMat.dead[,-bad.chrons]

dim(chronMat.living) # looks like we're down to 83 chronologies

# checking to make sure we have the same chrons in both matrices
colnames(chronMat.living) == colnames(chronMat.dead)

# let's store our remaining siteIDs in a vector
# we will probably need this to filter our climate data later
remaining.sites <- colnames(chronMat.living)

# now we're going to generate our lag matrices
# lets define a vector of NAs to insert at the beginning of each chron matrix
na.vec <- rep(NA, times = dim(chronMat.living)[2])

# insert it at the beginning of our living and dead chron matrix to push 
# everything down one and then removing the last row
chronMat.living.lag <- rbind(na.vec, chronMat.living)[-dim(chronMat.living)[1],]
chronMat.dead.lag <- rbind(na.vec, chronMat.dead)[-dim(chronMat.dead)[1],]

# replacing the rownames in our lag matrix with the rownames in our unlagged one
rownames(chronMat.living.lag) <- row.names(chronMat.living)
rownames(chronMat.dead.lag) <- row.names(chronMat.dead)

# for the model, we want log rwis, so let's take the log of all of our data in 
# the chron mats and lag chron mats

# start with empty matrices with appropriate dimentions
logRWI.living <- matrix(nrow = nyears, ncol = length(remaining.sites))
lagRWI.living <- matrix(nrow = nyears, ncol = length(remaining.sites))
logRWI.dead <- matrix(nrow = nyears, ncol = length(remaining.sites))
lagRWI.dead <- matrix(nrow = nyears, ncol = length(remaining.sites))

# now let's fill in those matrices
for(i in 1:length(remaining.sites)){ # loop through sites (cols)
  for(j in 1:nyears){ # loop through years (rows)
    logRWI.living[j,i] <- log(chronMat.living[j,i])
    lagRWI.living[j,i] <- log(chronMat.living.lag[j,i])
    logRWI.dead[j,i] <- log(chronMat.dead[j,i])
    lagRWI.dead[j,i] <- log(chronMat.dead.lag[j,i])
  }
}

# let's get the climate data

# precipitation
pre <- read.csv("pre_2017_11_28.csv")
# making an empty nyears*12*nsites array to store the precip data
# an array of matrices, with each matrix containing the precip data for a site
# matrix columns as months, rows as years
ppt.array <- array(dim = c(length(unique(pre$year)), 12, length(remaining.sites)))

# fill in the array
# loop through sites
for(i in 1:length(remaining.sites)){
  ppt <- pre %>% # grab precip data
    filter(site_num == remaining.sites[i]) %>% # filter by site
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
temp.array <- array(dim = c(length(unique(tmp$year)), 12, length(remaining.sites)))

for(i in 1:length(remaining.sites)){
  temp <- tmp %>%
    filter(site_num == remaining.sites[i]) %>% 
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
for(i in 1:length(remaining.sites)){
  ppt.ave[i] <- as.integer(mean(ppt.array[,,i], na.rm = T))
  temp.ave[i] <- as.integer(mean(temp.array[,,i], na.rm = T))
}

# lets compile our data lists
# living chron data
SAM_dat_living_chron <- list(Nv=2, Nlag=5, Nyears=nyears,
                             Nsites = length(remaining.sites),
                   xave=cbind(ppt.ave, temp.ave),
                   block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                             13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                             25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                             31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                             35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                     .Dim=c(12,5)),
                   BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                   LogWidth = logRWI.living,
                   AR1 = lagRWI.living,
                   ppt = ppt.array,
                   tave = temp.array)
save(SAM_dat_living_chron, file = paste0("./rboone/data/SAM_dat_living_chron.R"))

# dead chron data
SAM_dat_dead_chron <- list(Nv=2, Nlag=5, Nyears=nyears,
                             Nsites = length(remaining.sites),
                             xave=cbind(ppt.ave, temp.ave),
                             block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                       25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                       31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                       35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                               .Dim=c(12,5)),
                             BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                             LogWidth = logRWI.dead,
                             AR1 = lagRWI.dead,
                             ppt = ppt.array,
                             tave = temp.array)
save(SAM_dat_dead_chron, file = paste0("./rboone/data/SAM_dat_dead_chron.R"))

# I'm also going to make data lists with values saved for deltaX
# this is useful for generating inits if you are stupid and lazy like me

# living chron data with deltaX
SAM_dat_living_chron_with_deltaX <- list(Nv=2, Nlag=5, Nyears=nyears,
                             Nsites = length(remaining.sites),
                             xave=cbind(ppt.ave, temp.ave),
                             block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                       13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                       25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                       31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                       35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                               .Dim=c(12,5)),
                             deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                             BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                             LogWidth = logRWI.living,
                             AR1 = lagRWI.living,
                             ppt = ppt.array,
                             tave = temp.array)
save(SAM_dat_living_chron_with_deltaX, file = paste0("./rboone/data/SAM_dat_living_chron_with_deltaX.R"))

# dead chron data with deltax
SAM_dat_dead_chron_with_deltaX <- list(Nv=2, Nlag=5, Nyears=nyears,
                           Nsites = length(remaining.sites),
                           xave=cbind(ppt.ave, temp.ave),
                           block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                     13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                     25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                     31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                     35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                             .Dim=c(12,5)),
                           deltaX = structure(.Data = rep(1, times = 76), .Dim = c(38, 2)),
                           BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                           LogWidth = logRWI.dead,
                           AR1 = lagRWI.dead,
                           ppt = ppt.array,
                           tave = temp.array)
save(SAM_dat_dead_chron_with_deltaX, file = paste0("./rboone/data/SAM_dat_dead_chron_with_deltaX.R"))
