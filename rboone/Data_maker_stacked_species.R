# Script for compiling raw ring width data into chronology-level list compatible
# with RJAGS
# Author: Rohan David Boonce

# Rohan home PC working directory
setwd("G:/My Drive/SAM_DeadLive_trees")

# Load libraries
library(dplR)
library(tidyverse)
library(R.utils)

SiteSpecies <- read.csv("SiteSpeciesTable.csv")

site.tbl <- as.data.frame(cbind(site = as.numeric(names(table(SiteSpecies$Site))),
                                num.sp = as.vector(table(SiteSpecies$Site)))) %>%
  filter(num.sp == 1)

# load data file with status (living/dead) data for each tree
status <- read.csv("TreeTable2.csv")%>%
  select(TreeNum, statusCode, latitude)

# load data file with raw ring width data, joining status data by coreID
raw <- read.csv("RingWidths.csv")%>%
  select(SiteNum, Year, TreeNum, Width, YearNum, Age, Species) %>%
  left_join(status, by = "TreeNum")

# changing colnames
new_colnames = c("SiteNum", "Year", "CoreID", "Width",
                 "yrID", "Age", "Species", "statusCode", "lat")
colnames(raw)<- new_colnames

# let's get the climate data

# precipitation
pre <- read.csv("pre_2017_11_28.csv")
tmp <- read.csv("tmp_2017_11_28.csv")

# get some values
(nsites <- length(unique(raw$SiteNum)))
(nyears <- length(unique(pre$year)))
(yearmin <- min(pre$year))
(yearmax <- max(pre$year))
(nspecies <- length(unique(raw$Species)))

# filter out only living trees
raw.living <- filter(raw, statusCode == 1)
# filter out only dead trees
raw.dead <- filter(raw, statusCode == 2)

# looping through all sites
for(i in unique(raw.living$SiteNum)){
  # first I'm going to pull out the species and latitude at site i
  site.sp <- raw.living$Species[raw.living$SiteNum == i][1]
  site.lat <- raw.living$lat[raw.living$SiteNum == i][1]
  # now I'm going to get the data into a dplR-compatible format
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
    cbind(Year = Year - yearmin + 1) %>% # put year back on
    filter(samp.depth >= 5) %>% # filter for sample depth
    mutate(logRWI = log(xxxstd)) %>% # take the log of the standardized ring widths
    mutate(lagRWI = lag(logRWI)) %>% # make a lagged logRWI column for autoregressive term
    select(-samp.depth, -xxxstd) # remove sample depth and original rwi columns
  
  # add a column with site ID and column with species ID
  site.crn <- cbind(site.crn, SiteID = rep(i, times = length(site.crn$logRWI)),
                    SpeciesID = rep(site.sp, times = length(site.crn$logRWI)),
                    lat = rep(site.lat, times = length(site.crn$logRWI))) 
  
  # make an object to stack all of the chronologies
  if(i == 1){
    chronMat.living <- site.crn
  }else{
    chronMat.living <- rbind(chronMat.living, site.crn)
  }
}
chronMat.living <- chronMat.living %>%
  filter(lat > 0) %>%
  filter(SiteID %in% site.tbl$site == TRUE)

# generate year.start and year.end vectors to store the first and last year 
# represented in each chronology
year.start.living <- c()
year.end.living <- c()
for(i in 1:length(unique(chronMat.living$SiteID))){
  year.start.living[i] <- min(chronMat.living$Year[chronMat.living$SiteID == i])
  year.end.living[i] <- max(chronMat.living$Year[chronMat.living$SiteID == i])
}


# looping through all sites
for(i in unique(raw.dead$SiteNum)){
  # first I'm going to pull out the species and latitude at site i
  site.sp <- raw.dead$Species[raw.dead$SiteNum == i][1]
  site.lat <- raw.dead$lat[raw.dead$SiteNum == i][1]
  # now I'm going to get the data into a dplR-compatible format
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
    cbind(Year = Year - yearmin + 1) %>% # put year back on
    filter(samp.depth >= 5) %>% # filter for sample depth
    mutate(logRWI = log(xxxstd)) %>% # take the log of the standardized ring widths
    mutate(lagRWI = lag(logRWI)) %>% # make a lagged logRWI column for autoregressive term
    select(-samp.depth, -xxxstd) # remove sample depth and original rwi columns
  
  # add a column with site ID and column with species ID
  site.crn <- cbind(site.crn, SiteID = rep(i, times = length(site.crn$logRWI)),
                    SpeciesID = rep(site.sp, times = length(site.crn$logRWI)),
                    lat = rep(site.lat, times = length(site.crn$logRWI))) 
  
  # make an object to stack all of the chronologies
  if(i == 1){
    chronMat.dead <- site.crn
  }else{
    chronMat.dead <- rbind(chronMat.dead, site.crn)
  } 
}
chronMat.dead <- chronMat.dead %>%
  filter(lat > 0) %>%
  filter(SiteID %in% site.tbl$site == TRUE)

# Now let's remove sites from both dataset so we only have those present in both
# vector containing only the sites present in both the living and dead matrices
remaining.sites <- intersect(chronMat.living$SiteID, chronMat.dead$SiteID)
# filtering the matrices
chronMat.living <- chronMat.living[chronMat.living$SiteID %in% remaining.sites == TRUE,]
chronMat.dead <- chronMat.dead[chronMat.dead$SiteID %in% remaining.sites == TRUE,]

# Now we have to renumber the sites to go from 1 the number of remaining sites
# create a vector of all possible sites (1:116)
sites <- seq(1,nsites,1)
# start by isolating the removed sites
removed.sites <- setdiff(sites, remaining.sites)

# now we're going to loop through all the sites in sites and insert zeros where
# there is a missing site. This will move all the sites up one each time we 
# encouter a missing site. This will make sense later.
for(i in 1:nsites){
  if(i %in% removed.sites){
    sites <- insert(sites, i, 0)
  }
}

# now we replace the site numbers in out chronMats by referenceing our sites 
# vector at the position of the original site number... maybe it still doesn't
# make sense, but trust me it works
for(i in 1:length(chronMat.living$SiteID)){
  chronMat.living$SiteID[i] <- sites[chronMat.living$SiteID[i]]
}

for(i in 1:length(chronMat.dead$SiteID)){
  chronMat.dead$SiteID[i] <- sites[chronMat.dead$SiteID[i]]
}

# generate year.start and year.end vectors to store the first and last year 
# represented in each chronology
year.start.living <- c()
year.end.living <- c()
for(i in 1:length(unique(chronMat.living$SiteID))){
  year.start.living[i] <- min(chronMat.living$Year[chronMat.living$SiteID == i])
  year.end.living[i] <- max(chronMat.living$Year[chronMat.living$SiteID == i])
}
year.start.dead <- c()
year.end.dead <- c()
for(i in 1:length(unique(chronMat.dead$SiteID))){
  year.start.dead[i] <- min(chronMat.dead$Year[chronMat.dead$SiteID == i])
  year.end.dead[i] <- max(chronMat.dead$Year[chronMat.dead$SiteID == i])
}


# saving old and new site numbers to convert back later
site.conversion <- cbind("original" = remaining.sites, "new" = unique(chronMat.living$SiteID))
write.csv(site.conversion, "rboone/site_conversion_table.csv")

# now we do the same thing to renumber the species IDs
remaining.species <- intersect(chronMat.living$SpeciesID, chronMat.dead$SpeciesID)
species <- seq(1,nspecies,1)
removed.species <- setdiff(species, remaining.species)
n = 0
for(i in 1:length(species)){
  if(i %in% removed.species){
    species <- insert(species, i, 0)
  }
}

for(i in 1:length(chronMat.living$SpeciesID)){
  chronMat.living$SpeciesID[i] <- species[chronMat.living$SpeciesID[i]]
}

for(i in 1:length(chronMat.dead$SpeciesID)){
  chronMat.dead$SpeciesID[i] <- species[chronMat.dead$SpeciesID[i]]
}

species.conversion <- cbind("original" = remaining.species,
                            "new" = unique(chronMat.living$SpeciesID))
write.csv(species.conversion, "rboone/species_conversion_table.csv")

site.species <- unique(chronMat.living[,c('SiteID','SpeciesID')])[,'SpeciesID']

# define new nyears, nsites, and nspecies after filtering
(nsites <- length(remaining.sites))
(nspecies <- length(remaining.species))

# now we process the climate data

# filtering out the sites that were removed from the data set
pre <- pre[pre$site_num %in% remaining.sites,]

# reassigning site numbers using the same logic as before
for(i in 1:dim(pre)[1]){
  pre$site_num[i] <- sites[pre$site_num[i]]
}
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

# filtering out the sites that were removed from the data set
tmp <- tmp[tmp$site_num %in% remaining.sites,]
# reassigning site numbers using the same logic as before
for(i in 1:dim(tmp)[1]){
  tmp$site_num[i] <- sites[tmp$site_num[i]]
}
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

(nyears <- length(unique(union(pre$year, tmp$year))))

# lets compile our data lists
# living chron data
SAM_dat_living_chron <- list(Nv=2, Nlag=5,
                             Nyears= nyears,
                             Nring = length(chronMat.living$Year),
                             Nsites = nsites, Nspecies = nspecies,
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
                             SpeciesID = chronMat.living$SpeciesID,
                             species.site = site.species,
                             year.start = year.start.living,
                             year.end = year.end.living,
                             ppt = ppt.array,
                             tave = temp.array)
save(SAM_dat_living_chron, file = paste0("./rboone/data/SAM_dat_living_chron.R"))

SAM_dat_dead_chron <- list(Nv=2, Nlag=5,
                           Nyears = nyears,
                           Nring = length(chronMat.dead$Year),
                           Nsites = nsites,  Nspecies = nspecies,
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
                             SpeciesID = chronMat.dead$SpeciesID,
                             species.site = site.species,
                             year.start = year.start.dead,
                             year.end = year.end.dead,
                             ppt = ppt.array,
                             tave = temp.array)
save(SAM_dat_dead_chron, file = paste0("./rboone/data/SAM_dat_dead_chron.R"))

# I'm also going to make data lists with values saved for deltaX
# this is useful for generating inits if you are stupid and lazy like me

# living chron data with deltaX
SAM_dat_living_chron_with_deltaX <- list(Nv=2, Nlag=5,
                                         Nyears= nyears,
                                         Nring = length(chronMat.living$Year),
                                         Nsites = nsites,  Nspecies = nspecies,
                                         xave=cbind(ppt.ave, temp.ave),
                                         block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                                   13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                                   25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                                   31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                                   35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                                           .Dim=c(12,5)),
                                         deltaX = structure(.Data = rep(1, times = 1748), .Dim = c(38, 2, 23)),
                                         BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                                         LogWidth = chronMat.living$logRWI,
                                         AR1 = chronMat.living$lagRWI,
                                         Year = chronMat.living$Year,
                                         SiteID = chronMat.living$SiteID,
                                         SpeciesID = chronMat.living$SpeciesID,
                                         species.site = site.species,
                                         year.start = year.start.living,
                                         year.end = year.end.living,
                                         ppt = ppt.array,
                                         tave = temp.array)
save(SAM_dat_living_chron_with_deltaX, file = paste0("./rboone/data/SAM_dat_living_chron_with_deltaX.R"))

# dead chron data with deltax
SAM_dat_dead_chron_with_deltaX <- list(Nv=2, Nlag=5,
                                       Nyears= nyears,
                                       Nring = length(chronMat.dead$Year),
                                       Nsites = nsites,  Nspecies = nspecies,
                                       xave=cbind(ppt.ave, temp.ave),
                                       block = structure(.Data=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                                                 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
                                                                 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
                                                                 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34, 34,
                                                                 35, 35, 35, 36, 36, 36, 37, 37, 37, 38, 38, 38),
                                                         .Dim=c(12,5)),
                                       deltaX = structure(.Data = rep(1, times = 1748), .Dim = c(38, 2, 23)),
                                       BlockSize = c(1, 1, 2, 3, 3), Nblocks = 38,
                                       LogWidth = chronMat.dead$logRWI,
                                       AR1 = chronMat.dead$lagRWI,
                                       Year = chronMat.dead$Year,
                                       SiteID = chronMat.dead$SiteID,
                                       SpeciesID = chronMat.dead$SpeciesID,
                                       species.site = site.species,
                                       year.start = year.start.dead,
                                       year.end = year.end.dead,
                                       ppt = ppt.array,
                                       tave = temp.array)
save(SAM_dat_dead_chron_with_deltaX, file = paste0("./rboone/data/SAM_dat_dead_chron_with_deltaX.R"))

