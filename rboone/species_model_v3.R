# Tree-ring memory model
# data files for Pinus Edulis (PIED) growing near Montrose, CO:
# ring widths (mm): BUGS_data_widths_Montrose_PIED
# total monthly precipitation (mm): BUGS_data_PPT_montrose
# mean daily temperature (C): BUGS_data_Tave_montose
# sample sizes and time-scale info: BUGS_data_sample_sizes
# In the PPT and Tave data, the variables are tave[r,c], where r (row)
# represents a year (as indicated by Year.tave or Year.ppt), and c is the 
# month (site = 1 for montrose datasets).

# Model code. Note that this model set "antecedent importance" weights = 0 in year 1 for Oct-Dec (post-growing season),
# and uses varying time periods (for weigths) for each year such that yrs 1-2 = 12 blocks of 1 month, 
# year 3 = 6 blocks of 2 months, yrs 4-5 = 4 blocks of 3 months. Total # weights = 2*12 + 1*6 + 2*4 = 38.
model{
  #for(i in 1:Nobs){
  # Reorganize tree-ring data into a matrix (years x cores)
  # Note: this is not necessary, but it makes the code more efficient for this 
  # example application.
  #width.matrix[yrID[i],CoreID[i]] <- Width[i]
  #age[yrID[i],CoreID[i]] <- Age[i]
  #}
  # for(y in 1:Nyears){
  # 	for(c in 1:Ncores){
  # 		# log transform ring width (mm) data; because there are some missing rings
  # 		# (width = 0), add one to width value:
  # 		LogWidth[y,c] <- log(width.matrix[y,c]+1)
  # 		}
  # Data not used in ring width and climate data. If using newer versions of OpenBUGS,
  # can read-in unused data with-out error messages; if calling model from R (and running
  # via OpenBUGS or JAGS), can delete this line, and only supply data that are actually
  # used via R script.
  # nn[y] <- Year[y] + Year.tave[y] + Year.ppt[y] + xave[1]
  #	}
  
  # Likelihood and mean model (first year serves as the "initial condition" if
  # ring-widths are modeled as function of previous ring-width, so only
  # define likelihood for 2nd year of observation and beyond):
  for(r in 1:Nring){
    # Likelihood for log ring-width data:			
    LogWidth[r] ~ dnorm(mu.LogWidth[r], tau)
    # Replicated data for evaluating model fit:
    LogWidth.rep[r] ~ dnorm(mu.LogWidth[r], tau)
    # Mean model that incorporates age effect, antecedent climate effects (antX), and
    # previous ring-width, which is centered around a value that is close to the 
    # observed mean ring width of 0.5:
    #covariate centering: mean age = 50
    # KO: I've modified the indexing so that we "grab" Status and SiteID for the core that corresponds to ring r,
    # KO: Thus, need to provide Status and SiteID at the "level" of the individual cores (so the lenght of each should be 
    # the same as the total number of cores, not the total number of rings).
    mu.LogWidth[r] <- a[1,CoreID[r]] +
      a[2,CoreID[r]]*antX[Year[r],1,Status[CoreID[r]],SiteID[CoreID[r]]] + 
      a[3,CoreID[r]]*antX[Year[r],2,Status[CoreID[r]],SiteID[CoreID[r]]] + 
      a[4,CoreID[r]]*antX[Year[r],1,Status[CoreID[r]],SiteID[CoreID[r]]]*antX[Year[r],2,Status[CoreID[r]],SiteID[CoreID[r]]] + 
      a[5,CoreID[r]]*(AR1[r] - log(1+1))
    
    #data likelihood
    # if no missing values, this is not necessary
    # KO: Updating this since each core is missing at least one value, corresponding to the "ring" 
    # KO: before the first "observed" ring. This should probably also have it's own precision. 
    # KO: Not sure how much this will affect the model, but will be more realistic.
    AR1[r]~dnorm(mu.AR1[Year[r],SiteID[CoreID[r]]], tau.AR1)
  }
  
  # Compute antecedent climate variables for climate variable v and time "block" into 
  # the past t (t = 1 is the current year):
  #vary these by status as well
  for(v in 1:Nv){
    for(s in 1:2){ # status
      for(c in 1:Nsites){
        for(j in 1:Nblocks){
          # Assign a dirichlet(1,1,...,1) prior to the importance weights using the "delta-trick"
          # (as per the relationship between the dirichlet and gamma distributions). The "weightX"
          # terms are intermediate quantities that we generally don't care to make inferences 
          # about.
          deltaX[v,s,c,j] ~ dgamma(alpha[v,s,j],1)
          weightX[v,s,c,j] <- deltaX[v,s,c,j]/sum(deltaX[v,s,c,])
        }
        # Compute importance weights of interest:
        for(t in 1:Nlag){
          # Compute yearly weights (these describe the relative importance of climate
          # conditions (for each variable v) occuring at different years into the past 
          # (again, t = 1 is current year):
          yr.w[t,v,s,c] <- sum(weight[,t,v,s,c])
          
          # Define monthly importance weights for every month x year into past combo:
          for(m in 1:12){
            # Unnormalized monthly weights; set monthly weight = 0 if current year (t = 1)
            # and m > 9.5 (this months Oct, Nov, Dec); divide by BlockSize to account for
            # time-periods of different lengths:
            # set weight equal to zero past the growing season
            delta[m,t,v,s,c] <- (deltaX[v,s,c,block[m,t]]/BlockSize[t])*(1-equals(t,1)*step(m-9.5))
            # use equals and step in place of if statements in jags
            # Here, "weight" is the importance weight for each month m given year t, times 
            # the importance of past year t, for climate variable v:
            weight[m,t,v,s,c] <- delta[m,t,v,s,c]/sumD[v,s,c]
            # Reorder the weights such that the indexing for weightOrdered[j,v] corresponds to
            # j = 1, 2, 3, ..., 12*Nlag is the weight for the most recent month (j = 1; Dec of ring year),
            # previous month (j = 2; Nov of ring year), ...., to the last month (j = 12*Nlag; Jan of prior
            # year Nlag):
            weightOrdered[(t-1)*12 + (12-m+1),v,s,c] <- weight[m,t,v,s,c]
            # Loop through all years in the tree-ring datasets such that i = 1 corresponds
            # to the first year of data (in this example, 1910). In this exmaple, the climate variables 
            # (ppt and Tave) begin at year 1901 (9 years earlier than ring widths, thus, Nlag can't be
            # greater than 9).
            for(i in year.start[c]:year.end[c]){
            # for(i in 1:Nyears){
              # Compute the weighted climate variable for past year t and month m; v = 1
              # corresponds to precip, v = 2 corresponds to temperature. Set value for site
              # such that grab correct column in climate data (in this example, must set site = 1):
              antX1[i,m,t,v,s,c] <- weight[m,t,v,s,c]*(equals(v,1)*ppt[(i)-t+1, m, c] + 
                                                       equals(v,2)*tave[(i)-t+1, m, c])
            }
          }
          # Compute sum of deltas needed for computing importance weights:
          sumD1[t,v,s,c] <- sum(delta[,t,v,s,c])
        } 
        # Compute sum of deltas needed for computing importance weights:	
        sumD[v,s,c] <- sum(sumD1[,v,s,c])	
        for(i in year.start[c]:year.end[c]){
        # for(i in 1:Nyears){
          for(t in 1:Nlag){
            # Compute the antecedent variable by summing the weighted climate 
            # variables. First sum over past months:
            ant.sum1[i,t,v,s,c] <- sum(antX1[i,,t,v,s,c])
          }
          # Now sum over past years, and center the covariate about the empirical sample mean
          # (xmean, read-in with data):
          ant.sum2[i,v,s,c] <- sum(ant.sum1[i,,v,s,c])
          # Center antecedent climate variable around the overall, montly mean for each site-level climate variable:
          antX[i,v,s,c] <- ant.sum2[i,v,s,c] - xave[c,v]
        }
      }
    }
  } 
  
  # KO: Moved this plot of code for hierarchical weight priors up here, to following the
  # corresponding blocks of code
  # Priors for population (species) level importance weights.
    for(v in 1:Nv){ # precip and temperature
    for(s in 1:2){ # dead or alive
      # Since the first 3 "time blocks" are eventually assigned weight = 0 (since post-growing season),
      # just set their alpha's = 1, and don't use these alphas to compute the expected, pop-level weights.
      for(j in 1:3){
        alpha[v,s,j] <- 1
        }
      for(j in 4:Nblocks){ # blocks of time
        # Uniform prior for the Dirichlet parametes, alphas:
        alpha[v,s,j] ~ dunif(1,100)
        # Compute expected weight at population level; this is the "temporary" weight for
        # each block, and we need to "match-up" the blocks to specific year-lags and months.
        # Don't use the alpha values for the first 3 blocks as the correspond to post-growing season,
        # and their importance weights are automatially set to 0.
        EwX[v,s,j] <- alpha[v,s,j]/sum(alpha[v,s,4:Nblocks])
        for(m in 1:12){
          for(t in 1:Nlag){
            # Distributes the expected weight to associated month and lag-year, and sets
            # weights in Oct, Nov, Dec of current year (t = 1) to zero.
            Ew[m,t,v,s] <- (EwX[v,s,block[m,t]]/BlockSize[t])*(1-equals(t,1)*step(m-9.5))
            # Ordered pop-level weights:
            EwOrdered[(t-1)*12 + (12-m+1),v,s] <- Ew[m,t,v,s]
            }
          }
      }
    }
  }
  
  # Compute cumulative monthly weights:
  for(v in 1:Nv){
    for(s in 1:2){
      for(c in 1:Nsites){
        for(t in 1:(12*Nlag)){ 
          cum.weight[t,v,s,c] <- sum(weightOrdered[1:t,v,s,c])
        }
        for(m in 1:12){
          for(t in 1:Nlag){
            # Compute month GIVEN year weights (i.e., within each year, these
            # monthly weights (mgv.weights) add to one):
            mgv.weight[m,t,v,s,c] <- weight[m,t,v,s,c]/yr.w[t,v,s,c]
          } 
        }  
      }
    }
  }
  
    # Compute cumulative monthly weights (for pop-level):
  for(v in 1:Nv){
    for(s in 1:2){
        for(t in 1:(12*Nlag)){ 
          # Cumulative weights:
          cum.Ew[t,v,s] <- sum(EwOrdered[1:t,v,s])
          # Yearly weights:
          yr.Ew[t,v,s] <- sum(Ew[,t,v,s])
        }
        for(m in 1:12){
          for(t in 1:Nlag){
            # Compute month GIVEN year weights (i.e., within each year, these
            # monthly weights add to one):
            mgv.Ew[m,t,v,s] <- Ew[m,t,v,s]/yr.Ew[t,v,s]
          } 
        }  
      }
    }
  }
  
  # Assign hierarhical priors to the core-level effects in the mean model,
  # and assign relatively non-informative priors to population-level
  # parameters:	
  for(k in 1:5){ # loop through effects parameters.
    for(r in 1:Ncores){  # loop through cores.
      # core-level hierarchical prior:
      # KO: Status and SiteID should be vectors of length Ncores:
      a[k,r] ~ dnorm(mu.a[k,Status[r],SiteID[r]],tau.a[k,Status[r]])
    }
    # for (c in 1:Ncores) PUT IN CORE LOOOP
    for(s in 1:2){ # loop through status (dead/live)
      #site loop
      for(c in 1:Nsites){
        mu.a[k,s,c] ~ dnorm(mu.mu.a[k,s],tau.mu.a[k,s])
      }
      #  **use folded sig/folded cauchy for hierarchical variance priors**
      # KO: This is currently just using uniform priors for the sigs, which may be fine.
      mu.mu.a[k,s] ~ dnorm(0,0.0001)
      tau.mu.a[k,s] <- pow(sig.mu.a[k,s],-2)
      sig.mu.a[k,s] ~ dunif(0,100)
      tau.a[k,s] <- pow(sig.a[k,s],-2)
      sig.a[k,s] ~ dunif(0,100)
      # mu.a[k, st] ~ dnorm(0,0.0001)
    }
    # Compute pairwise differences to compare dead vs living parameters:
    # KO: Need to add a comment to define status (e.g., is status = 1 for dead or for living?).
    for(c in 1:Nsites){
      # Site level parameters:
      mu.a.diff[k,c] <- mu.a[k,1,c] - mu.a[k,2,c]
    }
    # Population or species-level:
    mu.mu.a.diff[k] <- mu.mu.a[k,1] - mu.mu.a[k,2]
  }
  
# Priors for the mean associated with the AR1 likelihood (for potential missing data):
for(c in 1:Nsites){
  for(y in 1:Nyears){
    mu.AR1[y,c] ~ dnorm(mu.mu.AR1, tau.muAR1)
  }
  }
  # KO: Prior for overall mean; may want to change the lower/upper limits on the uniform
# prior to better reflect the possible range of values for the AR1 variable.
  mu.mu.AR1[c] ~ dunif(-2,2)

  # Priors for standard deviations / precisions associated with the likelihoods:
  # KO: Will need to provide inits for sig.AR1 (use similar to what is use for sig)
  sig ~ dunif(0,100)
  tau <- pow(sig,-2)	
  sig.AR1 ~ dunif(0,100)
  tau.AR1 <- pow(sig.AR1,-2)	
  sig.muAR1 ~ dunif(0,100)
  tau.muAR1 <- pow(sig.muAR1,-2)
  
}

