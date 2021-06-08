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
    mu.LogWidth[r] <- a[1,CoreID[r]] +
      a[2,CoreID[r]]*antX[Year[r],1,SiteID[r],Status[r]] + 
      a[3,CoreID[r]]*antX[Year[r],2,SiteID[r],Status[r]] + 
      a[4,CoreID[r]]*antX[Year[r],1,SiteID[r],Status[r]]*antX[Year[r],2,SiteID[r],Status[r]] + 
      a[5,CoreID[r]]*(AR1[r] - log(1+1))
    
    #data likelihood
    # if no missing values, this is not necessary
    AR1[r]~dnorm(0, tau) #check with somebody about this please
  }
  
  
  # Compute antecedent climate variables for climate variable v and time "block" into 
  # the past t (t = 1 is the current year):
  #vary these by status as well
  for(v in 1:Nv){
    for(s in 1:2){ # status
      for(j in 1:Nblocks){
        # Assign a dirichlet(1,1,...,1) prior to the importance weights using the "delta-trick"
        # (as per the relationship between the dirichlet and gamma distributions). The "weightX"
        # terms are intermediate quantities that we generally don't care to make inferences 
        # about.
        deltaX[j,v,s] ~ dgamma(1,1)
        weightX[j,v,s] <- deltaX[j,v,s]/sum(deltaX[,v,s])
      }
      
      # Compute importance weights of interest:
      for(t in 1:Nlag){
        # Compute yearly weights (these describe the relative importance of climate
        # conditions (for each variable v) occuring at different years into the past 
        # (again, t = 1 is current year):
        yr.w[t,v,s] <- sum(weight[,t,v,s])
        
        # Define monthly importance weights for every month x year into past combo:
        for(m in 1:12){
          # Unnormalized monthly weights; set monthly weight = 0 if current year (t = 1)
          # and m > 9.5 (this months Oct, Nov, Dec); divide by BlockSize to account for
          # time-periods of different lengths:
          # set weight equal to zero past the growing season
          delta[m,t,v,s] <- (deltaX[block[m,t],v,s]/BlockSize[t])*(1-equals(t,1)*step(m-9.5)) 
          # use equals and step in place of if statements in jags
          # Here, "weight" is the importance weight for each month m given year t, times 
          # the importance of past year t, for climate variable v:
          weight[m,t,v,s] <- delta[m,t,v,s]/sumD[v,s]
          # Reorder the weights such that the indexing for weightOrdered[j,v] corresponds to
          # j = 1, 2, 3, ..., 12*Nlag is the weight for the most recent month (j = 1; Dec of ring year),
          # previous month (j = 2; Nov of ring year), ...., to the last month (j = 12*Nlag; Jan of prior
          # year Nlag):
          weightOrdered[(t-1)*12 + (12-m+1),v,s] <- weight[m,t,v,s]
          # Loop through all years in the tree-ring datasets such that i = 1 corresponds
          # to the first year of data (in this example, 1910). In this exmaple, the climate variables 
          # (ppt and Tave) begin at year 1901 (9 years earlier than ring widths, thus, Nlag can't be
          # greater than 9).
          for(c in 1:Nsites){
            # for(i in 5:Nyears){ 
            for(i in year.start[c]:year.end[c]){
              # Compute the weighted climate variable for past year t and month m; v = 1
              # corresponds to precip, v = 2 corresponds to temperature. Set value for site
              # such that grab correct column in climate data (in this example, must set site = 1):
              antX1[i,m,t,v,c,s] <- weight[m,t,v,s]*(equals(v,1)*ppt[(i)-t+1, m, c] + 
                                                       equals(v,2)*tave[(i)-t+1, m, c])
            }
          }	
        }
        # Compute sum of deltas needed for computing importance weights:
        sumD1[t,v,s] <- sum(delta[,t,v,s])
      } 
      # Compute sum of deltas needed for computing importance weights:	
      sumD[v,s] <- sum(sumD1[,v,s])	
      
      for(c in 1:Nsites){
        # for(i in 5:Nyears){
        for(i in year.start[c]:year.end[c]){
          for(t in 1:Nlag){
            # Compute the antecedent variable by summing the weighted climate 
            # variables. First sum over past months:
            ant.sum1[i,t,v,c,s] <- sum(antX1[i,,t,v,c,s])
          }
          # Now sum over past years, and center the covariate about the empirical sample mean
          # (xmean, read-in with data):
          ant.sum2[i,v,c,s] <- sum(ant.sum1[i,,v,c,s])
          antX[i,v,c,s] <- ant.sum2[i,v,c,s] - xave[c,v]
        }
      }
    }
  } 
  
  # Compute cumulative monthly weights:
  for(v in 1:Nv){
    for(s in 1:2){
      for(t in 1:(12*Nlag)){ 
        cum.weight[t,v,s] <- sum(weightOrdered[1:t,v,s])
      }
      for(m in 1:12){
        for(t in 1:Nlag){
          # Compute month GIVEN year weights (i.e., within each year, these
          # monthly weights (alpha) add to one):
          alpha[m,t,v,s] <- weight[m,t,v,s]/yr.w[t,v,s]
        } 
      } 
    }
  }
  
  # Assign hierarhical priors to the core-level effects in the mean model,
  # and assign relatively non-informative priors to population-level
  # parameters:	
  for(k in 1:5){
    for(r in 1:Ncores){
      # core-level hierarchical prior:
      a[k,r] ~ dnorm(mu.a[k,Status[r],SiteID[r]],tau.a[k,Status[r]])
    }
    # for (c in 1:Ncores) PUT IN CORE LOOOP
    for(s in 1:2){
      #site loop
      for(c in 1:Nsites){
        mu.a[k,s,c] ~ dnorm(mu.mu.a[k,s],tau.mu.a[k,s])
      }
      #  **use folded sig/folded cauchy for hierarchical variance priors**
      mu.mu.a[k,s] ~ dnorm(0,0.0001)
      tau.mu.a[k,s] <- pow(sig.mu.a[k,s],-2)
      sig.mu.a[k,s] ~ dunif(0,100)
      tau.a[k,s] <- pow(sig.a[k,s],-2)
      sig.a[k,s] ~ dunif(0,100)
      # mu.a[k, st] ~ dnorm(0,0.0001)
    }
    for(c in 1:Nsites){
      mu.a.diff[k,c] <- mu.a[k,1,c] - mu.a[k,2,c]
    }
    mu.mu.a.diff[k] <- mu.mu.a[k,1] - mu.mu.a[k,2]
  }
  sig ~ dunif(0,100)
  tau <- pow(sig,-2)	
  
}