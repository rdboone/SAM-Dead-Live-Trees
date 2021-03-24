# Tree-ring memory model
# data files for Pinus Edulis (PIED) growing near Montrose, CO:
# ring widths (mm): BUGS_data_widths_Montrose_PIED
# total monthly precipitation (mm): BUGS_data_PPT_montrose
# mean daily temperature (C): BUGS_data_Tave_montose
# sample sizes and time-scale info: BUGS_data_sample_sizes
# In the PPT and Tave data, the variables are tave[r,c], where r (row)
# represents a year (as indicated by Year.tave or Year.ppt), and c is the 
# month (site = 1 for montrose datasets).

# Model code. Note that this model set "antecedent importance" weights = 0 in year 1 for Oct-Dec (post-growing season), and uses varying time periods (for weigths) for each year such that yrs 1-2 = 12 blocks of 1 month, year 3 = 6 blocks of 2 months, yrs 4-5 = 4 blocks of 3 months. Total # weights = 2*12 + 1*6 + 2*4 = 38.
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
	for(c in 1:Nsites){
		for(y in 1:Nyears){
			# Likelihood for log ring-width data:			
			LogWidth[y,c] ~ dnorm(mu.LogWidth[y,c], tau)
			# Replicated data for evaluating model fit:
			LogWidth.rep[y,c] ~ dnorm(mu.LogWidth[y,c], tau)
			# Mean model that incorporates age effect, antecedent climate effects (antX), and
			# previous ring-width, which is centered around a value that is close to the 
			# observed mean ring width of 0.5:
			#covariate centering: mean age = 50
			mu.LogWidth[y,c] <- a[1,c] + a[2,c]*antX[y,1,c] + 
						a[3,c]*antX[y,2,c] + a[4,c]*antX[y,1,c]*antX[y,2,c] + 
						a[5,c]*(AR1[y,c] - log(1+1))
			
			#data likelihood
			AR1[y, c]~dnorm(0, tau) #check with somebody about this please
			}
	}
  

# Compute antecedent climate variables for climate variable v and time "block" into 
# the past t (t = 1 is the current year):
  #vary these by status as well
for(v in 1:Nv){
	for(j in 1:Nblocks){
		# Assign a dirichlet(1,1,...,1) prior to the importance weights using the "delta-trick"
		# (as per the relationship between the dirichlet and gamma distributions). The "weightX"
		# terms are intermediate quantities that we generally don't care to make inferences 
		# about.
		deltaX[j,v] ~ dgamma(1,1)
		weightX[j,v] <- deltaX[j,v]/sum(deltaX[,v])
		}

# Compute importance weights of interest:
  for(t in 1:Nlag){
	# Compute yearly weights (these describe the relative importance of climate
	# conditions (for each variable v) occuring at different years into the past 
	# (again, t = 1 is current year):
	yr.w[t,v] <- sum(weight[,t,v])
	
	# Define monthly importance weights for every month x year into past combo:
	for(m in 1:12){
		# Unnormalized monthly weights; set monthly weight = 0 if current year (t = 1)
		# and m > 9.5 (this months Oct, Nov, Dec); divide by BlockSize to account for
		# time-periods of different lengths:
		delta[m,t,v] <- (deltaX[block[t,m],v]/BlockSize[t])*(1-equals(t,1)*step(m-9.5))
		# Here, "weight" is the importance weight for each month m given year t, times 
		# the importance of past year t, for climate variable v:
		weight[m,t,v] <- delta[m,t,v]/sumD[v]
		# Reorder the weights such that the indexing for weightOrdered[j,v] corresponds to
		# j = 1, 2, 3, ..., 12*Nlag is the weight for the most recent month (j = 1; Dec of ring year),
		# previous month (j = 2; Nov of ring year), ...., to the last month (j = 12*Nlag; Jan of prior
		# year Nlag):
		weightOrdered[(t-1)*12 + (12-m+1),v] <- weight[m,t,v]
		# Loop through all years in the tree-ring datasets such that i = 1 corresponds
		# to the first year of data (in this example, 1910). In this exmaple, the climate variables 
		# (ppt and Tave) begin at year 1901 (9 years earlier than ring widths, thus, Nlag can't be
		# greater than 9).
		for(i in 1:Nyears){
		  for(s in 1:Nsites){
		    # Compute the weighted climate variable for past year t and month m; v = 1
		    # corresponds to precip, v = 2 corresponds to temperature. Set value for site
		    # such that grab correct column in climate data (in this example, must set site = 1):
		    antX1[i,m,t,v,s] <- weight[m,t,v]*(equals(v,1)*ppt[(i+4)-t+1, m, s] + 
		                                         equals(v,2)*tave[(i+4)-t+1, m, s])
		    }
			}	
		}
		# Compute sum of deltas needed for computing importance weights:
		sumD1[t,v] <- sum(delta[,t,v])
	} 
	# Compute sum of deltas needed for computing importance weights:	
	sumD[v] <- sum(sumD1[,v])	

  for(s in 1:Nsites){
    for(i in 1:Nyears){
      for(t in 1:Nlag){
        # Compute the antecedent variable by summing the weighted climate 
        # variables. First sum over past months:
        ant.sum1[i,t,v,s] <- sum(antX1[i,,t,v,s])
      }
      # Now sum over past years, and center the covariate about the empirical sample mean
      # (xmean, read-in with data):
      ant.sum2[i,v,s] <- sum(ant.sum1[i,,v,s])
      antX[i,v,s] <- ant.sum2[i,v,s] - xave[s,v]
    }
  }
	
} 

# Compute cumulative monthly weights:
for(v in 1:Nv){
	for(t in 1:(12*Nlag)){
		cum.weight[t,v] <- sum(weightOrdered[1:t,v])
		}
	for(m in 1:12){
		for(t in 1:Nlag){
			# Compute month GIVEN year weights (i.e., within each year, these
			# monthly weights (alpha) add to one):
			alpha[m,t,v] <- weight[m,t,v]/yr.w[t,v]
			} 
		} 
	}
	
# Assign hierarhical priors to the core-level effects in the mean model,
# and assign relatively non-informative priors to population-level
# parameters:	
	for(k in 1:5){
		for(c in 1:Nsites){
			# Core-level hierarchical prior:
			a[k,c] ~ dnorm(mu.a[k],tau.a[k])
		  # a[k,c] ~ dnorm(mu.a[k, status[c]],tau.a[k, status[c]])
		}
		tau.a[k] <- pow(sig.a[k],-2)
		sig.a[k] ~ dunif(0,100)
		# for(st in 1:2){}
		# tau.a[k] <- pow(sig.a[k],-2)
		# sig.a[k] ~ dunif(0,100)
		# 
		}
	for(k in 1:5){
	  # for(st in 1:2){}
		mu.a[k] ~ dnorm(0,0.0001)
	  # mu.a[k, st] ~ dnorm(0,0.0001)
		}
	sig ~ dunif(0,100)
	tau <- pow(sig,-2)	

}