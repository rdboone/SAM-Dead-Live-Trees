###############################################################################
#
# A function to extract values from a given iteration for each chain in an MCMC
# object.  The main use of this fuction is to extract the final state of a 
# model for use in initializing the modle for further update.
#
# Written by Michael Fell on Monday, May 8, 2017
#
###############################################################################

intfind <- function(mcmcin, iteration=0){
	# If mcmc.list convert to mcmc
	if(is.mcmc.list(mcmcin)==TRUE){
		mcmcin <- mcmc(data=mcmcin, thin=1)
	}
	
	# Get the number of chains
	n.chains <- length(mcmcin)
	
	# get variable names from a list
	var.names <- colnames(mcmcin[[1]])
	var.dims <- dim(mcmcin[[1]])
	if(iteration==0){
		iteration <- var.dims[1]
	}
	
	if(sum(var.names=="deviance")>0){
		var.names <- var.names[-which(var.names=="deviance")]
	}
	
	# Get names and dimension of each variable since the output is a table
	var.names2 <- apply(X=as.matrix(var.names), MARGIN=c(1), FUN=strsplit, split="\\x5B", perl=TRUE)
	var.info <- matrix(nrow=length(var.names), ncol=5)
	var.info[,5] <- var.names

	# var names and dimentions in sep. columns
	for(i in 1:length(var.names2)){
		var.info[i,1:length(unlist(var.names2[[i]]))] <- unlist(var.names2[[i]])
	}
	# Determin if a variable is a scaler (0) or not
	var.info[which(is.na(var.info[,2])),4] <- 0
	var.info[which(!is.na(var.info[,2])),4] <- 1
	
	# seperate the dimensions and remove commas
	for(i in 1:dim(var.info)[1]){
		if(var.info[i,4]==1){
			var.val <- unlist(strsplit(x=var.info[i,2], split="\\x5D", perl=TRUE))
			var.val <- unlist(strsplit(x=unlist(var.val), split=",", perl=TRUE))
			
			if(length(var.val)==1){
				var.info[i,2] <- var.val
			}else if(length(var.val==2)){
				var.info[i,2:3] <- var.val
				var.info[i,4] <- 2
			}else{
				stop("a variable has more than 2 dimensions and this function thus can't work")
			}
		}
	}

	var.names <- unique(var.info[,1])
	
	# create a list to contain the new initials
	initsoutall <- list()

	# put in variables and their values
	for(i in 1:n.chains){
		initsout <- list()
		for(j in 1:length(var.names)){
			var.dim <- which(var.info[,1]==var.names[j])
			if(var.info[var.dim[1],4]==0){
				initsout[[j]] <- mcmcin[[i]][iteration, which(colnames(mcmcin[[i]])==var.info[var.dim,5])]
			}else if(var.info[var.dim[1],4]==1){
				# Create a vector of the correct length
				temp <- numeric(max(var.info[var.dim,2]))
				# Fill the vector
				for(k in 1:length(temp)){
					temp[k] <- mcmcin[[i]][iteration, which(colnames(mcmcin[[i]])==var.info[var.dim[k],5])]
				}
				initsout[[j]] <- assign(x=paste0(var.names[j]), value=temp)
			}else if(var.info[var.dim[1],4]==2){
				# create an appropriate two dimensional variable
				temp <- matrix(nrow=max(as.numeric(var.info[var.dim,2])), ncol=max(as.numeric(var.info[var.dim,3])))
				# Fill the matrix
				for(k in 1:length(var.dim)){
					temp[as.numeric(var.info[var.dim[k],2]), as.numeric(var.info[var.dim[k],3])] <- mcmcin[[i]][iteration, which(colnames(mcmcin[[i]])==var.info[var.dim[k],5])]
				}
				initsout[[j]] <- assign(x=paste0(var.names[j]), value=temp)
			}else{
				stop("your code likely has variables over 2 dimensions! I was not smart enough to solve this dimensionality problem in this version.")
			}
		}
		names(initsout) <- var.names
		initsoutall[[i]] <- initsout
	}
	
	return(initsoutall)
}