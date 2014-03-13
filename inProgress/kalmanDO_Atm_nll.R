
KFnllDO <- function(Params, DO, Aitch, PAR, Chla, Temp, Zmix, Press){
	
	#!Pseudocode #1: Initial guesses for B, C, and Q t
	Beta <- 1
	LightGPPCoef <- Params[1] 
	TempRCoef <- Params[2]
	
	Sea <- matrix(c(LightGPPCoef,TempRCoef,1),nrow=1,ncol=3,byrow=TRUE) #Collecting coefficients to be fit into a matrix.  

	Queue <- Params[3]#Variance of the process error
	
	# Kt <- KO2(Temp,Freq,Wind)
	#     RefT <- SatdConc(Temp,Press) #I removed multiplying by .942 because it seems that this has already been done in the sonde's calculations (when the sondes are calibrated in the bubble bath, they go to about 94.2%)
	
	Ewe <- matrix(nrow=length(DO), ncol=3)
	Ewe[,1] <- PAR
	Ewe[,2] <- Temp
	#Ewe[,3] <- (Kt * (RefT-DO))/Zmix #maybe not best to use DO here
	
		
	
	#!Pseudocode #2: Set first true value equal to first observation
	Alpha <- DO[1]#Let's give this model some starting values
	
	#!Pseudocode #3: Set process covariance, P, equal to Q
	Pea <- Queue #starting value
	
	NegLogLikes <- c(0,rep(NA,(length(DO)-1)))#Fill in all likelihood values with 0's
		
	#!Pseudocode #4: Starting with 2nd time step, build a Time Series of Alpha and P
	for(i in 2:length(DO)){
	Ewe[i,3] <- ((Kt[i-1] * (RefT[i-1]-Alpha))/Zmix[i-1])
		
		#!Pseudocode #4a: Predictions
		Alpha <- Beta*Alpha + Sea[1,1]*Ewe[i,1] + Sea[1,2]*log(Ewe[i,2]) + Sea[1,3]*Ewe[i,3] 
		Pea <- (Beta^2)*Pea*(Kt[i]^2) + Queue 
		
		#!Pseudocode #4b: Update Predictions
		Eta <- DO[i] - Alpha
		Eff <- Pea + Aitch

		
		Alpha <- Alpha + Pea/Eff*Eta
		Pea <- Pea - Pea*Pea/Eff 
		
		#!Pseudocode #5: Calculate the NLL
		NegLogLikes[i] <- .5*log(2*pi) + .5*log(Eff)+ .5*Eta*Eta/Eff
		}#End cycling through data// series
		
		
	NLL <- sum(NegLogLikes)#Sum up all likelihoods for the total for this model/ set of parameters
	return(NLL)#NLL is the item which nlm() will attemp to minimize, thus yielding the 'most likely' combinations of parameter values
	}#End function
	