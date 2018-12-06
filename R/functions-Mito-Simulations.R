########################################################
#  Antagonistic coevolution between the sex-chromosomes
#
#  Necessary functions for stochastic simulation of 
#  invasion of Y-linked male-benefit alleles and  
#  MITOCHONDRIAL compensatory alleles. 
#
#  Author: Colin Olito
#
#  NOTES: 
#          


###############
# DEPENDENCIES
###############


######################
# Necessary functions  
######################


#' Frequency of mutant Y chromosome among male gametes
#'
#' @title Haplotype frequencies among gametes
#' @param qy     Adult frequency of mutant Y (of length = 25)
#' @param W_y.g  Fitness expressions for male gametes (i.e., fitness effects associated with the wild-type and mutant Y chromosomes)
#' @export
qy.g  <-  function(qy, W_y.g) {
	(qy * W_y.g) / (qy * W_y.g + (1 - qy))
}


#' @title Threshold y frequency for invasion of compensatory allele 
#' @param sm     Selection coefficient for mutant y chromosome
#' @param ho     Dominance coefficient for matings between wild-type 
#'               compensatory allele females and mutant males
#' @param so     Selection coefficient for matings between wild-type 
#'               compensatory allele females and mutant males
#' @param hc     Dominance coefficient for matings between mutant 
#'               compensatory allele females and mutant males
#' @param sc     Selection coefficient for matings between mutant 
#'               compensatory allele females and mutant males
#' @export
qyTilde.Mito  <-  function(sm, so, sc) {
	delta  <-  (sm - so)
	sc / (sc + (sm - delta)*(1 + sm))
}



#' Frequency of mutant Y chromosome among offspring
#'
#' @title Haplotype frequencies among gametes
#' @param qy     Adult frequency of mutant Y
#' @param W_y.o  Fitness expressions for male gametes (i.e., fitness effects associated with the wild-type and mutant Y chromosomes)
#' @param qm     Adult frequency of mutant m compensatory allele
#' @export
qy.o  <-  function(qy, W_y.o, qm) {
	(qy * ((1 - qm)*W_y.o[1] + qm*W_y.o[2])) / ((qy * ((1 - qm)*W_y.o[1] + qm*W_y.o[2])) + (1 - qy))
}


#' Run a single stochastic simulation for the invasion of a male-beneficial
#' mutant Y chromosome into a population initially fixed for the wild-type 
#' compensatory allele (i.e., the non-compensatory allele) 
#' 
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param N         Population size
#' @param sm        Selection favouring mutant y male gametes
#' @param delta     Difference between sm and selection coefficient for selection against offspring 
#'                  resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                  (delta = sm - so)
#' @param ho        dominance coefficient for so
#' @param sc        Cost of compensation selection coefficient against offspring resulting from matings
#'                  between wild-type fathers and mothers with mutant comensatory allele
#' @param hc        Dominance coefficient for sc
#' @param yInvade   Logical. If TRUE, simulates invasion from single copy of mutant y into population 
#'                  fixed for wild-type compensatory allele
#' @param m.init    Initial frequency of mutant m allele at mitochondrial comensatory locus. Defaults to 0.
#'                  be set to any vector of 3 frequencies that sum to 1
#' @param mInvade   Logical. If TRUE, simulates invasion from single copy of mutant mitochondrial compensatory 
#'                  allele into population initially fixed for the mutant y chromosome.
#' @param qy.init   Initial frequency of mutant y. Defaults to 1, but can be set to any frequency.
#' @param saveTrajectories  Save evolutionary trajectories of allele frequencies? Used for troubleshooting.
#' @export#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
yMitoInvadeFwdSim  <-  function(N = N, sm = sm, delta = delta, sc = sc,
								yInvade = TRUE, qm.init = 0,
								mInvade = FALSE, qy.init = 1, 
								saveTrajectories = FALSE, ...) {
	# Pre-emptive warnings
	if(!yInvade & !mInvade) {
		stop('Setting saveTrajectories = FALSE makes no sense if neither mutant allele is invading')
	}

	# Storage for result objects
	qyInv        <-  NA
	qyInvTime    <-  NA
	qyTildeTime  <-  NA
	qmInv        <-  NA

	# Define Fitness expressions
	so     <-  (sm - delta)
	W_y.g  <-  (1 + sm)
	W_o    <-  c((1 - so), 1)
	W_c    <-  c(1, (1 - sc))

	# Define threshold frequency for establishment of mutant Y
	qcrit.y  <-  8/(N*delta)
	qcrit.m  <-  8/(N*so)
	qyTilde  <-  qyTilde.Mito(sm = sm, so = so, sc = sc)

	# set initial genotype frequencies
		if(yInvade) {
			qy.init  <-  2/N
			qm       <-  qm.init
		}
		if(mInvade) {
			qy.init  <-  qy.init
			qm       <-  2/N
		}
		if(yInvade & mInvade) {
			qy.init  <-  qy.init
			Fii      <-  Fii.init
		}

	# Run simulation, saving frequencies for each generation
	if(saveTrajectories) {

		# Storage structures for individual simulation data
		qy.t   <-  rep(0, times=(4*N+1))
		qm.t   <-  rep(0, times=(4*N+1))
		E.qy   <-  rep(0, times=(4*N+1))
		E.qm   <-  rep(0, times=(4*N+1))
	
		# Initial frequency of mutant Y chromosome
		qy.t[1]  <-  qy.init
		E.qy[1]  <-  qy.t[1]
		E.qm[1]  <-  qm
		
		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t[gen] != 0 & qm.t[gen] < 1) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			E.qy_g  <-  qy.g(qy = qy.t[gen], W_y.g = W_y.g)
#			qy_g    <-  sum(rbinom(n = 1, size = N/2, prob = E.qy_g)) / (N/2)
			qy_g    <-  E.qy_g
			# 2) Mating and selection
			O  <-  matrix(c(((1 - qm)*W_c[1]), (qm*W_c[2]),
							((1 - qm)*W_o[1]), (qm*W_o[2])
							), nrow=2, byrow=TRUE)
			O     <- (O*matrix(c((1 - qy_g), (1 - qy_g), 
									  qy_g,       qy_g), nrow=2, byrow=TRUE))
			Wbar  <-  sum(O)
			# 4) Expected frequencies in offspring, after selection
			E.O     <-  O/Wbar
			E.qm[gen+1]     <-  colSums(E.O)[2]
			E.qy[gen+1]     <-  sum(E.O[2,])
			# 5) Realized frequencies in adults
#			if(mInvade) {
#				qy.t[gen+1]  <-  qy.init
#			} else{ 
				qy.t[gen+1]  <-  rbinom(n = 1, size = N/2, prob = E.qy[gen+1])/(N/2) 
#			}
			qm               <-  rbinom(n = 1, size = N/2, prob = E.qm[gen+1])/(N/2)
			qm.t[gen+1]      <-  qm

			gen  <-  gen+1
		}

		# Remove unnecessary 0's
		qy.t  <-  qy.t[qy.t != 0]
		E.qy  <-  E.qy[1:length(qy.t)]
		qm.t  <-  qm.t[1:length(qy.t)]
		E.qm  <-  E.qm[1:length(qy.t)]

		# Have the mutant alleles reached threshold frequencies for establishment?
		# When?

		if(any(qy.t >= qcrit.y)) {
			qyInv        <-  1
			qyTildeTime  <-  gen[qy.t >= qyTilde][1]
		}
		if(any(qm.t >= qcrit.m)) {
			qmInv      <-  1
		}
	
		# Save  simulation data
		res  <-  list(
					"yInvade"      =  yInvade,
					"mInvade"      =  mInvade,
					"delta"        =  delta,
					"qcrit.y"      =  qcrit.y,
					"qcrit.m"      =  qcrit.m,
					"qyFreq"       =  qy.t[1:gen-1],
					"E.qy"         =  E.qy[1:gen-1],
					"qmFreq"       =  qm.t[1:gen-1],
					"E.qm"         =  E.qm[1:gen-1],
					"Wbar"         =  Wbar[1:gen-1],
					"nGen"         =  gen,
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qmInv"        =  qmInv
 					)
	} 

	if(!saveTrajectories) {

		# hotfix to ensure while loop works when 
		# qy.init > qcrit.y
		if(mInvade) {
			qcrit.y  <-  1.1
			qnull.m  <-  0
		} else {qnull.m  <-  1}

		# Initial frequencies
		qy.t   <-  qy.init
		qm.t   <-  qm
		E.qy   <-  qy.t
		E.qm   <-  qm.t

		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t != 0 & qy.t < qcrit.y & qm.t < qcrit.m & qm.t != qnull.m) {

#browser()
			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			E.qy_g  <-  qy.g(qy = qy.t, W_y.g = W_y.g)
#			qy_g    <-  sum(rbinom(n = 1, size = N/2, prob = E.qy_g)) / (N/2)
			qy_g    <-  E.qy_g
			# 2) Mating and selection
			O  <-  matrix(c(((1 - qm)*W_c[1]), (qm*W_c[2]),
							((1 - qm)*W_o[1]), (qm*W_o[2])
							), nrow=2, byrow=TRUE)
			O     <- (O*matrix(c((1 - qy_g), (1 - qy_g), 
									  qy_g,       qy_g), nrow=2, byrow=TRUE))
			Wbar  <-  sum(O)
			# 4) Expected frequencies in offspring, after selection
			E.O     <-  O/Wbar
			E.qm    <-  colSums(E.O)[2]
			E.qy    <-  sum(E.O[2,])
			# 5) Realized frequencies in adults
#			if(mInvade) {
#				qy.t  <-  1
#			} else{
				qy.t  <-  rbinom(1, N/2, E.qy)/(N/2) 
#			}
			qm       <-  rbinom(1, N/2, E.qm)/(N/2)
			qm.t     <-  qm

			gen  <-  gen+1

			# Has the inversion reached threshold frequency for establishment (pcrit)? 
			# When did it first reach pcrit?
			if(qy.t >= qcrit.y) {
				qyInv        <-  1
				qyTildeTime  <-  gen[qy.t >= qyTilde][1]
			}
			if(qm.t >= qcrit.m) {
				qmInv      <-  1
			}
		}	

		# Save  simulation data
		res  <-  list(
					"yInvade"      =  yInvade,
					"mInvade"      =  mInvade,
					"delta"        =  delta,
					"qcrit.y"      =  qcrit.y,
					"qcrit.m"      =  qcrit.m,
					"qyFreq"       =  qy.t,
					"E.qy"         =  E.qy,
					"qmFreq"       =  qm.t,
					"E.qm"         =  E.qm,
					"Wbar"         =  Wbar,
					"nGen"         =  (gen-1),
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qmInv"        =  qmInv
					)
	} 
return(res)

}
