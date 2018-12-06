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
#' @param Fm     2 x 2 Matrix of adult male genotypic frequencies
#' @param W_y.g  Fitness expressions for male gametes (i.e., fitness effects associated with the wild-type and mutant Y chromosomes)
#' @export
maleGameteSel  <-  function(Fm, W_y.g) {
	(Fm * c(1, W_y.g, 1, W_y.g)) / sum((Fm * c(1, W_y.g, 1, W_y.g)))
}

offspringSel  <-  function(Fm.g, Ff, W_c, W_o) {
	cbind(rowSums(Fm.g),rowSums(Fm.g))*rbind(Ff,Ff)*rbind(c(W_c[1], W_c[2]), c(W_o[1], W_o[2]))
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
#' mutant Y chromosomes into a population initially fixed for the wild-type 
#' compensatory allele OR a mutant comensatory mutation into a population
#' initially fixed for the mutant y chromosome.
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
	if(!yInvade & !mInvade | yInvade & mInvade) {
		stop('Check the settings for yInvade and mInvade. If you are interested
		      in simulating invasion of both types of mutants, consider using 
		      the yMitoCoEvolCycle() function')
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
	# genotypes are a vector, turned into a matrix
	# males: Y/M, Y/m, 
	#		 y/M, y/m
	# females: M, m
		if(yInvade & !mInvade) {
			qy.init  <-  2/N
			Fm  <-  matrix(c((1 - qy.init), 0,
								  qy.init,  0), 
						   nrow=2, byrow=TRUE)
			Ff  <-  c(1, 0)
		}
		if(!yInvade & mInvade) {
			qm.init  <-  2/N
			Fm  <-  matrix(c(0,  0,
							 1,  0), 
						   nrow=2, byrow=TRUE)
			Ff  <-  c((1 - qm.init), qm.init)
		}
		if(yInvade & mInvade) {
			Fm  <-  matrix(c((1 - qy.init), 0,
								  qy.init,  0), 
						   nrow=2, byrow=TRUE)
			Ff  <-  c((1 - qm.init), qm.init)
		}

	# Run simulation, saving frequencies for each generation
	if(saveTrajectories) {

		# Storage structures for individual simulation data
		qy.t   <-  rep(0, times=(4*N+1))
		qm.t   <-  rep(0, times=(4*N+1))
		E.qy   <-  rep(0, times=(4*N+1))
		E.qm   <-  rep(0, times=(4*N+1))
	
		# Initial frequency of mutant Y chromosome and compensatory allele
		qy.t[1]  <-  rowSums(Fm)[2]
		E.qy[1]  <-  rowSums(Fm)[2]
		qm.t[1]  <-  Ff[2]
		E.qm[1]  <-  Ff[2]
		
		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t[gen] != 0 & qm.t[gen] < 1) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			Fm.g    <-  maleGameteSel(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
	# males: Y/M, Y/m, 
	#		 y/M, y/m
	# females: M, m
			O     <-  offspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			Wbar  <-  sum(O)
			# 4) Expected frequencies in offspring, after selection
			E.O     <-  O/Wbar
			E.qy[gen+1]     <-  sum(E.O[2,])
			E.qm[gen+1]     <-  sum(E.O[,2])
			# 5) Realized frequencies in adults
			Fm  <-  matrix(as.vector(rmultinom(n=1, size=N/2, prob=E.O))/(N/2), ncol=2)
			qy.t[gen+1]  <-  sum(Fm[2,])
			qm.t[gen+1]  <-  rbinom(n = 1, size = N/2, prob = E.qm[gen+1])/(N/2)
			Ff  <-  c((1 - qm.t[gen+1]), qm.t[gen+1])

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
		qy.t  <-  rowSums(Fm)[2]
		E.qy  <-  rowSums(Fm)[2]
		qm.t  <-  Ff[2]
		E.qm  <-  Ff[2]

		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t != 0 & qy.t < qcrit.y & qm.t < qcrit.m & qm.t != qnull.m) {

#browser()
			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			Fm.g    <-  maleGameteSel(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
	# males: Y/M, Y/m, 
	#		 y/M, y/m
	# females: M, m
			O     <-  offspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			Wbar  <-  sum(O)
			# 4) Expected frequencies in offspring, after selection
			E.O   <-  O/Wbar
			print(E.)
			E.qy  <-  sum(E.O[2,])
			E.qm  <-  sum(E.O[,2])
			# 5) Realized frequencies in adults
			Fm    <-  matrix(as.vector(rmultinom(n=1, size=N/2, prob=E.O))/(N/2), ncol=2)
			qy.t  <-  sum(Fm[2,])
			qm.t  <-  rbinom(n = 1, size = N/2, prob = E.qm)/(N/2)
			Ff    <-  c((1 - qm.t), qm.t)

			gen  <-  gen+1

			# Has the inversion reached threshold frequency for establishment (pcrit)? 
			# When did it first reach pcrit?
			if(qy.t >= qcrit.y) {
				qyInv        <-  1
				qyTildeTime  <-  gen
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












#' Run a single stochastic simulation of a complete co-evolutionary cycle
#' between the mutant Y chromosomes and Mitochondrial compensatory locus
#' A cycle is initit
#' 
#' @title Find the deterministic equilibeium genotype frequencies prior to introducing the inversion
#' @param N         Population size
#' @param sm        Selection favouring mutant y male gametes
#' @param delta     Difference between sm and selection coefficient for selection against offspring 
#'                  resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                  (delta = sm - so)
#' @param sc        Cost of compensation selection coefficient against offspring resulting from matings
#'                  between wild-type fathers and mothers with mutant comensatory allele
#' @param uy        Per-generation mutation rate to male-beneficial antagonistic y chromosome (Y --> y) .
#' @param qy.init   Initial frequency of mutant y chromosome. Defaults to "singleCopy": q.init = 2/N
#' @param um        Per-generation mutation rate to the compensatory allele at the mitochondrial locus
#'                   (M --> m) .
#' @param qm.init   Initial frequency of mutant m allele at mitochondrial comensatory locus. 
#'                   Defaults to 0.
#' @export#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
yMitoCoEvolCycle  <-  function(N = N, sm = sm, delta = delta, sc = sc,
							   uy = 10e-5, qy.init = "singleCopy",
							   um = 10e-5, qm.init = 0, ...) {
	# Pre-emptive warnings
#	if(!yInvade & !mInvade | yInvade & mInvade) {
#		stop('Check the settings for yInvade and mInvade. If you are interested
#		      in simulating invasion of both types of mutants, consider using 
#		      the yMitoCoEvolCycle() function')
#	}

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
	if(qy.init == "singleCopy") {
		qy.init  <-  2/N
	}
	Fm  <-  matrix(c((1 - qy.init), 0,
						  qy.init,  0), 
				   nrow=2, byrow=TRUE)
	Ff  <-  c((1 - qm.init), qm.init)

	# Storage structures for individual simulation data
	qy.t   <-  c()
	qm.t   <-  c()
	E.qy   <-  c()
	E.qm   <-  c()

	# Initial frequency of mutant Y chromosome and compensatory allele
	qy.t[1]  <-  rowSums(Fm)[2]
	E.qy[1]  <-  rowSums(Fm)[2]
	qm.t[1]  <-  Ff[2]
	E.qm[1]  <-  Ff[2]
		
	## Start forward simulation with new mutant Y chromosome
	gen  <-  1
	while(qy.t[gen] < 1 | qm.t[gen] < 1) {

		## Step through recursions:
		# 1) Mutation
		mutateY  <-  runif(1) <= uy
		mutateM  <-  runif(1) <= um
		if(mutateY & qy.t[gen] < 1) {
			mutant  <-  c(1:4)[as.vector(rmultinom(n=1, size=1, prob=Fm)) == 1]
			if(any(mutant == c(1,3))) {
				Fm  <-  as.vector(Fm)
				Fm[mutant]  <-  Fm[mutant] - (2/N)
				Fm[mutant + 1]  <-  Fm[mutant + 1] + (2/N)
				Fm  <-  matrix(Fm, nrow=2)
			}
		}
		if(mutateM) {
			mutant  <-  rbinom(n=1, size=1, prob=Ff[2])
			if(mutant == 0) {
				Ff[1]  <-  Ff[1] - (2/N)
				Ff[2]  <-  Ff[2] + (2/N)
			}
		}
		# 2) Calculate frequency of mutant Y among male gametes
		Fm.g    <-  maleGameteSel(Fm=Fm, W_y.g=W_y.g)
		# 3) Mating and selection
		O     <-  offspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
		Wbar  <-  sum(O)
		# 4) Expected frequencies in offspring, after selection
		E.O     <-  O/Wbar
		if(any(E.O < 0)) {
			browser()
		}
		E.qy[gen+1]  <-  sum(E.O[2,])
		E.qm[gen+1]  <-  sum(E.O[,2])
		# 5) Realized frequencies in adults
		Fm           <-  matrix(as.vector(rmultinom(n=1, size=N/2, prob=E.O))/(N/2), ncol=2)
		qy.t[gen+1]  <-  sum(Fm[2,])
		qm.t[gen+1]  <-  rbinom(n = 1, size = N/2, prob = E.qm[gen+1])/(N/2)
		Ff           <-  c((1 - qm.t[gen+1]), qm.t[gen+1])

		gen  <-  gen+1
	}

	# When did mutant alleles establish and fix, and when was qTilde crossed?
	Time         <-  seq_along(qy.t)
	yInvTime     <-  Time[qy.t >= qcrit.y][1]
	yFixTime     <-  Time[qy.t == 1][1]
	qyTildeTime  <-  Time[qy.t >= qyTilde][1]
	qmAtqyTilde  <-  qm.t[qyTildeTime]
	mInvTime     <-  Time[qm.t >= qcrit.m][1]
	mFixTime     <-  gen
	
	# Save  simulation data
	res  <-  list(
				"sm"           =  sm,
				"delta"        =  delta,
				"sc"           =  sc,
				"qcrit.y"      =  qcrit.y,
				"qcrit.m"      =  qcrit.m,
				"qy.t"         =  qy.t[1:gen-1],
				"E.qy"         =  E.qy[1:gen-1],
				"qm.t"         =  qm.t[1:gen-1],
				"E.qm"         =  E.qm[1:gen-1],
				"yInvTime"     =  yInvTime,
				"yFixTime"     =  yFixTime,
				"qyTildeTime"  =  qyTildeTime,
				"qmAtqyTilde"  =  qmAtqyTilde,
				"mInvTime"     =  mInvTime,
				"mFixTime"     =  mFixTime
				)
return(res)
}




#' Run replicate stochastic simulation to generate a data set for plotting
#' the probability of establishment/fixation for invaison of mutant y and
#' mitochondrial compensatory mutations 
#' 
#' @title makeDataYMitoInv
#' @param N          Population size
#' @param sm         Selection favouring mutant y male gametes
#' @param delta.vals Difference between sm and selection coefficient for selection against offspring 
#'                   resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                   (delta = sm - so)
#' @param sc.vals    Cost of compensation selection coefficient against offspring resulting from matings
#'                   between wild-type fathers and mothers with mutant comensatory allele
#' @export#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
makeDataYMitoInv  <-  function(N = N, sm = sm, delta = delta.vals, sc = sc.vals) {

	for(i in 1:length(delta.vals)) {
		for(j in 1:length(sc.vals)) {
			
		}
	}
}