########################################################
#  Antagonistic coevolution between the sex-chromosomes
#
#  Necessary functions for stochastic simulation of 
#  invasion of Y-linked male-benefit alleles and  
#  AUTOSOMAL compensatory alleles. 
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

#' Haplotype frequencies of compensatory alleles among gametes
#'
#' @title Haplotype frequencies among gametes
#' @param Fii Vector of adult genotypic frequencies (of length = 3)
#' @export
x.A  <-  function(Fii=Fii) {
	(2*Fii[1] + Fii[2])/(2*sum(Fii))
} 
x.a  <-  function(Fii=Fii) {
	(2*Fii[3] + Fii[2])/(2*sum(Fii))
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
qyTilde.Auto  <-  function(sm, ho, so, hc, sc) {
	delta  <-  (sm - so)
	(hc*sc) / (hc*sc + (sm - delta)*(1 - ho)*(1 + sm))
}

#' Offspring frequencies after random mating
#'
#' @title Offspring frequencies after random mating
#' @param xi Vector of haplotype frequencies among gametes (of length = 5)
#' @export
offFreqCo  <-  function(xi) {
	O  <-  c(x.A^2, (2*x.A*x.a), x.a^2)
	O
}


#' Frequency of mutant Y chromosome among offspring
#'
#' @title Haplotype frequencies among gametes
#' @param qy     Adult frequency of mutant Y (of length = 25)
#' @param W_y.g  Fitness expressions for male gametes (i.e., fitness effects associated with the wild-type and mutant Y chromosomes)
#' @param Fii    Adult genotype frequencies at the Autosomal compensatory locus (vector of length 3)
#' @export
qy.o  <-  function(qy, W_y.o, Fii) {
	(qy * (Fii[1]*W_y.o[1] + Fii[2]*W_y.o[2] + Fii[3]*W_y.o[3])) / ((qy * (Fii[1]*W_y.o[1] + Fii[2]*W_y.o[2] + Fii[3]*W_y.o[3])) + (1 - qy))
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
#' @param Fii.init  Initial genotype frequencies at comensatory locus. Defaults to c(1,0,0), but can 
#'                  be set to any vector of 3 frequencies that sum to 1
#' @param aInvade   Logical. If TRUE, simulates invasion from single copy of mutant a compensatory 
#'                  allele into population initially fixed for the mutant y chromosome.
#' @param qy.init   Initial frequency of mutant y. Defaults to 1, but can be set to any frequency.
#' @param saveTrajectories  Save evolutionary trajectories of allele frequencies? Used for troubleshooting.
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
yAutoInvadeFwdSim  <-  function(N = N, sm = sm, delta = delta, ho = ho, sc = sc, hc = hc, 
								yInvade = TRUE, Fii.init = Fii.init,
								aInvade = FALSE, qy.init = 1, 
								saveTrajectories = FALSE, ...) {
	# Pre-emptive warnings
	if(!yInvade & !aInvade) {
		stop('Setting saveTrajectories = FALSE makes no sense if neither mutant allele is invading')
	}

	# Storage for result objects
	qyInv        <-  NA
	qyInvTime    <-  NA
	qyTildeTime  <-  NA
	qaInv        <-  NA

	# Define Fitness expressions
	so     <-  (sm - delta)
	W_y.g  <-  (1 + sm)
	W_o    <-  c((1 - so), (1 - ho*so), 1)
	W_c    <-  c(1, (1 - hc*sc), (1 - sc))

	# Define threshold frequency for establishment of mutant Y
	qcrit.y  <-  8/(N*delta)
	qcrit.a  <-  2/(N*so)
	qyTilde  <-  qyTilde.Auto(sm = sm, ho = ho, so = so, hc = hc, sc = sc)

	# set initial genotype frequencies
		if(yInvade) {
			qy.init  <-  2/N
			Fii      <-  Fii.init
		}
		if(aInvade) {
			qy.init  <-  qy.init
			Fii      <-  c((1 - (1/(N))), (1/(N)), 0)
		}
		if(yInvade & aInvade) {
			qy.init  <-  qy.init
			Fii      <-  Fii.init
		}

	# Run simulation, saving frequencies for each generation
	if(saveTrajectories) {

		# Storage structures for individual simulation data
		qy.t   <-  rep(0, times=(4*N+1))
		qa.t   <-  rep(0, times=(4*N+1))
		E.qy   <-  rep(0, times=(4*N+1))
		E.qa   <-  rep(0, times=(4*N+1))
	
		# Initial frequency of mutant Y chromosome
		qy.t[1]  <-  qy.init
		E.qy[1]  <-  qy.t[1]
		E.qa[1]  <-  Fii[3] + Fii[2]/2
		
		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t[gen] != 0 & qa.t[gen] < 1) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			E.qy_g  <-  qy.g(qy = qy.t[gen], W_y.g = W_y.g)
#			qy_g    <-  sum(rbinom(n = 1, size = N/2, prob = E.qy_g)) / (N/2)
			qy_g    <-  E.qy_g
			# 2) Mating and selection
			O  <-  matrix(c((Fii[1]*W_c[1] + (Fii[2]/2)*W_c[2])*(Fii[1] + Fii[2]/2), 
							(Fii[1]*W_c[1] + (Fii[2]/2)*W_c[2])*(Fii[3] + Fii[2]/2) +
							 (Fii[3]*W_c[3] + (Fii[2]/2)*W_c[2])*(Fii[1] + Fii[2]/2), 
							(Fii[3]*W_c[3] + (Fii[2]/2)*W_c[2])*(Fii[3] + Fii[2]/2),
							(Fii[1]*W_o[1] + (Fii[2]/2)*W_o[2])*(Fii[1] + Fii[2]/2), 
							(Fii[1]*W_o[1] + (Fii[2]/2)*W_o[2])*(Fii[3] + Fii[2]/2) +
							 (Fii[3]*W_o[3] + (Fii[2]/2)*W_o[2])*(Fii[1] + Fii[2]/2), 
							(Fii[3]*W_o[3] + (Fii[2]/2)*W_o[2])*(Fii[3] + Fii[2]/2)
							), nrow=2, byrow=TRUE)
			O     <- (O*matrix(c((1 - qy_g), (1 - qy_g), (1 - qy_g), 
									  qy_g,       qy_g,       qy_g), nrow=2, byrow=TRUE))
			Wbar  <-  sum(O)
			# 4) Expected frequencies in offspring, after selection
			E.O     <-  O/Wbar
			E.Fii    <-  colSums(E.O)
			E.qy[gen+1]     <-  sum(E.O[2,])
			E.qa[gen+1]     <-  E.Fii[3] + E.Fii[2]/2
			# 5) Realized frequencies in adults
			if(aInvade) {
				qy.t[gen+1]  <-  qy.init
			} else{ qy.t[gen+1]  <-  sum(rbinom(n = 1, size = N/2, prob = E.qy[gen+1]))/(N/2) }
			Fii          <-  as.vector(rmultinom(1, N, E.Fii)/(N))
			qa.t[gen+1]	 <-  Fii[3] + Fii[2]/2

			gen  <-  gen+1
		}

		# Remove unnecessary 0's
		qy.t  <-  qy.t[qy.t != 0]
		E.qy  <-  E.qy[1:length(qy.t)]
		qa.t  <-  qa.t[1:length(qy.t)]
		E.qa  <-  E.qa[1:length(qy.t)]

		# Have the mutant alleles reached threshold frequencies for establishment?
		# When?

		if(any(qy.t >= qcrit.y)) {
			qyInv        <-  1
			qyTildeTime  <-  gen[qy.t >= qyTilde][1]
		}
		if(any(qa.t >= qcrit.a)) {
			qaInv      <-  1
		}
	
		# Save  simulation data
		res  <-  list(
					"yInvade"      =  yInvade,
					"aInvade"      =  aInvade,
					"delta"        =  delta,
					"qcrit.y"      =  qcrit.y,
					"qcrit.a"      =  qcrit.a,
					"qyFreq"       =  qy.t[1:gen-1],
					"E.qy"         =  E.qy[1:gen-1],
					"qaFreq"       =  qa.t[1:gen-1],
					"E.qa"         =  E.qa[1:gen-1],
					"Wbar"         =  Wbar[1:gen-1],
					"nGen"         =  gen,
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qaInv"        =  qaInv
 					)
	} 

	if(!saveTrajectories) {

		# hotfix to ensure while loop works when 
		# qy.init > qcrit.y
		if(aInvade) {
			qcrit.y  <-  1.1
			qnull.a  <-  0
		} else {qnull.a  <-  1}

		# Initial frequencies
		qy.t   <-  qy.init
		qa.t   <-  Fii[3] + Fii[2]/2
		E.qy   <-  qy.t
		E.qa   <-  qa.t

		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t != 0 & qy.t < qcrit.y & qa.t < qcrit.a & qa.t != qnull.a) {

#browser()
			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			E.qy_g  <-  qy.g(qy = qy.t, W_y.g = W_y.g)
#			qy_g    <-  sum(rbinom(n = 1, size = N/2, prob = E.qy_g)) / (N/2)
			qy_g    <-  E.qy_g
			# 2) Mating and selection
			O  <-  matrix(c((Fii[1]*W_c[1] + (Fii[2]/2)*W_c[2])*(Fii[1] + Fii[2]/2), 
							(Fii[1]*W_c[1] + (Fii[2]/2)*W_c[2])*(Fii[3] + Fii[2]/2) +
							 (Fii[3]*W_c[3] + (Fii[2]/2)*W_c[2])*(Fii[1] + Fii[2]/2), 
							(Fii[3]*W_c[3] + (Fii[2]/2)*W_c[2])*(Fii[3] + Fii[2]/2),
							(Fii[1]*W_o[1] + (Fii[2]/2)*W_o[2])*(Fii[1] + Fii[2]/2), 
							(Fii[1]*W_o[1] + (Fii[2]/2)*W_o[2])*(Fii[3] + Fii[2]/2) +
							 (Fii[3]*W_o[3] + (Fii[2]/2)*W_o[2])*(Fii[1] + Fii[2]/2), 
							(Fii[3]*W_o[3] + (Fii[2]/2)*W_o[2])*(Fii[3] + Fii[2]/2)
							), nrow=2, byrow=TRUE)
			O     <- (O*matrix(c((1 - qy_g), (1 - qy_g), (1 - qy_g), 
									  qy_g,       qy_g,       qy_g), nrow=2, byrow=TRUE))
			Wbar  <-  sum(O)
			# 4) Expected frequencies in offspring, after selection
			E.O     <-  O/Wbar
			E.Fii    <-  colSums(E.O)
			E.qy     <-  sum(E.O[2,])
			E.qa     <-  E.Fii[3] + E.Fii[2]/2
			# 5) Realized frequencies in adults
			if(aInvade) {
				qy.t  <-  1
			} else{
				qy.t  <-  sum(rbinom(1, N/2, E.qy))/(N/2) 
			}
			Fii       <-  as.vector(rmultinom(1, N, E.Fii)/(N))
			qa.t	  <-  Fii[3] + Fii[2]/2

			gen  <-  gen+1

			# Has the inversion reached threshold frequency for establishment (pcrit)? 
			# When did it first reach pcrit?
			if(qy.t >= qcrit.y) {
				qyInv        <-  1
				qyTildeTime  <-  gen[qy.t >= qyTilde][1]
			}
			if(qa.t >= qcrit.a) {
				qaInv      <-  1
			}
		}	

		# Save  simulation data
		res  <-  list(
					"yInvade"      =  yInvade,
					"aInvade"      =  aInvade,
					"delta"        =  delta,
					"qcrit.y"      =  qcrit.y,
					"qcrit.a"      =  qcrit.a,
					"qyFreq"       =  qy.t,
					"E.qy"         =  E.qy,
					"qaFreq"       =  qa.t,
					"E.qa"         =  E.qa,
					"Wbar"         =  Wbar,
					"nGen"         =  (gen-1),
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qaInv"        =  qaInv
					)
	} 
return(res)

}


#' Introduce new mutant inverion genotype
#'
#' @title Introduce new mutant inverion genotype
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#' 					 chosen randomly, given initial genotypic frequencies (Fii.init). 
#' 					 See params for runReplicateAutoInvSims().
#' @param Fii.init   Initial genotypic frequencies, from which to calculate probability of new 
#' 					 mutant inversion occurring
#' @param N			 Population size
introduceInversion  <-  function(newMutant, Fii.init, N) {

	# Toggle
	specifyNewMutant  <-  is.numeric(newMutant)

	# Choose mutant genotype randomly
	if(!specifyNewMutant) {

		probNewMutant     <-  c(Fii.init[c(4,9,14,16:18)], Fii.init[19]*2)/sum(Fii.init[c(4,9,14,16:18)], Fii.init[19]*2)
		newMut            <-  c(4,9,14,16:19)[as.vector(rmultinom(1,1,probNewMutant)) == 1]

		# Subtract new mutant individual from frequency of old genotype
		Fii.init[newMut]  <-  Fii.init[newMut] - 1/N
	}

	# Specify mutant genotype
	if(specifyNewMutant) {
		# Subtract new mutant individual from frequency of old genotype
		newMut  <-  newMutant
		Fii.init[newMutant]  <-  Fii.init[newMutant] - 1/N
	}

	# Add mutant individual to frequency of new inversion genotype
	if(newMut == 4 | newMut == 9 | newMut == 14)
		Fii.init[newMut + 1]  <-  1/N
	if(newMut == 16 | newMut == 17 | newMut == 18)
		Fii.init[newMut + 5]  <-  1/N

	# if inversion occurs on abab genotype, choose randomly whether it occurs on
	# the maternally or paternally inherited chromosome 
	if(newMut == 19) {
		if(runif(1) >= 1/2) {
			Fii.init[newMut + 1]  <-  1/N
		}
		else Fii.init[newMut + 5]  <-  1/N
	}

	Fii.init
}


#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population 
#'
#' @title Wright-Fisher forward simulation of genotypic frequencies (default parameter values in parentheses)
#' @param nReps  Numer of replicate simulations. With no deleterious mutations, and introducing 
#' 				 a single copy of the inversion, it takes 1,600,000 replicate simulations to 
#'				 get 10,000 where the inversion successfully establishes.
#' @param N      Effective population size
#' @param m      Migration rate for locally maladaptive alleles (m =  0.01)
#' @param s      Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02)
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param n      Number of loci at which deleterious mutations may occur.
#' @param u      Mutation rate (default value of u = 1e-6).
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Defalut value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 				 to the loci involved in local adaptation. 
#' @param newMutant  Switch to choose whether to specify new mutant genotypes, or if they are 
#' 					 chosen randomly, given initial genotypic frequencies (Fii.init). If inversion
#' 					 genotypes are to be chosen randomly, set newMutant = 'random'. If they are to 
#' 					 be specified, set newMutant to one of the following values: c(4,9,14,16:19), 
#' 					 but consider which genotype each of these correspond to.
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSim`
#' @export
#' @author Colin Olito.
runReplicateAutoInvSims  <-  function(nReps = 1000, N = 500, m = 0.01, s = 0.1, h = 1/2, r = 0.1, 
									  n = 100, u = 1e-5, h.del = 0, s.del = 1, noDel = FALSE,
									  newMutant = 'random', saveTrajectories = FALSE) {

	##  Preemptive Warnings
	if(any(c(N,m,s,h,r,n,u,h.del) < 0) | any(c(m,s,h,r,u,h.del) > 1) | r > 0.5)
		stop('The chosen parameter values fall outside of reasonable parameter space')

	specifyNewMutant  <-  is.numeric(newMutant)
	if(specifyNewMutant & all(newMutant != c(4,9,14,16:19)))
		stop('If specifying the genotype of new inversion mutants, newMutant must take 
			  one of the following values: 4,9,14,16:19')

	if(!specifyNewMutant & newMutant != 'random')
		stop('If the genotype of new inversion mutants is being chosen randomly, 
			  the parameter newMutant must equal random')

	try({
		 if(m >= s )
			stop('Warning: migration is stronger than than selection, 
				  adaptive alleles will be swamped by maladaptive migrants')
	}, silent=FALSE)

	try({
		 if(nReps > 1000 & saveTrajectories)
			stop('Warning: You have chosen to save evolutionary trajectories 
				  for a large number of replicate simulations. Thiss will be
				  memory intensive. Consider setting saveTrajectories = FALSE')
	}, silent=FALSE)

	##  Define Fitness Expressions for determining eq. frequencies in absence of inversion
	W.init  <-  c(1,          (1 + h*s),         (1 + h*s),         (1 + h*s)^2,        0,
			     (1 + h*s),   (1 + s),           (1 + h*s)^2,       (1 + h*s)*(1 + s),  0,
			     (1 + h*s),   (1 + h*s)^2,       (1 + s),           (1 + s)*(1 + h*s),  0,
			     (1 + h*s)^2, (1 + h*s)*(1 + s), (1 + s)*(1 + h*s), (1 + s)^2,          0,
			     0,           0,                 0,                 0,                  0)
  		
 	## Find deterministic equilibrium frequencies in absence of inversion  
	Fii.init  <-  findEqFreqs(W=W.init, m=m, r=r, threshold=1e-7)

	# Introduce rare mutant inversion
	Fii.init  <-  introduceInversion(newMutant=newMutant, Fii.init=Fii.init, N=N)

	# Storage structures for replicate simulation data
	finalInvFreq    <-  rep(0, times=nReps)
	finalE.InvFreq  <-  rep(0, times=nReps)
	finalW.mean     <-  rep(0, times=nReps)
	nGen            <-  rep(0, times=nReps)
	invEst          <-  rep(0, times=nReps)
	invEstTime      <-  rep(0, times=nReps)
	nDels           <-  rep(0, times=nReps)

	if(saveTrajectories) {
		replicateTraj  <-  c()
		InvFreqTraj    <-  c()
		E.InvFreqTraj  <-  c()
		W.meanTraj     <-  c()
	} 

	# Replicate simulation loop
#	print('Running Wright-Fisher Forward Simulations')
	pb   <-  txtProgressBar(min=0, max=nReps, style=3)
	setTxtProgressBar(pb, 0)
	for(i in 1:nReps) {

	## Sample stationary distribution of deleterious alleles
	delMutFreq  <-  rejectionSampler(n=n, Ne=N, u=u, s=s.del, h=h.del)
	n.del       <-  sum(delMutFreq > runif(n=n))

	# Define fitness expressions, including fitness effects of deleterious mutations
	if(noDel) {
		s.del  <-  0
	}
	W  <-  c(1,                                  (1 + h*s),                                 (1 + h*s),                                 (1 + h*s)^2,                       (1 + h*s)^2*(1 - h.del*s.del)^n.del,
			(1 + h*s),                           (1 + s),                                   (1 + h*s)^2,                               (1 + h*s)*(1 + s),                 (1 + h*s)*(1 + s)*(1 - h.del*s.del)^n.del,
			(1 + h*s),                           (1 + h*s)^2,                               (1 + s),                                   (1 + s)*(1 + h*s),                 (1 + s)*(1 + h*s)*(1 - h.del*s.del)^n.del,
			(1 + h*s)^2,                         (1 + h*s)*(1 + s),                         (1 + s)*(1 + h*s),                         (1 + s)^2,                         (1 + s)^2*(1 - h.del*s.del)^n.del,
			(1 + h*s)^2*(1 - h.del*s.del)^n.del, (1 + h*s)*(1 + s)*(1 - h.del*s.del)^n.del, (1 + s)*(1 + h*s)*(1 - h.del*s.del)^n.del, (1 + s)^2*(1 - h.del*s.del)^n.del, (1 + s)^2*(1 - s.del)^n.del)

		## RUN SIMULATION
		repRes  <-  autoInvFwdSim(Fii.init=Fii.init, N=N, W=W, m=m, r=r, saveTrajectories=saveTrajectories)

		# save results for each replicate
		finalInvFreq[i]    <-  repRes$InvFreq[length(repRes$InvFreq)]
		finalE.InvFreq[i]  <-  repRes$E.InvFreq[length(repRes$E.InvFreq)]
		finalW.mean[i]     <-  repRes$W.mean
		nGen[i]            <-  repRes$nGen
		invEst[i]          <-  repRes$invEst
		invEstTime[i]      <-  repRes$invEstTime
		nDels[i]           <-  n.del

		if(saveTrajectories) {
			replicateTraj  <-  c(replicateTraj, rep(i, times=length(repRes$InvFreq)))
			InvFreqTraj    <-  c(InvFreqTraj, repRes$InvFreq)
			E.InvFreqTraj  <-  c(E.InvFreqTraj, repRes$E.InvFreq)
			W.meanTraj     <-  c(W.meanTraj, repRes$W.mean)
			} 

	setTxtProgressBar(pb, i)
	}

	# Save results and return results as a list
	results.df  <-  data.frame(
							   "finalInvFreq"    =  finalInvFreq,
							   "finalE.InvFreq"  =  finalE.InvFreq,
							   "finalW.mean"     =  finalW.mean,
							   "nGen"            =  nGen,
							   "invEst"          =  invEst,
							   "invEstTime"      =  invEstTime,
							   "nDels"           =  nDels
							   )
	if(saveTrajectories) {
		traj.df  <-  data.frame(
								"replicateTraj"  =  replicateTraj,
								"InvFreqTraj"    =  InvFreqTraj,
								"E.InvFreqTraj"  =  E.InvFreqTraj,
								"W.meanTraj"     =  W.meanTraj
								)
	} else {
		traj.df  <-  NULL
	}
	res  <-  list(
				  "results.df"  =  results.df,
				  "traj.df"     =  traj.df
				  )
	return(res)
}


#' Wrapper function to run replicate forward simulations for invasion
#' of autosomal inversions in a Wright-Fisher population 
#' USING DESIGNATED PARAMETER VALUES
#'
#' @title Run replicate Wright-Fisher forward simulations for autosomal inversion under different parameter values 
#' @param N.vals Desired population sizes
#' @param m.vals desired migration rates for locally maladaptive alleles (m =  0.01)
#' @param s.del.vals  Desired selection coefficients for deleterious mutations (default value of s = 1).
#' @param s      Selective advantage of locally adaptive alleles over migrant alleles (s = 0.02)
#' @param h      Dominance coefficient for locally adaptive alleles relative to migrant alleles (h = 0.5)
#' @param r      Recombination rate among the two loci involved in local adaptation (r = 0.1).
#' @param n      Number of loci at which deleterious mutations may occur.
#' @param u      Mutation rate (default value of u = 1e-6).
#' @param h.del  Dominance of deleterious mutations (default value of h = 0).
#' @param noDel  Omit fitness effects of deleterious mutations? Defalut value of FALSE assumes that 
#' 				 selection against deleterious mutations is twice as strong as selection favouring
#' 				 the locally adaptive alleles. Setting noDel = TRUE runs the W-F simulations as if
#' 				 there were no delterious mutations segregating in the population that are linked
#' 				 to the loci involved in local adaptation. 
#' @param saveTrajectories  Save evolutionary trajectories of inversion frequencies? Setting this 
#' 							to TRUE can become extremely memory intensive if you are running many
#' 							replicate simulations (see warning).
#' 							otherwise
#' @seealso `offFreq`, `findEqFreqs`, `autoInvFwdSim`
#' @export
#' @author Colin Olito.
makeFastReplicateAutoSexSpecInvSimsData(nReps = 100, N = 20000, h = 1/2, 
										m.vals = c(0.0005, 0.001), m.deltas = NULL,
										s.vals = c(0.001, 0.05), s.deltas = NULL, 
										r.vals = seq(from = 0, to = 0.5, by = 0.025),
										n = 100, u = 1e-5, h.del = 0, noDel = FALSE, 
										fastSim = TRUE, newMutant=c("random","random"), saveTrajectories = FALSE)
	# Simulate deleterious mutations that are either 
	# 1) recessive lethals OR
	# 2) recessive experiencing purifying selection
	#    that is twice as strong as the selective 
	#    advantage of the locally adaptive allels  
	s.del.vals = c(0, 1, 2*s)


	# create empty data frame with same structure as we are going to need
	data  <-  data.frame(matrix(ncol=10, nrow=0))

	# Convenience variables to monitor progress
	prog  <-  0
	tot   <-  length(r.vals)*length(m.vals)*length(s.del.vals)

	# Loop over parameter values we want to explore 
	for(j in 1:length(r.vals)) {
		for(k in 1:length(m.vals)) {
			for(l in 1:length(s.del.vals)) {

				# Display progress in terminal
				prog  <-  prog + 1
				cat("\n",paste('Running simulations for parameter set ', prog, "/", tot),"\n")

				# Run simulations  
				res  <-  runReplicateAutoInvSims(nReps = nReps, N = N, m = m.vals[k], s = s, h = h, r = r.vals[j], 
												 n = n, u = u, h.del = h.del, s.del = s.del.vals[l], noDel = FALSE, 
												 newMutant = newMutant)

				# Save data 
				rs      <-  rep(r.vals[j], times=nrow(res$results.df))
				ms      <-  rep(m.vals[k], times=nrow(res$results.df))
				s.dels  <-  rep(s.del.vals[l], times=nrow(res$results.df))

				# Append to data frame
				df      <-  cbind(res$results.df, rs, ms, s.dels)
				data    <-  rbind(data, df)
				rm(df)

			}
		}
	}

	# Include constant variables in data frame
	ss      <-  rep(s, times=nrow(data))
	hs      <-  rep(h, times=nrow(data))
	Ns      <-  rep(N, times=nrow(data))
	us      <-  rep(u, times=nrow(data))
	h.dels  <-  rep(h.del, times=nrow(data))
	data    <-  cbind(data, ss, hs, Ns, us, h.dels)
	colnames(data)  <-  c("finalInvFreq","finalE.InvFreq","finalW.mean",
						  "nGen","invEst","invEstTime","nDels","r","m",
						  "s.dels","s","h","N","u","h.del")

	# create file name
	filename  <-  paste("./output/data/simResults/auto-InvSimsData", "_s", s, "_h", h, "_N", N, "_n", n, "_u", u, "_hDel", h.del, ".csv", sep="")

	# export data as .csv to ./output/data
	write.csv(data, file=filename, row.names = FALSE)

	#  Return results in case user wants it
	return(data)
	
}