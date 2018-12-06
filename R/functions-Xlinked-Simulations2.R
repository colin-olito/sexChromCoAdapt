########################################################
#  Antagonistic coevolution between the sex-chromosomes
#
#  Necessary functions for stochastic simulation of 
#  invasion of Y-linked male-benefit alleles and  
#  X-LINKED compensatory alleles. 
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
#' @param Fm     2 x 3 Matrix of adult male genotypic frequencies
#' @param W_y.g  Fitness expressions for male gametes (i.e., fitness effects associated with the wild-type and mutant Y chromosomes)
#' @export
maleGameteSel  <-  function(Fm, W_y.g) {
	(Fm * c(1, W_y.g, 1, W_y.g)) / sum((Fm * c(1, W_y.g, 1, W_y.g)))
}


#' Frequency of offspring genotypes after mating and selection
#'
#' @title Haplotype frequencies among gametes
#' @param Fm.g   2 x 3 Matrix of adult male genotypic frequencies at mating
#' @param Ff     Vector with length = 3 of adult female genotypic frequencies at mating
#' @param W_c    Fitness expressions for offspring survival resulting from matings with 
#' 				 wild-type males (vector with length = 3)
#' @param W_o    Fitness expressions for offspring survival resulting from matings with 
#' 				 mutant males (vector with length = 3)
#' @export
	# males: c(XY, Xy, xY, xy)
	# females: c(XX, Xx, xx)
maleOffspringSel  <-  function(Fm.g, Ff, W_c, W_o) {
				Fm_o  <-  c((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(Fm.g[1] + Fm.g[3]), 
							(Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(Fm.g[2] + Fm.g[4]),
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*(Fm.g[1] + Fm.g[3]),
							(Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*(Fm.g[2] + Fm.g[4])
						 )
				Fm_o/sum(Fm_o)
}
femaleOffspringSel  <-  function(Fm.g, Ff, W_c, W_o) {
				Ff_o  <-  c((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*Fm.g[1] + (Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*Fm.g[2], 
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*Fm.g[1] + (Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*Fm.g[3] + 
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*Fm.g[2] + (Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*Fm.g[4], 
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*Fm.g[3] + (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*Fm.g[4]
							)
				Ff_o/sum(Ff_o)
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

qyTilde.X  <-  function(sm, ho, so, hc, sc) {
	-((2*hc*sc - hc*sc*sm + (sm^2) + 2*so - 2*ho*so + hc*sc*so + sm*so - 
	 3*ho*sm*so + hc*sc*sm*so - (sm^2)*so - ho*(sm^2)*so + ho*(so^2) + 2*ho*sm*(so^2) + ho*(sm^2)*(so^2) - 
		sqrt(8*hc*sc*(sm*(so - 1) + so)*(-hc*sc + (-2 + ho)*so + sm*(1 - 2*so + ho*so)) + ((2 + ho*(-2 + so))*so + hc*sc*(2 + sm*(-1 + so) + so) + (sm^2)*(-1 + so)*(-1 + ho*so) + sm*so*(1 + ho*(-3 + 2*so)))^2)) / 
			(2*(sm*(-1 + so) + so)*(-hc*sc + (-2 + ho)*so + sm*(1 - 2*so + ho*so))))
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
#' @param xInvade   Logical. If TRUE, simulates invasion from single copy of mutant a compensatory 
#'                  allele into population initially fixed for the mutant y chromosome.
#' @param qy.init   Initial frequency of mutant y. Defaults to 1, but can be set to any frequency.
#' @param saveTrajectories  Save evolutionary trajectories of allele frequencies? Used for troubleshooting.
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
yXLinkedInvadeFwdSim  <-  function(N = N, sm = sm, delta = delta, ho = ho, sc = sc, hc = hc, 
								   yInvade = TRUE, Fii.init = Fii.init,
								   xInvade = FALSE, qy.init = 1, 
								   saveTrajectories = FALSE, ...) {
	# Pre-emptive warnings
	if(!yInvade & !xInvade) {
		stop('Setting saveTrajectories = FALSE makes no sense if neither mutant allele is invading')
	}

	# Storage for result objects
	qyInv        <-  NA
	qyInvTime    <-  NA
	qyTildeTime  <-  NA
	qxInv        <-  NA

	# Define Fitness expressions
	so     <-  (sm - delta)
	W_y.g  <-  (1 + sm)
	W_o    <-  c((1 - so), (1 - ho*so), 1)
	W_c    <-  c(1, (1 - hc*sc), (1 - sc))

	# Define threshold frequency for establishment of mutant Y
	qcrit.y  <-  8/(N*delta)
	qcrit.x  <-  8/(3*N*so)
	qyTilde  <-  qyTilde.X(sm = sm, ho = ho, so = so, hc = hc, sc = sc)

	# set initial genotype frequencies
	# males: c(XY, Xy, xY, xy)
	# females: c(XX, Xx, xx)
		if(yInvade & !xInvade) {
			qy.init  <-  2/N
			Fm  <-  c((1 - qy.init), qy.init, 0, 0)
			Ff  <-  c(1, 0, 0)
		}
		if(!yInvade & xInvade) {
			fMut     <-  rbinom(1, 1, 2/3)
			mMut     <-  (1 - fMut)
			qx.init  <-  2/N
			Fm  <-  c(0, (1 - (mMut*qx.init)), 0, mMut*qx.init)
			Ff  <-  c((1 - (fMut*qx.init)), fMut*qx.init, 0)
		}
		if(yInvade & xInvade) {
			qy.init  <-  qy.init
			qx.init  <-  Fii.init
			fMut     <-  rbinom(1, 1, 2/3)
			mMut     <-  (1 - fMut)
			dblMut   <-  rbinom(1, 1, qy.init)
			Fm       <-  c((1 - (qy.init + (1-dblMut)*mMut*qx.init)), 
						   qy.init - (mMut*dblMut*qx.init),  
						   (mMut*(1-dblMut)*qx.init), 
						   (mMut*dblMut*qx.init))
			Ff  <-  c((1 - fMut*qx.init), fMut*qx.init, 0)
		}

	# Run simulation, saving frequencies for each generation
	if(saveTrajectories) {

		# Storage structures for individual simulation data
		qy.t   <-  rep(0, times=(4*N+1))
		qx.t   <-  rep(0, times=(4*N+1))
		E.qy   <-  rep(0, times=(4*N+1))
		E.qx   <-  rep(0, times=(4*N+1))
	
		# Initial frequency of mutant Y chromosome and compensatory allele
		qy.t[1]  <-  sum(Fm[c(2,4)])
		E.qy[1]  <-  sum(Fm[c(2,4)])
		qx.t[1]  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3
		E.qx[1]  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3
		
		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t[gen] != 0 & qx.t[gen] < 1) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			Fm.g  <-  maleGameteSel(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
			E.Fm     <-  maleOffspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			E.Ff     <-  femaleOffspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			# 4) Expected frequencies in offspring, after selection
			E.qy[gen+1]     <-  sum(E.Fm[c(2,4)])
			E.qx[gen+1]     <-  (sum(E.Fm[3:4]) + 2*sum(E.Ff[3] + E.Ff[2]/2))/3
			# 5) Realized frequencies in adults
			Fm           <-  as.vector(rmultinom(n=1, size=N/2, prob=E.Fm)/(N/2))
			qy.t[gen+1]  <-  sum(Fm[c(2,4)])
			Ff           <-  as.vector(rmultinom(n = 1, size = N/2, prob = E.Ff)/(N/2))
			qx.t[gen+1]  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3

			gen  <-  gen+1
		}

		# Remove unnecessary 0's
		qy.t  <-  qy.t[qy.t != 0]
		E.qy  <-  E.qy[1:length(qy.t)]
		qx.t  <-  qx.t[1:length(qy.t)]
		E.qx  <-  E.qx[1:length(qy.t)]

		# Have the mutant alleles reached threshold frequencies for establishment?
		# When?

		if(any(qy.t >= qcrit.y)) {
			qyInv        <-  1
			qyTildeTime  <-  gen[qy.t >= qyTilde][1]
		}
		if(any(qx.t >= qcrit.x)) {
			qxInv      <-  1
		}
	
		# Save  simulation data
		res  <-  list(
					"yInvade"      =  yInvade,
					"xInvade"      =  xInvade,
					"delta"        =  delta,
					"qcrit.y"      =  qcrit.y,
					"qcrit.x"      =  qcrit.x,
					"qyFreq"       =  qy.t[1:gen-1],
					"E.qy"         =  E.qy[1:gen-1],
					"qxFreq"       =  qx.t[1:gen-1],
					"E.qx"         =  E.qx[1:gen-1],
					"nGen"         =  gen,
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qxInv"        =  qxInv
					)
	} 

	if(!saveTrajectories) {

		# hotfix to ensure while loop works when 
		# qy.init > qcrit.y
		if(xInvade) {
			qcrit.y  <-  1.1
			qnull.a  <-  0
		} else {qnull.a  <-  1}

		# Initial frequency of mutant Y chromosome and compensatory allele
		qy.t  <-  sum(Fm[c(2,4)])
		E.qy  <-  sum(Fm[c(2,4)])
		qx.t  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3
		E.qx  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3
		

		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t != 0 & qy.t < qcrit.y & qx.t < qcrit.x & qx.t != qnull.a) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			Fm.g  <-  maleGameteSel(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
			E.Fm     <-  maleOffspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			E.Ff     <-  femaleOffspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			# 4) Expected frequencies in offspring, after selection
			E.qy  <-  sum(E.Fm[c(2,4)])
			E.qx  <-  (sum(E.Fm[3:4]) + 2*sum(E.Ff[3] + E.Ff[2]/2))/3
			# 5) Realized frequencies in adults
			Fm    <-  as.vector(rmultinom(n=1, size=N/2, prob=E.Fm)/(N/2))
			qy.t  <-  sum(Fm[c(2,4)])
			Ff    <-  as.vector(rmultinom(n = 1, size = N/2, prob = E.Ff)/(N/2))
			qx.t  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3

			gen  <-  gen+1

			# Has the inversion reached threshold frequency for establishment (pcrit)? 
			# When did it first reach pcrit?
			if(qy.t >= qcrit.y) {
				qyInv        <-  1
				qyTildeTime  <-  gen[qy.t >= qyTilde][1]
			}
			if(qx.t >= qcrit.x) {
				qxInv      <-  1
			}
		}	

		# Save  simulation data
		res  <-  list(
					"yInvade"      =  yInvade,
					"xInvade"      =  xInvade,
					"delta"        =  delta,
					"qcrit.y"      =  qcrit.y,
					"qcrit.x"      =  qcrit.x,
					"qyFreq"       =  qy.t,
					"E.qy"         =  E.qy,
					"qxFreq"       =  qx.t,
					"E.qx"         =  E.qx,
					"nGen"         =  (gen-1),
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qxInv"        =  qxInv
					)
	} 
return(res)

}










#' Run a single stochastic simulation of a complete co-evolutionary cycle
#' between the mutant Y chromosomes and AUTOSOMAL compensatory locus
#' By default, a coevolutionary cycle is initiated by a single-copy mutation
#' at the male-beneficial Y-linked locus (this can be changed by changing qy.init)
#' and completes when the mutant compensatory allele (a) fixes in the population.
#' Both alleles experience recurrent one-way mutation (we assume no back-mutation
#' at either locus) 
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
#' @param uy        Per-generation mutation rate to male-beneficial antagonistic y chromosome (Y --> y) .
#' @param qy.init   Initial frequency of mutant y chromosome. Defaults to "singleCopy": q.init = 2/N
#' @param ua        Per-generation mutation rate to the compensatory allele at the AUTOSOMAL locus
#'                   (A --> a) .
#' @param qa.init   Initial frequency of mutant m allele at mitochondrial comensatory locus. 
#'                   Defaults to 0.
#' @export
#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
yXLinkedCoEvolCycle  <-  function(N = N, sm = sm, ho = ho, delta = delta, hc = hc, sc = sc,
							   uy = 1e-3, qy.init = "singleCopy",
							   ux = 3e-3, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25), ...) {

	# Define Fitness expressions
	so     <-  (sm - delta)
	W_y.g  <-  (1 + sm)
	W_o    <-  c((1 - so), (1 - ho*so), 1)
	W_c    <-  c(1, (1 - hc*sc), (1 - sc))

	# Define threshold frequency for establishment of mutant Y
	qcrit.y  <-  8/(N*delta)
	qcrit.x  <-  8/(3*N*so)
	qyTilde  <-  qyTilde.X(sm = sm, ho = ho, so = so, hc = hc, sc = sc)

	# set initial genotype frequencies
	if(qy.init == "singleCopy") {
		qy.init  <-  2/N
	}
	Fm  <-  c((1 - qy.init), qy.init, 0, 0)
	if(qx.init != 0) {
		Ff  <-  Ff.init
	} else {Ff  <-  c(1,0,0)}
	
	# Storage structures for individual simulation data
	qy.t  <-  c()
	qx.t  <-  c()
	E.qy  <-  c()
	E.qx  <-  c()

	# Initial frequency of mutant Y chromosome and compensatory allele
	qy.t[1]  <-  sum(Fm[c(2,4)])
	E.qy[1]  <-  sum(Fm[c(2,4)])
	qx.t[1]  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3
	E.qx[1]  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3
		
	## Start forward simulation with new mutant Y chromosome
	gen  <-  1
	while(qy.t[gen] < 1 | qx.t[gen] < 1) {

		## Step through recursions:
		# 1) Mutation
		mutateY  <-  runif(1) <= uy
		mutateX  <-  runif(1) <= ux
		if(mutateY & qy.t[gen] < 1) {
			mutant  <-  c(1:4)[as.vector(rmultinom(n=1, size=1, prob=Fm)) == 1]
			if(any(mutant == c(1,3))) {
				Fm[mutant]      <-  Fm[mutant] - (2/N)
				Fm[mutant + 1]  <-  Fm[mutant + 1] + (2/N)
			}
		}
		if(mutateX) {
			fMut     <-  rbinom(1, 1, 2/3)
			if(fMut == 1) {
				prob    <-  Ff*c(1,1/2,0)/sum(Ff*c(1,1/2,0))
				mutant  <-  c(1:3)[as.vector(rmultinom(n=1, size=1, prob=prob)) == 1]
				Ff[mutant]      <-  Ff[mutant] - (2/N)
				Ff[mutant + 1]  <-  Ff[mutant + 1] + (2/N)
			}
			if(fMut == 0) {
				mutant  <-  c(1:4)[as.vector(rmultinom(n=1, size=1, prob=Fm)) == 1]
				if(any(mutant == c(1,2))) {
					Fm[mutant]      <-  Fm[mutant] - (2/N)
					Fm[mutant + 2]  <-  Fm[mutant + 1] + (2/N)
				}
			}
		}
		# 2) Calculate frequency of mutant Y among male gametes
		Fm.g  <-  maleGameteSel(Fm=Fm, W_y.g=W_y.g)
		# 3) Mating and selection
		E.Fm     <-  maleOffspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
		E.Ff     <-  femaleOffspringSel(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
		# 4) Expected frequencies in offspring, after selection
		E.qy[gen+1]     <-  sum(E.Fm[c(2,4)])
		E.qx[gen+1]     <-  (sum(E.Fm[3:4]) + 2*sum(E.Ff[3] + E.Ff[2]/2))/3
		# 5) Realized frequencies in adults
		Fm           <-  as.vector(rmultinom(n=1, size=N/2, prob=E.Fm)/(N/2))
		qy.t[gen+1]  <-  sum(Fm[c(2,4)])
		Ff           <-  as.vector(rmultinom(n = 1, size = N/2, prob = E.Ff)/(N/2))
		qx.t[gen+1]  <-  (sum(Fm[3:4]) + 2*sum(Ff[3] + Ff[2]/2))/3

		gen  <-  gen+1
	}

	# When did mutant alleles establish and fix, and when was qTilde crossed?
	Time         <-  seq_along(qy.t)
	yInvTime     <-  Time[qy.t >= qcrit.y][1]
	yFixTime     <-  Time[qy.t == 1][1]
	qyTildeTime  <-  Time[qy.t >= qyTilde][1]
	qxAtqyTilde  <-  qx.t[qyTildeTime]
	xInvTime     <-  Time[qx.t >= qcrit.x][1]
	xFixTime     <-  gen
	
	# Save  simulation data
	res  <-  list(
				"sm"           =  sm,
				"ho"           =  ho,
				"delta"        =  delta,
				"hc"           =  hc,
				"sc"           =  sc,
				"qcrit.y"      =  qcrit.y,
				"qcrit.x"      =  qcrit.x,
				"qy.t"         =  qy.t[1:gen-1],
				"E.qy"         =  E.qy[1:gen-1],
				"qx.t"         =  qx.t[1:gen-1],
				"E.qx"         =  E.qx[1:gen-1],
				"yInvTime"     =  yInvTime,
				"yFixTime"     =  yFixTime,
				"qyTildeTime"  =  qyTildeTime,
				"qxAtqyTilde"  =  qxAtqyTilde,
				"xInvTime"     =  xInvTime,
				"xFixTime"     =  xFixTime
				)
return(res)
}



#' Run replicate stochastic simulation to generate a data set for plotting
#' the probability of establishment/fixation for invaison of mutant y and
#' X-LINKED compensatory mutations 
#' 
#' @title makeDataXLinkedInv
#' @param N          Population size
#' @param sm         Selection favouring mutant y male gametes
#' @param delta.vals Difference between sm and selection coefficient for selection against offspring 
#'                   resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                   (delta = sm - so)
#' @param sc.vals    Cost of compensation selection coefficient against offspring resulting from matings
#'                   between wild-type fathers and mothers with mutant comensatory allele
#' @export#' @seealso `offFreq`, `findEqFreqs`, `x.1`, ...
#' @author Colin Olito
makeDataYAutoInv  <-  function(reps = reps, N = N, sm = sm, 
							   ho = ho.vals, delta = delta, 
							   hc = hc.vals, sc = sc) {

	len.h  <-  length(ho.vals)

	pInvade_y  <-  c()
	pInvade_x  <-  c()
	qCrit_y    <-  c()
	qCrit_x    <-  c()
	qyTilde    <-  c()

	for(i in 1:len.h) {
		nyInvade  <-  0
		naInvade  <-  0
		for(n in 1:reps) {
			yInvData  <-  yXLinkedInvadeFwdSim(N = N, sm = sm, 
											   ho = ho.vals[i], delta = delta, 
											   hc = hc.vals[i], sc = sc, 
											   yInvade = TRUE, Fii.init = Fii.init,
											   aInvade = FALSE, qy.init = 1, 
											   saveTrajectories = FALSE)
			xInvData  <-  yXLinkedInvadeFwdSim(N = N, sm = sm, 
											   ho = ho.vals[i], delta = delta, 
											   hc = hc.vals[i], sc = sc, 
											   yInvade = TRUE, Fii.init = Fii.init,
											   aInvade = FALSE, qy.init = 1, 
											   saveTrajectories = FALSE)
			nyInvade  <-  nyInvade + as.numeric(is.numeric(yInvData$qyInv))
			nxInvade  <-  nxInvade + as.numeric(is.numeric(xInvData$qaInv))
		}
		pInvade_y[i]  <-  nyInvade/reps
		pInvade_x[i]  <-  nxInvade/reps
		qCrit_y[i]    <-  yInvData$qcrit.y
		qCrit_x[i]    <-  xInvData$qcrit.x
		qyTilde[i]    <-  xInvData$qyTilde
	}
	# compile data as dataframe
	data  <-  data.frame(
						 "sm"         = rep(sm, times=len.h),
						 "delta"      = rep(delta, times=len.h),
						 "sc"         = rep(sc, times=len.h),
						 "ho"         = ho.vals,
						 "hc"         =  hc.vals,
						 "pInvade_y"  = pInvade_y,
						 "pInvade_x"  =  pInvade_x,
						 "qCrit_y"    =  qCrit_y,
						 "qCrit_x"    =  qCrit_x,
						 "qyTilde"    =  qyTilde
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYXLinkedInv", "_sm", sm, "_delta", delta, "_sc", sc, "_reps", reps, ".csv", sep="")	
	write.csv(data, file=filename)
}
