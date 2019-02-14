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
#' @param Fm     2 x 3 Matrix of adult male genotypic frequencies
#' @param W_y.g  Fitness expressions for male gametes (i.e., fitness effects associated with the wild-type and mutant Y chromosomes)
#' @export
maleGameteSelAuto  <-  function(Fm, W_y.g) {
	(Fm * rep(c(1, W_y.g), 3)) / sum((Fm * rep(c(1, W_y.g), 3)))
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
offspringSelAuto  <-  function(Fm.g, Ff, W_c, W_o) {
				O  <-  matrix(c((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(Fm.g[1,1] + Fm.g[1,2]/2), 
							(Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(Fm.g[1,3] + Fm.g[1,2]/2) +
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*(Fm.g[1,1] + Fm.g[1,2]/2), 
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*(Fm.g[1,3] + Fm.g[1,2]/2),
							(Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(Fm.g[2,1] + Fm.g[2,2]/2), 
							(Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(Fm.g[2,3] + Fm.g[2,2]/2) +
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*(Fm.g[2,1] + Fm.g[2,2]/2), 
							(Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*(Fm.g[2,3] + Fm.g[2,2]/2)
							), nrow=2, byrow=TRUE)
				O/sum(O)
}

# male genotypes:   c(Y|AA, Y|Aa, Y|aa, 
#					  y|AA, y|Aa, y|aa)
# female genotypes: c(  AA,   Aa,   aa)

offspringMutSelAuto  <-  function(Fm.g, Ff, W_c, W_o, uy, ua) {
				O  <-  matrix(c( 

							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ua)) * ((Fm.g[1,1] + (Fm.g[1,2]/2))*(1 - uy)*(1 - ua)), 

							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ua)) * ((Fm.g[1,1] + (Fm.g[1,2]/2))*(1 - uy)*ua) + 
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ua)       * ((Fm.g[1,1] + (Fm.g[1,2]/2))*(1 - uy)*(1 - ua)) + 
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ua)) *  (Fm.g[1,3] + (Fm.g[1,2]/2))*(1 - uy) +
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])           * ((Fm.g[1,1] + (Fm.g[1,2]/2))*(1 - ua))*(1 - uy), 

							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ua) * ((Fm.g[1,1] + (Fm.g[1,2]/2))*(1 - uy)*ua) + 
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ua) *  (Fm.g[1,3] + (Fm.g[1,2]/2))*(1 - uy) +
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])     * ((Fm.g[1,1] + (Fm.g[1,2]/2))*(1 - uy)*ua) +  
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])     *  (Fm.g[1,3] + (Fm.g[1,2]/2))*(1 - uy),
							####
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ua)) * ((Fm.g[1,1] + (Fm.g[1,2]/2))*uy*(1 - ua)) + 
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ua)) * ((Fm.g[2,1] + (Fm.g[2,2]/2))*(1 - ua)), 

							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ua)) * ((Fm.g[1,1] + (Fm.g[1,2]/2))*uy*ua) + 
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ua)       * ((Fm.g[1,1] + (Fm.g[1,2]/2))*uy*(1 - ua))  + 
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ua)) *  (Fm.g[1,3] + (Fm.g[1,2]/2))*uy +
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])           * ((Fm.g[1,1] + (Fm.g[1,2]/2))*uy*(1 - ua)) + 
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ua)) * ((Fm.g[2,1] + (Fm.g[2,2]/2))*ua) + 
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*ua)       * ((Fm.g[2,1] + (Fm.g[2,2]/2))*(1 - ua)) + 
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ua)) *  (Fm.g[2,3] + (Fm.g[2,2]/2)) +
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])           * ((Fm.g[2,1] + (Fm.g[2,2]/2))*(1 - ua)),

							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ua) * ((Fm.g[1,1] + (Fm.g[1,2]/2))*uy*ua) + 
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ua) * ((Fm.g[1,3] + (Fm.g[1,2]/2))*uy) +
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])     * ((Fm.g[1,1] + (Fm.g[1,2]/2))*uy*ua) +  
							 (Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])     * ((Fm.g[1,3] + (Fm.g[1,2]/2))*uy) +
							
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*ua) * ((Fm.g[2,1] + (Fm.g[2,2]/2))*ua) + 
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*ua) *  (Fm.g[2,3] + (Fm.g[2,2]/2)) +
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])     * ((Fm.g[2,1] + (Fm.g[2,2]/2))*ua) +  
							(Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])      *  (Fm.g[2,3] + (Fm.g[2,2]/2))
							), nrow=2, byrow=TRUE)
				O/sum(O)
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
#' 	# males: c(Y/AA, Y/Aa, Y/aa,
	#		   y/AA, y/Aa, y/aa)
	# females: c(AA, Aa, aa)

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
		if(yInvade & !aInvade) {
			qy.init  <-  2/N
			Fm  <-  matrix(c((1 - qy.init), 0, 0,
								  qy.init,  0, 0), 
						   nrow=2, byrow=TRUE)
			Ff  <-  c(1, 0, 0)
		}
		if(!yInvade & aInvade) {
			fMut     <-  rbinom(1, 1, 1/2)
			mMut     <-  (1 - fMut)
			qa.init  <-  2/N
			Fm  <-  matrix(c(0, 0, 0,
							 (1 - (mMut*qa.init)), mMut*qa.init, 0), 
						   nrow=2, byrow=TRUE)
			Ff  <-  c((1 - (fMut*qa.init)), fMut*qa.init, 0)
		}
		if(yInvade & aInvade) {
			qy.init  <-  qy.init
			qx.init  <-  Fii.init
			fMut     <-  rbinom(1, 1, 1/2)
			mMut     <-  (1 - fMut)
			dblMut     <-  rbinom(1, 1, qy.init)
			Fm  <-  matrix(c((1 - (qy.init + (1-dblMut)*mMut*qa.init)), (mMut*(1-dblMut)*qa.init), 0,
								  qy.init - (mMut*dblMut*qa.init),  (mMut*dblMut*qa.init), 0), 
						   nrow=2, byrow=TRUE)
			Ff  <-  c((1 - fMut*qa.init), fMut*qa.init, 0)
		}

	# Run simulation, saving frequencies for each generation
	if(saveTrajectories) {

		# Storage structures for individual simulation data
		qy.t   <-  rep(0, times=(4*N+1))
		qa.t   <-  rep(0, times=(4*N+1))
		E.qy   <-  rep(0, times=(4*N+1))
		E.qa   <-  rep(0, times=(4*N+1))
	
		# Initial frequency of mutant Y chromosome and compensatory allele
		qy.t[1]  <-  rowSums(Fm)[2]
		E.qy[1]  <-  rowSums(Fm)[2]
		Fii      <-  (colSums(Fm) + Ff)/2
		qa.t[1]  <-  Fii[3] + Fii[2]/2
		E.qa[1]  <-  Fii[3] + Fii[2]/2
		
		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t[gen] != 0 & qa.t[gen] < 1) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			Fm.g         <-  maleGameteSelAuto(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating, selection, mutation
			E.Fm         <-  offspringSelAuto(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			# 4) Expected frequencies in offspring, after selection
			E.Ff         <-  colSums(E.Fm)
			E.qy[gen+1]  <-  sum(E.Fm[2,])
			E.qa[gen+1]  <-  E.Ff[3] + E.Ff[2]/2
			# 5) Realized frequencies in adults
			Fm           <-  matrix(rmultinom(n=1, size=N/2, prob=E.Fm)/(N/2), ncol=3)
			qy.t[gen+1]  <-  sum(Fm[2,])
			Ff           <-  as.vector(rmultinom(n = 1, size = N/2, prob = E.Ff)/(N/2))
			Fii          <-  (colSums(Fm) + Ff)/2
			qa.t[gen+1]  <-  Fii[3] + Fii[2]/2

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

		# Initial frequency of mutant Y chromosome and compensatory allele
		qy.t  <-  rowSums(Fm)[2]
		E.qy  <-  rowSums(Fm)[2]
		Fii   <-  (colSums(Fm) + Ff)/2
		qa.t  <-  Fii[3] + Fii[2]/2
		E.qa  <-  Fii[3] + Fii[2]/2

		## Start forward simulation with new mutant Y chromosome
		gen  <-  1
		while(gen < (4*N) & qy.t != 0 & qy.t < qcrit.y & qa.t < qcrit.a & qa.t != qnull.a) {

			## Step through recursions:
			# 1) Calculate frequency of mutant Y among male gametes
			Fm.g  <-  maleGameteSelAuto(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
			E.Fm  <-  offspringSelAuto(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			# 4) Expected frequencies in offspring, after selection
			E.Ff  <-  colSums(E.Fm)
			E.qy  <-  sum(E.Fm[2,])
			E.qa  <-  E.Ff[3] + E.Ff[2]/2
			# 5) Realized frequencies in adults
			Fm    <-  matrix(rmultinom(n=1, size=N/2, prob=E.Fm)/(N/2), ncol=3)
			qy.t  <-  sum(Fm[2,])
			Ff    <-  as.vector(rmultinom(n = 1, size = N/2, prob = E.Ff)/(N/2))
			Fii   <-  (colSums(Fm) + Ff)/2
			qa.t  <-  Fii[3] + Fii[2]/2

			gen   <-  gen+1

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
					"nGen"         =  (gen-1),
					"qyInv"        =  qyInv,
					"qyTilde"      =  qyTilde,
					"qyTildeTime"  =  qyTildeTime,
					"qaInv"        =  qaInv
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
yAutoCoEvolCycle  <-  function(N = N, sm = sm, ho = ho, delta = delta, hc = hc, sc = sc,
							   uy = 1e-3, qy.init = "singleCopy",
							   ua = 2e-3, qa.init = 0, ...) {
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
	if(qy.init == "singleCopy") {
		qy.init  <-  2/N
	}
	Fm  <-  matrix(c((1 - qy.init), 0, 0,
					      qy.init,  0, 0), 
				   nrow=2, byrow=TRUE)
	if(qa.init != 0) {
		Ff  <-  c((1 - qa.init), 0, qa.init)
	} else {Ff  <-  c((1 - qa.init), qa.init, 0)}
	
	# Storage structures for individual simulation data
	qy.t  <-  c()
	qa.t  <-  c()
	E.qy  <-  c()
	E.qa  <-  c()

	# Initial frequency of mutant Y chromosome and compensatory allele
	qy.t[1]  <-  rowSums(Fm)[2]
	E.qy[1]  <-  rowSums(Fm)[2]
	Fii      <-  (colSums(Fm) + Ff)/2
	qa.t[1]  <-  Fii[3] + Fii[2]/2
	E.qa[1]  <-  Fii[3] + Fii[2]/2
		
	## Start forward simulation with new mutant Y chromosome
	gen  <-  1
	while(qy.t[gen] < 1 | qa.t[gen] < 1) {

		# 1) Calculate frequency of mutant Y among male gametes
		Fm.g         <-  maleGameteSelAuto(Fm=Fm, W_y.g=W_y.g)
		# 2) Expected frequencies in offspring after mating, selection, mutation
#		E.Fm         <-  offspringSelAuto(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
#		E.Fm         <-  offspringMutAuto(O=E.Fm, uy=uy, ua=ua)
		E.Fm         <-  offspringMutSelAuto(Fm.g = Fm.g, Ff=Ff, W_c=W_c, W_o=W_o, uy=uy, ua=ua)
		E.Ff         <-  colSums(E.Fm)
		E.qy[gen+1]  <-  rowSums(E.Fm)[2]
		E.qa[gen+1]  <-  E.Ff[3] + E.Ff[2]/2
		# 3) Realized frequencies in adults
		if(any(E.Fm < 0)) {
			browser()
		}
		Fm           <-  matrix(rmultinom(n=1, size=N/2, prob=E.Fm)/(N/2), ncol=3)
		qy.t[gen+1]  <-  rowSums(Fm)[2]
		Ff           <-  as.vector(rmultinom(n=1, size=N/2, prob = E.Ff)/(N/2))
		Fii          <-  (colSums(Fm) + Ff)/2
		qa.t[gen+1]  <-  Fii[3] + Fii[2]/2

		gen  <-  gen+1
	}

	# When did mutant alleles establish and fix, and when was qTilde crossed?
	Time         <-  seq_along(qy.t)
	yInvTime     <-  Time[qy.t >= qcrit.y][1]
	yFixTime     <-  Time[qy.t == 1][1]
	qyTildeTime  <-  Time[qy.t >= qyTilde][1]
	qaAtqyTilde  <-  qa.t[qyTildeTime]
	aInvTime     <-  Time[qa.t >= qcrit.a][1]
	aFixTime     <-  gen
	
	# Save  simulation data
	res  <-  list(
				"sm"           =  sm,
				"ho"           =  ho,
				"delta"        =  delta,
				"hc"           =  hc,
				"sc"           =  sc,
				"qcrit.y"      =  qcrit.y,
				"qcrit.a"      =  qcrit.a,
				"qyTilde"      =  qyTilde,
				"qy.t"         =  qy.t[1:gen-1],
				"E.qy"         =  E.qy[1:gen-1],
				"qa.t"         =  qa.t[1:gen-1],
				"E.qa"         =  E.qa[1:gen-1],
				"yInvTime"     =  yInvTime,
				"yFixTime"     =  yFixTime,
				"qyTildeTime"  =  qyTildeTime,
				"qaAtqyTilde"  =  qaAtqyTilde,
				"aInvTime"     =  aInvTime,
				"aFixTime"     =  aFixTime
				)
return(res)
}





#' Run replicate stochastic simulation to generate a data set for plotting
#' the probability of establishment/fixation for invaison of mutant y and
#' AUTOSOMAL compensatory mutations 
#' 
#' @title makeDataYAutoInv
#' @param reps       How many replicate simulations will be run
#' @param N          Population size
#' @param sm         Selection favouring mutant y male gametes
#' @param ho.vals    Values of dominace coefficients for selection against offspring resulting from 
#'                   matings between mutant males and wild-type females to be explored
#' @param delta      Difference between sm and selection coefficient for selection against offspring 
#'                   resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                   (delta = sm - so)
#' @param hc.vals    Values of dominace coefficients for selection against offspring resulting from 
#' 					 matings between wild-type males and mutant females to be explored
#' @param sc        Cost of compensation selection coefficient against offspring resulting from matings
#'                  between wild-type fathers and mothers with mutant comensatory allele
#' @export
#' @seealso 
#' @author Colin Olito
makeDataYAutoInv  <-  function(reps = reps, N = N, sm = sm, 
							   ho.vals = ho.vals, delta = delta, 
							   hc.vals = hc.vals, sc = sc) {

	len.h  <-  length(ho.vals)

	pInvade_y  <-  c()
	pInvade_a  <-  c()
	qCrit_y    <-  c()
	qCrit_a    <-  c()
	qyTilde    <-  c()

	print('Running Replicate Simulations for Autosomal Model')
	for(i in 1:len.h) {
		nyInvade  <-  0
		naInvade  <-  0
		cat("\n")
		print(paste("Running dominance value", i, "/",len.h))
		pb   <-  txtProgressBar(min=0, max=reps, style=3)
		setTxtProgressBar(pb, 0)

		for(n in 1:reps) {
			yInvData  <-  yAutoInvadeFwdSim(N = N, sm = sm, 
											ho = ho.vals[i], delta = delta, 
											hc = hc.vals[i], sc = sc, 
											yInvade = TRUE, Fii.init = Fii.init,
											aInvade = FALSE, qy.init = 1, 
											saveTrajectories = FALSE)
			aInvData  <-  yAutoInvadeFwdSim(N = N, sm = sm, 
											ho = ho.vals[i], delta = delta, 
											hc = hc.vals[i], sc = sc, 
											yInvade = FALSE, Fii.init = Fii.init,
											aInvade = TRUE, qy.init = 1, 
											saveTrajectories = FALSE)
			nyInvade  <-  nyInvade + as.numeric(is.numeric(yInvData$qyInv))
			naInvade  <-  naInvade + as.numeric(is.numeric(aInvData$qaInv))
			setTxtProgressBar(pb, n)
		}
		pInvade_y[i]  <-  nyInvade/reps
		pInvade_a[i]  <-  naInvade/reps
		qCrit_y[i]    <-  yInvData$qcrit.y
		qCrit_a[i]    <-  aInvData$qcrit.a
		qyTilde[i]    <-  aInvData$qyTilde
	}
	# compile data as dataframe
	data  <-  data.frame(
						 "sm"         =  rep(sm, times=len.h),
						 "delta"      =  rep(delta, times=len.h),
						 "sc"         =  rep(sc, times=len.h),
						 "ho"         =  ho.vals,
						 "hc"         =  hc.vals,
						 "pInvade_y"  =  pInvade_y,
						 "pInvade_a"  =  pInvade_a,
						 "qCrit_y"    =  qCrit_y,
						 "qCrit_a"    =  qCrit_a,
						 "qyTilde"    =  qyTilde
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYAutoInv", "_sm", sm, "_delta", delta, "_sc", sc, "_N", N, "_reps", reps, ".csv", sep="")	
	write.csv(data, file=filename)
}





#' Run replicate stochastic simulation to generate a data set for plotting
#' the Time To fixation for invasion of mutant y and
#' X-LINKED compensatory mutations 
#' 
#' @title makeDataXLinkedInv
#' @title makeDataYAutoInv
#' @param reps       How many replicate simulations will be run
#' @param N          Population size
#' @param sm         Selection favouring mutant y male gametes
#' @param ho.vals    Values of dominace coefficients for selection against offspring resulting from 
#'                   matings between mutant males and wild-type females to be explored
#' @param delta      Difference between sm and selection coefficient for selection against offspring 
#'                   resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                   (delta = sm - so)
#' @param hc.vals    Values of dominace coefficients for selection against offspring resulting from 
#' 					 matings between wild-type males and mutant females to be explored
#' @param sc        Cost of compensation selection coefficient against offspring resulting from matings
#'                  between wild-type fathers and mothers with mutant comensatory allele
#' @export
#' @seealso
#' @author Colin Olito
makeDataYATimeFix  <-  function(reps = reps, N = N, sm = sm, 
									  ho.vals = ho.vals, delta = delta, 
									  hc.vals = hc.vals, sc = sc,
									  uy = 1e-3, qy.init = "singleCopy",
									  ua = 2e-3, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25)) {

	# length of dominance gradient being explored
	len.h  <-  length(ho.vals)

	# Define output variables
	tInvade_y    <-  c()
	tFix_y       <-  c()
	tInvade_a    <-  c()
	tFix_a       <-  c()
	qCrit_y      <-  c()
	qCrit_a      <-  c()
	qyTilde      <-  c()
	qyTildeTime  <-  c()
	deltaInv     <-  c()
	deltaFix     <-  c()
	deltaInvFix  <-  c()
	tCycle       <-  c()

	# Dominance gradient loop
	print('Running Y-A Coevolutionary Cycle Simulations')
	for(i in 1:len.h) {
		# Print progress to terminal
		cat("\n")
		print(paste("Running dominance value", i, "/",len.h))
		pb   <-  txtProgressBar(min=0, max=reps, style=3)
		setTxtProgressBar(pb, 0)
		
		# dummy variables for replicate simulation output
		tInvade_yTemp     <-  c()
		tFix_yTemp        <-  c()
		tInvade_aTemp     <-  c()
		tFix_aTemp        <-  c()
		deltaInv_Temp     <-  c()
		deltaFix_Temp     <-  c()
		deltaInvFix_Temp  <-  c()
		tCycle_Temp       <-  c()

		# Replicate simulation loop
		for(n in 1:reps) {
			coEvolCycleData  <-  yAutoCoEvolCycle(N = N, sm = sm, ho = ho.vals[i], 
													delta = delta, hc = hc.vals[i], sc = sc,
													uy = uy, qy.init = "singleCopy",
													ua = ua, qa.init = 0)

			tInvade_yTemp[n]  <-  coEvolCycleData$yInvTime
			tFix_yTemp[n]     <-  coEvolCycleData$yFixTime
			tInvade_aTemp[n]  <-  coEvolCycleData$aInvTime
			tFix_aTemp[n]     <-  coEvolCycleData$aFixTime
			deltaInv_Temp[n]     <-  tInvade_aTemp[n] - tInvade_yTemp[n]
			deltaFix_Temp[n]     <-  tFix_aTemp[n] - tFix_yTemp[n]
			deltaInvFix_Temp[n]  <-  tFix_aTemp[n] - tInvade_yTemp[n]
			
			# Back-calculate total time for complete coevlutionary cycles
			# from generation when compensatory mutation that ultimately fixed occurred
			yMutTime  <-  tInvade_yTemp[n]
			qy  <-  coEvolCycleData$qy.t[yMutTime]
			while(qy > 0) {
				yMutTime  <-  yMutTime - 1
				if(yMutTime == 0) {
					break
				}
				qy  <-  coEvolCycleData$qy.t[yMutTime]
			}
			tCycle_Temp[n]         <-  tFix_aTemp[n] - (yMutTime + 1)

			setTxtProgressBar(pb, n)
		}

		# output data for each level of dominance gradient, calculate necessary means
		tInvade_y[i]    <-  mean(tInvade_yTemp)
		tFix_y[i]       <-  mean(tFix_yTemp)
		tInvade_a[i]    <-  mean(tInvade_aTemp)
		tFix_a[i]       <-  mean(tFix_aTemp)
		qCrit_y[i]      <-  coEvolCycleData$qcrit.y
		qCrit_a[i]      <-  coEvolCycleData$qcrit.a
		qyTilde[i]      <-  coEvolCycleData$qyTilde
		qyTildeTime[i]  <-  coEvolCycleData$qyTildeTime
		deltaInv[i]     <-  mean(deltaInv_Temp)
		deltaFix[i]     <-  mean(deltaFix_Temp)
		deltaInvFix[i]  <-  mean(deltaInvFix_Temp)
		tCycle[i]       <-  mean(tCycle_Temp)

	}

	# compile output data as dataframe
	data  <-  data.frame(
						 "sm"           =  rep(sm,    times=len.h),
						 "delta"        =  rep(delta, times=len.h),
						 "sc"           =  rep(sc,    times=len.h),
						 "ho"           =  ho.vals,
						 "hc"           =  hc.vals,
						 "tInvade_y"    =  tInvade_y,
						 "tFix_y"       =  tFix_y,
						 "tInvade_a"    =  tInvade_a,
						 "tFix_a"       =  tFix_a,
						 "deltaInv"     =  deltaInv,
						 "deltaFix"     =  deltaFix,
						 "deltaInvFix"  =  deltaInvFix,
						 "tCycle"       =  tCycle,
						 "qCrit_y"      =  qCrit_y,
						 "qCrit_a"      =  qCrit_a,
						 "qyTilde"      =  qyTilde,
						 "qyTildeTime"  =  qyTildeTime
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYATimeFix", "_sm", sm, "_delta", delta, "_sc", sc, "_N", N, "_reps", reps, ".csv", sep="")	
	write.csv(data, file=filename)
}








#' Run replicate stochastic simulation to generate a data set for plotting
#' the Time To fixation for invasion of mutant y and
#' X-LINKED compensatory mutations 
#' 
#' @title makeDataXLinkedInv
#' @title makeDataYAutoInv
#' @param reps       How many replicate simulations will be run
#' @param N          Population size
#' @param sm         Selection favouring mutant y male gametes
#' @param ho.vals    Values of dominace coefficients for selection against offspring resulting from 
#'                   matings between mutant males and wild-type females to be explored
#' @param delta      Difference between sm and selection coefficient for selection against offspring 
#'                   resulting from matings between mutant y fathers and wild-type compensatory mothers.
#'                   (delta = sm - so)
#' @param hc.vals    Values of dominace coefficients for selection against offspring resulting from 
#' 					 matings between wild-type males and mutant females to be explored
#' @param sc        Cost of compensation selection coefficient against offspring resulting from matings
#'                  between wild-type fathers and mothers with mutant comensatory allele
#' @export
#' @seealso
#' @author Colin Olito
makeDataYATimeFixMutGrad  <-  function(reps = reps, N = N, sm = sm, 
									  ho = ho, delta = delta, 
									  hc = hc, sc = sc,
									  uy = 1e-3, qy.init = "singleCopy",
									  ua.vals = ua.vals, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25)) {

	# length of mutation rate gradient being explored
	len.ua  <-  length(ua.vals)

	# Define output variables
	tInvade_y    <-  c()
	tFix_y       <-  c()
	tInvade_a    <-  c()
	tFix_a       <-  c()
	qCrit_y      <-  c()
	qCrit_a      <-  c()
	qyTilde      <-  c()
	qyTildeTime  <-  c()
	deltaInv     <-  c()
	deltaFix     <-  c()
	deltaInvFix  <-  c()
	tCycle       <-  c()

	# Mutation rate gradient loop
	print('Running Y-A Coevolutionary Cycle Simulations')
	for(i in 1:len.ua) {
		# Print progress to terminal
		cat("\n")
		print(paste("Running ua value", i, "/",len.ua))
		pb   <-  txtProgressBar(min=0, max=reps, style=3)
		setTxtProgressBar(pb, 0)
		
		# dummy variables for replicate simulation output
		tInvade_yTemp     <-  c()
		tFix_yTemp        <-  c()
		tInvade_aTemp     <-  c()
		tFix_aTemp        <-  c()
		deltaInv_Temp     <-  c()
		deltaFix_Temp     <-  c()
		deltaInvFix_Temp  <-  c()
		tCycle_Temp       <-  c()

		# Replicate simulation loop
		for(n in 1:reps) {
			coEvolCycleData  <-  yAutoCoEvolCycle(N = N, sm = sm, ho = ho, 
													delta = delta, hc = hc, sc = sc,
													uy = uy, qy.init = "singleCopy",
													ua = ua.vals[i], qa.init = 0)

			tInvade_yTemp[n]     <-  coEvolCycleData$yInvTime
			tFix_yTemp[n]        <-  coEvolCycleData$yFixTime
			tInvade_aTemp[n]     <-  coEvolCycleData$aInvTime
			tFix_aTemp[n]        <-  coEvolCycleData$aFixTime
			deltaInv_Temp[n]     <-  tInvade_aTemp[n] - tInvade_yTemp[n]
			deltaFix_Temp[n]     <-  tFix_aTemp[n] - tFix_yTemp[n]
			deltaInvFix_Temp[n]  <-  tFix_aTemp[n] - tInvade_yTemp[n]
			
			# Back-calculate total time for complete coevlutionary cycles
			# from generation when compensatory mutation that ultimately fixed occurred
			yMutTime  <-  tInvade_yTemp[n]
			qy        <-  coEvolCycleData$qy.t[yMutTime]
			while(qy > 0) {
				yMutTime  <-  yMutTime - 1
				if(yMutTime == 0) {
					break
				}
				qy  <-  coEvolCycleData$qy.t[yMutTime]
			}
			tCycle_Temp[n]  <-  tFix_aTemp[n] - (yMutTime + 1)

			setTxtProgressBar(pb, n)
		}

		# output data for each level of mutation rate gradient, calculate necessary means
		tInvade_y[i]    <-  mean(tInvade_yTemp)
		tFix_y[i]       <-  mean(tFix_yTemp)
		tInvade_a[i]    <-  mean(tInvade_aTemp)
		tFix_a[i]       <-  mean(tFix_aTemp)
		qCrit_y[i]      <-  coEvolCycleData$qcrit.y
		qCrit_a[i]      <-  coEvolCycleData$qcrit.a
		qyTilde[i]      <-  coEvolCycleData$qyTilde
		qyTildeTime[i]  <-  coEvolCycleData$qyTildeTime
		deltaInv[i]     <-  mean(deltaInv_Temp)
		deltaFix[i]     <-  mean(deltaFix_Temp)
		deltaInvFix[i]  <-  mean(deltaInvFix_Temp)
		tCycle[i]       <-  mean(tCycle_Temp)

	}

	# compile output data as dataframe
	data  <-  data.frame(
						 "sm"           =  rep(sm,    times=len.ua),
						 "delta"        =  rep(delta, times=len.ua),
						 "sc"           =  rep(sc,    times=len.ua),
						 "ho"           =  rep(ho,    times=len.ua),
						 "hc"           =  rep(hc,    times=len.ua),
						 "uy"           =  rep(uy,    times=len.ua),
						 "ua"           =  ua.vals,
						 "tInvade_y"    =  tInvade_y,
						 "tFix_y"       =  tFix_y,
						 "tInvade_a"    =  tInvade_a,
						 "tFix_a"       =  tFix_a,
						 "deltaInv"     =  deltaInv,
						 "deltaFix"     =  deltaFix,
						 "deltaInvFix"  =  deltaInvFix,
						 "tCycle"       =  tCycle,
						 "qCrit_y"      =  qCrit_y,
						 "qCrit_a"      =  qCrit_a,
						 "qyTilde"      =  qyTilde,
						 "qyTildeTime"  =  qyTildeTime
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYATimeFixMutGrad2", "_sm", sm, "_delta", delta, "_ho", ho, "_sc", sc, "_hc", hc, "_N", N, "_reps", reps, ".csv", sep="")	
	write.csv(data, file=filename)
}
