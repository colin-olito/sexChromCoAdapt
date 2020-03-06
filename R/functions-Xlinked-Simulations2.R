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
maleGameteSelXLinked  <-  function(Fm, W_y.g) {
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
maleOffspringSelXLinked  <-  function(Fm.g, Ff, W_c, W_o) {
				Fm_o  <-  c((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(Fm.g[1] + Fm.g[3]), 
							(Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(Fm.g[2] + Fm.g[4]),
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*(Fm.g[1] + Fm.g[3]),
							(Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*(Fm.g[2] + Fm.g[4])
						 )
				Fm_o/sum(Fm_o)
}
maleOffspringMutXLinked  <-  function(O_m, uy, ux) {
						  c(O_m[1]*(1 - uy)*(1 - ux),
							O_m[1]*uy*(1 - ux) + O_m[2]*(1 - ux),
							O_m[1]*(1 - uy)*ux + O_m[3]*(1 - uy),
							O_m[1]*uy*ux + O_m[3]*uy + O_m[4]
						 )
}
femaleOffspringSelXLinked  <-  function(Fm.g, Ff, W_c, W_o) {
				Ff_o  <-  c((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*Fm.g[1] + (Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*Fm.g[2], 
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*Fm.g[1] + (Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*Fm.g[3] + 
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*Fm.g[2] + (Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*Fm.g[4], 
							(Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*Fm.g[3] + (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*Fm.g[4]
							)
				Ff_o/sum(Ff_o)
}
femaleOffspringMutXLinked  <-  function(O_f, ux) {
						  c(O_f[1]*(1 - ux),
							O_f[1]*ux + (O_f[2]/2) + (O_f[2]/2)*(1 - ux),
							(O_f[2]/2)*ux + O_f[3]
							)
}






	# males: c(XY, Xy, xY, xy)
	# females: c(XX, Xx, xx)
maleOffspringMutSelXLinked  <-  function(Fm.g, Ff, W_c, W_o, uy, ux) {
				Fm_o  <-  c(
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ux)) * ((Fm.g[1] + Fm.g[3])*(1 - uy)), 
							
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ux)) * ((Fm.g[1] + Fm.g[3])*uy) +
							((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ux)) * (Fm.g[2] + Fm.g[4]),
							
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ux)       * ((Fm.g[1] + Fm.g[3])*(1 - uy)) +
							((Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2])*(1 - ux)) * ((Fm.g[1] + Fm.g[3])*(1 - uy)),
							
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ux) * ((Fm.g[1] + Fm.g[3])*uy) +
							((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ux) * (Fm.g[2] + Fm.g[4]) + 
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])     * ((Fm.g[1] + Fm.g[3])*uy) +
							 (Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])     * (Fm.g[2] + Fm.g[4])
						 )
				Fm_o/sum(Fm_o)
}

	# males: c(XY, Xy, xY, xy)
	# females: c(XX, Xx, xx)


femaleOffspringMutSelXLinked  <-  function(Fm.g, Ff, W_c, W_o, ux) {
				Ff_o  <-  c(
							(((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ux)) * (Fm.g[1]*(1 - ux))) + 
							(((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ux)) * (Fm.g[2]*(1 - ux))), 

							(((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ux)       * (Fm.g[1]*(1 - ux))) + 
							(((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ux)) * (Fm.g[1]*ux)) + 
							(((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*ux)       * (Fm.g[2]*(1 - ux))) +
							(((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ux)) * (Fm.g[2]*ux)) +
							((Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2]) * (Fm.g[1]*(1 - ux))) + ((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*(1 - ux))*Fm.g[3] + 
							((Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2]) * (Fm.g[2]*(1 - ux))) + ((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*(1 - ux))*Fm.g[4], 

							(((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ux) * (Fm.g[1]*ux)) + 
							(((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*ux) * (Fm.g[2]*ux)) +
							((Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2]) * (Fm.g[1]*ux)) + ((Ff[1]*W_c[1] + (Ff[2]/2)*W_c[2])*ux)*Fm.g[3] + 
							((Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2]) * (Fm.g[2]*ux)) + ((Ff[1]*W_o[1] + (Ff[2]/2)*W_o[2])*ux)*Fm.g[4] +
							 ((Ff[3]*W_c[3] + (Ff[2]/2)*W_c[2]) * Fm.g[3]) + ((Ff[3]*W_o[3] + (Ff[2]/2)*W_o[2])*Fm.g[4])
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
			Fm.g         <-  maleGameteSelXLinked(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
			E.Fm         <-  maleOffspringSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			E.Ff         <-  femaleOffspringSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			# 4) Expected frequencies in offspring, after selection
			E.qy[gen+1]  <-  sum(E.Fm[c(2,4)])
			E.qx[gen+1]  <-  (sum(E.Fm[3:4]) + 2*sum(E.Ff[3] + E.Ff[2]/2))/3
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
			Fm.g  <-  maleGameteSelXLinked(Fm=Fm, W_y.g=W_y.g)
			# 2) Mating and selection
			E.Fm  <-  maleOffspringSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
			E.Ff  <-  femaleOffspringSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
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

	# males: c(XY, Xy, xY, xy)
	# females: c(XX, Xx, xx)

		# 1) Calculate frequency of mutant Y among male gametes
		Fm.g         <-  maleGameteSelXLinked(Fm=Fm, W_y.g=W_y.g)
		# 2) Mating, mutation, selection
#		E.Fm         <-  maleOffspringSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
#		E.Fm         <-  maleOffspringMutXLinked(O_m=E.Fm, uy=uy, ux=ux) 
		E.Fm         <-  maleOffspringMutSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o, uy=uy, ux=ux)
#		E.Ff         <-  femaleOffspringSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o)
#		E.Ff         <-  femaleOffspringMutXLinked(O_f=E.Ff, ux=ux)
		E.Ff         <-  femaleOffspringMutSelXLinked(Fm.g=Fm.g, Ff=Ff, W_c=W_c, W_o=W_o, ux=ux)
		# 3) Expected frequencies in offspring, after selection
		E.qy[gen+1]  <-  sum(E.Fm[c(2,4)])
		E.qx[gen+1]  <-  (sum(E.Fm[3:4]) + 2*sum(E.Ff[3] + E.Ff[2]/2))/3
		# 4) Realized frequencies in adults
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
				"qyTilde"      =  qyTilde,
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
makeDataYXLinkedInv  <-  function(reps = reps, N = N, sm = sm, 
								  ho.vals = ho.vals, delta = delta, 
								  hc.vals = hc.vals, sc = sc) {

	len.h  <-  length(ho.vals)

	pInvade_y  <-  c()
	pInvade_x  <-  c()
	qCrit_y    <-  c()
	qCrit_x    <-  c()
	qyTilde    <-  c()

	print('Running Replicate Simulations for X-linked Model')
	for(i in 1:len.h) {
	cat("\n")
	print(paste("Running dominance value", i, "/",len.h))
	pb   <-  txtProgressBar(min=0, max=reps, style=3)
	setTxtProgressBar(pb, 0)
		nyInvade  <-  0
		nxInvade  <-  0
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
											   yInvade = FALSE, Fii.init = Fii.init,
											   xInvade = TRUE, qy.init = 1, 
											   saveTrajectories = FALSE)
			nyInvade  <-  nyInvade + as.numeric(is.numeric(yInvData$qyInv))
			nxInvade  <-  nxInvade + as.numeric(is.numeric(xInvData$qxInv))
			setTxtProgressBar(pb, n)
		}
		pInvade_y[i]  <-  nyInvade/reps
		pInvade_x[i]  <-  nxInvade/reps
		qCrit_y[i]    <-  yInvData$qcrit.y
		qCrit_x[i]    <-  xInvData$qcrit.x
		qyTilde[i]    <-  xInvData$qyTilde
	}
	# compile data as dataframe
	data  <-  data.frame(
						 "sm"         =  rep(sm,    times=len.h),
						 "delta"      =  rep(delta, times=len.h),
						 "sc"         =  rep(sc,    times=len.h),
						 "ho"         =  ho.vals,
						 "hc"         =  hc.vals,
						 "pInvade_y"  =  pInvade_y,
						 "pInvade_x"  =  pInvade_x,
						 "qCrit_y"    =  qCrit_y,
						 "qCrit_x"    =  qCrit_x,
						 "qyTilde"    =  qyTilde
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYXLinkedInv", "_sm", sm, "_delta", delta, "_sc", sc, "_N", N, "_reps", reps, ".csv", sep="")	
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
makeDataYXTimeFix  <-  function(reps = reps, N = N, sm = sm, 
									  ho.vals = ho.vals, delta = delta, 
									  hc.vals = hc.vals, sc = sc,
									  uy = 1e-3, qy.init = "singleCopy",
									  ux = 1.5e-3, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25)) {

	# length of dominance gradient being explored
	len.h  <-  length(ho.vals)

	# Define output variables
	tInvade_y    <-  c()
	tFix_y       <-  c()
	tInvade_x    <-  c()
	tFix_x       <-  c()
	qCrit_y      <-  c()
	qCrit_x      <-  c()
	qyTilde      <-  c()
	qyTildeTime  <-  c()
	deltaInv     <-  c()
	deltaFix     <-  c()
	deltaInvFix  <-  c()
	tCycle       <-  c()

	# Dominance gradient loop
	print('Running Y-X Coevolutionary Cycle Simulations')
	for(i in 1:len.h) {
		# Print progress to terminal
		cat("\n")
		print(paste("Running dominance value", i, "/",len.h))
		pb   <-  txtProgressBar(min=0, max=reps, style=3)
		setTxtProgressBar(pb, 0)
		
		# dummy variables for replicate simulation output
		tInvade_yTemp     <-  c()
		tFix_yTemp        <-  c()
		tInvade_xTemp     <-  c()
		tFix_xTemp        <-  c()
		deltaInv_Temp     <-  c()
		deltaFix_Temp     <-  c()
		deltaInvFix_Temp  <-  c()
		tCycle_Temp       <-  c()

		# Replicate simulation loop
		for(n in 1:reps) {
			coEvolCycleData  <-  yXLinkedCoEvolCycle(N = N, sm = sm, ho = ho.vals[i], 
													 delta = delta, hc = hc.vals[i], sc = sc,
													 uy = uy, qy.init = "singleCopy",
													 ux = ux, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25))
			tInvade_yTemp[n]     <-  coEvolCycleData$yInvTime
			tFix_yTemp[n]        <-  coEvolCycleData$yFixTime
			tInvade_xTemp[n]     <-  coEvolCycleData$xInvTime
			tFix_xTemp[n]        <-  coEvolCycleData$xFixTime
			deltaInv_Temp[n]     <-  tInvade_xTemp[n] - tInvade_yTemp[n]
			deltaFix_Temp[n]     <-  tFix_xTemp[n] - tFix_yTemp[n]
			deltaInvFix_Temp[n]  <-  tFix_xTemp[n] - tInvade_yTemp[n]
			
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
			tCycle_Temp[n]         <-  tFix_xTemp[n] - (yMutTime + 1)
			setTxtProgressBar(pb, n)			
		}

		# output data for each level of dominance gradient, calculate necessary means
		tInvade_y[i]    <-  mean(tInvade_yTemp)
		tFix_y[i]       <-  mean(tFix_yTemp)
		tInvade_x[i]    <-  mean(tInvade_xTemp)
		tFix_x[i]       <-  mean(tFix_xTemp)
		qCrit_y[i]      <-  coEvolCycleData$qcrit.y
		qCrit_x[i]      <-  coEvolCycleData$qcrit.x
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
						 "tInvade_x"    =  tInvade_x,
						 "tFix_x"       =  tFix_x,
						 "deltaInv"     =  deltaInv,
						 "deltaFix"     =  deltaFix,
						 "deltaInvFix"  =  deltaInvFix,
						 "tCycle"       =  tCycle,
						 "qCrit_y"      =  qCrit_y,
						 "qCrit_x"      =  qCrit_x,
						 "qyTilde"      =  qyTilde,
						 "qyTildeTime"  =  qyTildeTime
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYXTimeFix", "_sm", sm, "_delta", delta, "_sc", sc, "_N", N, "_reps", reps, ".csv", sep="")	
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
makeDataYXTimeFixMutGrad  <-  function(reps = reps, N = N, sm = sm, 
									  ho = ho, delta = delta, 
									  hc = hc, sc = sc,
									  uy = 1e-3, qy.init = "singleCopy",
									  ux.vals = ux.vals, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25)) {

	# length of mutation rate gradient being explored
	len.ux  <-  length(ux.vals)

    # Calculate weighted mutation rate from
    # relative size of X vs. Autosomal genome in 
    # D. melanogaster (see Adams et al. 2000)
    MbA  <-  23 + 5.4 + 11 + 21.4 + 24.4 + 8.2 + 8.2 + 28 + 3.1 + 1.2
    MbX  <-  20 + 21.8
    MbY  <-  40.9
    weightMuA    <-  MbA / (MbA + MbX)
    weightMuX    <-  MbX / (MbA + MbX)

	ux.vals  <-  ux.vals * weightMuX

	# Define output variables
	tInvade_y    <-  c()
	tFix_y       <-  c()
	tInvade_x    <-  c()
	tFix_x       <-  c()
	qCrit_y      <-  c()
	qCrit_x      <-  c()
	qyTilde      <-  c()
	qyTildeTime  <-  c()
	deltaInv     <-  c()
	deltaFix     <-  c()
	deltaInvFix  <-  c()
	tCycle       <-  c()

	# Mutation rate gradient loop
	print('Running Y-X Coevolutionary Cycle Simulations')
	for(i in 1:len.ux) {

		# Print progress to terminal
		cat("\n")
		print(paste("Running u_x value", i, "/",len.ux))
		pb   <-  txtProgressBar(min=0, max=reps, style=3)
		setTxtProgressBar(pb, 0)
		
		# dummy variables for replicate simulation output
		tInvade_yTemp     <-  c()
		tFix_yTemp        <-  c()
		tInvade_xTemp     <-  c()
		tFix_xTemp        <-  c()
		deltaInv_Temp     <-  c()
		deltaFix_Temp     <-  c()
		deltaInvFix_Temp  <-  c()
		tCycle_Temp       <-  c()

		# Replicate simulation loop
		for(n in 1:reps) {
			coEvolCycleData  <-  yXLinkedCoEvolCycle(N = N, sm = sm, ho = ho, 
													 delta = delta, hc = hc, sc = sc,
													 uy = uy, qy.init = "singleCopy",
													 ux = ux.vals[i], qx.init = 0, Ff.init = c(0.25, 0.5, 0.25))
			tInvade_yTemp[n]     <-  coEvolCycleData$yInvTime
			tFix_yTemp[n]        <-  coEvolCycleData$yFixTime
			tInvade_xTemp[n]     <-  coEvolCycleData$xInvTime
			tFix_xTemp[n]        <-  coEvolCycleData$xFixTime
			deltaInv_Temp[n]     <-  tInvade_xTemp[n] - tInvade_yTemp[n]
			deltaFix_Temp[n]     <-  tFix_xTemp[n] - tFix_yTemp[n]
			deltaInvFix_Temp[n]  <-  tFix_xTemp[n] - tInvade_yTemp[n]
			
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
			tCycle_Temp[n]         <-  tFix_xTemp[n] - (yMutTime + 1)
			setTxtProgressBar(pb, n)			
		}

		# output data for each value of mutation rate gradient, calculate necessary means
		tInvade_y[i]    <-  mean(tInvade_yTemp)
		tFix_y[i]       <-  mean(tFix_yTemp)
		tInvade_x[i]    <-  mean(tInvade_xTemp)
		tFix_x[i]       <-  mean(tFix_xTemp)
		qCrit_y[i]      <-  coEvolCycleData$qcrit.y
		qCrit_x[i]      <-  coEvolCycleData$qcrit.x
		qyTilde[i]      <-  coEvolCycleData$qyTilde
		qyTildeTime[i]  <-  coEvolCycleData$qyTildeTime
		deltaInv[i]     <-  mean(deltaInv_Temp)
		deltaFix[i]     <-  mean(deltaFix_Temp)
		deltaInvFix[i]  <-  mean(deltaInvFix_Temp)
		tCycle[i]       <-  mean(tCycle_Temp)
	}
	# compile data as dataframe
	data  <-  data.frame(
						 "sm"           =  rep(sm,    times=len.ux),
						 "delta"        =  rep(delta, times=len.ux),
						 "sc"           =  rep(sc,    times=len.ux),
						 "ho"           =  rep(ho, times=len.ux),
						 "hc"           =  rep(hc, times=len.ux),
						 "uy"           =  rep(uy, times=len.ux),
						 "ux"           =  ux.vals,
						 "tInvade_y"    =  tInvade_y,
						 "tFix_y"       =  tFix_y,
						 "tInvade_x"    =  tInvade_x,
						 "tFix_x"       =  tFix_x,
						 "deltaInv"     =  deltaInv,
						 "deltaFix"     =  deltaFix,
						 "deltaInvFix"  =  deltaInvFix,
						 "tCycle"       =  tCycle,
						 "qCrit_y"      =  qCrit_y,
						 "qCrit_x"      =  qCrit_x,
						 "qyTilde"      =  qyTilde,
						 "qyTildeTime"  =  qyTildeTime
						)
	# Write data to file
	filename  <-  paste("./output/data/simData/dataYXTimeFixMutGrad2_weighted", "_sm", sm, "_delta", delta, "_ho", ho, "_sc", sc, "_hc", hc, "_N", N, "_reps", reps, ".csv", sep="")	
	write.csv(data, file=filename)
}







