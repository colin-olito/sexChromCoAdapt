#####################################################
#  Coadaptation between the sex-chromosomes
#
#  Necessary functions for analyses
#
#  Author: Colin Olito
#
#  NOTES:  
#          



###############################################################
## General Functions 
###############################################################

rm(list=ls())



#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


####################
#  Simulations to check approximation for change in frequency of Y-linked
#  sperm-competition allele

wgbar  <-  function(qy, wgYs) {
	((1 - qy)*wgYs[1]) + (qy*wgYs[2])
}
FgY  <-  function(qy, wgYs) {
	((1 - qy)*wgYs[1])/wgbar(qy = qy, wgYs = wgYs)
} 
Fgy  <-  function(qy, wgYs) {
	(qy*wgYs[2])/wgbar(qy = qy, wgYs = wgYs)
}

# Recursions
wobar  <-  function(qy, wgYs, Fij, woijY, woijy) {
	(FgY(qy=qy, wgYs=wgYs)*((Fij[1]*woijY[1]) + (Fij[2]*woijY[2]) + (Fij[3]*woijY[3]))) + 
	(Fgy(qy=qy, wgYs=wgYs)*((Fij[1]*woijy[1]) + (Fij[2]*woijy[2]) + (Fij[3]*woijy[3])))
}
FoY  <-  function(qy, wgYs, Fij, woijY, woijy) {
	(FgY(qy=qy, wgYs=wgYs)*((Fij[1]*woijY[1]) + (Fij[2]*woijY[2]) + (Fij[3]*woijY[3]))) / wobar(qy=qy, wgYs=wgYs, Fij=Fij, woijY=woijY, woijy=woijy)
}
Foy  <-  function(qy, wgYs, Fij, woijY, woijy) {
	(Fgy(qy=qy, wgYs=wgYs)*((Fij[1]*woijy[1]) + (Fij[2]*woijy[2]) + (Fij[3]*woijy[3]))) / wobar(qy=qy, wgYs=wgYs, Fij=Fij, woijY=woijY, woijy=woijy)
}
qyPr  <-  function(qy, wgYs, Fij, woijY, woijy) {
	Foy(qy=qy, wgYs=wgYs, Fij=Fij, woijY=woijY, woijy=woijy)
}
deltaqy  <-  function(qy, wgYs, Fij, woijY, woijy) {
	qyPr(qy=qy, wgYs=wgYs, Fij=Fij, woijY=woijY, woijy=woijy) - qy
}
qyt  <-  function(qy0, sel, t) {
	(qy0*exp(sel* t))/(1 - qy0 + qy0*exp(sel*t))
}


###################################
#  Simulations to check approximations for 
#  invasion of autosomal compensatory allele
#  when mutant y is fixed
FAAPr  <-  function(Fij, woijy){
	((Fij[2] + 2*Fij[1])*(Fij[2]*woijy[2] + 2*Fij[1]*woijy[1])) / (4*(Fij[1] + Fij[2] + Fij[3])*(Fij[3]*woijy[3] + Fij[2]*woijy[2] + Fij[1]*woijy[1]))
}
FAaPr  <-  function(Fij, woijy){
	(Fij[2]*(Fij[2] + Fij[1])*woijy[2] + Fij[3]*Fij[2]*(woijy[3] + woijy[2]) + Fij[2]*Fij[1]*woijy[1] + 2*Fij[3]*Fij[1]*(woijy[3] + woijy[1])) / (2*(Fij[1] + Fij[2] + Fij[3])*(Fij[3]*woijy[3] + Fij[2]*woijy[2] + Fij[1]*woijy[1]))
}
FaaPr  <-  function(Fij, woijy){
	((2*Fij[3] + Fij[2])*(2*Fij[3]*woijy[3] + Fij[2]*woijy[2])) / (4*(Fij[1] + Fij[2] + Fij[3])*(Fij[3]*woijy[3] + Fij[2]*woijy[2] + Fij[1]*woijy[1]))
}

qPr.std  <-  function(qa, woijy) {
	(qa*(1 - qa)*woijy[2]+(qa^2)*woijy[3]) / (((1-qa)^2)*woijy[1] + 2*qa*(1 - qa)*woijy[2]+(qa^2)*woijy[3])
}

qPr.App.std  <-  function(qa, ho,so) {
	(qa*(2 + (-1 + qa)*(1 + ho - qa + 2*ho*qa)*so))/(2 + 2*(-1 + qa)*(1 + (-1 + 2*ho)*qa)*so)
}

avSelRetard  <-  function(q, s) {
	f.ret  <-  ((1-q)^2)*(1/8) + 2*q*(1-q)*(3/8) + (q^2)*(9/12)
	s*f.ret
}

woijyRetard  <-  function(q, h, s) {
	c((1 - ((1-q)^2)*(1/8)*s), (1 - (2*q*(1-q)*(3/8)*h*s)), (q^2)*1)
}



###################################
#  Simulations to check approximations for 
#  invasion of X-linked compensatory allele
#  when mutant y is fixed
FXyPr  <-  function(Fiy, FXij, woijy){
	((FXij[2]*woijy[2])/2 + FXij[1]*woijy[1]) / (FXij[3]*woijy[3] + FXij[2]*woijy[1] + FXij[1]*woijy[1])
}
FxyPr  <-  function(Fiy, FXij, woijy){
	(FXij[3]*woijy[3] + (FXij[2]*woijy[2])/2) / (FXij[3]*woijy[3] + FXij[2]*woijy[1] + FXij[1]*woijy[1])
}

FXXPr  <-  function(Fiy, FXij, woijy){
	(Fiy[1]*(FXij[2]*woijy[2] + 2*FXij[1]*woijy[1])) / (2*(Fiy[2] + Fiy[1])*(FXij[3]*woijy[3] + FXij[2]*woijy[2] + FXij[1]*woijy[1]))
}
FXxPr  <-  function(Fiy, FXij, woijy){
	(2*FXij[3]*Fiy[1]*woijy[3] + FXij[2]*(Fiy[2] + Fiy[1])*woijy[2] + 2*FXij[1]*Fiy[2]*woijy[1]) / (2*(Fiy[2] + Fiy[1])*(FXij[3]*woijy[3] + FXij[2]*woijy[2] + FXij[1]*woijy[1]))
}
FxxPr  <-  function(Fiy, FXij, woijy){
	(Fiy[2]*(2*FXij[3]*woijy[3] + FXij[2]*woijy[2])) / (2*(Fiy[2] + Fiy[1])*(FXij[3]*woijy[3] + FXij[2]*woijy[2] + FXij[1]*woijy[1]))
}

qPr.XLinked  <-  function(qx, woim, woijf) {
	(1/3)*(2*(qx*(1-qx)*((1-qx)*(woijf[2] - woijf[1]) + qx*(woijf[3] - woijf[2]))) + 
		((qx*(1-qx)*(woim[2] - woim[1]))/((1 - qx)*woim[1] + qx*woim[2])))
}



###################################
#  Finite population allele trajectories
#  for invading alleles

MSApproxGivenFix  <-  function(N, q0, s, t) {
	(q0*exp(s*t)) / (1 - q0 + q0*exp(s*t))
}
deltaqyGivenFix  <-  function(N, qy, s) {
	s*qy*(1 - qy)*(1/tanh(((N*s*qy)/2)))
}
deltaqyGivenLoss  <-  function(N, qy, s) {
	s*qy*(1 - qy)*(1/tanh(((N*s*(qy - 1))/2)))
}

###################################
#  Genetic variance functions
oneLocVar  <-  function(q, s) {
	q*(1-q)*s^2
}

###################################
#  Threshold frequencies of y at 
#  which compensatory mutations are
#  favoured
qTildeA  <-  function(sm, ho, delta, hc, sc) {
	(hc*sc) / (hc*sc + (sm - delta)*(1 - ho)*(1 + sm))
}
qTildeM  <-  function(sm, delta, sc) {
	sc / (sc*(1 + sm)*(sm - delta))
}
b  <-  function(sm,delta,sc) {
	sm*(2 + sm*(2 - sm*(1 - sm))) - delta*(2 + sm + (sm^2) + (2*(sm^3))) + (delta^2)*(1 + sm^2) + sc*(2 + sm^2 - delta*(1 + sm))
}
d  <-  function(sm,delta,sc) {
	(sc + sm + (3*sm^2) - 3*delta*(1 + sm))*(sm^2 - delta*(1 + sm))
}
qTildeX  <-  function(sm, ho, delta, hc, sc) {
	(b(sm=sm, delta=delta, sc=sc) - sqrt(b(sm=sm, delta=delta, sc=sc)^2 - 8*sc*d(sm=sm, delta=delta, sc=sc))) / (2*d(sm=sm, delta=delta, sc=sc))
}


