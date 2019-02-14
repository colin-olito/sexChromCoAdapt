########################################################
#  Antagonistic coevolution between the sex-chromosomes
#
#  Functions to run stochastic simulation of 
#  invasion of Y-linked male-benefit alleles and  
#  MITOCHONDRIAL, AUTOSOMAL, and X-LINKED compensatory
#  mutations. Functions create output data frames for
#  plotting
#
#  Author: Colin Olito
#
#  NOTES: 
#          

rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')
source('R/functions-Mito-Simulations2.R')
source('R/functions-Auto-Simulations2.R')
source('R/functions-Xlinked-Simulations2.R')

######################################
# Equal dominance (ho = hc = h)
reps     <-  500000
N        <-  1000
sm       <-  0.1
delta    <-  0.095
sc       <-  0.005
ho.vals  <-  seq(from = 0, to = 1, by = 0.1)
hc.vals  <-  seq(from = 0, to = 1, by = 0.1)

makeDataYAutoInv(reps = reps, N = N, sm = sm, 
				 ho = ho.vals, delta = delta, 
				 hc = hc.vals, sc = sc)

makeDataYXLinkedInv(reps = reps, N = N, sm = sm, 
					ho.vals = ho.vals, delta = delta, 
					hc.vals = hc.vals, sc = sc)

delta.vals    <-  c(0.005, 0.05, 0.095)
for(i in 1:length(delta.vals)) {
	makeDataYAutoInv(reps = reps, N = N, sm = sm, 
					 ho = ho.vals, delta = delta.vals[i], 
					 hc = hc.vals, sc = sc)

	makeDataYXLinkedInv(reps = reps, N = N, sm = sm, 
						ho.vals = ho.vals, delta = delta.vals[i], 
						hc.vals = hc.vals, sc = sc)	
}

# Re-run sims for delta = 0.095 with higher reps to clean up plot
reps     <-  5000000
N        <-  1000
sm       <-  0.1
delta    <-  0.095
sc       <-  0.005
ho.vals  <-  seq(from = 0, to = 1, by = 0.1)
hc.vals  <-  seq(from = 0, to = 1, by = 0.1)

makeDataYAutoInv(reps = reps, N = N, sm = sm, 
				 ho = ho.vals, delta = delta, 
				 hc = hc.vals, sc = sc)

makeDataYXLinkedInv(reps = reps, N = N, sm = sm, 
					ho.vals = ho.vals, delta = delta, 
					hc.vals = hc.vals, sc = sc)

### Simulations for time to invasion/fixation
###  exploring dominance gradient
reps     <-  50000
N        <-  1000
sm       <-  0.1
delta    <-  0.01
sc       <-  0.01
ho.vals  <-  seq(from = 0, to = 1, by = 0.1)
hc.vals  <-  seq(from = 0, to = 1, by = 0.1)
uy       <-  1e-3
ux       <-  1.5e-3
ua       <-  2e-3
makeDataYXTimeFix(reps = reps, N = N, sm = sm, 
				  ho.vals = ho.vals, delta = delta, 
				  hc.vals = hc.vals, sc = sc,
				  uy = uy, qy.init = "singleCopy",
				  ux = ux, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25))

makeDataYATimeFix(reps = reps, N = N, sm = sm, 
				  ho.vals = ho.vals, delta = delta, 
				  hc.vals = hc.vals, sc = sc,
				  uy = uy, qy.init = "singleCopy",
				  ua = ua, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25))



### Simulations for time to invasion/fixation
###  exploring Mutation rate gradient

reps     <-  500
N        <-  1000
sm       <-  0.1
delta    <-  0.05
sc       <-  0.005
ho       <-  0.5
hc       <-  0.5
uy       <-  1e-3
u.vals   <-  seq(1e-7,uy,len=11)
#u.vals   <-  log10space(b=1, d1=-7, d2=-3, n=10)
#seq(u.vals[7], u.vals[10], len=4)
#u.vals  <-  c(u.vals[1:6], seq(u.vals[7], u.vals[10], len=4))

ux.vals  <-  u.vals
ua.vals  <-  u.vals


h.vals  <-  c(0.1, 0.5, 0.9)
sm.vals  <-  c(0.1)
delta.vals  <-  c(0.05, 0.05,
				  0.095)
sc.vals  <-  c(0.05, 0.005,
			   0.005)

for(i in 1:length(h.vals)) {
	for(j in 1:length(sm.vals)) {
		for(k in 1:length(delta.vals)) {
			makeDataYXTimeFixMutGrad(reps = reps, N = N, sm = sm.vals[j], 
				  ho = h.vals[i], delta = delta.vals[k], 
				  hc = h.vals[i], sc = sc.vals[k],
				  uy = uy, qy.init = "singleCopy",
				  ux = ux.vals, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25))

			makeDataYATimeFixMutGrad(reps = reps, N = N, sm = sm.vals[j], 
				  ho = h.vals[i], delta = delta.vals[k], 
				  hc = h.vals[i], sc = sc.vals[k],
				  uy = uy, qy.init = "singleCopy",
				  ua = ua.vals, qx.init = 0, Ff.init = c(0.25, 0.5, 0.25))


		}
	}
}
