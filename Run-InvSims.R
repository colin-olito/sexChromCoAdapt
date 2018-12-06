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
source('R/functions-Auto-Simulations2.R')

######################################
# Equal dominance (ho = hc = h)
reps     <-  30000
N        <-  1000
sm       <-  0.1
delta    <-  0.05
sc       <-  0.05
ho.vals  <-  seq(from = 0, to = 1, by = 0.05)
hc.vals  <-  seq(from = 0, to = 1, by = 0.05)

makeDataYAutoInv(reps = reps, N = N, sm = sm, 
				 ho = ho.vals, delta = delta, 
				 hc = hc.vals, sc = sc)

rm(list=ls())
source('R/functions-analyses.R')
source('R/functions-Xlinked-Simulations2.R')

reps     <-  30000
N        <-  1000
sm       <-  0.1
delta    <-  0.05
sc       <-  0.05
ho.vals  <-  seq(from = 0, to = 1, by = 0.05)
hc.vals  <-  seq(from = 0, to = 1, by = 0.05)

makeDataYAutoInv(reps = reps, N = N, sm = sm, 
				 ho = ho.vals, delta = delta, 
				 hc = hc.vals, sc = sc)
