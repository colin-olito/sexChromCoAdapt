#####################################################
#  2-locus SA with partial selfing
#
#  Functions to generate Figures for: 
#    Article title goes here...
#  
#  Author: Colin Olito
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
#          

rm(list=ls())

library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)

#source('paths.R')
source('R/functions-analyses.R')
source('R/functions-figures.R')


###############
# PAPER FIGURES
###############

toPdf(theoryFig(), figPath(name='theoryFig.pdf'), width=10.5, height=7)
embed_fonts(figPath(name='theoryFig.pdf'))

source('R/functions-figures.R')
toPdf(theoryFig2(), figPath(name='theoryFig2.pdf'), width=15, height=5)
embed_fonts(figPath(name='theoryFig2.pdf'))

# toPdf(Fig.1wk(), figPath(name='Fig1wk.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig1wk.pdf'))

# toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig2.pdf'))

# toPdf(Fig.2wk(), figPath(name='Fig2wk.pdf'), width=7, height=7)
# embed_fonts(figPath(name='Fig2wk.pdf'))

# toPdf(recSimFig_add(), figPath(name='recSimFig_add.pdf'), width=7, height=7)
# embed_fonts(figPath(name='recSimFig_add.pdf'))

# toPdf(recSimFig_domRev(), figPath(name='recSimFig_domRev.pdf'), width=7, height=7)
# embed_fonts(figPath(name='recSimFig_domRev.pdf'))


##############################
# SUPPLEMENTARY FIGURES
##############################

###############
# OTHER FIGURES
###############

toPdf(relInvasionProbFigLegend(), figPath(name='relInvasionProbTestFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='relInvasionProbTestFig.pdf'))

toPdf(relCycleTimeFig2(), figPath(name='relCycleTimeTestFig.pdf'), width=7, height=7)
embed_fonts(figPath(name='relCycleTimeTestFig.pdf'))
