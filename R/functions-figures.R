###############
# DEPENDENCIES
###############

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


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


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}



##############################################################
##############################################################
##  Final figures for paper

theoryFig  <-  function() {

    # Set plot layout
    layout.mat <- matrix(c(1,1,2,
                           1,1,3), nrow=2, ncol=3, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    # Calculate qTilda values for plotting
    sms  <-  seq(0,0.1, by=0.001)

    qTildeA_y1.add  <-  qTildeA(sm = sms, ho=1/2, so=0.0025, hc=1/2, sc=0.005)
    qTildeX_y1.add  <-  qTildeX(sm = sms, ho=1/2, so=0.0025, hc=1/2, sc=0.005)
    qTildeA_y2.add  <-  qTildeA(sm = sms, ho=1/2, so=0.005, hc=1/2, sc=0.005)
    qTildeX_y2.add  <-  qTildeX(sm = sms, ho=1/2, so=0.005, hc=1/2, sc=0.005)
    qTildeA_y3.add  <-  qTildeA(sm = sms, ho=1/2, so=0.025, hc=1/2, sc=0.005)
    qTildeX_y3.add  <-  qTildeX(sm = sms, ho=1/2, so=0.025, hc=1/2, sc=0.005)


# Panel A: qTilde_y examples
    par(omi=rep(0.5, 4), mar = c(3,3,4,4), bty='o', xaxt='s', yaxt='s')
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        lines(qTildeA_y3.add ~ sms, lwd=3, lty=1, col='grey70')
        lines(qTildeX_y3.add ~ sms, lwd=3, lty=2, col='grey70')
        lines(qTildeA_y2.add ~ sms, lwd=3, lty=1, col='grey50')
        lines(qTildeX_y2.add ~ sms, lwd=3, lty=2, col='grey50')
        lines(qTildeA_y1.add ~ sms, lwd=3, lty=1, col='#252525')
        lines(qTildeX_y1.add ~ sms, lwd=3, lty=2, col='#252525')
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.125, 0.5, expression(paste(tilde(italic(q))[italic(y)])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
#        proportionalLabel(-0.15, 0.5, expression(paste(tilde(italic(q))[italic(y)]^italic(A)/tilde(italic(q))[italic(y)]^italic(X))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.1, expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.032, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(s[o])~"="~"0.005")),
                          expression(paste(italic(s[o])~"="~"0.0025")),
                          expression(paste(italic(s[o])~"="~"0.0015"))),
              lty     =  1,
              lwd     =  3,
              col     =  c('#252525', 'grey50', 'grey70'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )
        legend(
              x       =  usr[2]*0.75,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              lty     =  c(1,2),
              lwd     =  3,
              seg.len =  3,
              col     =  '#252525',
              cex     =  1.2,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )

# Panel B: Relative Probability of invasion
    # Import Data
    aDat = './output/data/simData/dataYAutoInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    xDat = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    aDat1 = './output/data/simData/dataYAutoInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    xDat1 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    aDat2 = './output/data/simData/dataYAutoInv_sm0.1_delta0.095_sc0.005_N1000_reps5e+06.csv'
    xDat2 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.095_sc0.005_N1000_reps5e+06.csv'
    autoDat     <-  read.csv(aDat, header=TRUE)
    XLinkedDat  <-  read.csv(xDat, header=TRUE)
    autoDat1     <-  read.csv(aDat1, header=TRUE)
    XLinkedDat1  <-  read.csv(xDat1, header=TRUE)
    autoDat2     <-  read.csv(aDat2, header=TRUE)
    XLinkedDat2  <-  read.csv(xDat2, header=TRUE)

    # Calculate Relative Invasion Probabilities
    Py   <-  autoDat$pInvade_y / XLinkedDat$pInvade_y
    Pxa  <-  XLinkedDat$pInvade_x / autoDat$pInvade_a
    Py1   <-  autoDat1$pInvade_y / XLinkedDat1$pInvade_y
    Pxa1  <-  XLinkedDat1$pInvade_x / autoDat1$pInvade_a
    Py2  <-  autoDat2$pInvade_y / XLinkedDat2$pInvade_y
    Pxa2 <-  XLinkedDat2$pInvade_x / autoDat2$pInvade_a

        # Make the plot
#        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.9,1.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
#        points(Py ~ autoDat$ho, pch=21, ylim=c(0,2))
#        points(Py2 ~ autoDat2$ho, pch=21, ylim=c(0,2))
        abline(h=1,lty=2, lwd=2)
        points(Pxa  ~ autoDat$ho, pch=21, bg='grey60')
        points(Pxa1  ~ autoDat$ho, pch=21, bg='grey80')
        points(Pxa2 ~ autoDat$ho, pch=21, bg=transparentColor('#252525', opacity=0.85))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.25, 0.5, expression(paste(Pi[X]/Pi[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(h[i]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.99,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(italic(delta)~"="~"0.005")),
                          expression(paste(italic(delta)~"="~"0.05")),
                          expression(paste(italic(delta)~"="~"0.095"))),
              pch     =  21,
              pt.bg   =  c('grey60', 'grey80', transparentColor('#252525', opacity=0.85)),
              cex     =  0.75,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )



# Panel C: Time to complete coevolutionary cycle
    YA1         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"
    YX1         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"
    YA2         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.095_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"
    YX2         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.095_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"
    aData       <-  read.csv(YA1, header=TRUE)
    xData       <-  read.csv(YX1, header=TRUE)
    aData2      <-  read.csv(YA2, header=TRUE)
    xData2      <-  read.csv(YX2, header=TRUE)
    tCycle      <-  xData$tCycle / aData$tCycle
    tCycle2     <-  xData2$tCycle / aData2$tCycle
    relMutRate  <-  aData$ua / aData$uy

        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.75,1.2), ylab='', xlab='', cex.lab=1.2)
#        plot(NA, axes=FALSE, type='n', main='',xlim = c(min(log10(relMutRate)),max(log10(relMutRate))), ylim = c(0.5,1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        abline(h=1,lty=2, lwd=2)
        points(tCycle ~ (relMutRate), pch=21, bg="#252525", col='black', ylim=c(0,2))
        points(tCycle2 ~ (relMutRate), pch=21, bg="grey80", col='black', ylim=c(0,2))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.25, 0.5, expression(paste(T["cycle"]^X/T["cycle"]^A)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(mu[i])/italic(mu[y]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*1,
              y       =  usr[4]*0.6,
              legend  =  c(
                          expression(paste(italic(s[m])~"="~"0.1;"~italic(delta)~"="~"0.050;"~italic(s[c])~"="~"0.005;"~italic(N)~"="~"1000;")),
                          expression(paste(italic(s[m])~"="~"0.1;"~italic(delta)~"="~"0.090;"~italic(s[c])~"="~"0.000;"~italic(N)~"="~"1000;"))),
              pch     =  21,
              col     =  'black',
              pt.bg   =  c('#252525', 'grey80'),
              cex     =  0.7,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )

}




theoryFig2  <-  function() {

    # Set plot layout
#    layout.mat <- matrix(c(1,1,2,
#                           1,1,3), nrow=2, ncol=3, byrow=TRUE)
    layout.mat <- matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    # Calculate qTilda values for plotting
    sms  <-  seq(0,0.1, by=0.001)

    qTildeA_y1.add  <-  qTildeA(sm = sms, ho=1/2, so=0.0025, hc=1/2, sc=0.005)
    qTildeX_y1.add  <-  qTildeX(sm = sms, ho=1/2, so=0.0025, hc=1/2, sc=0.005)
    qTildeA_y2.add  <-  qTildeA(sm = sms, ho=1/2, so=0.005, hc=1/2, sc=0.005)
    qTildeX_y2.add  <-  qTildeX(sm = sms, ho=1/2, so=0.005, hc=1/2, sc=0.005)
    qTildeA_y3.add  <-  qTildeA(sm = sms, ho=1/2, so=0.025, hc=1/2, sc=0.005)
    qTildeX_y3.add  <-  qTildeX(sm = sms, ho=1/2, so=0.025, hc=1/2, sc=0.005)


# Panel A: qTilde_y examples
    par(omi=c(0.5,0.5,0.5,0.5), mar = c(3.5,4,3.5,4), bty='o', xaxt='s', yaxt='s')
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        lines(qTildeA_y3.add ~ sms, lwd=3, lty=1, col='grey70')
        lines(qTildeX_y3.add ~ sms, lwd=3, lty=2, col='grey70')
        lines(qTildeA_y2.add ~ sms, lwd=3, lty=1, col='grey50')
        lines(qTildeX_y2.add ~ sms, lwd=3, lty=2, col='grey50')
        lines(qTildeA_y1.add ~ sms, lwd=3, lty=1, col='#252525')
        lines(qTildeX_y1.add ~ sms, lwd=3, lty=2, col='#252525')
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.2, 0.5, expression(paste(tilde(italic(q))[italic(y)])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
#        proportionalLabel(-0.15, 0.5, expression(paste(tilde(italic(q))[italic(y)]^italic(A)/tilde(italic(q))[italic(y)]^italic(X))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.2, expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(delta)~"="~"0.095")),
                          expression(paste(italic(delta)~"="~"0.0975")),
                          expression(paste(italic(delta)~"="~"0.0985"))),
#              legend  =  c(
#                          expression(paste(italic(s[o])~"="~"0.005")),
#                          expression(paste(italic(s[o])~"="~"0.0025")),
#                          expression(paste(italic(s[o])~"="~"0.0015"))),
              lty     =  1,
              lwd     =  3,
              col     =  c('#252525', 'grey50', 'grey70'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )
        legend(
              x       =  usr[2]*0.6,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              lty     =  c(1,2),
              lwd     =  3,
              seg.len =  3,
              col     =  '#252525',
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )

# Panel B: Relative Probability of invasion
    # Import Data
    aDat = './output/data/simData/dataYAutoInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    xDat = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    aDat1 = './output/data/simData/dataYAutoInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    xDat1 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    aDat2 = './output/data/simData/dataYAutoInv_sm0.1_delta0.095_sc0.005_N1000_reps7500000.csv'
    xDat2 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.095_sc0.005_N1000_reps7500000.csv'
    autoDat     <-  read.csv(aDat, header=TRUE)
    XLinkedDat  <-  read.csv(xDat, header=TRUE)
    autoDat1     <-  read.csv(aDat1, header=TRUE)
    XLinkedDat1  <-  read.csv(xDat1, header=TRUE)
    autoDat2     <-  read.csv(aDat2, header=TRUE)
    XLinkedDat2  <-  read.csv(xDat2, header=TRUE)

    # Calculate Relative Invasion Probabilities
    Py   <-  autoDat$pInvade_y / XLinkedDat$pInvade_y
    Pxa  <-  XLinkedDat$pInvade_x / autoDat$pInvade_a
    Py1   <-  autoDat1$pInvade_y / XLinkedDat1$pInvade_y
    Pxa1  <-  XLinkedDat1$pInvade_x / autoDat1$pInvade_a
    Py2  <-  autoDat2$pInvade_y / XLinkedDat2$pInvade_y
    Pxa2 <-  XLinkedDat2$pInvade_x / autoDat2$pInvade_a

        # Make the plot
#        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.9,1.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
#        points(Py ~ autoDat$ho, pch=21, ylim=c(0,2))
#        points(Py2 ~ autoDat2$ho, pch=21, ylim=c(0,2))
        abline(h=1,lty=2, lwd=2)
        points(Pxa  ~ autoDat$ho, pch=21, bg='grey80', cex=1.5)
        points(Pxa1  ~ autoDat$ho, pch=21, bg='grey60', cex=1.5)
        points(Pxa2 ~ autoDat$ho, pch=21, bg=transparentColor('#252525', opacity=0.85), cex=1.5)
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.2, 0.5, expression(paste(Pi[X]/Pi[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.2, expression(paste(italic(h))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.98,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(delta)~"="~"0.005")),
                          expression(paste(italic(delta)~"="~"0.05")),
                          expression(paste(italic(delta)~"="~"0.095"))),
              pch     =  21,
              pt.bg   =  c('grey80', 'grey60', transparentColor('#252525', opacity=0.85)),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )



# Panel C: Time to complete coevolutionary cycle
    YA1          <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.005_hc0.1_N1000_reps500.csv"
    YX1          <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.05_hc0.1_N1000_reps500.csv"
    YA2         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.05_hc0.1_N1000_reps500.csv"
    YX2         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.005_hc0.1_N1000_reps500.csv"
    YA3         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"
    YX3         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.05_hc0.5_N1000_reps500.csv"
    YA4         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.05_hc0.5_N1000_reps500.csv"
    YX4         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"

    YA5         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.005_hc0.9_N1000_reps500.csv"
    YX5         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.05_hc0.9_N1000_reps500.csv"
    YA6         <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.05_hc0.9_N1000_reps500.csv"
    YX6         <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.005_hc0.9_N1000_reps500.csv"

    aData1      <-  read.csv(YA1, header=TRUE)
    xData1      <-  read.csv(YX1, header=TRUE)
    aData2      <-  read.csv(YA2, header=TRUE)
    xData2      <-  read.csv(YX2, header=TRUE)
    aData3      <-  read.csv(YA3, header=TRUE)
    xData3      <-  read.csv(YX3, header=TRUE)
    aData4      <-  read.csv(YA4, header=TRUE)
    xData4      <-  read.csv(YX4, header=TRUE)
    aData5      <-  read.csv(YA5, header=TRUE)
    xData5      <-  read.csv(YX5, header=TRUE)
    aData6      <-  read.csv(YA6, header=TRUE)
    xData6      <-  read.csv(YX6, header=TRUE)
    tCycle1     <-  xData1$tCycle / aData1$tCycle
    tCycle2     <-  xData2$tCycle / aData2$tCycle
    tCycle3     <-  xData3$tCycle / aData3$tCycle
    tCycle4     <-  xData4$tCycle / aData4$tCycle
    tCycle5     <-  xData5$tCycle / aData5$tCycle
    tCycle6     <-  xData6$tCycle / aData6$tCycle
    relMutRate  <-  aData1$ua / aData1$uy

        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.7,1.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        abline(h=1,lty=2, lwd=2)
        points(tCycle1 ~ (relMutRate), pch=22, bg="#252525", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle2 ~ (relMutRate), pch=22, bg="grey80", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle3 ~ (relMutRate), pch=21, bg="#252525", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle4 ~ (relMutRate), pch=21, bg="grey80", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle5 ~ (relMutRate), pch=24, bg="#252525", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle6 ~ (relMutRate), pch=24, bg="grey80", col='black', cex=1.5, type='b', ylim=c(0,2))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.2, 0.5, expression(paste(T[x]/T[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.2, expression(paste(italic(mu[i])/italic(mu[y]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.65,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(s[c])~"="~"0.05")),
                          expression(paste(italic(s[c])~"="~"0.005"))),
              pch     =  21,
              col     =  'black',
              pt.bg   =  c('#252525', 'grey80'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )

        legend(
              x       =  usr[2]*0.98,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(h[o])~"="~italic(h[c])~"="~"0.9")),
                          expression(paste(italic(h[o])~"="~italic(h[c])~"="~"0.5")),
                          expression(paste(italic(h[o])~"="~italic(h[c])~"="~"0.1"))),
              pch     =  c(24,21,22),
              col     =  'black',
              pt.bg   =  NA,
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )


}





theoryFig2vert  <-  function() {

    # Set plot layout
#    layout.mat <- matrix(c(1,1,2,
#                           1,1,3), nrow=2, ncol=3, byrow=TRUE)
    layout.mat <- matrix(c(1,2,3), nrow=3, ncol=1, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    # Calculate qTilda values for plotting
    sms  <-  seq(0,0.1, by=0.001)

    qTildeA_y1.add  <-  qTildeA(sm = sms, ho=1/2, so=0.0025, hc=1/2, sc=0.005)
    qTildeX_y1.add  <-  qTildeX(sm = sms, ho=1/2, so=0.0025, hc=1/2, sc=0.005)
    qTildeA_y2.add  <-  qTildeA(sm = sms, ho=1/2, so=0.005, hc=1/2, sc=0.005)
    qTildeX_y2.add  <-  qTildeX(sm = sms, ho=1/2, so=0.005, hc=1/2, sc=0.005)
    qTildeA_y3.add  <-  qTildeA(sm = sms, ho=1/2, so=0.025, hc=1/2, sc=0.005)
    qTildeX_y3.add  <-  qTildeX(sm = sms, ho=1/2, so=0.025, hc=1/2, sc=0.005)


# Panel A: qTilde_y examples
    par(omi=c(0.5,0.5,0.5,0.5), mar = c(4,4,4,3), bty='o', xaxt='s', yaxt='s')
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        lines(qTildeA_y3.add ~ sms, lwd=3, lty=1, col='grey70')
        lines(qTildeX_y3.add ~ sms, lwd=3, lty=2, col='grey70')
        lines(qTildeA_y2.add ~ sms, lwd=3, lty=1, col='grey50')
        lines(qTildeX_y2.add ~ sms, lwd=3, lty=2, col='grey50')
        lines(qTildeA_y1.add ~ sms, lwd=3, lty=1, col='#252525')
        lines(qTildeX_y1.add ~ sms, lwd=3, lty=2, col='#252525')
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.2, 0.5, expression(paste(tilde(italic(q))[italic(y)])), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
#        proportionalLabel(-0.15, 0.5, expression(paste(tilde(italic(q))[italic(y)]^italic(A)/tilde(italic(q))[italic(y)]^italic(X))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.2, expression(paste(italic(s[m]))), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(A))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.95,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(delta)~"="~"0.095")),
                          expression(paste(italic(delta)~"="~"0.0975")),
                          expression(paste(italic(delta)~"="~"0.0985"))),
#              legend  =  c(
#                          expression(paste(italic(s[o])~"="~"0.005")),
#                          expression(paste(italic(s[o])~"="~"0.0025")),
#                          expression(paste(italic(s[o])~"="~"0.0015"))),
              lty     =  1,
              lwd     =  3,
              col     =  c('#252525', 'grey50', 'grey70'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )
        legend(
              x       =  usr[2]*0.6,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste("Autosomal")),
                          expression(paste("X-linked"))),
              lty     =  c(1,2),
              lwd     =  3,
              seg.len =  3,
              col     =  '#252525',
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )

# Panel B: Relative Probability of invasion
    # Import Data
    aDat = './output/data/simData/dataYAutoInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    xDat = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    aDat1 = './output/data/simData/dataYAutoInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    xDat1 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    aDat2 = './output/data/simData/dataYAutoInv_sm0.1_delta0.095_sc0.005_N1000_reps7500000.csv'
    xDat2 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.095_sc0.005_N1000_reps7500000.csv'
    autoDat     <-  read.csv(aDat, header=TRUE)
    XLinkedDat  <-  read.csv(xDat, header=TRUE)
    autoDat1     <-  read.csv(aDat1, header=TRUE)
    XLinkedDat1  <-  read.csv(xDat1, header=TRUE)
    autoDat2     <-  read.csv(aDat2, header=TRUE)
    XLinkedDat2  <-  read.csv(xDat2, header=TRUE)

    # Calculate Relative Invasion Probabilities
    Py   <-  autoDat$pInvade_y / XLinkedDat$pInvade_y
    Pxa  <-  XLinkedDat$pInvade_x / autoDat$pInvade_a
    Py1   <-  autoDat1$pInvade_y / XLinkedDat1$pInvade_y
    Pxa1  <-  XLinkedDat1$pInvade_x / autoDat1$pInvade_a
    Py2  <-  autoDat2$pInvade_y / XLinkedDat2$pInvade_y
    Pxa2 <-  XLinkedDat2$pInvade_x / autoDat2$pInvade_a

        # Make the plot
#        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.9,1.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
#        points(Py ~ autoDat$ho, pch=21, ylim=c(0,2))
#        points(Py2 ~ autoDat2$ho, pch=21, ylim=c(0,2))
        abline(h=1,lty=2, lwd=2)
        points(Pxa  ~ autoDat$ho, pch=21, bg='grey80', cex=1.5)
        points(Pxa1  ~ autoDat$ho, pch=21, bg='grey60', cex=1.5)
        points(Pxa2 ~ autoDat$ho, pch=21, bg=transparentColor('#252525', opacity=0.85), cex=1.5)
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.2, 0.5, expression(paste(Pi[X]/Pi[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.2, expression(paste(italic(h))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(B))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.98,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(delta)~"="~"0.005")),
                          expression(paste(italic(delta)~"="~"0.05")),
                          expression(paste(italic(delta)~"="~"0.095"))),
              pch     =  21,
              pt.bg   =  c('grey80', 'grey60', transparentColor('#252525', opacity=0.85)),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )



# Panel C: Time to complete coevolutionary cycle
    YA1      <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.005_hc0.1_N1000_reps500.csv"
    YX1      <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.05_hc0.1_N1000_reps500.csv"
    YA2      <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.05_hc0.1_N1000_reps500.csv"
    YX2      <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.1_sc0.005_hc0.1_N1000_reps500.csv"
    YA3      <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"
    YX3      <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.05_hc0.5_N1000_reps500.csv"
    YA4      <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.05_hc0.5_N1000_reps500.csv"
    YX4      <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps500.csv"

    YA5      <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.005_hc0.9_N1000_reps500.csv"
    YX5      <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.05_hc0.9_N1000_reps500.csv"
    YA6      <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.05_hc0.9_N1000_reps500.csv"
    YX6      <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.9_sc0.005_hc0.9_N1000_reps500.csv"

    aData1   <-  read.csv(YA1, header=TRUE)
    xData1   <-  read.csv(YX1, header=TRUE)
    aData2   <-  read.csv(YA2, header=TRUE)
    xData2   <-  read.csv(YX2, header=TRUE)
    aData3   <-  read.csv(YA3, header=TRUE)
    xData3   <-  read.csv(YX3, header=TRUE)
    aData4   <-  read.csv(YA4, header=TRUE)
    xData4   <-  read.csv(YX4, header=TRUE)
    aData5   <-  read.csv(YA5, header=TRUE)
    xData5   <-  read.csv(YX5, header=TRUE)
    aData6   <-  read.csv(YA6, header=TRUE)
    xData6   <-  read.csv(YX6, header=TRUE)
    tCycle1  <-  xData1$tCycle / aData1$tCycle
    tCycle2  <-  xData2$tCycle / aData2$tCycle
    tCycle3  <-  xData3$tCycle / aData3$tCycle
    tCycle4  <-  xData4$tCycle / aData4$tCycle
    tCycle5  <-  xData5$tCycle / aData5$tCycle
    tCycle6  <-  xData6$tCycle / aData6$tCycle

    relMutRate  <-  aData1$ua / aData1$uy

        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.7,1.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        abline(h=1,lty=2, lwd=2)
        points(tCycle1 ~ (relMutRate), pch=22, bg="#252525", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle2 ~ (relMutRate), pch=22, bg="grey80", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle3 ~ (relMutRate), pch=21, bg="#252525", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle4 ~ (relMutRate), pch=21, bg="grey80", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle5 ~ (relMutRate), pch=24, bg="#252525", col='black', cex=1.5, type='b', ylim=c(0,2))
        points(tCycle6 ~ (relMutRate), pch=24, bg="grey80", col='black', cex=1.5, type='b', ylim=c(0,2))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.2, 0.5, expression(paste(T[X]/T[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.2, expression(paste(italic(mu[i])/italic(mu[y]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.075, expression(paste(bold(C))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.65,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(s[c])~"="~"0.05")),
                          expression(paste(italic(s[c])~"="~"0.005"))),
              pch     =  21,
              col     =  'black',
              pt.bg   =  c('#252525', 'grey80'),
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )

        legend(
              x       =  usr[2]*0.98,
              y       =  usr[4]*0.98,
              legend  =  c(
                          expression(paste(italic(h[o])~"="~italic(h[c])~"="~"0.9")),
                          expression(paste(italic(h[o])~"="~italic(h[c])~"="~"0.5")),
                          expression(paste(italic(h[o])~"="~italic(h[c])~"="~"0.1"))),
              pch     =  c(24,21,22),
              col     =  'black',
              pt.bg   =  NA,
              cex     =  1,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )


}

##############################################################
##############################################################
##  Functions for exploratory figures

relInvasionProbFig  <-  function(aDat, xDat) {

    # Import Data
    autoDat     <-  read.csv(aDat, header=TRUE)
    XLinkedDat  <-  read.csv(xDat, header=TRUE)

    # Calculate Relative Invasion Probabilities
    Py  <-  autoDat$pInvade_y / XLinkedDat$pInvade_y
    Pxa  <-  XLinkedDat$pInvade_x / autoDat$pInvade_a

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.75,1.5), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        points(Py ~ autoDat$ho, pch=21, ylim=c(0,2))
        abline(h=1,lty=2, lwd=2)
        points(Pxa ~ autoDat$ho, pch=21, bg='grey80')
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
#        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.15, 0.5, expression(paste(Pi[X]/Pi[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.15, expression(paste(italic(h))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ', 0)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

}

relInvasionProbFigLegend  <-  function() {

    # Import Data
    aDat = './output/data/simData/dataYAutoInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    xDat = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.005_sc0.005_N1000_reps5e+05.csv'
    aDat1 = './output/data/simData/dataYAutoInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    xDat1 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.05_sc0.005_N1000_reps5e+05.csv'
    aDat2 = './output/data/simData/dataYAutoInv_sm0.1_delta0.095_sc0.005_N1000_reps1e+06.csv'
    xDat2 = './output/data/simData/dataYXLinkedInv_sm0.1_delta0.095_sc0.005_N1000_reps1e+06.csv'
    autoDat     <-  read.csv(aDat, header=TRUE)
    XLinkedDat  <-  read.csv(xDat, header=TRUE)
    autoDat1     <-  read.csv(aDat1, header=TRUE)
    XLinkedDat1  <-  read.csv(xDat1, header=TRUE)
    autoDat2     <-  read.csv(aDat2, header=TRUE)
    XLinkedDat2  <-  read.csv(xDat2, header=TRUE)

    # Calculate Relative Invasion Probabilities
    Py   <-  autoDat$pInvade_y / XLinkedDat$pInvade_y
    Pxa  <-  XLinkedDat$pInvade_x / autoDat$pInvade_a
    Py1   <-  autoDat1$pInvade_y / XLinkedDat1$pInvade_y
    Pxa1  <-  XLinkedDat1$pInvade_x / autoDat1$pInvade_a
    Py2  <-  autoDat2$pInvade_y / XLinkedDat2$pInvade_y
    Pxa2 <-  XLinkedDat2$pInvade_x / autoDat2$pInvade_a

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.9,1.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
#        points(Py ~ autoDat$ho, pch=21, ylim=c(0,2))
#        points(Py2 ~ autoDat2$ho, pch=21, ylim=c(0,2))
        abline(h=1,lty=2, lwd=2)
        points(Pxa  ~ autoDat$ho, pch=21, bg='grey60')
        points(Pxa1  ~ autoDat$ho, pch=21, bg='grey80')
        points(Pxa2 ~ autoDat$ho, pch=21, bg=transparentColor('#252525', opacity=0.85))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.15, 0.5, expression(paste(Pi[X]/Pi[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.15, expression(paste(italic(h[i]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ', 0)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*0.99,
              y       =  usr[4],
              legend  =  c(
                          expression(paste(italic(s[m])~"="~"0.1;"~italic(delta)~"="~"0.050;"~italic(s[c])~"="~"0.010;"~italic(N)~"="~"10,000;")),
                          expression(paste(italic(s[m])~"="~"0.1;"~italic(delta)~"="~"0.095;"~italic(s[c])~"="~"0.005;"~italic(N)~"="~"10,000;"))),
              pch     =  21,
              pt.bg   =  c('grey80',transparentColor('#252525', opacity=0.85)),
              cex     =  0.75,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )
}


relCycleTimeFig  <-  function(aDat = './output/data/simData/dataYATimeFix_sm0.1_delta0.05_sc0.05_N1000_reps1000.csv',
                              xDat = './output/data/simData/dataYXTimeFix_sm0.1_delta0.05_sc0.05_N1000_reps1000.csv') {

    # Import Data
    aData  <-  read.csv(aDat, header=TRUE)
    xData  <-  read.csv(xDat, header=TRUE)

    # Calculate additional variables for Plotting
    DeltaInv     <-  aData$deltaInv    / xData$deltaInv
    DeltaFix     <-  aData$deltaFix    / xData$deltaFix
    DeltaInvFix  <-  aData$deltaInvFix / xData$deltaInvFix
    tCycle       <-  aData$tCycle      / xData$tCycle

    # relative size of X vs. Autosomal genome in 
    # D. melanogaster (see Adams et al. 2000)
    MbA  <-  23 + 5.4 + 11 + 21.4 + 24.4 + 8.2 + 8.2 + 28 + 3.1 + 1.2
    MbX  <-  20+21.8
    MbY  <-  40.9

        # Make the plot
#        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.5,1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        points(Py ~ aData$ho, pch=21, bg='red', lwd=2, col=2, ylim=c(0,2))
        abline(h=1,lty=2, lwd=2)
        points(Pax ~ aData$ho, pch=21, bg='grey80', lwd=2, col=1, ylim=c(1/2,1))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
#        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(Pi[A]/Pi[X])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(h))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ', 0)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

}


relCycleTimeFig2  <-  function() {

    # Import Data
YA1  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.05_hc0.5_N1000_reps300.csv"
YX1  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.05_hc0.5_N1000_reps300.csv"
YA2  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps300.csv"
YX2  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0.005_hc0.5_N1000_reps300.csv"
YA3  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0_hc0.5_N1000_reps300.csv"
YX3  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.05_ho0.5_sc0_hc0.5_N1000_reps300.csv"
YA4  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.09_ho0.5_sc0.005_hc0.5_N1000_reps300.csv"
YX4  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.09_ho0.5_sc0.005_hc0.5_N1000_reps300.csv"
YA5  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.09_ho0.5_sc0_hc0.5_N1000_reps300.csv"
YX5  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.09_ho0.5_sc0_hc0.5_N1000_reps300.csv"
YA6  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.095_ho0.5_sc0.005_hc0.5_N1000_reps300.csv"
YX6  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.095_ho0.5_sc0.005_hc0.5_N1000_reps300.csv"
YA7  <-  "./output/data/simData/dataYATimeFixMutGrad2_sm0.1_delta0.095_ho0.5_sc0_hc0.5_N1000_reps300.csv"
YX7  <-  "./output/data/simData/dataYXTimeFixMutGrad2_sm0.1_delta0.095_ho0.5_sc0_hc0.5_N1000_reps300.csv"

    aDat   <-  YA3
    xDat   <-  YX3
    aDat2  <-  YA6
    xDat2  <-  YX6
    aData  <-  read.csv(aDat, header=TRUE)
    xData  <-  read.csv(xDat, header=TRUE)
    aData2  <-  read.csv(aDat2, header=TRUE)
    xData2  <-  read.csv(xDat2, header=TRUE)

    # Calculate additional variables for Plotting
    DeltaInv      <-     xData$deltaInv / aData$deltaInv
    DeltaInv2     <-     xData2$deltaInv / aData2$deltaInv
    DeltaFix      <-     xData$deltaFix / aData$deltaFix
    DeltaFix2     <-     xData2$deltaFix / aData2$deltaFix
    DeltaInvFix   <-  xData$deltaInvFix / aData$deltaInvFix
    DeltaInvFix2  <-  xData2$deltaInvFix / aData2$deltaInvFix
    tCycle        <-       xData$tCycle / aData$tCycle
    tCycle2       <-       xData2$tCycle / aData2$tCycle

    # relative size of X vs. Autosomal genome in 
    # D. melanogaster (see Adams et al. 2000)
    MbA  <-  23 + 5.4 + 11 + 21.4 + 24.4 + 8.2 + 8.2 + 28 + 3.1 + 1.2
    MbX  <-  20 + 21.8
    MbY  <-  40.9
    relMuA  <-  MbA / (MbA + MbX)
    relMuX  <-  MbX / (MbA + MbX)
    relMutRateA  <-  (aData$ua * relMuA) / aData$uy
    relMutRateX  <-  (aData$ux * relMuA) / aData$uy

    # Set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

        # Make the plot
# Panel 1: Time to complete coevolutionary cycle
        par(omi=rep(0.5, 4), mar = c(3,3,4,3), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.75,1.2), ylab='', xlab='', cex.lab=1.2)
#        plot(NA, axes=FALSE, type='n', main='',xlim = c(min(log10(relMutRate)),max(log10(relMutRate))), ylim = c(0.5,1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        abline(h=1,lty=2, lwd=2)
        points(tCycle ~ (relMutRate), pch=21, bg=transparentColor('dodgerblue',opacity=0.7), col='dodgerblue', ylim=c(0,2))
        points(tCycle2 ~ (relMutRate), pch=21, bg=transparentColor('tomato',opacity=0.85), col='tomato', ylim=c(0,2))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
#        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(T[X]/T[A])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(mu[i])/italic(mu[y]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste("Time to complete coevolutionary cycle")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        #legend
        legend(
              x       =  usr[2]*1,
              y       =  usr[4]*0.6,
              legend  =  c(
                          expression(paste(italic(s[m])~"="~"0.1;"~italic(delta)~"="~"0.050;"~italic(s[c])~"="~"0.005;"~italic(N)~"="~"1000;")),
                          expression(paste(italic(s[m])~"="~"0.1;"~italic(delta)~"="~"0.090;"~italic(s[c])~"="~"0.000;"~italic(N)~"="~"1000;"))),
              pch     =  21,
              col     =  c(transparentColor('dodgerblue', opacity=0.85),transparentColor('tomato', opacity=0.85)),
              pt.bg   =  c(transparentColor('dodgerblue', opacity=0.85),transparentColor('tomato', opacity=0.85)),
              cex     =  0.7,
              xjust   =  1,
              yjust   =  1,
              bty     =  'n',
              border  =  NA
        )


# Panel 3: Time between y and a/x fixation
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.5,1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        abline(h=1,lty=2, lwd=2)
        points(DeltaFix ~ relMutRate, pch=21, bg=transparentColor('dodgerblue',opacity=0.7), col='dodgerblue', ylim=c(0,2))
        points(DeltaFix2 ~ relMutRate, pch=21, bg=transparentColor('tomato',opacity=0.85), col='tomato', ylim=c(0,2))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
#        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(Delta["fix,X"]/Delta["fix,A"])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(mu[i])/italic(mu[y]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste("Time between y and a/x fixation")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

# Panel 4: Time from y invasion and a/x fixation
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0.5,1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # plot data
        abline(h=1,lty=2, lwd=2)
        points(DeltaInvFix ~ relMutRate, pch=21, bg=transparentColor('dodgerblue',opacity=0.7), col='dodgerblue', ylim=c(0,2))
        points(DeltaInvFix2 ~ relMutRate, pch=21, bg=transparentColor('tomato',opacity=0.85), col='tomato', ylim=c(0,2))
        # axes and labels
        axis(1, las=1)
        axis(2, las=1)
#        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(Delta["inv-fix,X"]/Delta["inv-fix,A"])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(mu[i])/italic(mu[y]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste("Time from y invasion and a/x fixation")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

}