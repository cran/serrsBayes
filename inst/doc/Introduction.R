## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(1234)

## ----data, message=FALSE, fig.width=6-----------------------------------------
library(serrsBayes)
data("lsTamra", package = "serrsBayes")
wavenumbers <- lsTamra$wavenumbers
spectra <- lsTamra$spectra
plot(wavenumbers, spectra[1,], type='l', col=4, main="Raman Spectrum of TAMRA+DNA",
     xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)")
nWL <- length(wavenumbers)

## ----priors-------------------------------------------------------------------
peakLocations <- c(615, 631, 664, 673, 702, 705, 771, 819, 895, 923, 1014,
                   1047, 1049, 1084, 1125, 1175, 1192, 1273, 1291, 1307, 1351, 1388, 1390,
                   1419, 1458, 1505, 1530, 1577, 1601, 1615, 1652, 1716)
nPK <- length(peakLocations)
priors <- list(loc.mu=peakLocations, loc.sd=rep(25,nPK), scaG.mu=log(10) - (0.34^2)/2,
               scaG.sd=0.34, scaL.mu=log(10) - (0.4^2)/2, scaL.sd=0.4, noise.nu=5, noise.sd=20,
               bl.smooth=1, bl.knots=121, beta.exp=80)

## ----SMC, results='hide'------------------------------------------------------
data("result", package = "serrsBayes")
if(!exists("result")) {
  t1 <- Sys.time()
  result <- fitVoigtPeaksSMC(wavenumbers, spectra, priors, mcSteps=50, minPart = 7900, npart = 8000)
  result$time <- Sys.time() - t1
  save(result, file="Figure 2/result.rda")
}

## ----time---------------------------------------------------------------------
print(paste(format(result$time,digits=4)," for",length(result$ess),"SMC iterations."))

## ----ess, fig.show='hold', fig.width=3, fig.height=3--------------------------
plot.ts(result$ess, ylab="ESS", main="Effective Sample Size",
        xlab="SMC iteration", ylim=c(0,max(result$ess)))
abline(h=max(result$ess)/2, col=4, lty=2)
abline(h=0,lty=2)
plot.ts(result$accept, ylab="accept", main="M-H Acceptance Rate",
        xlab="SMC iteration", ylim=c(0,max(result$accept)))
abline(h=0.234, col=4, lty=2)
abline(h=0,lty=2)
plot.ts(result$times/60, ylab="time (min)", main="Elapsed Time", xlab="SMC iteration")
plot.ts(result$kappa, ylab=expression(kappa), main="Likelihood Tempering",
        xlab="SMC iteration")
abline(h=0,lty=2)
abline(h=1,lty=3,col=4)

## ----baseline, fig.show='hold'------------------------------------------------
samp.idx <- sample.int(length(result$weights), 50, prob=result$weights)
samp.mat <- resid.mat <- matrix(0,nrow=length(samp.idx), ncol=nWL)
samp.sigi <- samp.lambda <- numeric(length=nrow(samp.mat))
samp.loc <- colMeans(result$location)
spectra <- as.matrix(spectra)
plot(wavenumbers, spectra[1,], type='l', xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)")
for (pt in 1:length(samp.idx)) {
  k <- samp.idx[pt]
  samp.mat[pt,] <- mixedVoigt(result$location[k,], result$scale_G[k,],
                              result$scale_L[k,], result$beta[k,], wavenumbers)
  samp.sigi[pt] <- result$sigma[k]
  samp.lambda[pt] <- result$lambda[k]
  
  Obsi <- spectra[1,] - samp.mat[pt,]
  g0_Cal <- length(Obsi) * samp.lambda[pt] * result$priors$bl.precision
  gi_Cal <- crossprod(result$priors$bl.basis) + g0_Cal
  mi_Cal <- as.vector(solve(gi_Cal, crossprod(result$priors$bl.basis, Obsi)))

  bl.est <- result$priors$bl.basis %*% mi_Cal # smoothed residuals = estimated basline
  lines(wavenumbers, bl.est, col=2)
  lines(wavenumbers, bl.est + samp.mat[pt,], col=4)
  resid.mat[pt,] <- Obsi - bl.est[,1]
}
title(main="Baseline for TAMRA")
rug(samp.loc, col=4,lwd=2)

plot(range(wavenumbers), range(samp.mat), type='n', xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)")
abline(h=0,lty=2)
for (pt in 1:length(samp.idx)) {
  lines(wavenumbers, samp.mat[pt,], col=4)
}
title(main="Spectral Signature")
rug(samp.loc,col=4,lwd=2)

## -----------------------------------------------------------------------------
result$voigt <- result$FWHM <- matrix(nrow=nrow(result$beta), ncol=ncol(result$beta))
for (k in 1:nrow(result$beta)) {
  result$voigt[k,] <- getVoigtParam(result$scale_G[k,], result$scale_L[k,])
  f_G <- result$scale_G[k,]
  f_L <- result$scale_L[k,]
  result$FWHM[k,] <- 0.5346*f_L + sqrt(0.2166*f_L^2 + f_G^2)
}

## ----confint, warning=FALSE, echo=FALSE, results='asis'-----------------------
result$voigt <- result$FWHM <- matrix(nrow=nrow(result$beta), ncol=ncol(result$beta))
for (k in 1:nrow(result$beta)) {
  result$voigt[k,] <- getVoigtParam(result$scale_G[k,], result$scale_L[k,])
  f_G <- 2*result$scale_G[k,]*sqrt(2*log(2))
  f_L <- 2*result$scale_L[k,]
  result$FWHM[k,] <- 0.5346*f_L + sqrt(0.2166*f_L^2 + f_G^2)
}

voigtCI <- apply(result$voigt, 2, function(x) quantile(x,probs=c(0.025,0.975)))
fmtVoigtCI <- format(voigtCI,digits=2)
fwhmCI <- apply(result$FWHM, 2, function(x) quantile(x,probs=c(0.025,0.975)))
fmtFWHMci <- format(fwhmCI,digits=2)
ampCI <- apply(result$beta, 2, function(x) quantile(x,probs=c(0.025,0.975)))
fmtAmpCI <- format(ampCI,digits=1,scientific=FALSE)
locCI <- apply(result$location, 2, function(x) quantile(x,probs=c(0.025,0.975)))
fmtLocCI <- format(locCI,digits=3)

library(knitr)
tabCI <- cbind(paste0("[", fmtLocCI[1,], "; ", fmtLocCI[2,], "]"),
  paste0("[", fmtAmpCI[1,], "; ", fmtAmpCI[2,], "]"),
  paste0("[", fmtFWHMci[1,], "; ", fmtFWHMci[2,], "]"),
  paste0("[", fmtVoigtCI[1,], "; ", fmtVoigtCI[2,], "]"))
colnames(tabCI) <- c("Location (cm-1)", "Amplitude", "FWHM (cm-1)", "Voigt")
kable(tabCI, caption="95% highest posterior density intervals for pseudo-Voigt peaks",
      align = 'rrrr')

