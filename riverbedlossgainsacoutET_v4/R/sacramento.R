## hydromad: Hydrological Modelling and Analysis of Data
##
## Copyright (c) Felix Andrews <felix@nfrac.org>
##

#' Sacramento Soil Moisture Accounting model
#'
#' Sacramento Soil Moisture Accounting model.  Developed by the US National
#' Weather Service.
#'
#' This description of the model is given by Burnash (1995):
#'
#' \dQuote{The moisture accounting system utilized in the Sacramento Catchment
#' Model is a carefully structured representation of the catchment's soil
#' moisture storage system. It is based on using simple approximations of many
#' of those soil moisture processes which have been reported in the hydrologic
#' literature. The authors have organised these approximations in a manner
#' which would allow the determination of many catchment characteristics from
#' carefully selected portions of the catchment's hydrologic record. Inasmuch
#' as many of the catchment characteristics are related to the soil moisture
#' capabilities of the catchment, an intelligent application of the model start
#' with a good understanding of the three basic types of soil moisture which
#' can potentially influence catchment runoff conditions. These soil moisture
#' types are: (1) Hygroscopic Water, (2) Tension Water and (3) Free Water. }
#'
#' [...]
#'
#' \dQuote{Streamflow as computed by the Sacramento Catchment Model is the
#' result of processing precipiatation through an algorithm representing the
#' uppermost soil mantle identified as the upper zone and a deeper portion of
#' the soil mantle or lower zone. The algorithm computes runoff in five basic
#' forms. These are (1) direct runoff from permanant and temporary impervious
#' areas, (2) surface runoff due to precipitation occurring at a rate faster
#' than percolation and interflow can take place when both upper zone storages
#' are full, (3) interflow resulting from the lateral drainage of a temporary
#' free water storage, (4) supplemental base flow, and (5) primary base flow.}
#' (Burnash, 1995)
#'
#' The default parameter ranges were taken from Blasone et. al. (2008).
#'
#' Note that the Sacramento model potentially suffers from numerical
#' instabilities, which can be seen for example as discontinuities in output
#' and derivatives of outputs (see Hendrickson et al. 1988). Ideally, the
#' underlying differential equations of the model would be solved using a
#' numerically robust timestepping scheme (see Clark & Kavetski 2010). The
#' hydromad package makes use of an existing implementation. To help remedy the
#' numerical instability, the argument \code{min_ninc} has been added, which
#' defines the minimum number of inner loops used within each timestep. The
#' user is encouraged to test the effect of increasing \code{min_ninc} on their
#' dataset.
#'
#' @name sacramento
#' @aliases sacramento.sim
#' @param DATA time-series-like object with columns \code{P} (precipitation,
#' mm) and \code{E} (potential evapo-transpiration, mm, scaled by
#' \code{etmult}).
#' @param uztwm Upper zone tension water maximum capacity (mm).
#' @param uzfwm Upper zone free water maximum capacity (mm).
#' @param uzk Lateral drainage rate of upper zone free water expressed as a
#' fraction of contents per day.
#' @param pctim The fraction of the catchment which produces impervious runoff
#' during low flow conditions.
#' @param adimp The additional fraction of the catchment which exhibits
#' impervious characteristics when the catchment's tension water requirements
#' are met.
#' @param zperc Maximum percolation (from upper zone free water into the lower
#' zone) rate coefficient.
#' @param rexp An exponent determining the rate of change of the percolation
#' rate with changing lower zone water contents.
#' @param lztwm Lower zone tension water maximum capacity (mm).
#' @param lzfsm Lower zone supplemental free water maximum capacity (mm).
#' @param lzfpm Lower zone primary free water maximum capacity (mm).
#' @param lzsk Lateral drainage rate of lower zone supplemental free water
#' expressed as a fraction of contents per day.
#' @param lzpk Lateral drainage rate of lower zone primary free water expressed
#' as a fraction of contents per day.
#' @param pfree Direct percolation fraction from upper to lower zone free water
#' (the percentage of percolated water which is available to the lower zone
#' free water aquifers before all lower zone tension water deficiencies are
#' satisfied).
#' @param etmult Multiplier applied to \code{DATA$E} to estimate potential
#' evapotranspiration.
#' @param dt Length of each time step in days.
#' @param uztwc_0 Initial upper zone tension water contents as proportion of
#' \code{uztwm}
#' @param uzfwc_0 Initial upper zone free water content as proportion of
#' \code{uzfwm}
#' @param lztwc_0 Initial lower zone tension water content as proportion of
#' \code{lztwm}
#' @param lzfsc_0 Initial lower zone free water secondary as proportion of
#' \code{lzfsm}
#' @param lzfpc_0 Initial lower zone free water primary as proportion of
#' \code{lzfpm}
#' @param adimc_0 Initial additional impervious flow store, as proportion of
#' \code{uztwm+lztwm}
#' @param min_ninc Minimum number of inner iterations. This is a simple attempt
#' to improve numerical stability. See Details.
#' @param return_state to return time series of each state variable and flow
#' component
#' @return the simulated effective rainfall (\dQuote{total channel inflow}), a
#' time series of the same length as the input series.
#'
#' if \code{return_state=TRUE}, a list with components: \item{uztwc}{Upper zone
#' tension water content} \item{uzfwc}{Upper zone free water content}
#' \item{lztwc}{Lower zone tension water content} \item{lzfsc}{Lower zone free
#' secondary water content} \item{lzfpc}{Lower zone free primary water content}
#' \item{adimc}{Tension water contents of the additional impervious area}
#' \item{sett}{Cumulative total evapotranspiration} \item{se1}{Cumulative
#' evapotranspiration from upper zone tension water} \item{se3}{Cumulative
#' evapotranspiration from lower zone tension water} \item{se4}{Cumulative
#' evapotranspiration} \item{se5}{Cumulative evapotranspiration from riparian
#' zone} \item{roimp}{Runoff from impervious area} \item{sdro}{Six hour sum of
#' runoff (?)} \item{ssur}{Surface runoff} \item{sif}{Interflow}
#' \item{bfp}{Primary baseflow} \item{bfs}{Secondary baseflow}
#' \item{bfcc}{Channel baseflow (bfp+bfs)}
#' @author Felix Andrews \email{felix@@nfrac.org} and Joseph Guillaume, based
#' on code from the University of Arizona MOSCEM project
#' @seealso \code{\link{hydromad}(sma = "sacramento")} to work with models as
#' objects (recommended).
#' @references Burnash, R.J.C (1995). The NWS River Forecast System --
#' Catchment Modeling.  In: Vijay P. Singh (ed.), \emph{Computer models of
#' watershed hydrology.} Revised edition, Highlands Ranch, Colo. : Water
#' Resources Publications, c1995.  \url{http://www.wrpllc.com/books/cmwh.html}.
#'
#' Blasone, R., J.A. Vrugt, H. Madsen, D. Rosbjerg, B.A. Robinson, G.A.
#' Zyvoloski (2008). Generalized likelihood uncertainty estimation (GLUE) using
#' adaptive Markov Chain Monte Carlo sampling. \emph{Advances in Water
#' Resources} 31, pp. 630-648.
#'
#' Hendrickson, Jene' D., Soroosh Sorooshian, and Larry E. Brazil (1988)
#' Comparison of Newton-Type and Direct Search Algorithms for Calibration of
#' Conceptual Rainfall-Runoff Models. \emph{Water Resources Research} 24 (5):
#' 691-700.  \url{http://dx.doi.org/10.1029/WR024i005p00691}
#'
#' Clark, Martyn P., and Dmitri Kavetski (2010) Ancient Numerical Daemons of
#' Conceptual Hydrological Modeling: 1. Fidelity and Efficiency of Time
#' Stepping Schemes.” Water Resources Research 46 (10).
#' \url{http://dx.doi.org/10.1029/2009WR008894}
#' @keywords models
#' @examples
#'
#' ## view default parameter ranges:
#' str(hydromad.options("sacramento"))
#'
#' data(HydroTestData)
#' mod0 <- hydromad(HydroTestData, sma = "sacramento")
#' mod0
#'
#' ## simulate with some arbitrary parameter values
#' set.seed(2)
#' mod1 <- simulate(update(mod0, etmult = 0.01), 1,
#'   sampletype =
#'     "random"
#' )[[1]]
#'
#' testQ <- predict(mod1, return_state = TRUE)
#' xyplot(window(cbind(HydroTestData[, 1:2], sacramento = testQ), start = 100))
#' mod1
#'
#' ## show effect of increase/decrease in each parameter
#' parRanges <- hydromad.getOption("sacramento")
#' parsims <- mapply(
#'   val = parRanges, nm = names(parRanges),
#'   FUN = function(val, nm) {
#'     lopar <- min(val)
#'     hipar <- max(val)
#'     names(lopar) <- names(hipar) <- nm
#'     fitted(runlist(
#'       decrease = update(mod1, newpars = lopar),
#'       increase = update(mod1, newpars = hipar)
#'     ))
#'   }, SIMPLIFY = FALSE
#' )
#'
#' xyplot.list(parsims,
#'   superpose = TRUE, layout = c(1, NA),
#'   strip = FALSE, strip.left = TRUE,
#'   main = "Simple parameter perturbation example"
#' ) +
#'   latticeExtra::layer(panel.lines(fitted(mod1), col = "grey", lwd = 2))
#' useDynLib hydromad sma_sac
#' useDynLib hydromad sma_sac_state
#' export
sacramento.sim <-
  function(DATA,
           uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
           lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree,
           etmult = 1, dt = 1,
           uztwc_0 = 0.5, uzfwc_0 = 0.5,
           lztwc_0 = 0.5, lzfsc_0 = 0.5, lzfpc_0 = 0.5,
           adimc_0 = 0.5, min_ninc = 20,
           return_state = FALSE,
           state_S_ts=NA,
           state_S2_ts=NA,
           state_S3_ts=NA) {
    stopifnot(c("P", "E") %in% colnames(DATA))
    ## check values
    stopifnot(uztwm >= 0)
    stopifnot(uzfwm >= 0)
    stopifnot(uzk >= 0)
    stopifnot(0 <= pctim && pctim <= 1)
    stopifnot(adimp >= 0)
    stopifnot(zperc >= 0)
    stopifnot(lztwm >= 0)
    stopifnot(lzfsm >= 0)
    stopifnot(lzfpm >= 0)
    stopifnot(lzsk >= 0)
    stopifnot(lzpk >= 0)
    stopifnot(pfree >= 0)
    stopifnot(etmult >= 0)
    stopifnot(dt >= 0)
    

    xpar <-
      c(
        uztwm, uzfwm, uzk, pctim, adimp, zperc, rexp,
        lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree
      )

    P <- DATA[, "P"]
    E <- DATA[, "E"]
    ## skip over missing values
    bad <- is.na(P) | is.na(E)
    P[bad] <- 0
    E[bad] <- 0
    ## NOTE: there is no pure.R.code implementation
    if(is.na(state_S_ts[1])){
      state_S_ts<-as.double(rep(-1e6,NROW(DATA)-1))
      use_state_S_ts<-as.integer(0)
    } else {
      stopifnot(length(state_S_ts) == (NROW(DATA)-1))
      if(any(state_S_ts[!(state_S_ts==-1e6)]<0)) stop("shouldn't be UZTWC < 0")
      state_S_ts<-as.double(state_S_ts)
      use_state_S_ts<-as.integer(1)
    }
    if(is.na(state_S2_ts[1])){
      state_S2_ts<-as.double(rep(-1e6,NROW(DATA)-1))
      use_state_S2_ts<-as.integer(0)
    } else {
      stopifnot(length(state_S2_ts) == (NROW(DATA)-1))
      if(any(state_S2_ts[!(state_S2_ts==-1e6)]<0)) stop("shouldn't be UZFWC < 0")
      state_S2_ts<-as.double(state_S2_ts)
      use_state_S2_ts<-as.integer(1)
    }
    if(is.na(state_S3_ts[1])){
      state_S3_ts<-as.double(rep(-1e6,NROW(DATA)-1))
      use_state_S3_ts<-as.integer(0)
    } else {
      stopifnot(length(state_S3_ts) == (NROW(DATA)-1))
      if(any(state_S3_ts[!(state_S3_ts==-1e6)]<0)) stop("shouldn't be LZFPC < 0")
      state_S3_ts<-as.double(state_S3_ts)
      use_state_S3_ts<-as.integer(1)
    }
    if (return_state) {
      states <- .C(sma_sac_state,
        as.double(P),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(xpar),
        as.double(etmult),
        as.double(dt),
        U = double(NROW(DATA)),
        uztwc = double(NROW(DATA)),
        uzfwc = double(NROW(DATA)),
        lztwc = double(NROW(DATA)),
        lzfsc = double(NROW(DATA)),
        lzfpc = double(NROW(DATA)),
        adimc = double(NROW(DATA)),
        sett = double(NROW(DATA)),
        se1 = double(NROW(DATA)),
        se3 = double(NROW(DATA)),
        se4 = double(NROW(DATA)),
        se5 = double(NROW(DATA)),
        roimp = double(NROW(DATA)),
        sdro = double(NROW(DATA)),
        ssur = double(NROW(DATA)),
        sif = double(NROW(DATA)),
        bfp = double(NROW(DATA)),
        bfs = double(NROW(DATA)),
        bfcc = double(NROW(DATA)),
        as.double(uztwc_0), as.double(uzfwc_0),
        as.double(lztwc_0), as.double(lzfsc_0), as.double(lzfpc_0),
        as.double(adimc_0),
        as.integer(min_ninc),
        state_S_ts=state_S_ts,
        use_state_S_ts=use_state_S_ts,
        state_S2_ts=state_S2_ts,
        use_state_S2_ts=use_state_S2_ts,
        state_S3_ts=state_S3_ts,
        use_state_S3_ts=use_state_S3_ts,
        NAOK = FALSE #PACKAGE = "hydromad"
      )
      for (i in 7:25) attributes(states[[i]]) <- attributes(P)
      ans <- do.call(cbind, states[7:25])
      if(any(is.infinite(ans[,1])) | any(is.na(ans[,1]))){
        ans<-NA
      }
      return(ans)
    } else {
      U <- .C(sma_sac,
        as.double(P),
        as.double(E),
        as.integer(NROW(DATA)),
        as.double(xpar),
        as.double(etmult),
        as.double(dt),
        U = double(NROW(DATA)),
        as.double(uztwc_0), as.double(uzfwc_0),
        as.double(lztwc_0), as.double(lzfsc_0), as.double(lzfpc_0),
        as.double(adimc_0),
        as.integer(min_ninc),
        state_S_ts=state_S_ts,
        use_state_S_ts=use_state_S_ts,
        state_S2_ts=state_S2_ts,
        use_state_S2_ts=use_state_S2_ts,
        state_S3_ts=state_S3_ts,
        use_state_S3_ts=use_state_S3_ts,
        NAOK = FALSE #, PACKAGE = "hydromad"
      )$U
      ## make it a time series object again
      attributes(U) <- attributes(P)
      ## re-insert missing values
      U[bad] <- NA
      if(any(is.infinite(U)) | any(is.na(U))){
        U<-NA
      }
      return(U)
    }
  }

sma_sac<-"sma_sac"
sma_sac_state<-"sma_sac_state"

sacramento.ranges <- function() {
  list(
    uztwm = c(1, 150),
    uzfwm = c(1, 150),
    uzk = c(0.1, 0.5),
    pctim = c(0.000001, 0.1),
    adimp = c(0, 0.4),
    zperc = c(1, 250),
    rexp = c(0, 5),
    lztwm = c(1, 500),
    lzfsm = c(1, 1000),
    lzfpm = c(1, 1000),
    lzsk = c(0.01, 0.25),
    lzpk = c(0.0001, 0.25),
    pfree = c(0, 0.6)
  )
}

#' @import utils
utils::globalVariables(c("sma_sac_state"))
