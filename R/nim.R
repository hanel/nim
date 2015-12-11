#' nim: A package for fitting and assessment of non-stationary index floods models
#'
#' The nim package provides functions for fitting, uncertainty and goodness-of-fit assessment of Non-stationary Index-flood Models (nims) with smooth or parametric trends in the parameters of the GEV model together with some classes and methods for easy manipulation and assessment of results.
#' @author Martin Hanel \email{hanel@@fzp.czu.cz}
#' @references Hanel, M., Buishand, T. A., Ferro, C. A. T. (2009)  A non-stationary index-flood model for precipitation extremes in transient Regional Climate Model simulations, Journal of Geophysical Research., 114(D15107). \href{http://onlinelibrary.wiley.com/doi/10.1029/2009JD011712/abstract}{link}
#' @docType package
#' @name nim-package
NULL
#> NULL


#' Fitting the non-stationary index-flood models
#'
#' @param formula one or two sided formula describing the model, see Details
#' @param data data.frame including the covariate and the extremes to be fitted
#' @param method method used for the parameter optimization see ?optim
#' @param tol minimum improvement of the log-likelihood between the iteration steps, if the difference between two succesive iteration steps is lower, the iteration ends
#' @param ... other parameters passed to \code{optim(control = list(...))}
#'
#' @return Object of class \link{nim} - i.e. simply
#' @export
#' @section Details:
#' The function is a simple wrapper for fitting non-stationary index-flood models with parametric (currently only linear) and smoothing trends in the parameters of the GEV model. Stationary index-flood models can be fitted also. For the data (\code{data.frame}) that are used within the function, the attribute \code{extremes} (see \link{extremes}) has to be set. The trend can be either smooth (\link{s}) or parametric (\link{p})n.
#' @section See also:
#' \link{sample}, \link{fit}, ...
#' @examples
#' data('precip_max')
#' head(precip_max)
#' # indicate which columns contain the extremes to be fitted
#' extremes(precip_max) = 2:ncol(precip_max)
#'
#' # stationary model
#' nim( ~1, data = precip_max)
#'
#' # smooth trend in xi, with bandwidth h = 0.2
#' nim(xi ~ s(YR, h = 0.2), data = precip_max)
#'
#' # smooth trend in all parameters, default bandwidth
#' nim(xi + g + k ~ s(YR), data = precip_max)
#'
#' # linear trend in xi and gamma
#' nim(xi + g ~ p(YR), data = precip_max)
nim = function(formula =   ~ 1, data, method = 'Nelder-Mead', tol = 0.1, ...){

  ext = attr(data, 'extremes')
  data = data.frame(data)
  extremes(data) = ext
  estSmoothReg = function(){

    if (h > 1) h = h / length(xcov)
    r = list()
    for (i in eval_at){

      W = w((xcov - i)/ (d * h))
      win = dmx[W>0,]
      xpr = list(t0 = i, t = xcov[W>0])

      .xt = xpr$t - xpr$t0
      .xi = THETA$.xi(xpr)

      params = c(xi = mean(THETA$REG$XI[W > 0]), xi1 = 0.000000001, g = mean(THETA$REG$G[W > 0]), g1 = 0.000000001, k = mean(THETA$REG$K[W > 0]))
      start = params[sort(charmatch( c(pars, paste0(pars, '1')), names(params) ))]
      fp = function(par)cReg(p = par, xpar = .xt, xi = .xi, g = THETA$REG$G[1], k = THETA$REG$K[1])

      f = cffx(data = as.matrix(win), fpar = fp, start = start, xpar = .xt, method = method, w = W[W>0])

      params[names(f$estimate)] = f$estimate
      r[[as.character(i)]] = params
    }

    params = list(XI = approx(x = eval_at, xout = THETA$REG[[cvrt]], y = sapply(r, function(x)x[1]))$y, G =approx(x = eval_at, xout = THETA$REG[[cvrt]], y = sapply(r, function(x)x[3]))$y, K =approx(x = eval_at, xout = THETA$REG[[cvrt]], y = sapply(r, function(x)x[5]))$y)
    return(params)
  }

  estLinReg = function(){

    .xi = THETA$.xi()
    params = c(xi = mean(THETA$REG$XI), xi1 = 0, g = mean(THETA$REG$G), g1 = 0, k = mean(THETA$REG$K), k1 = 0)
    start = params[sort(charmatch( c(pars, paste0(pars, '1')), names(params) ))]
    fp = function(par)cReg(p = par, xpar = xcov, xi = .xi, g = THETA$REG$G[1], k = THETA$REG$K[1])
    f = cffx(data = as.matrix(dmx), fpar = fp, start = start, xpar = xcov, method = method)

    params[names(f$estimate)] = f$estimate
    list(XI = params['xi'] + xcov * params['xi1'], G = params['g'] + xcov * params['g1'], K = params['k'] + xcov * params['k1'])

  }

  type = formula[[length(formula)]][[1]]
  vars = all.vars(formula)
  cvrt = vars[!grepl('xi|g|k', vars)]
  pars = vars[!(vars %in% cvrt)]
  if (length(cvrt)==0) cvrt = 'I'
  xcov = if (type==1) (rep(1, nrow(data))) else (data[[cvrt]])
  dmx = extremes(data)#data[,!grepl(cvrt, colnames(data))]

  if (all(c('xi', 'g') %in% pars) & !('k' %in% pars)) cReg = cppREG_constK
  if (all(c('xi') %in% pars) & all(!(c('g', 'k') %in% pars))) cReg = cppREG_constGK

  if (type == 's'){
    estReg = estSmoothReg
    if (all(c('xi', 'g', 'k') %in% pars)) cReg = cppREG
    h = eval(formula[[3]])$h
    eval_at = eval(formula[[3]])$eval_at
    if (is.function(eval_at) | is.character(eval_at)) eval_at = 1:length(xcov)
    eval_at = c(1, eval_at, length(xcov))
    eval_at = eval_at[!duplicated(eval_at)]
    eval_at = xcov[eval_at]
    w = eval(formula[[3]])$w
    d = diff(range(xcov))
  }

    if (type == 'p'){
    estReg = estLinReg
    if (all(c('xi', 'g', 'k') %in% pars)) cReg = cppREG_lin
  }

  if (type == 1){
    estReg = function()THETA$REG
  }

  ini = lgev(as.matrix(dmx))
  THETA = list(XI = ini$location, .xi = function(xpar = list(t = THETA$REG[[1]])){matrix(THETA$XI, nrow = length(xpar$t), ncol = length(THETA$XI), byrow = TRUE)}, REG = data.frame(COV = xcov, XI = 1, G = log(ini$dispersion), K = (ini$shape)))
  names(THETA$REG)[1] = cvrt
  names(THETA$XI) = names(dmx)
  dev0 = Inf
  dev1 = nll.gev(par = THETA, fpar = lpar, dat = as.matrix(dmx))

  while (dev1 - dev0 < -tol){

   message('\ndeviance: ', dev1, '\t DIFF: ', dev1-dev0)
   dev0 = dev1

   message( '\testimating at-site xi')

   for (i in 1:ncol(dmx)){
     THETA$XI[i] = cffx(data = as.matrix(dmx[, i]), fpar = function(par){cppXI0(p = par, xi = THETA$REG$XI, g = THETA$REG$G, k = THETA$REG$K) }, start = c(THETA$XI[i]), xpar = xcov, method = 'Brent', lower = min(dmx), upper = max(dmx))$estimate
   }

   message( '\testimating regional parameters')

   if ((!('g' %in% pars) & !('k' %in% pars)) | (type == 1) ){

      est = cffx(data = as.matrix(dmx), fpar = function(par)lpar_gk(par, THETA), start = c(THETA$REG$G[1], THETA$REG$K[1]), xpar = xcov, method = method)$estimate
      THETA$REG$G = est[1]
      THETA$REG$K = est[2]
   }

   if (('g' %in% pars) & !('k' %in% pars) & ('xi' %in% pars)){

     THETA$REG$K = cffx(data = as.matrix(dmx), fpar = function(par)lpar_k(par, THETA), start = THETA$REG$K[1], xpar = xcov, method = 'Brent', lower = -.5, upper = .5)$estimate

   }

   est = estReg()
   THETA$REG$XI = est$XI
   THETA$REG$G = est$G
   THETA$REG$K = est$K

   dev1 = nll.gev(par = THETA, fpar = lpar, dat = as.matrix(dmx))

  }
  extremes(data) = ext
 out = structure(.Data = THETA,  logLik = -dev1, call = match.call(), class = 'nim', data = data)#, resid = provideResid(data, THETA), observed = data, h = h, call = match.call(), class = 'nprgev')
 structure(out, resid = provideResid(out))
 # out
}

