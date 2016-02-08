#' Specification of form of the trend within the nim formula.
#'
#' The function specifying the trend are not meant to be used outside the formula of the \link{nim} function.
#'
#' @param ... The name of the covariate within the nim formula
#' @param h Width of the smoothing window (bandwidth) as a fraction of the number of the time steps included in the window (for \code{h <= 1}) or as a number of time-steps (for \code{h>1})
#' @param eval_at Either numerical indices of the time series where to calculate the local regression estimates, all (defaults) or 'all' for evaluation at each time step. The outcome (smoothed parameter) are linearly approximated in order to provide value for each time step.
#' @param w The kernel function - defaults to Epanechnikov kernel. Could be specified as a function of u: (-Inf, Inf) -> (0, 1), typically with largest values for u = 0.
#'
#' @return Returns the environment of the function.
#' @export
#' @name trend-spec
#' @section See also:
#' \link{nim}
#' @examples
#' data('precip_max')
#' head(precip_max)
#'
#' # indicate which columns contain the extremes to be fitted
#' extremes(precip_max) = 2:ncol(precip_max)
#'
#' # smooth trend in xi, with bandwidth h = 0.2
#' nim(xi ~ s(YR, h = 0.2), data = precip_max)
#'
#' # linear trend in xi and gamma
#' nim(xi + g ~ p(YR), data = precip_max)

s = function(..., h = .1, eval_at = all, w = function(u){pmax(1 - u^2, rep(0, length(u)))}) {
  environment()
  }

#' @export
#' @name trend-spec
p = function(...) {
  environment()
  }


#' Title
#'
#' @param data
#' @param value
#'
#' @return Either nothing,
#' @export
#'
#' @examples
"extremes<-" = function(data, value){
  structure(data, extremes = value)
}

extremes = function(...){
  UseMethod('extremes')
}

extremes.data.frame = function(data){
  if (is.null(attr(data, 'extremes'))) stop('Extremes not specified.')
  data[, attr(data, 'extremes'), drop = FALSE]
}

extremes.nim = function(nim){
  extremes(attr(nim, 'data'))
}

model_info = function(...){
  UseMethod('model_info')
}


model_info.formula = function(formula){
  type = formula[[length(formula)]][[1]]
  vars = all.vars(formula)
  cvrt = vars[!grepl('xi|g|k', vars)]
  pars = vars[!(vars %in% cvrt)]
  if (length(cvrt)==0) cvrt = 'I'
  if (type == 's'){
    h = eval(formula[[3]])$h
    #  eval_at = eval(formula[[3]])$eval_at
    #  if (is.function(eval_at) | is.character(eval_at)) eval_at = 1:length(xcov)
    #  eval_at = xcov[eval_at]
    w = eval(formula[[3]])$w
  }
  as.list(environment())
}

model_info.nim = function(nim){
  model_info.formula(attr(nim, 'call')[['formula']])
}

print.nim = function(nim){
  nfo = model_info(nim)
  if (nfo$type == 's') type = paste('smooth in', paste0(nfo$pars, collapse = ', '))
  if (nfo$type == 'p') type = paste('parametric in', paste0(nfo$pars, collapse = ', '))
  if (nfo$type == 1) type = 'stationary'
  cat("- \n")
  cat("(Non)-stationary index-flood model\n")
  cat("- \n")
  cat("trend:\t", type)
  cat("\nlocal parameters\n")
  print(nim$.xi()[1, ])

  cat("\nregional parameters\n")
  print(apply(nim$REG, 2, range))
}

logLik.nim = function(nim){
  attr(nim, 'logLik')
}

# TDD
# fitted.nprgev = function(x){
#   attr(x, 'fitted')
# }

# TDD
# AD = function(x, ...){
#   UseMethod('AD', x)
# }

residuals.nim = function(nim, type = 'Gumbel'){
  std = function(xx)qnorm(exp(-exp(- xx )))
  r = attr(nim, 'resid')
  nfo = model_info(nim)
  rr = switch (type,
               'Gumbel' = r,
               'Normal' = data.table(r[[nfo$cvrt]], r[, sapply(.SD, std), .SDcols = 2:ncol(r)])
  )
  rr = data.frame(rr)
  extremes(rr) = attr(attr(nim, 'data'), 'extremes')
  return(rr)

}

resid.nim = function(nim, type = 'Gumbel'){
  residuals(nim, type)
}

sample = function(...){
  UseMethod('sample')
}

sample.default = base::sample


#' Create nims object.
#'
#' The function concatenates nim objects (typically the result of bootstrapping) into a nims object. Various helper functions are provided for comparison of various characteristics of individual nims.
#'
#' @param ...
#'
#' @return `nims` object
#' @export
#'
#' @examples
#'
nims = function(...){

  if (all(sapply(list(...), class )=='nim')) {n = structure(list(...), names = as.character(match.call()[1:length(list(...)) + 1]))} else {
  if (is.list(...)) n = structure(...)}

  nfo = sapply(n, model_info)
  if (length(unique((nfo['cvrt', ])))!=1) stop('All covariates must have same names.') else (cvrt = nfo[['cvrt', 1]])
  structure(n, class = 'nims', cvrt = cvrt)
}

lengths = function(...){
  UseMethod('lengths')
}

lengths.default = base::lengths

lengths.nims = function(nims, sites = FALSE){
  if (!sites) return(sapply(nims, function(x)nrow(x$REG)))
  sapply(nims, function(x)length(x$XI))
}

extract_from_nims = function(nims, FUN, return = 'data.table', melt = FALSE, sites = FALSE, ...){
  l = lapply(nims, FUN)
  if (return == 'list') return(l)
  if (return == 'data.table'){
    dt = data.table(NIM = data.table(names(nims), lengths(nims, sites = sites))[, rep(V1, each = V2), by = V1][, V1], do.call(rbind, l))
    if (melt) return(melt(dt, id.vars = if (!sites) {c('NIM', attr(nims, 'cvrt'))} else {'NIM'}, ...)) else return(dt)
  }
  stop('`return` must be one of "data.table" or "list".')
}

residuals.nims = function(nims, return = 'data.table', melt = FALSE, ...){
  extract_from_nims(nims, resid, return, melt, ...)
}

regional = function(...){
  UseMethod('regional')
}

regional.nim = function(nim, transform_G = TRUE){
  o = nim$REG
  if (transform_G) o$G = exp(o$G)
  data.table(o)
}

regional.nims = function(nims, return = 'data.table', melt = FALSE, ...){
  extract_from_nims(nims, regional, return = return, melt = melt, ...)
}

atsite = function(...){
  UseMethod('atsite')
}

atsite.nim = function(nim){
  nim$XI
}

atsite.nims = function(nims, return = 'data.table', melt = FALSE, ...){
  extract_from_nims(nims, atsite, return = return, melt = melt, sites = TRUE, ...)
}

params2data.nims = function(nims, return = 'data.table', melt = FALSE, ...){
  extract_from_nims(nims, params2data, return = return, melt = melt, ...)
}

AD.nims = function(nims, return = 'data.table', melt = FALSE, ...){
  e = extract_from_nims(nims, AD, return = return, melt = FALSE, sites = TRUE, ...)
  setnames(e, 'V2', 'ad')
  if (!melt) return(dcast(e, NIM ~ ID)) else (return(e))
}
