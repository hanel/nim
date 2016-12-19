lgev = function(data,era=0,K=NA){

  require(lmom)

  mi=vector('numeric',dim(data)[2])
  ga=0
  ks=0
  if (era==0) (era=dim(data)[1])

 # L=Lmoments(data[1:era,],4)
  L = t(apply(data[1:era, , drop = FALSE], 2, function(x) samlmu(x, ratios = FALSE)))
  t3=L[,3]/L[,2]
  cc=2/(3+t3)-log(2)/log(3)
  if (is.na(K)) (k=7.8590*cc+2.9554*cc^2) else (k=K)
  a=(L[,2]*k)/((1-2^(-k))*gamma(1+k))
  mi=L[,1]-a*(1-gamma(1+k))/k
  GA=a/mi
  KK=-k
  Data=c(data[1:era,])


  #L=Lmoments(Data,4)
  L=samlmu(unlist(data), ratios = FALSE)
  t3=L[3]/L[2]
  cc=2/(3+t3)-log(2)/log(3)
  if (is.na(K)) (k=7.8590*cc+2.9554*cc^2) else (k=K)
  a=(L[2]*k)/((1-2^(-k))*gamma(1+k))
  mmi=L[1]-a*(1-gamma(1+k))/k
  ks=-k
  sig=a
  ga=sig/mmi

  if(ga<0) ga = sig # more clever solution?

  list(location=mi,j.location=mmi,dispersion=ga,shape=ks,scale=sig,Dispersion=GA,Shape=KK)
}

nll.gev = function(par,fpar,xpar,dat=NA,opar=NA,theta=NA,narm=T) {
  #browser()
  if (any(!is.na(dat))) (data=dat)
  if (any(is.na(opar))) (pmat <- fpar(par, xpar)) else (pmat <- fpar(par, xpar,opar,theta))
  loc <- pmat$loc
  scale <- pmat$scale
  shape <- pmat$shape
  if(any(scale <= 0)) return(1e+20)
  gumbel <- (abs(shape) < 1e-06)
  #browser()
  y <- (data - loc) / scale
  z <- 1 + shape * y
  if(any(z <= 0, na.rm = TRUE)) return(1e+20)
  nll <- (1 + 1 / shape) * log(z) + z^(-1 / shape)
  nll[gumbel] <- y[gumbel] + exp(-y[gumbel])
  if (narm) (out=sum(nll + log(scale), na.rm = TRUE)) else ({
    if (any(is.na(nll+log(scale)))) (out=1e+20) else(out=sum(nll + log(scale), na.rm = TRUE))
  })
}

prob = function(x){
  (rank(x)-.3) / (length(x)+.4)
}


params2data = function(...){
  UseMethod('params2data')
}

params2data.nim = function(nim, data = NULL){

  if (is.null(data)) data = attr(nim, 'data')
  if (is.null(data)) stop('Data must be provided.')
  i = attr(attr(nim, 'data'), 'extremes')
  nfo = model_info(nim)

  if (nfo$type==1) {
    data = data.frame(I = 1:nrow(data), extremes(data))
    extremes(data) = 2:ncol(data)
    nim$REG[['I']] = 1:nrow(data)
    }

  data = data[, names(data) %in% c(names(extremes(data)), nfo$cvrt)]
  mx = data.table(melt(data, id.vars = nfo$cvrt, variable.name = 'ID'))
  #setnames(mx, nfo$cvrt, 'COV')

  xi = data.table(COV = data[[nfo$cvrt]], nim$.xi())
  setnames(xi, names(xi), names(data))
 # setnames(xi, nfo$cvrt, 'COV')
  mxi = melt(xi, id.vars = nfo$cvrt, variable.name = 'ID', value.name = 'XI0')
  mxi = mxi[data.table(nim$REG), on = nfo$cvrt]
  mxi[, X:=XI0 * XI]
  m = mx[mxi, on = c(nfo$cvrt, 'ID')]
  m = m[, .(eval(parse(text = nfo$cvrt)), ID, value, XI = X, G, K, S = X * exp(G))]
  setnames(m, 1, nfo$cvrt)
  m
}


provideResid = function(nim, data = NULL){

  nfo = model_info(nim)
  m = params2data(nim, data)
  if (nfo$type == 1) {
    m[, I:=1:.N, by = ID]
  }
  m[, RESID:= (1 / K) * log ( 1 + (K / exp(G) * (value/XI - 1)  ) ) ]
  s = m[, .(eval(parse(text = nfo$cvrt)), ID, RESID)]
  setnames(s, 1, nfo$cvrt)
  res = dcast.data.table(s, eval(parse(text = nfo$cvrt)) ~ ID, value.var = 'RESID')
  setnames(res, 'nfo', nfo$cvrt)
  res
}

#' Create bootstrap sample based on fitted nim
#'
#' @param nim fitted model
#' @param length number of bootstrap samples
#' @param type one of \code{parametric_average_cor} (the default), \code{parametric} and \code{nonparametric}
#'
#' @return list of bootstrap samples
#' @export sample.nim
#'
#' @examples
sample.nim = function(nim, length = 1, type = 'parametric_average_cor', impute_NA = TRUE){

  nfo = model_info(nim)
  out = list()
  obs = attr(nim, 'data')
  if (nfo$cvrt=='I') {
    obs[[names(obs)[-attr(obs, 'extremes')]]] = 1
    names(obs)[-attr(obs, 'extremes')] = 'I'
  }#obs$I = 1
  #attr(obs, 'extremes') = c(attr(obs, 'extremes'), -which(colnames(obs)=='I'))
  rsd = copy(obs)

 # if (nfo$cvrt=='I') obs$I = 1

  if (type %in% c('parametric', 'parametric_average_cor') ){
    #cc = cov(resid(nim, type = 'Normal')[, 2:ncol(obs)])
    cc = cov(resid(nim, type = 'Normal')[, attr(obs, 'extremes'), drop = FALSE], use = 'pairwise.complete.obs')

    if (type == 'parametric_average_cor'){
      co = cov2cor(cc)
      co[] = mean(co[upper.tri(co)], na.rm = TRUE)
      diag(co) = 1
      sdMat = diag(sqrt(diag(cc)))
      cc = sdMat %*% co %*% t(sdMat)
    }

    sam = function(...){
      sresid = data.table(rmvnorm(nrow(nim$REG), sigma = cc, method = 'chol'))
      #sresid = data.frame(COV = obs[[nfo$cvrt]], sresid[, sapply(.SD, function(xx) (-log(-log(pnorm(xx) ))))])
      sresid = sresid[, sapply(.SD, function(xx) (-log(-log(pnorm(xx) ))))]
      rsd[, attr(obs, 'extremes')] = sresid

      #names(sresid) = c(nfo$cvrt, names(extremes(obs)))#names(obs)
      #extremes(sresid) = attr(obs, 'extremes')
      m = params2data(nim, rsd)
      m[, value := XI * (1 + exp(G) * ( (exp(K * value) -1) / K ) )]
      s = m[, .(eval(parse(text = nfo$cvrt)), ID, value)]
      setnames(s, 1, nfo$cvrt)
      res = dcast.data.table(s, eval(parse(text = nfo$cvrt)) ~ ID, value.var = 'value')
      setnames(res, 'nfo', nfo$cvrt)
      #res # TDD - what to do with NAs ???!!!
      if (impute_NA) {data.table(res[, 1, with = FALSE], as.matrix(res[, -1, with = FALSE]) * (extremes(nim)/extremes(nim)))} else {res} ## simplest way - copy the NA structure from data
    }
  }

  # TDD opravit - zobecnit pro stacionarni model
  if (type == 'nonparametric'){

    re = data.table(extremes(resid(nim, type = 'Gumbel')))
    #setkeyv(re, nfo$cvrt)
    sam = function(...){
      sresid = re[sample(1:nrow(re), nrow(re), replace = TRUE), ]
      rsd[, attr(obs, 'extremes')] = sresid

      #sresid$COV = nim$REG$COV
      #sresid = data.frame(sresid)
      #names(sresid) = names(obs)
      #extremes(sresid) = attr(obs, 'extremes')
      m = params2data(nim, rsd)
      m[, value:=XI * (1 + exp(G) * ( (exp(K * value) -1) / K ) )]
      res = dcast.data.table(m[, .(eval(parse(text = nfo$cvrt)), ID, value)], eval(parse(text = nfo$cvrt)) ~ ID, value.var = 'value')
      setnames(res, 'nfo', nfo$cvrt)
      res
    }

  }

  out = mapply(sam, 1:length, SIMPLIFY = FALSE)
  names(out) = paste0('BSP_', 1:length)

  out = lapply(out, function(x){
    x = data.frame(x)
    extremes(x) = attr(obs, 'extremes')
    return(x)})

  c(list(BSP_0 = obs), out)
}

ad = function(sgv){
  sgv = sgv[!is.na(sgv)]
  nr = length(sgv)
  u=sort(exp(-exp(-(sgv))))
  -nr-1/nr*sum((2*1:nr-1)*log(u)+(2*nr-2*1:nr+1)*log(1-u))
}

fit = function(nim, smp, verbose = FALSE, mc.cores = 2){#, pullData = TRUE){

  verb = if (!verbose) suppressMessages else identity
  RES = list()
  nfo = model_info(nim)
  cl = attr(nim, 'call')
 # RES[['BSP_0']] = nim

  # for (i in 1:length(smp)){
  #   message('sample ', i)
  #   cl$data = smp[[i]]
  #   verb({RES[[paste0('BSP_', i)]] = eval(cl)})
  # }
  message('Progress reporting of forked tasks is not supported in RStudio. To enable', immediate. = TRUE)
  r = mclapply(1:length(smp), function(i){
       message('sample ', i)
       cl$data = smp[[i]]
       eval(cl)
  }, mc.cores = 8)

  names(r) = paste0('BSP_', 1:length(smp))

  structure(c(RES, r), class = 'nims')

  # res = list()
  # res$results = RES
  # res[[nfo$cvrt]] = RES[[1]]$REG[[nfo$cvrt]]
  # res$ID = names(extremes(nim))
  # res$residuals = function(melt = FALSE, type = 'Gumbel', ...){
  #   p = data.table(BSP = paste0('BSP_', rep(1:length(res$results) - 1, each = length(res[[nfo$cvrt]]))), do.call(rbind, lapply(res$results, resid, type = type)))
  #   if (melt) return(melt(p, id.vars = c('BSP', nfo$cvrt), variable.name = 'ID', ...)) else return(p)
  # }
  # res$REG = function(melt = FALSE, ...){
  #   p =  data.table(BSP = paste0('BSP_', rep(1:length(res$results) - 1, each = length(res[[nfo$cvrt]]))), do.call(rbind, lapply(res$results, function(x)x$REG)))
  #   if (melt) return(melt(p, id.vars = c('BSP', nfo$cvrt), variable.name = 'PAR', ...)) else return(p)
  # }
  # res$atsite = function(melt = FALSE, ...){
  #   p = data.table(BSP = paste0('BSP_', 1:length(res$results)-1), do.call(rbind, lapply(res$results, function(x)x$XI)))
  #   setnames(p, 2:ncol(p), res$ID)
  #   if (melt) return(melt(p, id.vars = c('BSP'), variable.name = 'ID', value.name = 'XI', ...)) else return(p)
  # }
  # res$params2data = function(melt = FALSE, ...){
  #   p = data.table(BSP = paste0('BSP_', 1:length(res$results)-1), do.call(rbind, lapply(res$results, params2data)))
  #   if (melt) return(melt(p, id.vars = c('BSP', nfo$cvrt, 'ID', 'value'), variable.name = 'PAR', ...)) else return(p)
  # }
  # class(res) = c('fitted_sample')
  # return(res)
}

AD = function(...){
  UseMethod('AD')
}

AD.nim = function(nim){
  data.table(ID = names(extremes(nim)), apply(extremes(resid(nim)), 2, nim::ad))
}

AD.fitted_sample = function(fit){
  res = data.frame(t(sapply(f$results, AD)))
  res = data.frame(BSP = rownames(res), res, row.names = NULL)
  res
}



# TDD
# provideFitted = function(data, THETA){
#   mx = melt(data, id.vars = 1, variable.name = 'INDC')
#   mx[, p:= prob(value), by = INDC]
#   mu = data.table(YR = cmx$YR, THETA$.mu())
#   setnames(mu, names(mu), names(cmx))
#   mmu = melt(mu, id.vars = 1, variable.name = 'INDC', value.name = 'MU0')
#   setkey(mmu, YR)
#   setkey(THETA$REG, YR)
#   muu = mmu[THETA$REG]
#   muu[, MU:=MU0 * M]
#   setkey(muu, YR, INDC)
#   setkey(mx, YR, INDC)
#   m = mx[muu]
#   m[, FITTED:= evd::qgev(p, loc = MU, scale = MU * exp(G), shape = K), by = 1:nrow(m)]
#   dcast.data.table(m[, .(YR, INDC, FITTED)], YR ~ INDC, value.var = 'FITTED')
# }


# TDD
# predict.nprgev = function(obj, mx){
#   mmx = melt(mx, id.var = 'YR', variable.name = 'INDC')
#   setkey(obj$REG, YR)
#   setkey(mmx, YR)
#   mmx = obj$REG[mmx]
#   mmx[, sval:=value/M]
#   mmx[, MU0:= lgev(matrix(sval))$location, by = INDC]
#   mmx[, p:=prob(value), by = INDC]#evd::pgev(value, loc = MU0 * M, scale = (MU0 * M) * exp(G), shape = K), by = 1:nrow(mmx)]
#   mmx[, PRED:= evd::qgev(p, loc = MU0 * M, scale = (MU0 * M) * exp(G), shape = K), by = 1:nrow(mmx)]
#
#   ll = mmx[, nll.gev.plain(loc = MU0 * M, scale = (MU0 * M) * exp(G), shape = K, value = value)]
#   aa = mmx[,  ad((1 / K) * log ( 1 + (K / exp(G) * (value/(MU0 * M) - 1)  ) ))]
#   rmse = mmx[, mean(sqrt((PRED-value)^2)) ]
#   out = data.table(dcast.data.table(mmx[, .(YR, INDC, PRED)], YR ~ INDC, value.var = 'PRED'))
#   structure(.Data = list(out), logLik = ll, AD = aa, RMSE = rmse, class = 'predict.nprgev')
# }


nll.gev.plain = function(loc, scale, shape, value){

  if(any(scale <= 0)) return(1e+20)
  gumbel <- (abs(shape) < 1e-06)
  y <- (value - loc) / scale
  z <- 1 + shape * y
  if(any(z <= 0, na.rm = TRUE)) return(1e+20)
  nll <- (1 + 1 / shape) * log(z) + z^(-1 / shape)
  nll[gumbel] <- y[gumbel] + exp(-y[gumbel])
  if (any(is.na(nll+log(scale)))) (out=1e+20) else(out=sum(nll + log(scale), na.rm = TRUE))
  out

}



# getREG = function(x, ...){
#   UseMethod('getREG', x)
# }
#

# getAtSite = function(x, ...){
#   UseMethod('getAtSite', x)
# }

# getParDat = function(x, ...){
#   UseMethod('getParDat', x)
# }

# getREG.boot_nprgev = function(x, melt = FALSE, ...){
#   p =  data.table(BSP = paste0('BSP_', rep(1:length(x$results) - 1, each = length(x$YR))), do.call(rbind, lapply(x$results, function(y)y$REG)))
#   if (melt) return(melt(p, id.vars = c('BSP', 'YR'), ...)) else return(p)
# }

# residuals.boot_nprgev = function(x, melt = FALSE, type = 'Gumbel', ...){
#   p = data.table(BSP = paste0('BSP_', rep(1:length(x$results) - 1, each = length(x$YR))), do.call(rbind, lapply(x$results, resid, type = type)))
#   if (melt) return(melt(p, id.vars = c('BSP', 'YR'), ...)) else return(p)
# }

# resid.boot_nprgev = function(x, ...){
#   residuals(x, ...)
# }

# getAtSite.boot_nprgev = function(x, melt = FALSE, ...){
#   p = data.table(BSP = paste0('BSP_', 1:length(x$results)-1), do.call(rbind, lapply(x$results, function(xx)xx$MU)))
#   setnames(p, 2:ncol(p), x$ID)
#   if (melt) return(melt(p, id.vars = c('BSP'), ...)) else return(p)
# }

# getParDat = function(x, melt = FALSE, ...){
#   p = data.table(BSP = paste0('BSP_', 1:length(x$results)-1), do.call(rbind, lapply(x$results, params2data)))
#   if (melt) return(melt(p, id.vars = c('BSP', 'YR', 'INDC', 'value'), value.name = 'parVal', ...)) else return(p)
# }


extremity = function(nim, type = 'fitted', return_period = FALSE){
  d = attr(nim, 'data')
  r = resid(nim)
  switch(type,
    'fitted' = {
      o = apply(extremes(r), 2, function(x) evd::pgev(x))
    },
    'empirical' = {
      o = apply(extremes(r), 2, function(x) prob(x))
    })
  if (return_period) o = apply(o, 2, function(x) (1/(1-x)))
  d[, attr(d, 'extremes')] = o
  d
}

# quantile.nim = function(nim, p = NULL, T = c(2, 5, 10, 50), at_site = TRUE){
#
#   cl = match.call()
#   qO = Vectorize(function(p){
#     outer(XI , (1 - rg$G/rg$K * (1 - ( - log(p) ) ^ (-rg$K) )))
#   }  )
#
#   XI = if (at_site == TRUE) (nim$XI) else (1)
#   if (is.null(p)) p = 1 - 1 / T
#   rg = nim$REG
#   rg$G = exp(rg$G)
#   res = data.frame(qO(p))
#   #browser()
#   names(res) = if (is.null(cl$p)) (paste0(eval(cl$T), 'yr')) else (paste0(eval(cl$p), '%'))
#   res = data.frame(nim$REG[[1]], res, check.names = FALSE)
#   names(res)[1] = model_info(nim)$cvrt
#   res
# }

quantile.nim = function(nim, p = NULL, T = c(2, 5, 10, 50), at_site = TRUE){

  cl = match.call()
  qO = Vectorize(function(p){
    r = data.frame(t(outer(XI , (1 - rg$G/rg$K * (1 - ( - log(p) ) ^ (-rg$K) )))))
    if (model_info(nim)$cvrt == 'I') {
      r = data.frame(I = 1, r[1, ])
      names(r) = c(names(r)[1], names(XI))#if (!at_site) {names(r)[2] = 'regional'}
      r
    } else {
      r = data.frame(nim$REG[[1]], r, check.names = FALSE)
      names(r)[1] = model_info(nim)$cvrt
      r
    }
  }, SIMPLIFY = FALSE  )

  XI = if (at_site == TRUE) (nim$XI) else (structure(1, names = 'regional') )
  if (is.null(p)) p = 1 - 1 / T
  rg = nim$REG
  rg$G = exp(rg$G)
  res = qO(p)
  #  browser()
  names(res) = if (is.null(cl$p)) (paste0(eval(cl$T), 'yr')) else (paste0(eval(cl$p)*100, '%'))
  res = rbindlist(res, idcol = 'q')
  #res = data.frame(nim$REG[[1]], res, check.names = FALSE)
  #names(res)[1] = model_info(nim)$cvrt
  res
}

quantile.nims = function(nims, p = NULL, T = c(2, 5, 10, 50), at_site = TRUE){
  r = lapply(nims, function(x) melt(quantile(x, p, T, at_site), id.vars = 1:2))
  r = rbindlist(r, idcol = 'NIM')
  r[, p := as.double(gsub('%', '', q))/100]
  r[, T := 1/(1-p)]
  return(copy(r))
}

detrend = function(nim, wrt = 1){
  nfo = model_info(nim)
  rescale = function(x, xi, g, k){
    xi * (1 + g * (exp(k * x) -1) / k)
  }
  r = resid(nim)
  r = params2data(nim, data = r)
  r[, c('XI', 'G', 'K') := list(XI[wrt], G[wrt], K[wrt]), by = ID]
  r[, rsc:= rescale(value, XI, exp(G), K)]
  res = dcast.data.table(r[, .(eval(parse(text = nfo$cvrt)), ID, rsc)], eval(parse(text = nfo$cvrt)) ~ ID, value.var = 'rsc')
  setnames(res, 'nfo', nfo$cvrt)
  data.frame(res)
  }


## TDD growth curve
