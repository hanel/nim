lpar = function(p, ...){
  loc = p$.xi() * p$REG$XI
  scale = loc * exp(p$REG$G)
  shape = (loc^0) * c(p$REG$K)
  list(loc = loc, scale = scale, shape = shape)
}

lpar00 = function(p, ...){
  loc = p$.xi() * p$REG$XI
  scale = loc * exp(p$REG$G)
  shape = (loc^0) * c(p$REG$K)
  list(loc = loc, scale = scale, shape = shape)
}

# estimate k
lpar_k = function(p, THETA, ...){
  loc = THETA$.xi() * THETA$REG$XI
  scale = loc * exp(THETA$REG$G)
  shape = p * loc^0
  list(loc = loc, scale = scale, shape = shape)
}

lpar_gk = function(p, THETA, ...){
  loc = THETA$.xi() * THETA$REG$XI
  scale = loc * exp(p[1])
  shape = p[2] * loc^0
  list(loc = loc, scale = scale, shape = shape)
}








# estimate mu0
lpar_mu0 = function(p, ...){
  loc = p[1] * THETA$REG$M
  scale = loc * exp(THETA$REG$G)
  shape = THETA$REG$K
  list(loc = loc, scale = scale, shape = shape)
}

# estimate regional
lpar_REG = function(p, xpar){
  loc = THETA$.mu(xpar) * (p[1] + p[2] * (xpar$t - xpar$t0))
  scale = loc * exp(p[3] + p[4] * (xpar$t - xpar$t0))
  shape = THETA$REG$K[1] * loc^0 #p[5] + p[6] * (xpar$t - xpar$t0) * loc^0
  list(loc = loc, scale = scale, shape = shape)
}

lpar_REG2 = function(p, xpar){
  loc = .mu * (p[1] + p[2] * (.xt))
  scale = loc * exp(p[3] + p[4] * (.xt))
  shape = p[5] + p[6] * (.xt) * loc^0
  list(loc = loc, scale = scale, shape = shape)
}


# give all
fpar00 <- function(p, ...) {
  loc = p$.mu() * p$REG$M
  scale = loc * exp(p$REG$G)
  shape = p$REG$K * (loc^0)
  list(loc = loc, scale = scale, shape = shape)
}



ffx = function(data, start, fpar, xpar, std.err = TRUE, w = 1, ...) {

  nll.gev <- function(par) {
    pmat <- fpar(par, xpar)
    loc <- pmat$loc
    scale <- pmat$scale
    shape <- pmat$shape
    if(any(scale <= 0)) return(1e+20)
    gumbel <- (abs(shape) < 1e-06)
    y <- (data - loc) / scale
    z <- 1 + shape * y
    #		if(any(abs(shape)> .3))return(1e+20)
    if(any(z <= 0, na.rm = TRUE)) return(1e+20)
    nll <- (1 + 1 / shape) * log(z) + z^(-1 / shape)
    nll[gumbel] <- y[gumbel] + exp(-y[gumbel])
    sum(w * (nll + log(scale)), na.rm = TRUE)
  }

  call <- match.call()
  opt <- nlminb(start, nll.gev, control=list(...))
  gev <- fpar(opt$par, xpar)
  #cat(opt$convergence,'\n',opt$counts,'\n')
  out <- list(estimate = opt$par, deviance = 2 * opt$obj, convergence = opt$con, counts = opt$iter, message = opt$message, loc = gev$loc, scale = gev$scale, shape = gev$shape)
  #out = opt$par
  structure(c(out, call = call), class = "evd")
}



cffx = function(data, start, fpar, xpar, std.err = TRUE, w = rep(1, nrow(data)), method = 'Nelder-Mead', lower = -Inf, upper = Inf, ...) {

  nll.gev <- function(par) {
    nll_gev(par, fpar, data, w)
  }

  call <- match.call()
  #opt <- nlminb(start, nll.gev, control=list(...))
  opt <- optim(start, nll.gev, control=list(...), method = method, lower = lower, upper = upper)
  gev <- fpar(opt$par)
  #cat(opt$convergence,'\n',opt$counts,'\n')
  #out <- list(estimate = opt$par, deviance = 2 * opt$obj, convergence = opt$con, counts = opt$iter, message = opt$message, loc = gev$loc, scale = gev$scale, shape = gev$shape, eval = opt$eval)
  out <- list(estimate = opt$par, deviance = 2 * opt$value, convergence = opt$con, counts = opt$coun, message = opt$message, loc = gev$loc, scale = gev$scale, shape = gev$shape)
  #out = opt$par
  structure(c(out, call = call), class = "evd")
}
