Mort2Dsmooth_poisson = function (x, y, Z, offset, W, overdispersion = FALSE, ndx = c(floor(length(x)/5), 
                                                                             floor(length(y)/5)), deg = c(3, 3), pord = c(2, 2), lambdas = NULL, 
                         df = NULL, method = 1, coefstart = NULL, control = list()) 
{
  if (missing(W)) {
    W <- matrix(1, length(x), length(y))
  }
  if (missing(x)) {
    x <- 1:nrow(Z)
  }
  if (missing(y)) {
    y <- 1:ncol(Z)
  }
  m <- length(x)
  n <- length(y)
  if (missing(offset)) 
    offset <- matrix(0, m, n)
  if (length(offset) == 1) 
    offset <- matrix(offset, m, n)
  
  check <- Mort2Dsmooth_checker(x = x, y = y, Z = Z, offset = offset, 
                                W = W, overdispersion = overdispersion, ndx = ndx, deg = deg, 
                                pord = pord, lambdas = lambdas, df = df, method = method, 
                                coefstart = coefstart, control = control)
  x <- check$x
  y <- check$y
  Z <- check$Z
  m <- check$m
  n <- check$n
  offset <- check$offset
  offsetINIT <- check$offsetINIT
  wei <- check$W
  over <- check$overdispersion
  ndx <- check$ndx
  deg <- check$deg
  lambdas <- check$lambdas
  df <- check$df
  a.init <- check$coefstart
  MET <- check$method
  MON <- check$control$MON
  TOL1 <- check$control$TOL1
  TOL2 <- check$control$TOL2
  RANGEx <- check$control$RANGEx
  RANGEy <- check$control$RANGEy
  MAX.IT <- check$control$MAX.IT
  call <- match.call()
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  Bx <- MortSmooth_bbase(x, xmin, xmax, ndx[1], deg[1])
  nbx <- ncol(Bx)
  yl <- min(y)
  yr <- max(y)
  ymax <- yr + 0.01 * (yr - yl)
  ymin <- yl - 0.01 * (yr - yl)
  By <- MortSmooth_bbase(y, ymin, ymax, ndx[2], deg[2])
  nby <- ncol(By)
  Bx1 <- kronecker(matrix(1, ncol = nbx, nrow = 1), Bx)
  Bx2 <- kronecker(Bx, matrix(1, ncol = nbx, nrow = 1))
  RTBx <- Bx1 * Bx2
  By1 <- kronecker(matrix(1, ncol = nby, nrow = 1), By)
  By2 <- kronecker(By, matrix(1, ncol = nby, nrow = 1))
  RTBy <- By1 * By2
  Dx <- diff(diag(nbx), diff = pord[1])
  Dy <- diff(diag(nby), diff = pord[2])
  Px <- kronecker(diag(nby), t(Dx) %*% Dx)
  Py <- kronecker(t(Dy) %*% Dy, diag(nbx))
  # wei = W
  if (is.null(a.init)) {
    Z[is.na(Z)] <- 0
    xx <- rep(x, n)
    yy <- rep(y, each = m)
    fit0 <- glm(round(c(Z)) ~ xx + yy + offset(c(offset)), 
                family = poisson(link = "log"), weights = c(wei))
    etaGLM <- matrix(log(fit0$fitted) - c(offset), m, n)
    eta0 <- log((Z)) - offset
    eta0[wei == 0] <- etaGLM[wei == 0]
    BBx <- solve(t(Bx) %*% Bx + diag(nbx) * 1e-06, t(Bx))
    BBy <- solve(t(By) %*% By + diag(nby) * 1e-06, t(By))
    a.init <- MortSmooth_BcoefB(BBx, BBy, eta0)
  }
  else {
    a.init = a.init
  }
  psi2 <- 1
  if (MET == 1 | MET == 2) {
    if (over) {
      tol.over <- 10
      i.over <- 0
      while (tol.over > 0.001 && i.over < 5) {
        i.over <- i.over + 1
        lambdas.hat <- Mort2Dsmooth_optimize(x = x, y = y, 
                                             Z = Z, offset = offset, wei = wei, psi2 = psi2, 
                                             Bx = Bx, By = By, nbx = nbx, nby = nby, RTBx = RTBx, 
                                             RTBy = RTBy, Px = Px, Py = Py, a.init = a.init, 
                                             MON = MON, TOL1 = TOL1, TOL2 = TOL2, RANGEx = RANGEx, 
                                             RANGEy = RANGEy, MAX.IT = MAX.IT, MET = MET)
        FIT <- Mort2Dsmooth_estimate(x = x, y = y, Z = Z, 
                                     offset = offset, wei = wei, psi2 = psi2, Bx = Bx, 
                                     By = By, nbx = nbx, nby = nby, RTBx = RTBx, 
                                     RTBy = RTBy, lambdas = lambdas.hat, Px = Px, 
                                     Py = Py, a.init = a.init, MON = MON, TOL1 = TOL1, 
                                     MAX.IT = MAX.IT)
        psi2.old <- psi2
        psi2 <- FIT$dev/(m * n - FIT$df)
        tol.over <- abs(psi2 - psi2.old)/abs(psi2)
      }
    }
    else {
      lambdas.hat <- Mort2Dsmooth_optimize(x = x, y = y, 
                                           Z = Z, offset = offset, wei = wei, psi2 = psi2, 
                                           Bx = Bx, By = By, nbx = nbx, nby = nby, RTBx = RTBx, 
                                           RTBy = RTBy, Px = Px, Py = Py, a.init = a.init, 
                                           MON = MON, TOL1 = TOL1, TOL2 = TOL2, RANGEx = RANGEx, 
                                           RANGEy = RANGEy, MAX.IT = MAX.IT, MET = MET)
      FIT <- Mort2Dsmooth_estimate(x = x, y = y, Z = Z, 
                                   offset = offset, wei = wei, psi2 = psi2, Bx = Bx, 
                                   By = By, nbx = nbx, nby = nby, RTBx = RTBx, RTBy = RTBy, 
                                   lambdas = lambdas.hat, Px = Px, Py = Py, a.init = a.init, 
                                   MON = MON, TOL1 = TOL1, MAX.IT = MAX.IT)
      psi2 <- FIT$dev/(m * n - FIT$df)
    }
    if (log10(lambdas.hat[1]) >= log10(RANGEx[2]) | log10(lambdas.hat[1]) <= 
        log10(RANGEx[1])) {
      warning(paste("optimal lambda for x at the edge of its grid."))
    }
    if (log10(lambdas.hat[2]) >= log10(RANGEy[2]) | log10(lambdas.hat[2]) <= 
        log10(RANGEy[1])) {
      warning(paste("optimal lambda for y at the edge of its grid."))
    }
  }
  if (MET == 3) {
    lambdas.hat <- lambdas
    FIT <- Mort2Dsmooth_estimate(x = x, y = y, Z = Z, offset = offset, 
                                 wei = wei, psi2 = psi2, Bx = Bx, By = By, nbx = nbx, 
                                 nby = nby, RTBx = RTBx, RTBy = RTBy, lambdas = lambdas.hat, 
                                 Px = Px, Py = Py, a.init = a.init, MON = MON, TOL1 = TOL1, 
                                 MAX.IT = MAX.IT)
    psi2 <- FIT$dev/(m * n - FIT$df)
  }
  if (MET == 4) {
    Mort2Dsmooth_opt_df <- function(X) {
      FIT <- Mort2Dsmooth_estimate(x = x, y = y, Z = Z, 
                                   offset = offset, wei = wei, psi2 = psi2, Bx = Bx, 
                                   By = By, nbx = nbx, nby = nby, RTBx = RTBx, RTBy = RTBy, 
                                   lambdas = c(X[1], X[2]), Px = Px, Py = Py, a.init = a.init, 
                                   MON = MON, TOL1 = TOL1, MAX.IT = MAX.IT)
      return(abs(FIT$df - df))
    }
    by.lambda.x <- length(seq(RANGEx[1], RANGEx[2], by = TOL2))
    by.lambda.y <- length(seq(RANGEy[1], RANGEy[2], by = TOL2))
    by.lambda <- max(by.lambda.x, by.lambda.y)
    l.med.x <- median(log10(RANGEx))
    l.med.y <- median(log10(RANGEy))
    lambdas.hat <- cleversearch(fn = Mort2Dsmooth_opt_df, 
                                lower = c(RANGEx[1], RANGEy[1]), upper = c(RANGEx[2], 
                                                                           RANGEy[2]), startvalue = c(l.med.x, l.med.y), 
                                ngrid = by.lambda, logscale = TRUE, verbose = FALSE)[[1]]
    if (log10(lambdas.hat[1]) >= log10(RANGEx[2]) | log10(lambdas.hat[1]) <= 
        log10(RANGEx[1])) {
      warning(paste("optimal lambda for x at the edge of the grid."))
    }
    if (log10(lambdas.hat[2]) >= log10(RANGEy[2]) | log10(lambdas.hat[2]) <= 
        log10(RANGEy[1])) {
      warning(paste("optimal lambda for y at the edge of the grid."))
    }
    FIT <- Mort2Dsmooth_estimate(x = x, y = y, Z = Z, offset = offset, 
                                 wei = wei, psi2 = psi2, Bx = Bx, By = By, nbx = nbx, 
                                 nby = nby, RTBx = RTBx, RTBy = RTBy, lambdas = lambdas.hat, 
                                 Px = Px, Py = Py, a.init = a.init, MON = MON, TOL1 = TOL1, 
                                 MAX.IT = MAX.IT)
    psi2 <- FIT$dev/(m * n - FIT$df)
  }
  aic <- FIT$aic
  bic <- FIT$bic
  df <- FIT$df
  dev <- FIT$dev
  coef <- FIT$a
  psi2 <- psi2
  h <- FIT$h
  tolerance <- FIT$tol
  eta.hat <- matrix(MortSmooth_BcoefB(Bx, By, coef), m, n, 
                    dimnames = list(x, y))
  fitted.values <- matrix(exp(offset + eta.hat), m, n, dimnames = list(x, 
                                                                       y))
  lambdas.hat <- lambdas.hat
  res0 <- sign(Z - fitted.values)
  res1 <- log(ifelse(Z == 0, 1, Z/fitted.values))
  res2 <- Z - fitted.values
  res <- res0 * sqrt(2 * (Z * res1 - res2))
  res <- matrix(res, m, n, dimnames = list(x, y))
  wei <- matrix(wei, m, n, dimnames = list(x, y))
  object <- list(call = call, m = m, n = n, tolerance = tolerance, 
                 residuals = res, aic = aic, bic = bic, lev = h, df = df, 
                 dev = dev, psi2 = psi2, lambdas = lambdas.hat, ndx = ndx, 
                 deg = deg, pord = pord, x = x, y = y, Z = Z, offset = offsetINIT, 
                 W = wei, Bx = Bx, By = By, fitted.values = fitted.values, 
                 linear.predictors = eta.hat + offset, coefficients = coef, 
                 logmortality = eta.hat)
  class(object) <- "Mort2Dsmooth"
  object
}



Mort2Dsmooth_checker = function (x, y, Z, offset, W, overdispersion, ndx, deg, pord, 
                               lambdas, df, method, coefstart, control)
{
  if (missing(x)) {
    x <- 1:nrow(Z)
  }
  if (missing(y)) {
    y <- 1:ncol(Z)
  }
  offsetINIT <- offset
  m <- length(x)
  n <- length(y)
  whioff <- which(is.infinite(offset))
  whiwei <- which(W == 0)
  if (any(!whioff %in% whiwei)) {
    stop("weights different from zero associated with infinitive offset values")
  }
  offset[c(whioff, whiwei)] <- 100
  if (length(x) != nrow(Z)) 
    stop("length of x must be equal to number of rows in Z")
  if (length(y) != ncol(Z)) 
    stop("length of y must be equal to number of columns in Z")
  if (dim(Z)[1] != m | dim(Z)[2] != n) 
    stop("Argument arrays of wrong length")
  if (dim(offset)[1] != m | dim(offset)[2] != n) 
    stop("Argument arrays of wrong length")
  if (deg[1] < 1 | deg[1] >= 10 | deg[2] < 1 | deg[2] >= 10) 
    stop("Wrong values for deg")
  if (pord[1] <= 0 | pord[1] >= 5 | pord[2] <= 0 | pord[2] >= 
      5) 
    stop("Wrong value for pord")
  if (ndx[1] < 2 | ndx[1] >= floor(m * 0.9) | ndx[2] < 2 | 
      ndx[2] >= floor(n * 0.8)) 
    stop("Wrong value for ndx")
  coefstart.check <- is.null(coefstart)
  if (!coefstart.check) {
    if (nrow(coefstart) != (ndx[1] + deg[1]) | ncol(coefstart) != 
        (ndx[2] + deg[2])) {
      stop("coefstart must be a ndx[1]+deg[1] times ndx[2]+deg[2] matrix")
    }
  }
  if (method != 1 & method != 2 & method != 3 & method != 4) 
    stop("Wrong value for method")
  lambdas.check <- is.null(lambdas)
  df.check <- is.null(df)
  MET <- NULL
  if (lambdas.check & df.check & method == 1) {
    MET = 1
  }
  if (lambdas.check & df.check & method == 2) {
    MET = 2
  }
  if (lambdas.check & df.check & method == 3) {
    stop("with method 3, provide lambdas")
  }
  if (lambdas.check & df.check & method == 4) {
    stop("with method 4, provide df")
  }
  if (lambdas.check & !df.check & method == 1) {
    stop("df and method 1 cannot be chosen together")
  }
  if (lambdas.check & !df.check & method == 2) {
    stop("df and method 2 cannot be chosen together")
  }
  if (lambdas.check & !df.check & method == 3) {
    stop("df and method 3 cannot be chosen together")
  }
  if (lambdas.check & !df.check & method == 4) {
    MET = 4
    warning("Isotropic smoothing is applied", call. = FALSE)
  }
  if (!lambdas.check & df.check & method == 1) {
    stop("lambdas and method 1 cannot be chosen together")
  }
  if (!lambdas.check & df.check & method == 2) {
    stop("lambdas and method 2 cannot be chosen together")
  }
  if (!lambdas.check & df.check & method == 3) {
    MET = 3
  }
  if (!lambdas.check & df.check & method == 4) {
    stop("lambdas and method 4 cannot be chosen together")
  }
  if (!lambdas.check & !df.check) {
    stop("lambdas and df cannot be chosen together")
  }
  if (!lambdas.check && lambdas < 0) 
    stop("lambdas must be positive")
  if (!df.check && df < (pord[1] + pord[2])) 
    stop("df must be larger than the sum of pord")
  if (!df.check && df > ((ndx[1] + deg[1]) * (ndx[2] + deg[2]))) 
    stop("df must be smaller than the product of (ndx+deg) values")
  if (!df.check & length(df) != 1) 
    stop("df must be length 1")
  if (!lambdas.check & length(lambdas) != 1 & length(lambdas) != 
      2) 
    stop("lambda must be length 1 or 2")
  if (!lambdas.check & length(lambdas) == 1) {
    lambdas <- rep(lambdas, 2)
    warning("Isotropic smoothing is applied", call. = FALSE)
  }
  con <- list(MON = FALSE, TOL1 = 1e-06, TOL2 = 0.5, RANGEx = c(10^-4, 
                                                                10^6), RANGEy = c(10^-4, 10^6), MAX.IT = 50)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  if (nrow(W) != m | ncol(W) != n) {
    stop("dimensions of W and Z must be equal")
  }
  if (any(W == 0)) {
    warning("Interpolation and/or extrapolation is taking place", 
            call. = FALSE)
  }
  if (overdispersion & method == 3) 
    warning("given method 3, overdispersion is computed a posteriori")
  if (overdispersion & method == 4) 
    warning("given method 4, overdispersion is computed a posteriori")
  if (min(W) < 0) {
    warning(paste("At least one weight entry is negative"))
  }
  llist <- list(x = x, y = y, Z = Z, offset = offset, W = W, 
                m = m, n = n, offsetINIT = offsetINIT, overdispersion = overdispersion, 
                ndx = ndx, deg = deg, pord = pord, lambdas = lambdas, 
                df = df, method = method, coefstart = coefstart, control = con)
  llist
}

MortSmooth_bbase = function (x, xl, xr, ndx, deg) 
{
  dx <- (xr - xl)/ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, MortSmooth_tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  B
}

MortSmooth_BcoefB = function (X1, X2, mat) 
{
  BcoefB <- X1 %*% mat %*% t(X2)
  return(BcoefB)
}


Mort2Dsmooth_optimize = function (x, y, Z, offset, wei, psi2, Bx, By, nbx, nby, RTBx, 
                                  RTBy, Px, Py, a.init, MON, TOL1, TOL2, RANGEx, RANGEy, MAX.IT, 
                                  MET) 
{
  print("test")
  Mort2Dsmooth_opt_ic <- function(X) {
    FIT <- Mort2Dsmooth_estimate(x = x, y = y, Z = Z, offset = offset, 
                                 psi2 = psi2, wei = wei, Bx = Bx, By = By, nbx = nbx, 
                                 nby = nby, RTBx = RTBx, RTBy = RTBy, lambdas = c(X[1], 
                                                                                  X[2]), Px = Px, Py = Py, a.init = a.init, MON = MON, 
                                 TOL1 = TOL1, MAX.IT = MAX.IT)
    return(ifelse(MET == 2, FIT$aic, FIT$bic))
  }
  lambdas.x0 <- seq(log10(RANGEx[1]), log10(RANGEx[2]), TOL2 * 
                      4)
  lambdas.y0 <- seq(log10(RANGEx[1]), log10(RANGEx[2]), TOL2 * 
                      4)
  by.lambdas.x0 <- length(lambdas.x0)
  by.lambdas.y0 <- length(lambdas.y0)
  by.lambda.0 <- max(by.lambdas.x0, by.lambdas.y0)
  l.st.x0 <- median(lambdas.x0)
  l.st.y0 <- median(lambdas.y0)
  lambdas.hat0 <- cleversearch(fn = Mort2Dsmooth_opt_ic, lower = c(log10(RANGEx[1]), 
                                                                   log10(RANGEy[1])), upper = c(log10(RANGEx[2]), log10(RANGEy[2])), 
                               startvalue = c(l.st.x0, l.st.y0), ngrid = by.lambda.0, 
                               logscale = TRUE, verbose = FALSE)[[1]]
  l.st.x <- log10(lambdas.hat0[1])
  l.st.y <- log10(lambdas.hat0[2])
  min.l.x <- max(l.st.x - TOL2 * 4, log10(RANGEx[1]))
  max.l.x <- min(l.st.x + TOL2 * 4, log10(RANGEx[2]))
  min.l.y <- max(l.st.y - TOL2 * 4, log10(RANGEy[1]))
  max.l.y <- min(l.st.y + TOL2 * 4, log10(RANGEy[2]))
  lambdas.x <- seq(min.l.x, max.l.x, TOL2)
  lambdas.y <- seq(min.l.y, max.l.y, TOL2)
  by.lambdas.x <- length(lambdas.x)
  by.lambdas.y <- length(lambdas.y)
  by.lambda <- max(by.lambdas.x, by.lambdas.y)
  lambdas.hat <- cleversearch(fn = Mort2Dsmooth_opt_ic, lower = c(min.l.x, 
                                                                  min.l.y), upper = c(max.l.x, max.l.y), startvalue = c(l.st.x, 
                                                                                                                        l.st.y), ngrid = by.lambda, logscale = TRUE, verbose = FALSE)[[1]]
  return(lambdas.hat)
}


Mort2Dsmooth_estimate = function (x, y, Z, offset, psi2, wei, Bx, By, nbx, nby, RTBx, 
                                  RTBy, lambdas, Px, Py, a.init, MON, TOL1, MAX.IT) 
{
  P <- (lambdas[1] * Px) + (lambdas[2] * Py)
  tol <- 1
  i <- 0
  a <- a.init
  a.old <- 10
  if (MON) {
    cat("lambda.x =", lambdas[1], "\n")
    cat("lambda.y =", lambdas[2], "\n")
    cat("Iter         tol", "\n")
  }
  while (tol > TOL1 && i < MAX.IT) {
    i <- i + 1
    a <- Mort2Dsmooth_update(x = x, y = y, Z = Z, offset = offset, 
                             psi2 = psi2, wei = wei, Bx = Bx, By = By, nbx = nbx, 
                             nby = nby, RTBx = RTBx, RTBy = RTBy, P = P, a = a)
    tol <- max(abs(a - a.old)/abs(a))
    a.old <- a
    if (MON) {
      cat(i, "      ", tol, "\n")
    }
  }
  if (i > (MAX.IT - 1)) {
    warning(paste("parameter estimates did NOT converge in", 
                  MAX.IT, "iterations. Increase MAX.IT in control."))
  }
  eta <- MortSmooth_BcoefB(Bx, By, a)
  mu <- exp(offset + eta)
  W <- mu
  z <- eta + (1/mu) * (Z - mu)
  z[which(wei == 0)] <- 0
  WW <- wei * W
  BWB <- MortSmooth_BWB(RTBx, RTBy, nbx, nby, WW)
  BWBpP <- BWB + psi2 * P
  BWz <- MortSmooth_BcoefB(t(Bx), t(By), (WW * z))
  a0 <- solve(BWBpP, c(BWz))
  a <- matrix(a0, nrow = nbx)
  H <- solve(BWBpP, BWB)
  h <- diag(H)
  Z1 <- Z
  Z1[Z == 0] <- 10^(-4)

  dev <- 2 * (sum(wei * (Z1 * log(Z1/mu)) - (Z1 - mu), na.rm = TRUE))/psi2    # Dev = -2 (logLH of the model - logLH of the saturated model)
  df <- sum(h)


  # aic <- dev/psi2 + 2 * df  # I do not understand this definition, let's change it:
  logLH =  (-sum(mu) + sum(Z1 * log(mu)) - sum(lfactorial(Z1)))
  PenalizedlogLH = logLH - (1/2) * (t(a0) %*% P %*% (a0)) 

  aic = -2 * PenalizedlogLH/psi2 + 2 * df
  
  
  
  bic <- dev/psi2 + log(sum(wei)) * df  #
  
  llist <- list(a = a, h = h, df = df, aic = aic, bic = bic, 
                dev = dev, tol = tol, BWB = BWB, P = P)
  llist
}


Mort2Dsmooth_update = function (x, y, Z, offset, psi2, wei, Bx, By, nbx, nby, RTBx, 
                                RTBy, P, a) 
{
  eta <- MortSmooth_BcoefB(Bx, By, a)
  mu <- exp(offset + eta)
  W <- mu
  z <- eta + (1/mu) * (Z - mu)
  z[which(wei == 0)] <- 0
  WW <- wei * W
  BWB <- MortSmooth_BWB(RTBx, RTBy, nbx, nby, WW)
  BWz <- MortSmooth_BcoefB(t(Bx), t(By), (WW * z))
  a0 <- solve(BWB + psi2 * P, c(BWz))
  a <- matrix(a0, nrow = nbx)
  a
}

