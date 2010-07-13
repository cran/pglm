lnl.binomial <- function(param, y, X, id, model, link, rn, gradient = FALSE, hessian = FALSE,
                         start.sigma = FALSE, opposite = FALSE, direction = rep(0, length(param)),
                         initial.value = NULL, steptol = 1E-01){
  opposite <- ifelse(opposite, -1, +1)
  if (link == 'probit'){
    F <- pnorm
    f <- dnorm
    e <- function(x) - x * f(x)
  }
  if (link == 'logit'){
    F <- function(x) exp(x) / (1 + exp(x))
    f <- function(x) F(x) * (1 - F(x))
    e <- function(x) (exp(x) - exp(3*x))/(1+exp(x))^4
  }
  mills <- function(x) f(x) / F(x)
  millsp <- function(x) e(x)/F(x) - (f(x)/F(x))^2

  K <- ncol(X)
  n <- length(unique(id))
  N <- length(y)
  q <- 2 * y - 1
  step <- 2
  
  repeat{
    step <- step / 2
    if (step <  steptol) break
    beta <- param[1:K] + step * direction[1:K]
    bX <- as.numeric(crossprod(t(X), beta))
    
    if (start.sigma){
      mu <- -  tapply(q * mills(q * bX), id, sum) / tapply(millsp(q * bX), id, sum)
      return(sqrt(2) * sd(mu))
    }

    if (model == "random"){
      sigma <- param[K + 1L] + step * direction[K + 1L]
      Pitr <- lapply(rn$nodes, function(x) F( q * (bX + sigma * x)))
      Pir <- lapply(Pitr, function(x) tapply(x, id, prod))
      Li <- Reduce("+", mapply("*", Pir, rn$weights, SIMPLIFY = FALSE)) / sqrt(pi)
    }
    if (model == "pooling") Li <- F(q * bX)
    lnl <- opposite * sum(log(Li))
    if (is.null(initial.value) || lnl <= initial.value) break
  }    
  if (gradient){
    if (model == "random"){
      gitr <- mapply(function(w, x, p)
                     q * w * cbind(rep(1, N), x) *
                     as.numeric(p[as.character(id)]) *
                     mills( q * (bX + sigma * x)),
                     rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE)
      g <- Reduce("+", gitr)
      gradi <- cbind(g[, 1] * X, g[, 2]) /
        as.numeric(Li[as.character(id)])/ sqrt(pi)
    }
    if (model == "pooling") gradi <- q * mills(q * bX) * X

    attr(lnl, 'gradient') <- opposite * gradi
  }
  if (hessian){
    if (model == "pooling")
      H <- crossprod(millsp(q * bX) * X, X)
    if (model == "random"){
      Hr <- mapply(
                   function(w, v, p){
                     p <- p / Li
                     P <- p[as.character(id)]
                     p <- as.numeric(p)
                     P <- as.numeric(P)
                     z <- q * (bX + sigma * v)
                     gi <- q * mills(z) * cbind(X, v)
                     gi <- apply(gi, 2, tapply, id, sum)
                     H1 <- crossprod(p * gi, gi)
                     H2 <- crossprod(P * millsp(z) * cbind(X, v), cbind(X, v))
                     (H1 + H2) * w
                   },
                   rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE)
      H <- Reduce("+", Hr) / sqrt(pi) - crossprod(apply(gradi, 2, tapply, id, sum))
    }                  
    attr(lnl, 'hessian') <-   opposite * H
  }
  if (step < steptol) lnl <- NULL
  else{
    attr(lnl, 'step') <- step
  }
  lnl
}  
