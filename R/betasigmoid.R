
#' Logit and Sigmoid function
#'
#' @param x numerical vectors
#' @param infl position of the inflection point of the sigmoid
#' @param scale the width of the sigmoid
#'
#' @export
logit <- function(x, infl=0, scale=1, oneminus = FALSE){
  # infl + scale * log(x/(1-x))
  if(scale < 0){
    scale <- -scale
    oneminus <- ! oneminus
  }
  qlogis(x, location=infl, scale=scale, log.p = FALSE, lower.tail = ! oneminus)
}

#' @rdname logit
#' @export
sigmoid <- function(x, infl=0, scale=1, log=FALSE, oneminus = FALSE){
  # 1/(1+ exp(-(x-infl)/scale))
  if(scale < 0){
    scale <- -scale
    oneminus <- ! oneminus
  }
  plogis(x, location=infl, scale=scale, log.p =log, lower.tail = ! oneminus)
}


#' Beta-Sigmoid probability distribution
#'
#' The distribution is obtained as a transformation of the Beta distribution with a sigmoid function. It
#' is useful for modelling the distribution of a mean if all that is known is a number of success and
#' failures and the sigmoid distribution that assignes the success and failures:
#' \deqn{m \sim \mathrm{Bernoulli}(\mathrm{sigmoid(\mu)})}{b ~ Bernoulli(sigmoid(µ))}
#' where a is number of times that m is 1 and b is the number of times it is 0.
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param a,b non-negative parameters of the Beta-Sigmoid distribution which can be understood as number of success and failures.
#'    They correspond to the `shape1` and `shape2` parameters of the Beta distribution. Default: 1
#' @param infl position of the inflection point of the sigmoid. Default: 0
#' @param scale positive number that describes how broad the sigmoid is. Default: 1
#' @param mean,sd mean and standard deviation of the distribution. It is an alternative
#'    parametrization to \code{infl} and \code{scale} only use either of the two.
#' @param log,log.p logical; if `TRUE`, probabilities are given as log(p)
#' @param lower.tail logical; if `TRUE` (default), probabilities are `P[X < x]`, otherwise `P[X > x]`
#' \deqn{x \sim \mathrm{BetaSigmoid}(\alpha, \beta)}{x ~ BetaSigmoid(a, b, i, s)}
#' @details
#' The probability density function of the Beta-Sigmoid distribution is
#' \deqn{f_{BS_{\alpha, \beta, i, s}}(x) = \frac{\Gamma(\alpha + \beta)}{\Gamma(\alpha) \Gamma(\beta)}
#'        \Big(\frac{1}{1+e^{-x}}\Big)^\alpha \Big(1  - \frac{1}{1+e^{-(x - i) / s}}\Big)^\beta }
#' @examples
#' dbetasigmoid(4, a=3, b=1, infl=3.5, scale=3)
#' pbetasigmoid(-60, a=0.3, b=10, infl=3.5, scale=5)
#' qbetasigmoid(0.8, a=3, b=3, infl=0, scale=0.5)
#' rbetasigmoid(3, a=19, b=39)
#' @export
dbetasigmoid <- function(x,  a=1, b=1, infl=0, scale=1,
                         mean = infl + scale * (digamma(a) - digamma(b)),
                         sd =  scale * sqrt(trigamma(a) + trigamma(b)),
                         log=FALSE){

  # Alternative parameterization using mu and sd
  if (!missing(infl) && !missing(mean)) {
    if (abs((infl + scale * (digamma(a) - digamma(b))) - mu) < 1e-15)
      warning("specify 'infl' or 'mean' but not both")
    else stop("specify 'infl' or 'mean' but not both")
  }
  if (!missing(scale) && !missing(sd)) {
    if (abs(( scale * sqrt(trigamma(a) + trigamma(b))) - sd) < 1e-15)
      warning("specify 'scale' or 'sd' but not both")
    else stop("specify 'scale' or 'sd' but not both")
  }

  scale <- sd / sqrt((trigamma(a) + trigamma(b)))
  infl <- mean - scale * (digamma(a) - digamma(b))
  # End Alternative parameterization


  sigx <- sigmoid(x, infl=infl, scale=scale, log=TRUE)
  sig1mx <- sigmoid(x, infl=infl, scale=scale, log=TRUE, oneminus = TRUE)
  rval <- sigx * (a - 1) + sig1mx * (b - 1) - lbeta(a, b) +
    sigx + sig1mx - log(abs(scale))
  if(! log){
    exp(rval)
  }else{
   rval
  }
}

#' @rdname dbetasigmoid
#' @export
rbetasigmoid <- function(n, a=1, b=1, infl=0, scale=1,
                         mean = infl + scale * (digamma(a) - digamma(b)),
                         sd =  scale * sqrt(trigamma(a) + trigamma(b))){
  # Alternative parameterization using mu and sd
  if (!missing(infl) && !missing(mean)) {
    if (abs((infl + scale * (digamma(a) - digamma(b))) - mu) < 1e-15)
      warning("specify 'infl' or 'mean' but not both")
    else stop("specify 'infl' or 'mean' but not both")
  }
  if (!missing(scale) && !missing(sd)) {
    if (abs(( scale * sqrt(trigamma(a) + trigamma(b))) - sd) < 1e-15)
      warning("specify 'scale' or 'sd' but not both")
    else stop("specify 'scale' or 'sd' but not both")
  }

  scale <- sd / sqrt((trigamma(a) + trigamma(b)))
  infl <- mean - scale * (digamma(a) - digamma(b))
  # End Alternative parameterization
  logit(stats::rbeta(n, shape1=a, shape2=b)) * scale + infl
}

#' @rdname dbetasigmoid
#' @export
pbetasigmoid <- function(q, a=1, b=1, infl=0, scale=1,
                         mean = infl + scale * (digamma(a) - digamma(b)),
                         sd =  scale * sqrt(trigamma(a) + trigamma(b)),
                         lower.tail=TRUE, log.p=FALSE){
  # Alternative parameterization using mu and sd
  if (!missing(infl) && !missing(mean)) {
    if (abs((infl + scale * (digamma(a) - digamma(b))) - mu) < 1e-15)
      warning("specify 'infl' or 'mean' but not both")
    else stop("specify 'infl' or 'mean' but not both")
  }
  if (!missing(scale) && !missing(sd)) {
    if (abs(( scale * sqrt(trigamma(a) + trigamma(b))) - sd) < 1e-15)
      warning("specify 'scale' or 'sd' but not both")
    else stop("specify 'scale' or 'sd' but not both")
  }

  scale <- sd / sqrt((trigamma(a) + trigamma(b)))
  infl <- mean - scale * (digamma(a) - digamma(b))
  # End Alternative parameterization
  stats::pbeta(sigmoid(q, infl=infl, scale=scale), shape1=a, shape2=b, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname dbetasigmoid
#' @export
qbetasigmoid <- function(p, a=1, b=1, infl=0, scale=1,
                         mean = infl + scale * (digamma(a) - digamma(b)),
                         sd =  scale * sqrt(trigamma(a) + trigamma(b)),
                         lower.tail=TRUE, log.p=FALSE){
  # Alternative parameterization using mu and sd
  if (!missing(infl) && !missing(mean)) {
    if (abs((infl + scale * (digamma(a) - digamma(b))) - mu) < 1e-15)
      warning("specify 'infl' or 'mean' but not both")
    else stop("specify 'infl' or 'mean' but not both")
  }
  if (!missing(scale) && !missing(sd)) {
    if (abs(( scale * sqrt(trigamma(a) + trigamma(b))) - sd) < 1e-15)
      warning("specify 'scale' or 'sd' but not both")
    else stop("specify 'scale' or 'sd' but not both")
  }

  scale <- sd / sqrt((trigamma(a) + trigamma(b)))
  infl <- mean - scale * (digamma(a) - digamma(b))
  # End Alternative parameterization
  logit(stats::qbeta(p, shape1=a, shape2=b, lower.tail=lower.tail, log.p=log.p), infl=infl, scale=scale)
}

mode_betasigmoid <- function(a=1, b=1, infl=0, scale = 1){
  infl + scale * (log(a) - log(b))
}

mean_betasigmoid <- function(a=1, b=1, infl=0, scale = 1){
   infl + scale * (digamma(a) - digamma(b))
}

var_betasigmoid <- function(a=1, b=1, infl=0, scale = 1){
  scale^2 * (trigamma(a) + trigamma(b))
}


