\name{HealthIns}
\docType{data}
\alias{HealthIns}
\title{Health Insurance and Doctor Visits}
\description{
  an unbalanced panel of 5908 individuals for 1 to 5 years

  \emph{number of observations} :  20186

  \emph{observation} :  Individuals
  
  \emph{country} :  United States

}
\usage{data(HealthIns)}
\format{A dataframe containing :
  \describe{
    \item{id}{the individual index,}
    \item{year}{the year,}
    \item{med}{medical expenses,}
    \item{mdu}{number of face-to face medical visits,}
    \item{coins}{coinsurance rate,}
    \item{disease}{count of chronic diseases,}
    \item{sex}{a factor with level \code{'male'} and \code{'female'},}
    \item{age}{the age,}
    \item{size}{the size of the family,}
    \item{child}{a factor with levels \code{'no'} and \code{'yes'}.}
  }
}

\source{
  Rand Health Insurance Experiment
}
\references{
  
Manning W. G., Newhouse J. P.,  Duan N., Keeler E. B. and Leibowitz
  A. (1987) : \dQuote{Health Insurance and the Demand for Medical Care:
  Evidence from a Randomized Experiment}, \emph{American Economic
  Review}, 77(3), pp. 251-277.

Deb P. and P. K. Trivedi (2002) : \dQuote{The structure of demand for
  medical care: latent class versus two-part models}, \emph{Journal of
    Health Economics}, 21, pp. 601-625.
  
Cameron C. and P. K. Trivedi (2008) : \emph{Microeconometrics Using
  Stata}, Stata Press.
}

\keyword{datasets}
