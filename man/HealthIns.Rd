\name{HealthIns}
\docType{data}
\alias{HealthIns}
\title{Health Insurance and Doctor Visits }
\description{
   observations of 401 individuals

 \emph{number of observations} :  20186

   \emph{country} :  United States

\emph{economic topic} : Health Economics

\emph{econometrics topic} : censored dependant variable

}
\usage{data(HealthIns)}
\format{A time serie containing :
  \describe{
    \item{id}{the individual index}
    \item{year}{the year}
    \item{med}{medical expenses}
    \item{mdu}{number of face-to face medical visits}
    \item{coins}{coinsurance rate}
    \item{disease}{count of chronic diseases}
    \item{sex}{a factor with level \code{'male'} and \code{'female'}}
    \item{age}{the age}
    \item{size}{the size of the family}
    \item{child}{a factor with levels \code{'no'} and \code{'yes'}}
  }
}
\source{
    Manning,  W. G.,  J. P.   Newhouse,   N.   Duan,  E. B.    Keeler and  A.  Leibowitz (1987) \dQuote{Health Insurance and the Demand for Medical Care: Evidence from a Randomized Experiment}, \emph{American Economic Review}, \bold{77(3)}, 251-277.

    Deeb P. , and  P.K.  Trivedi (2002) \dQuote{The structure of demand for medical care: latent class versus two-part models}, \emph{Journal of Health Economics}, \bold{21}, 601-625..

}
\references{
      \url{http://cameron.econ.ucdavis.edu/musbook/mus.html}.
}

\keyword{datasets}
