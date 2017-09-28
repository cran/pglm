\name{Fairness}
\docType{data}
\alias{Fairness}
\title{Perveived Fairness of Rules for Allocating Seats in Trains and Parking Spaces }
\description{
   observations of 401 individuals

 \emph{number of observations} :  5614

   \emph{country} :  France

\emph{economic topic} : public economics

\emph{econometrics topic} : ordered response

}
\usage{data(Fairness)}
\format{A dataframe containing :
  \describe{
    \item{id}{the individual index}
    \item{answer}{a factor with levels 0 (very unfair), 1 (essentially unfair), 2 (essentially fair) and 3 (very fair)}
    \item{good}{one of \code{'tgv'} (French fast train) and \code{'Parking'}}
    \item{rule}{the allocation rule, a factor with levels \code{'peak'}, \code{'admin'}, \code{'lottery'}, \code{'addsupply'}, \code{'queuing'}, \code{'moral'} and \code{'compensation'}}
    \item{driving}{does the individual has the driving license ?}
    \item{education}{does the individual has a diploma ?}
    \item{recurring}{does the allocation problem is reccuring ?}
  }
}
\source{
    provided by the authors.

}
\references{
    Charles Raux, Stephanie Souche and Yves Croissant (2009) \dQuote{How fair is pricing perceived to be? An empirical study}, \emph{Public Choice}, \bold{139(1)}, 227-240.
}
\keyword{datasets}
