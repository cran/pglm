\name{Fairness}
\docType{data}
\alias{Fairness}
\title{Perceived Fairness of Different Rules for Allocating Seats in
  Trains and Parking Spaces}
\description{
  401 individuals who answer 14 questions

  \emph{number of observations} :  5614
  
  \emph{observation} :  Individuals
  
  \emph{country} :  France

}
\usage{data(Fairness)}
\format{A dataframe containing :
  \describe{
    \item{id}{the individual index,}
    \item{answer}{a factor with levels 0 (very unfair), 1 (essentially
      unfair), 2 (essentially fair) and 3 (very fair),}
    \item{good}{one of \code{'tgv'} (French fast train) and \code{'Parking'},}
    \item{rule}{the allocation rule, a factor with levels \code{'peak'},
    \code{'admin'}, \code{'lottery'}, \code{'addsupply'},
    \code{'queuing'}, \code{'moral'} and \code{'compensation'},}
    \item{driving}{does the individual has the driving license ?}
    \item{education}{does the individual has a diploma ?}
    \item{recurring}{does the allocation problem is reccuring ?}
  }
}
\source{
Raux C., Souche S. and Croissant Y., 2009, How fair is pricing to be? An
empirical study, \emph{Public Choice}, vol. 139, 1, pp.227-240.
}
\references{
Raux C., Souche S. and Croissant Y., 2009, How fair is pricing to be? An
empirical study, \emph{Public Choice}, vol. 139, 1, pp.227-240.
}
\keyword{datasets}
