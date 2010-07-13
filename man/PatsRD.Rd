\name{PatsRD}
\docType{data}
\alias{PatsRD}
\title{Dynamic Relation Between Patents and R\&D }
\description{
  a panel of 346 observations  from 1970  to 1979 

 \emph{number of observations} :  3460

 \emph{observation} :  production units

   \emph{country} :  United States

}
\usage{data(PatsRD)}
\format{A dataframe containing :
  \describe{
    \item{cusip}{ Compustat's identifying number for the firm (Committee
      on Uniform Security Identification Procedures number)}
    \item{year}{year}
    
    \item{ardsic}{a two-digit code for the applied R&D industrial
    classification (roughly that in Bound, Cummins, Griliches, Hall, and
    Jaffe, in the Griliches R&D, Patents, and Productivityvolume)}
    \item{scisect}{is the firm in the scientific sector ?}
    \item{logk}{the logarithm of the book value of capital in 1972.}
    \item{sumpat}{the sum of patents applied for between 1972-1979. }
    \item{logr}{ the logarithm of R&D spending during the year (in 1972 dollars)}
    \item{pat}{the number of patents applied for during the year that
      were eventually granted}
  }
}
\source{
  Hall, Bronwyn , Zvi Griliches and Jerry Hausman (1986) \dQuote{Patents
  and R&D: Is There a Lag?}, \emph{International Economic Review},
  \bold{27}, 265-283.
}

\references{
  Cameron, A.C.  and Trivedi P.K. (1998) \emph{Regression analysis of
count data}, Cambridge University Press,
\url{http://cameron.econ.ucdavis.edu/racd/racddata.html}, chapter 9.

Cameron, A.C.  and P.K.  Trivedi (2005) \emph{Microeconometrics :
methods and applications}, Cambridge, pp. 792--5.
}
\keyword{datasets}
