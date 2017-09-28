\name{PatentsRDUS}
\docType{data}
\alias{PatentsRDUS}
\title{Dynamic Relation Between Patents and R\&d }
\description{
  yearly observations of 346 production units

 \emph{number of observations} :  3460

   \emph{country} :  United States

\emph{economic topic} : industrial economics

\emph{econometrics topic} : count data

}
\usage{data(PatentsRDUS)}
\format{A dataframe containing :
  \describe{
    \item{cusip}{compustat's identifying number for the firm}
    \item{year}{year}
    \item{ardssic}{a two-digit code for the applied R&D industrial classification}
    \item{scisect}{is the firm in the scientific sector ?}
    \item{capital72}{book value of capital in 1972}
    \item{sumpat}{the sum of patents applied for between 1972-1979}
    \item{rd}{ R&D spending during the year (in 1972 dollars)}
    \item{patents}{the number of patents applied for during the year that were eventually granted}
  }
}
\source{
    Hall, Browyn,  Zvi   Griliches and  Jerry  Hausman (1986) \dQuote{Patents and R and D: Is there a Lag?}, \emph{International Economic Review}, \bold{27}, 265-283.
}
\references{
      \url{http://cameron.econ.ucdavis.edu/racd/racddata.html}, chapter 9..
}

\keyword{datasets}
