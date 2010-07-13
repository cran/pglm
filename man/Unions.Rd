\name{Unions}
\docType{data}
\alias{Unions}
\title{Unions and Wages}
\description{
  a balanced panel of 545 individuals for 8 years  1980  to 1987

  \emph{number of observations} :  4360

  \emph{observation} :  Individuals
  
  \emph{country} :  United States

}
\usage{data(Unions)}
\format{A dataframe containing :
  \describe{
    \item{id}{the individual index,}
    \item{year}{the year,}
    \item{age}{the age,}
    \item{exper}{the experience, computed as age - 6 - schooling,}
    \item{union}{does the wage is set by collective bargaining,}
    \item{married}{is the individual married ?}
    \item{ethn}{one of \code{'black'}, \code{'hisp'} and \code{'other'},}
    \item{disability}{does the individual has health disability ?}
    \item{rural}{does the individual lives in a rural area ?}
    \item{region}{the region, one of \code{'NorthEast'},
      \code{'NothernCentral'}, \code{'South'} and \code{'other'},}
    \item{wage}{hourly wage in US dollars,}
    \item{sector}{one of \code{'agricultural'}, \code{'mining'},
      \code{'construction'}, \code{'trade'}, \code{'transportation'},
      \code{'finance'},  \code{'businessrepair'}, \code{'personalservice'},
      \code{'entertainment'}, \code{'manufacturing'},
      \code{'pro.rel.service'}, \code{'pub.admin'},}
    \item{occ}{one of \code{'proftech'}, \code{'manoffpro'}, \code{'sales'},
      \code{'clerical'}, \code{'craftfor'}, \code{'operative'},
      \code{'laborfarm'}, \code{'farmlabor'}, \code{'service'}.}
  }
}

\source{
National Longitudinal Survey (Youth Sample)
}
\references{
Vella F. and Verbeek M. (1998) Whose wages do unions raise ? A dynamic
  model of unionism and wage rate determination for young men,
  \emph{Journal of Applied Economics}, 13, pp.163-183.
}

\keyword{datasets}
