\name{promin}
\alias{promin}
\title{Promin rotation}
\description{
This is a rotation criterion, developed by Lorenzo-Seva (1999), in which oblique target rotation (\code{\link{tarrotob}}) is applied using the Weighted Varimax solution (\code{\link{wvarim}}) as the target matrix. 
}
\usage{
promin(F1, nrs = 20)
}
\arguments{
  \item{F1}{Matrix to be rotated}
  \item{nrs}{Number of random starts}
}
\value{
\item{Th }{Transformation matrix to the pattern}
\item{loadings }{Rotated matrix}
\item{U}{Transformation matrix to the structure}
}
\references{
\cite{Lorenzo-Seva, U. (1999). Promin: A method for oblique factor rotation. Multivariate Behavioral Research, 34(3), 347-365.}
}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\examples{
Px <- matrix(rnorm(36),12,3)
print(Px)

Px_r <- promin(Px)
print(Px_r$loadings)
}
\keyword{multivariate}
\keyword{regression}
