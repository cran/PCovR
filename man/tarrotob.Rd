\name{tarrotob}
\alias{tarrotob}
\title{Oblique target rotation}
\description{
Oblique target rotation
}
\usage{
tarrotob(F1, W)
}
\arguments{
  \item{F1}{Matrix to be rotated}
  \item{W}{Target binary matrix}
}
\value{
\item{T }{Rotation matrix}
\item{A }{Rotated matrix}
}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\examples{
Px <- matrix(rnorm(36),12,3)
print(Px)

W <- matrix(rbinom(36,1,.4),12,3)
Px_r <- tarrotob(Px,W)
print(Px_r$A)
}
\keyword{multivariate}
\keyword{regression}
