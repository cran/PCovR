\name{ErrorRatio}
\alias{ErrorRatio}
\title{
Error variance ratio
}
\description{
Estimating the ratio of the error variance of the predictor block, versus the error variance of the criterion block.
}
\usage{
ErrorRatio(X, Y, Rmin = 1, Rmax = ncol(X)/3, prepX="stand",prepY="stand")
}
\arguments{
  \item{X}{Dataframe containing predictor scores}
  \item{Y}{Dataframe containing criterion scores}
  \item{Rmin}{Lowest number of components considered}
  \item{Rmax}{Highest number of components considered}
  \item{prepX}{Preprocessing of predictor scores: standardizing (\kbd{stand}) or centering data (\kbd{cent})}
  \item{prepY}{Preprocessing of criterion scores: standardizing (\kbd{stand}) or centering data (\kbd{cent})}
}
\details{
An estimate for the error variance of \kbd{X} can be obtained by applying principal component analysis to \kbd{X} and determining the optimal number of components through a scree test; the estimate equals the associated percentage of unexplained variance. The estimate for the error variance of \kbd{Y} boils down to the percentage of unexplained variance when \kbd{Y} is regressed on \kbd{X}.  This approach for estimating  and  was based on the work of Wilderjans, Ceulemans, Van Mechelen, and Van den Berg (2011).
}
\value{
The returned value is the estimated error variance of \kbd{X}, divided by the estimated error variance of \kbd{Y}.
}
\references{
\cite{Wilderjans, T. F., Ceulemans, E., Van Mechelen, I., & Van den Berg, R. A. (2011). Simultaneous analysis of coupled data matrices subject to different amounts of noise. British Journal of Mathematical and Statistical Psychology , 64, 277-290.}
}
\author{
Marlies Vervloet
}
\examples{
data(psychiatrists)
ratio <- ErrorRatio(psychiatrists$X,psychiatrists$Y)
}
\keyword{regression}