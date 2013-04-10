\name{SortLoadings}
\alias{SortLoadings}
\title{Sorting a component loading matrix}
\description{
A loading matrix indicates how predictors that have been reduced to components - e.g., in principal covariates regression (De Jong & Kiers, 1992) - relate to these components. Usually, components are interpreted by looking at what the predictors with a clear non-zero loading have in common.
To make this easier, this function changes the order of the predictors presented in a loading matrix, so that the firstly, the predictors with clear non-zero loadings on the first component (with decreasing loadings) are presented, then the predictors with clear non-zero loadings on the second component, etc.
}
\usage{
SortLoadings(Px, cutoff = .20)
}
\arguments{
  \item{Px}{Dataframe that contains component loadings (components x predictors)}
  \item{cutoff}{Values higher than this will be considered as a clear non-zero loading}
}
\value{
\code{sortloads} returns a dataframe with the same dimensions and labels as the original loading matrix, but with the columns (referring to the predictors) presented in a different order.
}
\references{
\cite{De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: Part I. Theory. Chemometrics and Intelligent Laboratory Systems , 155-164.}
}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\seealso{
\code{\link{pcovr}}
}
\examples{
# Compute loading matrix of alexithymia dataset
data(alexithymia)
results <- pcovr(alexithymia$X,alexithymia$Y)
Px <- results$Px
print(Px)

# Sort loading matrix
sorted_Px <- SortLoadings(results$Px)
print(sorted_Px)
}
\keyword{multivariate}
\keyword{regression}
