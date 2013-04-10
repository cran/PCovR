\name{PCovR-package}
\alias{PCovR-package}
\alias{PCovR}
\docType{package}
\title{Principal Covariates Regression}
\description{Analyzing regression data with many and/or highly collinear predictor variables, by simultaneously reducing the predictor variables to a limited number of components and regressing the criterion variables on these components. Several rotation options are provided in this package, as well as model selection options.}
\details{
\tabular{ll}{
Package: \tab \pkg{PCovR}\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-04-10\cr
License: \tab GPL (>= 2)\cr
}
This package contains the function \code{\link{pcovr}}, which runs a full PCovR analysis of a data set and provides several preprocessing, model selection, and rotation options. This function calls the function \code{\link{pcovr_est}}, which estimates the PCovR parameters given a specific weigthing parameter value and a particular number of components. This function was originally written in MATLAB by De Jong & Kiers (1992). Moreover, an illustrative data set is included, \code{\link{alexithymia}}.
}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\references{\cite{S. de Jong, H.A.L. Kiers, Principal covariates regression: Part I. Theory, Chemom. intell. lab. syst 14 (1992) 155-164.}}
\keyword{ package }
\seealso{
\code{\link{pcovr}}

\code{\link{pcovr_est}}

\code{\link{alexithymia}}

\code{\link{SortLoadings}}
}
\examples{
data(alexithymia)
results <- pcovr(alexithymia$X, alexithymia$Y)
summary(results)
plot(results)
print(t(SortLoadings(results$Px)))
}