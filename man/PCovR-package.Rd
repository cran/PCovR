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
Version: \tab 2.6\cr
Date: \tab 2015-05-05\cr
License: \tab GPL (>= 2)\cr
}
This package contains the function \code{\link{pcovr}}, which runs a full PCovR analysis of a data set and provides several preprocessing, model selection, and rotation options. This function calls the function \code{\link{pcovr_est}}, which estimates the PCovR parameters given a specific weigthing parameter value and a particular number of components. This function was originally written in MATLAB by De Jong & Kiers (1992). Two illustrative data sets are included: \code{\link{alexithymia}} and \code{\link{psychiatrists}}.
}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\references{\cite{S. de Jong, H.A.L. Kiers, Principal covariates regression: Part I. Theory, Chemom. intell. lab. syst 14 (1992) 155-164.}

\cite{Marlies Vervloet, Henk A. Kiers, Wim Van den Noortgate, Eva Ceulemans (2015). PCovR: An R Package for Principal Covariates Regression. Journal of Statistical Software, 65(8), 1-14. URL http://www.jstatsoft.org/v65/i08/.}}
\keyword{ package }
\seealso{
\code{\link{pcovr}}

\code{\link{pcovr_est}}

\code{\link{alexithymia}}

\code{\link{psychiatrists}}
}
\examples{
data(alexithymia)
results <- pcovr(alexithymia$X, alexithymia$Y)
summary(results)
plot(results)
}