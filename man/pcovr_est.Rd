\name{pcovr_est}
\alias{pcovr_est}
\title{Estimation of Principal Covariates Regression parameters, given a prespecified weighting value and number of components}
\description{Analyzing regression data with many and/or highly collinear predictor variables, by simultaneously reducing the predictor variables to a given number of components and regressing the criterion variables on these components. A weighting parameter value is specified that determines the extent to which both aspects influence the solution. Cross-validation (Hastie, Tibshirani & Friedman, 2001) options are included.}
\usage{
pcovr_est(X, Y, r, a, cross = FALSE, fold = "LeaveOneOut")
}
\arguments{
  \item{X}{Matrix containing predictor scores (observations x predictors)}
  \item{Y}{Matrix containing criterion scores (observations x criteria)}
  \item{r}{The desired number of components}
  \item{a}{The desired weighting parameter value}
  \item{cross}{Logical. If \kbd{TRUE} cross-validation is performed}
  \item{fold}{Value of \emph{k} when performing \emph{k}-fold cross-validation. By default, leave-one-out cross-validation is performed.}
}
\value{
  \item{W}{Component weights matrix (predictors x components)}
  \item{B}{Regression weights for predictors (predictors x criteria)}
  \item{Rx2}{Proportion of explained variance in \env{X}}
  \item{Ry2}{Proportion of explained variance in \env{Y}}
  \item{Te}{Component score matrix (observations x components)}
  \item{Px}{Loading matrix of components (components x predictors)}
  \item{Py}{Regression weights matrix (components x criteria)}
  \item{Qy2}{Cross-validation fit}
}
\references{\cite{De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: Part I. Theory. Chemometrics and Intelligent Laboratory Systems , 155-164.}

\cite{Hastie, T., Tibshirani, R., & Friedman, J. (2001). The elements of statistical learning: Data mining, inference and prediction. New York: Springer.}}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\examples{
# standardize data
data(alexithymia)
Z <- nrm2(scale(alexithymia, scale = FALSE))*122^(1/2)
X <- Z[,1:20]
Y <- array(Z[,21:22],c(122,2))

# pcovr-analysis
results <- pcovr_est(X, Y, r=2, a=.90)
str(results)
}
\keyword{multivariate}
\keyword{regression}
