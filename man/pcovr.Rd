\name{pcovr}
\alias{pcovr.default}
\alias{pcovr}
\alias{plot.pcovr}
\alias{print.pcovr}
\alias{summary.pcovr}
\title{Full Principal covariates regression analysis of a specific data set
}
\description{
Application of a PCovR analysis consists of the following steps: preprocessing the data, running PCovR analyses with different numbers of components and/or weighting parameter values, performing model selection, and rotating the retained solution for easier interpretation.
}
\usage{
pcovr(X,Y,modsel="seq",Rmin=1,Rmax=ncol(X)/3,R=NULL,weight=NULL,rot="varimax", 
target=NULL, prepX="stand",prepY="stand", ratio="estimation", fold="LeaveOneOut",
zeroloads=ncol(X))

\method{plot}{pcovr}(x,cpal=NULL,lpal=NULL,\dots)
}
\arguments{
  \item{X}{Dataframe containing predictor scores}
  \item{Y}{Dataframe containing criterion scores}
  \item{modsel}{Model selection procedure (\kbd{seq}, \kbd{seqRcv}, \kbd{seqAcv} or \kbd{sim})
}
  \item{Rmin}{Lowest number of components considered }
  \item{Rmax}{Highest number of components considered}
  \item{R}{Number of components (overrules \kbd{Rmin} and \kbd{Rmax})}
  \item{weight}{Weighting values considered}
  \item{rot}{Rotation criterion (\kbd{varimax}, \kbd{quartimin}, \kbd{targetT}, \kbd{targetQ}, \code{\link{wvarim}} or \code{\link{promin}})}
  \item{target}{Target matrix for target rotation (components x predictor variables)}
  \item{prepX}{Preprocessing of predictor scores: standardizing (\kbd{stand}) or centering data (\kbd{cent})}
  \item{prepY}{Preprocessing of criterion scores: standardizing (\kbd{stand}) or centering data (\kbd{cent})}
  \item{ratio}{Ratio of the estimated error variances of the predictor block and the criterion block}
  \item{fold}{Value of \emph{k} when performing \emph{k}-fold cross-validation. By default, leave-one-out cross-validation is performed.}
  \item{zeroloads}{Number of near-zero loadings of the target for \kbd{simplimax} rotation}
  \item{x}{An object of the type produced by \kbd{pcovr}}
  \item{cpal}{Vector of \code{\link{colors}} used for model selection plots}
  \item{lpal}{Vector of line types used for model selection plots}
  \item{\dots}{Further graphical arguments}
}
\details{
\subsection{Preprocessing}{
The PCovR package includes two preprocessing options, which can be applied to \kbd{X} and/or \kbd{Y}. Specifically, it is possible to only center the data (\kbd{prepX="cent"}, \kbd{prepY="cent"}). However, the default option is to standardize the data (\kbd{prepX="stand"}, \kbd{prepY="stand"}), which implies that \kbd{X} and/or \kbd{Y} are centered and normalized (i.e., each variable has a mean of zero and a standard deviation of one).}
\subsection{Model selection}{\subsection{Sequential procedure}{
The fastest and therefore default model selection setting (\kbd{modsel="seq"}) implies a sequential procedure in which the weighting value is determined on the basis of maximum likelihood principles (Vervloet, Van den Noortgate, Van Deun, & Ceulemans, 2013), but taking the  weighting values entered by the user (i.e., specified with the parameter \kbd{weight}) into account. Specifically, if the weighting value does not equal one of those values, the entered weighting value that is closest to the maximum likelihood weighting value (in absolute sense) is used. Note that the default error variance ratio is estimated with the function \code{\link{ErrorRatio}}, but can be specified otherwise with the parameter \kbd{ratio}. However, this is only possible for datasets with more observations than predictor variables. Among all models with the selected weighting value and a number of components between \kbd{Rmin} and \kbd{Rmax}, the solution is picked that has the highest \emph{st} value (Cattell, 1966; Wilderjans, Ceulemans, & Meers, 2012). However, models for which the fit is less than 1\% better than the fit of a less complex model are excluded. Note that the  assessment of the optimal number of components can be overruled, in case one is only interested in the solutions with a particular number of components. In particular, when specifying the input parameter \kbd{R}, \kbd{Rmin} and \kbd{Rmax} will be ignored, and the specified number of components will be used when running the analysis and determining the weighting value.

The package also provides two sequential procedures that incorporate a cross-validation step (\kbd{modsel="seqRcv"} and \kbd{modsel="seqAcv"}). \kbd{seqRcv} also starts with the selection of the weighting value based on maximum likelihood principles, but in the next step, the number of components is determined using leave-one-out cross-validation. \kbd{seqAcv} is identical to the default procedure, but has an extra step: after the selection of the number of components, leave-one-out cross-validation is applied to choose the weighting value.}
\subsection{Simultaneous procedure}{The simultaneous procedure (\kbd{modsel="sim"}) performs leave-one-out cross-validation for all considered weighting values (\kbd{weight}; by default, 100 values between .01 and 1) and all numbers of components between \kbd{Rmin} (default: 1) and \kbd{Rmax} (default: number of predictors divided by 3). The weighting parameter value and number of components that maximize the cross-validation fit are retained. Note that the parameter \kbd{fold} can be used to alter the number of roughly equal-sized parts in which the data are split for cross-validation (Hastie, Tibshirani, & Friedman, 2001). The default value of \kbd{fold} is \kbd{"LeaveOneOut"}, implying that the data is split in \emph{N} (number of observations) parts. }
}
\subsection{Interpreting the component matrices}{The rotation criteria that are implemented in the PCovR package are \kbd{varimax}, \kbd{quartimin}, \kbd{targetT}, \kbd{targetQ}, \code{\link{wvarim}} and \code{\link{promin}}. One can also request the original solution by typing \kbd{rot="none"}. Target rotation (Browne, 1972) orthogonally rotates the loading matrix towards a target matrix (\kbd{target}) that is specified by the user. 
Note that Simplimax requires the specification of a number of zero elements. By default, \kbd{zeroloads} equals the number of predictors.   

The interpretation of the obtained solution usually starts with the interpretation of the loading matrix. Specifically, the components are labeled by considering what the predictors that have the highest loadings (in absolute sense), have in common. Given these labels, the regression weights can be interpreted.}}
\value{
  \code{pcovr} returns a list that contains the following objects (note that some objects can be empty, depending on the model selection settings used) :
  \item{Px }{Loading matrix (components x predictor variables)}
  \item{Py }{Regression weights matrix (components x criterion variables)}
  \item{Te }{Component score matrix (observations x components)}
  \item{W }{Component weights matrix (predictor variables x components)}
  \item{Rx2 }{Proportion of explained variance in \env{X}}
  \item{Ry2 }{Proportion of explained variance in \env{Y}}
  \item{Qy2 }{Cross-validation fit as a function of weighting parameter and number of components (weighting parameter x number of components)}
  \item{VAFsum}{Weighted sum of the variance accounted for in \env{X} and in \env{Y} as a function of number of components (1 x number of components)}
  \item{alpha}{Selected value of the weighting parameter}
  \item{R}{Selected number of components}
  \item{modsel}{Model selection procedure that was used}
  \item{rot}{Rotation criterion that was used}
  \item{prepX}{Method of preprocessing that was used for the predictor scores}
  \item{prepY}{Method of preprocessing that was used for the criterion scores}
  \item{Rvalues}{Numbers of components that were considered}
  \item{Alphavalues}{Weighting parameter values that were considered}
}
\references{
\cite{Browne, M. W. (1972). Oblique rotation to a partially specified target. British Journal of Mathematical and Statistical Psychology , 25 (2), 207-212.}

\cite{Cattell, R. B. (1966). The scree test for the number of factors. Multivariate behavioral research , 1 (2), 245-276.}

\cite{De Jong, S., & Kiers, H. A. (1992). Principal covariates regression: Part I. Theory. Chemometrics and Intelligent Laboratory Systems , 155-164.}

\cite{Hastie, T., Tibshirani, R., & Friedman, J. (2001). The elements of statistical learning: Data mining, inference and prediction. New York: Springer.}

\cite{Vervloet, M., Van Deun, K., Van den Noortgate, W., & Ceulemans, E. (2013). On the selection of the weighting parameter value in Principal Covariates Regression. Chemometrics and Intelligent Laboratory Systems.}

\cite{Marlies Vervloet, Henk A. Kiers, Wim Van den Noortgate, Eva Ceulemans (2015). PCovR: An R Package for Principal Covariates Regression. Journal of Statistical Software, 65(8), 1-14. URL http://www.jstatsoft.org/v65/i08/.}

\cite{Wilderjans, T. F., Ceulemans, E., & Meers, K. (2012). CHull: A generic convex-hull-based model selection method. Behavior research methods .}}
\author{Marlies Vervloet (\email{marlies.vervloet@ppw.kuleuven.be})}
\examples{
data(alexithymia)
results <- pcovr(alexithymia$X, alexithymia$Y)
summary(results)
plot(results)
}
\keyword{regression}