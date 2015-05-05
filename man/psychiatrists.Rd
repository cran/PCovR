\name{psychiatrists}
\alias{psychiatrists}
\docType{data}
\title{
Effect of psychiatric symptoms on toxicomania, schizophrenia, depression and anxiety disorder
}
\description{
The data contain the scores of 30 Belgian psychiatric patients on 23 psychiatric symptoms and 4 psychiatric disorders (toxicomania, schizophrenia, depression, and anxiety disorder). Each score is the sumscore of the binary symptom and disorder scores that were given by 15 different psychiatrists. The data can be used to examine the extent to which the degree of toxicomania, schizophrenia, depression and anxiety disorder, can be predicted by the 23 psychiatric symptoms.
}
\usage{data(psychiatrists)}
\format{
  The format is:
List of 2

\subsection{$ X:'data.frame':  30 obs. of  23 variables:}{
  \describe{
    \item{\code{desorganised_speech}}{}
    \item{\code{agitation}}{}
    \item{\code{hallucinations}}{}
    \item{\code{inappropriate}}{}
    \item{\code{desorientation}}{}
    \item{\code{depression}}{}
    \item{\code{fear}}{}
    \item{\code{suicide}}{}
    \item{\code{somatic_concern}}{}
    \item{\code{narcotics}}{}
    \item{\code{antisocial}}{}
    \item{\code{retardation}}{}
    \item{\code{social_isolation}}{}
    \item{\code{routine}}{}
    \item{\code{alcohol}}{}
    \item{\code{negativism}}{}
    \item{\code{denial}}{}
    \item{\code{grandeur}}{}
    \item{\code{suspicion}}{}
    \item{\code{intellectual_obstruction}}{}
    \item{\code{impulse_control}}{}
    \item{\code{social_leveling}}{}
    \item{\code{occupational_dysfunction}}{}
    }}

\subsection{$ Y:'data.frame':  30 obs. of  4 variables:}{
  \describe{
    \item{\code{toxicomania}}{}
    \item{\code{schizophrenia}}{}
    \item{\code{depression}}{}
    \item{\code{anxiety_disorder}}{}
  }}
  }
\references{
\cite{Van Mechelen, I., & De Boeck, P. (1990). Projection of a binary criterion into a model of hierarchical classes. Psychometrika , 55 (4), 677-694.}

\cite{Marlies Vervloet, Henk A. Kiers, Wim Van den Noortgate, Eva Ceulemans (2015). PCovR: An R Package for Principal Covariates Regression. Journal of Statistical Software, 65(8), 1-14. URL http://www.jstatsoft.org/v65/i08/.}
}
\examples{
data(psychiatrists)
str(psychiatrists)
}
\keyword{datasets}
