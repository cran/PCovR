\name{psychiatrists}
\alias{psychiatrists}
\docType{data}
\title{
Effect of psychiatric symptoms on toxicomania, schizophrenia, depression and anxiety disorder
}
\description{
The data contain the scores of 30 Belgian psychiatric patients on 23 psychiatric symptoms.The scores on the 23 psychiatric symptoms are the sums of the scores that were given to each patient by 15 different psychiatrists.
These data can be used to see which components underly the symptoms, and how these components can predict the scores on 4 psychiatric disorders: toxicomania, schizophrenia, depression and anexiety disorder. Also the scores on the disorders are the sums of scores that were given by the 15 psychiatrists.
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
\examples{
data(psychiatrists)
str(psychiatrists)
}
\keyword{datasets}
