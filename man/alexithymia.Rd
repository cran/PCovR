\name{alexithymia}
\alias{alexithymia}
\docType{data}
\title{
Effect of alexithymia on depression and self-esteem
}
\description{
The data contain the scores of 122 Belgian psychology students on the 20-item Toronto Alexithymia Scale (TAS-20; Bagby, Parker, & Taylor, 1994), which measures the inability to recognize and verbalize emotions, the Center for Epidemiological Studies Depression Scale (CES-D; Radloff, 1977), and the Rosenberg Self-Esteem Scale (RSE; Rosenberg, 1989). These data can be used to examine the extent to which the degree of depressive symptomatology (measured by the total CES-D score), and the degree of self-esteem (measured by the total RSE-score), can be predicted by the separate items of the TAS-20. We investigate the individual items because Bankier, Aigner and Bach (2001) emphasize that alexithymia is a multidimensional construct and authors disagree about the number and nature of the dimensions.}
\usage{data(alexithymia)}
\format{
  
List of 2
\subsection{$ X:'data.frame':  122 obs. of  20 variables:}{
  \describe{
    \item{\code{confused}}{I am often confused about what emotion I am feeling}
    \item{\code{right words}}{It is difficult for me to find the right words for my feelings}
    \item{\code{sensations}}{I have physical sensations that even doctors don't understand}
    \item{\code{describe}}{I am able to describe my feelings easily}
    \item{\code{analyze problems}}{I prefer to analyze problems rather than just describe them}
    \item{\code{upset}}{When I am upset, I don't know if I am sad, frightened, or angry}
    \item{\code{puzzled}}{I am often puzzled by sensations in my body}
    \item{\code{let happen}}{I prefer to just let things happen rather than to understand why they turned out that way}
    \item{\code{let happen}}{I have feelings that I can't quite identify}
    \item{\code{essential}}{Being in touch with emotions is essential}
    \item{\code{feel about people}}{I find it hard to describe how I feel about people}
    \item{\code{describe more}}{People tell me to describe my feelings more}
    \item{\code{going on}}{I don't know what's going on inside me}
    \item{\code{why angry}}{I often don't know why I am angry}
    \item{\code{daily activities}}{I prefer talking to people about their daily activities rather than their feelings}
    \item{\code{entertainment}}{I prefer to watch "light" entertainment shows rather than psychological dramas}
    \item{\code{reveal feelings}}{It is difficult for me to reveal my innermost feelings, even to close friends}
    \item{\code{close}}{I can feel close  to someone, even in moments of silence}
    \item{\code{useful}}{I find examination of my feelings useful in solving personal problems}
    \item{\code{hidden meanings}}{Looking for hidden meanings in movies or plays distracts from their enjoyment}
  }}
  \subsection{$ Y:'data.frame':  122 obs. of  2 variables:}{
  \describe{
    \item{\code{CES-D}}{Degree of depressive symptomatology}
    \item{\code{RSE}}{Degree of self-esteem}
  }}
}
\references{
\cite{Bagby, R. M., Parker, J. D., & Taylor, G. J. (1994). The twenty-item Toronto Alexithymia Scale: Item selection and cross-validation of the factor structure. Journal of Psychosomatic Research , 38 (1), 23-32.}

\cite{Bankier, B., Aigner, M., & Bach, M. (2001). Alexithymia in DSM-IV Disorder: Comparative Evaluation of Somatoform Disorder, Panic Disorder, Obsessive-Compulsive Disorder, and Depression. Psychosomatics , 42 (3), 235-240.}

\cite{Radloff, R. S. (1977). The CES-D scale: A self-report depression scale for research in the general population. Applied psychological measurement , 1 (3), 385-401.}

\cite{Rosenberg, M. (1989). Society and the adolescent self-image. Middletown: Wesleyan University Press.}
}
\examples{
data(alexithymia)
str(alexithymia)
}
\keyword{datasets}
