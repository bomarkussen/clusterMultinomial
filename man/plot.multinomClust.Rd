\name{plot.multinomClust}
\alias{plot.multinomClust}
\title{Plot of multinomial clusters}
\description{Some plots}
\usage{
\method{plot}{multinomClust}(x, what = "dendrogram", add.names = TRUE, number.clusters = NULL, variables = NULL,
pval.correction = "holm", significance.level = 0.05, \dots)
}
\arguments{
  \item{x}{Fit from \code{\link{cluster}}.}
  \item{what}{Which plots to make? Possibilities are \code{dendrogram}, \code{pvalue} and \code{histogram}.}
  \item{add.names}{Should observation names be added in dendrogram? Also possible to give a character vector of names.}
  \item{number.clusters}{How many clusters should be used in histograms? If \code{NULL}, then taken from \code{obj}.}
  \item{variables}{Which variables to display in histograms. Either as characters or as numeric, where the latter gives variable numbers. If \code{NULL} then all variables are used.}
  \item{pval.correction}{Which correction for multiple testing should be used? Same as in \code{\link{p.adjust}}.}
  \item{significance.level}{Significance level to be used for color coding.}
  \item{\dots}{Not used.}
}
\details{To be done}
\value{Either \code{NULL} for dendrogram, which give base R graphics as side effect. Otherwise a \code{ggplot2} object.}
\author{Bo Markussen (bomar@math.ku.dk)}
\examples{
data(SEScluster)
plot(SEScluster,"histogram",variables=1:5)
}
\keyword{ ~graphs }
\keyword{ ~methods }
