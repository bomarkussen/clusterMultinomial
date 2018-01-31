\name{cluster}
\alias{cluster}
\title{Clustering of multinomial data}
\description{Agglomerative clustering of observations from a collection of multinomial variables.
}
\usage{
cluster(Y_data,X_data=NULL,logit.model=NULL,ordering=NULL,method="pvalue",verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_data}{data.frame with responses, which all should be factors}
  \item{X_data}{data.frame with possible predictors}
  \item{logit.model}{formula for logistic regression separating the observations}
  \item{ordering}{numeric used for cluster ordering}
  \item{method}{If \code{pvalue}, then p-values are used for selection. Otherwise the likelihood ratio is used.}
  \item{verbose}{output of algorithm tracking}
}
\details{To be added}
\value{
Object of class \code{multinomClust}, which is a LIST with 19 elements. Among these
  \item{Y_data }{Data frame with response data}
  \item{X_data }{Data frame with covariate data}
  \item{logit.model }{Model formula for logistic model}
  \item{category.names }{Names of variables in \code{Y_data}}
  \item{clusters.sequence }{List of list of clusters }
  \item{p.list }{Vector with p-values for agglomerative clustering}
  \item{number.clusters }{Integer with estimated number of clusters}
}
\author{Bo Markussen (bomar@math.ku.dk)}
\examples{
\dontrun{
data(SES)
SEScluster <- cluster(SES[,-(1:5)],SES[,1:5],.~GL,SES$dpi)
}

# Plot
data(SEScluster)
plot(SEScluster)
}
\keyword{ ~models }
\keyword{ ~multivariate }
\keyword{ ~cluster }
