\name{predict.multinomClust}
\alias{predict.multinomClust}
\title{Prediction of probabilities for cluster membership}
\description{Either within data used for clustering, or on new data.}
\usage{
\method{predict}{multinomClust}(object, number.clusters = NULL,
new_Y_data = NULL, new_X_data = NULL, \dots)
}
\arguments{
  \item{object}{Object of class \code{multinomClust} as fitted by \code{\link{cluster}}.}
  \item{number.clusters}{Number of used cluster. If \code{NULL}, then number of clusters taken from \code{obj}.}
  \item{new_Y_data}{Possible data frame with new response data.}
  \item{new_X_data}{Possible data frame with new predictor data.}
  \item{\dots}{Not used}
}
\details{To de done.}
\value{Output given as a matrix. Details to be provided later.}
\author{Bo Markussen (bomar@math.ku.dk)}
\examples{
\dontrun{
data(SEScluster)
predict(SEScluster)
}
}
\keyword{ ~cluster }
