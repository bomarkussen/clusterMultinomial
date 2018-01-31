\name{SEScluster}
\alias{SEScluster}
\docType{data}
\title{Multinomial clustering of Standard Event System dataset}
\description{The clustering suggests 16 clusters}
\usage{data("SEScluster")}
\format{Object of class \code{multinomClust}.}
\examples{
\dontrun{
data(SES)
SEScluster <- cluster(SES[,-(1:5)],SES[,1:5],.~GL,SES$dpi)
}

# Plot
data(SEScluster)
plot(SEScluster)
}
\keyword{datasets}
