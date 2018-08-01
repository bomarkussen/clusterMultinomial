\name{SES}
\alias{SES}
\docType{data}
\title{Standard Event System on 100 pig embryo}
\description{Collected by Kristian Kamstrup}
\usage{data("SES")}
\format{
  A data frame with 100 observations on the following 26 variables.
  \describe{
    \item{\code{Sow}}{a factor with sow id}
    \item{\code{Embryo}}{a factor with pig id}
    \item{\code{Carnegie.stage}}{Carnegie stage}
    \item{\code{dpi}}{days post insamination}
    \item{\code{GL}}{greatest length}
    \item{\samp{Neural tube}}{a factor with levels \code{Primitive streak} \code{Neural folds closure} \code{Anterior neuropore closed} \code{Neuropores closed}}
    \item{\code{Somites}}{a factor with levels \code{6-10 pairs} \code{11-15 pairs} \code{16-20 pairs} \code{21-25 pairs} \code{26-30 pairs} \code{31-35 pairs} \code{36-40 pairs} \code{41-45 pairs} \code{46-50 pairs} \code{>50 pairs}}
    \item{\code{Head}}{a factor with levels \code{Head bulbus} \code{Anterior cephalic projection} \code{Head projection disappeared}}
    \item{\code{Nose}}{a factor with levels \code{Olfactory pit} \code{External nares}}
    \item{\code{Ear}}{a factor with levels \code{Otic pit} \code{Otic vesicle} \code{Pinna fold}}
    \item{\code{Eye}}{a factor with levels \code{Optic vesicle} \code{Lens vesicle} \code{Optic fissure} \code{Contour lens/iris} \code{Pupil forms}}
    \item{\code{Ribs}}{a factor with levels \code{Rib primordia}}
    \item{\code{Heart}}{a factor with levels \code{Ventricle bulbus} \code{Thoracal bulbus disappeared} \code{Ventricle:liver 50:50} \code{Liver biggest}}
    \item{\code{Tail}}{a factor with levels \code{Tail bud} \code{Tail}}
    \item{\code{Forelimbs}}{a factor with levels \code{Ridge} \code{Bud} \code{Elongated} \code{Apical Ectodamal Ridge} \code{Elbow} \code{Paddle} \code{Digital plate}}
    \item{\code{Hindlimbs}}{a factor with levels \code{Apical Ectodamal Ridge} \code{Paddle} \code{Digital plate}}
    \item{\samp{Maxillary process}}{a factor with levels \code{Bud} \code{Posterior eye} \code{Midline eye} \code{Anterior lens} \code{Anterior eye} \code{Frontonasal fuse}}
    \item{\samp{Mandibular process}}{a factor with levels \code{Arch bud} \code{Posterior eye} \code{Posterior lens} \code{Midline eye} \code{Anterior lens} \code{Anterior eye} \code{Level frontonasal} \code{Occlusion point} \code{Insisur}}
    \item{\samp{Pharyngeal arches}}{a factor with levels \code{2nd arch} \code{3rd arch} \code{4th arch}}
    \item{\samp{Pharyngeal slits}}{a factor with levels \code{1st slit} \code{2nd slit} \code{3rd slit} \code{Slits closed}}
    \item{\code{Neck}}{a factor with levels \code{Cervical flexure 90?} \code{Cervical flexure disappeared} \code{45? angle}}
    \item{\samp{Eye lids}}{a factor with levels \code{Lower lid} \code{Begun overgrow} \code{Ventral lens} \code{Dorsal and ventral lens}}
    \item{\samp{Forefoot}}{a factor with levels \code{Digital grooves} \code{Digital serration} \code{Finger}}
    \item{\samp{CS}}{a factor with levels \code{Cervical sinus}}
    \item{\samp{Mammary line}}{a factor with levels \code{Posterior ear} \code{Mid ear} \code{Anterior ear} \code{1/3 from ear to eye} \code{2/3 from ear to eye} \code{Posterior eye} \code{Mid eye}}
    \item{\samp{Hindfoot}}{a factor with levels \code{Digital serration} \code{Finger}}
  }
}
\references{Kristians PhD Thesis}
\examples{
data(SES)
plot(GL~dpi,col=Sow,data=SES)
}
\keyword{datasets}
