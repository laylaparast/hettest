\name{mu1swy.h}
\alias{mu1swy.h}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates conditional mean given S and W
}
\description{
Nonparametric estimation of the conditional expectation of the outcome given S and W
}
\usage{
mu1swy.h(s.use, w.use, y.use, s.apply, w.grd, h.paper3, h.paper2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s.use}{
surrogate marker that is being used to estimate conditional mean (prior study)
}
  \item{w.use}{
baseline covariate that is being used to estimate conditional mean (prior study)
}
  \item{y.use}{
outcome that is being used to estimate conditional mean (prior study)
}
  \item{s.apply}{
surrogate marker for which the conditional mean is being applied to (current study)
}
  \item{w.grd}{
baseline covariate for which the conditional mean is being applied to (current study)
}
  \item{h.paper3}{
bandwidth
}
  \item{h.paper2}{
bandwidth}
}

\value{
estimated conditional expectation of the outcome given S and W
}

\author{
Layla Parast
}
\keyword{internal}