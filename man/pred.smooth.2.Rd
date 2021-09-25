\name{pred.smooth.2}
\alias{pred.smooth.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provides conditional mean estimate
}
\description{
Nonparametric estimate of the conditional mean: Y given X, applied to a new value X*}
\usage{
pred.smooth.2(kernel.use, kernel.apply, bw, outcome)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernel.use}{
values to condition on, X 
}
  \item{kernel.apply}{
new value for which the predicted outcome is estimated, X*
}
  \item{bw}{
bandwidth
}
  \item{outcome}{
outcome, Y
}
}
\value{
Predicted outcome for X* based on conditional mean estimate
}

\author{
Layla Parast
}
\keyword{internal}

