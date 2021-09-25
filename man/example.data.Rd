\name{example.data}
\alias{example.data}
\docType{data}
\title{
Example data
}
\description{
Example data
}
\usage{data("example.data")}
\format{
 A list with 9 elements:
  \describe{
  	  \item{\code{w0.p}}{the baseline covariate in the control group in the prior study}
    \item{\code{s0.p}}{the surrogate marker in the control group in the prior study}
    \item{\code{y0.p}}{the primary outcome in the control group in the prior study}
  	 \item{\code{w1}}{the baseline covariate in the treatment group in the current study}
    \item{\code{w0}}{the baseline covariate in the control group in the current study}
     \item{\code{s1}}{the surrogate marker in the treatment group in the current study}
    \item{\code{s0}}{the surrogate marker in the control group in the current study}
      \item{\code{y1}}{the primary outcome in the treatment group in the current study}
    \item{\code{y0}}{the primary outcome in the control group in the current study}
}
}
\examples{
data(example.data)
names(example.data)
}

