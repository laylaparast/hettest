\name{delta.e.estimate}
\alias{delta.e.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests for a treatment effect on the primary outcome using surrogate marker information, ignoring potential heterogeneity}
\description{
Nonparametric test for a treatment effect on the primary outcome using surrogate marker information, ignoring potential heterogeneity. This test borrows information from a prior study about the relationship between the surrogate and the primary outcome to test for a treatment effect in the current study.}
\usage{
delta.e.estimate(sone = NULL, szero = NULL, szerop, yzerop, extrapolate = TRUE, 
mat = NULL, n1 = NULL, n0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sone}{
surrogate marker in the treated group in the current study}
  \item{szero}{
surrogate marker in the control group in the current study
}
  \item{szerop}{
surrogate marker in the control group in the prior study
}
  \item{yzerop}{
primary outcome in the control group in the prior study
}
  \item{extrapolate}{
TRUE or FALSE; extrapolate for values outside of the support in the prior study
}
  \item{mat}{
for the current study, the user can either provide sone and szero or can provide a vector, mat, where the first n1 values are the surrogate marker in the treated group in the current study, and the remaining values are the surrogate marker in the control group in the current study
}
  \item{n1}{
sample size of treated group in the current study; only needed if mat is provided instead of sone and szero
}
  \item{n0}{
sample size of control group in the current study; only needed if mat is provided instead of sone and szero
}
}
\value{
\item{delta.e}{estimated treatment effect using surrogate marker information}
\item{sd.e}{estimated standard error of treatment effect estimate}
\item{test.statistic.e}{test statistic for treatment effect}
\item{p.value.e}{p-value for test statistic}
}
\references{
Parast, Cai, and Tian (2021+). Using a Surrogate with Heterogeneous Utility to Test for a Treatment Effect. 
}
\author{
Layla Parast
}

\examples{
data(example.data)
delta.e.estimate(sone = example.data$s1, szero = example.data$s0, szerop = example.data$s0.p, 
yzerop = example.data$y0.p)
}
\keyword{nonparametric}
\keyword{test}