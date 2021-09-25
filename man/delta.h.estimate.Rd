\name{delta.h.estimate}
\alias{delta.h.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests for a treatment effect on the primary outcome using surrogate marker information, accounting for heterogeneity}
\description{
Nonparametric test for a treatment effect on the primary outcome using surrogate marker information, accounting for heterogeneity in the utility of the surrogate. This test borrows information from a prior study about the relationship between the surrogate and the primary outcome and the baseline covariate to test for a treatment effect in the current study.}
\usage{
delta.h.estimate(sone = NULL, szero = NULL, wone = NULL, wzero = NULL, szerop, 
wzerop, yzerop, extrapolate = TRUE, mat = NULL, n1 = NULL, n0 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{sone}{
surrogate marker in the treated group in the current study}
  \item{szero}{
surrogate marker in the control group in the current study
}
 \item{wone}{
baseline covariate in the treated group in the current study}
  \item{wzero}{
baseline covariate in the control group in the current study
}
  \item{szerop}{
surrogate marker in the control group in the prior study
}
  \item{wzerop}{
baseline covariate in the control group in the prior study
}
  \item{yzerop}{
primary outcome in the control group in the prior study
}
  \item{extrapolate}{
TRUE or FALSE; extrapolate for values outside of the support in the prior study
}
  \item{mat}{
for the current study, the user can either provide sone, szero, wone, wzero or can provide a vector, mat, where the first n1 values are the surrogate marker in the treated group in the current study, the second n0 values are the surrogate marker in the control group in the current study, the next n1 values are the baseline covariate in the treated group in the current study, the next n0 values are the baseline covariate in the control group in the current study
}
  \item{n1}{
sample size of treated group in the current study; only needed if mat is provided instead of sone, szero, wone, wzero
}
  \item{n0}{
sample size of control group in the current study; only needed if mat is provided instead of sone, szero, wone, wzero
}
}
\value{
\item{delta.h}{estimated treatment effect using surrogate marker information, account for heterogeneity}
\item{sd.h}{estimated standard error of treatment effect estimate}
\item{test.statistic.h}{test statistic for treatment effect}
\item{p.value.h}{p-value for test statistic}
}
\references{
Parast, Cai, and Tian (2021+). Using a Surrogate with Heterogeneous Utility to Test for a Treatment Effect. 
}
\author{
Layla Parast
}

\examples{
data(example.data)
delta.h.estimate(sone = example.data$s1, szero = example.data$s0, wone = example.data$w1,
 wzero = example.data$w0, szerop = example.data$s0.p, wzerop = example.data$w0.p, 
 yzerop = example.data$y0.p)
}
\keyword{nonparametric}
\keyword{test}