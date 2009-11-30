\name{NameNeedle-package}
\Rdversion{1.1}
\alias{NameNeedle-package}
\alias{NameNeedle}
\docType{package}
\title{
  Implements Needleman-Wunsch simple global alignment algorithm
}
\description{
  The Needleman-Wunsch simple gap algorithm was one of the first
  methods introduced for global alignment of biological sequences.  The
  same algorithm can be used to match cell line names or sample names
  fromtwo related data sets; we provide examples in the documentaiton,
  using data that accompanies this package.
  
  While the \pkg{NameNeedle} package can be used for biological sequence
  alignment, the \pkg{Biostrings} package from Bioconductor contains
  much more sophisticated tools for that problem.
}
\details{
\tabular{ll}{
Package: \tab NameNeedle\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-11-18\cr
License: \tab GPL Version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\author{
Kevin R. Coombes

Maintainer: Kevin R. Coombes <kcoombes@mdanderson.org>
}
\keyword{ package }
\examples{
needles("HCC123", "HCC-123")
}
