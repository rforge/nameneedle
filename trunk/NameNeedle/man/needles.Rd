\name{needles}
\Rdversion{1.1}
\alias{needles}
\alias{needleScores}
\alias{defaultNeedleParams}
\title{
  Needleman-Wunsch global alignment algorithm
}
\description{
  This package contains a simple implementation of the Needleman-Wunsch
  global alignment algorithm.
}
\usage{
needles(pattern, subject, params = defaultNeedleParams)
needleScores(pattern, subjects, params = defaultNeedleParams)
defaultNeedleParams
}
\arguments{
  \item{pattern}{
    A character string (i.e., a character vector of length 1).
}
  \item{subject}{
    A character string (i.e., a character vector of length 1).
}
  \item{subjects}{
    A character vector.
}
  \item{params}{
    A list containing four required components. The default values are
    specified by the object \code{defaultNeedleParams}, which contains
    the following values:
\preformatted{   
 $ MATCH   : num 1
 $ MISMATCH: num -1
 $ GAP     : num -1
 $ GAPCHAR : chr "*"
}
}
}
\details{
  The Needleman-Wunsch global alignment algorithm was one of the first
  algorithms used to align DNA, RNA, or protein sequences.  The basic
  algorithm uses dynamic programming to find an optimal alignment
  between two sequences, with parameters that specify penalties for
  mismatches and gaps and a reward for exact matches.  More elaborate
  algorithms (not implemented here) make use of matrices with different
  penalties depending on different kinds of mismatches.  The version
  implemented here is based on the perl implementation in the first
  sectionn of Chapter 3 of the book BLAST, by Korf, Yandell, and Bedell. 
}
\value{
  The \command{needles} function returns a list with five components
\item{score}{The raw alignment score.}
\item{align1}{The final (optimal) alignment for the \code{pattern}.}
\item{align2}{The final (optimal) alignment for the \code{subject}.}
\item{sm}{The score matrix.}
\item{dm}{The backtrace matrix.}

  The \command{needleScores} function returns a numeric vector the same
  length as the \code{subjects} argument, with each entry equal to the
  corresponding raw alignment score.
}
\references{
Needleman SB, Wunsch CD: A general method applicable to the search for
similarities in the amino acid sequence of two proteins. J Mol Biol
1970, 48(3):443--453.

Korf I, Yandel M, Bedell J.  BLAST. O'Reilly Media, 2003.
}
\author{
Kevin R. Coombes <kcoombes@mdanderson.org>
}
\seealso{
  The \pkg{Biostrings} package from Bioconductor used to contain a
  function called \code{needwunQS} that provided a
  simple gap implementation of Needleman-Wunsch, similar to the one
  presented here. That function ha been deprecated in favor of a more
  elaborate interface called \code{\link[Biostrings]{pairwiseAlignment}}
  that incorporates a variety of other alignment methods in addition.
  While \code{pairwiseAlignment} is much more useful for applications to
  biological sequences, it is serious overkill for the application we
  have in mind for matching cell line or other sample names.
}
\examples{
data(cellLineNames)
myParam <- defaultNeedleParams
myParam$MATCH <- 2
myParam$MISMATCH <- -2
needles(sf2Names[2], illuNames[1], myParam)
scores <- needleScores(sf2Names[6], illuNames, myParam)
w <- which(scores==max(scores))
w
sf2Names[6]

needles(sf2Names[6], illuNames[w], myParam)
}
\keyword{ character }
