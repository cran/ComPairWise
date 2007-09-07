\name{dots.to.bases}
\alias{dots.to.bases}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert alignment formats from dot-match to all bases}
\description{
dot.to.bases takes an alignment in which bases that match the first taxon are shown with a dot or other character, and replaces them with the base.
}
\usage{
dots.to.bases(alignment, matchchar = ".")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{An alignment object}
  \item{matchchar}{Character designating matched base}
}
\value{
A new alignment object, whose elements are
  \item{nb}{Number of aligned sequences}
  \item{nam}{Character vector of the names of the sequences}
  \item{seq}{Character vector of the sequences}
  \item{com}{currently always NA; could be a comment}
}
\author{TER}
\note{Designed as an internal function, and is called by read.nexus, but works as a standalone.  }
\seealso{\code{\link{read.nexus}}, \code{\link{read.phylip}}, \code{\link{aln.to.matrix}}}
\examples{

## The function is currently defined as
function (alignment, matchchar = ".") 
{
    dotsout <- function(aln.column) {
        aln.column[aln.column == matchchar] <- aln.column[1]
        return(aln.column)
    }
    aln.new <- aln.to.matrix(alignment, F)
    aln.new <- apply(aln.new, 2, dotsout)
    aln.new <- matrix.to.aln(aln.new)
    alignment$seq <- aln.new$seq
    return(alignment)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
