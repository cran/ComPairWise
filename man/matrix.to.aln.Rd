\name{matrix.to.aln}
\alias{matrix.to.aln}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert a data matrix to an alignment}
\description{
This is the inverse of \code{\link{aln.to.matrix}}.  It reads a character matrix and turns it into an alignment object.
}
\usage{
matrix.to.aln(matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix}{A matrix of aligned data, usually but not necessarily DNA.  If taxon names are the row names of the matrix, they will be preserved.  If the matrix has no row names, the output taxon names will be NULL.}
}
\value{
  ~Describe the value returned
An alignment object, whose elements are
  \item{nb}{Number of aligned sequences}
  \item{nam}{Character vector of the names of the sequences (NULL if none included in the input matrix)}
  \item{seq}{Character vector of the sequences}
  \item{com}{currently always NA; could be a comment}
}
\author{TER}

\seealso{ \code{\link{aln.to.matrix}}}
\examples{
## The function is currently defined as
function (matrix) 
{
    x <- unname(apply(matrix, 1, paste, collapse = ""))
    aln <- list(nb = nrow(matrix), nam = rownames(matrix), seq = x, 
        com = NA)
    class(aln) <- "alignment"
    return(aln)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
