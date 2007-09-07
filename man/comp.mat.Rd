\name{comp.mat}
\alias{comp.mat}
\title{Pairwise comparison of alignment matrices}
\description{
Internal ComPairWise function.  Probably not useful as a standalone.
}
\usage{
comp.mat(align1, align2, cols1, cols2, ref, i, j, name1, name2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{align1}{First alignment, in array form}
  \item{align2}{Second alignment, in array form}
  \item{cols1}{Number of alignment positions in alignment 1}
  \item{cols2}{Number of alignment positions in alignment 2}
  \item{ref}{Name or index of the reference alignment}
  \item{i}{Index of the first alignment}
  \item{j}{Index of the second alignment}
  \item{name1}{Original name of the first alignment}
  \item{name2}{Original name of the second alignment}
}
\details{
comp.mat is the key to pairwise comparisons in ComPairWise.  It determines which of align1 and align2 should be treated as the reference in this comparison, by seeing whether either of them is actually ref.  It then compares them column by column.
}
\value{
A logical vector of length cols1 or cols2, depending on which was used as the reference alignment.  TRUE for identical columns, FALSE for non-identical columns, NA for all gap/missing columns.
}
\author{TER}
\seealso{\code{\link{comp.all}} for multiple (non-pairwise) comparisons}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (align1, align2, cols1, cols2, ref, i, j, name1, name2) 
{
    all.na <- function(column) all(is.na(column))
    if (ref == "longest") {
        if (cols2 > cols1) {
            ref.align <- align2
            other.align <- align1
        }
        else {
            ref.align <- align1
            other.align <- align2
        }
    }
    else {
        if (ref == name1 || ref == i) {
            ref.align <- align1
            other.align <- align2
        }
        if (ref == name2 || ref == j) {
            ref.align <- align2
            other.align <- align1
        }
    }
    ref.align <- ref.align[, which(apply(ref.align, 2, all.na) == 
        FALSE)]
    other.align <- other.align[, which(apply(other.align, 2, 
        all.na) == FALSE)]
    comp.col <- function(column) {
        row <- which(column > 0)[1]
        if (!is.na(row)) {
            col <- which(other.align[row, ] == column[row])
            comp <- identical(column, other.align[, col])
        }
        else {
            comp <- NA
        }
    }
    col.ident <- apply(ref.align, 2, comp.col)
  }
}

\keyword{internal}