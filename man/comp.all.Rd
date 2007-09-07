\name{comp.all}
\alias{comp.all}
\title{Groupwise comparison of alignment matrices}
\description{
comp.all is an internal ComPairWise function.  It compares a group (greater than 2) of matrices all at once, to determine which columns are identical in all matrices and which aren't.
}
\usage{
comp.all(align, n, ref.index)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{align}{a three-dimensional array of n taxa x character alignments}
  \item{n}{how many matrices (should be the same as dim(align)[3]}
  \item{ref.index}{integer; which alignment is being treated as the reference alignment?}
}
\details{
comp.all is the analog to \code{\link{comp.mat}} for non-pairwise comparisons when n > 3.  It is an internal function and isn't designed to be called as a standalone.
}
\value{
A logical vector of length ncols(align[,,ref.index]), i.e. the number of columns in the reference alignment.  TRUE for identical columns, FALSE for non-identical columns, NA for all gap/missing columns.
}
\author{ TER}
\seealso{ \code{\link{comp.mat}} for pairwise comparisons}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (align, n, ref.index) 
{
    comp.cols <- function(column) {
        if (any(!is.na(column))) {
            ref.row <- which(column > 0)[1]
            if (!is.na(ref.row)) {
                comp.vec <- c()
                for (i in 1:n) {
                  col <- which(align[ref.row, , i] == column[ref.row])
                  comp.vec <- c(comp.vec, identical(column, align[, 
                    col, i]))
                }
                comp <- all(comp.vec)
            }
            else {
                comp <- NA
            }
        }
        else {
            comp <- NA
        }
    }
    ref.align <- align[, , ref.index]
    col.ident <- apply(align[, , ref.index], 2, comp.cols)
  }
}

\keyword{internal}