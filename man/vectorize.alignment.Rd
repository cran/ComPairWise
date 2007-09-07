\name{vectorize.alignment}
\alias{vectorize.alignment}
\title{Turn DNA/RNA alignment into numbers}
\description{
Takes an alignment of DNA/RNA, and converts each data line into integers from 1 to nbases.  Internal ComPairWise function.
}
\usage{
vectorize.alignment(dna.mat, all.bases, non.bases)
}
\arguments{
  \item{dna.mat}{data matrix; rows are taxa, columns are characters.}
  \item{all.bases}{Character vector containing everything that should be considered a base (including ambiguities).}
  \item{non.bases}{Character vector containing everything that should be considered a gap or missing data.}
}
\details{
vectorize.alignment is an internal ComPairWise function.  It turns each sequence of bases into integers, from 1 to the number of bases, with 0 for gap or missing data.

A "vectorized" alignment maintains information on the position of the base in the sequence, but not of the character state.
}
\value{
A matrix the same size and shape as dna.mat (the input), but integers rather than characters.
}
\author{TER}
\note{
Not designed for standalone use, but should work.
}
\examples{
## The function is currently defined as
function (dna.mat, all.bases, non.bases) 
{
    vect.row <- function(row) {
        baseno <- 1
        dna.vect.row <- matrix(nrow = 1, ncol = length(row))
        for (j in 1:length(row)) {
            if (row[j] \%in\% all.bases) {
                dna.vect.row[j] <- baseno
                baseno <- baseno + 1
            }
            else {
                if (row[j] \%in\% non.bases) 
                  dna.vect.row[j] <- 0
                else {
                  dna.vect.row[j] <- -1
                  writeLines(paste("\t\t\t", as.character(row[j]), 
                    "in ref column", as.character(j)))
                }
            }
            dna.vect.row <- dna.vect.row
        }
    }
    dna.vect <- t(apply(dna.mat, 1, vect.row))
  }
}
\keyword{internal}