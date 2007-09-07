\name{aln.to.matrix}
\alias{aln.to.matrix}
\title{Convert an alignment to a taxa-by-characters matrix}
\description{
Convert an alignment object to a matrix of individual characters.  
}
\usage{
aln.to.matrix(alignment, taxa = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{An object of class alignment, usually the output of various alignment-reading functions in seqinr or ComPairWise such as read.nexus and read.phylip.  Or, any object with a component named seq containing data as a character vector, and (if taxa =TRUE) a component named nam containing the taxon names as a character vector.}
  \item{taxa}{Logical; should taxon names be returned as row names of the matrix?}
}
\details{
aln.to.matrix uses \code{\link{strsplit}} to split sequences into single-character strings.  It was designed for DNA, but should work with any datatype whose character states are single characters (e.g. "1" is fine, but "10" will get split into "1" and "0").  If taxa = TRUE the matrix will retain the taxon names as row names; otherwise the matrix will have no names.
}
\value{
A matrix with rows = taxa and columns = characters.  If taxa = TRUE the row names will be taxon names.

}
\author{TER}
\note{
Designed for objects of class "alignment", but should work with anything as long as there's a component named seq and (if taxa = TRUE) a component named nam.}
\examples{

data(sample.aln)
sample.mat <- aln.to.matrix(sample.aln, taxa = TRUE)

## The function is currently defined as
function (alignment, taxa = FALSE) 
{
    dna.only <- alignment$seq
    dna.mat <- t(sapply(strsplit(dna.only, split = ""), as.matrix))
    if (taxa) 
        rownames(dna.mat) <- alignment$nam
    return(dna.mat)
}

}
\keyword{manip}