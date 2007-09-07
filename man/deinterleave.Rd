\name{deinterleave}
\alias{deinterleave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert interleaved alignments to sequential}
\description{
Takes an interleaved alignment and turns it into a sequential one.
}
\usage{
deinterleave(alignment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{An alignment object, presumably one with every taxon split into more than one sequence}
}
\details{
This is a workaround function for reading alignments.  It re-assembles interleaved alignments, which otherwise would come out with every line of sequence as a different taxon.

deinterleave requires that all the parts of one sequence have the same taxon name, and that otherwise no taxon names be duplicate.  It concatenates sequences based only on their taxon names.
}
\value{
  ~Describe the value returned
An alignment object with the input sequences reassembled.
  \item{nb}{Number of sequences in the alignment}
  \item{nam}{Names of the sequences}
  \item{seq}{The sequences themselves}
  \item{com}{Currently NA; comment}
}
\author{TER}
\note{ 
Designed as an internal function but works as a standalone.  Called by \code{\link{read.nexus}}
}
\seealso{\code{\link{read.nexus}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (alignment) 
{
    taxa <- unique(alignment$nam)
    deint <- function(taxon, alignment) paste(alignment$seq[alignment$nam == 
        taxon], sep = "", collapse = "")
    x <- unlist(lapply(taxa, deint, alignment))
    newaln <- list(nb = length(taxa), nam = taxa, seq = x, com = NA)
    class(newaln) <- "alignment"
    return(newaln)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
