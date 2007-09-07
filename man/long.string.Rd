\name{long.string}
\alias{long.string}
\title{Make the +-o string}
\description{
Internal function in ComPairWise.  Makes a string of "+", "-", and "o" describing which columns in the alignments were the same, different, or all gap/all missing.
}
\usage{
long.string(col.ident)
}
\arguments{
  \item{col.ident}{Logical vector; output of \code{\link{comp.mat}}}
}
\details{
Generates a single character string in which identical columns are represented with "+", different columns with "-", and all gap/all missing columns with "o".  Designed to be copied/pasted into a text alignment.
}
\value{
A single string, as long as the reference alignment.
}
\author{TER}

\examples{

x<-c(TRUE, FALSE, TRUE, TRUE, NA, FALSE)
y<-long.string(x); y	#should return "+-++o-"


## The function is currently defined as
function (col.ident) 
{
    ident.string <- col.ident
    ident.string[which(col.ident)] <- "+"
    ident.string[which(!col.ident)] <- "-"
    ident.string[which(is.na(col.ident))] <- "o"
    id.s <- paste(ident.string, sep = "", collapse = "")
    return(id.s)
  }
}

\keyword{manip}