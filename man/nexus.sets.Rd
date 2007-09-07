\name{nexus.sets}
\alias{nexus.sets}
\title{Create NEXUS sets block}
\description{
Creates a NEXUS-format sets block defining charsets of alignment columns that are the same and different.  Designed as internal ComPairWise function, although useable in other contexts.
}
\usage{
nexus.sets(ident.set, diff.set)
}
\arguments{
  \item{ident.set}{numeric vector; indices of elements to include in identical set}
  \item{diff.set}{numeric vector; indices of elements to include in different set}
}
\details{
Given two vectors, generates a character set of each.  Designed to be mutually exclusive, but this isn't a requirement.  Output is a multi-line character string, designed to be printed with writeLines, that can be used as-is in any NEXUS application that supports the sets block, or in pieces in any NEXUS application that supports character sets in another block.
}
\value{
Multi-line character string
}
\author{TER}

\examples{
set1<-c(1,2,3,8,9,10)
set2<-c(4,5,6,7)
x<-nexus.sets(set1,set2)
writeLines(x)

## The function is currently defined as
function (col.ident, ident.set, diff.set) 
{
    if (!is.na(ident.set[1])) 
        nexus.ident.set.line <- paste("charset ident_align=", 
            paste(ident.set, sep = "", collapse = " "), ";", 
            sep = "")
    else nexus.ident.set.line <- ""
    if (!is.na(diff.set[1])) 
        nexus.diff.set.line <- paste("charset diff_align=", paste(diff.set, 
            sep = "", collapse = " "), ";", sep = "")
    else nexus.diff.set.line <- ""
    nexus.sets.block <- paste("begin sets;", "\n", nexus.ident.set.line, 
        "\n", nexus.diff.set.line, "\n", "end;", "\n", sep = "")
    return(nexus.sets.block)
  }
}
\keyword{utilities}