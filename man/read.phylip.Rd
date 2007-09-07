\name{read.phylip}
\alias{read.phylip}
\title{Read a NEXUS-format character alignment and return it as a matrix}
\description{
Reads a sequential (non-interleaved) PHYLIP-format alignment file, and returns it as a matrix with or without the taxon names as row names.  Designed for DNA but should work for anything as long as character states are single characters.
}
\usage{
read.phylip(filename)
}
\arguments{
  \item{filename}{Name of the file to be read.}
}
\details{
read.phylip reads the file specified by filename, using scan.  It skips the first line (which should not contain data) and assumes that the rest of the file contains data.  

It will try to guess how to parse taxon names by looking for whitespace in the data.  If it finds whitespace, it will guess that the taxon names are what comes before the first space.  If it finds no whitespace it will assume that the first ten characters of each line are the taxon name.

}
\value{
An object of class alignment, of the same format as alignments in the package seqinr.  It is a list with the following components:

  \item{nb}{The number of taxa}
  \item{nam}{The taxon names as a character vector}
  \item{seq}{The data as a character vector (data for each taxon is a single string)}
  \item{com}{Currently always NA; present for compatibility with other formats that might have comments}
}
\author{TER}
\note{
Will fail if it guesses incorrectly about where the character data start; the first data column must be either directly after the first whitespace (space or tab) in the line or the 11th character in the line.

For some PHYLIP files the phylip reader in seqinr works better.
}
\seealso{\code{\link{check.format}} for automated format checking, \code{\link{read.nexus}} for NEXUS, \code{\link[seqinr]{read.alignment}} for some other formats}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

\dontrun{read.phylip("sample.phy")		#could run this if the sample file sample.phy were in the current working directory}

## The function is currently defined as
function (filename) 
{
    phy.in <- scan(file = filename, skip = 1, sep = "\n", what = "c", 
        strip.white = T, multi.line = T, quiet = T)
    if (length(grep("[[:blank:]]", phy.in)) > 0) {
        dna.only <- sub("^\\w+( |\t)+", "", phy.in)
        extr.first <- function(vector) vector[1]
        taxon.labels <- sapply(strsplit(grep("^\\w+( |\t)+", 
            phy.in, value = TRUE), split = "[[:blank:]]"), extr.first)
    }
    else {
        dna.only <- substr(phy.in, 11, nchar(phy.in))
        taxon.labels <- substr(phy.in, 1, 10)
    }
    dna.only <- gsub("[[:blank:]]+", "", dna.only)
    phy.align <- list(nb = length(dna.only), nam = taxon.labels, 
        seq = dna.only, com = NA)
    class(phy.align) <- "alignment"
    return(phy.align)
  }
}
\keyword{file}