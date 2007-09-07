\name{read.phylip}
\alias{read.phylip}
\title{Read a PHYLIP-format character alignment}
\description{
Reads a PHYLIP-format alignment file and returns it as an alignment object}
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

\dontrun{
oldwd <- getwd()
setwd(file.path(.Library, "ComPairWise", "examples"))
read.phylip("sample.phy")
setwd(oldwd)
rm(oldwd)
}}
\keyword{file}