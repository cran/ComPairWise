\name{read.nexus}
\alias{read.nexus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read a NEXUS-format character alignment}
\description{
Reads a NEXUS-format alignment file, and returns it as an alignment object.
}
\usage{
read.nexus(filename)
}
\arguments{
  \item{filename}{Character string; name of the file to be read.}
}
\details{
read.nexus reads the file specified by filename, using scan.  It goes through several manipulations with grep and gsub to remove comments, whitespace, and taxon names, and extract the actual data.  It returns the data as an alignment object.  Because it depends heavily on interpretations of regular expressions, it is possible for your locale to affect the results.

Data is assumed to be everything between the first instance of the word "matrix" and the next semicolon.  All other blocks are ignored, as is header information.

The alignment may be sequential or interleaved.  The function will decide which one based on the "interleave" option in the format command (if there is one).

In sequential format, each line in the matrix should be

<taxon-name-1><whitespace><all the data for taxon 1>

In interleaved format, all lines for a given taxon must have exactly the same taxon name.  

If "matchchar" is specified in the format command of the NEXUS file, the function will use the specified character to replace all matches with the base in the first taxon.

Whitespace separating taxon names and data can be any number or mixture of spaces and tabs.  The taxon name is assumed to start at the beginning of the line and cannot contain spaces or tabs.  Single-quote characters are OK but will be removed.  It might be wise to stay away from special characters and punctuation, although they should be fine.  Anything that's usually legal in NEXUS taxa should be OK here.

Whitespace in the sequences themselves is OK, and will be removed once the taxon names have been dealt with.

Comments are acceptable, but all comments (including those inside other comments) must be terminated.

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
May fail unpredictably if square brackets (which should indicate comments) aren't matched, particularly for unmatched brackets within other comments.  Parsing depends heavily on \code{\link{grep}} and \code{\link{gsub}} and on interpretation of regular expressions, so odd things may happen if these find unexpected characters.  This can be a problem for taxon names, or if scan misinterprets ends of lines.
}
\seealso{\code{\link{check.format}} for automatic guessing what the format is, \code{\link{read.phylip}} for phylip files, \code{\link[seqinr]{read.alignment}} for some other formats}
\examples{
\dontrun{
oldwd <- getwd()
setwd(file.path(.Library, "ComPairWise", "examples"))
nex.aln <- read.nexus("sample.nex")		
setwd(oldwd)
rm(oldwd)
}}
\keyword{file}