\name{check.format}
\alias{check.format}
\title{Determine DNA/RNA alignment format}
\description{
Determines whether an alignment file is NEXUS, PHYLIP, or neither, by looking at the first line of the file.  Internal function in ComPairWise.
}
\usage{
check.format(filename)
}

\details{
Used to try to figure out whether to call read.nexus or read.phylip, or to use an existing object.  It first checks whether an object with the name filename already exists and has class alignment; if so it assumes that's the object requested.

If there is no appropriate existing object, it uses scan to check the first line of the file filename.  If the first line contains "\#NEXUS" or "begin", it concludes that the file is a nexus file.  If the first line contains two numbers, it concludes that it's a phylip file.  

Designed as an internal function for ComPairWise, but could be useful in other contexts.
}
\arguments{
  \item{filename}{String.  Can be either the name of a file, or the name (as a string) of an existing object.}
}
\value{
A single string: "nexus" if the name seems to refer to a file in NEXUS format; "phylip" if the name seems to refer to a file in PHYLIP format; "object" if the name seems to refer to an existing object of class "alignment".  Otherwise, exits with an error (and doesn't return anything).
}
\author{TER}

\note{
check.format depends mostly on what's in the first line of the file.  It will usually fail to identify the file format if this is something unexpected, including a blank line.}
\examples{
## The function is currently defined as
function (filename) 
{
    options(warn = -1)
    on.exit(options(warn = 0))
    if (exists(filename)) {
        if (class(get(filename)) == "alignment") 
            aln.format <- "object"
    }
    else {
        is.nexus <- length(grep("#nexus|begin", scan(file = filename, 
            what = "c", nlines = 1, quiet = T), ignore.case = T)) > 
            0
        if (is.nexus) 
            is.phylip <- FALSE
        if (!is.nexus) 
            is.phylip <- (!is.na(as.numeric(scan(file = filename, 
                nlines = 1, quiet = T, what = "c"))) && length(scan(file = filename, 
                nlines = 1, quiet = T)) == 2)
        if (!is.nexus && !is.phylip) {
            aln.format <- "neither"
            stop(paste(filename, "isn't a NEXUS or a PHYLIP file!  Exiting."), 
                call. = FALSE)
        }
        if (is.nexus) 
            aln.format <- "nexus"
        if (is.phylip) 
            aln.format <- "phylip"
    }
    return(aln.format)
  }
}
\keyword{file}