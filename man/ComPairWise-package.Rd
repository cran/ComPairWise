\name{ComPairWise-package}
\alias{ComPairWise-package}
\alias{ComPairWise}
\docType{package}
\title{
Compare phylogenetic or population genetic DNA/RNA alignments
}
\description{
ComPairWise compares alternate alignments of a single DNA/RNA data set, and determines which columns are identically aligned.  It also contains functions for reading NEXUS and PHYLIP-format alignment files. 
}
\details{
\tabular{ll}{
Package: \tab ComPairWise\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-09-07\cr
License: \tab GNU GPL 2.0 or greater\cr
}
Contents:

\item{cpw}{The basic function.  Various command-line options are available, all of which are optional.  If run with nothing on the command line, the program will assume you want to compare two alignments and print output to the screen, and will prompt for all other necessary information.  }

\item{read.nexus}{Read a nexus file and generate an object of class alignment.}

\item{read.phylip}{Read a phylip file and generate an object of class alignment.}

\item{aln.to.matrix}{Take the data from an alignment object and convert it into a taxa-by-characters matrix, with or without rownames.}

The package also includes code to read NEXUS and PHYLIP alignments, which can be used as stand-alone functions.
}
\author{
Trina E. Roberts

Maintainer: Trina E. Roberts <trina.roberts@uaf.edu>
}
\references{
Roberts, T.E. and L.E. Olson.  ComPairWise...
}
\examples{

#The simplest case: 2 alignments, no options changed:
#To run this, move the sample files from the package's data directory to the current working directory, or specify the path to the data directory
\dontrun{cpw("sample.nex", "sample1.nex")}



#other functions
#To run these, move the sample files from the package's data directory to the current working directory, or specify the path to the data directory
\dontrun{nex.aln <- read.nexus("sample.nex")}

\dontrun{phy.aln <- read.phylip("sample.phy")}

data(sample.aln)
aln1.mat <- aln.to.matrix(sample.aln)

}
\keyword{package}