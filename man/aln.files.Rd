\name{aln.files}
\alias{aln.files}
\title{Prompt user for filenames if necessary, and combine all names into a vector.}
\description{
aln.files is an internal ComPairWise function.  It takes input filenames or objectnames either from the command line or (if not enough were given on the command line) by prompting the user, and combines them into a vector.
}
\usage{
aln.files(name1, name2, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name1}{String; first filename from the command line.}
  \item{name2}{String; second filename from the command line.}
  \item{n}{Integer; number of filenames needed, from the command line}
}
\details{
The behavior of aln.files depends on how much input the user gives when calling the main ComPairWise function, \code{\link{cpw}}.

If no filenames are given on the command line, file1 and file2 will both be "", and aln.files will prompt the user for n filenames.  It will not prompt for any filename != "".  If n = 2 and both filenames are given on the command line, it will not prompt at all.

If a prompt is offered and no filename is entered, the function exits with a message.

Once it has all the filenames, it returns them as a character vector.

This is an internal ComPairWise function.
}
\value{
A character vector of names of objects and/or files}
\author{TER}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

\dontrun{aln.files("test.nex", "", 2)	#prompts for the missing filename, and returns a character vector}

## The function is currently defined as
function(name1, name2, n){
	name.vector<-vector(mode="character", length=n)
	name.vector[1]<-name1
	name.vector[2]<-name2
	for (i in 1:n) {
		if (name.vector[i]=="") name.vector[i]<-readline(paste("Filename for alignment ", i, "? ", sep=""))
		if (name.vector[i]=="") stop("Exiting--you didn't enter a filename!", call.=FALSE)
	}
return(name.vector)
  }
}
\keyword{internal}