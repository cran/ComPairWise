\name{cpw}
\alias{cpw}
\title{Compare alternate alignments of phylogenetic or population genetic data}
\description{
cpw is the main (outer) function of ComPairWise.  It takes 2 or more DNA/RNA alignment files and compares them, column by column.  It outputs information about which columns are identical and which aren't.
}
\usage{
cpw(name1 = "", name2 = "", ref = "", n = 2, out = "screen", outfile = "", graph = "screen", gr.file = "", pairwise = TRUE, onegraph = TRUE, nexus.out = TRUE, keep = FALSE, gap = "-", miss = "?", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name1}{String.  Name of the first alignment}
  \item{name2}{String.  Name of the second alignment}
  \item{ref}{Integer or string.  Which alignment should be treated as the reference?}
  \item{n}{Integer.  How many alignments?}
  \item{out}{String.  Where should the text output go?}
  \item{outfile}{String.  What should the output file (if any) be named?}
  \item{graph}{String or FALSE.  Where should graphics output go?}
  \item{gr.file}{String.  What should the graphics output file (if any) be named?}
  \item{pairwise}{Logical.  If n > 2, do you want pairwise comparisons or a single group comparison of alignments?}
  \item{onegraph}{Logical.  If n > 2 and pairwise = TRUE, do you want graphical output on one graph or multiple graphs?}
  \item{nexus.out}{Logical.  Do you want a NEXUS-format sets block of identical and different columns?}
  \item{keep}{Logical.  Should R do a global assignment of certain variables, thus leaving them in the workspace after the function is done?}
  \item{gap}{String.  What is the gap character?}
  \item{miss}{String.  What is the missing data character?}
  \item{\dots}{Currently does nothing}
}
\details{
cpw is the function implementing ComPairWise, pairwise or group comparison of alternative alignments of a give DNA/RNA sequence data set.  The method is described in Roberts & Olson.

Command-line options have been designed so that users with no R experience should be able to get the output they want easily.  

Some options do interact, and the program will ignore conflicting instructions if necessary.

A more complete description of the options and arguments:

\item{\option{name1, name2}}{Default: "".  If left blank (the default), the program will prompt for these later.  For repetitive or batch runs, they can be specified here.

They can be either valid filenames, or the names (as strings) of R objects of class alignment, e.g. those created by a call to \code{\link[seqinr]{read.alignment}}.  Objects and filenames can be freely mixed.}

\item{\option{ref}}{Default: "".  If an integer, the index of the alignment to treat as the reference, e.g. 1 for the first alignment.  If a string, the name of the alignment to treat as the reference.  If this is left blank (the default) the program will prompt for it later on; if nothing is entered at the prompt, it will use the longest alignment as the reference.}

\item{\option{n}}{Default: 2.  How many alignments are being compared?}

\item{\option{out}}{Default: "screen".  Valid options are "screen", "file", and "" (no text output).  Interacts with \option{outfile}, and will be automatically set to "file" if \option{outfile} is given.  If \option{out} = "file" and \option{outfile} = "", the filename will be "compairwise.out".}

\item{\option{outfile}}{Default: "".  Name of file for text output.  Interacts with \option{out}; if \option{outfile} is specified, it will force output to go to a file no matter what was specified for \option{out}.  This often means that only one or the other needs to be set.}

\item{\option{graph}}{Default: "screen".  Valid options are "screen", "" or "none" or FALSE (for no graphing), and options for various types of files: "jpeg", "jpg", "ps", "postscript", "eps", "png", "bmp", "pdf".  \option{graph} = "screen" will make R open one or more graphics devices (depending on \option{pairwise} and \option{onegraph}).  For file output it will open a graphics device depending on the file type desired.  File types not supported by your version of R will, obviously, not work; similarly, it would be wise to set \option{graph} = "" on a system that doesn't support graphics at all.} 

\item{\option{gr.file}}{Default: "".  Base filename for graphics output, if /option{graph} is set to send graphics to a file.  A file extension (which depends on the file type) will be appended to this.  If \option{onegraph} = FALSE and multiple graphs are generated, the indices of the alignments in each comparison will also be appended.}

\item{\option{pairwise}}{Default: TRUE.  If n > 2, you have three options for how to do comparisons, depending on the settings for \option{ref} and \option{pairwise}.  If \option{pairwise} = FALSE, cpw will do only a group comparison to determine which columns are identical in all alignments.  If \option{pairwise} = TRUE, it will do all (n * (n - 1))/2 pairwise comparisons if no reference alignment is specified, or all n - 1 comparisons to the reference alignment if one is specified.  If \option{pairwise} = TRUE, \option{onegraph} controls whether you get one graph or many.}

\item{\option{onegraph}}{Default: TRUE.  If multiple comparisons are being done, should graphical output all be on one graph, or separate graphs?}

\item{\option{nexus.out}}{Default: TRUE.  Should the text output include a NEXUS-format sets block with a character set for identical and non-identical columns?}

\item{\option{keep}}{Default: FALSE.  If \option{keep} = TRUE, cpw will make global assignments of some variables with the <<- operator, thus leaving them in the workspace after the function exits.  Advanced R users might want to use this option for further analysis of their alignments.  However, many of the same objects can be kept by assigning the value of the function cpw() instead, and this is a better way if all you want is summary information.}

\item{\option{gap}}{Default: "-".  A string (one or more values) telling cpw what values indicate an alignment gap.  Any non-base character is acceptable.}

\item{\option{miss}}{Default: "?".  A string (one or more values) telling cpw what values indicate missing data.  Any non-base character is acceptable.  Note however that "N" is treated as an ambiguous base and should not be used for missing data, which is treated quite differently!}

\item{\option{\dots}}{Currently does nothing.}


cpw implements the alignment comparison described in ComPairWise.  For every column of a reference alignment, it determines whether an identically-aligned column exists in an alternate alignment.

As well as the value (which is returned invisibly if no assignment is being made), the program will print some text output to the screen or a file (or neither), and will make a sketchy graph of the results to the screen or a file (or neither).

}
\value{
A list, returned invisibly if no assignment is being made.  For each comparison that was done:
  \item{totals}{Summary list:
  	\item{cols\_align...}{Number of columns in each alignment}
  	\item{cols\_compared}{Number of columns compared (depends on ref)}
  	\item{n\_ident}{Number of columns that were identical}
  	\item{n\_diff}{Number of columns that were different}
  	\item{n\_gap\_miss}{Number of columns that were all gap/missing}}
  \item{diff}{Indices of columns that were different}
  \item{ident}{Indices of columns that were identical}
  \item{lab}{String of +,-,o for identical, different, and all gap/missing columns, respectively}
  \item{id.num}{Numeric vector of 1 for identical columns, 0 for different, NA for all gap/missing}
  \item{id.log}{Logical vector of TRUE for identical columns, FALSE for different, NA for all gap/missing}
}
\references{Roberts, T.E. and L.E. Olson.  ComPairWise...}
\author{Trina E. Roberts, \email{trina.roberts@uaf.edu} }

\examples{

#The simplest case: 2 alignments, no options changed:
#The example below sets the working directory so it can find the example files

\dontrun{	
oldwd <- getwd()
setwd(file.path(.Library, "ComPairWise", "examples"))
cpw("sample.nex", "sample2.nex")
setwd(oldwd)
rm(oldwd)
}
}
\keyword{misc}