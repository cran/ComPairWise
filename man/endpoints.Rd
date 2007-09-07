\name{endpoints}
\alias{endpoints}
\title{Find the breakpoints between identical and different alignment columns}
\description{
  Finds identical columns followed by different columns, and different columns followed by identical columns.  Internal function to find breakpoints for annotating axes.
}
\usage{
endpoints(ident.set, diff.set)
}
\arguments{
  \item{ident.set}{Numeric indices of identical columns}
  \item{diff.set}{Numeric indices of non-identical columns}
}
\value{
A vector of integers.   
}
\author{TER}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

identical<-c(1:10,12,15:20)
different<-c(11,13:14,21:25)
a<-endpoints(identical,different);a	#should return c(10,11,12,14,20)

## The function is currently defined as
function (ident.set, diff.set) 
{
    end.points <- c()
    len.all <- length(ident.set) + length(diff.set)
    end.points <- which((1:(len.all - 1) \%in\% ident.set & 2:(len.all) \%in\% 
        diff.set) | (1:(len.all - 1) \%in\% diff.set & 2:(len.all) \%in\% 
        ident.set))
  }
}
\keyword{dplot}
