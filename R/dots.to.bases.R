`dots.to.bases` <-
function (alignment, matchchar = ".") 
{
    dotsout <- function(aln.column) {
        aln.column[aln.column == matchchar] <- aln.column[1]
        return(aln.column)
    }
    aln.new <- aln.to.matrix(alignment, F)
    aln.new <- apply(aln.new, 2, dotsout)
    aln.new <- matrix.to.aln(aln.new)
    alignment$seq <- aln.new$seq
    return(alignment)
}
