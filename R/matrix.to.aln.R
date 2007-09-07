`matrix.to.aln` <-
function (matrix) 
{
    x <- unname(apply(matrix, 1, paste, collapse = ""))
    aln <- list(nb = nrow(matrix), nam = rownames(matrix), seq = x, 
        com = NA)
    class(aln) <- "alignment"
    return(aln)
}
