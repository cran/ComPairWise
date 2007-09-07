`aln.to.matrix` <-
function (alignment, taxa = FALSE) 
{
    dna.only <- alignment$seq
    dna.mat <- t(sapply(strsplit(dna.only, split = ""), as.matrix))
    if (taxa) 
        rownames(dna.mat) <- alignment$nam
    return(dna.mat)
}
