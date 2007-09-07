`deinterleave` <-
function (alignment) 
{
    taxa <- unique(alignment$nam)
    deint <- function(taxon, alignment) paste(alignment$seq[alignment$nam == 
        taxon], sep = "", collapse = "")
    x <- unlist(lapply(taxa, deint, alignment))
    newaln <- list(nb = length(taxa), nam = taxa, seq = x, com = NA)
    class(newaln) <- "alignment"
    return(newaln)
}
