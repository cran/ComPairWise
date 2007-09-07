`read.phylip` <-
function (filename) 
{
    phy.in <- scan(file = filename, skip = 1, sep = "\n", what = "c", 
        strip.white = T, multi.line = T, quiet = T)
    if (length(grep("[[:blank:]]", phy.in)) > 0) {
        dna.only <- sub("^[^[:blank:]]+[[:blank:]]+", "", phy.in)
        extr.first <- function(vector) vector[1]
        taxon.labels <- sapply(strsplit(grep("^[^[:blank:]]+[[:blank:]]+", 
            phy.in, value = TRUE), split = "[[:blank:]]"), extr.first)
    }
    else {
        dna.only <- substr(phy.in, 11, nchar(phy.in))
        taxon.labels <- substr(phy.in, 1, 10)
    }
    dna.only <- gsub("[[:blank:]]+", "", dna.only)
    phy.align <- list(nb = length(dna.only), nam = taxon.labels, 
        seq = dna.only, com = NA)
    class(phy.align) <- "alignment"
    return(phy.align)
}
