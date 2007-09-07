`read.nexus` <-
function (filename) 
{
    nex.in <- scan(file = filename, sep = "\n", what = "c", strip.white = T, 
        multi.line = F, quiet = T)
    while (length(grep("\\[[^]\\[]*]", nex.in)) > 0) nex.in <- gsub("\\[[^]\\[]*]", 
        "", nex.in)
    comm.start <- grep("\\[", nex.in)
    comm.end <- grep("]", nex.in)
    if (length(comm.start) > 0 && length(comm.end) > 0) {
        for (i in 1:(length(comm.start))) {
            if (comm.end[i] == comm.start[i] + 1) {
                nex.in[comm.start[i]] <- gsub("\\[.*$", "", nex.in[comm.start[i]])
                nex.in[comm.end[i]] <- gsub("^.*]", "", nex.in[comm.end[i]])
            }
            else {
                for (j in (comm.start[i] + 1):(comm.end[i] - 
                  1)) {
                  nex.in[j] <- ""
                }
                nex.in[comm.start[i]] <- gsub("\\[.*$", "", nex.in[comm.start[i]])
                nex.in[comm.end[i]] <- gsub("^.*]", "", nex.in[comm.end[i]])
            }
        }
    }
    nex.in <- gsub("'", "", nex.in)
    nex.in <- nex.in[grep("^.+$", nex.in)]
    tax1 <- (grep("matrix", nex.in, ignore.case = T)) + 1
    semicolons <- grep(";", nex.in, ignore.case = T)
    taxlast <- semicolons[semicolons > tax1][1] - 1
    nex.header <- nex.in[1:(grep("matrix", nex.in, ignore.case = T))[1]]
    format.start <- grep("format", nex.header, ignore.case = T)
    nex.format <- unlist(strsplit(gsub("=|;", " ", nex.header[format.start:semicolons[semicolons >= 
        format.start][1]]), split = " "))
    if ("interleave" %in% nex.format) {
        if (grep("interleave", nex.format) == length(nex.format)) {
            interleave <- TRUE
        }
        else {
            if (nex.format[grep("interleave", nex.format) + 1] == 
                "no") 
                interleave <- FALSE
            else interleave <- TRUE
        }
    }
    else interleave <- FALSE
    dna.only <- sub("^[^[:blank:]]+[[:blank:]]+", "", nex.in[tax1:taxlast])
    dna.only <- gsub("[[:blank:]]+", "", dna.only)
    extr.first <- function(vector) vector[1]
    taxon.labels <- sapply(strsplit(grep("^[^[:blank:]]+[[:blank:]]+", 
        nex.in[tax1:taxlast], value = TRUE), split = "[[:blank:]]"), 
        extr.first)
    nex.align <- list(nb = length(dna.only), nam = taxon.labels, 
        seq = dna.only, com = NA)
    class(nex.align) <- "alignment"
    if ("matchchar" %in% nex.format) 
        nex.align <- dots.to.bases(nex.align, matchchar = nex.format[grep("matchchar", 
            nex.format, ignore.case = T)[1] + 1])
    if (interleave) 
        nex.align <- deinterleave(nex.align)
    return(nex.align)
}
