`check.format` <-
function (filename) 
{
    options(warn = -1)
    on.exit(options(warn = 0))
    if (exists(filename)) {
        if (class(get(filename)) == "alignment") 
            aln.format <- "object"
    }
    else {
        is.nexus <- length(grep("#nexus|begin", scan(file = filename, 
            what = "c", nlines = 1, quiet = T), ignore.case = T)) > 
            0
        if (is.nexus) 
            is.phylip <- FALSE
        if (!is.nexus) 
            is.phylip <- (!is.na(as.numeric(scan(file = filename, 
                nlines = 1, quiet = T, what = "c"))) && length(scan(file = filename, 
                nlines = 1, quiet = T)) == 2)
        if (!is.nexus && !is.phylip) {
            aln.format <- "neither"
            stop(paste(filename, "isn't a NEXUS or a PHYLIP file!  Exiting."), 
                call. = FALSE)
        }
        if (is.nexus) 
            aln.format <- "nexus"
        if (is.phylip) 
            aln.format <- "phylip"
    }
    return(aln.format)
}
