`aln.files` <-
function (name1, name2, n) 
{
    name.vector <- vector(mode = "character", length = n)
    name.vector[1] <- name1
    name.vector[2] <- name2
    for (i in 1:n) {
        if (name.vector[i] == "") 
            name.vector[i] <- readline(paste("Filename for alignment ", 
                i, "? ", sep = ""))
        if (name.vector[i] == "") 
            stop("Exiting--you didn't enter a filename!", call. = FALSE)
    }
    return(name.vector)
}
