`nexus.sets` <-
function (ident.set, diff.set) 
{
    if (!is.na(ident.set[1])) 
        nexus.ident.set.line <- paste("charset ident_align=", 
            paste(ident.set, sep = "", collapse = " "), ";", 
            sep = "")
    else nexus.ident.set.line <- ""
    if (!is.na(diff.set[1])) 
        nexus.diff.set.line <- paste("charset diff_align=", paste(diff.set, 
            sep = "", collapse = " "), ";", sep = "")
    else nexus.diff.set.line <- ""
    nexus.sets.block <- paste("begin sets;", "\n", nexus.ident.set.line, 
        "\n", nexus.diff.set.line, "\n", "end;", "\n", sep = "")
    return(nexus.sets.block)
}
