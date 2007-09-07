`long.string` <-
function (col.ident) 
{
    ident.string <- col.ident
    ident.string[which(col.ident)] <- "+"
    ident.string[which(!col.ident)] <- "-"
    ident.string[which(is.na(col.ident))] <- "o"
    id.s <- paste(ident.string, sep = "", collapse = "")
    return(id.s)
}
