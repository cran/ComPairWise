`endpoints` <-
function (ident.set, diff.set) 
{
    end.points <- c()
    len.all <- length(ident.set) + length(diff.set)
    end.points <- which((1:(len.all - 1) %in% ident.set & 2:(len.all) %in% 
        diff.set) | (1:(len.all - 1) %in% diff.set & 2:(len.all) %in% 
        ident.set))
}
