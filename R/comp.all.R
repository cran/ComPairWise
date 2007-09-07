`comp.all` <-
function (align, n, ref.index) 
{
    comp.cols <- function(column) {
        if (any(!is.na(column))) {
            ref.row <- which(column > 0)[1]
            if (!is.na(ref.row)) {
                comp.vec <- c()
                for (i in 1:n) {
                  col <- which(align[ref.row, , i] == column[ref.row])
                  comp.vec <- c(comp.vec, identical(column, align[, 
                    col, i]))
                }
                comp <- all(comp.vec)
            }
            else {
                comp <- NA
            }
        }
        else {
            comp <- NA
        }
    }
    ref.align <- align[, , ref.index]
    col.ident <- apply(align[, , ref.index], 2, comp.cols)
}
