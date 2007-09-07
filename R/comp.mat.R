`comp.mat` <-
function (align1, align2, cols1, cols2, ref, i, j, name1, name2) 
{
    all.na <- function(column) all(is.na(column))
    if (ref == "longest") {
        if (cols2 > cols1) {
            ref.align <- align2
            other.align <- align1
        }
        else {
            ref.align <- align1
            other.align <- align2
        }
    }
    else {
        if (ref == name1 || ref == i) {
            ref.align <- align1
            other.align <- align2
        }
        if (ref == name2 || ref == j) {
            ref.align <- align2
            other.align <- align1
        }
    }
    ref.align <- ref.align[, which(apply(ref.align, 2, all.na) == 
        FALSE)]
    other.align <- other.align[, which(apply(other.align, 2, 
        all.na) == FALSE)]
    comp.col <- function(column) {
        row <- which(column > 0)[1]
        if (!is.na(row)) {
            col <- which(other.align[row, ] == column[row])
            comp <- identical(column, other.align[, col])
        }
        else {
            comp <- NA
        }
    }
    col.ident <- apply(ref.align, 2, comp.col)
}
