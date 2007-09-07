`vectorize.alignment` <-
function (dna.mat, all.bases, non.bases) 
{
    vect.row <- function(row) {
        baseno <- 1
        dna.vect.row <- matrix(nrow = 1, ncol = length(row))
        for (j in 1:length(row)) {
            if (row[j] %in% all.bases) {
                dna.vect.row[j] <- baseno
                baseno <- baseno + 1
            }
            else {
                if (row[j] %in% non.bases) 
                  dna.vect.row[j] <- 0
                else {
                  dna.vect.row[j] <- -1
                  writeLines(paste("\t\t\t", as.character(row[j]), 
                    "in ref column", as.character(j)))
                }
            }
            dna.vect.row <- dna.vect.row
        }
    }
    dna.vect <- t(apply(dna.mat, 1, vect.row))
}
