`cpw` <-
function (name1 = "", name2 = "", ref = "", n = 2, out = "screen", 
    outfile = "", graph = "screen", gr.file = "", pairwise = TRUE, 
    onegraph = TRUE, nexus.out = TRUE, keep = FALSE, gap = "-", 
    miss = "?", ...) 
{
    cpw.out <- c()
    if (graph == "" || graph == "none" || graph == FALSE) {
        graph <- FALSE
    }
    else {
        gr.type <- graph
        graph <- TRUE
    }
    if ((out == "screen" || out == "") && outfile != "") 
        out <- "file"
    if (out == "file") {
        if (outfile == "") 
            outfile <- "compairwise.out"
        sink(outfile)
    }
    name <- aln.files(name1, name2, n)
    if (ref == "") {
        ref <- readline("If you want to define a reference alignment, which one? (or leave blank to use longest alignment as ref) ")
    }
    if ((ref %in% name || ref %in% (1:n)) == FALSE) 
        ref <- "longest"
    for (i in 1:n) {
        writeLines(paste("Alignment ", i, ": ", name[i], sep = ""))
    }
    writeLines(paste("Referencing to ", ref, "\n", sep = ""))
    formats <- sapply(name, check.format)
    bases <- c("A", "G", "C", "T", "a", "c", "g", "t", "U", "u")
    ambig.bases <- c("Y", "R", "K", "M", "W", "S", "B", "H", 
        "D", "V", "N", "y", "r", "k", "m", "w", "s", "b", "h", 
        "d", "v", "n")
    all.bases <- c(bases, ambig.bases)
    non.bases <- c(gap, miss)
    align.count <- c()
    for (i in 1:n) {
        if (formats[i] == "nexus") 
            aln <- read.nexus(name[i])
        if (formats[i] == "phylip") 
            aln <- read.phylip(name[i])
        if (formats[i] == "object") 
            aln <- get(name[i])
        dna.mat <- aln.to.matrix(aln, taxa = FALSE)
        if (all(dna.mat %in% c(all.bases, non.bases) == FALSE)) 
            writeLines(paste("\n", "Warning! Non-DNA, non-gap character(s) found in ", 
                name[i], " ...", sep = ""))
        align.new <- vectorize.alignment(dna.mat, all.bases, 
            non.bases)
        align.count[i] <- ncol(align.new)
        if (i == 1) {
            longest <- 1
            align <- array(c(align.new, rep(NA, times = (ncol(align.new) * 
                nrow(align.new) * (n - 1)))), c(nrow(align.new), 
                ncol(align.new), n))
        }
        else {
            if (ncol(align.new) > ncol(align[, , 1])) {
                longest <- i
                align.temp <- array(NA, c(nrow(align.new), ncol(align.new), 
                  n))
                for (j in 1:(i - 1)) {
                  align.temp[, , j] <- c(align[, , j], rep(NA, 
                    times = (ncol(align.temp[, , 1]) - ncol(align[, 
                      , 1])) * nrow(align.temp[, , 1])))
                }
                align <- align.temp
            }
            align[, , i] <- c(align.new, rep(NA, times = (ncol(align[, 
                , 1]) - ncol(align.new)) * nrow(align.new)))
        }
    }
    if (pairwise == TRUE) {
        if (ref == "longest") {
            i.set <- 1:(n - 1)
        }
        else {
            if (ref %in% 1:n) {
                i.set <- as.numeric(ref)
                j.set <- (1:n)
                j.set <- j.set[which(j.set != ref)]
            }
            else {
                if (ref %in% name) {
                  i.set <- which(name == ref)
                  j.set <- (1:n)
                  j.set <- j.set[which(j.set != ref)]
                }
            }
        }
        for (i in i.set) {
            if (ref == "longest") {
                j.set <- (i + 1):n
            }
            else {
                j.set <- (1:n)[which((1:n) != i.set)]
            }
            for (j in j.set) {
                col.ident <- comp.mat(align1 = align[, , i], 
                  align2 = align[, , j], cols1 = align.count[i], 
                  cols2 = align.count[j], ref, i = i, j = j, 
                  name1 = name[i], name2 = name[j])
                ident.set <- which(col.ident == TRUE)
                diff.set <- which(col.ident == FALSE)
                num.ident <- as.numeric(col.ident)
                end.points <- endpoints(ident.set, diff.set)
                totals <- data.frame(cols1 <- align.count[i], 
                  cols2 <- align.count[j], cols = length(col.ident), 
                  idents = length(which(col.ident)), nonidents = length(which(!col.ident)), 
                  all_gap_missing = length(which(is.na(col.ident))))
                colnames(totals) <- c("cols_align1", "cols_align2", 
                  "cols_compared", "n_ident", "n_diff", "n_gap_miss")
                id.s <- long.string(col.ident)
                writeLines(paste("\nComparing ", name[i], " and ", 
                  name[j], ".", sep = ""))
                writeLines(paste("\n", "As a string (+ matches, - doesn't, o is all gap/missing):", 
                  "\n", id.s, "\n", sep = ""))
                writeLines(paste("Columns in alignment 1:", "\n", 
                  totals$cols_align1, sep = ""))
                writeLines(paste("Columns in alignment 2:", "\n", 
                  totals$cols_align2, sep = ""))
                writeLines(paste("Columns considered:", "\n", 
                  totals$cols_compared, sep = ""))
                writeLines(paste("Identical columns:", "\n", 
                  totals$n_ident, sep = ""))
                writeLines(paste("Different columns:", "\n", 
                  totals$n_diff, sep = ""))
                writeLines(paste("All gap/missing columns in ref alignment:", 
                  "\n", totals$n_gap_miss, sep = ""))
                if (graph) {
                  if (gr.file == "") {
                    if (outfile == "") 
                      gr.file <- "compairwise"
                    else gr.file <- outfile
                  }
                  if (onegraph == FALSE) {
                    if (i == i.set[1] && j == j.set[1]) 
                      gr.file.in <- gr.file
                    gr.file <- paste(gr.file.in, as.character(i), 
                      as.character(j), sep = "")
                  }
                  if (onegraph == FALSE || (onegraph == TRUE && 
                    i == i.set[1] && j == j.set[1])) {
                    if (gr.type == "ps" || gr.type == "postscript" || 
                      gr.type == "PS") 
                      postscript(file = paste(gr.file, ".ps", 
                        sep = ""), width = 10, height = 4)
                    if (gr.type == "eps" || gr.type == "EPS") 
                      postscript(file = paste(gr.file, ".eps", 
                        sep = ""), onefile = FALSE, width = 10, 
                        height = 4)
                    if (gr.type == "jpeg" || gr.type == "jpg" || 
                      gr.type == "JPG" || gr.type == "JPEG") 
                      jpeg(filename = paste(gr.file, ".jpg", 
                        sep = ""), width = 960, height = 480)
                    if (gr.type == "pdf" || gr.type == "PDF") 
                      pdf(file = paste(gr.file, ".pdf", sep = ""), 
                        width = 10, height = 4)
                    if (gr.type == "png" || gr.type == "PNG") 
                      png(filename = paste(gr.file, ".png", sep = ""), 
                        width = 960, height = 480)
                    if (gr.type == "bmp" || gr.type == "BMP") 
                      bmp(filename = paste(gr.file, ".bmp", sep = ""), 
                        width = 960, height = 480)
                    plot(c(1, dim(align)[2]), c(0, 1.1), axes = F, 
                      xaxt = "n", yaxt = "n", bty = "L", type = "n", 
                      xlab = "position in reference alignment", 
                      ylab = "", font.lab = 4, cex.lab = 0.75)
                  }
                  if (onegraph && n > 2) 
                    linecol <- matrix(rainbow(n^2), n, n)[i, 
                      j]
                  else linecol <- "blue"
                  lines(num.ident, lwd = 1.5, col = linecol)
                  axis(side = 1, las = 1, lwd = 1, at = c(1, 
                    end.points, length(num.ident)), font = 1, 
                    cex.axis = 0.75)
                  axis(side = 2, lty = 0, at = c(0, 1), lwd = 0, 
                    labels = c("diff", "same"), font = 3, las = 1, 
                    pos = 1, cex.axis = 0.75)
                  if (!onegraph) 
                    title(main = paste("alignment", as.character(i), 
                      "vs alignment", as.character(j), sep = " "), 
                      cex.main = 1, font.main = 4)
                  if (gr.type != "screen" && !onegraph) 
                    dev.off()
                }
                if (nexus.out) {
                  nexus.sets.block <- nexus.sets(ident.set, diff.set)
                  writeLines(paste("NEXUS sets block:", "\n", 
                    nexus.sets.block, sep = ""))
                }
                cpw.out <- rbind(cpw.out, list(totals = totals, 
                  diff = diff.set, ident = ident.set, lab = id.s, 
                  id.num = num.ident, id.log = col.ident))
            }
        }
        if (graph) 
            if (gr.type != "screen" && onegraph) 
                dev.off()
    }
    if (pairwise == FALSE) {
        if (ref == "longest") 
            ref.index <- longest
        if (ref %in% 1:n) 
            ref.index <- as.numeric(ref)
        if (ref %in% name) 
            ref.index <- which(name == ref)
        col.ident <- comp.all(align, n, ref.index)
        ident.set <- which(col.ident == TRUE)
        diff.set <- which(col.ident == FALSE)
        num.ident <- as.numeric(col.ident)
        end.points <- endpoints(ident.set, diff.set)
        for (i in 1:n) writeLines(paste("Columns in alignment ", 
            as.character(i), ": ", as.character(align.count[i]), 
            sep = ""))
        totals <- data.frame(cols = length(col.ident), idents = length(which(col.ident)), 
            nonidents = length(which(!col.ident)), all_gap_missing = length(which(is.na(col.ident))))
        colnames(totals) <- c("cols_compared", "n_ident", "n_diff", 
            "n_gap_miss")
        id.s <- long.string(col.ident)
        writeLines(paste("\n", "As a string (+ matches, - doesn't):", 
            "\n", id.s, "\n", sep = ""))
        writeLines(paste("Columns considered:", "\n", totals$cols_compared, 
            sep = ""))
        writeLines(paste("Identical columns:", "\n", totals$n_ident, 
            sep = ""))
        writeLines(paste("Different columns:", "\n", totals$n_diff, 
            sep = ""))
        writeLines(paste("All gap/missing columns in ref alignment:", 
            "\n", totals$n_gap_miss, sep = ""))
        if (graph) {
            if (gr.file == "") {
                if (outfile == "") 
                  gr.file <- "compairwise"
                else gr.file <- outfile
            }
            if (gr.type == "ps" || gr.type == "postscript" | 
                gr.type == "PS") 
                postscript(file = paste(gr.file, ".ps", sep = ""), 
                  width = 10, height = 4)
            if (gr.type == "eps" || gr.type == "EPS") 
                postscript(file = paste(gr.file, ".eps", sep = ""), 
                  onefile = FALSE, width = 9, height = 4)
            if (gr.type == "jpeg" || gr.type == "jpg" || gr.type == 
                "JPG" || gr.type == "JPEG") 
                jpeg(filename = paste(gr.file, ".jpg", sep = ""), 
                  width = 960, height = 480)
            if (gr.type == "pdf" || gr.type == "PDF") 
                pdf(file = paste(gr.file, ".pdf", sep = ""), 
                  width = 10, height = 4)
            if (gr.type == "png" || gr.type == "PNG") 
                png(filename = paste(gr.file, ".png", sep = ""), 
                  width = 960, height = 480)
            if (gr.type == "bmp" || gr.type == "BMP") 
                bmp(filename = paste(gr.file, ".bmp", sep = ""), 
                  width = 960, height = 480)
            plot(c(1, dim(align)[2]), c(0, 1.1), axes = F, xaxt = "n", 
                yaxt = "n", bty = "L", type = "n", xlab = "position in reference alignment", 
                ylab = "", font.lab = 4, cex.lab = 0.75)
            lines(num.ident, lwd = 1.5, col = "blue")
            axis(side = 1, las = 1, lwd = 1, at = c(1, end.points, 
                length(num.ident)), font = 1, cex.axis = 0.75)
            axis(side = 2, lty = 0, at = c(0, 1), lwd = 0, labels = c("diff", 
                "same"), font = 3, las = 1, pos = 1, cex.axis = 0.75)
            if (gr.type != "screen") 
                dev.off()
        }
        if (nexus.out) {
            nexus.sets.block <- nexus.sets(ident.set, diff.set)
            writeLines(paste("NEXUS sets block:", "\n", nexus.sets.block, 
                sep = ""))
        }
        cpw.out <- rbind(cpw.out, list(totals = totals, diff = diff.set, 
            ident = ident.set, lab = id.s, id.num = num.ident, 
            id.log = col.ident))
    }
    if (keep) {
        totals <<- totals
        filenames <<- name
        diff.set <<- diff.set
        id.s <<- id.s
        ident.set <<- ident.set
        num.ident <<- num.ident
        col.ident <<- col.ident
        align <<- align
    }
    if (out == "file") 
        sink()
    invisible(cpw.out)
}
