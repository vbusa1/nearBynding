#' @title visualizeStereogene
#'
#' @description Creates a visual output of a single RNA structure context
#' relative to protein binding.
#'
#' @param dir_stereogene_output Directory of stereogene output. Default working
#' directory
#' @param context_file A single context file name for visualization with the
#' protein_file(s). File names must exclude extensions such as ".bedGraph".
#' Required.
#' @param protein_file A vector of at least one protein file name to be averaged
#' for visualization. File names must exclude extensions such as ".bedGraph".
#' All files in the list should be experimental or biological replicates.
#' Required.
#' @param protein_file_input A protein file name of background input to be
#' subtracted from protein_file signal. File name must exclude extension. Only
#' one input file is permitted. Optional.
#' @param x_lim A vector of two values denoting the lower and upper x axis
#' limits. Cannot exceed wSize/2 from write_config. Default (-100, 100)
#' @param y_lim A vector of two values denoting the lower and upper y axis
#' limits. Optional
#' @param out_file Name of output file, excluding extension. ".pdf" or ".jpeg"
#' will be added as relevant to the output file name. Default "out_file"
#' @param legend Whether a legend should be included with the output graph.
#' Default TRUE
#' @param heatmap Whether the output graph should be in the form of a heatmap
#' (TRUE) or of a line graph (FALSE). Default FALSE
#'
#' @export

visualizeStereogene <- function(dir_stereogene_output = ".",
                                context_file,
                                protein_file,
                                protein_file_input = NULL,
                                x_lim = c(-100, 100),
                                y_lim = NULL,
                                out_file = "out_file",
                                legend = TRUE,
                                heatmap = FALSE) {
    if (length(protein_file) < 1) {
        stop("Requires at least one protein file prefix to calculate distance")
    }
    if (length(protein_file) > 20) {
        stop("There are > 20 protein files input. This is likely in error")
    }
    if (length(protein_file_input) > 1) {
        stop("Only input one background track per protein.")
    }
    dist_1 <- NULL
    for (n in seq(length(protein_file))) {
        assign(paste0("dist_", n), read.table(paste0(dir_stereogene_output, "/",
                                            context_file, "~", protein_file[n],
                                            ".dist"), header = TRUE))
    }
    if (!is.null(protein_file_input)) {
        dist_input <- read.table(paste0(dir_stereogene_output, "/",
                                        context_file,
                                        "~", protein_file_input,
                                        ".dist"), header = TRUE)
    }
    dist <- as.data.frame(matrix(NA, ncol = (2 * length(protein_file)) +
                                                    1, nrow = nrow(dist_1)))
    colnames(dist) <- c("x", paste0("Fg", seq(length(protein_file))),
                        paste0("Bkg", seq(length(protein_file))))
    dist$x <- dist_1$x
    for (n in seq(length(protein_file))) {
        dist[, 1 + n] <- eval(parse(text = paste0("dist_", n)))$Fg
        dist[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_", n)))$Bkg
    }
    if (!is.null(protein_file_input)) {
        dist[, 2:(length(protein_file) + 1)] <- dist[,
                        2:(length(protein_file) + 1)] - dist_input$Fg
        dist[, (length(protein_file) + 2):ncol(dist)] <- dist[,
                        (length(protein_file) + 2):ncol(dist)] - dist_input$Bkg
    }
    if (length(protein_file) == 1) {
        dist$Fg <- dist$Fg1
        dist$Fg_se <- dist$Bkg_se <- 0
        dist$Bkg <- dist$Bkg1
    } else {
        dist$Fg <- rowMeans(dist[, 2:(length(protein_file) + 1)])
        dist$Fg_se <- rowSds(as.matrix(dist[, 2:(length(protein_file) +
                                        1)]))/sqrt(length(protein_file))
        dist$Bkg <- rowMeans(dist[, (length(protein_file) +
                                        2):((length(protein_file) * 2) + 1)])
        dist$Bkg_se <- rowSds(as.matrix(dist[, (length(protein_file) +
                                        2):((length(protein_file) * 2) + 1)]))/
            sqrt(length(protein_file))
    }
    if (heatmap == FALSE) {
        # create line plot
        out_file <- paste0(out_file, ".pdf")
        # save plot to pdf
        pdf(out_file, height = 4, width = 6)
        # create the plot
        old.par <- par(no.readonly = TRUE)
        par(oma = c(0, 0, 0, 0), mar = c(3, 3, 3, 1), mgp = c(1.6, 0.45, 0))
        if (is.null(y_lim)) {
            max_y <- max(dist$Fg) + 0.1
            min_y <- min(dist$Fg) - 0.1
            y_lim <- c(min_y, max_y)
        }
        plot(dist$x, dist$Bkg, type = "l", col = "black",
                xlim = x_lim, ylim = y_lim, main = NULL,
                xlab = "Distance", ylab = "Density x 100",
                cex.axis = 0.8, cex.lab = 1, cex.main = 1.2, lwd = 2)
        arrows(dist$x, dist$Bkg - dist$Bkg_se, dist$x,
                dist$Bkg + dist$Bkg_se, length = 0, angle = 90,
                code = 3, col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        abline(v = 0, col = "grey", lty = 2)
        lines(dist$x, dist$Fg, col = "blue", lwd = 2)
        arrows(dist$x, dist$Fg - dist$Fg_se, dist$x,
                dist$Fg + dist$Fg_se, length = 0, angle = 90,
                code = 3, col = rgb(blue = 1, red = 0, green = 0, alpha = 0.5))
        if (legend == TRUE) {
            legend("bottomright", legend = c("Background", "Signal"),
                    col = c("black", "blue"), lty = 1, cex = 0.8)
        }
        dev.off()
    } else {
        out_file <- paste0(out_file, ".jpeg")
        coords_1 <- which(dist$x == x_lim[1])
        coords_2 <- which(dist$x == x_lim[2])
        heatmap_plot <- data.frame(Signal = dist$Fg[coords_1:coords_2],
                                    Background = dist$Bkg[coords_1:coords_2])
        heatmap_plot <- as.matrix(t(heatmap_plot))
        zero <- (ncol(heatmap_plot) + 1)/2
        if (is.null(y_lim)) {
            max_y <- max(c(max(dist$Bkg), max(dist$Fg))) + 0.1
            min_y <- min(c(min(dist$Bkg), min(dist$Fg))) - 0.1
            y_lim <- c(-max(abs(max_y), abs(min_y)),
                        max(abs(max_y), abs(min_y)))
        }
        key <- TRUE
        if (legend == FALSE) {key <- FALSE}
        # make plot
        jpeg(out_file, height = 4, width = 15, units = "in", res = 100)
        heatmap.2(heatmap_plot, Rowv = NA, Colv = NA,
                col = redblue(256), dendrogram = "none",
                labCol = FALSE, trace = "none", density.info = "none",
                symkey = FALSE, key = key, colsep = c(zero - 1, zero),
                sepcolor = "grey", margins = c(1, 20),
                breaks = seq(y_lim[1], y_lim[2], length.out = 257))
        dev.off()
    }
}
