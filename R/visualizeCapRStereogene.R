#' @title visualizeCapRStereogene
#'
#' @description Creates a visual output of all CapR RNA structure contexts
#' relative to protein binding.
#'
#' @param dir_stereogene_output Directory of stereogene output. Default working
#' directory
#' @param CapR_prefix The prefix common to CapR output files of protein_file.
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
#' @return heatmap (JPEG) or line graph (PDF) image file
#'
#' @examples
#' ## pull example files
#'get_outfiles()
#'## heatmap
#'visualizeCapRStereogene(CapR_prefix = "chr4and5_3UTR",
#'                        protein_file = "chr4and5_liftOver",
#'                        heatmap = TRUE,
#'                        out_file = "all_contexts_heatmap",
#'                        x_lim = c(-500, 500))
#'## line graph
#'visualizeCapRStereogene(CapR_prefix = "chr4and5_3UTR",
#'                        protein_file = "chr4and5_liftOver",
#'                        x_lim = c(-500, 500),
#'                        out_file = "all_contexts_line",
#'                        y_lim = c(-18, 22))
#'
#' @importFrom utils read.table
#' @importFrom matrixStats rowSds
#' @importFrom graphics abline arrows lines par plot
#' @importFrom grDevices dev.off jpeg pdf rgb
#' @importFrom gplots heatmap.2 redblue
#'
#' @export

visualizeCapRStereogene <- function(dir_stereogene_output = ".",
                                    CapR_prefix,
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
    dist_bulge_1 <- NULL
    for (n in seq(length(protein_file))) {
        assign(paste0("dist_bulge_", n), read.table(paste0(
                                    dir_stereogene_output,
                                    "/", CapR_prefix, "_bulge_liftOver~",
                                    protein_file[n], ".dist"), header = TRUE))
        assign(paste0("dist_multibranch_", n), read.table(paste0(
                                    dir_stereogene_output,
                                    "/", CapR_prefix, "_multibranch_liftOver~",
                                    protein_file[n], ".dist"), header = TRUE))
        assign(paste0("dist_hairpin_", n), read.table(paste0(
                                    dir_stereogene_output,
                                    "/", CapR_prefix, "_hairpin_liftOver~",
                                    protein_file[n], ".dist"), header = TRUE))
        assign(paste0("dist_stem_", n), read.table(paste0(
                                    dir_stereogene_output,
                                    "/", CapR_prefix, "_stem_liftOver~",
                                    protein_file[n], ".dist"), header = TRUE))
        assign(paste0("dist_internal_", n), read.table(paste0(
                                    dir_stereogene_output,
                                    "/", CapR_prefix, "_internal_liftOver~",
                                    protein_file[n], ".dist"), header = TRUE))
        assign(paste0("dist_exterior_", n), read.table(paste0(
                                    dir_stereogene_output,
                                    "/", CapR_prefix, "_exterior_liftOver~",
                                    protein_file[n], ".dist"), header = TRUE))
    }
    if (!is.null(protein_file_input)) {
        dist_bulge_input <- read.table(paste0(dir_stereogene_output,
                                    "/", CapR_prefix, "_bulge_liftOver~",
                                    protein_file_input, ".dist"), header = TRUE)
        dist_multibranch_input <- read.table(paste0(dir_stereogene_output,
                                    "/", CapR_prefix, "_multibranch_liftOver~",
                                    protein_file_input, ".dist"), header = TRUE)
        dist_stem_input <- read.table(paste0(dir_stereogene_output,
                                    "/", CapR_prefix, "_stem_liftOver~",
                                    protein_file_input, ".dist"), header = TRUE)
        dist_hairpin_input <- read.table(paste0(dir_stereogene_output,
                                    "/", CapR_prefix, "_hairpin_liftOver~",
                                    protein_file_input, ".dist"), header = TRUE)
        dist_internal_input <- read.table(paste0(dir_stereogene_output,
                                    "/", CapR_prefix, "_internal_liftOver~",
                                    protein_file_input, ".dist"), header = TRUE)
        dist_exterior_input <- read.table(paste0(dir_stereogene_output,
                                    "/", CapR_prefix, "_exterior_liftOver~",
                                    protein_file_input, ".dist"), header = TRUE)
    }
    dist_bulge <- as.data.frame(matrix(NA,
                                    ncol = (2 * length(protein_file)) + 1,
                                    nrow = nrow(dist_bulge_1)))
    dist_stem <- as.data.frame(matrix(NA,
                                    ncol = (2 * length(protein_file)) + 1,
                                    nrow = nrow(dist_bulge_1)))
    dist_multibranch <- as.data.frame(matrix(NA,
                                    ncol = (2 * length(protein_file)) +1,
                                    nrow = nrow(dist_bulge_1)))
    dist_hairpin <- as.data.frame(matrix(NA,
                                    ncol = (2 * length(protein_file)) + 1,
                                    nrow = nrow(dist_bulge_1)))
    dist_exterior <- as.data.frame(matrix(NA,
                                    ncol = (2 * length(protein_file)) + 1,
                                    nrow = nrow(dist_bulge_1)))
    dist_internal <- as.data.frame(matrix(NA,
                                    ncol = (2 * length(protein_file)) + 1,
                                    nrow = nrow(dist_bulge_1)))
    colnames(dist_bulge) <- c("x", paste0("Fg", seq(length(protein_file))),
                                paste0("Bkg", seq(length(protein_file))))
    colnames(dist_stem) <- c("x", paste0("Fg", seq(length(protein_file))),
                                paste0("Bkg", seq(length(protein_file))))
    colnames(dist_internal) <- c("x", paste0("Fg", seq(length(protein_file))),
                                paste0("Bkg", seq(length(protein_file))))
    colnames(dist_multibranch) <- c("x", paste0("Fg",seq(length(protein_file))),
                                paste0("Bkg", seq(length(protein_file))))
    colnames(dist_exterior) <- c("x", paste0("Fg", seq(length(protein_file))),
                                paste0("Bkg", seq(length(protein_file))))
    colnames(dist_hairpin) <- c("x", paste0("Fg", seq(length(protein_file))),
                                paste0("Bkg", seq(length(protein_file))))
    dist_hairpin$x <-dist_exterior$x <-dist_multibranch$x <-dist_internal$x <-
        dist_stem$x <-dist_bulge$x <-dist_bulge_1$x
    for (n in seq(length(protein_file))) {
        dist_bulge[, 1 + n] <-
            eval(parse(text = paste0("dist_bulge_", n)))$Fg
        dist_stem[, 1 + n] <-
            eval(parse(text = paste0("dist_stem_",  n)))$Fg
        dist_hairpin[, 1 + n] <-
            eval(parse(text = paste0("dist_hairpin_", n)))$Fg
        dist_multibranch[, 1 + n] <-
            eval(parse(text = paste0("dist_multibranch_", n)))$Fg
        dist_exterior[, 1 + n] <-
            eval(parse(text = paste0("dist_exterior_", n)))$Fg
        dist_internal[, 1 + n] <-
            eval(parse(text = paste0("dist_internal_", n)))$Fg
        dist_bulge[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_bulge_", n)))$Bkg
        dist_stem[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_stem_", n)))$Bkg
        dist_hairpin[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_hairpin_", n)))$Bkg
        dist_multibranch[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_multibranch_", n)))$Bkg
        dist_exterior[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_exterior_", n)))$Bkg
        dist_internal[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_internal_", n)))$Bkg
    }
    if (!is.null(protein_file_input)) {
        dist_bulge[, 2:(length(protein_file) + 1)] <-
            dist_bulge[, 2:(length(protein_file) + 1)] - dist_bulge_input$Fg
        dist_bulge[, (length(protein_file) + 2):ncol(dist_bulge)] <-
            dist_bulge[, (length(protein_file) + 2):ncol(dist_bulge)] -
            dist_bulge_input$Bkg
        dist_stem[, 2:(length(protein_file) + 1)] <-
            dist_stem[, 2:(length(protein_file) + 1)] - dist_stem_input$Fg
        dist_stem[, (length(protein_file) + 2):ncol(dist_stem)] <-
            dist_stem[, (length(protein_file) + 2):ncol(dist_stem)] -
            dist_stem_input$Bkg
        dist_hairpin[, 2:(length(protein_file) + 1)] <-
            dist_hairpin[, 2:(length(protein_file) + 1)] -
            dist_hairpin_input$Fg
        dist_hairpin[, (length(protein_file) + 2):ncol(dist_hairpin)] <-
            dist_hairpin[, (length(protein_file) + 2):ncol(dist_hairpin)] -
            dist_hairpin_input$Bkg
        dist_multibranch[, 2:(length(protein_file) +  1)] <-
            dist_multibranch[,2:(length(protein_file) +  1)] -
            dist_multibranch_input$Fg
        dist_multibranch[, (length(protein_file) + 2):ncol(dist_multibranch)] <-
            dist_multibranch[,(length(protein_file)+2):ncol(dist_multibranch)] -
            dist_multibranch_input$Bkg
        dist_internal[, 2:(length(protein_file) + 1)] <-
            dist_internal[, 2:(length(protein_file) + 1)] -
            dist_internal_input$Fg
        dist_internal[, (length(protein_file) + 2):ncol(dist_internal)] <-
            dist_internal[, (length(protein_file) + 2):ncol(dist_internal)] -
            dist_internal_input$Bkg
        dist_exterior[, 2:(length(protein_file) + 1)] <-
            dist_exterior[, 2:(length(protein_file) + 1)] -
            dist_exterior_input$Fg
        dist_exterior[, (length(protein_file) + 2):ncol(dist_exterior)] <-
            dist_exterior[, (length(protein_file) + 2):ncol(dist_exterior)] -
            dist_exterior_input$Bkg
    }
    if (length(protein_file) == 1) {
        dist_bulge$Fg <- dist_bulge$Fg1
        dist_bulge$Fg_se <- dist_bulge$Bkg_se <- 0
        dist_bulge$Bkg <- dist_bulge$Bkg1
        dist_hairpin$Fg <- dist_hairpin$Fg1
        dist_hairpin$Fg_se <- dist_hairpin$Bkg_se <- 0
        dist_hairpin$Bkg <- dist_hairpin$Bkg1
        dist_internal$Fg <- dist_internal$Fg1
        dist_internal$Fg_se <- dist_internal$Bkg_se <- 0
        dist_internal$Bkg <- dist_internal$Bkg1
        dist_exterior$Fg <- dist_exterior$Fg1
        dist_exterior$Fg_se <- dist_exterior$Bkg_se <- 0
        dist_exterior$Bkg <- dist_exterior$Bkg1
        dist_stem$Fg <- dist_stem$Fg1
        dist_stem$Fg_se <- dist_stem$Bkg_se <- 0
        dist_stem$Bkg <- dist_stem$Bkg1
        dist_multibranch$Fg <- dist_multibranch$Fg1
        dist_multibranch$Fg_se <- dist_multibranch$Bkg_se <- 0
        dist_multibranch$Bkg <- dist_multibranch$Bkg1
    } else {
        dist_bulge$Fg <- rowMeans(dist_bulge[, 2:(length(protein_file) + 1)])
        dist_bulge$Fg_se <- rowSds(as.matrix(dist_bulge[,
                2:(length(protein_file) + 1)]))/sqrt(length(protein_file))
        dist_bulge$Bkg <- rowMeans(dist_bulge[,
                (length(protein_file) + 2):((length(protein_file) * 2) + 1)])
        dist_bulge$Bkg_se <- rowSds(as.matrix(dist_bulge[,
                (length(protein_file) + 2):((length(protein_file) *
                2) + 1)]))/sqrt(length(protein_file))
        dist_hairpin$Fg <- rowMeans(dist_hairpin[, 2:(length(protein_file) +
                                                                            1)])
        dist_hairpin$Fg_se <- rowSds(as.matrix(dist_hairpin[,
                2:(length(protein_file) + 1)]))/sqrt(length(protein_file))
        dist_hairpin$Bkg <- rowMeans(dist_hairpin[,
                (length(protein_file) + 2):((length(protein_file) * 2) + 1)])
        dist_hairpin$Bkg_se <- rowSds(as.matrix(dist_hairpin[,
                (length(protein_file) + 2):((length(protein_file) *
                2) + 1)]))/sqrt(length(protein_file))
        dist_stem$Fg <- rowMeans(dist_stem[, 2:(length(protein_file) + 1)])
        dist_stem$Fg_se <- rowSds(as.matrix(dist_stem[,
                2:(length(protein_file) + 1)]))/sqrt(length(protein_file))
        dist_stem$Bkg <- rowMeans(dist_stem[,
                (length(protein_file) + 2):((length(protein_file) * 2) + 1)])
        dist_stem$Bkg_se <- rowSds(as.matrix(dist_stem[,
                (length(protein_file) + 2):((length(protein_file) *
                2) + 1)]))/sqrt(length(protein_file))
        dist_multibranch$Fg<-rowMeans(dist_multibranch[,2:(length(protein_file)+
                                                                            1)])
        dist_multibranch$Fg_se <- rowSds(as.matrix(dist_multibranch[,
                2:(length(protein_file) + 1)]))/sqrt(length(protein_file))
        dist_multibranch$Bkg <- rowMeans(dist_multibranch[,
                (length(protein_file) + 2):((length(protein_file) * 2) + 1)])
        dist_multibranch$Bkg_se <- rowSds(as.matrix(dist_multibranch[,
                (length(protein_file) + 2):((length(protein_file) *
                2) + 1)]))/sqrt(length(protein_file))
        dist_internal$Fg <- rowMeans(dist_internal[, 2:(length(protein_file) +
                                                                            1)])
        dist_internal$Fg_se <- rowSds(as.matrix(dist_internal[,
                2:(length(protein_file) + 1)]))/sqrt(length(protein_file))
        dist_internal$Bkg <- rowMeans(dist_internal[,
                (length(protein_file) + 2):((length(protein_file) * 2) + 1)])
        dist_internal$Bkg_se <- rowSds(as.matrix(dist_internal[,
                (length(protein_file) + 2):((length(protein_file) *
                2) + 1)]))/sqrt(length(protein_file))
        dist_exterior$Fg <- rowMeans(dist_exterior[, 2:(length(protein_file) +
                                                                            1)])
        dist_exterior$Fg_se <- rowSds(as.matrix(dist_exterior[,
                2:(length(protein_file) + 1)]))/sqrt(length(protein_file))
        dist_exterior$Bkg <- rowMeans(dist_exterior[,
                (length(protein_file) + 2):((length(protein_file) * 2) + 1)])
        dist_exterior$Bkg_se <- rowSds(as.matrix(dist_exterior[,
                (length(protein_file) + 2):((length(protein_file) *
                2) + 1)]))/sqrt(length(protein_file))
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
            max_y <- max(c(max(dist_exterior$Fg), max(dist_internal$Fg),
                            max(dist_stem$Fg), max(dist_hairpin$Fg),
                            max(dist_multibranch$Fg), max(dist_bulge$Fg))) + 0.1
            min_y <- min(c(min(dist_exterior$Fg), min(dist_internal$Fg),
                            min(dist_stem$Fg), min(dist_hairpin$Fg),
                            min(dist_multibranch$Fg), min(dist_bulge$Fg))) - 0.1
            y_lim <- c(min_y, max_y)
        }
        plot(dist_bulge$x, dist_bulge$Bkg, type = "l",
            col = "black", xlim = x_lim, ylim = y_lim,
            main = NULL, xlab = "Distance", ylab = "Density x 100",
            cex.axis = 0.8, cex.lab = 1, cex.main = 1.2, lwd = 2)
        lines(dist_multibranch$x, dist_multibranch$Bkg, col = "black", lwd = 2)
        lines(dist_stem$x, dist_stem$Bkg, col = "black", lwd = 2)
        lines(dist_hairpin$x, dist_hairpin$Bkg, col = "black", lwd = 2)
        lines(dist_internal$x, dist_internal$Bkg, col = "black", lwd = 2)
        lines(dist_exterior$x, dist_exterior$Bkg, col = "black", lwd = 2)
        arrows(dist_bulge$x, dist_bulge$Bkg - dist_bulge$Bkg_se,
                dist_bulge$x, dist_bulge$Bkg + dist_bulge$Bkg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        arrows(dist_multibranch$x, dist_multibranch$Bkg -
                dist_multibranch$Bkg_se, dist_multibranch$x,
                dist_multibranch$Bkg + dist_multibranch$Bkg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        arrows(dist_stem$x, dist_stem$Bkg - dist_stem$Bkg_se,
                dist_stem$x, dist_stem$Bkg + dist_stem$Bkg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        arrows(dist_hairpin$x, dist_hairpin$Bkg - dist_hairpin$Bkg_se,
                dist_hairpin$x, dist_hairpin$Bkg + dist_hairpin$Bkg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        arrows(dist_internal$x, dist_internal$Bkg -
                dist_internal$Bkg_se, dist_internal$x,
                dist_internal$Bkg + dist_internal$Bkg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        arrows(dist_exterior$x, dist_exterior$Bkg -
                dist_exterior$Bkg_se, dist_exterior$x,
                dist_exterior$Bkg + dist_exterior$Bkg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 0, alpha = 0.5))
        abline(v = 0, col = "grey", lty = 2)
        lines(dist_bulge$x, dist_bulge$Fg, col = "blue", lwd = 2)
        lines(dist_multibranch$x, dist_multibranch$Fg, col = "red", lwd = 2)
        lines(dist_stem$x, dist_stem$Fg, col = "green", lwd = 2)
        lines(dist_hairpin$x, dist_hairpin$Fg, col = "purple", lwd = 2)
        lines(dist_internal$x, dist_internal$Fg, col = "orange", lwd = 2)
        lines(dist_exterior$x, dist_exterior$Fg, col = "cyan", lwd = 2)
        arrows(dist_bulge$x, dist_bulge$Fg - dist_bulge$Fg_se,
                dist_bulge$x, dist_bulge$Fg + dist_bulge$Fg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 1, red = 0, green = 0, alpha = 0.5))
        arrows(dist_multibranch$x, dist_multibranch$Fg -
                dist_multibranch$Fg_se, dist_multibranch$x,
                dist_multibranch$Fg + dist_multibranch$Fg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 1, green = 0, alpha = 0.5))
        arrows(dist_stem$x, dist_stem$Fg - dist_stem$Fg_se,
                dist_stem$x, dist_stem$Fg + dist_stem$Fg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 0, red = 0, green = 1, alpha = 0.5))
        arrows(dist_hairpin$x, dist_hairpin$Fg - dist_hairpin$Fg_se,
                dist_hairpin$x, dist_hairpin$Fg + dist_hairpin$Fg_se,
                length = 0, angle = 90, code = 3,
                col = rgb(blue = 1, red = 1, green = 0, alpha = 0.5))
        arrows(dist_internal$x, dist_internal$Fg -
                dist_internal$Fg_se, dist_internal$x, dist_internal$Fg +
                dist_internal$Fg_se, length = 0, angle = 90,
                code = 3, col = rgb(blue = 0, red = 1, green = 1, alpha = 0.5))
        arrows(dist_exterior$x, dist_exterior$Fg -
                dist_exterior$Fg_se, dist_exterior$x, dist_exterior$Fg +
                dist_exterior$Fg_se, length = 0, angle = 90,
                code = 3, col = rgb(blue = 1, red = 0, green = 1, alpha = 0.5))
        if (legend == TRUE) {
            legend("bottomright",
                legend = c("Background", "Bulge", "Multibranch", "Stem",
                            "Hairpin", "Internal", "Exterior"),
                col = c("black", "blue", "red", "green", "purple",
                        "orange", "cyan"), lty = 1, cex = 0.8)
        }
        dev.off()
    } else {
        out_file <- paste0(out_file, ".jpeg")
        coords_1 <- which(dist_stem$x == x_lim[1])
        coords_2 <- which(dist_stem$x == x_lim[2])
        heatmap_plot <- data.frame(stem = dist_stem$Fg[coords_1:coords_2],
                        hairpin = dist_hairpin$Fg[coords_1:coords_2],
                        bulge = dist_bulge$Fg[coords_1:coords_2],
                        internal = dist_internal$Fg[coords_1:coords_2],
                        multibranch = dist_multibranch$Fg[coords_1:coords_2],
                        exterior = dist_exterior$Fg[coords_1:coords_2])
        heatmap_plot <- as.matrix(t(heatmap_plot))
        zero <- (ncol(heatmap_plot) + 1)/2
        if (is.null(y_lim)) {
            max_y <- max(c(max(dist_exterior$Fg), max(dist_internal$Fg),
                            max(dist_stem$Fg), max(dist_hairpin$Fg),
                            max(dist_multibranch$Fg), max(dist_bulge$Fg))) + 0.1
            min_y <- min(c(min(dist_exterior$Fg), min(dist_internal$Fg),
                            min(dist_stem$Fg), min(dist_hairpin$Fg),
                            min(dist_multibranch$Fg), min(dist_bulge$Fg))) - 0.1
            y_lim <- c(-max(abs(max_y), abs(min_y)),
                        max(abs(max_y), abs(min_y)))
        }
        key <- TRUE
        if (legend == FALSE) {key <- FALSE}
        # make plot
        jpeg(out_file, height = 4, width = 8, units = "in", res = 300)
        heatmap.2(heatmap_plot, Rowv = NA, Colv = NA,
                    col = redblue(256), dendrogram = "none",
                    labCol = FALSE, trace = "none", density.info = "none",
                    symkey = FALSE, key = key, colsep = c(zero - 1, zero),
                    sepcolor = "grey", margins = c(1, 8),
                    breaks = seq(y_lim[1], y_lim[2], length.out = 257))
        dev.off()
    }
}
