#' @title bindingContextDistance
#'
#' @description Calculate the Wasserstein distance between two replicates' or
#' two proteins' binding contexts for CapR-generated RNA contexts.
#'
#' @param dir_stereogene_output Directory of Stereogene output for first
#' protein. Default current directory.
#' @param RNA_context Name of the RNA context file input to Stereogene. File
#' names must exclude extensions such as ".bedGraph". Requred
#' @param protein_file A vector of at least one protein file name to be
#' averaged for calculation of distance. File names must exclude extensions
#' such as ".bedGraph". All files in the list should be experimental/biological
#' replicates. Required.
#' @param protein_file_input A protein file name of background input to be
#' subtracted from protein_file signal. File name must exclude extension. Only
#' one input file is permitted. Optional.
#' @param dir_stereogene_output_2 Directory of Stereogene output for second
#' protein. Default dir_stereogene_output.
#' @param RNA_context_2 Name of the RNA context file input to Stereogene. File
#' names must exclude extensions such as ".bedGraph". Default same as
#' RNA_context.
#' @param protein_file_2 Similar to protein_file. A second vector of at least
#' one protein file name to be averaged for calculation of distance. File names
#' must exclude extensions such as ".bedGraph". All files in the list should be
#' experimental/biological replicates. Default same as protein_file
#' @param protein_file_input_2 Similar to protein_file_input. A second protein
#' file name of background input to be subtracted from protein_file_2 signal.
#' File name must exclude extension. Only one input file is permitted. Optional.
#' @param range A vector of two integers denoting the range upstream and
#' downstream of the center of protein binding to consider in the comparison.
#' Ranges that are too small miss the holistic binding context, while large
#' ranges amplify distal noise in the binding data. Cannot exceed wSize/2 from
#' write_config. Default c(-200, 200)
#'
#' @return Wasserstein distance between the two protein file sets provided for
#' the RNA structure context specified, minus the input binding signal if
#' applicable
#'
#' @note Either RNA_context_2 or protein_file_2 must be input. Otherwise, the
#' distance would be calculated between the same file and equal 0.
#' @note Wasserstein distance calculations are reciprocal, so it does not matter
#' which protein is first or second so long as replicates and input files
#' correspond to one another.
#'
#' @examples
#' ## pull example files
#' get_outfiles()
#' ## distance between stem and hairpin contexts
#' bindingContextDistance(RNA_context = "chr4and5_3UTR_stem_liftOver",
#'                        protein_file = "chr4and5_liftOver",
#'                        RNA_context_2 = "chr4and5_3UTR_hairpin_liftOver")
#'
#' ## distance between internal and hairpin contexts
#' bindingContextDistance(RNA_context = "chr4and5_3UTR_internal_liftOver",
#'                        protein_file = "chr4and5_liftOver",
#'                        RNA_context_2 = "chr4and5_3UTR_hairpin_liftOver")
#'
#' @importFrom utils read.table
#' @importFrom matrixStats rowSds
#' @importFrom magrittr '%>%'
#' @importFrom dplyr filter
#' @importFrom transport wasserstein1d
#'
#' @export

bindingContextDistance <- function(dir_stereogene_output = ".",
                                    RNA_context, protein_file,
                                    protein_file_input = NULL,
                                    dir_stereogene_output_2 = NULL,
                                    RNA_context_2 = NULL,
                                    protein_file_2 = NULL,
                                    protein_file_input_2 = NULL,
                                    range = c(-200, 200)) {
    if (length(protein_file) < 1) {
        stop("Requires at least one protein file prefix to calculate distance")
    }
    if (length(protein_file) > 20) {
        stop("There are > 20 protein files input. This is likely in error")
    }
    if (length(protein_file_input) > 1) {
        stop("Only input one background track per protein.")
    }
    if (is.null(protein_file_2) & is.null(RNA_context_2)) {
        stop("Requires either RNA_context_2 or protein_file_2")
    }
    if (length(protein_file_2) > 20) {
        stop("There are > 20 protein files input. This is likely in error")
    }
    if (length(protein_file_input_2) > 1) {
        stop("Only input one background track per protein.")
    }
    if (is.null(dir_stereogene_output_2)) {
        dir_stereogene_output_2 <- dir_stereogene_output
    }
    if (is.null(RNA_context_2)) {
        RNA_context_2 <- RNA_context
    }
    if (is.null(protein_file_2)) {
        protein_file_2 <- protein_file
    }
    dist_1 <- NULL
    second_dist_1 <- NULL
    for (n in seq(length(protein_file))) {
        assign(paste0("dist_", n), read.table(paste0(dir_stereogene_output,
            "/", RNA_context, "~", protein_file[n], ".dist"), header = TRUE) %>%
        dplyr::filter(range[1] <= .data$x, .data$x <= range[2]))
    }
    if (!is.null(protein_file_input)) {
        dist_input <- read.table(paste0(dir_stereogene_output,
            "/", RNA_context, "~", protein_file_input, ".dist"),
            header = TRUE) %>%
        dplyr::filter(range[1] <= .data$x, .data$x <= range[2])
    }
    dist <- as.data.frame(matrix(NA, ncol = (2 * length(protein_file)) + 1,
                                    nrow = nrow(dist_1)))
    colnames(dist) <- c("x", paste0("Fg", seq(length(protein_file))),
                        paste0("Bkg", seq(length(protein_file))))
    dist$x <- dist_1$x
    for (n in seq(length(protein_file))) {
        dist[, 1 + n] <- eval(parse(text = paste0("dist_", n)))$Fg
        dist[, 1 + n + length(protein_file)] <-
            eval(parse(text = paste0("dist_", n)))$Bkg
    }
    if (!is.null(protein_file_input)) {
        dist[, 2:(length(protein_file) + 1)] <-
            dist[, 2:(length(protein_file) + 1)] - dist_input$Fg
        dist[, (length(protein_file) + 1):ncol(dist)] <-
            dist[, (length(protein_file) + 1):ncol(dist)] - dist_input$Bkg
    }
    if (length(protein_file) > 1) {
        dist$Fg <- rowMeans(dist[, 2:(length(protein_file) + 1)])
        dist$Fg_se <- rowSds(as.matrix(dist[, 2:(length(protein_file) +
                        1)]))/sqrt(length(protein_file))
        dist$Bkg <- rowMeans(dist[, (length(protein_file) +1):(ncol(dist) - 2)])
        dist$Bkg_se <- rowSds(as.matrix(dist[, (length(protein_file) +
                        1):(ncol(dist) - 3)]))/sqrt(length(protein_file))
    } else {
        dist$Fg <- dist[, 2]
        dist$Fg_se <- 0
        dist$Bkg <- dist[, 3]
        dist$Bkg_se <- 0
    }
    for (n in seq(length(protein_file_2))) {
        assign(paste0("second_dist_", n), read.table(paste0(
            dir_stereogene_output_2, "/", RNA_context_2, "~", protein_file_2[n],
            ".dist"), header = TRUE) %>%
        dplyr::filter(range[1] <= .data$x, .data$x <= range[2]))
    }
    if (!is.null(protein_file_input_2)) {
        second_dist_input <- read.table(paste0(dir_stereogene_output_2,
        "/", RNA_context_2, "~", protein_file_input_2, ".dist"),
        header = TRUE) %>%
        dplyr::filter(range[1] <= .data$x, .data$x <= range[2])
    }
    second_dist <- as.data.frame(matrix(NA, ncol = (2 *
                    length(protein_file_2)) + 1, nrow = nrow(second_dist_1)))
    colnames(second_dist) <- c("x", paste0("Fg", seq(length(protein_file_2))),
                                paste0("Bkg", seq(length(protein_file_2))))
    second_dist$x <- second_dist_1$x
    for (n in seq(length(protein_file_2))) {
        second_dist[, 1 + n] <- eval(parse(text = paste0("second_dist_", n)))$Fg
        second_dist[, 1 + n + length(protein_file_2)] <-
            eval(parse(text = paste0("second_dist_", n)))$Bkg
    }
    if (!is.null(protein_file_input_2)) {
        second_dist[, 2:(length(protein_file_2) + 1)] <-
            second_dist[, 2:(length(protein_file_2) + 1)] - second_dist_input$Fg
        second_dist[, (length(protein_file_2) + 1):ncol(second_dist)] <-
            second_dist[, (length(protein_file_2) + 1):ncol(second_dist)] -
            second_dist_input$Bkg
    }
    if (length(protein_file) > 1) {
        second_dist$Fg <- rowMeans(second_dist[, 2:(length(protein_file) + 1)])
        second_dist$Fg_se <- rowSds(as.matrix(second_dist[,
                                                2:(length(protein_file) + 1)]))/
            sqrt(length(protein_file))
        second_dist$Bkg <- rowMeans(second_dist[, (length(protein_file) +
                                                    1):(ncol(second_dist) - 2)])
        second_dist$Bkg_se <- rowSds(as.matrix(second_dist[,
                        (length(protein_file) + 1):(ncol(second_dist) - 3)]))/
            sqrt(length(protein_file))
    } else {
        second_dist$Fg <- second_dist[, 2]
        second_dist$Fg_se <- 0
        second_dist$Bkg <- second_dist[, 3]
        second_dist$Bkg_se <- 0
    }
    wasserstein_distance <- suppressWarnings(wasserstein1d(dist$Fg,
                                            second_dist$Fg) %>% as.numeric())
    return(wasserstein_distance)
}
