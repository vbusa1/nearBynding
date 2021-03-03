#' @title symmetryCapR
#'
#' @description Calculate the symmetry of a binding context.
#'
#' @param dir_stereogene_output Directory of Stereogene output for first
#' protein. Default current directory.
#' @param CapR_prefix The prefix common to CapR output files of protein_file, if
#' applicable. Equivalent to output_prefix from runStereogeneOnCapR. Default ""
#' @param protein_file A vector of strings with at least one protein file name
#' to be averaged for calculation of distance. File names must exclude
#' extensions such as ".bedGraph". All files in the list should be
#' experimental/biological replicates. Required.
#' @param protein_file_input A protein file name of background input to be
#' subtracted from protein_file signal. File name must exclude extension. Only
#' one input file is permitted. Optional.
#' @param context The RNA structure context being interrogated. Acceptable
#' contexts include "all", which sums the distance of all six contexts, or any
#' of the contexts individually ("bulge", "hairpin", "stem", "exterior",
#' "multibranch", or "internal"). Default "all"
#' @param range A vector of two integers denoting the range upstream and
#' downstream of the center of protein binding to consider in the comparison.
#' Ranges that are too small miss the holistic binding context, while large
#' ranges amplify distal noise in the binding data. Cannot exceed wSize/2 from
#' write_config. Default c(-200, 200)
#'
#' @return Wasserstein distance between the two halves of the binding context,
#' with lower values suggesting greater symmetry.
#'
#'
#' @examples
#' ## load example StereoGene output
#' get_outfiles()
#'
#' ## This boring example compares a protein's binding with itself for all
#' ## contexts, therefore the distance is 0
#' symmetryCapR(CapR_prefix = "chr4and5_3UTR",
#'                 protein_file = "chr4and5_liftOver")
#'
#' @importFrom utils read.table
#' @importFrom matrixStats rowSds
#' @importFrom magrittr '%>%'
#' @importFrom dplyr filter
#' @importFrom transport wasserstein1d
#'
#' @export

symmetryCapR <- function(dir_stereogene_output = ".",
                                       CapR_prefix = "",
                                       protein_file,
                                       protein_file_input = NULL,
                                       context = "all",
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
    if (!(context %in% c("all", "bulge", "hairpin",
                         "stem", "exterior", "multibranch", "internal"))) {
        stop("This context is not available for CapR output analysis. Available
        contexts are \"all\", \"bulge\", \"hairpin\", \"stem\", \"exterior\",
        \"multibranch\",  or \"internal\"")
    }
    get_dist <- NULL
    dist_1 <- NULL
    second_dist_1 <- NULL
    get_dist <- function(context) {
        for (n in seq(length(protein_file))) {
            assign(paste0("dist_", n), read.table(paste0(dir_stereogene_output,
                                                         "/", CapR_prefix, "_", context, "_liftOver~",
                                                         protein_file[n], ".dist"), header = TRUE) %>%
                dplyr::filter(range[1] <= .data$x, .data$x <= 0))
            assign(paste0("second_dist_", n), read.table(paste0(dir_stereogene_output,
                                                         "/", CapR_prefix, "_", context, "_liftOver~",
                                                         protein_file[n], ".dist"), header = TRUE) %>%
                dplyr::filter(0 <= .data$x, .data$x <= range[2]))
        }
        if (!is.null(protein_file_input)) {
            input <- read.table(paste0(dir_stereogene_output,
                                            "/", CapR_prefix, "_", context, "_liftOver~",
                                            protein_file_input, ".dist"), header = TRUE)
            dist_input<- input%>%
                dplyr::filter(range[1] <= .data$x, .data$x <= 0)
            second_dist_input <- input %>%
                dplyr::filter(0 <= .data$x, .data$x <= range[2])
        }
        dist <- as.data.frame(matrix(NA,
                                     ncol = (length(protein_file)) + 1,
                                     nrow = nrow(dist_1)))
        colnames(dist) <- c("x", paste0("Fg", seq(length(protein_file))))
        dist$x <- dist_1$x
        second_dist <- as.data.frame(matrix(NA,
                                            ncol = (length(protein_file)) + 1,
                                            nrow = nrow(second_dist_1)))
        colnames(second_dist) <- c("x", paste0("Fg",
                                               seq(length(protein_file))))
        second_dist$x <- second_dist_1$x
        for (n in seq(length(protein_file))) {
            dist[, 1 + n] <- eval(parse(text = paste0("dist_", n)))$Fg
            second_dist[, 1 + n] <- eval(parse(text = paste0("second_dist_", n)))$Fg
        }
        if (!is.null(protein_file_input)) {
            dist[, 2:(length(protein_file) + 1)] <- dist[,
                                                         2:(length(protein_file) + 1)] - dist_input$Fg
            second_dist[, 2:(length(protein_file) + 1)] <- second_dist[,
                                                         2:(length(protein_file) + 1)] - second_dist_input$Fg
        }
        if (length(protein_file) > 1) {
            dist$Fg <- rowMeans(dist[, 2:(length(protein_file) + 1)])
            second_dist$Fg <- rowMeans(second_dist[, 2:(length(protein_file) + 1)])
        } else {
            dist$Fg <- dist[, 2]
        }
        second_dist$x<- -second_dist$x

        # scale tracks
        max<-max(abs(dist$Fg), abs(second_dist$Fg))
        dist$Fg<-dist$Fg/max
        second_dist$Fg<-second_dist$Fg/max

        wasserstein_distance <- suppressWarnings(wasserstein1d(dist$Fg,
                                                               second_dist$Fg) %>% as.numeric())
        return(wasserstein_distance)
    }
    if (context == "all") {
        distance <- sum(get_dist("bulge"), get_dist("hairpin"),
                        get_dist("internal"), get_dist("exterior"),
                        get_dist("stem"), get_dist("multibranch"))
    } else {
        distance <- get_dist(context)
    }
    return(distance)
}
