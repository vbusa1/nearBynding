.bedtools_run <- function(cmd) {
    system2("bedtools", cmd)
}
.is_bedtools_installed <- function() {
    code <- suppressWarnings(
        system2("bedtools", "--version", stdout = NULL, stderr = NULL)
    )
    code == 0
}

.sort_run <- function(cmd) {
    system2("sort", cmd)
}
.is_sort_installed <- function() {
    code <- suppressWarnings(
        system2("sort", "--version", stdout = NULL, stderr = NULL)
    )
    code == 0
}

.is_CapR_installed <- function() {
    code <- suppressWarnings(
        system2("CapR", stdout = NULL, stderr = NULL)
    )
    code == 0
}
.CapR_run <- function(in_file, out_file, max_dist) {
    system2("CapR", args = c(in_file, out_file, max_dist))
}

.StereoGene_run <- function(cmd) {
    system2("StereoGene", cmd)
}
.is_StereoGene_installed <- function() {
    code <- suppressWarnings(
        system2("StereoGene", "-h", stdout = NULL, stderr = NULL)
    )
    code == 0
}
