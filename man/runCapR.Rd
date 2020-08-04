\name{runCapR}
\alias{runCapR}
\title{
Run CapR}
\description{
Calls CapR on local machine and inputs .fa sequences for folding.}
\usage{
runCapR(in_file,
        out_file = NA,
        max_dist = 100)
}
\arguments{
    \item{in_file}{An .fa file like from ExtractTranscriptomeSequence that is a
    list of fasta sequences to be folded. Required}
    \item{out_file}{Name of the CapR output file of nucleotide folding
    probabilities. Default is in_file prefix.out}
    \item{max_dist}{Maximum distance between folded nucleotides in sequences.
    Recommeded between 50 and 100, with higher values yielding potentially more
    accurate results but dramatically increasing run time. Default 100.}
}
\value{
generates CapR outfile
}
\examples{
## make dummy file
write_fasta(paste0(sample(c("A", "T", "G", "C"), 50, replace = TRUE),
                    collapse = ""),
            "test",
            "test.fa")
## run CapR
\dontrun{
runCapR("test.fa")
}
}
