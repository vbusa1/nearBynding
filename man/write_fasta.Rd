\name{write_fasta}
\alias{write_fasta}
\title{
Write FASTA from sequences}
\description{
Writes a FASTA file from a vector of sequences and corresponding names}
\usage{
write_fasta(sequences,
            names,
            file.out)
}
\arguments{
  \item{sequences}{A vector of sequences to be written into a FASTA file}
  \item{names}{A vector of names corresponding to the sequences}
  \item{file.out}{Name of output FASTA file}
}
\value{
writes FASTA file into directory
}
\examples{
sequences<-c(paste0(sample(c("A", "T", "G", "C"), 20, replace = TRUE),
                    collapse = ""),
             paste0(sample(c("A", "T", "G", "C"), 20, replace = TRUE),
                    collapse = ""),
             paste0(sample(c("A", "T", "G", "C"), 20, replace = TRUE),
                    collapse = ""))
write_fasta(sequences,
            c("one", "two", "three"),
            "test.fa")
}
