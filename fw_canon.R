#!/usr/bin/env Rscript

# Max Brown 2021; Wellcome Sanger Institute

###############
## Libraries ##
###############

library(argparse)
library(data.table)

######################
## Argument Parsing ##
######################

parser <- ArgumentParser()

parser$add_argument("-t", "--tsv", help = "input TSV file from fasta_windows output.")

args <- parser$parse_args()

###############
## Functions ##
###############

# get the kmers from the colnames of out data table
get_kmers <- function(x) {
  colnms <- colnames(x)
  matches <- regmatches(colnms, gregexpr("A*|G*|C*|T*", colnms))
  matches_list <- lapply(matches, function(e) paste0(e, collapse = ""))

  kmers <- unlist(matches_list)
  kmers[kmers != ""]
}

# reverse complement a (vector of)
# string of DNA letters
revcomp <- function(x) {
  # kmer assignment vector
  kmers <- character(length(x))

  # loop over vector of kmers
  for (i in seq_len(length(x))) {
    # bases assignement vector
    vector <- character(length(kmers[1]))
    bases <- strsplit(x[i], "")[[1]]

    for (j in seq_len(length(bases))) {
      # complement
      switch(bases[j],
        "A" = {
          bases[j] <- "T"
        },
        "G" = {
          bases[j] <- "C"
        },
        "C" = {
          bases[j] <- "G"
        },
        "T" = {
          bases[j] <- "A"
        }
      )
      vector[j] <- bases[j]
    }
    # and reverse
    kmers[i] <- paste0(rev(vector), collapse = "")
  }
  return(kmers)
}

# from a full set of non-canonical kmers
# generate the set of canonical kmers
generate_canonical_kmers <- function(x) {
  kmers <- x
  rev_comp_kmers <- revcomp(x)

  rev_comp_kmers[!(kmers < rev_comp_kmers)]
}

# aggregate columns of our data with the same name
agg <- function(df) {
  do.call(cbind, lapply(split(as.list(df), names(df)), function(x) {
    Reduce(`+`, x)
  }))
}

##########
## Main ##
##########

File <- args$tsv

if (is.null(File)) {
  stop("[-]\tNo input file detected.")
}

dat <- fread(File)

# get kmers from our input data
kmers <- get_kmers(dat)

# rename columns
colnames(dat)[
  !grepl(pattern = "ID|description|start|end", x = colnames(dat))
] <-
  ifelse(kmers < revcomp(kmers), kmers, revcomp(kmers))

# aggregate data on equivalent colnames
# and set as DT
dat2 <- agg(dat)
dat2 <- setDT(as.data.frame(dat2))

# generate canonical kmers and lexicographically order them
canon_kmers <- generate_canonical_kmers(kmers)
canon_kmers_ordered <- canon_kmers[order(canon_kmers)]

# apply this to output
setcolorder(dat2, c(colnames(dat)[
  grepl(pattern = "ID|description|start|end", x = colnames(dat))
], canon_kmers_ordered))

# write to stdout
fwrite(x = dat2, file = "", sep = "\t")