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
parser$add_argument("-w", "--window",
  type = "integer", default = 1000,
  help = "Size of window used in fasta_windows [fixed at 1000]",
  metavar = "number"
)
parser$add_argument("-g", "--group",
  type = "integer", default = 10000,
  help = "Size of window to group into [default 10000]",
  metavar = "number"
)
parser$add_argument("-r", "--round",
  type = "integer", default = 3,
  help = "Round numbers to this many decimal places. Must be < 3 [default 3]",
  metavar = "number"
)
parser$add_argument("-f", "--function",
  type = "character", default = "mean",
  help = "Function to apply to groups [default mean: other options are Mode, median, var, sd]",
  metavar = "string"
)
parser$add_argument("-c", "--chromosomal",
  action = "store_true",
  help = "Aggregate statistics at the chromosomal level?"
)
args <- parser$parse_args()

File <- args$tsv
window_size <- args$window
group_size <- args$group
fun <- args$`function`
chromosomal <- args$chromosomal
r <- args$round

# remove scientific notation
options(scipen = 999)

# some parsing logic

if (is.null(File)) {
  stop("[-]\tNo input file detected.")
}

# keep just for future dev potentially
if (window_size != 1000) {
  stop("[-]\tWindow size should be 1000bp (1kb).")
}

# could be more flexible..?
group_size_options <- c(10000, 100000, 1000000)

if (group_size %in% group_size_options == FALSE) {
  stop("[-]\tGroup should be either 10000 (10kb), 100000 (100kb), or 1000000 (1Mb).")
}

fun_options <- c("mean", "Mode", "median", "var", "sd")

if (fun %in% fun_options == FALSE) {
  stop("[-]\tGrouping function should be either mean, Mode, median, var, or sd.")
}

if (r > 3) {
  stop("[-]\tr should be < 3, as the original output is capped at 3 decimal places.")
}

# dinky mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

###############
## Functions ##
###############

calc_aggregations <- function(file, window_size, group_size, fun = c("mean", "Mode", "median", "var", "sd"), chromosomal) {
  # read in file
  windows <- fread(file = file)
  # get aggregation function
  fun <- match.fun(fun)
  # fun_char <- as.character(substitute(fun))
  # column names will be consistent, so subset is okay here.
  column_names <- names(windows)[4:14]

  if (chromosomal) {
    chromosomal_table <- windows[, .(
      fun(get(column_names[1])),
      fun(get(column_names[2])),
      fun(get(column_names[3])),
      fun(get(column_names[4])),
      fun(get(column_names[5])),
      fun(get(column_names[6])),
      fun(get(column_names[7])),
      fun(get(column_names[8])),
      fun(get(column_names[9])),
      fun(get(column_names[10])),
      fun(get(column_names[11]))
    ),
    by = .(ID)
    ]

    setnames(chromosomal_table, 2:12, column_names)
    # do some rounding
    rounding <- names(chromosomal_table)[2:12]
    chromosomal_table[, (rounding) := round(.SD, r), .SDcols = rounding]

    return(chromosomal_table)
  }
  # add column of group size
  windows[, group := start %/% group_size]

  # aggregate using function
  grouped <- windows[, .(
    fun(get(column_names[1])),
    fun(get(column_names[2])),
    fun(get(column_names[3])),
    fun(get(column_names[4])),
    fun(get(column_names[5])),
    fun(get(column_names[6])),
    fun(get(column_names[7])),
    fun(get(column_names[8])),
    fun(get(column_names[9])),
    fun(get(column_names[10])),
    fun(get(column_names[11]))
  ),
  by = .(ID, group)
  ]

  # make bed like columns
  grouped[, `:=`(start = group * group_size, end = (group + 1) * group_size, group = NULL)]

  setcolorder(grouped, c(
    "ID",
    "start",
    "end",
    "V1",
    "V2",
    "V3",
    "V4",
    "V5",
    "V6",
    "V7",
    "V8",
    "V9",
    "V10",
    "V11"
  ))
  setnames(grouped, 4:14, column_names)

  rounding <- names(grouped)[4:14]
  grouped[, (rounding) := round(.SD, r), .SDcols = rounding]

  return(grouped)
}

##########
## Main ##
##########

res <- calc_aggregations(
  file = File,
  window_size = window_size,
  group_size = group_size,
  fun = fun,
  chromosomal = chromosomal
)

fwrite(x = res, file = "", sep = "\t")