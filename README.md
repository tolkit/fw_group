# fw_group

A tiny R script which takes the output of fasta_windows (<a href="https://github.com/tolkit/fasta_windows/tree/v2">version 2</a>) and aggregates the 1kb records. A small variety of grouping functions are supported.

## Usage

Run ` Rscript fw_group.R -h` to display the help.

```
usage: fw_group.R [-h] [-t TSV] [-w number] [-g number] [-f string] [-c]

optional arguments:
  -h, --help            show this help message and exit
  -t TSV, --tsv TSV     input TSV file from fasta_windows output.
  -w number, --window number
                        Size of window used in fasta_windows [fixed at 1000]
  -g number, --group number
                        Size of window to group into [default 10000]
  -f string, --function string
                        Function to apply to groups [default mean: other
                        options are Mode, median, var, sd]
  -c, --chromosomal     Aggregate statistics at the chromosomal level?
```

Currently only fasta windows output with 1kb windows is supported.

## Examples

`Rscript fw_group.R -t <TSV>`

Running on fasta windows output (1kb windows) without any flags will average statistics across 10kb windows.

`Rscript fw_group.R -t <TSV> -g 100000`

Will average statistics across 100kb windows. 1Mb (1000000) averaging also supported.

`Rscript fw_group.R -t <TSV> -f var`

Calculates variance of the statistics, instead of default mean.

`Rscript fw_group.R -t <TSV> -c`

Adding `-c` or `--chromosomal` averages across the contigs/scaffolds/chromsomes.

Note that these statistics will be different from running fasta windows at a larger window size.

## Requirements

- <a href="https://github.com/Rdatatable/data.table">data.table</a>
- <a href="https://github.com/trevorld/r-argparse>">argparse</a>

## TODO's

- Sort out decimal rounding.
