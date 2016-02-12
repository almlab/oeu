# oeu

## Requirements
You'll need the R packages:
- ggplot2
- dplyr
- optparse
- reshape2

You can install a package by running `install.packages("foo")` in R.

## Installation
These are standalone scripts and shouldn't require installation. You can
run them as if they were shell scripts using [`Rscript`](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html).

## Usage
The two R executables `oeu.R` and `make-plots.R` can call OTUs and make
plots of OTUs in each OEU.

Both scripts use optparse, so you can call, say, `oeu.R --help` to get
reminders about the order of arguments.

### `oeu.R`
This script takes an OTU table that:
- is normalized (i.e., it shows relative abundance in sample), 
- is preprocessed (i.e., replicated have been pooled, bad OTUs/samples have been removed),
- has OTU names in the first column (preferably under "OTU\_ID"), and
- is tab-separated.

This script produces a "groups file", which is a tab-separated file with
two columns: OTU ID and group number.

### `make-plots.R`
This script takes
- an OTU table, the same kind that fed to `oeu.R`, and
- a groups file, the same kind produced by `oeu.R`.

The script will produce a one plot per OEU in a directory you choose.
