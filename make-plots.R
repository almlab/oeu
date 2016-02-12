#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))

p <- OptionParser(usage="%prog [options] groups_dat otu_table")
p <- add_option(p, c("-p", "--plot_dir"), default="plots",
  help="directory in which to put all the plots (default: %default)", metavar="DIR")
args <- parse_args(p, positional_arguments=2)

plot.dir <- args$options$plot_dir
groups.fn <- args$args[1]
otu.fn <- args$args[2]

stat <- file.info(plot.dir)
if (is.na(stat) || !stat$isdir) {
  stop(sprintf("plot directory '%s' does not exist. maybe mkdir it?", plot.dir))
}


# read in the groups data and OTU table
group <- read.table(groups.fn, header=T, sep="\t")
otu <- read.table(otu.fn, header=T, sep="\t", row.names=1)

# rescale the rows (i.e., make profiles)
otu <- apply(otu, 1, function(x) x / sum(x)) %>% t %>% as.data.frame
otu$otu <- rownames(otu)

# create a data frame with entries otu, group, sample, value (= ra)
dat <- merge(group, otu) %>% melt(id.vars=c("otu", "group"), variable.name="sample")

# add the depths by parsing the sample names
dat$depth <- sapply(dat$sample %>% as.list, function(x) sub("M", "", x) %>% as.numeric)

# plot each piece of the output
for (g in 1:max(dat$group, na.rm=T)) {
  foo <- filter(dat, group == g)
  p <- dat %>% filter(group == g) %>%
    ggplot(aes(x=depth, y=value)) + geom_line(aes(group=otu)) +
    theme_minimal()
  fn <- sprintf("%s/plot_%02d.pdf", plot.dir, g)
  ggsave(fn, plot=p)
}
