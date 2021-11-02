#!/usr/bin/env Rscript
## transcript produced by Segtools 1.2.4

## Experimental R transcript
## You may not be able to run the R code in this file exactly as written.

segtools.r.dirname <-
            system2("python",
            c("-c", "'import segtools; print segtools.get_r_dirname()'"),
            stdout = TRUE)

source(file.path(segtools.r.dirname, 'common.R'))
source(file.path(segtools.r.dirname, 'length.R'))
save.length('length_distribution', 'length_distribution', 'length_distribution/length_distribution.tab', mnemonic_file = '', clobber = FALSE)
save.segment.sizes('length_distribution', 'segment_sizes', 'length_distribution/segment_sizes.tab', mnemonic_file = '', clobber = FALSE, show_segments = TRUE, show_bases = TRUE)
