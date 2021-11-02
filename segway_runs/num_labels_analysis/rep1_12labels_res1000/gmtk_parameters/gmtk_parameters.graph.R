#!/usr/bin/env Rscript
## transcript produced by Segtools 1.2.4

## Experimental R transcript
## You may not be able to run the R code in this file exactly as written.

segtools.r.dirname <-
            system2("python",
            c("-c", "'import segtools; print segtools.get_r_dirname()'"),
            stdout = TRUE)

source(file.path(segtools.r.dirname, 'common.R'))
source(file.path(segtools.r.dirname, 'transition.R'))
save.transition('gmtk_parameters', 'gmtk_parameters', 'params.params', clobber = FALSE, mnemonic_file = '', ddgram = FALSE, gmtk = TRUE)
