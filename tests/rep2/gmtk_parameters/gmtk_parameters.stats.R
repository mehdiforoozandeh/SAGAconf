#!/usr/bin/env Rscript
## transcript produced by Segtools 1.2.4

segtools.r.dirname <-
            system2("python",
            c("-c", "'import segtools; print segtools.get_r_dirname()'"),
            stdout = TRUE)

source(file.path(segtools.r.dirname, 'common.R'))
source(file.path(segtools.r.dirname, 'signal.R'))
source(file.path(segtools.r.dirname, 'track_statistics.R'))
save.track.stats('gmtk_parameters', 'gmtk_parameters.stats', 'params.params', mnemonic_file = '', translation_file = '', as_regex = FALSE, gmtk = TRUE, clobber = FALSE, label_order = list(), track_order = list())
