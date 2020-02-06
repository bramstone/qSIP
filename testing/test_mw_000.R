# Testing of the calculation of MW values between qSIP package and manual
# calculation code: 000

############################################
# qsip package
############################################

mdq <- specify_qsip(md,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    iso='18O',
                    iso_trt='tmt')

# calculate WADs
mdq <- calc_mw(mdq, separate_label=F, separate_light=F)

# convert to data frame, use c() to remove attributes
mw <- data.frame(OTU=names(mdq@qsip[['mw_label']]),
                 mw_label=c(mdq@qsip[['mw_label']]),
                 mw_light=c(mdq@qsip[['mw_light']]),
                 stringsAsFactors=F)

mw <- reshape(mw,
              idvar='OTU',
              varying=list(c('mw_label', 'mw_light')),
              times=c('18O', '16O'),
              timevar='tmt',
              v.names='mw',
              direction='long')
rownames(mw) <- NULL

mw_qsip <- mw



############################################
# manual calculation
############################################

# add treatment (isotope) data for grouping
mw <- merge(wads_man, unique(mdl[,c('RepID', 'tmt')]), all.x=T)

# split by treatment group and taxon
mw <- split(mw, interaction(mw$tmt, mw$OTU))

# average WADs, remove replicate data
mw <- lapply(mw, function(x) {
  x$wad <-  mean(x$wad, na.rm=T)
  x$RepID <- NULL
  x[1,]
})

# separate light and heavy MWs
mwl <- mw[grep('16O', names(mw))]
mwh <- mw[grep('18O', names(mw))]

# calculate MW-light
mwl <- lapply(mwl, function(x) {
  gc <- (1 / 0.083506) * (x$wad - 1.646057)
  x$mw <- (0.496 * gc) + 307.691
  x
})

# calculate MW-heavy
mwh <- Map(function(x,y) {
  x$mw <- (((x$wad - y$wad)/y$wad) + 1) * y$mw
  x
}, mwh, mwl)

# combine
mw <- rbind(do.call(rbind, mwl),
            do.call(rbind, mwh))
rownames(mw) <- NULL
mw_manual <- mw


############################################
# compare calculations
############################################

mw_000 <- merge(mw_qsip, mw_manual[,c('OTU', 'tmt', 'mw')],
                by=c('OTU', 'tmt'),
                suffixes=c('_qsip', '_manual'))
