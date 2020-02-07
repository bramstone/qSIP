# Testing of the calculation of MW values between qSIP package and manual
# calculation code: 100

############################################
# qsip package
############################################

mdq <- specify_qsip(md,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    rep_group='group',
                    iso='18O',
                    iso_trt='tmt')

# calculate WADs
mdq <- calc_mw(mdq, separate_label=F, separate_light=F)

mwh <- mat_to_df(mdq@qsip[['mw_label']], 'mw_label')
mwl <- mat_to_df(mdq@qsip[['mw_light']], 'mw_light')

mw <- merge(mwh, mwl, by=c('RepID', 'OTU'), all.x=T)

mw <- reshape(mw,
              idvar=c('RepID', 'OTU'),
              varying=list(c('mw_label', 'mw_light')),
              times=c('18O', '16O'),
              timevar='tmt',
              v.names='mw',
              direction='long')
rownames(mw) <- NULL

names(mw)[1] <- 'group'
mw$group <- sub('(.*)\\.(.*)', '\\2', mw$group)

mw_qsip <- mw



############################################
# manual calculation
############################################

# add treatment (isotope) and group data for grouping
mw <- merge(wads_man, unique(mdl[,c('RepID', 'tmt', 'group')]), all.x=T)

# split by treatment, group, and taxon
mw <- split(mw, interaction(mw$tmt, mw$group, mw$OTU))

# average WADs, remove replicate data
mw <- lapply(mw, function(x) {
  x$wad <-  mean(x$wad, na.rm=T)
  x$RepID <- NULL
  x[1,]
})

# recombine, then split by group and treatment
mw <- do.call(rbind, mw)
mw <- split(mw, interaction(mw$tmt, mw$group))

# separate light and heavy MWs
mwl <- mw[grep('16O', names(mw))]
mwh <- mw[grep('18O', names(mw))]

# calculate MW-light
mwl <- lapply(mwl, function(x) {
  gc <- (1 / 0.083506) * (x$wad - 1.646057)
  x$mw <- (0.496 * gc) + 307.691
  gc <- NULL
  x
})

# calculate MW-heavy across groups
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

mw_100 <- merge(mw_qsip, mw_manual[,!names(mw_manual) %in% 'wad'],
                by=c('OTU', 'group', 'tmt'),
                suffixes=c('_qsip', '_manual'))
