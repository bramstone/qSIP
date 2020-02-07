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
mdq <- calc_mw(mdq, separate_label=F, separate_light=T)

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

mw$group <- sub('^(\\w{1})(.*)', '\\1', mw$RepID)

mw_qsip <- mw



############################################
# manual calculation
############################################

# add treatment (isotope) and group data for grouping
mw <- merge(wads_man, unique(mdl[,c('RepID', 'tmt', 'group')]), all.x=T)

# separate light and heavy WADs
mwh <- mw[mw$tmt=='18O',]
mwl <- mw[mw$tmt=='16O',]


# split by group
mwh <- split(mwh, mwh$group)
mwl <- split(mwl, mwl$group)

# average heavy WADs, remove replicate data
mwh <- lapply(mwh, function(x) aggregate(wad ~ OTU, x, mean))

# merge light and heavy, disregarding treatment
mw <- Map(function(x,y) merge(x, y[,!names(y) %in% 'tmt'],
                              by='OTU',
                              all.y=T,
                              suffixes=c('_label', '_light')),
          mwh, mwl)

# calculate MWs
mw <- lapply(mw, function(x) {
  # MW-light
  x$gc <- (1 / 0.083506) * (x$wad_light - 1.646057)
  x$mw_light <- (0.496 * x$gc) + 307.691
  x$gc <- NULL
  # MW-label
  x$mw_label <- (((x$wad_label - x$wad_light) / x$wad_light) + 1) * x$mw_light
  x
})

# recombine
mw <- Map(function(x,y) {x$group <- y; x}, mw, names(mw))
mw <- do.call(rbind, mw)


mw <- reshape(mw,
              idvar=c('RepID', 'OTU', 'group'),
              varying=list(c('mw_label', 'mw_light')),
              drop=c('wad_label', 'wad_light'),
              times=c('18O', '16O'),
              timevar='tmt',
              v.names='mw',
              direction='long')
rownames(mw) <- NULL

mw_manual <- mw


############################################
# compare calculations
############################################

mw_101 <- merge(mw_qsip, mw_manual[,!names(mw_manual) %in% 'wad'],
                by=c('OTU', 'group', 'tmt'),
                suffixes=c('_qsip', '_manual'))
