# Testing of the calculation of MW values between qSIP package and manual
# calculation code: 110

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
mdq <- calc_mw(mdq, separate_label=T, separate_light=F)

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

# average heavy WADs, remove replicate data
mwl <- aggregate(wad ~ OTU + group, mwl, mean, na.rm=T)

# merge light and heavy, disregarding treatment
mw <- merge(mwh[,!names(mwh) %in% 'tmt'], mwl,
            by=c('OTU', 'group'),
            all.x=T,
            suffixes=c('_label', '_light'))

# calculate MWs
mw <- within(mw, {
  # MW-light
  gc <- (1 / 0.083506) * (wad_light - 1.646057)
  mw_light <- (0.496 * gc) + 307.691
  gc <- NULL
  # MW-label
  mw_label <- (((wad_label - wad_light) / wad_light) + 1) * mw_light
})

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

mw_110 <- merge(mw_qsip, mw_manual,
                by=c('OTU', 'RepID', 'group', 'tmt'),
                suffixes=c('_qsip', '_manual'))
