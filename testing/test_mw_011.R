# Testing of the calculation of MW values between qSIP package and manual
# calculation code: 011

############################################
# qsip package
############################################

mdq <- specify_qsip(md,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    rep_num='rep',
                    iso='18O',
                    iso_trt='tmt')

# calculate WADs
mdq <- calc_mw(mdq, separate_label=T, separate_light=T)

mwh <- mat_to_df(mdq@qsip[['mw_label']], 'mw_label')
mwl <- data.frame(OTU=names(mdq@qsip[['mw_light']]),
                  mw_light=mdq@qsip[['mw_light']],
                  stringsAsFactors=F)

# convert to data frame, use c() to remove attributes
mw <- merge(mwh, mwl, by='OTU', all.x=T)

mw <- reshape(mw,
              idvar=c('RepID', 'OTU'),
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

# separate light and heavy WADs
mwh <- mw[mw$tmt=='18O',]
mwl <- mw[mw$tmt=='16O',]

# average light WADs
mwl <- aggregate(wad ~ OTU, mwl, mean)

# merge back in, disregarding treatment
mw <- merge(mwh[,!names(mwh) %in% 'tmt'], mwl, by='OTU', all.y=T, suffixes=c('_label', '_light'))

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
              idvar=c('RepID', 'OTU'),
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

mw_010 <- merge(mw_qsip, mw_manual,
                by=c('OTU', 'RepID', 'tmt'),
                suffixes=c('_qsip', '_manual'))
