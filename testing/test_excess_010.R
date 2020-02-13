# calc code 010

mdq <- specify_qsip(md,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    iso='18O',
                    iso_trt='tmt')

# calculate atom excess
mdq <- calc_excess(mdq, separate_label=T, separate_light=F)

# convert to data frame, use c() to remove attributes
ae <- mat_to_df(mdq@qsip[['atom_excess']], 'excess')


ae_qsip <- ae


# manual calculation-------------------------------------------
ae <- mw_010[,!names(mw_010) %in% 'mw_qsip']

# separate MW-light and MW-label
aeh <- ae[ae$tmt=='18O',]
ael <- ae[ae$tmt=='16O',]

# DON'T average light MWs, they've already been calculated based on averages

# calculate MW-heavymax
adjust <- 12.07747
nat_abund <- 0.002010429
ael$mw_max <- (adjust + ael$mw_manual)

# merge, disregarding tmt
ae <- merge(aeh[,!names(aeh) %in% 'tmt'],
            ael[,!names(ael) %in% 'tmt'],
            by=c('RepID', 'OTU'),
            suffixes=c('_label', '_light'))

# calculate atom excess
ae <- within(ae, {
  excess <- ((mw_manual_label - mw_manual_light)/(mw_max - mw_manual_light)) * (1 - nat_abund)
})


# compare calculations-----------------------------------------
ae_010 <- merge(ae_qsip,
                ae[,c('RepID', 'OTU', 'excess')],
                by=c('RepID', 'OTU'),
                suffixes=c('_qsip', '_manual'))
