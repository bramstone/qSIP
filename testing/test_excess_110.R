# calc code 110

mdq <- specify_qsip(md,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    rep_group='group',
                    iso='18O',
                    iso_trt='tmt')

# calculate atom excess
mdq <- calc_excess(mdq, separate_label=T, separate_light=F)

# convert to data frame, use c() to remove attributes
ae <- mat_to_df(mdq@qsip[['atom_excess']], 'excess')

ae$group <- sub('^(\\w{1})(.*)', '\\1', ae$RepID)

ae_qsip <- ae


# manual calculation-------------------------------------------
# add treatment (isotope) and group data for grouping
ae <- merge(mw_110, unique(mdl[,c('RepID', 'tmt', 'group')]), all.x=T)

# separate light and heavy MWs
aeh <- ae[ae$tmt=='18O',]
ael <- ae[ae$tmt=='16O',]

# DON'T average light MWs, they've already been calculated based on averages

# calculate MW-heavymax
adjust <- 12.07747
nat_abund <- 0.002011429
ael$mw_max <- (adjust + ael$mw_manual)

# merge together, disregarding tmt
ae <- merge(aeh[,!names(aeh) %in% 'tmt'],
            ael[,!names(ael) %in% 'tmt'],
            by=c('RepID', 'group', 'OTU'),
            suffixes=c('_label', '_light'))

# calculate atom excess
ae <- within(ae, {
  excess <- ((mw_manual_label - mw_manual_light)/(mw_max - mw_manual_light)) * (1 - nat_abund)
})


# compare calculations-----------------------------------------
ae_110 <- merge(ae_qsip,
                ae[,c('RepID', 'group', 'OTU', 'excess')],
                by=c('RepID', 'group', 'OTU'),
                suffixes=c('_qsip', '_manual'))
