# calc code 100

mdq <- specify_qsip(mdq,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_group='group',
                    rep_id='RepID',
                    iso='18O',
                    iso_trt='tmt')

# calculate atom excess
mdq <- calc_excess(mdq, separate_label=F, separate_light=F)

# convert to data frame, use c() to remove attributes
ae <- mat_to_df(mdq@qsip[['atom_excess']], 'excess')

names(ae)[1] <- 'group'
ae$group <- sub('(.*)\\.(.*)', '\\2', ae$group)

ae_qsip <- ae


# manual calculation-------------------------------------------
# add treatment (isotope) and group data for grouping
ae <- merge(mw_100, unique(mdl[,c('RepID', 'tmt', 'group')]), all.x=T)

# average MWs by OTU, group, and treatment
ae <- aggregate(mw_manual ~ OTU + group + tmt, ae, mean, na.rm=T)

# separate light and heavy MWs
aeh <- ae[ae$tmt=='18O',]
ael <- ae[ae$tmt=='16O',]

# calculate MW-heavymax
adjust <- 12.07747
nat_abund <- 0.002011429
ael$mw_max <- (adjust + ael$mw_manual)

# merge, disregarding tmt
ae <- merge(aeh[,!names(aeh) %in% 'tmt'],
            ael[,!names(ael) %in% 'tmt'],
            by=c('group', 'OTU'),
            suffixes=c('_label', '_light'))

# calculate atom excess
ae <- within(ae, {
  excess <- ((mw_manual_label - mw_manual_light)/(mw_max - mw_manual_light)) * (1 - nat_abund)
})


# compare calculations-----------------------------------------
ae_100 <- merge(ae_qsip,
                ae[,c('group', 'OTU', 'excess')],
                by=c('group', 'OTU'),
                suffixes=c('_qsip', '_manual'))
