# calc code 011

# combine new replicate numbers with data
mdq <- md
mdq@sam_data$rep_number <- rep_matching_011$rep_number

mdq <- specify_qsip(mdq,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    rep_num='rep_number',
                    iso='18O',
                    iso_trt='tmt')

# calculate atom excess
mdq <- calc_excess(mdq, separate_label=T, separate_light=T)

# convert to data frame, use c() to remove attributes
ae <- mat_to_df(mdq@qsip[['atom_excess']], 'excess')


ae_qsip <- ae


# manual calculation-------------------------------------------
ae <- mw_011[,!names(mw_011) %in% 'mw_qsip']
ae <- merge(ae, unique(rep_matching_011[,c('RepID', 'tmt', 'rep_number')]), all.x=T)

# separate MW-light and MW-label
aeh <- ae[ae$tmt=='18O',]
ael <- ae[ae$tmt=='16O',]

# calculate MW-heavymax
adjust <- 12.07747
nat_abund <- 0.002011429
ael$mw_max <- (adjust + ael$mw_manual)

# merge, disregarding tmt
ae <- merge(aeh[,!names(aeh) %in% 'tmt'],
            ael[,!names(ael) %in% 'tmt'],
            by=c('rep_number', 'OTU'),
            suffixes=c('_label', '_light'))

# calculate atom excess
ae <- within(ae, {
  excess <- ((mw_manual_label - mw_manual_light)/(mw_max - mw_manual_light)) * (1 - nat_abund)
})


# compare calculations-----------------------------------------
ae_011 <- merge(ae_qsip,
                ae[,c('RepID_label', 'OTU', 'excess')],
                by.x=c('RepID', 'OTU'),
                by.y=c('RepID_label', 'OTU'),
                suffixes=c('_qsip', '_manual'))
