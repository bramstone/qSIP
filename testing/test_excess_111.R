# calc code 111

# combine new replicate numbers with data
mdq <- md
mdq@sam_data$rep_number_111 <- rep_matching$rep_number_111

mdq <- specify_qsip(mdq,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    rep_group='group',
                    rep_num='rep_number_111',
                    iso='18O',
                    iso_trt='tmt')

# calculate atom excess
mdq <- calc_excess(mdq, separate_label=T, separate_light=T)

# convert to data frame, use c() to remove attributes
ae <- mat_to_df(mdq@qsip[['atom_excess']], 'excess')

ae_qsip <- ae


# manual calculation-------------------------------------------
# add treatment (isotope) data, groups, and replicate numbers for grouping
ae <- merge(mw_111, unique(rep_matching[,c('RepID', 'tmt', 'group', 'rep_number_111')]), all.x=T)

# separate light and heavy WADs
aeh <- ae[ae$tmt=='18O',]
ael <- ae[ae$tmt=='16O',]

# merge light and heavy, matching by replicate number, OTU, and group
ae <- merge(aeh[,!names(aeh) %in% c('tmt', 'RepID')],
            ael[,!names(aeh) %in% c('tmt', 'RepID')],
            by=c('OTU', 'group', 'rep_number_111'),
            suffixes=c('_label', '_light'))

# calculate atom excess
adjust <- 12.07747
nat_abund <- 0.002011429
#
ae <- within(ae, {
  mw_max <- (adjust + mw_manual_light)
  excess <- ((mw_manual_label - mw_manual_light)/(mw_max - mw_manual_light)) * (1 - nat_abund)
})

# get RepID again
ae <- merge(unique(rep_matching[,c('RepID', 'tmt', 'group', 'rep_number_111')]),
            ae,
            all.y=T)
ae$rep_number_111 <- NULL

# compare calculations-----------------------------------------
ae_111 <- merge(ae_qsip,
                ae[,c('RepID', 'OTU', 'excess')],
                by=c('RepID', 'OTU'),
                suffixes=c('_qsip', '_manual'))
