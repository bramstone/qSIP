# Testing of the calculation of MW values between qSIP package and manual
# calculation code: 111

############################################
# qsip package
############################################

# first need to modify replicate number column so that they are unique in each group
rep_matching <- unique(md@sam_data[,c('RepID', 'tmt', 'group')])
rep_matching <- as(rep_matching, 'data.frame')
rep_matching <- split(rep_matching, interaction(rep_matching$tmt, rep_matching$group))
# sapply(rep_matching, nrow) # check row number - 6 in each 16O and 18O in each group
# add in numbers - NOTE: sample names are aligned alphabetically such that the numbers in each treatment correspond
rep_matching <- lapply(rep_matching, function(x) {x$rep_number_111 <- 1:nrow(x); x})
rep_matching <- do.call(rbind, rep_matching)
# get into same format as md@sam_data
rep_matching <- merge(as(md@sam_data, 'data.frame'), rep_matching, all.x=T)
rep_matching <- rep_matching[match(md@sam_data$sample.rep.fraction, rep_matching$sample.rep.fraction),]

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

# calculate WADs
mdq <- calc_mw(mdq, separate_label=T, separate_light=T)

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

mw_111 <- merge(mw_qsip, mw_manual,
                by=c('OTU', 'RepID', 'group', 'tmt'),
                suffixes=c('_qsip', '_manual'))
