# Testing of the calculation of MW values between qSIP package and manual
# calculation code: 011

############################################
# qsip package
############################################

# first need to modify replicate number column so that they are unique across all groups
rep_matching_011 <- unique(md@sam_data[,c('RepID', 'tmt')])
rep_matching_011 <- as(rep_matching_011, 'data.frame')
rep_matching_011 <- split(rep_matching_011, rep_matching_011$tmt)
# sapply(rep_matching_011, nrow) # check row number - 12 in both 16O and 18O - it matches
# add in numbers - NOTE: sample names are aligned alphabetically such that the numbers in each treatment correspond
rep_matching_011 <- lapply(rep_matching_011, function(x) {x$rep_number <- 1:nrow(x); x})
rep_matching_011 <- do.call(rbind, rep_matching_011)
# get into same format as md@sam_data
rep_matching_011 <- merge(as(md@sam_data, 'data.frame'), rep_matching_011, all.x=T)
rep_matching_011 <- rep_matching_011[match(md@sam_data$sample.rep.fraction, rep_matching_011$sample.rep.fraction),]

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

# calculate WADs
mdq <- calc_mw(mdq, separate_label=T, separate_light=T)

mwh <- mat_to_df(mdq@qsip[['mw_label']], 'mw')
mwl <- mat_to_df(mdq@qsip[['mw_light']], 'mw')

mwh$tmt <- '18O'
mwl$tmt <- '16O'

# convert to data frame, use c() to remove attributes
mw <- rbind(mwh, mwl)

mw_qsip <- mw


############################################
# manual calculation
############################################

# add treatment (isotope) data and replicate numbers for grouping
mw <- merge(wads_man, unique(rep_matching_011[,c('RepID', 'tmt', 'rep_number')]), all.x=T)

# separate light and heavy WADs
mwh <- mw[mw$tmt=='18O',]
mwl <- mw[mw$tmt=='16O',]

# merge back in, disregarding treatment
mw <- merge(mwh[,!names(mwh) %in% 'tmt'],
            mwl[,!names(mwl) %in% 'tmt'],
            by=c('OTU', 'rep_number'),
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
              idvar=c('rep_number', 'OTU'),
              varying=list(c('mw_label', 'mw_light')),
              drop=c('wad_label', 'wad_light', 'RepID_label', 'RepID_light'),
              times=c('18O', '16O'),
              timevar='tmt',
              v.names='mw',
              direction='long')
rownames(mw) <- NULL

# get RepID again
mw <- merge(unique(rep_matching_011[,c('RepID', 'tmt', 'rep_number')]),
            mw,
            all.y=T)
mw$rep_number <- NULL

mw_manual <- mw


############################################
# compare calculations
############################################

mw_011 <- merge(mw_qsip, mw_manual,
                by=c('OTU', 'RepID', 'tmt'),
                suffixes=c('_qsip', '_manual'))
