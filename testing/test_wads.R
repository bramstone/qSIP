# Testing of the calculation of WAD values between qSIP package and manual


############################################
# qsip package
############################################

mdq <- specify_qsip(md,
                    abund='qPCR.16S.copies.ul',
                    density='density.g.ml',
                    rep_id='RepID',
                    iso='18O',
                    iso_trt='tmt')

# calculate WADs
mdq <- calc_wad(mdq)

# convert to data frame
wads <- mat_to_df(mdq@qsip[['wad']], 'wad')
wads_qsip <- wads[!is.na(wads$wad),]


############################################
# manual calculation
############################################

# first, re-calculate relative abundances, since that's what the calc_wad function does
wads <- split(mdl, mdl$Sample)
wads <- lapply(wads, function(x) {x$Abundance <-  x$Abundance / sum(x$Abundance, na.rm=T); x})
wads <- do.call(rbind, wads)

# split by sample/replicate and taxon
wads <- split(wads, interaction(wads$RepID, wads$OTU))

# If any sample has NA for density, make those qPCR abundances NA as well
# calculate relative abundance of each taxa in each fraction
wads <- lapply(wads, function(x) {
  x$rel_abund <- x$Abundance * x$qPCR.16S.copies.ul
  x$rel_abund <- x$rel_abund / sum(x$rel_abund, na.rm=T)
  x
})

# multiply relative abundances by density, then sum
wads <- lapply(wads, function(x) {
  x$wad <- x$density.g.ml * x$rel_abund
  x$wad <- sum(x$wad, na.rm=T)
  x
})

# return single row, remove fraction info
wads <- lapply(wads, function(x) {
  x <- x[1,]
  x[,!names(x) %in% grep('fraction|Sample', names(x), value=T)]
})

# re-combine
wads <- do.call(rbind, wads)

# remove NA wads (where wad=0)
wads <- wads[wads$wad > 0,]
wads <- wads[,c('RepID', 'OTU', 'wad')]
rownames(wads) <- NULL
wads_man <- wads


############################################
# compare calculations
############################################

wads_qsip <- wads_qsip[order(wads_qsip$OTU, wads_qsip$RepID),]
wads_man <- wads_man[order(wads_man$OTU, wads_man$RepID),]

wads <- merge(wads_qsip, wads_man,
              by=c('RepID', 'OTU'),
              suffixes=c('_qsip', '_manual'))
