#' Calculation of weighted average density differences
#'
#' Calculates differences in weighted average densities across replicates due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#' @param filter Logical value specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#' @param correction Logical value indicating whether or not to apply tube-level correction to labeled WAD values.
#' @param offset_taxa Value from 0 to 1 indicating the percentage of the taxa to utilize for calculating offset correction values.
#'   Taxa are ordered by lowest difference in WAD values.
#'   Default is \code{0.1} indicating 10 percent of taxa with the lowest difference in WAD values.
#' @param separate_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, unlabeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight in the absence of isotope addition.
#' @param separate_label Logical value indicating whether or not WAD-label scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, labeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight as a result of isotope addition. The default is \code{TRUE}.
#' @param recalc Logical value indicating whether or not to recalculate WAD values or use existing values. Default is \code{TRUE}.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities have not been calculated
#'   beforehand, \code{calc_d_wad} will compute those first.
#'
#'   More details about tube-level corrections.
#'
#'   Unless a particular site is known to have significantly higher levels of naturally occurring isotopes, it is recommended to
#'   combine unlabeled (light) WAD scores in order to produce an average unlabeled WAD value that more accurately represents a
#'   taxon's natural WAD value.
#'
#'   The equation for the difference in weighted average density of taxon \emph{i}, \eqn{\Delta W_{i}}, is given by:
#'
#'   \deqn{\Delta W_{i} = W_{Lab,i} - W_{Light,i}}
#'
#'   Where \eqn{W_{Lab,i}} indicates the weighted average density of taxon \emph{i} in heavy isotope-labeled treatment(s)
#'   and \eqn{W_{Light,i}} indicates that taxon's corresponding weighted average density in light, or unlabeled treatment(s)
#'
#' @return \code{calc_d_wad} adds two S4 Matrix objects to the \code{data@@qsip@@.Data} slot, one for differences
#'   in weighted average density, and the other for weighted average density values of light treatments only (to be used in
#'   future calculations). The row and column specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}},
#'   meaning if taxa are listed on the table rows, they will in the resulting S4 Matrix objects
#'
#' @seealso \code{\link{calc_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'
#'  # Calculate weighted average density differences
#'
#' @references
#'  Hungate, Bruce, \emph{et al.} 2015. Quantitative microbial ecology through stable isotope probing.
#'  \emph{Applied and Environmental Microbiology} \strong{81}: 7570 - 7581.
#'
#' @export

calc_d_wad <- function(data, filter=FALSE, return_diffs=FALSE, correction=FALSE, offset_taxa=0.1, separate_light=FALSE, separate_label=TRUE, recalc=TRUE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if recalculation wanted, or if WAD values don't exist, calculate those first, this will also handle rep_id validity
  if(recalc | is.null(data@qsip[['wad']])) data <- calc_wad(data, filter=filter)
  #if(length(data@qsip@rep_group)==0) stop('Must specify replicate groupings with rep_group')
  if(length(data@qsip@iso_trt)==0) stop('Must specify treatment and controls with iso_trt')
  trt_levels <- unique(data@sam_data[[data@qsip@iso_trt]])
  trt_levels <- trt_levels[!is.na(trt_levels)]
  if(length(trt_levels) > 2) stop('More than two treatment levels present for this comparison')
  if(offset_taxa > 1 && correction) offset_taxa <- 1
  if(offset_taxa <= 0 && correction) stop('Must specify non-negative proportion of taxa generate offset values for WAD correction')
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip[['wad']]
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  n_taxa <- ncol(ft)
  tax_names <- colnames(ft)
  # manipulate data matrix and calculate
  # split by replicate groups, but keep track of light and heavy fractions
  ft <- valid_samples(data, ft, 'iso')
  iso_group <- ft[[2]]; ft <- ft[[1]]
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F, keep_names=1)
  # If there is no replicate grouping (i.e., all replicates in a treatment are grouped)...
  iso_group2 <- unique(iso_group[,!names(iso_group) %in% 'replicate']) # only get unique elements to match levels in ft
  iso_group2 <- iso_group2[match(names(ft), iso_group2$interaction),]
  # average light values across all groups
  if(!separate_light) {
    grouped_light <- ft[which(as.numeric(iso_group2$iso)==1)]
    grouped_light <- do.call(rbind, grouped_light)
    grouped_light <- colMeans(grouped_light, na.rm=T)
    grouped_light[is.nan(grouped_light)] <- NA
  }
  #
  # Determine which groups to keep and which to remove based on sample grouping
  # If grouping labeled values
  if(!separate_label) {
    # calculate average WAD per taxa for each replicate group
    ft <- base::lapply(ft, colMeans, na.rm=T)
    # remove any NaNs resulting from when a taxon is missing in all replicates
    ft <- base::lapply(ft, function(x) {x[is.nan(x)] <- NA; x})
    if(length(data@qsip@rep_group)==0) {
      # If there are no groups, replicates are designated either light or label
      keep_groups <- !logical(2)
    } else {
      # If there are groups, need to warn user if not all labeled groups have corresponding light group
      keep_groups <- !logical(nlevels(iso_group$interaction))
      for(i in 1:length(ft)) {
        # use numbers to reference non-labeled additions since they're isotope agnostic
        # any NA values result here when a taxa is completely missing from a heavy or light treatment in a replicate group
        which_light <- which(as.numeric(iso_group2$grouping)==i &
                               as.numeric(iso_group2$iso)==1)
        which_heavy <- which(as.numeric(iso_group2$grouping)==i &
                               as.numeric(iso_group2$iso)==2)
        # use grouped light values if specified
        # If there's no light OR no heavy treatment for a group of replicates, remove them
        # Only worry about unpaired light groups if separate_light==TRUE
        if((length(which_light)==0 && separate_light) || length(which_heavy)==0) {
          warning('Labeled or unlabeled isotope treatment missing in replicate group(s): ', names(ft)[i],
                  '\nRemoving sample(s): ', paste(as.character(iso_group[iso_group$grouping==names(ft)[i], 'replicate']), collapse=', '),
                  ' - from calculation', call.=FALSE)
          keep_groups[i] <- FALSE
          next
        }
      }
      # end for-loop
    }
  # If keeping labeled replicates separate
  } else if(separate_label) {
    # If replicates are not grouped, don't worry about unpaired groups, worry about unpaired replicates
    if(length(data@qsip@rep_group)==0) {
      if(separate_light) {
        # check if number of light replicates equals matches number of labeled replicates
        # if these don't match, user needs to group light values, or both label and light
      }
      keep_groups <- !logical(2)
    } else { # use a for-loop to check labeled and light replicates in each group
      n_reps <- iso_group[as.numeric(iso_group$iso)==2,]
      n_reps <- table(n_reps$grouping)
      # For each repliate group: identify which elements of ft are light and which are heavy, then get difference
      keep_groups <- !logical(length(n_reps))
      for(i in 1:length(ft)) {
        # use numbers to reference non-labeled additions since they're isotope agnostic
        # any NA values result here when a taxa is completely missing from a heavy or light treatment in a replicate group
        which_light <- which(as.numeric(iso_group2$grouping)==i &
                               as.numeric(iso_group2$iso)==1)
        which_heavy <- which(as.numeric(iso_group2$grouping)==i &
                               as.numeric(iso_group2$iso)==2)
        # If there's no light OR no heavy treatment for a group of replicates, remove them
        # Only worry about unpaired light groups if separate_light==TRUE
        if((length(which_light)==0 && separate_light) || length(which_heavy)==0) {
          warning('Labeled or unlabeled isotope treatment missing in replicate group(s): ', names(ft)[i],
                  '\nRemoving sample(s): ', paste(as.character(iso_group[iso_group$grouping==names(ft)[i], 'replicate']), collapse=', '),
                  ' - from calculation', call.=FALSE)
          keep_groups[i] <- FALSE
          next
        }
      }
      # end for-loop
    }
  }
  # organize and add new data as S4 matrix
  # return weighted average densities of light calcs only
  ft <- ft[iso_group2$grouping %in% levels(iso_group2$grouping)[keep_groups]] # remove unpaired & dropped groups
  iso_group2 <- iso_group2[iso_group2$grouping %in% levels(iso_group2$grouping)[keep_groups],] # remove unpaired & dropped groups
  wl <- ft[which(as.numeric(iso_group2$iso)==1)]
  wh <- ft[which(as.numeric(iso_group2$iso)==2)]
  if(length(data@qsip@rep_group)!=0) {
    names(wl) <- sub('\\w{3}\\.', '', names(wl))
    names(wh) <- sub('\\w{3}\\.', '', names(wh))
  }
  # Apply tube-level correction?
  if(correction) {
    if(length(data@qsip@rep_group)!=0) light <- grouped_light else light <- do.call(rbind, wl)
    # only utilize taxa that are present in every replicate
    #pa <- split_data(data, pa, iso_group$interaction, grouping_w_phylosip=F)
    #pa <- pa[which(as.numeric(iso_group2$iso)==2)]
    #pa <- base::lapply(pa, function(x) colSums(x) / nrow(x)) # creates frequency code 0 - 1
    #names(pa) <- names(wl)
    shift <- base::lapply(wh, '-', light)
    #shift <- base::Map(function(x, y) x[y==1], shift, pa)
    # sort by lowest diff WAD
    shift <- base::lapply(shift, sort, decreasing=F)
    # calculate median of lowest x% of diff WADs, default is 10%
    shift <- base::lapply(shift, function(x) median(x[1:floor(offset_taxa * length(x))], na.rm=T))
    # subtract shift from labeled WAD values
    wh <- base::Map('-', wh, shift)
  }
  wh <- base::lapply(wh, function(x) {x[is.nan(x)] <- NA; x})
  if(separate_label) {
    wh <- do.call(rbind, wh)
  }
  if(!separate_light) wl <- grouped_light
  data <- collate_results(data, wh, tax_names=tax_names, 'wad_label', sparse=TRUE)
  data <- collate_results(data, wl, tax_names=tax_names, 'wad_light', sparse=TRUE)
  return(data)
}
