#' Grouping of weighted average fraction
#'
#' Groups weighted average fraction values of microbial taxa by isotope treatment
#'
#' @param data Data as a \code{phyloseq} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#' @param correction Logical value indicating whether or not to apply tube-level correction to labeled WAF values.
#' @param offset_taxa Value from 0 to 1 indicating the percentage of the taxa to utilize for calculating offset correction values.
#'   Taxa are ordered by lowest difference in WAD values.
#'   Default is \code{0.1} indicating 10 percent of taxa with the lowest difference in WAD values.
#' @param separate_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, unlabeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight in the absence of isotope addition.
#' @param separate_label Logical value indicating whether or not WAD-label scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, labeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight as a result of isotope addition. The default is \code{TRUE}.
#' @param global_light Logical value indicating whether or not to use WAD-light scores that are global averages (\emph{i.e.,} averaged across
#'   all samples rather than averaged across any specified replicate groups). The default is \code{FALSE}.
#' @param rel_abund Logical value specifying if relative abundances of taxa are to be calculated prior to calculations. The default is \code{TRUE}.
#'   This parameter is passed to \code{calc_wad}.
#' @param recalc Logical value indicating whether or not to recalculate WAD, WAD heavy, and WAD light values or use existing values. Default is \code{TRUE}.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average fractions
#'   have not been calculated beforehand, \code{calc_d_waf} will compute those first.
#'
#'   Explain consequences of different grouping actions on results.
#'
#' @return \code{calc_d_waf} adds three S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of weighted average fraction for each taxon at each group of replicates in the labeled and unlabeled groups. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#'   If timepoint is specified through \code{data@@qsip@@timepoint}, then \code{calc_d_waf} will group molecular weights by timepoint as a grouping factor.
#'   If both \code{data@@qsip@@timepoint} and \code{data@@qsip@@rep_group} are specified, then \code{calc_d_waf} will group weighted average fractions by the
#'   interaction of these terms. There is no option to change this behavior except by removing the specification for \code{data@@qsip@@timepoint}.
#'   This is because calculating weighted average fraction for separate timepoints is a necessity in \code{\link{calc_pop}}
#'
#' @seealso \code{\link{calc_waf}}
#'
#' @examples
#'  # Load in example data
#'
#'
#' @export

calc_d_waf <- function(data, filter=FALSE, correction=FALSE, offset_taxa=0.1, separate_light=FALSE, separate_label=TRUE,
                    global_light=FALSE, recalc=TRUE, rel_abund=TRUE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(separate_label && separate_light && length(data@qsip@rep_num)==0) {
    stop('Must specify replicate matching with data@qsip@rep_num', call.=FALSE)
  }
  # if WAF values don't exist, or if recalculation wanted, calculate those first
  # this will also handle rep_id validity (through calc_waf)
  if(recalc | is.null(data@qsip[['waf']])) {
    data <- calc_waf(data, filter=filter, rel_abund=rel_abund)
  }
  # extract WAF values and convert to S3 matrix
  waf <- data@qsip[['waf']]
  waf <- as(waf, 'matrix')
  if(phyloseq::taxa_are_rows(data)) waf <- t(waf)
  n_taxa <- ncol(waf)
  tax_names <- colnames(waf)
  #
  # if timepoint specified in data but not group, timepoint becomes the grouping
  if(length(data@qsip@timepoint)==1 && length(data@qsip@rep_group)==0) {
    data@qsip@rep_group <- data@qsip@timepoint
  # else if timepoint AND group specified, group:timepoint interaction become the grouping
  } else if(length(data@qsip@timepoint)==1 && length(data@qsip@rep_group)==1) {
    orig_group <- data@qsip@rep_group
    data@sam_data$timepoint.rep_group <- interaction(data@sam_data[[data@qsip@rep_group]],
                                                     data@sam_data[[data@qsip@timepoint]])
  }
  # split matrix by replicate, remove samples with NA for isotope trt, and keep track of light and heavy fractions
  waf <- valid_samples(data, waf, 'iso')
  iso_group <- waf[[2]]; waf <- waf[[1]]
  # separate labeled and unlabeled samples
  wh <- waf[as.numeric(iso_group$iso)==2,]
  wl <- waf[as.numeric(iso_group$iso)==1,]
  iso_h <- droplevels(iso_group[as.numeric(iso_group$iso)==2,])
  iso_l <- droplevels(iso_group[as.numeric(iso_group$iso)==1,])
  #
  # ----------------
  # Different methods of grouping data produce slightly different calculations
  # These are coded as a 3-digit binary code based on 3 grouping criteria
  #   - Replicate separated by group:  0/1 (no/yes)
  #   - Labeled densities separated:   0/1
  #   - Unlabeled densities separated: 0/1
  # Example: 010 means calculate MW with separate labeled values, and averaged unlabeled values
  # Example: 110 means calculate MW with separate labeled values, and unlabeled values averaged in each group
  # ----------------
  #
  # if there is no replicate grouping
  if(length(data@qsip@rep_group)==0) {
    if(!separate_label) {
      #
      wh <- colMeans(wh, na.rm=TRUE)
      wh[is.nan(wh)] <- NA
      # ...........................................
      if(!separate_light) { # CODE 000
        #
        wl <- colMeans(wl, na.rm=TRUE)
        wl[is.nan(wl)] <- NA
        # WAF correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        # ...........................................
      } else if(separate_light) { # CODE 001
        #
        # WAF correction
        if(correction) {
          shift <- sweep(wl, 2, wh, function(wL, wH) wH - wL)
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        shift <- sweep(wl, 2, wh, function(wL, wH) wH - wL)
        #
      }
    } else if(separate_label) {
      # ...........................................
      if(!separate_light) { # CODE 010
        #
        wl <- colMeans(wl, na.rm=T)
        wl[is.nan(wl)] <- NA
        # WAF correction
        if(correction) {
          shift <- sweep(wh, 2, wl, function(wH, wL) wH - wL)
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        shift <- sweep(wh, 2, wl, '-')
        #
        # ...........................................
      } else if(separate_light) { # CODE 011
        #
        # evaluate that individual samples align for comparison
        # remove comparisons with missing labeled samples, replace unlabeled samples with global unlabeled average
        good_vals <- match_reps(data, wh, wl, iso_group, rep_group=FALSE)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        # WAF correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        #
      }
    }
  # if there is replicate grouping
  } else if(length(data@qsip@rep_group)==1) {
    #
    wh <- split_data(data, wh, iso_h$interaction, grouping_w_phylosip=FALSE, keep_names=1)
    wl <- split_data(data, wl, iso_l$interaction, grouping_w_phylosip=FALSE, keep_names=1)
    #
    if(!separate_label) {
      #
      wh <- base::lapply(wh, colMeans, na.rm=TRUE)
      wh <- base::lapply(wh, function(x) {x[is.nan(x)] <- NA; x})
      # ...........................................
      if(!separate_light) { # CODE 100
        #
        wl <- base::lapply(wl, colMeans, na.rm=TRUE)
        wl <- base::lapply(wl, function(x) {x[is.nan(x)] <- NA; x})
        #
        # evaluate that grouped samples align for comparison
        # separately, global light will just replace ALL unlabeled WADs with the global average
        good_vals <- match_groups(data, wh, wl, iso_group, global_light=global_light)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        wh <- do.call(rbind, wh)
        wl <- do.call(rbind, wl)
        # WAF correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        #
        # ...........................................
      } else if(separate_light) { # CODE 101
        #
        # evaluate that there are data in each group comparison
        # remove comparisons with missing labeled samples, replace unlabeled samples with global unlabeled average
        # separately, global light will just replace ALL unlabeled WADs with the global average
        good_vals <- match_groups(data, wh, wl, iso_group, global_light=global_light)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        # WAF correction
        if(correction) {
          shift <- base::Map(function(x, y) sweep(x, 2, y, function(wL, wH) wH - wL),
                             wl, wh)
          shift <- base::lapply(shift, sort, decreasing=FALSE)
          shift <- base::lapply(shift, function(x) median(x[1:floor(offset_taxa * length(x))], na.rm=T))
          wh <- base::lapply(wh, '-', shift)
        }
        #
      }
    } else if(separate_label) {
      # ...........................................
      if(!separate_light) { # CODE 110
        #
        wl <- base::lapply(wl, colMeans, na.rm=TRUE)
        wl <- base::lapply(wl, function(x) {x[is.nan(x)] <- NA; x})
        #
        # evaluate that there are data in each group comparison
        # remove comparisons with missing labeled samples, replace unlabeled groups with global average
        # separately, global light will just replace ALL unlabeled WADs with the global average
        good_vals <- match_groups(data, wh, wl, iso_group, global_light=global_light)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        # WAD correction
        if(correction) {
          shift <- base::Map(function(x, y) sweep(x, 2, y, '-'), wh, wl)
          shift <- base::lapply(shift, sort, decreasing=FALSE)
          shift <- base::lapply(shift, function(x) median(x[1:floor(offset_taxa * length(x))], na.rm=T))
          wh <- base::lapply(wh, '-', shift)
          #
        }
        # ...........................................
      } else if(separate_light) { # CODE 111
        #
        # evaluate that individual samples align for comparison (RETURNS MATRICES)
        # remove comparisons with missing labeled samples, replace unlabeled samples with average unlabeled
        # separately, global light will just replace ALL unlabeled WADs with the global average
        good_vals <- match_reps(data, wh, wl, iso_group, rep_group=TRUE, global_light=global_light)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        # WAD correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        shift <- wh - wl
        #
      }
    }
  }
  # organize and add new data as S4 matrices
  data <- collate_results(data, wh, tax_names=tax_names, 'waf_label', sparse=TRUE)
  data <- collate_results(data, wl, tax_names=tax_names, 'waf_light', sparse=TRUE)
  # add attributes summarizing the calculation method
  output_attr <- c(rep_group=as.logical(length(data@qsip@rep_group)),
                   sep_label=separate_label,
                   sep_light=separate_light)
  attributes(data@qsip[['waf_label']])$calc_method <- output_attr
  attributes(data@qsip[['waf_light']])$calc_method <- output_attr
  #
  # if timepoint used as temporary grouping, remove the assignment
  if(length(data@qsip@timepoint)==1 && !exists('orig_group')) {
    data@qsip@rep_group <- character()
    # else remove timepoint:group interaction, and replace with original grouping
  } else if(length(data@qsip@timepoint)==1 && exists('orig_group')) {
    data@sam_data$timepoint.rep_group <- NULL
    data@qsip@rep_group <- orig_group
  }
  return(data)
}
