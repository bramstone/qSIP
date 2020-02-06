#' Calculation of molecular weights
#'
#' Calculates molecular weights of microbial taxa due to isotope incorporation
#'
#' @param data Data as a \code{phyloseq} object
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
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
#' @param global_light Logical value indicating whether or not to use WAD-light scores that are global averages (\emph{i.e.,} averaged across
#'   all samples rather than averaged across any specified replicate groups). The default is \code{FALSE}.
#' @param rel_abund Logical value specifying if relative abundances of taxa are to be calculated prior to calculations. The default is \code{TRUE}.
#'   This parameter is passed to \code{calc_wad}.
#' @param recalc Logical value indicating whether or not to recalculate WAD, WAD heavy, and WAD light values or use existing values. Default is \code{TRUE}.
#'
#' @details Some details about proper isotope control-treatment factoring. If weighted average densities or the change in weighted average densities
#'   have not been calculated beforehand, \code{calc_mw} will compute those first.
#'
#'   The equations for calculating the molecular weights of taxon \emph{i}, designated \eqn{M_{Lab,i}} for labeled and \eqn{M_{Light,i}} for
#'   unlabeled, are:
#'
#'   \deqn{M_{Light,i} = 0.496 \cdot G_{i} + 307.691}
#'   \deqn{M_{Lab,i} = \left( \frac{\Delta W}{W_{Light,i}} + 1 \right) \cdot M_{Light,i}}
#'
#'   Where
#'
#'   \deqn{G_{i} = \frac{1}{0.083506} \cdot (W_{Light,i} - 1.646057)}
#'   Which indicates the GC content of taxon \emph{i} based on the density of its DNA when unlabeled
#'
#'   Explain consequences of different grouping actions on results.
#'
#' @return \code{calc_mw} adds three S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of molecular weights for each taxon at each group of replicates in the labeled and unlabeled groups. The row and column
#'   specifications will mirror those of the \code{phylosip}'s \code{\link{otu_table}}, meaning if taxa are listed on the table rows,
#'   they will in the resulting S4 Matrix class.
#'
#' @seealso \code{\link{calc_wad}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate molecular weights
#'
#' @references
#'  Hungate, Bruce, \emph{et al.} 2015. Quantitative microbial ecology through stable isotope probing.
#'  \emph{Applied and Environmental Microbiology} \strong{81}: 7570 - 7581.
#'
#' @export

calc_mw <- function(data, filter=FALSE, correction=FALSE, offset_taxa=0.1, separate_light=FALSE, separate_label=TRUE,
                    global_light=FALSE, recalc=TRUE, rel_abund=TRUE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  if(separate_label && separate_light && length(data@qsip@rep_num)==0) {
    stop('Must specify replicate matching with data@qsip@rep_num', call.=FALSE)
  }
  # if WAD values don't exist, or if recalculation wanted, calculate those first
  # this will also handle rep_id validity (through calc_wad)
  if(recalc | is.null(data@qsip[['wad']])) {
    data <- calc_wad(data, filter=filter, rel_abund=rel_abund)
  }
  # extract WAD values and convert to S3 matrix
  ft <- data@qsip[['wad']]
  ft <- as(ft, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  n_taxa <- ncol(ft)
  tax_names <- colnames(ft)
  # split matrix by replicate, remove samples with NA for isotope trt, and keep track of light and heavy fractions
  ft <- valid_samples(data, ft, 'iso')
  iso_group <- ft[[2]]; ft <- ft[[1]]
  # separate labeled and unlabeled samples
  wh <- ft[as.numeric(iso_group$iso)==2,]
  wl <- ft[as.numeric(iso_group$iso)==1,]
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
        # WAD correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        # calculate GC content of each taxa
        gc <- (1 / 0.083506) * (wl - 1.646057)
        # calculate mol. weight of taxa
        mw_l <- (0.496 * gc) + 307.691
        mw_h <- (((wh - wl)/wl) + 1) * mw_l
        # ...........................................
      } else if(separate_light) { # CODE 001
        #
        # calculate GC content of each taxa
        gc <- (1 / 0.083506) * (wl - 1.646057)
        # calculate mol. weight of taxa
        mw_l <- (0.496 * gc) + 307.691
        mw_h <- sweep(wl, 2, wh, function(wL, wH) (((wH - wL)/wL) + 1))
        mw_h <- mw_h * mw_l
        shift <- sweep(wl, 2, wh, function(wL, wH) wH - wL)
        #
      }
    } else if(separate_label) {
      # ...........................................
      if(!separate_light) { # CODE 010
        #
        wl <- colMeans(wl, na.rm=T)
        wl[is.nan(wl)] <- NA
        # WAD correction
        if(correction) {
          shift <- sweep(wh, 2, wl, function(wH, wL) wH - wL)
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        # calculate GC content of each taxa
        gc <- (1 / 0.083506) * (wl - 1.646057)
        # calculate mol. weight of taxa
        mw_l <- (0.496 * gc) + 307.691
        mw_h <- sweep(wh, 2, wl, function(wH, wL) (((wH - wL)/wL) + 1))
        mw_h <- sweep(mw_h, 2, mw_l, '*')
        shift <- sweep(wh, 2, wl, '-')
        # ...........................................
      } else if(separate_light) { # CODE 011
        #
        # evaluate that individual samples align for comparison
        # remove comparisons with missing labeled samples, replace unlabeled samples with global unlabeled average
        good_vals <- match_reps(data, wh, wl, iso_group, rep_group=FALSE)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        # WAD correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        # calculate GC content of each taxa
        gc <- (1 / 0.083506) * (wl - 1.646057)
        # calculate mol. weight of taxa
        mw_l <- (0.496 * gc) + 307.691
        mw_h <- (((wh - wl)/wl) + 1) * mw_l
      }
    }
  }
  # if there is replicate grouping
  if(length(data@qsip@rep_group)==1) {
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
        # remove comparisons with missing labeled samples, replace unlabeled samples with global unlabeled average
        good_vals <- match_groups(data, wh, wl, iso_group)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        wh <- do.call(rbind, wh)
        wl <- do.call(rbind, wl)
        # WAD correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        # calculate GC content of each taxa
        gc <- (1 / 0.083506) * (wl - 1.646057)
        # calculate mol. weight of taxa
        mw_l <- (0.496 * gc) + 307.691
        mw_h <- (((wh - wl)/wl) + 1) * mw_l
        # ...........................................
      } else if(separate_light) { # CODE 101
        #
        wh <- base::lapply(wh, colMeans, na.rm=TRUE)
        wh <- base::lapply(wh, function(x) {x[is.nan(x)] <- NA; x})
        #
        # evaluate that there are data in each group comparison
        # remove comparisons with missing labeled samples, replace unlabeled samples with global unlabeled average
        good_vals <- match_groups(data, wh, wl, iso_group)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        # WAD correction
        if(correction) {
          shift <- base::Map(function(x, y) sweep(x, 2, y, function(wL, wH) wH - wL),
                             wl, wh)
          shift <- base::lapply(shift, sort, decreasing=FALSE)
          shift <- base::lapply(shift, function(x) median(x[1:floor(offset_taxa * length(x))], na.rm=T))
          wh <- base::lapply(wh, '-', shift)
        }
        # calculate GC content of each taxa
        gc <- base::lapply(wl, function(wL) (1 / 0.083506) * (wL - 1.646057))
        # calculate mol. weight of taxa
        mw_l <- base::lapply(gc, function(GC) (0.496 * GC) + 307.691)
        mw_h <- base::Map(function(x, y) sweep(x, 2, y, function(wL, wH) (((wH - wL)/wL) + 1)),
                          wl, wh)
        mw_h <- base::Map('*', mw_h, mw_l)
        shift <- base::Map(function(x, y) sweep(x, 2, y, function(wL, wH) wH - wL), wl, wh)
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
        good_vals <- match_groups(data, wh, wl, iso_group)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        #
        # WAD correction
        if(correction) {
          shift <- base::Map(function(x, y) sweep(x, 2, y, '-'), wh, wl)
          shift <- base::lapply(shift, sort, decreasing=FALSE)
          shift <- base::lapply(shift, function(x) median(x[1:floor(offset_taxa * length(x))], na.rm=T))
          wh <- base::lapply(wh, '-', shift)
        }
        # calculate GC content of each taxa
        gc <- base::lapply(wl, function(wL) (1 / 0.083506) * (wL - 1.646057))
        # calculate mol. weight of taxa
        mw_l <- base::lapply(gc, function(GC) (0.496 * GC) + 307.691)
        mw_h <- base::Map(function(x, y) sweep(x, 2, y,  function(wH, wL) (((wH - wL)/wL) + 1)),
                          wh, wl)
        mw_h <- base::Map('*', mw_h, mw_l)
        shift <- base::Map(function(x, y) sweep(x, 2, y, '-'), wh, wl)
        # ...........................................
      } else if(separate_light) { # CODE 111
        #
        # evaluate that individual samples align for comparison
        # remove comparisons with missing labeled samples, replace unlabeled samples with average unlabeled
        good_vals <- match_reps(data, wh, wl, iso_group, rep_group=TRUE)
        wh <- good_vals[[1]]
        wl <- good_vals[[2]]
        if(global_light) {
          global_wl <- do.call(rbind, wl)
          global_wl <- colMeans(global_wl, na.rm=TRUE)
          global_wl[is.nan(global_wl)] <- NA
          wl[grepl('light_avg', rownames(wl)),] <- global_wl
        }
        #
        wh <- do.call(rbind, wh)
        wl <- do.call(rbind, wl)
        # WAD correction
        if(correction) {
          shift <- wh - wl
          shift <- sort(shift, decreasing=FALSE)
          shift <- median(shift[1:floor(offset_taxa * length(shift))], na.rm=T)
          wh <- wh - shift
        }
        # calculate GC content of each taxa
        gc <- (1 / 0.083506) * (wl - 1.646057)
        # calculate mol. weight of taxa
        mw_l <- (0.496 * gc) + 307.691
        mw_h <- (((wh - wl)/wl) + 1) * mw_l
        shift <- wh - wl
      }
    }
  }
  # organize and add new data as S4 matrices
  data <- collate_results(data, mw_h, tax_names=tax_names, 'mw_label', sparse=TRUE)
  data <- collate_results(data, mw_l, tax_names=tax_names, 'mw_light', sparse=TRUE)
  # add attributes summarizing the calculation method
  output_attr <- c(rep_group=as.logical(length(data@qsip@rep_group)),
                   sep_label=separate_label,
                   sep_light=separate_light)
  attributes(data@qsip[['mw_label']])$calc_method <- output_attr
  attributes(data@qsip[['mw_light']])$calc_method <- output_attr
  return(data)
}
