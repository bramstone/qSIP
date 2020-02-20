#' Calculation of population dynamics
#'
#' Calculates population growth and death rates
#'
#' @param data Data as a \code{phyloseq} object
#' @param ci_method Character value indicating how to calculate confidence intervals of stable isotope atom excess.
#'   Options are \code{bootstrap} or \code{none} (see \code{details} below for discussion on their differences).
#'   The default is \code{none} indicating that no confidence intervals will be calculated.
#' @param ci Numeric value from 0 to 1 indicating the width of the confidence interval for bootsrapped atom excess values.
#' @param iters Number of (subsampling) iterations to conduct to calculate confidence intervals. Default is \code{999}.
#' @param filter Logical vector specifying whether or not to filter taxa from the weighted average density calculation.
#'   This will require \code{data} to have a filter applied with \code{\link{filter_qsip}}.
#' @param growth_model Character vector specifying whether growth rates should be calculated using exponential or linear growth.
#'   Default is \code{exponential}.
#' @param mu Numeric value from 0 to 1 indicating the expected incorporation of labeled O into DNA during replication.
#'   Default is \code{0.60}.
#' @param correction Logical value indicating whether or not to apply tube-level correction to labeled WAD values.
#' @param offset_taxa Value from 0 to 1 indicating the percentage of the taxa to utilize for calculating offset correction values.
#'   Taxa are ordered by lowest difference in WAD values.
#'   Default is \code{0.1} indicating 10 percent of taxa with the lowest difference in WAD values.
#' @param max_label Numeric value indicating the maximum possible isotope labeling in an experiment.
#'   Keeping the value at \code{1} will ensure that the maximum possible atom excess value of 1 corresponds to complete updake of the isotope.
#'   Recommended for experiments with lower atom percent enrichment treatments (see Details).
#' @param separate_light Logical value indicating whether or not WAD-light scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, unlabeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight in the absence of isotope addition.
#' @param separate_label Logical value indicating whether or not WAD-label scores should be averaged across all replicate groups or not.
#'   If \code{FALSE}, labeled WAD scores across all replicate groups will be averaged, creating a single molecular weight score per taxon
#'   representing it's genetic molecular weight as a result of isotope addition. The default is \code{TRUE}, resulting in no averaging across replicates.
#' @param global_light Logical value indicating whether or not to use WAD-light scores that are global averages (\emph{i.e.,} averaged across
#'   all samples rather than averaged across any specified replicate groups). The default is \code{FALSE}.
#' @param separate_time Logical value indicating whether or not per-capita rates for individual replicates should be calculated using
#'   abundances that match in sample origin (\code{TRUE}) or using abundances that have been averaged across
#'   each group of replicates (\code{FALSE}, the default). Requires that replicate matches have been recorded and specified in the \code{@@rep_num} slot.
#' @param sequential_time Logical value indicating whether or not mortality (turnover) rates should be calculated by comparing the estimated proportion of
#'   remaining unlabeled population size at time \code{t} to the unlabeled population at time \code{t-1} (\code{TRUE}) or to population size at the initial
#'   timepoint (\code{FALSE}).
#' @param rm_light_abund Logical value indicating whether to remove the abundance of taxa from light replicates (\code{TRUE}) and average abundances at
#'   time \emph{t} only from labeled replicates. The alternative (value of \code{FALSE}) is to average abundances at time \emph{t} using all replicates,
#'   labeled and unlabeled which is the default action.
#' @param rel_abund Logical value specifying if relative abundances of taxa are to be calculated prior to calculations. The default is \code{TRUE}.
#'   This parameter is passed to \code{calc_wad}.
#' @param recalc Logical value indicating whether or not to recalculate WAD and molecular weight values or use existing values. Default is \code{TRUE}.
#'   Using bootstrapped calculations will automatically recalculate all values.
#'
#' @details Some details about proper isotope control-treatment factoring and timepoint specification. If weighted average densities or the
#'   change in weighted average densities have not been calculated beforehand, \code{calc_pop} will compute those first.
#'
#'   Timepoint should be in units of days, so that birth will be new 16S copies d-1 and death will be loss (\emph{i.e.,} turnover) of light 16S copies d-1.
#'   Use of different time increments will yield growth rates (e.g. per hour), but must be appropriate for the frequency of sampling.
#'
#'   The \code{recalc} argument is necessary to support the bootstrap subsampling implementation, which re-draws from the table of WAD
#'   values to create a bootstrap resampled WAD table. As such, automatic recalculation of WAD values is inappropriate within each
#'   bootstrap iteration. Typically, users should not set \code{recalc} to \code{FALSE}. Note that the true, observed WAD and
#'   molecular weight values will be returned following completion of bootstrapping.
#'
#'   Setting \code{max_label < 1} will return birth and death values higher than would otherwise be returned. A \code{max_label} value of 1 indicates
#'   that the maximum molecular weight of an organism respresents its molecular weight under the assumption of complete, or 100\%, isotope incorporation.
#'   For various reasons, complete isotope incorporation will be impossible. However, a \code{max_label} value less than 1 will indicate molecular weight
#'   where 0 indicates no isotope incorporation and 1 indicates the highest possible incorporation, as constrained by atom percent enrichment
#'   provided in the experiment. For example, an experiment enriching soil with 18-O at 50\% atom enrichment will want to specify \code{max_label=0.5}
#'   and an atom excess fraction value of 1 in this case corresponds to an organism that has succeeded in incorporating 18-O into it's nucleic acids
#'   at 50\%.
#'
#' @return \code{calc_pop} adds two S4 Matrix class objects (which more efficiently stores sparse matrix data) to the \code{data@@qsip@@.Data} slot
#'   of population birth rates for each taxon at each group of replicates. The row and column specifications will mirror those of the \code{phylosip}'s
#'   \code{\link{otu_table}}, meaning if taxa are listed on the table rows, they will in the resulting S4 Matrix class.
#'
#'   Note that the bootstrap method produces a \emph{single} bootstrapped median (and matching confidence intervals) for groups of replicates,
#'   either grouped by isotope treatment alone, or also by some other grouping factor (if \code{data@@qsip@@rep_group} is specified).
#'   Using no bootstrap value allows separate enrichment values to be attained for each replicate, if \code{separate_label=TRUE}.
#'
#'
#' @seealso \code{\link{calc_mw}}
#'
#' @examples
#'  # Load in example data
#'
#'  # Calculate population fluxes
#'
#' @references
#'  Hungate, Bruce, \emph{et al.} 2015. Quantitative microbial ecology through stable isotope probing.
#'  \emph{Applied and Environmental Microbiology} \strong{81}: 7570 - 7581.
#'
#'  Koch, Benjamin, \emph{et al.} 2018. Estimating taxon-specific population dynamics in diverse microbial communities.
#'  \emph{Ecosphere} \strong{9}: e02090.
#'
#' @export

# NOTE: MAX_LABEL IS NOT CURRENTLY IMPLEMENTED IN THE CALCULATIONS. NEED TO LOOK AT BEST WAY TO DO THIS.
calc_pop <- function(data, ci_method=c('none', 'bootstrap'), ci=.95, iters=999, filter=FALSE, growth_model=c('exponential', 'linear'),
                     mu=0.6, correction=FALSE, offset_taxa=0.1, max_label=1, separate_label=FALSE, separate_light=TRUE, global_light=FALSE,
                     separate_time=FALSE, sequential_time=FALSE, rm_light_abund=FALSE, rel_abund=TRUE, recalc=TRUE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object', call.=FALSE)
  ci_method <- match.arg(tolower(ci_method))
  growth_model <- match.arg(tolower(growth_model))
  if(data@qsip@iso!='18O') stop('Must use 18O-labeled treatment to calculate population change', call.=FALSE)
  if(length(data@qsip@timepoint)==0) stop('Must specify different sample times with timepoint', call.=FALSE)
  times <- data@sam_data[[data@qsip@timepoint]]
  times <- times[!is.na(times)]
  if(nlevels(times)==1 || length(unique(times))==1) stop('Only one timepoint present in the data - cannot calculate population change', call.=FALSE)
  rm(times)
  if(separate_time && length(data@qsip@rep_num)==0) stop('Must specify replicate (sample origin) matches with rep_num', call.=FALSE)
  #
  # -------------------------------------------------------------
  # no CI and resampling
  #
  if(ci_method=='none') {
    # if recalculation wanted, do that
    # this will also handle rep_id validity (through calc_wad) and rep_group/iso_trt validity (through calc_d_wad)
    if(recalc | is.null(data@qsip[['mw_label']])) {
      data <- suppressWarnings(calc_mw(data,
                                       filter=filter,
                                       correction=correction,
                                       offset_taxa=offset_taxa,
                                       separate_light=separate_light,
                                       separate_label=separate_label,
                                       global_light=global_light,
                                       rel_abund=rel_abund,
                                       recalc=TRUE))
    }
    # transform sequencing abundances to 16S copy numbers
    # returns feature table (as matrix) with taxa as columns, samples as rows
    nt <- copy_no(data, rel_abund=rel_abund)
    n_taxa <- ncol(nt)
    if(filter) {
      tax_names <- data@qsip@filter
      nt <- nt[,colnames(nt) %in% tax_names]
    }
    tax_names <- colnames(nt)
    # calculate per-taxon total 16S copy abundance for each sample (i.e., sum over fractions)
    nt <- split_data(data, nt, data@qsip@rep_id)
    nt <- lapply(nt, colSums, na.rm=T)
    nt <- do.call(rbind, nt)
    # remove samples with NA for rep_id, timepoint, or rep_group
    nt <- valid_samples(data, nt, 'time')
    time_group <- nt[[2]]; nt <- nt[[1]]
    # identify labeling in samples
    iso_reps <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    time_group <- merge(time_group, iso_reps[,c('replicate', 'iso')], all.x=TRUE)
    time_group <- time_group[match(rownames(nt), time_group$replicate),]
    # identify lowest timepoint (doesn't have to be 0) and separate by that
    nt_0 <- nt[as.numeric(time_group$time)==1,]
    nt_t <- nt[as.numeric(time_group$time) > 1,]
    time_0 <- droplevels(time_group[as.numeric(time_group$time)==1,])
    time_t <- droplevels(time_group[as.numeric(time_group$time) > 1,])
    #
    # extract MW-labeled and convert to S3 matrix with taxa as columns
    mw_h <- data@qsip[['mw_label']]
    mw_l <- data@qsip[['mw_light']]
    if(phyloseq::taxa_are_rows(data)) {
      if(is.matrix(mw_h)) mw_h <- t(mw_h)
      if(is.matrix(mw_l)) mw_l <- t(mw_l)
    }
    if(is.matrix(mw_h)) tax_names <- colnames(mw_h) else tax_names <- names(mw_h)
    # create MW heavy max
    mw_max <- (12.07747 * mu) + mw_l
    #
    # ----------------
    # Different methods of grouping data produce slightly different calculations
    # These are coded as a 5-digit binary code based on 3 grouping criteria
    #   - Replicate separated by group:  0/1 (no/yes)
    #   - Labeled densities separated:   0/1
    #   - Unlabeled densities separated: 0/1
    #   - Time=0 abundances separated:   0/1
    #   - Sequential timepoint calcs:    0/1 (no - calc from time=0 / yes)
    # Example: 10101 means calculate pop flux with seprate groups, separate labeled values, averaged unlabeled values,
    #   separate time 0 abundances, and pop rates from the second timepoint will be calculated using abundances from the first timepoint
    # ----------------
    #
    # if there is no replicate grouping
    if(length(data@qsip@rep_group)==0) {
      # w/o replicate grouping, only split is based on timepoint
      unique_time_t <- unique(time_t[,c('time', 'interaction')])
      #
      if(!separate_label) {
        if(!separate_light) {
          if(!separate_time) {  # CODE 0000*
            #
            # average all abundances
            nt_0 <- colMeans(nt_0, na.rm=T)
            if(rm_light_abund) {
              nt_t <- nt_t[as.numeric(time_t$iso)==2,]
              time_t <- time_t[as.numeric(time_t$iso)==2,]
            }
            nt_t <- split_data(data, nt_t, time_t$time, grouping_w_phylosip=FALSE)
            nt_t <- base::lapply(nt_t, colMeans, na.rm=T)
            nt_t <- do.call(rbind, nt_t)
            #
            # length of incubation times?
            incubate_t <- unique_time_t$time
            incubate_t <- as.numeric(as.character(incubate_t))
            #
            # calculate unlabeled abundances at time t
            nl_t <- nt_t * ((mw_max - mw_h) / (mw_max - mw_l))
            #
            # add different timepoints to time 0 for sequential death rates or just repeat time 0?
            if(sequential_time) { # CODE 00001
              nt_0 <- rbind(nt_0, nl_t[-nrow(nl_t),])
            } else { # CODE 00000
              nt_00 <- vector('list', nrow(nl_t))
              nt_00 <- base::lapply(nt_00, function(x) nt_0)
              nt_0 <- do.call(rbind, nt_00)
            }
            #
          } else if(separate_time) {  # CODE 0001*
            #
            # remove light abundances or average light and labeled abundances?
            if(rm_light_abund) {
              nt_t <- nt_t[as.numeric(time_t$iso)==2,]
              time_t <- time_t[as.numeric(time_t$iso)==2,]
            } else {
              reps <- unique(as(data@sam_data[,c(data@qsip@rep_id, data@qsip@rep_num)], 'data.frame'))
              names(reps) <- c('replicate', 'replicate_num')
              time_t <- merge(time_t, reps, all.x=TRUE)
              time_t <- time_t[match(rownames(nt_t), time_t$replicate),]
              nt_t <- split_data(data, nt_t, time_t$rep_num, grouping_w_phylosip=FALSE, keep_names=1)
              nt_t <- base::lapply(nt_t, colMeans, na.rm=TRUE)
              nt_t <- do.call(rbind, nt_t)
            }
            # split time t abundances, and MWs by timepoint
            nt_t <- split_data(data, nt_t, time_t$time, grouping_w_phylosip=FALSE, keep_names=1)
            mw_h <- split_data(data, mw_h, unique_time_t$time, grouping_w_phylosip=FALSE)
            mw_l <- split_data(data, mw_l, unique_time_t$time, grouping_w_phylosip=FALSE)
            mw_max <- split_data(data, mw_max, unique_time_t$time, grouping_w_phylosip=FALSE)
            #
            # length of incubation times?
            incubate_t <- time_t$time
            incubate_t <- as.numeric(as.character(incubate_t))
            #
            # calculate unlabeled abundances at time t
            enrich <- base::Map(function(mwH, mwL, mwMax) (mwMax - mwH) / (mwMax - mwL),
                                mw_h, mw_l, mw_max)
            nl_t <- base::Map(function(Nt, enrich) sweep(Nt, 2, enrich, '*'), nt_t, enrich)
            #
            # add different timepoints to time 0 for sequential death rates or just repeat time 0?
            if(sequential_time) { # CODE 00011
              nt_0 <- rbind(nt_0, nl_t[[-length(nl_t)]])
            } else { # CODE 00010
              nt_00 <- vector('list', length(nl_t))
              nt_00 <- base::lapply(nt_00, function(x) nt_0)
              nt_0 <- do.call(rbind, nt_00)
            }
            #
            nl_t <- do.call(rbind, nl_t)
            nt_t <- do.call(rbind, nt_t)
          }
          #
        } else if(separate_light) {
          if(!separate_time) {  # CODE 0001*
            #
            # average all abundances
            nt_0 <- colMeans(nt_0, na.rm=T)
            if(rm_light_abund) {
              nt_t <- nt_t[as.numeric(time_t$iso)==2,]
              time_t2 <- time_t[as.numeric(time_t$iso)==2,]
            }
            time_t1 <- time_t[as.numeric(time_t$iso)==1,]
            # split time t abundances, and MWs by timepoint
            nt_t <- split_data(data, nt_t, time_t2$time, grouping_w_phylosip=FALSE)
            nt_t <- base::lapply(nt_t, colMeans, na.rm=T)
            mw_h <- split_data(data, mw_h, time_t1$time, grouping_w_phylosip=FALSE)
            mw_l <- split_data(data, mw_l, time_t1$time, grouping_w_phylosip=FALSE)
            mw_max <- split_data(data, mw_max, time_t1$time, grouping_w_phylosip=FALSE)
            #
            # length of incubation times (focusing only on 16O)
            incubate_t <- time_t1$time
            incubate_t <- as.numeric(as.character(incubate_t))
            #
            # calculate unlabeled abundances at time t
            enrich <- base::Map(function(mwH, mwL, mwMax) (mwMax - mwH) / (mwMax - mwL),
                                mw_h, mw_l, mw_max)
            nl_t <- base::Map(function(Nt, enrich) sweep(enrich, 2, Nt, '*'), nt_t, enrich)
            #
            # create abundance matrix of same dimensions as nl_t
            mat_0 <- base::lapply(nl_t, function(x) {x[] <- 0; x})
            mat_nt <- base::Map(function(mat_0, Nt) sweep(mat_0, 2, Nt, '+'), mat_0, nt_t)
            mat_n0 <- base::lapply(mat_0, function(x) sweep(x, 2, nt_0, '+'))
            #
            nl_t <- do.call(rbind, nl_t)
            nt_t <- do.call(rbind, mat_nt)
            nt_0 <- do.call(rbind, mat_n0)
            rm(mat_0)
            #
            # add different timepoints to time 0 for sequential death rates or just repeat time 0?
            if(sequential_time) { # CODE 00101
              nt_0 <- nt_0[1:nrow(time_0),]
              all_but_last_t <- as.numeric(time_t1$time) < max(as.numeric(time_t1$time))
              nt_0 <- rbind(nt_0, nl_t[all_but_last_t,])
            } # else: CODE 00100
            #
          } else if(separate_time) {  # CODE 0011*
            #
            # remove light abundances or average light and labeled abundances?
            if(rm_light_abund) {
              nt_t <- nt_t[as.numeric(time_t$iso)==2,]
              time_t <- time_t[as.numeric(time_t$iso)==1,]
            } else {
              reps <- unique(as(data@sam_data[,c(data@qsip@rep_id, data@qsip@rep_num)], 'data.frame'))
              names(reps) <- c('replicate', 'rep_num')
              time_t <- merge(time_t, reps, all.x=TRUE)
              time_t <- time_t[match(rownames(nt_t), time_t$replicate),]
              nt_t <- split_data(data, nt_t, interaction(time_t$rep_num, time_t$time), grouping_w_phylosip=FALSE, keep_names=1)
              nt_t <- base::lapply(nt_t, colMeans, na.rm=TRUE)
              nt_t <- do.call(rbind, nt_t)
              time_t <- time_t[as.numeric(time_t$iso)==1,]
            }
            # split time t abundances, and MWs by timepoint
            nt_t <- split_data(data, nt_t, time_t$time, grouping_w_phylosip=FALSE)
            mw_h <- split_data(data, mw_h, time_t$time, grouping_w_phylosip=FALSE, keep_names=1)
            mw_l <- split_data(data, mw_l, time_t$time, grouping_w_phylosip=FALSE)
            mw_max <- split_data(data, mw_max, time_t$time, grouping_w_phylosip=FALSE)
            #
            # length of incubation times?
            incubate_t <- time_t$time
            incubate_t <- as.numeric(as.character(incubate_t))
            #
            # calculate unlabeled abundances at time t
            enrich <- base::Map(function(mwH, mwL, mwMax) (mwMax - mwH) / (mwMax - mwL),
                                mw_h, mw_l, mw_max)
            nl_t <- base::Map('*', enrich, nt_t)
            #
            # add different timepoints to time 0 for sequential death rates or just repeat time 0?
            if(sequential_time) { # CODE 00111
              nt_0 <- rbind(nt_0, nl_t[[-length(nl_t)]])
            } else { # CODE 00110
              nt_00 <- vector('list', length(nl_t))
              nt_00 <- base::lapply(nt_00, function(x) nt_0)
              nt_0 <- do.call(rbind, nt_00)
            }
            #
            nl_t <- do.call(rbind, nl_t)
            nt_t <- do.call(rbind, nt_t)
            #
          }
        }
        #
      } else if(separate_label) {
        if(!separate_light) {
          if(!separate_time) {  # CODE 0100*
          } else if(separate_time) {  # CODE 0101*
          }
        } else if(separate_light) {
          if(!separate_time) {  # CODE 0110*
          } else if(separate_time) {  # CODE 0111*
          }
        }
      }
    } else if(length(data@qsip@rep_group)==1) {
      if(!separate_label) {
        if(!separate_light) {
          if(!separate_time) {  # CODE 1000*
          } else if(separate_time) {  # CODE 1001*
          }
        } else if(separate_light) {
          if(!separate_time) {  # CODE 1010*
          } else if(separate_time) {  # CODE 1011*
          }
        }
      } else if(separate_label) {
        if(!separate_light) {
          if(!separate_time) {  # CODE 1100*
          } else if(separate_time) {  # CODE 1101*
          }
        } else if(separate_light) {
          if(!separate_time) {  # CODE 1110*
          } else if(separate_time) {  # CODE 1111*
          }
        }
      }
    }
    #
    # Abundance matrices should be same dimensions
    if(!isTRUE(all.equal(dim(nt_0), dim(nl_t)))) {
      stop('Time 0 and time t abundance data should be same size for groupings specified but are not',
           '\nCheck that @rep_num is correct or if @rep_group should be specified', call.=FALSE)
    }
    #
    # calculate birth and death rates
    if(growth_model=='exponential') {
      b <- log(nt_t / nl_t) * (1/incubate_t)
      d <- log(nl_t / nt_0) * (1/incubate_t)
    } else if(growth_model=='linear') {
      b <- (nt_t - nl_t) * (1/incubate_t)
      d <- (nl_t - nt_0) * (1/incubate_t)
    }



    # separate samples based on timepoint, keeping only valid samples
    # if matching replicates, re-order and match according to replicate number
    if(separate_time && length(data@qsip@rep_num)==1) {
      ft <- valid_samples(data, ft, 'time')
      time_group <- ft[[2]]; ft <- ft[[1]]
    } else{
      ft <- valid_samples(data, ft, 'time')
      time_group <- ft[[2]]; ft <- ft[[1]]
    }
    # remove light samples from abundance calcs
    if(rm_light_abund) {
      iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
      light_group <- iso_group[as.numeric(iso_group$iso)==1,]
      ft <- ft[!rownames(ft) %in% light_group$replicate,]
    }
    # separate by time or time:group interaction
    time_group <- time_group[match(rownames(ft), time_group$replicate),]
    ft <- split_data(data, ft, time_group$interaction, grouping_w_phylosip=FALSE, keep_names=1)
    time_group2 <- unique(time_group[,!names(time_group) %in% c('replicate', 'replicate_num')]) # only get unique elements to match levels in ft
    # if keeping labeled 16S abundances and MWs separate ft will be in a list
    if(separate_label) {
      ft <- ft[match(time_group2$interaction, names(ft))]
      nt_0 <- ft[as.numeric(time_group2$time)==1]
      # average time 0 only if not matching abundances to each replicate across incubation
      if(!separate_time) {
        nt_0 <- base::lapply(nt_0, colMeans, na.rm=TRUE)
      }
      ft[as.numeric(time_group2$time)==1] <- nt_0
      ft <- base::lapply(ft, function(x) {x[x==0] <- NA; x})
    # if averaging labeled 16S abundances and MWs ft will be in matrix
    } else {
      ft <- base::lapply(ft, colMeans, na.rm=T)
      ft <- do.call(cbind, ft)
      ft[ft==0] <- NA
      ft <- ft[,match(time_group2$interaction, colnames(ft))] # re-order columns to match time_group2$interaction
    }
    # Does a time 0 exist in these data?
    if(!any(time_group2$time==0)) {
      warning('No timepoints designated as time 0; using time ',
              levels(time_group2$time)[1],
              ' as time before isotope addition', call.=FALSE)
    }
    # get 16S copy numbers for different timepoints
    n_t_names <- paste0('n_t_',levels(time_group2$time))
    if(separate_label) { # get n_t values from list
      for(time in 1:nlevels(time_group2$time)) {
        nt_t <- ft[time_group2$time==levels(time_group2$time)[time]]
        assign(n_t_names[time], t(do.call(rbind, nt_t)))
      }; rm(nt_t)
    } else { # get n_t values from matrix
      for(time in 1:nlevels(time_group2$time)) {
        assign(n_t_names[time], ft[,time_group2$time==levels(time_group2$time)[time]])
      }
    }
    # extract MW-labeled and convert to S3 matrix with taxa as ROWS (opposite all other calcs)
    mw_h <- data@qsip[['mw_label']]
    mw_l <- data@qsip[['mw_light']]
    if(!phyloseq::taxa_are_rows(data)) mw_h <- t(mw_h)
    if(!phyloseq::taxa_are_rows(data) && is.matrix(mw_l)) mw_l <- t(mw_l)
    # calculate proportion in light fraction (N_light) at any time after 0
    n_l_names <- paste0('n_l_', levels(time_group2$time))
    for(time in 2:nlevels(time_group2$time)) {
      if(length(data@qsip@rep_group)==0) {
        # if there is no grouping separate from timepoints, go by the column timepoints
        mw_h_t <- mw_h[,time - 1]
        if(is.null(dim(mw_l))) mw_l_t <- mw_l else mw_l_t <- as.matrix(mw_l)[,time - 1]
        # else break mw_label and mw_light into list separated by timepoint
      } else {
        if(separate_label) {
          time_group_t <- time_group[as.numeric(time_group$time) > 1,]
          time_group_t <- time_group_t[match(colnames(mw_h), time_group_t$replicate),]
        }
        else if(!separate_label){
          time_group_t <- time_group2[as.numeric(time_group2$time) > 1,]
        }
        time_group_t$time <- factor(time_group_t$time)
        mw_h_t <- split_data(data, t(mw_h), time_group_t$time, grouping_w_phylosip=FALSE, keep_names=1)
        mw_h_t <- t(mw_h_t[[time - 1]])
        if(is.null(dim(mw_l))) {
          mw_l_t <- mw_l
        } else {
            mw_l_t <- suppressWarnings(split_data(data, t(mw_l), time_group_t$time, grouping_w_phylosip=FALSE))
            mw_l_t <- t(mw_l_t[[time - 1]])
        }
      }
      # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
      mw_max <- (12.07747 * mu) + mw_l_t
      if(separate_label)  mw_h_t <- mw_h_t[,match(colnames(get(n_t_names[time])), colnames(mw_h_t))]
      # calculate abundances
      if(isTRUE(all.equal(dim(mw_max), dim(mw_h_t)))) {
        n <- ((mw_max - mw_h_t)/(mw_max - mw_l_t)) * get(n_t_names[time])
      } else {
        if(is.null(dim(mw_h_t))) mw_h_t <- as.matrix(mw_h_t)
        if(is.null(dim(mw_l_t))) mw_l_t <- as.matrix(mw_l_t)
        num <- sweep(mw_h_t, 1, mw_max) * -1 # mw_max - mw_h
        denom <- sweep(mw_l_t, 1, mw_max) * -1 # mw_max - mw_l
        n <- sweep(num, 1, denom, '/') * get(n_t_names[time]) # MW_proportion * N_t
      }
      if(!separate_label) colnames(n) <- colnames(get(n_t_names[time]))
      # remove abundances less than 0 (occurs when labeled MWs are heavier than heavymax)
      n[n < 0] <- NA
      assign(n_l_names[time], n)
    }; suppressWarnings(rm(n, mw_h_t, mw_l_t, time_group_t, mw_max))
    # calculate birth and death rate for each timepoint after 0
    b_names <- paste0('b_', levels(time_group2$time))
    d_names <- paste0('d_', levels(time_group2$time))
    for(time in 2:nlevels(time_group2$time)) {
      incubate_time <- as.numeric(levels(time_group2$time)[time])
      # n_col_t <- ncol(get(n_t_names[time]))
      # n_col_0 <- ncol(get(n_t_names[1]))
      # if match_replicate=F, separate labeled samples which will have different dimensions than unlabeled
      if(!match_replicate) {
      # if(n_col_t > n_col_0) {
        sam_names_t <- colnames(get(n_t_names[time]))
        group_repeat <- time_group[match(sam_names_t, time_group$replicate), 'grouping']
        time_group_0 <- time_group2[time_group2$time==levels(time_group2$time)[1],]
        group_repeat <- as.character(time_group_0[match(group_repeat, time_group2$grouping), 'interaction'])
        assign(n_t_names[1], get(n_t_names[1])[,group_repeat])
      }
      if(growth_model=='exponential') {
        b <- get(n_t_names[time]) / get(n_l_names[time])
        d <- get(n_l_names[time]) / get(n_t_names[1])
        b <- log(b) / incubate_time
        d <- log(d) / incubate_time
      } else {
        b <- get(n_t_names[time]) - get(n_l_names[time])
        d <- get(n_l_names[time]) - get(n_t_names[1])
        b <- (b / get(n_t_names[1])) / incubate_time
        d <- (d / get(n_t_names[1])) / incubate_time
      }
      colnames(b) <- colnames(d) <- colnames(get(n_t_names[time]))
      assign(b_names[time], b)
      assign(d_names[time], d)
    }; suppressWarnings(rm(b, d, n_col_t, n_col_0, sam_names_t, group_repeat, time_group_0))
    # if more than two timepoints (0, and t), combine resulting matrices, but skip time 0
    if(length(n_t_names) > 2) {
      b <- t(do.call(cbind, mget(b_names[2:length(b_names)])))
      d <- t(do.call(cbind, mget(d_names[2:length(d_names)])))
    } else {
      b <- get(b_names[2])
      d <- get(d_names[2])
    }
    # organize and add new data as S4 matrices
    data <- collate_results(data, t(b), tax_names=tax_names, 'birth_rate', sparse=T)
    data <- collate_results(data, t(d), tax_names=tax_names, 'death_rate', sparse=T)
    return(data)
  #
  # -------------------------------------------------------------
  # CI values obtained through bootstrap subsampling (will need to bootstrap abundances AND WADs/MWs)
  #
  } else if(ci_method=='bootstrap') {
    # calculate and create subsampling criteria for WADs......................................
    data <- suppressWarnings(calc_wad(data, filter=filter, rel_abund=rel_abund))
    wads <- as(data@qsip[['wad']], 'matrix')
    if(phyloseq::taxa_are_rows(data)) wads <- t(wads)
    n_taxa <- ncol(wads)
    tax_names <- colnames(wads)
    # keep only valid samples
    wads <- valid_samples(data, wads, 'iso', quiet=TRUE)
    iso_group <- wads[[2]]; wads <- wads[[1]]
    # split by replicate groups
    sam_names_wads <- rownames(wads)
    wads <- split_data(data, wads, iso_group$interaction, grouping_w_phylosip=FALSE)
    # how many samples in each group to subsample WADS with?
    subsample_n <- base::lapply(wads, nrow)
    subsample_wads <- base::lapply(subsample_n,
                              function(x) sample.int(x, size=iters*x, replace=TRUE))
    subsample_wads <- base::mapply(matrix,
                                   subsample_wads,
                                   nrow=subsample_n,
                                   byrow=F, SIMPLIFY=FALSE)
    # calculate and create subsampling criteria for 16S gene copy number......................
    # transform sequencing abundances to 16S copy numbers (taxa as columns)
    ft <- copy_no(data)
    if(filter) ft <- ft[,colnames(ft) %in% tax_names]
    # calculate per-taxon 16S copy abundance for each sample
    ft <- split_data(data, ft, data@qsip@rep_id)
    ft <- lapply(ft, colSums, na.rm=T)
    ft <- do.call(rbind, ft)
    # separate abundances based on timepoint, keeping only valid samples
    ft <- valid_samples(data, ft, 'time')
    time_group <- ft[[2]]; ft <- ft[[1]]
    sam_names <- time_group$replicate
    # remove light samples from abundance calcs
    if(rm_light_abund) {
      iso_group_ft <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
      light_group <- iso_group_ft[as.numeric(iso_group_ft$iso)==1,]
      ft <- ft[!rownames(ft) %in% light_group$replicate,]
    }
    # separate by time or time:group interaction
    time_group <- time_group[match(rownames(ft), time_group$replicate),]
    time_group <- droplevels(time_group)
    ft <- split_data(data, ft, time_group$interaction, grouping_w_phylosip=F)
    # how many samples in each group to subsample with?
    subsample_n <- base::lapply(ft, nrow)
    subsample <- base::lapply(subsample_n,
                              function(x) sample.int(x, size=iters*x, replace=TRUE))
    subsample <- base::mapply(matrix,
                              subsample,
                              nrow=subsample_n,
                              byrow=F, SIMPLIFY=FALSE)
    # collect output in matrices (each column is a pop matrix from that iterations' subsampling)
    if(length(data@qsip@rep_group)==0) {
      n_groups <- 1
      boot_rnames <-  expand.grid(tax_names,
                                  levels(time_group$time)[2:nlevels(time_group$time)],
                                  stringsAsFactors=FALSE)
      boot_rnames <- interaction(boot_rnames[,1], boot_rnames[,2], sep=':')
    } else {
      n_groups <- nlevels(time_group$grouping)
      boot_rnames <- expand.grid(tax_names,
                                 levels(time_group$grouping),
                                 levels(time_group$time)[2:nlevels(time_group$time)],
                                 stringsAsFactors=FALSE)
      boot_rnames <- interaction(boot_rnames[,1], boot_rnames[,2], boot_rnames[,3], sep=':')
    }
    n_timepoints <- nlevels(time_group$time) - 1
    boot_collect_b <- matrix(0, nrow=n_taxa * n_groups * n_timepoints, ncol=iters)
    rownames(boot_collect_b) <- boot_rnames
    boot_collect_d <- boot_collect_b
    rm(boot_rnames)
    #
    for(i in 1:iters) {
      # subsample WADs
      subsample_i_wads <- lapply(subsample_wads, function(x) x[,i])
      wads_i <- mapply(function(x, y) x[y,,drop=FALSE], wads, subsample_i_wads, SIMPLIFY=FALSE)
      wads_i <- recombine_in_order(wads_i, iso_group, n_taxa)
      rownames(wads_i) <- sam_names_wads
      # calc diff_WADs, MWs, and N values
      data <- suppressWarnings(collate_results(data, wads_i, tax_names=tax_names, 'wad', sparse=TRUE))
      data <- suppressWarnings(calc_mw(data,
                                       separate_label=FALSE,
                                       separate_light=FALSE,
                                       recalc=FALSE))
      mw_h <- data@qsip[['mw_label']]
      mw_h <- as(mw_h, 'matrix')
      mw_l <- data@qsip[['mw_light']]
      if(!is.null(dim(mw_l))) mw_l <- as(mw_l, 'matrix')
      if(!phyloseq::taxa_are_rows(data)) mw_h <- t(mw_h)
      if(!phyloseq::taxa_are_rows(data) && is.matrix(mw_l)) mw_l <- t(mw_l)
      # subsample abundances
      # calculate per-taxon average 16S copy abundance for each group:time interaction point
      subsample_i <- lapply(subsample, function(x) x[,i])
      ft_i <- mapply(function(x, y) x[y,,drop=FALSE], ft, subsample_i, SIMPLIFY=FALSE)
      ft_i <- lapply(ft_i, colMeans, na.rm=T)
      ft_i <- t(recombine_in_order(ft_i, time_group, n_taxa, condensed_grouping=TRUE))
      colnames(ft_i) <- time_group$interaction[!duplicated(time_group$interaction)]
      ft_i[ft_i==0] <- NA
      # get per-taxon 16S copy numbers for different timepoints
      time_group2 <- unique(time_group[,!names(time_group) %in% 'replicate']) # only get unique elements to match levels in ft
      ft_i <- ft_i[,match(time_group2$interaction, colnames(ft_i))] # re-order columns to match time_group2$interaction
      n_t_names <- paste0('n_t_',levels(time_group2$time))
      # t represents different timepoints
      for(t in 1:nlevels(time_group2$time)) {
        assign(n_t_names[t],
               ft_i[,time_group2$time==levels(time_group2$time)[t]])
      }
      # calculate proportion in light fraction (N_light) at any time after 0
      n_l_names <- paste0('n_l_', levels(time_group2$time))
      for(time in 2:nlevels(time_group2$time)) {
        if(length(data@qsip@rep_group)==0) {
          # if there is no grouping separate from timepoints, go by the column timepoints
          mw_h_t <- mw_h[,time - 1]
          mw_l_t <- as.matrix(mw_l)[,time - 1]
          # else break mw_label and mw_light into list separated by timepoint
        } else {
          time_group_t <- time_group2[as.numeric(time_group2$time) > 1,]
          time_group_t$time <- factor(time_group_t$time)
          mw_h_t <- split_data(data, t(mw_h), time_group_t$time, grouping_w_phylosip=FALSE)
          mw_l_t <- suppressWarnings(split_data(data, t(mw_l), time_group_t$time, grouping_w_phylosip=FALSE))
          mw_h_t <- t(mw_h_t[[time - 1]])
          mw_l_t <- t(mw_l_t[[time - 1]])
        }
        # calculate mol. weight heavy max (i.e., what is maximum possible labeling)
        mw_max <- (12.07747 * mu) + mw_l_t
        if(all(dim(mw_max)==dim(mw_h_t))) {
          n <- ((mw_max - mw_h_t)/(mw_max - mw_l_t)) * get(n_t_names[time])
        } else {
          num <- sweep(mw_h_t, 1, mw_max) * -1
          denom <- sweep(mw_l_t, 1, mw_max) * -1
          n <- sweep(num, 1, denom, '/') * get(n_t_names[time])
        }
        colnames(n) <- colnames(get(n_t_names[time]))
        # remove abundances less than 0 (occurs when labeled MWs are heavier than heavymax)
        n[n < 0] <- NA
        assign(n_l_names[time], n)
      }; suppressWarnings(rm(n, mw_h_t, mw_l_t, time_group_t, mw_max))
      # calculate birth and death rate for each timepoint after 0
      b_names <- paste0('b_', levels(time_group2$time))
      d_names <- paste0('d_', levels(time_group2$time))
      for(time in 2:nlevels(time_group2$time)) {
        incubate_time <- as.numeric(levels(time_group2$time)[time])
        if(growth_model=='exponential') {
          b <- get(n_t_names[time]) / get(n_l_names[time])
          d <- get(n_l_names[time]) / get(n_t_names[1])
          b <- log(b) / incubate_time
          d <- log(d) / incubate_time
        } else {
          b <- get(n_t_names[time]) - get(n_l_names[time])
          d <- get(n_l_names[time]) - get(n_t_names[1])
          b <- (b / get(n_t_names[1])) / incubate_time
          d <- (d / get(n_t_names[1])) / incubate_time
        }
        assign(b_names[time], b)
        assign(d_names[time], d)
      }; rm(b,d)
      # if more than two timepoints (0, and t), combine resulting matrices, but skip time 0
      if(length(n_t_names) > 2) {
        b <- do.call(cbind, mget(b_names[2:length(b_names)]))
        d <- do.call(cbind, mget(d_names[2:length(d_names)]))
      } else {
        b <- get(b_names[2])
        d <- get(d_names[2])
      }
      # organize and add data as single columns in bootstrap output matrices
      boot_collect_b[,i] <- c(b)
      boot_collect_d[,i] <- c(d)
    }
    # END OF BOOTSTRAP ITERATIONS
    #
    # clean workspace
    suppressWarnings(rm(ft_i, wads_i, subsample, subsample_wads,
       subsample_n, subsample_i, subsample_i_wads, b, d,
       b_names, d_names, n_l_names, n_t_names, ft, wads, mw_lab, mw_l))
    # summarize birth, death, flux across iterations (lower CI, median, upper CI)
    ci_birth <- summarize_ci(boot_collect_b, ci,
                             grouping=time_group,
                             ncols=n_taxa,
                             data=data)
    ci_death <- summarize_ci(boot_collect_d, ci,
                             grouping=time_group,
                             ncols=n_taxa,
                             data=data)
    rm(boot_collect_b, boot_collect_d)
    # collate results
    objects <- c('ci_birth', 'ci_death')
    metric <- c('birth_rate', 'death_rate')
    ci_level <- c('ci_l', '', 'ci_u')
    for(i in 1:2) {
      for(j in 1:2) {
        data <- collate_results(data,
                                get(objects[i])[[j]],
                                tax_names=tax_names,
                                metric=sub('_{1}$', '', paste(metric[i], ci_level[j], sep='_')),
                                sparse=TRUE)
      }
    }
    # recalculate WAD, diff_WAD, and MW values (they've been replaced by bootstrapped versions)
    data <- suppressWarnings(calc_mw(data,
                                     filter=filter,
                                     correction=correction,
                                     offset_taxa=offset_taxa,
                                     separate_light=separate_light,
                                     separate_label=separate_label,
                                     recalc=TRUE))
    return(data)
  }
}
