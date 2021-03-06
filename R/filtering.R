#' Specify multiple taxa filting options for qSIP calculations
#'
#' Creates a data frame of per-replicate and per-fraction minimum incidence frequencies to retain microbial taxa in the dataset
#'
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param rm_combn Optional character vector specifying any particular combinations of replicate and fraction frequency to remove.
#'   Replicate and frequency combinations should be specified by separation with \code{:} (\emph{e.g.}, \code{'3:12'})
#' @param hard Single-length numeric value indicating which filtering option is to be made a hard filtering criteria.
#'   Valid values are \code{1:nrow(data)}
#' @param soft Single-length numeric value indicating which filtering option is to be made a soft filtering criteria.
#'   Valid values are \code{1:nrow(data)}
#' @param rm_low_freq Logical values indicating whether or not a given filtering specification should trust replicates where a taxon
#'   fails to meet the minimum filtering threshold, even if that taxon occurs frequently in other replicates in a group
#'
#' @details \code{create_filters} should be used to specify the desired frequencies of microbial taxa \emph{prior} to calculation of
#'   atom excess fraction / atom percent excess. Filters must also be specified so that diagnostic functions can be run as well.
#'   Because \code{create_filters} uses \code{\link[base]{expand.grid}}, all possible combinations of the first two arguments will
#'   be used to construct the resulting data frame If certain combinations are not desirable, the \code{rm_combn} argument can be
#'   supplied to remove them.
#'
#'   Hard and soft filtering parameters indicate whether taxa should be removed completely if they fail to satisfy the minimum
#'   frequency threshold across all cases (hard) or if they should be removed only from invdividual comparisons if they fail to
#'   meet the frequency threshold in that comparison. For example, in a data set with two replicate groups (representing some two
#'   distinct biological or ecological units), a taxon meets the minimum frequency requirment in group one (meaning it occurs in both
#'   the unlabeled and labeled treatments of group one) but not group two (it is not frequent enough in the labeled treatment).
#'   A hard filter would remove it from both group one and two while a soft filter would remove it only from the group two comparisons.
#'
#'   Specifying \code{rm_low_freq=TRUE} for one or several options means that a taxon will be kept only in those replicates where
#'   it meets or exceeds the fraction threshold. It has been observed that in replicates where a taxon is infrequent, those weighted
#'   average density values are suspect, even if it is frequent in other replicates. We \emph{strongly} recommend that \code{rm_low_freq}
#'   be \code{TRUE}, as this should reduce the occurrence of negative shifts in density.
#'
#' @return \code{create_filters} produces a data frane of replicate frequencies and within-replicate-fraction frequencies to be
#'   investigated in downstream analyses. In addition the logical columns \code{hard} and \code{soft} indicate which frequencies to
#'   be utilized for hard or soft cut-offs.
#'
#' @examples
#' # Only filter on samples
#'  create_filters(2:3)
#'
#' # Filter on samples and fractions
#' create_filters(2:3, 10:12)
#'
#' # remove certain combinations
#' create_filters(2:3, 10:13, c('3:13', '2:10'))
#'
#' @export

create_filters <- function(replicate=0, fraction=0, rm_combn=character(), hard=0, soft=0, rm_low_freq=TRUE) {
  if(length(hard) > 1 || length(soft) > 1) warning('Only one value should be indicated for hard or soft cut-offs. Using first value(s)')
  filters <- expand.grid(replicate, fraction)
  names(filters) <- c('rep_freq', 'frac_freq')
  filters <- filters[order(filters$rep_freq),]
  if(!missing(rm_combn)) {
    rm_combn <- strsplit(rm_combn, ':')
    rm_combn <- do.call(rbind, rm_combn)
    storage.mode(rm_combn) <- 'integer'
    for(i in 1:nrow(rm_combn)) {
      filters <- filters[filters$rep_freq!=rm_combn[i,1] |
                           filters$frac_freq!=rm_combn[i,2],]
    }
  }
  filters$hard <- filters$soft <- logical(nrow(filters))
  filters$hard[hard[1]] <- TRUE
  filters$soft[soft[1]] <- TRUE
  filters$rm_low_freq <- rm_low_freq
  rownames(filters) <- NULL
  return(filters)
}

#' Creates filtering criteria
#'
#' Returns binary vector specifying whether taxa from a phylosip object satisfy minimum frequency specifications
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#'
#' @details \code{impose_filter} is primarily utilized within other functions and imposes a hard filter on the data.
#'
#' @return Returns a
#'
#' @seealso \code{\link{create_filters}}
#'
#' @examples
#' # Filter a tax table

impose_filter <- function(data, replicate=0, fraction=0) {
  # extract feature table and convert to matrix with taxa as columns
  ft <- as(data@otu_table, 'matrix')
  if(phyloseq::taxa_are_rows(data)) ft <- t(ft)
  # make presence-absence
  ft <- ceiling(ft / max(ft))
  storage.mode(ft) <- 'integer'
  # split by replicates
  ft <- split_data(data, ft, data@qsip@rep_id)
  # calculate within replicate (i.e., fraction frequency)
  ft <- lapply(ft, colSums, na.rm=T)
  # combine, apply filter
  #ft <- lapply(ft, function(x) ifelse(x >= fraction, 1, 0))
  ft <- do.call(rbind, ft)
  ft <- ifelse(ft >= fraction, 1, 0)
  # split by replicate group
  iso_group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
  if(length(data@qsip@timepoint) > 0) {
    time_group <- time_grouping(data, data@qsip@timepoint, data@qsip@rep_id, data@qsip@rep_group)
    time_group <- time_group[match(iso_group$replicate, time_group$replicate),]
    iso_group$interaction <- factor(paste(iso_group$interaction, time_group$interaction, sep=':'))
  }
  # Drop any rows (probably NA) that don't appear in ft rownames
  iso_group <- iso_group[match(rownames(ft), iso_group$replicate),]
  ft <- ft[!is.na(iso_group$iso),]
  iso_group <- iso_group[!is.na(iso_group$iso),]
  ft <- split_data(data, ft, iso_group$interaction, grouping_w_phylosip=F)
  # combine, apply filter
  ft <- lapply(ft, colSums, na.rm=T)
  #ft <- lapply(ft, function(x) ifelse(x >= replicate, 1, 0))
  ft <- do.call(rbind, ft)
  ft <- ifelse(ft >= replicate, 1, 0)
  # add in taxa names and calculate presence or absence after filters
  colnames(ft) <- phyloseq::taxa_names(data)
  output <- colSums(ft)
  output <- ifelse(output > 0, 1L, 0L)
}

#' Explore different filtering levels
#'
#' Generates an S4 \code{filter} class object which stores filtering statistics
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param filters Optional data frame specifying filtering levels to explore.
#'   The default is to utilize specifications from \code{data@@qsip@@filter}.
#'   If \code{NULL}, \code{explore_filters} will generate statistics for all possible frequencies at both the replicate and fraction level.
#'
#' @details Some text here
#'
#' @return Some text here
#'
#' @seealso \code{\link{create_filters}}
#'
#' @examples
#' # Filter a tax table
#'
#' @export

explore_filters <- function(data, filters=data@qsip@filter_levels) {
  # if no filters provided, will use all
  if(is.null(filters)) {
    # number of fractions in data?
    fractions <- table(data@sam_data[[data@qsip@rep_id]])
    # group by isotope treatment, replicate ID, replicate group
    group <- iso_grouping(data, data@qsip@iso_trt, data@qsip@rep_id, data@qsip@rep_group)
    group <- group[match(names(fractions), iso_group$replicate),]
    replicates <- split(fractions, group$interaction) # will drop NA isotope values
    # number of replicates in data?
    relicates <- sapply(replicates, length)
    filters <- create_filters(1:max(replicates), 1:max(fractions))
  }
  # create blank matrix (taxa as rows, number of filter combos as columns)
  output <- matrix(0L, nrow=phyloseq::ntaxa(data),
                   ncol=nrow(filters),
                   dimnames=list(phyloseq::tax_names(data0),
                                 interaction(filters, sep=':')))
  # for every combination, run impose filters, stick output to column i of matrix
  for(i in 1:nrow(filters)) {
    output[,i] <- impose_filter(data,
                                replicate=filters$rep_freq[i],
                                fraction=filters$frac_freq[i])
  }
  # convert to S4 Matrix class
  output <- Matrix::Matrix(output, sparse=T)
  # add other statistics into output
  # NEED TO CREATE NEW S4 CLASS FOR FILTERING STATISTICS
  # NEED TO UPDATE SET_CLASS AND SET_INIT_METHODS SCRIPTS
  # maybe one of the slots could be an extension of the input filtering data frame
  #  with columns for the proportion (or number) of remaining taxa and another for remaining samples
}

#' Filter taxa from phylosip
#'
#' Filters taxa from qSIP portion of phyloseq object
#'
#' @param data \code{Phylosip}-class object to pull feature taxa table from.
#' @param replicate Numeric vector specifying the minimum frequency of occurrence of microbial taxa across replicates.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the replicate level
#' @param fraction Numeric vector specifying the minimum frequency of occurrence of microbial taxa across fractions within a sample.
#'   Keeping the default value of \code{0} will apply no frequency threshold at the fraction level
#' @param filter_phyloseq Logical value indicating whether or not to filter taxa from \code{data@@otu_table} and \code{data@@tax_table}.
#'
#' @details \code{filter_qsip} Has twofold actions. First, it produces a character vector of taxa names that have passed the specified
#'   filtering threshold and which are stored in \code{data@@qsip@@filter}. Additionally, any qsip-related data that have been previously
#'   generated and stored in \code{data@@qsip@@.Data} will be filtered to only retain those taxa.
#'
#'   It is not recommended to filter taxa from the \code{phyloseq}-inhereted portions of the data (\emph{e.g.}, from \code{data@@otu_table},
#'   and \code{data@@tax_table}), simply because it will be difficult to call this information again if needed.
#'
#' @return Returns a character vector in \code{data@@qsip@@filter} of taxa names that have passed the specified filtering threshold.
#'
#' @seealso \code{\link{create_filters}}, \code{\link{explore_filters}}
#'
#' @examples
#' # Filter a phyloseq object
#'
#' @export

filter_qsip <- function(data, replicate=0, fraction=0, filter_phyloseq=FALSE) {
  if(is(data)[1]!='phylosip') stop('Must provide phylosip object')
  # if user supplies no frequencies, apply from @qsip@filter_levels
  filter_levels <- data@qsip@filter_levels
  if(missing(replicate)) {
    if(any(filter_levels$hard)) {
      replicate <- filter_levels$rep_freq[which(filter_levels$hard==TRUE)[1]]
    }
  }
  if(missing(fraction)) {
    if(any(filter_levels$hard)) {
      fraction <- filter_levels$frac_freq[which(filter_levels$hard==TRUE)[1]]
    }
  }
  tax_filter <- impose_filter(data, replicate=replicate, fraction=fraction)
  tax_filter <- names(tax_filter[tax_filter > 0])
  if(length(data@qsip) > 0) {
    types <- sapply(data@qsip@.Data, function(x) class(x)[1])
    # remove taxa from each qSIP-related output object in the data@qsip@.Data list
    for(i in 1:length(data@qsip)) {
      x <- data@qsip@.Data[[1]]
      if(types[i]=='dgCMatrix') {
        if(phyloseq::taxa_are_rows(data)) {
          x <- x[rownames(x) %in% tax_filter,]
        }
        else {
          x <- x[,colnames(x) %in% tax_filter]
        }
      } else if(types=='numeric') {
        x <- x[names(x) %in% tax_filter]
      }
      data@qsip@.Data[[1]] <- x
    }
  }
  if(filter_phyloseq) {
    data <- phyloseq::filter_taxa(data, function(x) phyloseq::taxa_names(x@otu_table) %in% tax_filter, TRUE)
  }
  data@qsip@filter <- tax_filter
  return(data)
}
