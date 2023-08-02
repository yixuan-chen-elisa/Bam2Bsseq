#' Construct BSseq objects from BED files
#'
#' @details These functions read in nanopore data, modified BED files, to Bsseq objects
#'
#' @param M matrix, methylation data
#' @param Cov matrix, coverage data
#' @param filter matrix, ambiguous methylation modification data
#' @param coef matrix, smoothing estimates
#' @param se.coef matrix, smoothing standard errors
#' @param trans function, smoothing transformation
#' @param parameters list, smoothing parameters
#' @param pData data frame, phenotypic data with samples as rows and variables as columns
#' @param gr GRanges, genomic positions
#' @param pos vector, locations
#' @param chr vector, chromosomes
#' @param sampleNames vector, sample names
#' @param rmZeroCov logical, should genomic locations with zero coverage in all samples be removed
#' @param files vector, BED files.
#'
#' @return BSseq object
#'
#' @import bsseq
#' @export

BSseq2 <- function(M = NULL,Cov = NULL,filter = NULL,coef = NULL,se.coef = NULL,
                   trans = NULL, parameters = NULL, pData = NULL, gr = NULL,
                   pos = NULL, chr = NULL, sampleNames = NULL,
                   rmZeroCov = FALSE) {

  # Argument checks ----------------------------------------------------------

  # Process assays.
  # NOTE: Nothing to do for 'coef', and 'se.coef'.
  if (is.null(M) || is.null(Cov)) {
    stop("Need 'M' and 'Cov'.")
  }
  # Process 'trans' and 'parameters'.
  if (is.null(trans)) {
    trans <- function() NULL
    environment(trans) <- emptyenv()
  }
  if (is.null(parameters)) {
    parameters <- list()
  }
  # Process 'sampleNames' and 'pData'.
  if (is.null(sampleNames)) {
    if (is.null(pData)) {
      # BSseq object will have no colnames.
      pData <- make_zero_col_DFrame(ncol(M))
    } else {
      # BSseq object will have 'sampleNames' as colnames.
      pData <- DataFrame(row.names = sampleNames)
    }
  } else {
    if (is.null(pData)) {
      # BSseq object will have 'sampleNames' as colnames.
      pData <- DataFrame(row.names = sampleNames)
    } else {
      if (is.null(rownames(pData))) {
        rownames(pData) <- sampleNames
      } else {
        stopifnot(identical(rownames(pData), sampleNames))
      }
    }
  }
  # Process 'gr', 'pos', and 'chr'.
  if (is.null(gr)) {
    if (is.null(pos) || is.null(chr)) {
      stop("Need 'pos' and 'chr' if 'gr' not supplied.")
    }
    gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1L))
  }
  if (!is(gr, "GRanges")) {
    stop("'gr' needs to be a GRanges.")
  }
  # Process 'rmZeroCov'.
  stopifnot(isTRUEorFALSE(rmZeroCov))

  # Collapse duplicate loci --------------------------------------------------

  is_duplicated <- duplicated(gr)
  if (any(is_duplicated)) {
    warning("Detected duplicate loci. Collapsing counts in 'M' and 'Cov' ",
            "at these positions.")
    if (!is.null(coef) || !is.null(se.coef)) {
      stop("Cannot collapse when 'coef' or 'se.coef' are non-NULL.")
    }
    loci <- gr[!is_duplicated]
    ol <- findOverlaps(gr, loci, type = "equal")
    M <- rowsum(x = M, group = subjectHits(ol), reorder = FALSE)
    rownames(M) <- NULL
    Cov <- rowsum(x = Cov, group = subjectHits(ol), reorder = FALSE)
    rownames(Cov) <- NULL
    if (!is.null(filter)){
      filter <- rowsum(x = filter, group = subjectHits(ol), reorder = FALSE)
      rownames(filter) <- NULL}
  } else {
    loci <- gr
  }

  # Optionally, remove positions with zero coverage --------------------------

  if (rmZeroCov) {
    loci_with_zero_cov <- rowAlls(Cov, value = 0)
    if (any(loci_with_zero_cov)) {
      loci_with_nonzero_cov <- !loci_with_zero_cov
      loci <- loci[loci_with_nonzero_cov]
      M <- M[loci_with_nonzero_cov, , drop = FALSE]
      Cov <- Cov[loci_with_nonzero_cov, , drop = FALSE]
      coef <- coef[loci_with_nonzero_cov, , drop = FALSE]
      se.coef <- se.coef[loci_with_nonzero_cov, , drop = FALSE]
      filter <- filter[loci_with_nonzero_cov, , drop = FALSE]
    }
  }

  # Construct BSseq object ---------------------------------------------------

  assays <- SimpleList(M = M, Cov = Cov, coef = coef, se.coef = se.coef,
                       filter = filter)
  assays <- assays[!S4Vectors:::sapply_isNULL(assays)]
  se <- SummarizedExperiment(
    assays = assays,
    rowRanges = loci,
    colData = pData)
  bsseq:::.BSseq(se, trans = trans, parameters = parameters)
}


create_bsseq_filter = function(files, pData = NULL){
  gr_list = list()
  sampleNames = sub("\\.bed$","",basename(files))
  if (!is.null(pData)){
    rownames(pData) <- sampleNames
  }

  for (i in seq_along(files)){
    data = read.table(files[i],header = FALSE, sep="\t",
                      stringsAsFactors=FALSE, quote="")
    gr = GRanges(seqnames = data$V1,
                 ranges = IRanges(start = data$V2+1,end = data$V3))
    mcols(gr)$m = data$V13
    mcols(gr)$cov = data$V12 + data$V13
    mcols(gr)$filter = data$V14
    names(gr) = sampleNames[i]
    gr_list[[sampleNames[i]]] = gr
  }

  overlap_gr = Reduce(subsetByOverlaps, gr_list)

  m_cov_list = lapply(gr_list, function(gr){
    overlap_data = gr[gr %over% overlap_gr]
    data.frame(m = overlap_data$m, cov = overlap_data$cov,
               filter = overlap_data$filter)})

  m = do.call(cbind,lapply(m_cov_list,`[[`, "m"))
  cov = do.call(cbind, lapply(m_cov_list, `[[`, "cov"))
  filter = do.call(cbind, lapply(m_cov_list, `[[`, "filter"))

  bsseq_obj = BSseq2(M = as.matrix(m), Cov = as.matrix(cov),
                     filter = as.matrix(filter),
                     coef = NULL,se.coef = NULL,
                     pos = start(overlap_gr),trans = NULL,
                     parameters = NULL, pData = pData, gr = NULL,
                     chr = as.vector(seqnames(overlap_gr)),
                     sampleNames = sampleNames,rmZeroCov = FALSE)

  return(bsseq_obj)
}

