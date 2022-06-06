#' @importFrom tidyr separate
#' @importFrom dplyr select mutate
#' @importFrom stringr str_split
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
ranges_lines <- function(dat, tool, extend_left = 1000, extend_right = 1000) {
  chrs <- genes <- strands <- ids <- vector("character", length = nrow(dat))
  lefts <- rights <- vector("numeric", length = nrow(dat))
  coordses <- vector("list", length = nrow(dat))

  dat <- switch(
    tool,
    rmats = dat,
    cash = dat |>
      tidyr::separate(Location, c("chr", "coord1", "coord2"), sep = "[:-]"),
    leafcutter = dat |>
      tidyr::separate(cluster, c("chr", "clus2"), sep = "\\:") |>
      dplyr::select(!clus2),
    spladder = dat,
    majiq = dat,
    suppa = dat |> dplyr::mutate(chr = str_extract(coord, "chr[0-9XY]*")),
    psichomics = dat |>
      tidyr::separate(event, c("event", "coord"), sep = "\\_[+-]\\_")
  )

  colnames_coord <- switch(
    tool,
    rmats = c(
      "exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES",
      "downstreamEE", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base",
      "2ndExonEnd", "longExonStart_0base", "longExonEnd", "shortES", "shortEE",
      "flankingES", "flankingEE", "riExonStart_0base", "riExonEnd"
    ),
    cash = c("coord1", "coord2"),
    leafcutter = c("myintron"),
    spladder = c(
      "exon_pre_start", "exon_pre_end", "exon_start", "exon_end",
      "exon_aft_start", "exon_aft_end", "exon1_start", "exon1_end",
      "exon2_start", "exon2_end", "exon_const_start", "exon_const_end",
      "exon_alt1_start", "exon_alt1_end", "exon_alt2_start", "exon_alt2_end",
      "intron_start", "intron_end", "exon_starts", "exon_ends"
    ),
    majiq = c("exons_coords", "junctions_coords", "ir_coords"),
    suppa = c("coord"),
    psichomics = c("coord")
  )
  split <- switch(
    tool,
    rmats = F, cash = F, leafcutter = T, spladder = F, majiq = T, suppa = T,
    psichomics = T
  )
  split_by <- switch(
    tool,
    rmats = NULL, cash = NULL, leafcutter = "[\\:\\,]", spladder = NULL,
    majiq = "[;-]", suppa = "[:-]", psichomics = "\\_"
  )

  for (i in 1:nrow(dat)) {
    chrs[i] <- dat$chr[i]
    genes[i] <- dat$gene_symbol[i]
    # strands[i] <- dat$strand[i]
    ids[i] <- dat$myID[i]
    if (split) {
      tmp <- unlist(stringr::str_split(dat[i, colnames_coord], split_by))
    } else {
      tmp <- dat[i, colnames_coord]
    }
    coordses[[i]] <- tmp |>
      as.numeric() |>
      suppressWarnings() |>
      na.omit()
    lefts[i] <- min(coordses[[i]])
    rights[i] <- max(coordses[[i]])
  }
  ranges <- GenomicRanges::GRanges(
    chrs, IRanges::IRanges(lefts - extend_left, rights + extend_right)
  ) # , strand = strands

  list(ranges = ranges, lines = coordses, genes = genes, ids = ids)
}

#' @import  trackViewer
plot_tracks <- function(
  dir_plot, plot_basics, ranges_lines, txdb, annotation_db, ...
) {
  myCreatedir(dir_plot)
  ranges <- ranges_lines$ranges
  lines <- ranges_lines$lines
  genes <- ranges_lines$genes
  ids <- ranges_lines$ids

  bams_1 <- plot_basics$bams_1
  bams_2 <- plot_basics$bams_2
  # lt_1 <- plot_basics$library_types_1
  # lt_2 <- plot_basics$library_types_2

  list_bams <- vector("list", length = length(c(bams_1, bams_2))) |>
    setNames(names(c(bams_1, bams_2)))
  for (s in names(list_bams)) {
    list_bams[[s]] <- trackViewer::importBam(
      c(bams_1, bams_2)[s], ranges = ranges
    ) # , pairs = c(lt_1, lt_2)[s]
  }

  for (i in 1:length(ranges)) {
    message(glue("ploting {i} out of {length(ranges)}"))
    trs <- trackViewer::geneModelFromTxdb(txdb, db, gr = ranges[i])
    optSty <- trackViewer::optimizeStyle(trackViewer::trackList(list_bams, trs))
    trackList <- optSty$tracks
    viewerStyle <- optSty$style
    trackViewer::setTrackViewerStyleParam(viewerStyle, "xaxis", T)
    trackViewer::setTrackViewerStyleParam(
      viewerStyle, "margin", c(0.15, 0.15, 0.05, 0.05)
    )
    for (b in names(bams_1)) {
      trackViewer::setTrackStyleParam(
        trackList[[b]], "color", c("red", "black")
      )
      trackViewer::setTrackYaxisParam(trackList[[b]], "main", F)
    }
    for (b in names(bams_2)) {
      trackViewer::setTrackStyleParam(
        trackList[[b]], "color", c("blue", "black")
      )
      trackViewer::setTrackYaxisParam(trackList[[b]], "main", F)
    }
    pdf(
      glue::glue("{dir_plot}/{i}_{genes[i]}_{ids[i]}.pdf"),
      width = 10, height = 9
    )
    vp <- trackViewer::viewTracks(
      trackList, gr = ranges[i], viewerStyle = viewerStyle,
      autoOptimizeStyle = T
    )
    trackViewer::addGuideLine(lines[[i]], vp = vp)
    dev.off()
  }
}

getPlotBasics <- function(asp, plot_samples, ...) {
  if (class(plot_samples) == "numeric") {
    bams_1 <- sample(asp@basics$bams_1, plot_samples)
    bams_2 <- sample(asp@basics$bams_2, plot_samples)
  } else if (class(plot_samples) == "character") {
    bams_1 <- na.omit(asp@basics$bams_1[plot_samples]) # |> as.character()
    bams_2 <- na.omit(asp@basics$bams_2[plot_samples]) # |> as.character()
    if (all(is.na(bams_1))) bams_1 <- NULL
    if (all(is.na(bams_2))) bams_2 <- NULL
  }
  # samples_1 <- names(bams_1)
  # samples_2 <- names(bams_2)
  # if (all(is.na(samples_1))) samples_1 <- NULL
  # if (all(is.na(samples_2))) samples_2 <- NULL
  # library_types_1 <- asp@sampletable$library_types[match(samples_1, asp@sampletable$samples)]
  # library_types_2 <- asp@sampletable$library_types[match(samples_2, asp@sampletable$samples)]
  # if (all(is.na(library_types_1))) library_types_1 <- NULL
  # if (all(is.na(library_types_2))) library_types_2 <- NULL
  # lt_1 <- vector("logical", length = length(library_types_1)) |>
  #   setNames(samples_1)
  # lt_2 <- vector("logical", length = length(library_types_2)) |>
  #   setNames(samples_2)
  # for (i in 1:length(lt_1)) lt_1[i] <- ifelse(library_types_1[i] == "paired-end", T, F)
  # for (i in 1:length(lt_2)) lt_2[i] <- ifelse(library_types_2[i] == "paired-end", T, F)
  # if (all(is.na(lt_1))) lt_1 <- NULL
  # if (all(is.na(lt_2))) lt_2 <- NULL

  list(bams_1 = bams_1, bams_2 = bams_2) # , library_types_1 = lt_1, library_types_2 = lt_2
}




