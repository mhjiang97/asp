library(trackViewer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(glue)

gtf <- "~/doc/reference/gtf/gencode.v32.annotation.gtf"
txdb <- makeTxDbFromGFF(gtf)

sig_rmats <- read.delim("~/projects/dux4/analysis/mgy/rnaseq/as/sig_rmats.txt")
ri_rmats <- sig_rmats[sig_rmats$as_type == "RI", ]

for (i in 1:nrow(ri_rmats)) {
  chr <- ri_rmats[i, 4]
  strand <- ri_rmats[i, 5]
  coords <- as.numeric(ri_rmats[i, c(8:11, 36, 37)])
  left <- min(coords)
  right <- max(coords)
  gene <- ri_rmats[i, 3]

  dux4_1 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_1/DUX4_1.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  dux4_2 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_2/DUX4_2.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  dux4_3 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_3/DUX4_3.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_1 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_1/VEC_1.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_2 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_2/VEC_2.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_3 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_3/VEC_3.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )

  trs <- geneModelFromTxdb(
    txdb, # TxDb.Hsapiens.UCSC.hg38.knownGene
    org.Hs.eg.db,
    gr = GRanges(chr, IRanges(left - 1000, right + 1000), strand = strand)
  )
  viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.07, .02, .02, .02))
  pdf(
    glue("~/projects/dux4/analysis/mgy/rnaseq/figs/track/rmats/{i}_{gene}.pdf"),
    width = 9.5, height = 8.5
  )
  vp <- viewTracks(
    trackList(
      dux4_1, dux4_2, dux4_3, vec_1, vec_2, vec_3,
      trs
    ),
    gr = GRanges(chr, IRanges(left - 1000, right + 1000), strand = strand),
    viewerStyle = viewerStyle,
    autoOptimizeStyle = TRUE
  )
  addGuideLine(coords, vp = vp)
  dev.off()
}


sig_whippet <- read.delim("~/projects/dux4/analysis/mgy/rnaseq/as/sig_whippet.txt")
ri_whippet <- sig_whippet[sig_whippet$Type == "RI", ]
tx2gene <- read.delim("~/doc/tx2gene/tx2gene_v32.txt", header = F)
ri_whippet <- dplyr::left_join(ri_whippet, tx2gene[match(unique(tx2gene[, 2]), tx2gene[, 2]), 2:3], by = c("Gene" = "V2"))

for (i in 1:nrow(ri_whippet)) {
  chr <- strsplit(ri_whippet[i, 3], ":")[[1]][1]
  strand <- ri_whippet[i, 4]
  coords <- as.numeric(strsplit(strsplit(ri_whippet[i, 3], ":")[[1]][2], "-")[[1]])
  left <- min(coords)
  right <- max(coords)
  gene <- ri_whippet[i, 12]

  dux4_1 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_1/DUX4_1.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  dux4_2 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_2/DUX4_2.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  dux4_3 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_3/DUX4_3.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_1 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_1/VEC_1.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_2 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_2/VEC_2.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_3 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_3/VEC_3.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )

  trs <- geneModelFromTxdb(
    txdb,
    org.Hs.eg.db,
    gr = GRanges(chr, IRanges(left - 1000, right + 1000), strand = strand)
  )
  viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.07, .02, .02, .02))
  pdf(
    glue("~/projects/dux4/analysis/mgy/rnaseq/figs/track/whippet/{i}_{gene}.pdf"),
    width = 10, height = 7
  )
  vp <- viewTracks(
    trackList(
      dux4_1, dux4_2, dux4_3, vec_1, vec_2, vec_3,
      trs
    ),
    gr = GRanges(chr, IRanges(left - 1000, right + 1000), strand = strand),
    viewerStyle = viewerStyle,
    autoOptimizeStyle = TRUE
  )
  addGuideLine(coords, vp = vp)
  dev.off()
}


sig_majiq <- read.delim("~/projects/dux4/analysis/mgy/rnaseq/as/sig_majiq.txt")

for (i in 2:nrow(sig_majiq)) {
  chr <- sig_majiq$seqid[i]
  strand <- sig_majiq$strand[i]
  coords1 <- strsplit(sig_majiq$junctions_coords[i], "[;-]")[[1]]
  coords2 <- strsplit(sig_majiq$exons_coords[i], "[;-]")[[1]]
  coords3 <- strsplit(sig_majiq$ir_coords[i], "[;-]")[[1]]
  coords <- unique(as.numeric(c(coords1, coords2, coords3)))
  if (any(is.na(coords))) coords <- coords[!is.na(coords)]
  left <- min(coords)
  right <- max(coords)
  gene <- sig_majiq$gene_name[i]

  dux4_1 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_1/DUX4_1.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  dux4_2 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_2/DUX4_2.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  dux4_3 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_3/DUX4_3.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_1 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_1/VEC_1.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_2 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_2/VEC_2.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )
  vec_3 <- importBam(
    "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_3/VEC_3.sortedByCoord.bam",
    ranges = GRanges(chr, IRanges(left - 1000, right + 1000))
  )

  trs <- geneModelFromTxdb(
    txdb,
    org.Hs.eg.db,
    gr = GRanges(chr, IRanges(left - 1000, right + 1000), strand = strand)
  )
  viewerStyle <- trackViewerStyle()
  setTrackViewerStyleParam(viewerStyle, "xaxis", T)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.07, .02, .02, .02))
  pdf(
    glue("~/projects/dux4/analysis/mgy/rnaseq/figs/track/majiq/{i}_{gene}.pdf"),
    width = 10, height = 7
  )
  vp <- viewTracks(
    trackList(
      dux4_1, dux4_2, dux4_3, vec_1, vec_2, vec_3,
      trs
    ),
    gr = GRanges(chr, IRanges(left - 1000, right + 1000), strand = strand),
    viewerStyle = viewerStyle,
    autoOptimizeStyle = TRUE
  )
  addGuideLine(coords, vp = vp)
  dev.off()
}






