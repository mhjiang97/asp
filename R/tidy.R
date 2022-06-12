#' @importFrom dplyr arrange desc filter
querySig <- function(dat, ps) {
  p_adj <- ps$p_adj
  p_value <- ps$p_value
  delta_psi <- ps$delta_psi
  statistic_col <- ps$statistic_col
  statistic_threshold <- ps$statistic_threshold
  filter_direction <- ps$filter_direction
  effect_size_col <- ps$effect_size_col
  effect_size_threshod <- ps$effect_size_threshod

  if (!is.null(statistic_col)) {
    dat[[statistic_col]] <- as.numeric(dat[[statistic_col]])
    if (filter_direction[1] == ">") dat_sig <- dat[dat[[statistic_col]] >= statistic_threshold, ] |> dplyr::arrange(dplyr::desc(.data[[statistic_col]]))
    if (filter_direction[1] == "<") dat_sig <- dat[dat[[statistic_col]] <= statistic_threshold, ] |> dplyr::arrange(.data[[statistic_col]])
  } else {
    if (is.null(p_value)) {
      dat$padj <- as.numeric(dat$padj)
      dat_sig <- dplyr::filter(dat, padj <= p_adj) |> dplyr::arrange(padj)
    } else {
      dat$p <- as.numeric(dat$p)
      dat_sig <- dplyr::filter(dat, p <= p_value) |> dplyr::arrange(p)
    }
  }

  ## always filter according effect size after p value ##
  if (!is.null(effect_size_col)) {
    dat_sig[[effect_size_col]] <- as.numeric(dat_sig[[effect_size_col]])
    dat_sig <- dat_sig[abs(dat_sig[[effect_size_col]]) >= effect_size_threshod, ] |> dplyr::arrange(dplyr::desc(.data[[effect_size_col]]))
  } else {
    dat_sig$dpsi <- as.numeric(dat_sig$dpsi)
    dat_sig <- dplyr::filter(dat_sig, abs(dpsi) >= delta_psi) |> dplyr::arrange(dplyr::desc(dpsi))
  }

  dat_sig
}

#' @importFrom vroom vroom
#' @importFrom utils read.table
myVroom <- function(file, na_append = NULL, delim = "\t", comment = "", col_names = T) {
  na_origin <- c("", "NA")
  myna <- c(na_origin, na_append)

  n_col <- utils::read.table(
    file, header = T, sep = "\t", fill = T, na.strings = c("NA", na_append), nrows = 10
  ) |>
    ncol()

  dat <- vroom::vroom(
    file, delim = delim, col_types = paste0(rep("c", n_col), collapse = ""),
    na = myna, comment = comment, col_names = col_names
  )

  dat
}

#' @importFrom glue glue
#' @importFrom stats setNames
#' @importFrom futile.logger flog.warn flog.error
findReadFiles <- function(types_all, types_event, dir, patterns) {
  et <- intersect(types_event, types_all)
  if ("all" %in% types_event) et <- types_all

  files_read <- vector("character", length(et)) |> stats::setNames(et)
  for (e in names(files_read)) {
    tmp <- list.files(path = dir, pattern = patterns[e], full.names = T)
    if (length(tmp) == 0) {
      futile.logger::flog.warn(glue::glue("The file of type '{e}' can't be found"))
      tmp <- ""
    }
    files_read[e] <- tmp
  }

  if (any(files_read == "")) files_read <- files_read[-which(files_read == "")]
  if (length(files_read) == 0) return(futile.logger::flog.error("Can't find any files to read in"))

  files_read
}

#' @importFrom glue glue
#' @importFrom dplyr mutate select starts_with rename bind_rows
#' @importFrom tibble remove_rownames column_to_rownames rownames_to_column as_tibble
read_rmats <- function(asp) {
  dir <- glue::glue("{asp@dir_out}/rMATS/")
  read_type <- "JC"
  event_type <- "all"
  types_all <- c("SE", "MXE", "A3SS", "A5SS", "RI")
  patterns <- glue::glue("{types_all}.*\\.JC\\.txt.*") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )

  list_rmats <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_rmats)) {
    f <- files_read[e]
    list_rmats[[e]] <- f |>
      myVroom(na_append = c("nan", "NaN")) |>
      dplyr::mutate(as_type = e) |>
      tibble::remove_rownames() |>
      tibble::column_to_rownames("ID...1") |>
      dplyr::select(!dplyr::starts_with("ID")) |>
      tibble::rownames_to_column("ID") |>
      dplyr::rename(
        padj = FDR, dpsi = IncLevelDifference,
        gene_id = GeneID, gene_symbol = geneSymbol,
        p = PValue, psis_1 = IncLevel1, psis_2 = IncLevel2
      ) |>
      tibble::as_tibble()
  }

  rmats <- dplyr::bind_rows(list_rmats) |>
    dplyr::mutate(myID = paste0(as_type, ": ", ID))
  rmats
}

#' @importFrom glue glue
#' @importFrom dplyr rename mutate case_when all_of rowwise ungroup left_join
#' @importFrom tidyr separate
read_cash <- function(asp) {
  MergePval <- "G"
  if (!is.null(asp@parameters$MergePval)) MergePval <- asp@parameters$MergePval
  file <- glue::glue(
    "{asp@dir_out}/CASH/cash_{MergePval}.{asp@basics$condition_1}vs{asp@basics$condition_2}.alldiff.txt"
  )
  condition_1 <- asp@basics$condition_1
  condition_2 <- asp@basics$condition_2

  cash <- myVroom(file, na_append = c("nan", "NaN")) |>
    dplyr::rename(
      gene_symbol = AccID, dpsi = delta_PSI, p = `P-Value`, padj = FDR
    ) |>
    dplyr::mutate(
      as_type = dplyr::case_when(
        SplicingType == "IR" ~ "RI",
        SplicingType == "Cassette" ~ "SE",
        SplicingType == "Cassette_multi" ~ "SME",

        SplicingType == "MXE" ~ "MXE",
        SplicingType == "A3SS" ~ "A3SS",
        SplicingType == "A5SS" ~ "A5SS",
        SplicingType == "AltStart" ~ "AFE",
        SplicingType == "AltEnd" ~ "ALE"
      )
    ) |>
    tidyr::separate(
      dplyr::all_of(glue::glue("{condition_1}_Junc_Inclusive::Exclusive")),
      c("inclusive_1", "exlclusive_1"),
      "::"
    ) |>
    tidyr::separate(
      dplyr::all_of(glue::glue("{condition_2}_Junc_Inclusive::Exclusive")),
      c("inclusive_2", "exlclusive_2"),
      "::"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      average_psi_1 = as.numeric(inclusive_1) / (as.numeric(inclusive_1) + as.numeric(exlclusive_1)),
      average_psi_2 = as.numeric(inclusive_2) / (as.numeric(inclusive_2) + as.numeric(exlclusive_2))
    ) |>
    dplyr::ungroup() |>
    dplyr::left_join(asp@tx2gene_convert, by = "gene_symbol") |>
    dplyr::mutate(myID = paste0(gene_symbol, ": ", Location))

  cash
}

#' @importFrom glue glue
#' @importFrom dplyr select left_join mutate case_when rename all_of
#' @importFrom tidyr separate_rows
read_majiq <- function(asp) {
  threshold <- 0.2
  if (!is.null(asp@parameters$threshold)) threshold <- asp@parameters$threshold
  show_all <- F
  if (!is.null(asp@parameters$show_all)) show_all <- asp@parameters$show_all
  condition_1 <- asp@basics$condition_1
  condition_2 <- asp@basics$condition_2

  delta_psi_file <- glue::glue(
    "{asp@dir_out}/MAJIQ/deltapsi/{condition_1}_{condition_2}.deltapsi.tsv"
  )
  if (show_all) {
    voila_file <- glue::glue(
      "{asp@dir_out}/MAJIQ/voila/{condition_1}_{condition_2}_{threshold}_showall.tsv"
    )
  } else {
    voila_file <- glue::glue(
      "{asp@dir_out}/MAJIQ/voila/{condition_1}_{condition_2}_{threshold}.tsv"
    )
  }

  delta_psi <- myVroom(delta_psi_file) |>
    dplyr::select(`LSV ID`, A5SS, A3SS, ES, `Num. Junctions`, `Num. Exons`)
  voila <- myVroom(voila_file, comment = "#") |>
    dplyr::left_join(delta_psi, by = c("lsv_id" = "LSV ID")) |>
    dplyr::mutate(
      as_type = dplyr::case_when(
        as.numeric(`Num. Junctions`) == 1 & ES == "False" & A5SS == "False" & A3SS == "False" & !is.na(`ir_coords`) ~ "RI",
        as.numeric(`Num. Junctions`) == 2 & as.numeric(`Num. Exons`) == 3 ~ "SE",
        as.numeric(`Num. Junctions`) == 4 & as.numeric(`Num. Exons`) == 4 ~ "MXE",
        as.numeric(`Num. Junctions`) %in% c(2, 3) & A3SS == "True" ~ "A3SS",
        as.numeric(`Num. Junctions`) %in% c(2, 3) & A5SS == "True" ~ "A5SS",
        T ~ "OTHER"
      )
    ) |>
    dplyr::rename(
      gene_symbol = gene_name,
      chr = seqid,
      dpsi = mean_dpsi_per_lsv_junction,
      p = probability_changing,
      average_psi_1 = dplyr::all_of(glue::glue("{condition_1}_mean_psi")),
      average_psi_2 = dplyr::all_of(glue::glue("{condition_2}_mean_psi"))
    ) |>
    dplyr::mutate(myID = lsv_id)

  if (T) {
    voila_flatten <- voila |>
      tidyr::separate_rows(
        dpsi, p, probability_non_changing, average_psi_1, average_psi_2,
        sep = ";"
      )
  }

  list(voila = voila, voila_flatten = voila_flatten)
}

#' @importFrom glue glue
#' @importFrom dplyr rowwise mutate ungroup group_by distinct select left_join rename
#' @importFrom stringr str_split
#' @importFrom tidyr separate_rows
read_leafcutter <- function(asp) {
  file_significance <- glue::glue("{asp@dir_out}/LeafCutter/ds/leafcutter_cluster_significance.txt")
  file_effect_sizes <- glue::glue("{asp@dir_out}/LeafCutter/ds/leafcutter_effect_sizes.txt")
  condition_1 <- asp@basics$condition_1
  condition_2 <- asp@basics$condition_2

  leafcutter_size <- myVroom(file_effect_sizes) |>
    dplyr::rowwise() |>
    dplyr::mutate(clus = stringr::str_split(intron, ":")[[1]][4]) |>
    dplyr::ungroup() |>
    dplyr::group_by(clus) |>
    dplyr::mutate(
      average_psi_1 = paste(.data[[condition_1]], collapse = ","),
      average_psi_2 = paste(.data[[condition_2]], collapse = ","),
      dpsi = paste(deltapsi, collapse = ","),
      myintron = paste(intron, collapse = ",")
    ) |>
    dplyr::ungroup() |>
    dplyr::distinct(clus, .keep_all = T) |>
    dplyr::select(!c(logef, intron, s1, s2, deltapsi))

  leafcutter_result <- myVroom(file_significance, na_append = c("nan", "NaN")) |>
    dplyr::rowwise() |>
    dplyr::mutate(clus = stringr::str_split(cluster, ":")[[1]][2]) |>
    dplyr::ungroup() |>
    tidyr::separate_rows(genes, sep = ",")

  leafcutter <- dplyr::left_join(leafcutter_result, leafcutter_size, by = "clus") |>
    dplyr::rename(padj = p.adjust, gene_symbol = genes) |>
    dplyr::mutate(myID = cluster, as_type = "OTHER")

  leafcutter_flatten <- leafcutter  |>
    tidyr::separate_rows(average_psi_1, average_psi_2, dpsi, myintron, sep = ",")

  list(leafcutter = leafcutter, leafcutter_flatten = leafcutter_flatten)
}

#' @importFrom glue glue
#' @importFrom stringr str_replace_all
#' @importFrom dplyr rowwise mutate ungroup rename case_when left_join bind_rows
read_spladder <- function(asp) {
  dir_test <- glue::glue("{asp@dir_out}/SplAdder/testing_{asp@basics$condition_1}_vs_{asp@basics$condition_2}")
  dir_event <- glue::glue("{asp@dir_out}/SplAdder/")
  event_type <- "all"
  confidence_level <- 3
  if (!is.null(asp@parameters$confidence_level)) confidence_level <- asp@parameters$confidence_level
  types_all <- c("exon_skip", "mutex_exons", "alt_3prime", "alt_5prime", "intron_retention", "mult_exon_skip")

  patterns_event <- glue::glue(
    "merge_graphs_{types_all}_C{confidence_level}.confirmed.txt"
  ) |>
    setNames(types_all)
  files_read_event <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir_event, patterns = patterns_event
  )

  list_spladder_event <- vector("list", length(files_read_event)) |>
    setNames(names(files_read_event))
  for (e in names(list_spladder_event)) {
    f <- files_read_event[e]
    list_spladder_event[[e]] <- myVroom(f, na_append = c("nan", "NaN", "inf", "-inf")) |>
      dplyr::rowwise() |>
      dplyr::mutate(event_id = stringr::str_replace_all(event_id, "\\.", "\\_")) |>
      dplyr::ungroup()
  }

  patterns <- glue::glue(".*extended_C{confidence_level}_{types_all}\\.tsv.*") |>
    as.character() |>
    setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir_test, patterns = patterns
  )

  list_spladder <- vector("list", length(files_read)) |>
    setNames(names(files_read))
  for (e in names(list_spladder)) {
    f <- files_read[e]
    list_spladder[[e]] <- myVroom(f, na_append = c("nan", "NaN", "inf", "-inf")) |>
      dplyr::mutate(as_type = e) |>
      dplyr::rename(
        padj = p_val_adj, gene_id = gene, p = p_val
      ) |>
      dplyr::mutate(
        as_type = dplyr::case_when(
          as_type == "exon_skip" ~ "SE",
          as_type == "mutex_exons" ~ "MXE",
          as_type == "alt_3prime" ~ "A3SS",
          as_type == "alt_5prime" ~ "A5SS",
          as_type == "intron_retention" ~ "RI",
          as_type == "mult_exon_skip" ~ "SME"
        )
      )
  }

  for (e in names(list_spladder)) {
    list_spladder[[e]] <- list_spladder[[e]] |>
      dplyr::left_join(
        list_spladder_event[[e]],
        by = c("gene_id" = "gene_name", "event_id" = "event_id")
      )
  }

  spladder <- dplyr::bind_rows(list_spladder) |>
    dplyr::left_join(asp@tx2gene_convert, by = "gene_id") |>
    dplyr::mutate(myID = event_id) |>
    dplyr::rename(chr = contig)

  spladder
}

#' @importFrom glue glue
#' @importFrom dplyr rename select left_join mutate
read_bandits <- function(asp) {
  transcript_result_file <- glue::glue("{asp@dir_out}/BANDITS/results_transcript.txt")
  gene_result_file <- glue::glue("{asp@dir_out}/BANDITS/results_gene.txt")
  bandits <- myVroom(transcript_result_file) |>
    dplyr::rename(
      gene_id = Gene_id, transcript_id = Transcript_id,
      p = p.values, padj = adj.p.values
    )

  if (!is.null(gene_result_file)) {
    bandits_gene <- myVroom(gene_result_file) |>
      dplyr::select(Gene_id, DTU_measure)

    bandits <- bandits |>
      dplyr::left_join(
        bandits_gene, by = c("gene_id" = "Gene_id")
      )
  }
  bandits <- bandits |>
    dplyr::left_join(asp@tx2gene_convert, by = "gene_id") |>
    dplyr::mutate(myID = paste0(gene_id, transcript_id), as_type = "OTHER")

  bandits
}

#' @importFrom glue glue
#' @importFrom dplyr mutate case_when bind_rows left_join
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom tidyr separate drop_na
#' @importFrom utils read.delim
#' @importFrom stats setNames
read_suppa <- function(asp) {
  dir <- glue::glue("{asp@dir_out}/SUPPA/ds/")
  event_type <- "all"
  gene_correction <- T
  types_all <- c("SE", "MX", "A3", "A5", "RI", "AF", "AL")
  colnames_suppa <- c("dpsi", "p")
  if (gene_correction) colnames_suppa <- c("dpsi", "padj")

  patterns <- glue::glue(".*\\.{types_all}\\.dpsi.*") |>
    as.character() |>
    stats::setNames(types_all)
  files_read <- findReadFiles(
    types_all = types_all, types_event = event_type, dir = dir, patterns = patterns
  )

  list_suppa <- vector("list", length(files_read)) |>
    stats::setNames(names(files_read))
  for (e in names(list_suppa)) {
    f <- files_read[e]
    list_suppa[[e]] <- utils::read.delim(f) |>
      stats::setNames(colnames_suppa) |>
      dplyr::mutate(as_type = e) |>
      tibble::rownames_to_column("event") |>
      tidyr::separate(event, c("gene_id", "coord"), sep = ";") |>
      dplyr::mutate(
        as_type = dplyr::case_when(
          as_type == "MX" ~ "MXE",
          as_type == "A3" ~ "A3SS",
          as_type == "A5" ~ "A5SS",

          as_type == "SE" ~ "SE",
          as_type == "RI" ~ "RI",
          as_type == "AF" ~ "AFE",
          as_type == "AL" ~ "ALE"
        )
      ) |>
      tibble::as_tibble()
  }

  suppa <- dplyr::bind_rows(list_suppa) |>
    dplyr::left_join(asp@tx2gene_convert, by = "gene_id") |>
    dplyr::mutate(myID = coord) |>
    tidyr::drop_na(dpsi)

  suppa
}

#' @importFrom glue glue
#' @importFrom dplyr rename mutate case_when left_join select bind_rows
read_psichomics <- function(asp) {
  result_file <- glue::glue("{asp@dir_out}/psichomics/results.txt")
  psichomics <- myVroom(result_file) |>
    dplyr::rename(
      event = `...1`, chr = Chromosome, strand = Strand, gene_symbol = Gene,
      p = `T-test p-value`, padj = `T-test p-value (BH adjusted)`, dpsi = "\u2206 Median",
      average_psi_1 = `Median (s1)`, average_psi_2 = `Median (s2)`
    ) |>
    dplyr::mutate(
      as_type = dplyr::case_when(
        `Event type` == "Alternative 3' splice site (A3SS)" ~ "A3SS",
        `Event type` == "Alternative 5' splice site (A5SS)" ~ "A5SS",
        `Event type` == "Alternative first exon (AFE)" ~ "AFE",
        `Event type` == "Alternative last exon (ALE)" ~ "ALE",
        `Event type` == "Mutually exclusive exon (MXE)" ~ "MXE",
        `Event type` == "Skipped exon (SE)" ~ "SE"
      ),
      chr = paste0("chr", chr)
    ) |>
    dplyr::left_join(asp@tx2gene_convert, by = "gene_symbol")

  psichomics_symbol <- psichomics[!grepl("ENSG", psichomics$gene_symbol), ]
  psichomics_id <- psichomics[grepl("ENSG", psichomics$gene_symbol), ] |>
    dplyr::select(!gene_id) |>
    dplyr::rename(gene_id = gene_symbol) |>
    dplyr::left_join(asp@tx2gene_convert, by = "gene_id")
  psichomics <- dplyr::bind_rows(psichomics_symbol, psichomics_id) |>
    dplyr::mutate(myID = event)

  psichomics
}




