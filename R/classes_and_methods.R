#' Class \code{ASP}
#' @name ASP-class
#' @description  All information stored in ASP object.
#'    You can use \code{creatASP} to create an ASP object.
#'    In this package, most of the functions will use
#'    ASP object as input, and return a modified ASP as well.
#' @slot sampletable data.frame. An experiment design table and
#'     these columns are mandatory:
#'     \itemize{
#'         \item samples
#'         \item conditions
#'         \item files_bam
#'         \item dirs_salmon
#'         \item read_lengths
#'         \item library_types
#'         \item strandedness
#'     }
#' @slot tx2gene data.frame. At least three types of
#'     information for each transcript should be provided:
#'     \itemize{
#'         \item gene_id
#'         \item gene_symbol
#'         \item transcript_id
#'     }
#' @slot tx2gene_convert data.frame. Automatically generated from tx2gene.
#' @slot gtf character. The path to your gtf file.
#' @slot gff character. The path to your gff3 file.
#' @slot fa character. The path to your genome fasta file.
#' @slot fa_transcript character. The path to your transcriptome fasta file.
#' @slot genome_version character. hg38 or hg19.
#' @slot dir_out character. The parent directory for all the next procedures.
#' @slot nproc integer. The number of threads.
#' @slot novel boolean. If detect novel splicing or not.
#' @slot write_log boolean. If write the stdout and stderr to a log file
#'     beginning with an underscore.
#' @slot parallel boolean. If run different tools or multiple samples
#'     in parallel.
#' @slot np integer. How many samples will be processed at the same time.
#' @slot conda_path character. The path to your conda or vritual environment.
#'     Different tools can rely on different paths.
#' @slot basics list. Automatically generated based on the sampletable.
#' @slot cmds list. Storing commands of different tools to run or to check.
#' @slot results list. Storing results outputted by different tools.
#' @slot sig_results list. Results fulfilling significant thresholds.
#' @slot parameters list.
#' @slot intersection list. Results intersections.
#' @slot txdb TxDb. Generated from the gtf file.
#'     \code{\link[GenomicFeatures]{makeTxDbFromGFF}}.
#' @slot annotation_db OrgDb. Bioconductor annotation data package.
#'     for example: \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}.
#' @exportClass ASP
#' @export ASP
#' @return NULL
ASP <- setClass(
  "ASP", slots = c(
    sampletable = "data.frame",
    tx2gene = "data.frame",
    tx2gene_convert = "data.frame",
    gtf = "character",
    gff = "character",
    fa = "character",
    fa_transcript = "character",
    genome_version = "character",
    dir_out = "character",
    nproc = "numeric",
    novel = "logical",
    write_log = "logical",
    parallel = "logical",
    np = "numeric",
    conda_path = "character",

    basics = "list",
    cmds = "list",
    results = "list",
    sig_results = "list",
    parameters = "list",
    intersection = "list",
    txdb = "ANY",
    annotation_db = "ANY"
  )
)

setMethod("initialize", "ASP", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  if (all(is.na(.Object@tx2gene))) .Object@tx2gene <- data.frame(gene_id = NULL, gene_symbol = NULL, transcript_id = NULL)
  if (all(is.na(.Object@dir_out))) .Object@dir_out <- "./"
  if (all(is.na(.Object@nproc))) .Object@nproc <- 5
  if (all(is.na(.Object@novel))) .Object@novel <- F
  if (all(is.na(.Object@write_log))) .Object@write_log <- T
  if (all(is.na(.Object@parallel))) .Object@parallel <- T
  if (all(is.na(.Object@np))) .Object@np <- 5

  if (all(is.na(.Object@gtf))) .Object@gtf <- ""
  if (all(is.na(.Object@gff))) .Object@gff <- ""
  if (all(is.na(.Object@fa))) .Object@fa <- ""
  if (all(is.na(.Object@fa_transcript))) .Object@fa_transcript <- ""
  if (all(is.na(.Object@genome_version))) .Object@genome_version <- ""
  if (all(is.na(.Object@conda_path))) .Object@conda_path <- ""
  .Object@basics <- getBasics(.Object@sampletable)
  # .Object@tx2gene_convert <- .Object@tx2gene |>
  #   dplyr::distinct(gene_id, .keep_all = T) |>
  #   dplyr::select(gene_id, gene_symbol)

  .Object
})

#' @importFrom glue glue
validAspObject <- function(object) {
  mycols_needed <- c(
    "samples", "conditions", "read_lengths", "library_types", "strandedness"
  )
  mycols_dispensable <- c(
    "files_bam_tophat", "files_bam", "files_fq", "dirs_salmon"
  )
  mycols_atleast1 <- c(
    "files_bam_tophat", "files_bam", "files_fq", "dirs_salmon"
  )
  for (mycol in mycols_dispensable) {
    if (!mycol %in% colnames(object@sampletable)) {
      warning(glue::glue("You didn't provide ASP {mycol}, so some alternative splicing algorithm cannot be implemented."))
    }
  }
  if (intersect(c("paired-end", "single-end"), unique(object@sampletable$library_types)) == 0) {
    errorCondition(
      "A wrong sampletable provided. Please check your library types!\n
      Only paired-end and/or single-end are supported!"
    )
  }
  if (nrow(object@tx2gene) == 0 || sum(c("gene_id", "gene_symbol", "transcript_id") %in% colnames(object@tx2gene)) != 3) {
    warning(
      "tx2gene is needed and should contain at least three following columns: gene_id, gene_symbol, transcript_id.
      Automatically generated then."
    )
  }

  if (
    length(intersect(mycols_needed, colnames(object@sampletable))) != length(mycols_needed) ||
    all(is.na(colnames(object@sampletable))) ||
    length(intersect(mycols_atleast1, colnames(object@sampletable))) < 1
  ) {errorCondition("A wrong sampletable provided!")}
}
setValidity("ASP", validAspObject)

#' Create an \code{\linkS4class{ASP}} object
#' @description This function is about how to build an \code{\linkS4class{ASP}}
#'    object. An \code{\linkS4class{ASP}} object is the base for the whole
#'    analysing workflow.
#' @name createASP
#' @param ... paramters pass to \code{\linkS4class{ASP}}.
#' @importFrom rtracklayer import
#' @importFrom tibble as_tibble
#' @importFrom dplyr distinct select
#' @export
#' @return An \code{\linkS4class{ASP}} object.
createASP <- function(...) {
  object <- ASP(...)

  if (nrow(object@tx2gene) == 0 || sum(c("gene_id", "gene_symbol", "transcript_id") %in% colnames(object@tx2gene)) != 3) {
    tmp <- rtracklayer::import(asp@gtf)
    tmp <- tmp[tmp$type == "transcript"]
    chr <- as.character(tmp@seqnames)
    gene_id <- tmp$gene_id
    gene_symbol <- tmp$gene_name
    transcript_id <- tmp$transcript_id

    object@tx2gene <- data.frame(
      chr = chr, gene_id = gene_id, gene_symbol = gene_symbol, transcript_id = transcript_id
    ) |> tibble::as_tibble()
  }

  object@tx2gene_convert <- object@tx2gene |>
    dplyr::distinct(gene_id, .keep_all = T) |>
    dplyr::select(gene_id, gene_symbol)

  object
}

#' Generate commands for AS tools
#' @description This function is about to generate desire commands to run or to
#'     check.
#' @name cmd_asp
#' @param asp An \code{\linkS4class{ASP}} object.
#' @param conda_env The path to environment each tools depend on.
#'     Different tools can rely on different envs.
#' @param tx2gene,dir_out,nproc,gtf,gff,fa,fa_transcript,genome_version
#'     These are same
#'     as the flags in \code{\linkS4class{ASP}}, if you want to change any of
#'     them, then provide a new one.
#' @param parallel,np,basics,sampletable,novel,conda_path,write_log
#'     These are same
#'     as the flags in \code{\linkS4class{ASP}}, if you want to change any of
#'     them, then provide a new one.
#' @param tool Any combination of "rmats", "cash", "leafcutter", "majiq",
#'     "suppa", "bandits", "psichomics", and "spladder".
#' @param ... paramters pass to \code{\link{rmats}}, \code{\link{cash}},
#'     \code{\link{leafcutter}}, \code{\link{majiq}}, \code{\link{suppa}},
#'     \code{\link{bandits}}, \code{\link{psichomics}},
#'     and \code{\link{spladder}}.
#' @importFrom dplyr distinct select
#' @export
#' @return An \code{\linkS4class{ASP}} object with commands of each workflows
#'     saved in cmds slot.
cmd_asp <- function(
  asp, conda_env = NULL, tool, tx2gene = NULL,
  dir_out = NULL, nproc = NULL, gtf = NULL, gff = NULL, fa = NULL,
  fa_transcript = NULL, genome_version = NULL, parallel = NULL, np = NULL,
  basics = NULL, sampletable = NULL, novel = NULL, conda_path = NULL, write_log = NULL, ...
) {
  # asp@parameters_addition <- list(...)

  if (is.null(dir_out)) {
    dir_out <- asp@dir_out
  } else {
    asp@dir_out <- dir_out
  }
  if (is.null(nproc)) nproc <- asp@nproc
  if (is.null(gtf)) gtf <- asp@gtf
  if (is.null(gff)) gff <- asp@gff
  if (is.null(fa)) fa <- asp@fa
  if (is.null(fa_transcript)) fa_transcript <- asp@fa_transcript
  if (is.null(genome_version)) genome_version <- asp@genome_version
  if (is.null(parallel)) parallel <- asp@parallel
  if (is.null(np)) np <- asp@np
  if (is.null(basics)) basics <- asp@basics
  if (is.null(sampletable)) {
    sampletable <- asp@sampletable
  } else {
    asp@sampletable <- sampletable
    asp@basics <- getBasics(asp@sampletable)
  }
  if (is.null(tx2gene)) {
    tx2gene <- asp@tx2gene
  } else {
    asp@tx2gene <- tx2gene
    asp@tx2gene_convert <- asp@tx2gene |>
      dplyr::distinct(gene_id, .keep_all = T) |>
      dplyr::select(gene_id, gene_symbol)
  }
  if (is.null(novel)) novel <- asp@novel
  if (is.null(conda_path)) conda_path <- asp@conda_path
  if (is.null(write_log)) write_log <- asp@write_log

  # for (t in tool) asp@parameters[[t]] <- list(...)
  asp@parameters <- list(...)

  if (!is.null(conda_env)) {
    if (length(conda_env) != length(tool)) {
      if (is.null(names(conda_env))) {
        conda_env <- rep(conda_env[1], length(tool))
      } else {
        ce <- rep(NA, length(tool)) |> setNames(tool)
        ce[names(conda_env)] <- conda_env
        conda_env <- ce
      }
    }
    if (is.null(names(conda_env))) names(conda_env) <- tool
  }

  if (!is.null(conda_path)) {
    if (length(conda_path) != length(tool)) {
      if (is.null(names(conda_path))) {
        conda_path <- rep(conda_path[1], length(tool))
      } else {
        cp <- rep(asp@conda_path, length(tool)) |> setNames(tool)
        cp[names(conda_path)] <- conda_path
        conda_path <- cp
      }
    }
    if (is.null(names(conda_path))) names(conda_path) <- tool
  }
  # for (i in length(conda_path)) {
  #   if (is.na(conda_path[i])) conda_path[i] <- asp@conda_path
  # }

  if ("rmats" %in% tool) {
    asp@cmds[["rmats"]] <- rmats(
      dir_out = dir_out, basics = basics, nproc = nproc, gtf = gtf, novel = novel,
      conda_path = conda_path["rmats"], conda_env = conda_env["rmats"],
      write_log = write_log, ...
    )
  }
  if ("cash" %in% tool) {
    asp@cmds[["cash"]] <- cash(
      dir_out = dir_out, basics = basics, gtf = gtf, novel = novel,
      conda_path = conda_path["cash"], conda_env = conda_env["cash"],
      write_log = write_log, ...
    )
  }
  if ("majiq" %in% tool) {
    asp@cmds[["majiq"]] <- majiq(
      dir_out = dir_out, basics = basics, sampletable = sampletable, nproc = nproc,
      gff = gff, novel = novel, fa = fa,
      conda_path = conda_path["majiq"], conda_env = conda_env["majiq"],
      write_log = write_log, ...
    )
  }
  if ("leafcutter" %in% tool) {
    asp@cmds[["leafcutter"]] <- leafcutter(
      dir_out = dir_out, sampletable = sampletable, nproc = nproc, np = np, gtf = gtf,
      conda_path = conda_path["leafcutter"], conda_env = conda_env["leafcutter"],
      write_log = write_log, parallel = parallel, ...
    )
  }
  if ("spladder" %in% tool) {
    asp@cmds[["spladder"]] <- spladder(
      dir_out = dir_out, basics = basics, gtf = gtf, novel = novel,
      conda_path = conda_path["spladder"], conda_env = conda_env["spladder"],
      write_log = write_log, ...
    )
  }
  if ("bandits" %in% tool) {
    asp@cmds[["bandits"]] <- bandits(
      dir_out = dir_out, sampletable = sampletable, basics = basics, nproc = nproc,
      conda_path = conda_path["bandits"], conda_env = conda_env["bandits"],
      write_log = write_log, tx2gene = tx2gene |> dplyr::select(gene_id, transcript_id),
      ...
    )
  }
  if ("suppa" %in% tool) {
    asp@cmds[["suppa"]] <- suppa(
      dir_out = dir_out, gtf = gtf, np = np, parallel = parallel,
      conda_path = conda_path["suppa"], conda_env = conda_env["suppa"],
      write_log = write_log, ...
    )
  }
  if ("psichomics" %in% tool) {
    asp@cmds[["psichomics"]] <- psichomics(
      dir_out = dir_out, basics = basics, sampletable = sampletable,
      conda_path = conda_path["psichomics"], conda_env = conda_env["psichomics"],
      novel = novel, write_log = write_log, ...
    )
  }

  asp
}

#' Run the command generated by \code{\link{cmd_asp}} and/or
#'     read in the results
#' @description This function is about to run each workflow and read in results.
#' @name run_asp
#' @param asp An \code{\linkS4class{ASP}} object.
#' @param tool Any combination of "rmats", "cash", "leafcutter", "majiq",
#'     "suppa", "bandits", "psichomics", and "spladder".
#' @param run Default is TRUE. If you just want to check the command, set to
#'     FALSE, and the whole workflow will be printed.
#' @param read_results Default is TRUE.
#'     If read in results after workflows are done.
#' @param parallel Default is TRUE. If run multiple workflows simultaneously.
#'     The function \code{\link[future]{future}} is imposed.
#' @param block Default is TRUE. If wait for the end of running.
#' @import future
#' @export
#' @return An \code{\linkS4class{ASP}} object with results saved in results slot.
run_asp <- function(asp, tool, run = T, read_results = T, parallel = T, block = T) {
  if (run && length(tool) > 1) future::plan(multisession)

  if (run == F && read_results == T && parallel == T) {
    stop("please set parallel to FALSE when reading in results only!")
  }
  if (read_results && block == F && parallel) {
    warning(
      "you didn't wait for the finish but wanted to read in results,
      read_results has been set to FALSE automatically"
    )
    read_results <- F
  }
  if (length(tool) <= 1) {
    message("the number of tools you want to run is 1, and reset parallel to FALSE")
    parallel <- F
  }

  myreadin <- c(
    "rmats" = read_rmats,
    "cash" = read_cash,
    "leafcutter" = read_leafcutter,
    "majiq" = read_majiq,
    "spladder" = read_spladder,
    "bandits" = read_bandits,
    "suppa" = read_suppa,
    "psichomics" = read_psichomics
  )
  mydirs <- list(
    "rmats" = "",
    "cash" = "",
    "leafcutter" = c("ds", "juncs", "cluster"),
    "majiq" = c("build", "psi", "deltapsi", "voila"),
    "spladder" = "",
    "bandits" = "",
    "suppa" = c("tpm", "psi", "ds"),
    "psichomics" = ""
  )

  f <- vector("list", length = length(tool)) |> setNames(tool)

  for (t in tool) {
    cmds <- asp@cmds[[t]]
    message(glue::glue("{t} is running"))
    if (run) myCreatedir(glue::glue("{asp@dir_out}/{getOption('asp_tools')[t]}/{mydirs[[t]]}"))

    if (parallel) {
      cmd_in_R <- glue::glue("env_[t] %<-% {run_cmds(cmds, run)}", .open = "[", .close = "]")
      eval(parse(text = cmd_in_R))
      cmd_in_R <- glue::glue("f[[t]] <- future::futureOf(env_{t})")
      eval(parse(text = cmd_in_R))
    } else {
      run_cmds(cmds, run)
    }
  }

  ##### make run_asp blocked #####
  if (parallel && block && run) {
    tool_running <- tool
    while (length(tool_running) != 0) {
      for (t in tool_running) {
        if (future::resolved(f[[t]])) {
          eval(parse(text = glue::glue("env_{t}")))
          tool_running <- setdiff(tool_running, t)
        }
      }
    }
  }

  if (read_results && block) {
    asp@results <- asp@results[!names(asp@results) %in% tool]
    for (t in tool) {
      message(glue::glue("read in {t} results"))
      asp@results[[t]] <- myreadin[[t]](asp)
    }
  }

  asp
}

#' Query significant results
#' @description This function is about to extract the significant results.
#' @name sig_asp
#' @param tool Any combination of "rmats", "cash", "leafcutter", "majiq",
#'     "suppa", "bandits", "psichomics", and "spladder".
#' @param p_adj The threshold for the padj column.
#' @param p_value The threshold for the p column.
#' @param delta_psi The threshold for the dpsi column.
#' @param effect_size_col If results have no dpsi column, provide an
#'     alternative.
#' @param effect_size_threshod The threshold for the alternative
#'     effect size column.
#' @param filter_direction One of \code{c(">", "<")} indicating the significance
#'     direction.
#' @param statistic_col If results have no padj or p column, provide an
#'     alternative.
#' @param statistic_threshold The threshold for the alternative
#'     statistics column.
#' @export
#' @return An \code{\linkS4class{ASP}} object with significant results saved
#'     in sig_results slot.
sig_asp <- function(
  asp, tool, p_adj = 0.05, p_value = NULL, delta_psi = 0,
  statistic_col = NULL, statistic_threshold = 0, filter_direction = c(">", "<"),
  effect_size_col = NULL, effect_size_threshod = 0
) {
  ps <- list(
    p_adj = p_adj,
    p_value = p_value,
    delta_psi = delta_psi,
    statistic_col = statistic_col,
    statistic_threshold = statistic_threshold,
    filter_direction = filter_direction,
    effect_size_col = effect_size_col,
    effect_size_threshod = effect_size_threshod
  )

  for (t in tool) {
    dat <- asp@results[[t]]
    if (t == "majiq") dat <- asp@results[[t]]$voila_flatten
    if (t == "leafcutter") dat <- asp@results[[t]]$leafcutter_flatten
    asp@sig_results[[t]] <- querySig(dat = dat, ps)
  }

  asp
}

#' Prepare for plotting tracks of significant differential AS genes
#' @description This function is to load a txdb and an annotation database.
#' @name plot_pre_asp
#' @param asp An \code{\linkS4class{ASP}} object.
#' @param gtf If it's NULL, the gtf stored in \code{\linkS4class{ASP}} will be
#'     used.
#' @param annotation_db OrgDb.
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @export
#' @return An \code{\linkS4class{ASP}} object.
plot_pre_asp <- function(asp, gtf = NULL, annotation_db = NULL) {
  if (is.null(gtf)) gtf <- asp@gtf
  asp@txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
  if (is.null(annotation_db)) stop("please provide an annotation db, perhaps org.Hs.eg.db")
  asp@annotation_db <- annotation_db

  asp
}

#' Plot tracks of bam files to visualise significant differential AS genes
#' @description This function is to plot tracks
#'     using \code{\link[trackViewer]{trackViewer-package}}.
#' @name plot_asp
#' @param asp An \code{\linkS4class{ASP}} object.
#' @param tool Any combination of "rmats", "cash", "leafcutter", "majiq",
#'     "suppa", "bandits", "psichomics", and "spladder".
#' @param plot_samples Can be a integer (sampling in each condition)
#'     or a character vector (sample names).
#' @param dat If it's NULL, plot \code{top_n}.
#' @param top_n Default is 10.
#' @param dir_plot If it's NULL, pdf files will be deposited under the subfolder
#'     \emph{track} of \emph{dir_out} stored in \code{\linkS4class{ASP}}.
#' @param ... Accept two additional flags: \strong{extend_left} and
#'     \strong{extend_right}, mean how many based should be extended around the
#'     AS events.
#' @importFrom glue glue
#' @export
#' @return An invisible \code{\linkS4class{ASP}} object.
plot_asp <- function(asp, tool, plot_samples = NULL, dat = NULL, dir_plot = NULL, top_n = 10, ...) {
  basics_plot <- getPlotBasics(asp = asp, plot_samples = plot_samples, ...)

  for (t in tool) {
    if (is.null(dat)) {
      d <- asp@sig_results[[t]][1:top_n, ]
    } else {
      d <- dat[[t]]
    }
    ranges_lines <- ranges_lines(dat = d, tool = t, ...)

    if (is.null(dir_plot)) dir_plot <- glue::glue("{asp@dir_out}/tracks/{getOption('asp_tools')[t]}/")
    plot_tracks(
      dir_plot = dir_plot, txdb = asp@txdb, annotation_db = asp@annotation_db,
      plot_basics = basics_plot, ranges_lines = ranges_lines, ...
    )
    if (is.null(dat)) nm <- glue::glue("{asp@dir_out}/tracks/{getOption('asp_tools')[t]}/dat_plotted_{top_n}.txt")
    if (!is.null(dat)) nm <- glue::glue("{asp@dir_out}/tracks/{getOption('asp_tools')[t]}/dat_plotted_custom.txt")
    write.table(d, nm, sep = "\t", quote = F)
  }

  invisible(asp)
}

#' Intersect the results
#' @description This function is to retrieve intersections of significant
#'     results of different workflows and to visualize it.
#' @name intersect_asp
#' @param asp An \code{\linkS4class{ASP}} object.
#' @param tool Any combination of "rmats", "cash", "leafcutter", "majiq",
#'     "suppa", "bandits", "psichomics", and "spladder".
#'     Or use "all" to retrieve the intersections of all significant results
#'     stored in your \code{\linkS4class{ASP}} object.
#' @param type One of "SE", "MXE", "A3SS", "A5SS", "RI", "SME", "ALE" and "AFE".
#'     Or use "all" to retrieve the intersections of all AS types.
#' @param width pdf width.
#' @param height pdf height.
#' @param dir_plot If it's NULL, pdf files will be deposited under the
#'     \emph{dir_out} folder stored in \code{\linkS4class{ASP}}.
#' @importFrom dplyr %>% filter select mutate
#' @importFrom yyplot ggvenn
#' @importFrom ggplot2 guides theme element_blank theme_void
#' @importFrom ggsci scale_fill_npg
#' @importFrom tibble column_to_rownames as_tibble rownames_to_column
#' @importFrom UpSetR upset
#' @importFrom ggplotify as.ggplot
#' @importFrom ggimage geom_subview
#' @importFrom glue glue
#' @export
#' @return An \code{\linkS4class{ASP}} object with intersection tibble
#'     stored in the intersection slot.
intersect_asp <- function(asp, tools = "all", type = "all", dir_plot = NULL, width = 10, height = 7) {
  type <- type[1]
  results <- asp@sig_results
  if (tools != "all") results <- results[tools]
  if (type != "all") {
    results <- lapply(results, function(x) {
      tmp <- dplyr::filter(x, as_type %in% type)
      if (nrow(tmp) == 0) tmp <- NULL
      tmp
    })
    results <- results[!sapply(results, is.null)]
  }

  ds_gene_all <- sapply(results, function(x) {x$gene_symbol}) |>
    unlist() |>
    unique() |>
    na.omit()
  dat_ds_gene <- data.frame(matrix(0, nrow = length(ds_gene_all), ncol = length(results)))
  colnames(dat_ds_gene) <- names(results)
  dat_ds_gene$ds_gene <- ds_gene_all
  dat_ds_gene <- dplyr::select(dat_ds_gene, ds_gene, everything())
  for (t in names(results)) {
    ds_tbl <- table(results[[t]]$gene_symbol) %>%
      {
        .[match(dat_ds_gene$ds_gene, names(.))] -> .
        .
      } %>%
      as.numeric()
    dat_ds_gene[, t] <- ds_tbl
  }
  dat_ds_gene <- dat_ds_gene |> tibble::column_to_rownames("ds_gene")
  dat_ds_gene[is.na(dat_ds_gene)] <- 0
  dat_intersection <- dat_ds_gene
  dat_ds_gene[dat_ds_gene != 0] <- 1
  ##### venn plot #####
  yyplot::ggvenn(dat_ds_gene, alpha = 0.7) + ggplot2::guides(fill = "none") +
    ggsci::scale_fill_npg() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    ) -> p1
  ##### upset plot #####
  UpSetR::upset(
    dat_ds_gene, order.by = "freq", matrix.color = "#87CEFA", nsets = length(results), nintersects = NA,
    main.bar.color = "#424242", # sets.bar.color = color_sw[c(2, 1, 8, 4, 5, 6, 7, 3)],
    mainbar.y.label = "AS gene intersection", sets.x.label = "number of AS gene",
    point.size = 1.25, line.size = 0.6, mb.ratio = c(0.7, 0.3),
    att.pos = "top", shade.color = "black",
    shade.alpha = 0.1, matrix.dot.alpha = 0.5, text.scale = 0.8
  ) |>
    ggplotify::as.ggplot() -> p2
  pp <- p2 + ggimage::geom_subview(
    subview = p1 + ggplot2::theme_void(), x = .75, y = .7, w = .5, h = .5,
  )
  if (is.null(dir_plot)) dir_plot <- glue::glue("{asp@dir_out}/")
  pdf(
    glue::glue("{dir_plot}/intersection_{type}.pdf"),
    width = width, height = height
  )
  print(pp)
  dev.off()

  tmp <- rowSums(dat_ds_gene)
  asp@intersection[[type]] <- dat_intersection %>%
    dplyr::mutate(
      frequency = rowSums(.), frequency_unique = tmp
    ) |>
    tibble::rownames_to_column("gene") |>
    tibble::as_tibble()

  asp
}




