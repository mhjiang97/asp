#' rMATS
#' @description Generating commands of rMATS workflow.
#' @name rmats
#' @aliases rmats
#' @param dir_out,basics,nproc,gtf,novel,conda_path,conda_env,write_log
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param variable_read_length Default is FALSE.
#' @param paired_stats Default is FALSE.
#' @param cstat Default is 0.0001.
#' @param allow_clipping Default is TRUE.
#' @param ... nothing
#' @importFrom glue glue
#' @references \url{https://github.com/Xinglab/rmats-turbo}.
#' @return A list containing rMATS workflow.
rmats <- function(
  dir_out, basics, nproc, gtf, novel, conda_path, conda_env, write_log,
  variable_read_length = F, paired_stats = F, cstat = 0.0001, allow_clipping = T, ...
) {
  # myCreatedir(glue("{dir_out}/rMATS"))
  readLength <- basics$read_length_most
  library_type <- basics$library_type_most
  t <- ifelse(library_type == "paired-end", "paired", "single")
  strandedness <- basics$strandedness_most
  if (strandedness == "no") libType <- "fr-unstranded"
  if (strandedness == "yes") libType <- "fr-secondstrand"
  if (strandedness == "reverse") libType <- "fr-firststrand"
  b1 <- paste(path.expand(basics$bams_1), collapse = ",")
  b2 <- paste(path.expand(basics$bams_2), collapse = ",")

  cmds <- list()
  cmds[["generating input files"]] <- c(
    glue::glue('echo "{b1}\n" > {dir_out}/rMATS/{basics$condition_1}.txt'),
    glue::glue('echo "{b2}\n" > {dir_out}/rMATS/{basics$condition_2}.txt')
  )
  if (dir.exists(glue::glue("{dir_out}/rMATS/tmp/"))) cmds[["remove the tmp dir"]] <- glue::glue("rm -rf {dir_out}/rMATS/tmp/*")
  cmds[["running rMATS"]] <- glue::glue(
    "rmats.py \\
    --b1 {dir_out}/rMATS/{basics$condition_1}.txt \\
    --b2 {dir_out}/rMATS/{basics$condition_2}.txt \\
    --od {dir_out}/rMATS/ \\
    --gtf {gtf} \\
    --readLength {readLength} \\
    -t {t} \\
    --nthread {nproc} \\
    --cstat {cstat} \\
    --tstat {nproc} \\
    --libType {libType} \\
    --tmp {dir_out}/rMATS/tmp/"
  )
  if (!is.null(conda_env) && !is.na(conda_env)) cmds[["running rMATS"]] <- paste(glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[["running rMATS"]])
  if (novel) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--novelSS")
  if (variable_read_length) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--variable-read-length")
  if (paired_stats) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--paired-stats")
  if (allow_clipping) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--allow-clipping")
  if (write_log) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], glue::glue("1>{dir_out}/rMATS/_log 2>&1"))

  cmds
}

#' CASH
#' @description Generating commands of CASH workflow.
#' @name cash
#' @aliases cash
#' @param dir_out,basics,gtf,novel,conda_path,conda_env,write_log
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param cash_jar No default. The path to cash.jar.
#' @param MergePval Default is "G".
#' @param ram Default is "100g".
#' @param ... nothing
#' @importFrom glue glue
#' @references \url{https://soft.novelbio.com/cash/}.
#' @return A list containing CASH workflow.
cash <- function(
  dir_out, basics, gtf, novel, conda_path, conda_env, write_log,
  cash_jar, MergePval = "G", ram = "100g", ...
) {
  strandedness <- basics$strandedness_most
  if (strandedness == "no") StrandSpecific <- "NONE"
  if (strandedness == "yes") StrandSpecific <- "F"
  if (strandedness == "reverse") StrandSpecific <- "R"
  SpliceCons <- ifelse(novel, "True", "False")
  b1 <- paste(path.expand(basics$bams_1), collapse = ",")
  b2 <- paste(path.expand(basics$bams_1), collapse = ",")

  cmds <- list()

  cmds[["running cash"]] <- glue::glue(
    "java -jar -Xmx{ram} -Xms{ram} {cash_jar} \\
    --Case:{basics$condition_1} {b1} \\
    --Control:{basics$condition_2} {b2} \\
    --MergePval {MergePval} \\
    --GTF {gtf} \\
    --Output {dir_out}/CASH/cash_{MergePval} \\
    --SpliceCons {SpliceCons} \\
    --StrandSpecific {StrandSpecific}"
  )

  if (!is.null(conda_env) && !is.na(conda_env)) cmds[["running cash"]] <- paste(glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[["running cash"]])
  if (write_log) cmds[["running cash"]] <- paste(cmds[["running cash"]], glue::glue("1>{dir_out}/CASH/_log 2>&1"))

  cmds
}

#' MAJIQ
#' @description Generating commands of MAJIQ workflow.
#' @name majiq
#' @aliases majiq
#' @param
#' dir_out,basics,sampletable,nproc,gff,novel,fa,conda_path,conda_env,write_log
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param show_all Default is FALSE.
#' @param threshold Default is 0.2.
#' @param ... nothing
#' @importFrom glue glue
#' @references \url{https://majiq.biociphers.org}.
#' @return A list containing MAJIQ workflow.
majiq <- function(
  dir_out, basics, sampletable, nproc, gff, novel, fa, conda_path, conda_env, write_log,
  show_all = F, threshold = 0.2, ...
) {
  strandedness <- basics$strandedness_most
  if (strandedness == "no") strandness <- "None"
  if (strandedness == "yes") strandness <- "forward"
  if (strandedness == "reverse") strandness <- "reverse"
  bams_1 <- gsub("\\.bam", "", basics$bams_1)
  bams_2 <- gsub("\\.bam", "", basics$bams_2)
  b1 <- paste(basename(bams_1), collapse = ",")
  b2 <- paste(basename(bams_2), collapse = ",")
  bamdirs <- paste(path.expand(dirname(sampletable$files_bam)), collapse = ",")
  bams_exception <- sampletable$files_bam[sampletable$strandedness != strandedness]
  if (!all(is.na(bams_exception))) {
    b_e <- basename(gsub("\\.bam", "", bams_exception))
    s_e <- strandednesses_exception <-
      sampletable$strandedness[sampletable$strandedness != strandedness]
    for (i in 1:length(strandednesses_exception)) {
      if (strandednesses_exception[i] == "no") s_e[i] <- "None"
      if (strandednesses_exception[i] == "yes") s_e[i] <- "forward"
      if (strandednesses_exception[i] == "reverse") s_e[i] <- "reverse"
    }
    optional <- paste(glue("{b_e}=strandness:{s_e}"), collapse = ",\\n")
  }

  cmds <- list()

  cmds[["generating config"]] <- glue::glue(
    'echo "[info]\nbamdirs={bamdirs}\\
    \ngenome={fa}\nstrandness={strandness}\n\\
    [experiments]\n{basics$condition_1}={b1}\n{basics$condition_2}={b2}" > {dir_out}/MAJIQ/config'
  )
  if (!all(is.na(bams_exception))) {
    cmds[["generating config"]] <- glue::glue(
      'echo "[info]\nbamdirs={bamdirs}\\
      \ngenome={fa}\nstrandness={strandness}\n\\
      [experiments]\n{basics$condition_1}={b1}\n{basics$condition_2}={b2}\n\\
      [optional]\n{optional}" > {dir_out}/MAJIQ/config'
    )
  }

  cmds[["building majiq"]] <- glue::glue(
    "majiq build {gff} \\
    -c {dir_out}/MAJIQ/config \\
    -j {nproc} \\
    -o {dir_out}/MAJIQ/build/"
  )
  if (!novel) cmds[["building majiq"]] <- paste(cmds[["building majiq"]], "--disable-denovo --disable-denovo-ir")

  m1 <- gsub("\\.bam", "\\.majiq", basename(basics$bams_1))
  m2 <- gsub("\\.bam", "\\.majiq", basename(basics$bams_2))
  majiqs_1 <- glue::glue("{dir_out}/MAJIQ/build/{m1}") |> paste(collapse = " ")
  majiqs_2 <- glue::glue("{dir_out}/MAJIQ/build/{m2}") |> paste(collapse = " ")
  cmds[["calculating psi"]] <- glue::glue(
    "majiq psi \\
    {c(majiqs_1, majiqs_2)} \\
    -j {nproc} \\
    -o {dir_out}/MAJIQ/psi \\
    -n {c(basics$condition_1, basics$condition_2)}"
  )

  cmds[["calculating deltapsi"]] <- glue::glue(
    "majiq deltapsi \\
    -grp1 {majiqs_1} \\
    -grp2 {majiqs_2} \\
    -j {nproc} \\
    -o {dir_out}/MAJIQ/deltapsi \\
    -n {basics$condition_1} {basics$condition_2}"
  )


  cmds[["generating voila tsv"]] <- glue::glue(
    "voila tsv \\
    --threshold {threshold} \\
    {dir_out}/MAJIQ/build/splicegraph.sql \\
    {dir_out}/MAJIQ/deltapsi/{basics$condition_1}_{basics$condition_2}.deltapsi.voila \\
    -f {dir_out}/MAJIQ/voila/{basics$condition_1}_{basics$condition_2}_{threshold}.tsv"
  )

  if (show_all) {
    cmds[["generating voila tsv"]] <- glue::glue(
      "voila tsv \\
      --threshold {threshold} \\
      {dir_out}/MAJIQ/build/splicegraph.sql \\
      {dir_out}/MAJIQ/deltapsi/{basics$condition_1}_{basics$condition_2}.deltapsi.voila \\
      -f {dir_out}/MAJIQ/voila/{basics$condition_1}_{basics$condition_2}_{threshold}_showall.tsv \\
      --show-all"
    )
  }

  if (write_log) {
    cmds[["building majiq"]] <- paste(
      cmds[["building majiq"]], glue::glue("1>{dir_out}/MAJIQ/build/_log 2>&1")
    )
    cmds[["calculating psi"]] <- paste(
      cmds[["calculating psi"]], glue::glue("1>{dir_out}/MAJIQ/psi/_log.{c(basics$condition_1, basics$condition_2)} 2>&1")
    )
    cmds[["calculating deltapsi"]] <- paste(
      cmds[["calculating deltapsi"]], glue::glue("1>{dir_out}/MAJIQ/deltapsi/_log 2>&1")
    )
    cmds[["generating voila tsv"]] <- paste(
      cmds[["generating voila tsv"]], glue::glue("1>{dir_out}/MAJIQ/voila/_log 2>&1")
    )
  }

  if (!is.null(conda_env) && !is.na(conda_env)) {
    for (
      s in c(
        "building majiq", "calculating psi", "calculating deltapsi", "generating voila tsv"
      )
    ) {
      cmds[[s]] <- paste(
        glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[[s]]
      )
    }
  }

  cmds
}

#' LeafCutter
#' @description Generating commands of LeafCutter workflow.
#' @name leafcutter
#' @aliases leafcutter
#' @param
#' dir_out,sampletable,nproc,np,gtf,conda_path,conda_env,write_log,parallel
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param exon_file Any path is okay but should ends with .gz.
#'     Keep NULL will make the file be generated as
#'     \emph{dir_out/LeafCutter/exon_file.gz}.
#' @param leafcutter_dir No default. The path to
#'     your downloaded LeafCutter directory.
#' @param ... nothing
#' @importFrom glue glue
#' @importFrom utils write.table
#' @importFrom stats setNames
#' @importFrom futile.logger flog.info
#' @references \url{https://davidaknowles.github.io/leafcutter/}.
#' @return A list containing LeafCutter workflow.
leafcutter <- function(
  dir_out, sampletable, nproc, np, gtf, conda_path, conda_env, write_log, parallel,
  exon_file = NULL, leafcutter_dir, ...
) {
  if (is.null(exon_file)) exon_file <- glue::glue("{dir_out}/LeafCutter/exon_file.gz")

  cmds <- list()
  if (!file.exists(exon_file)) {
    myCreatedir(glue::glue("{dirname(exon_file)}"))
    cmds[["generating exon file"]] <- glue::glue(
      "Rscript {leafcutter_dir}/scripts/gtf_to_exons.R {gtf}.gz {exon_file}"
    )
  }

  strandedness <- sampletable$strandedness |>
    stats::setNames(sampletable$samples)
  s <- vector("character", nrow(sampletable)) |>
    stats::setNames(sampletable$samples)
  for (i in names(s)) {
    if (strandedness[i] == "no") s[i] <- "0"
    if (strandedness[i] == "yes") s[i] <- "2"
    if (strandedness[i] == "reverse") s[i] <- "1"
  }

  futile.logger::flog.info("Generating groups file used when analyzing differential intron excision")
  myCreatedir(glue::glue("{dir_out}/LeafCutter/ds/"))
  utils::write.table(
    sampletable[, c("samples", "conditions")],
    glue::glue("{dir_out}/LeafCutter/ds/groups_file.txt"),
    row.names = F, col.names = F, sep = "\t", quote = F
  )

  cmds[["converting bam to junc"]] <- glue::glue(
    "regtools junctions extract -s {s} {sampletable$files_bam} -t XS \\
    >{dir_out}/LeafCutter/juncs/{sampletable$samples}.junc && \\
    echo {dir_out}/LeafCutter/juncs/{sampletable$samples}.junc \\
    >>{dir_out}/LeafCutter/juncs/juncfiles.txt"
  )

  cmds[["intron clustering"]] <- glue::glue(
    "python2 {leafcutter_dir}/clustering/leafcutter_cluster_regtools.py \\
    -j {dir_out}/LeafCutter/juncs/juncfiles.txt \\
    -r {dir_out}/LeafCutter/cluster/ \\
    -o leafcutter \\
    -s \\
    -k"
  )

  cmds[["analyzing differential intron excision"]] <- glue::glue(
    "{leafcutter_dir}/scripts/leafcutter_ds.R \\
    --num_threads {nproc} \\
    --seed 1 \\
    --exon_file {exon_file} \\
    -o {dir_out}/LeafCutter/ds/leafcutter \\
    {dir_out}/LeafCutter/cluster/leafcutter_perind_numers.counts.gz \\
    {dir_out}/LeafCutter/ds/groups_file.txt"
  )

  if (min(table(sampletable$conditions)) < 5) cmds[["analyzing differential intron excision"]] <- paste(cmds[["analyzing differential intron excision"]], glue::glue("-i {min(table(sampletable$conditions))}"))
  if (min(table(sampletable$conditions)) < 3) cmds[["analyzing differential intron excision"]] <- paste(cmds[["analyzing differential intron excision"]], glue::glue("-g {min(table(sampletable$conditions))}"))

  if (write_log) {
    cmds[["intron clustering"]] <- paste(
      cmds[["intron clustering"]],
      glue::glue("1>{dir_out}/LeafCutter/cluster/_log 2>&1")
    )
    cmds[["analyzing differential intron excision"]] <- paste(
      cmds[["analyzing differential intron excision"]],
      glue::glue("1>{dir_out}/LeafCutter/ds/_log 2>&1")
    )
  }
  if (!is.null(conda_env) && !is.na(conda_env)) {
    for (
      s in c(
        "generating exon file", "converting bam to junc", "intron clustering", "analyzing differential intron excision"
      )
    ) {
      cmds[[s]] <- paste(
        glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[[s]]
      )
    }
  }
  if (parallel) cmds[["converting bam to junc"]] <- myParallel2(cmds[["converting bam to junc"]], n = np)

  cmds
}

#' SplAdder
#' @description Generating commands of SplAdder workflow.
#' @name spladder
#' @aliases spladder
#' @param dir_out,basics,gtf,novel,conda_path,conda_env,write_log
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param confidence_level Default is 3.
#' @param ... nothing
#' @importFrom glue glue
#' @references \url{https://spladder.readthedocs.io/en/latest/}.
#' @return A list containing SplAdder workflow.
spladder <- function(
  dir_out, basics, gtf, novel, conda_path, conda_env, write_log, parallel, np,
  confidence_level = 3, ...
) {
  readlen <- basics$read_length_most
  b1 <- paste(path.expand(basics$bams_1), collapse = ",")
  b2 <- paste(path.expand(basics$bams_2), collapse = ",")

  cmds <- list()

  if (T) {
    cmds[["generating single graph"]] <- glue::glue(
      "spladder build \\
      -o {dir_out}/SplAdder \\
      -a {gtf} \\
      -b {c(basics$bams_1, basics$bams_2)} \\
      --merge-strat single \\
      -n {readlen} \\
      -c {confidence_level} \\
      --no-extract-ase \\
      --validate-sg"
    )
    cmds[["merging graphs"]] <- glue::glue(
      "spladder build \\
      -o {dir_out}/SplAdder \\
      -a {gtf} \\
      -b {paste(b1, b2, sep = ',')} \\
      --merge-strat merge_graphs \\
      -n {readlen} \\
      -c {confidence_level} \\
      --no-extract-ase \\
      --validate-sg"
    )
    cmds[["quantifying"]] <- glue::glue(
      "spladder build \\
      -o {dir_out}/SplAdder \\
      -a {gtf} \\
      -b {c(basics$bams_1, basics$bams_2)} \\
      --merge-strat merge_graphs \\
      --quantify-graph \\
      --qmode single \\
      -n {readlen} \\
      -c {confidence_level} \\
      --no-extract-ase \\
      --validate-sg"
    )
    cmds[["collecting"]] <- glue::glue(
      "spladder build \\
      -o {dir_out}/SplAdder \\
      -a {gtf} \\
      -b {paste(b1, b2, sep = ',')} \\
      --merge-strat merge_graphs \\
      --quantify-graph \\
      --qmode collect \\
      -n {readlen} \\
      -c {confidence_level} \\
      --no-extract-ase \\
      --validate-sg"
    )
    cmds[["calling events"]] <- glue::glue(
      "spladder build \\
      -o {dir_out}/SplAdder \\
      -a {gtf} \\
      -b {paste(b1, b2, sep = ',')} \\
      --event-types exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons \\
      --output-txt \\
      -n {readlen} \\
      -c {confidence_level} \\
      --validate-sg"
    )
    if (!novel) {
      cmds[["generating single graph"]] <- paste(cmds[["generating single graph"]], "--no-insert-ir --no-insert-es --no-insert-ni")
      cmds[["merging graphs"]] <- paste(cmds[["merging graphs"]], "--no-insert-ir --no-insert-es --no-insert-ni")
      cmds[["quantifying"]] <- paste(cmds[["quantifying"]], "--no-insert-ir --no-insert-es --no-insert-ni")
      cmds[["collecting"]] <- paste(cmds[["collecting"]], "--no-insert-ir --no-insert-es --no-insert-ni")
      cmds[["calling events"]] <- paste(cmds[["calling events"]], "--no-insert-ir --no-insert-es --no-insert-ni")
    }
    if (write_log) {
      cmds[["generating single graph"]] <- paste(cmds[["generating single graph"]], glue::glue("1>{dir_out}/SplAdder/_log.build.{c(basics$samples_1, basics$samples_2)} 2>&1"))
      cmds[["merging graphs"]] <- paste(cmds[["merging graphs"]], glue::glue("1>{dir_out}/SplAdder/_log.merge_graphs 2>&1"))
      cmds[["quantifying"]] <- paste(cmds[["quantifying"]], glue::glue("1>{dir_out}/SplAdder/_log.quantification.{c(basics$samples_1, basics$samples_2)} 2>&1"))
      cmds[["collecting"]] <- paste(cmds[["collecting"]], glue::glue("1>{dir_out}/SplAdder/_log.collect 2>&1"))
      cmds[["calling events"]] <- paste(cmds[["calling events"]], glue::glue("1>{dir_out}/SplAdder/_log.event_calling 2>&1"))
    }
  }

  cmds[["testing alternative splicing events"]] <- glue::glue(
    "spladder test \\
    --readlen {readlen} \\
    --conditionA {b1} \\
    --conditionB {b2} \\
    --labelA {basics$condition_1} \\
    --labelB {basics$condition_2} \\
    --outdir {dir_out}/SplAdder"
  )
  if (write_log) cmds[["testing alternative splicing events"]] <- paste(cmds[["testing alternative splicing events"]], glue::glue("1>{dir_out}/SplAdder/_log.test 2>&1"))

  if (!is.null(conda_env) && !is.na(conda_env)) {
    for (
      s in c(
        "generating single graph", "merging graphs", "quantifying", "collecting", "calling events", "testing alternative splicing events"
      )
    ) {
      cmds[[s]] <- paste(
        glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[[s]]
      )
    }
  }

  if (parallel) {
    cmds[["generating single graph"]] <- myParallel2(cmds[["generating single graph"]], n = np)
    cmds[["quantifying"]] <- myParallel2(cmds[["quantifying"]], n = np)
  }

  cmds
}

#' BANDITS
#' @description Generating commands of BANDITS workflow.
#' @name bandits
#' @aliases bandits
#' @param
#' dir_out,sampletable,basics,nproc,conda_path,conda_env,write_log,tx2gene
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param infer_prior Default is TRUE.
#' @param max_genes_per_group Default is 50.
#' @param ... nothing
#' @importFrom glue glue
#' @importFrom utils write.table
#' @references
#' \url{https://bioconductor.org/packages/release/bioc/html/BANDITS.html}.
#' @return A list containing BANDITS workflow.
bandits <- function(
  dir_out, sampletable, basics, nproc, conda_path, conda_env, write_log, tx2gene,
  infer_prior = T, max_genes_per_group = 50, ...
) {
  myCreatedir(glue::glue("{dir_out}/BANDITS/"))
  write.table(sampletable, glue::glue("{dir_out}/BANDITS/sampleTable.txt"), quote = F, sep = "\t", row.names = F)
  cmds <- list()
  rscript <- glue::glue(
    '
    library(BANDITS)
    library(tximport)
    library(dplyr)
    library(glue)

    infer_prior <- {{infer_prior}}
    max_genes_per_group <- {{max_genes_per_group}}


    tx2gene <- read.delim("{{tx2gene}}", header = F) # |> setNames(c("gene_id", "transcript_id"))

    sampleTable <- read.delim("{{dir_out}}/BANDITS/sampleTable.txt")
    samples_design <- sampleTable |>
      dplyr::rename(sample_id = samples, group = conditions) |>
      dplyr::select(sample_id, group)

    quant_files <- vector("character", length = nrow(samples_design)) |>
      setNames(samples_design$sample_id)
    for (s in names(quant_files)) {
      quant_files[s] <- list.files(
        glue("{sampleTable$dirs_salmon[sampleTable$samples == s]}/"),
        pattern = "quant.sf*",
        full.names = T
      )
    }

    txi <- tximport(quant_files, type = "salmon", txOut = T)
    eff_len = eff_len_compute(x_eff_len = txi$length)

    equiv_classes_files <- vector("character", length = nrow(samples_design)) |>
      setNames(samples_design$sample_id)
    for (s in names(equiv_classes_files)) {
      equiv_classes_files[s] <- list.files(
        glue("{sampleTable$dirs_salmon[sampleTable$samples == s]}/aux_info/"),
        pattern = "eq_classes.txt*",
        full.names = T
      )
    }

    input_data <- create_data(
      salmon_or_kallisto = "salmon",
      gene_to_transcript = tx2gene,
      salmon_path_to_eq_classes = equiv_classes_files,
      eff_len = eff_len,
      n_cores = nrow(samples_design),
      max_genes_per_group = {{max_genes_per_group}}
    )

    precision <- NULL
    if (infer_prior) {
      precision <- prior_precision(
        gene_to_transcript = tx2gene,
        transcript_counts = txi$counts,
        n_cores = nrow(samples_design),
        transcripts_to_keep = NULL
      )
    }

    set.seed(1)
    results <- test_DTU(
      BANDITS_data = input_data,
      precision = precision$prior,
      samples_design = samples_design,
      group_col_name = "group",
      R = 10^4,
      burn_in = 2*10^3,
      n_cores = {{nproc}},
      gene_to_transcript = tx2gene
    )

    results_gene <- top_genes(results)
    results_transcript <- top_transcripts(results)

    write.table(
      results_gene, "{{dir_out}}/BANDITS/results_gene.txt",
      row.names = F, sep = "\\t", quote = F
    )
    write.table(
      results_transcript, "{{dir_out}}/BANDITS/results_transcript.txt",
      row.names = F, sep = "\\t", quote = F
    )
    save.image("{{dir_out}}/BANDITS/bandits.RData", compress = T)
  ',
  .open = "{{", .close = "}}"
)

  cmds[["generating the bandits rscript"]] <- glue::glue("echo '{rscript}' | sed 's/..//' > {dir_out}/BANDITS/BANDITS.R")
  cmds[["running bandits"]] <- glue::glue("Rscript {dir_out}/BANDITS/BANDITS.R")
  if (!is.null(conda_env) && !is.na(conda_env)) cmds[["running bandits"]] <- paste(glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[["running bandits"]])
  if (write_log) cmds[["running bandits"]] <- paste(cmds[["running bandits"]], glue::glue("1>{dir_out}/BANDITS/_log 2>&1"))

  cmds
}

#' SUPPA
#' @description Generating commands of SUPPA workflow.
#' @name suppa
#' @aliases suppa
#' @param dir_out,basics,np,gtf,conda_path,conda_env,parallel,write_log
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param dir_events Any directory is ok.
#'     The directory holding SUPPA ioe files. Keep NULL will make the ioe files
#'     be generated in the \emph{dir_out/SUPPA/events}.
#' @param ... nothing
#' @importFrom glue glue
#' @importFrom tximport tximport
#' @importFrom purrr flatten_chr
#' @importFrom utils write.table
#' @references \url{https://github.com/comprna/SUPPA}.
#' @return A list containing SUPPA workflow.
suppa <- function(
  dir_out, basics, gtf, np, conda_path, conda_env, parallel, write_log,
  dir_events = NULL, ...
) {
  if (is.null(dir_events)) dir_events <- glue::glue("{dir_out}/SUPPA/events/")
  events <- c("A3", "A5", "AF", "AL", "MX", "RI", "SE")

  cmds <- list()
  if (length(list.files(dir_events, pattern = glue::glue(".*\\.ioe"), full.names = T)) < length(events)) {
    myCreatedir(dir_events)
    cmds[["generating events"]] <- glue::glue(
      "suppa.py generateEvents -i {gtf} -o {dir_events}/events -f ioe -e SE SS MX RI FL"
    )
    run_cmds(cmds)
    cmds <- list()
  }

  ioe_suppa <- vector("character", length = length(events))
  names(ioe_suppa) <- events
  for (e in events) ioe_suppa[e] <- list.files(dir_events, pattern = glue::glue(".*{e}.*\\.ioe"), full.names = T)
  files_tpm <- c(basics$salmons_1, basics$salmons_2)
  txi <- tximport::tximport(files_tpm, type = "salmon", txOut = T)
  tpm_1 <- txi$abundance[, basics$samples_1]
  tpm_2 <- txi$abundance[, basics$samples_2]
  utils::write.table(tpm_1, glue::glue("{dir_out}/SUPPA/tpm/{basics$condition_1}.tpm"), quote = F, sep = "\t")
  utils::write.table(tpm_2, glue::glue("{dir_out}/SUPPA/tpm/{basics$condition_2}.tpm"), quote = F, sep = "\t")

  calculate_psi <- vector("list", length = length(events))
  for (e in events) {
    calculate_psi[[e]] <- glue::glue(
      "suppa.py psiPerEvent \\
      -i {ioe_suppa[e]} \\
      -e {dir_out}/SUPPA/tpm/{c(basics$condition_1, basics$condition_2)}.tpm \\
      -o {dir_out}/SUPPA/psi/{c(basics$condition_1, basics$condition_2)}.{e}"
    )
    if (write_log) calculate_psi[[e]] <- paste(calculate_psi[[e]], glue::glue("1>{dir_out}/SUPPA/psi/_log.{c(basics$condition_1, basics$condition_2)}.{e} 2>&1"))
  }
  cmds[["calculating psi"]] <- purrr::flatten_chr(calculate_psi)

  perform_ds <- vector("list", length = length(events))
  for (e in events) {
    perform_ds[[e]] <- glue::glue(
      "suppa.py diffSplice \\
      -m empirical \\
      -p {dir_out}/SUPPA/psi/{basics$condition_1}.{e}.psi {dir_out}/SUPPA/psi/{basics$condition_2}.{e}.psi \\
      -e {dir_out}/SUPPA/tpm/{basics$condition_1}.tpm {dir_out}/SUPPA/tpm/{basics$condition_2}.tpm \\
      -i {ioe_suppa[e]} \\
      --gene-correction \\
      -o {dir_out}/SUPPA/ds/{basics$condition_1}_{basics$condition_2}.{e}"
    )
    if (write_log) perform_ds[[e]] <- paste(perform_ds[[e]], glue::glue("1>{dir_out}/SUPPA/ds/_log.{e} 2>&1"))
  }
  cmds[["performing differential splicing analysis"]] <- purrr::flatten_chr(perform_ds)

  if (!is.null(conda_env) && !is.na(conda_env)) {
    for (
      s in c(
        "calculating psi", "performing differential splicing analysis"
      )
    ) {
      cmds[[s]] <- paste(
        glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[[s]]
      )
    }
  }
  if (parallel) {
    cmds[["calculating psi"]] <- myParallel2(cmds[["calculating psi"]], n = np)
    cmds[["performing differential splicing analysis"]] <- myParallel2(cmds[["performing differential splicing analysis"]], n = np)
  }

  cmds
}

#' psichomics
#' @description Generating commands of psichomics workflow.
#' @name psichomics
#' @aliases psichomics
#' @param dir_out,basics,sampletable,novel,conda_path,conda_env,write_log
#'     These are obtained from \code{\linkS4class{ASP}} object.
#' @param anno_suppa The directory holding SUPPA event files.
#' @param prefix_suppa No default. The prefix of SUPPA event files.
#' @param anno_rmats The directory holding rMATS event files.
#' @param anno_miso The directory holding MISO event files.
#' @param prefix_miso No default. The prefix of MISO event files.
#' @param anno_vast The directory holding VAST-TOOLS event files.
#' @param prefix_vast No default. The prefix of VAST-TOOLS event files.
#' @param types_miso Should be like this: "SE,MXE,A5SS,A3SS,RI".
#' @param ... nothing
#' @importFrom glue glue
#' @importFrom utils write.table
#' @references
#' \url{https://bioconductor.org/packages/release/bioc/html/psichomics.html}.
#' @return A list containing psichomics workflow.
psichomics <- function(
  dir_out, basics, sampletable, conda_path, conda_env, novel, write_log,
  anno_suppa = NULL, anno_rmats = NULL, anno_miso = NULL, anno_vast = NULL,
  prefix_suppa = NULL, prefix_miso = NULL, prefix_vast = NULL, types_miso = NULL,
  ...
) {
  myCreatedir(glue::glue("{dir_out}/psichomics/"))
  utils::write.table(
    sampletable, glue::glue("{dir_out}/psichomics/sampleTable.txt"),
    quote = F, sep = "\t", row.names = F
  )
  cmds <- list()
  rscript <- glue::glue(
    '
    library(psichomics)

    anno_suppa <- {{anno_suppa}}
    anno_rmats <- {{anno_rmats}}
    anno_miso <- {{anno_miso}}
    anno_vast <- {{anno_vast}}
    prefix_suppa <- {{prefix_suppa}}
    prefix_miso <- {{prefix_miso}}
    prefix_vast <- {{prefix_vast}}
    types_miso <- {{unlist(strsplit(types_miso, split = ","))}}

    sampleTable <- read.delim("{dir_out}/psichomics/sampleTable.txt")

    strandedness <- {{basics$strandedness_most}}
    if (strandedness == "no") sn <- "unstranded"
    if (strandedness == "yes") sn <- "stranded"
    if (strandedness == "reverse") sn <- "stranded (reverse)"

    files_genequant <- {{c(basics$quants_star_1, basic$quants_star_2)}}
    files_junctionquant <- {{c(basics$juncs_star_1, basics$juncs_star_2)}}

    if (!is.null(anno_suppa) && !is.null(prefix_suppa)) {
      suppa <- parseSuppaAnnotation(anno_suppa, genome = prefix_suppa)
    } else {suppa <- NULL}
    if (!is.null(anno_rmats)) {
      mats <- parseMatsAnnotation(anno_rmats, novelEvents = {{novel}})
    } else {mats <- NULL}
    if (!is.null(anno_miso) && !is.null(prefix_miso)) {
      miso <- parseMisoAnnotation(anno_miso, genome = prefix_miso, types = types_miso)
    } else {miso <- NULL}
    if (!is.null(anno_vast) && !is.null(prefix_vast)) {
      vast <- parseVastToolsAnnotation(anno_vast, genome = prefix_vast)
    } else {vast <- NULL}
    annot <- prepareAnnotationFromEvents(suppa, vast, mats, miso)
    saveRDS(annot, file = "{{dir_out}}/psichomics/annot.rds")

    prepareGeneQuant(
      files_genequant, strandedness = sn,
      output = "{{dir_out}}/psichomics/psichomics_gene_counts.txt"
    )
    prepareJunctionQuant(
      files_junctionquant, output = "{{dir_out}}/psichomics/psichomics_junctions.txt"
    )
    write.table(
      dplyr::select(sampleTable, samples, conditions) |> dplyr::rename(`Sample ID` = samples),
      "{{dir_out}}/psichomics/psichomics_metadata.txt",
      sep = "\t", quote = F, row.names = F
    )

    data <- loadLocalFiles("{{dir_out}}/psichomics/")
    sampleInfo <- data[[1]]$`Sample metadata`
    geneExpr <- data[[1]]$`Gene expression`
    junctionQuant <- data[[1]]$`Junction quantification`
    if (F) {
      filter <- filterGeneExpr(geneExpr)
    } else {
      filter <- NULL
    }
    geneExprNorm <- normaliseGeneExpression(geneExpr, geneFilter = filter)
    psi <- quantifySplicing(annot, junctionQuant)

    s1VSs2 <- list(
      sampleTable$samples[sampleTable$conditions == {basics$condition_1}],
      sampleTable$samples[sampleTable$conditions == {basics$condition_2}]
    ) |>
    setNames(unique(sampleTable$conditions))
    diffSplicing <- diffAnalyses(psi, s1VSs2)
    write.table(
      diffSplicing, "{{dir_out}}/psichomics/results_tmp.txt",
      quote = F, sep = "\t"
    )
    write.table(
      psi, "{{dir_out}}/psichomics/psi.txt",
      quote = F, sep = "\t"
    )
    save.image("{{dir_out}}/psichomics/psichomics.RData", compress = T)
    ', .open = "{{", .close = "}}"
  )

  cmds[["generating the psichomics rscript"]] <- glue::glue('echo "{rscript}" > {dir_out}/psichomics/psichomics.R')
  cmds[["running psichomics"]] <- glue::glue("Rscript {dir_out}/psichomics/psichomics.R")
  cmds[["tidying output"]] <- glue::glue("awk '{if(NR==1){print \"\t\"$0}else{print $0}}' {{dir_out}}/psichomics/results_tmp.txt > {{dir_out}}/psichomics/results.txt", .open = "{{", .close = "}}")
  if (!is.null(conda_env) && !is.na(conda_env)) cmds[["running psichomics"]] <- paste(glue::glue("source {conda_path}/bin/activate {conda_env} &&"), cmds[["running psichomics"]])
  if (write_log) cmds[["running psichomics"]] <- paste(cmds[["running psichomics"]], glue::glue("1>{dir_out}/psichomics/_log 2>&1"))

  cmds
}





