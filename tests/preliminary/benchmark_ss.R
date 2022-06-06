library(glue)
library(dplyr)
library(purrr)
library(vroom)
source("~/projects/pkgs/R/asp/R/utils.R")
##### samples #####
samples_1 <- paste0("sample", 1:4)
samples_2 <- paste0("sample", 5:8)

suffix_bam <- "SortedByCoord.bam"
dir_bam <- "~/projects/AS/data/sim/rsem/novel_splice_site/star/"
bams_1 <- sprintf("%s/%s/%s.%s", dir_bam, samples_1, samples_1, suffix_bam)
bams_2 <- sprintf("%s/%s/%s.%s", dir_bam, samples_2, samples_2, suffix_bam)


suffix_fq <- "fq.gz"
dir_fq <- "~/projects/AS/data/sim/rsem/paper_gb/sim_1/fq/"
fqs_1 <- paste(
  sprintf("%s/%s_R1.%s", dir_fq, samples_1, suffix_fq),
  sprintf("%s/%s_R2.%s", dir_fq, samples_1, suffix_fq),
  sep = ";"
)
fqs_2 <- paste(
  sprintf("%s/%s_R1.%s", dir_fq, samples_2, suffix_fq),
  sprintf("%s/%s_R2.%s", dir_fq, samples_2, suffix_fq),
  sep = ";"
)

##### generate a sample table #####
sampleTable <- data.frame(
  samples = c(samples_1, samples_2),
  conditions = c(rep("s1", 4), rep("s2", 4)),
  files_bam = c(bams_1, bams_2),
  files_fq = c(fqs_1, fqs_2),
  read_lengths = rep(101, 8),
  library_types = rep("paired-end", 8),
  strandedness = rep("no", 8)
)

##### required reference and annotation files #####
gtf <- "~/projects/AS/data/sim/rsem/novel_splice_site/gencode.v32.annotation.primary_assembly.protein_coding2.trimed_ss.gtf"
gff <- "~/projects/AS/data/sim/rsem/novel_splice_site/gencode.v32.annotation.protein_coding2.trimed_ss.gff3"
fa <- "~/doc/reference/fa/GRCh38.primary_assembly.genome.fa"
fa_transcript <- "~/projects/AS/data/sim/rsem/novel_splice_site/gencode.v32.transcripts_ss.fa"
genome_version <- "hg38"

##### other parameters #####
basics <- getBasics(sampleTable)

dir_out <- "~/projects/AS/analysis/novel_splice_site/"
nproc <- 15
novel <- T
write_log <- T
parallel <- T
np <- 5


##### rMATS #####
myCreatedir(glue("{dir_out}/rMATS"))

variable_read_length <- F
paired_stats <- F
cstat <- 0.05

rL <- table(sampleTable$read_lengths) %>% .[. == max(.)] %>% names() %>% as.numeric()
library_type <- table(sampleTable$library_types) %>% .[. == max(.)] %>% names()
tt <- ifelse(library_type == "paired-end", "paired", "single")
strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") lT <- "fr-unstranded"
if (strand == "yes") lT <- "fr-secondstrand"
if (strand == "reverse") lT <- "fr-firststrand"
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]]
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]]
b1 <- paste(path.expand(bams_1), collapse = ",")
b2 <- paste(path.expand(bams_2), collapse = ",")

cmds <- list()

cmds[["generating input files"]] <- c(
  glue('echo -e "{b1}\\n" >{dir_out}/rMATS/{conditions[1]}.txt'),
  glue('echo -e "{b2}\\n" >{dir_out}/rMATS/{conditions[2]}.txt')
)

if (dir.exists(glue("{dir_out}/rMATS/tmp/"))) cmds[["remove files in tmp dir"]] <- glue("rm -rf {dir_out}/rMATS/tmp/*")

cmds[["running rMATS"]] <- glue(
  "rmats.py \\
  --b1 {dir_out}/rMATS/{conditions[1]}.txt \\
  --b2 {dir_out}/rMATS/{conditions[2]}.txt \\
  --od {dir_out}/rMATS/ \\
  --gtf {gtf} \\
  --readLength {rL} \\
  -t {tt} \\
  --nthread {nproc} \\
  --cstat {cstat} \\
  --tstat {nproc} \\
  --libType {lT} \\
  --tmp {dir_out}/rMATS/tmp/"
)

if (novel) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--novelSS")
if (variable_read_length) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--variable-read-length")
if (paired_stats) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], "--paired-stats")
if (write_log) cmds[["running rMATS"]] <- paste(cmds[["running rMATS"]], glue("1>{dir_out}/rMATS/_log 2>&1"))


##### CASH #####
myCreatedir(glue("{dir_out}/CASH/"))

cash <- "~/projects/as/software/cash_v2.2.0/cash.jar"
MergePval <- "G"
ram <- "100g"

strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") SS <- "NONE"
if (strand == "yes") SS <- "F"
if (strand == "reverse") SS <- "R"
SC <- ifelse(novel, "True", "False")
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]]
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]]
b1 <- paste(path.expand(bams_1), collapse = ",")
b2 <- paste(path.expand(bams_2), collapse = ",")

cmds <- list()

cmds[["running CASH"]] <- glue(
  "java -jar -Xmx{ram} -Xms{ram} {cash} \\
  --Case:{conditions[1]} {b1} \\
  --Control:{conditions[2]} {b2} \\
  --MergePval {MergePval} \\
  --GTF {gtf} \\
  --Output {dir_out}/CASH/cash \\
  --SpliceCons {SC} \\
  --StrandSpecific {SS}"
)

if (write_log) cmds[["running CASH"]] <- paste(cmds[["running CASH"]], glue("1>{dir_out}/CASH/_log 2>&1"))


##### LeafCutter #####
myCreatedir(glue('{dir_out}/LeafCutter/{c("ds", "juncs", "cluster")}/'))

method_junc <- "regtools"
exon_file <- "~/doc/leafcutter/gencode_v32_exons.txt.gz"

strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") ss <- "0"
if (strand == "yes") ss <- "2"
if (strand == "reverse") ss <- "1"
script_dir <- "~/projects/as/software/leafcutter"

message("generating groups file used when analyzing differential intron excision")
write.table(
  sampleTable[, c("samples", "conditions")],
  glue("{dir_out}/LeafCutter/ds/groups_file.txt"),
  row.names = F, col.names = F, sep = "\t", quote = F
)

cmds <- list()

cmds[["converting bam to junc"]] <- glue(
  "regtools junctions extract -s {ss} {sampleTable$files_bam} \\
  >{dir_out}/LeafCutter/juncs/{sampleTable$samples}.junc && \\
  echo {dir_out}/LeafCutter/juncs/{sampleTable$samples}.junc \\
  >>{dir_out}/LeafCutter/juncs/juncfiles.txt"
)
if (parallel) cmds[["converting bam to junc"]] <- myParallel(cmds[["converting bam to junc"]])

cmds[["intron clustering"]] <- glue(
  "python {script_dir}/clustering/leafcutter_cluster_regtools.py \\
  -j {dir_out}/LeafCutter/juncs/juncfiles.txt \\
  -r {dir_out}/LeafCutter/cluster/ \\
  -o leafcutter \\
  -s \\
  -k"
)

cmds[["analyzing differential intron excision"]] <- glue(
  "{script_dir}/scripts/leafcutter_ds.R --num_threads {nproc} \\
  --exon_file {exon_file} \\
  {dir_out}/LeafCutter/cluster/leafcutter_perind_numers.counts.gz \\
  {dir_out}/LeafCutter/ds/groups_file.txt \\
  -o {dir_out}/LeafCutter/ds/leafcutter"
)

if (min(table(sampleTable$conditions)) < 5) cmds[["analyzing differential intron excision"]] <- paste(cmds[["analyzing differential intron excision"]], glue("-i {min(table(sampleTable$conditions))}"))
if (min(table(sampleTable$conditions)) < 3) cmds[["analyzing differential intron excision"]] <- paste(cmds[["analyzing differential intron excision"]], glue("-g {min(table(sampleTable$conditions))}"))

if (write_log) {
  cmds[["intron clustering"]] <- paste(
    cmds[["intron clustering"]],
    glue("1>{dir_out}/LeafCutter/cluster/_log 2>&1")
  )
  cmds[["analyzing differential intron excision"]] <- paste(
    cmds[["analyzing differential intron excision"]],
    glue("1>{dir_out}/LeafCutter/ds/_log 2>&1")
  )
}


##### MAJIQ #####
myCreatedir(glue('{dir_out}/MAJIQ/{c("psi", "deltapsi", "voila")}/'))

show_all <- F
threshold <- 0.05

rl <- max(sampleTable$read_lengths)
strand <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strand == "no") ss <- "None"
if (strand == "yes") ss <- "forward"
if (strand == "reverse") ss <- "reverse"
conditions <- unique(sampleTable$conditions)
bams_1 <- sampleTable$files_bam[sampleTable$conditions == conditions[1]] %>% gsub("\\.bam", "", .)
bams_2 <- sampleTable$files_bam[sampleTable$conditions == conditions[2]] %>% gsub("\\.bam", "", .)
b1 <- paste(basename(bams_1), collapse = ",")
b2 <- paste(basename(bams_2), collapse = ",")
bd <- paste(path.expand(dirname(sampleTable$files_bam)), collapse = ",")

cmds <- list()

cmds[["generating config"]] <- glue(
  'echo -e "[info]\\nreadlen={rl}\\nbamdirs={bd}\\
  \\ngenome={fa}\\nstrandness={ss}\\n\\
  [experiments]\\n{conditions[1]}={b1}\\n{conditions[2]}={b2}" > {dir_out}/MAJIQ/config'
)

if (!dir.exists(glue("{dir_out}/MAJIQ/build/"))) dir.create(glue("{dir_out}/MAJIQ/build/"), recursive = T)

cmds[["building majiq"]] <- glue(
  "majiq build {gff} \\
  -c {dir_out}/MAJIQ/config \\
  -j {nproc} \\
  -o {dir_out}/MAJIQ/build/"
)

if (!novel) cmds[["building majiq"]] <- paste(cmds[["building majiq"]], "--disable-denovo --disable-denovo-ir")

m1 <- basename(sampleTable$files_bam[sampleTable$conditions == conditions[1]]) %>%
  gsub("\\.bam", "\\.majiq", .)
m2 <- basename(sampleTable$files_bam[sampleTable$conditions == conditions[2]]) %>%
  gsub("\\.bam", "\\.majiq", .)
files_majiq1 <- glue("{dir_out}/MAJIQ/build/{m1}") %>% paste(collapse = " ")
files_majiq2 <- glue("{dir_out}/MAJIQ/build/{m2}") %>% paste(collapse = " ")
cmds[["calculating psi"]] <- glue(
  "majiq psi \\
  {c(files_majiq1, files_majiq2)} \\
  -j {nproc} \\
  -o {dir_out}/MAJIQ/psi \\
  -n {conditions}"
)

cmds[["calculating deltapsi"]] <- glue(
  "majiq deltapsi \\
  -grp1 {files_majiq1} \\
  -grp2 {files_majiq2} \\
  -j {nproc} \\
  -o {dir_out}/MAJIQ/deltapsi \\
  -n {conditions[1]} {conditions[2]}"
)

cmds[["generating voila tsv"]] <- glue(
  "voila tsv \\
  --threshold {threshold} \\
  {dir_out}/MAJIQ/build/splicegraph.sql \\
  {dir_out}/MAJIQ/deltapsi/{conditions[1]}_{conditions[2]}.deltapsi.voila \\
  -f {dir_out}/MAJIQ/voila/{conditions[1]}_{conditions[2]}_{threshold}.tsv"
)

if (show_all) {
  cmds[["generating voila tsv"]] <- glue(
    "voila tsv \\
    --threshold {threshold} \\
    {dir_out}/MAJIQ/build/splicegraph.sql \\
    {dir_out}/MAJIQ/deltapsi/{conditions[1]}_{conditions[2]}.deltapsi.voila \\
    -f {dir_out}/MAJIQ/voila/{conditions[1]}_{conditions[2]}_{threshold}_showall.tsv \\
    --show-all"
  )
}

if (write_log) {
  cmds[["building majiq"]] <- paste(
    cmds[["building majiq"]],
    glue("1>{dir_out}/MAJIQ/build/_log 2>&1")
  )
  cmds[["calculating psi"]] <- paste(
    cmds[["calculating psi"]],
    glue("1>{dir_out}/MAJIQ/psi/_log.{conditions} 2>&1")
  )
  cmds[["calculating deltapsi"]] <- paste(
    cmds[["calculating deltapsi"]],
    glue("1>{dir_out}/MAJIQ/deltapsi/_log 2>&1")
  )
  cmds[["generating voila tsv"]] <- paste(
    cmds[["generating voila tsv"]],
    glue("1>{dir_out}/MAJIQ/voila/_log 2>&1")
  )
}


##### Whippet #####
myCreatedir(glue('{dir_out}/Whippet/{c("index", "quant", "delta", "merge")}/'))

conditions <- unique(sampleTable$conditions)

prefix_index <- "~/projects/AS/analysis/novel_splice_site/Whippet/index/index"
bam_min_reads <- 2

if (length(list.files(dirname(prefix_index), pattern = basename(prefix_index))) != 2) {
  myCreatedir(dirname(prefix_index))
  cmds <- list()
  if (novel) {
    cmds[["merging bam"]] <- glue(
      "samtools merge -@ {nproc} -l 9 {dir_out}/Whippet/merge/merged.bam {paste(sampleTable$files_bam, collapse = ' ')}"
    )
    cmds[["removing duplicates"]] <- glue(
      "samtools rmdup -S {dir_out}/Whippet/merge/merged.bam {dir_out}/Whippet/merge/merged.rmdup.bam"
    )
    cmds[["indexing rmdup bam"]] <- glue(
      "samtools index -@ {nproc} {dir_out}/Whippet/merge/merged.rmdup.bam"
    )
    cmds[["generating index"]] <- glue(
      "whippet-index.jl \\
      --fasta {fa}.gz \\
      --gtf {gtf}.gz \\
      --suppress-low-tsl \\
      --bam-both-novel \\
      --bam-min-reads {bam_min_reads} \\
      --bam {dir_out}/Whippet/merge/merged.rmdup.bam \\
      -x {prefix_index}"
    )
  } else {
    cmds[["generating index"]] <- glue(
      "whippet-index.jl \\
      --fasta {fa}.gz \\
      --gtf {gtf}.gz \\
      --suppress-low-tsl \\
      -x {prefix_index}"
    )
  }
  run_cmds(cmds)
}

cmds <- list()
quantify_fq <- vector("list", length = nrow(sampleTable))
names(quantify_fq) <- sampleTable$samples
for (s in names(quantify_fq)) {
  tmp_fq <- strsplit(sampleTable$files_fq[sampleTable$samples == s], ";") %>% unlist() %>% paste(collapse = " ")
  quantify_fq[[s]] <- glue(
    "whippet-quant.jl \\
    {tmp_fq} \\
    -x {prefix_index} \\
    -o {dir_out}/Whippet/quant/{s}"
  )
  if (write_log) quantify_fq[[s]] <- paste(quantify_fq[[s]], glue("1>{dir_out}/Whippet/quant/_log.{s} 2>&1"))
}
cmds[["quantifying fastq"]] <- flatten_chr(quantify_fq)

cmds[["quantifying fastq"]] <- paste(paste(cmds[["quantifying fastq"]], collapse = " & "), "& wait")

files_quant1 <- vector("character", length = length(which(sampleTable$conditions == conditions[1]))) %>%
  setNames(sampleTable$samples[sampleTable$conditions == conditions[1]])
for (s in names(files_quant1)) {
  files_quant1[s] <- list.files(glue("{dir_out}/Whippet/quant/"), pattern = glue("{s}.psi.gz"), full.names = T)
}
files_quant2 <- vector("character", length = length(which(sampleTable$conditions == conditions[2]))) %>%
  setNames(sampleTable$samples[sampleTable$conditions == conditions[2]])
for (s in names(files_quant2)) {
  files_quant2[s] <- list.files(glue("{dir_out}/Whippet/quant/"), pattern = glue("{s}.psi.gz"), full.names = T)
}

cmds[["comparing psi"]] <- glue(
  "whippet-delta.jl \\
  -a {paste(files_quant1, collapse = ',')} \\
  -b {paste(files_quant2, collapse = ',')} \\
  -o {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}"
)
if (write_log) cmds[["comparing psi"]] <- paste(cmds[["comparing psi"]], glue("1>{dir_out}/Whippet/delta/_log 2>&1"))

cmds[["tidy output"]] <- glue(
  "zcat {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}.diff.gz | \\
  cut -f 1,2,3,4,5,6,7,8,9,10,11 > {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}.diff2 | \\
  gzip {dir_out}/Whippet/delta/{conditions[1]}_{conditions[2]}.diff2"
)


##### psichomics #####
myCreatedir(glue("{dir_out}/psichomics/"))

library(psichomics)
strandedness <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
if (strandedness == "no") sn <- "unstranded"
if (strandedness == "yes") sn <- "stranded"
if (strandedness == "reverse") sn <- "stranded (reverse)"

files_genequant <- glue(
  "~/projects/AS/data/sim/rsem/novel_splice_site/star/sample{1:8}/ReadsPerGene.out.tab"
) %>% setNames(sampleTable$samples)
files_junctionquant <- glue(
  "~/projects/AS/data/sim/rsem/novel_splice_site/star/sample{1:8}/SJ.out.tab"
) %>% setNames(sampleTable$samples)
anno_suppa <- "~/projects/AS/analysis/novel_splice_site/SUPPA/"
anno_rmats <- glue("{dir_out}/rMATS/")
anno_miso <- "~/projects/AS/analysis/novel_splice_site/MISO/index/commonshortest/"
# anno_vast <- "~/bin/vast-tools/VASTDB/Hs2/TEMPLATES/"
prefix_suppa <- "events"
prefix_miso <- "hg38"
# prefix_vast <- "Hs2"
fil <- F

suppa <- parseSuppaAnnotation(anno_suppa, genome = prefix_suppa)
mats <- parseMatsAnnotation(anno_rmats, novelEvents = novel)
miso <- parseMisoAnnotation(anno_miso, genome = prefix_miso, types = c("SE", "MXE", "A5SS", "A3SS", "RI"))
# vast <- parseVastToolsAnnotation(anno_vast, genome = prefix_vast)
annot <- prepareAnnotationFromEvents(suppa, mats, miso) # , vast
saveRDS(annot, file = glue("{dir_out}/psichomics/annot.rds"))

prepareGeneQuant(
  files_genequant[sampleTable$samples], strandedness = sn,
  output = glue("{dir_out}/psichomics/psichomics_gene_counts.txt")
)
prepareJunctionQuant(
  files_junctionquant, output = glue("{dir_out}/psichomics/psichomics_junctions.txt")
)
write.table(
  dplyr::select(sampleTable, samples, conditions) %>% dplyr::rename(`Sample ID` = samples),
  glue("{dir_out}/psichomics/psichomics_metadata.txt"),
  sep = "\t", quote = F, row.names = F
)

data <- loadLocalFiles(glue("{dir_out}/psichomics/"))
sampleInfo <- data[[1]]$`Sample metadata`
geneExpr <- data[[1]]$`Gene expression`
junctionQuant <- data[[1]]$`Junction quantification`
if (fil) {
  filter <- filterGeneExpr(geneExpr)
} else {
  filter <- NULL
}
geneExprNorm <- normaliseGeneExpression(geneExpr, geneFilter = filter)
psi <- quantifySplicing(annot, junctionQuant)

condition_1 <- unique(sampleTable$conditions)[1]
condition_2 <- unique(sampleTable$conditions)[2]
s1VSs2 <- list(
  sampleTable$samples[sampleTable$conditions == condition_1],
  sampleTable$samples[sampleTable$conditions == condition_2]
) %>% setNames(unique(sampleTable$conditions))
diffSplicing <- diffAnalyses(psi, s1VSs2)
write.table(
  diffSplicing, glue("{dir_out}/psichomics/results.txt"),
  quote = F, sep = "\t"
)
write.table(
  psi, glue("{dir_out}/psichomics/psi.txt"),
  quote = F, sep = "\t"
)
save.image(glue("{dir_out}/psichomics/psichomics.RData"), compress = T)


##### JuncBASE #####
## 0.05 ##
myCreatedir(glue("{dir_out}/JuncBASE/cufflinks/"))

dir_juncbase <- "~/bin/juncBASE-1.2-beta/"
db <- "juncbase2"

write.table(
  dplyr::select(sampleTable, samples, files_bam), glue("{dir_out}/JuncBASE/bams.txt"),
  sep = "\t", col.names = F, row.names = F, quote = F
)

library_type <- c()
for (i in 1:nrow(sampleTable)) {
  if (sampleTable$strandedness[i] == "no") library_type[i] <- "fr-unstranded"
  if (sampleTable$strandedness[i] == "yes") library_type[i] <- "fr-secondstrand"
  if (sampleTable$strandedness[i] == "reverse") library_type[i] <- "fr-firststrand"
}

cmds <- list()
if (novel) {
  cmds[["building a de-novo transcript database (1): cufflinks"]] <- glue(
    "cufflinks \\
    -p {nproc} \\
    -o {dir_out}/JuncBASE/cufflinks/{sampleTable$samples} \\
    -g {gtf} \\
    -u \\
    --library-type {library_type} \\
    {sampleTable$files_bam}"
  )
  cmds[["building a de-novo transcript database (2): cuffmerge"]] <- glue(
    "find {dir_out}/JuncBASE/cufflinks/ | grep transcripts.gtf > {dir_out}/JuncBASE/cufflinks/cufflinks_denovo.txt && \\
    mkdir -p {dir_out}/JuncBASE/tmp/ && \\
    python {dir_juncbase}/other_scripts/runCuffmerge.py \\
    -f {dir_out}/JuncBASE/cufflinks/cufflinks_denovo.txt \\
    --no_single_exons \\
    -g {gtf} \\
    -o {dir_out}/JuncBASE/cuffmerge/ \\
    -s {fa} \\
    --tmp_dir {dir_out}/JuncBASE/tmp/"
  )
  cmds[["building a de-novo transcript database (3): build DB"]] <- glue(
    "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/build_DB_FromGTF.py \\
    -g {dir_out}/JuncBASE/cuffmerge/merged.gtf \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/ \\
    --initialize && \\
    ~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/build_DB_FromGTF.py \\
    -g {dir_out}/JuncBASE/cuffmerge/merged.gtf \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/"
  )
} else {
  cmds[["building annotated AS DB in mysql"]] <- glue(
    "python {dir_juncbase}/build_DB_FromGTF.py \\
    -g {gtf} \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/ \\
    --initialize && \\
    python {dir_juncbase}/build_DB_FromGTF.py \\
    -g {gtf} \\
    -d {db} \\
    -f {fa} \\
    --sqlite_db_dir {dir_out}/JuncBASE/"
  )
}

cmds[["processing bam files"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/run_preProcess_by_chr_step1.py \\
  -i {dir_out}/JuncBASE/bams.txt \\
  -o {dir_out}/JuncBASE/input_files/ \\
  -p {nrow(sampleTable)} \\
  --preProcess_options \"\" \\
  --nice"
)

cmds[["disambiguating splice junction orietations"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/disambiguate_junctions.py \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -g {fa} \\
  --by_chr"
)

cmds[["identifying all junctions"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/preProcess_getASEventReadCounts_by_chr_step2.py \\
  -i {dir_out}/JuncBASE/input_files/ \\
  --by_chr"
)

cmds[["creating exon-intron junction count files"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/run_preProcess_step3_by_chr.py \\
  --input_dir {dir_out}/JuncBASE/input_files/ \\
  --num_processes {nrow(sampleTable)}"
)

cmds[["creating a pseudo all junction sample"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/createPseudoSample.py \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -s {sample(sampleTable$samples, 1)} \\
  --by_chr"
)

cmds[["identifying alternative splicing events and quantify events from each sample"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/run_getASEventReadCounts_multiSample.py \\
  -s {paste(sampleTable$samples, collapse = ',')} \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -o {dir_out}/JuncBASE/aseventreadcounts/ \\
  --txt_db1 {db} \\
  --txt_db2 {db} \\
  --txt_db3 {db} \\
  --sqlite_db_dir {dir_out}/JuncBASE/ \\
  --jcn_seq_len {(basics$read_length_most - 6) * 2} \\
  -p {nrow(sampleTable)} \\
  --nice \\
  --by_chr"
)

myCreatedir(glue("{dir_out}/JuncBASE/tables/"))

cmds[["creating tables of raw and length-normalized read counts of exclusion and inclusion isoforms"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/run_createAS_CountTables.py \\
  -d {dir_out}/JuncBASE/aseventreadcounts/ \\
  -i {dir_out}/JuncBASE/input_files/ \\
  -s {paste(sampleTable$samples, collapse = ',')} \\
  --num_processes {nrow(sampleTable)} \\
  --jcn_seq_len {(basics$read_length_most - 6) * 2} && \\
  ~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/combine_createAS_CountTables_by_chr.py \\
  -d {dir_out}/JuncBASE/aseventreadcounts/ \\
  -o {dir_out}/JuncBASE/tables/juncbase"
)

cmds[["differential splicing between two groups of samples"]] <- glue(
  "~/bin/anaconda3/envs/assw/bin/python2.7 {dir_juncbase}/compareSampleSets.py \\
  --in_prefix {dir_out}/JuncBASE/tables/juncbase \\
  --all_psi_output {dir_out}/JuncBASE/sample_set_comparison.txt \\
  --mt_correction BH \\
  --thresh 10 \\
  --delta_thresh 0.0 \\
  --sample_set1 {paste(basics$samples_1, collapse = ',')} \\
  --sample_set2 {paste(basics$samples_2, collapse = ',')} \\
  --html_dir {dir_out}/JuncBASE/html/ \\
  --pdf \\
  --html_out_sign_thresh 1.0"
)

