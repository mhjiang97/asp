pkgs <- c("glue", "dplyr", "org.Hs.eg.db", "vroom")
for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = T))

library(asp)
##### dux4 #####
bams_dux4 <- glue(
  "~/projects/dux4/analysis/mgy/rnaseq/star/DUX4_{c(1, 2, 3)}/DUX4_{c(1, 2, 3)}.sortedByCoord.bam"
)
bams_vec <- glue(
  "~/projects/dux4/analysis/mgy/rnaseq/star/VEC_{c(1, 2, 3)}/VEC_{c(1, 2, 3)}.sortedByCoord.bam"
)
fqs_dux4 <- glue(
  "~/projects/dux4/data/mgy/rnaseq/DUX4_{c(1, 2, 3)}.fq.gz"
)
fqs_vec <- glue(
  "~/projects/dux4/data/mgy/rnaseq/VEC_{c(1, 2, 3)}.fq.gz"
)
dirs_salmon_dux4 <- glue(
  "~/projects/dux4/analysis/mgy/rnaseq/salmon/DUX4_{c(1, 2, 3)}"
)
dirs_salmon_vec <- glue(
  "~/projects/dux4/analysis/mgy/rnaseq/salmon/VEC_{c(1, 2, 3)}"
)

sampleTable <- data.frame(
  samples = paste0(rep(c("DUX4", "VEC"), each = 3), "_", rep(c(1:3), 2)),
  conditions = rep(c("dux4", "vec"), each = 3),
  files_bam = c(bams_dux4, bams_vec),
  files_fq = c(fqs_dux4, fqs_vec),
  dirs_salmon = c(dirs_salmon_dux4, dirs_salmon_vec),
  read_lengths = rep(50, 6),
  library_types = rep("single-end", 6),
  strandedness = rep("no", 6)
)

asp <- ASP(
  sampletable = sampleTable,
  gtf = "~/doc/reference/gtf/gencode.v32.annotation.gtf",
  gff = "~/doc/reference/gtf/gff/gencode.v32.annotation.gff3",
  fa = "~/doc/reference/fa/GRCh38.p13.genome.fa",
  fa_transcript = "~/doc/reference/fa/gencode.v32.transcripts.fa",
  genome_version = "hg38",
  dir_out = "~/projects/dux4/analysis/mgy/rnaseq/as/",
  nproc = 20,
  novel = T,
  write_log = T,
  parallel = T,
  np = 5,
  conda_path = "~/miniconda3/"
)

asp <- cmd_asp(
  asp = asp, conda_env = "python3", tool = "rmats",
  variable_read_length = T, paired_stats = F, cstat = 0.01, allow_clipping = T
)
asp <- run_asp(asp = asp, tool = "rmats", run = T, read_results = T)
asp <- sig_asp(asp = asp, tool = "rmats", p_adj = 0.05, delta_psi = 0.1)

asp <- cmd_asp(
  asp = asp, tool = "cash",
  cash_jar = "~/projects/as/software/cash_v2.2.0/cash.jar", MergePval = "G", ram = "100g"
)
asp <- run_asp(asp = asp, tool = "cash", run = T, read_results = T)
asp <- sig_asp(asp = asp, tool = "cash", p_value = 0.1)

asp <- cmd_asp(
  asp = asp, conda_path = "~/env", conda_env = "", tool = "majiq",
  show_all = F, threshold = 0.1
)
asp <- run_asp(asp = asp, tool = "majiq", run = T, read_results = T)
asp <- sig_asp(asp = asp, tool = "majiq", p_value = 0.05)

asp <- cmd_asp(
  asp = asp, tool = "leafcutter", conda_env = "assw_3.8",
  exon_file = "~/doc/leafcutter/gencode_v32_exons.txt.gz",
  leafcutter_dir = "~/projects/as/software/leafcutter"
)
asp <- cmd_asp(
  asp = asp, tool = "bandits", conda_env = "r4.1",
  tx2gene = "~/doc/tx2gene/tx2gene_v32.txt"
)


##### new features #####
asp <- cmd_asp(
  asp = asp, conda_env = c("python3", NA, ""), conda_path = c(majiq = "~/env/"),
  tool = c("rmats", "cash", "majiq"),
  variable_read_length = T, paired_stats = F, cstat = 0.01, allow_clipping = T,
  cash_jar = "~/projects/as/software/cash_v2.2.0/cash.jar", MergePval = "G", ram = "100g",
  show_all = F, threshold = 0.1
)
asp <- run_asp(asp = asp, tool = c("rmats", "majiq", "cash"), run = T, read_results = T, parallel = T)


##### tidy #####
bams_s1 <- glue(
  "~/projects/AS/data/sim/rsem/paper_gb/sim_1/star/sample{1:4}/sample{1:4}.SortedByCoord.bam"
)
bams_s2 <- glue(
  "~/projects/AS/data/sim/rsem/paper_gb/sim_1/star/sample{5:8}/sample{5:8}.SortedByCoord.bam"
)
fqs_s1 <- paste(
  glue("~/projects/AS/data/sim/rsem/paper_gb/sim_1/fq/sample{1:4}_R1.fq.gz"),
  glue("~/projects/AS/data/sim/rsem/paper_gb/sim_1/fq/sample{1:4}_R2.fq.gz"),
  sep = ";"
)
fqs_s2 <- paste(
  glue("~/projects/AS/data/sim/rsem/paper_gb/sim_1/fq/sample{5:8}_R1.fq.gz"),
  glue("~/projects/AS/data/sim/rsem/paper_gb/sim_1/fq/sample{5:8}_R2.fq.gz"),
  sep = ";"
)
dirs_salmon_s1 <- glue(
  "~/projects/AS/data/sim/rsem/paper_gb/sim_1/salmon/sample{1:4}"
)
dirs_salmon_s2 <- glue(
  "~/projects/AS/data/sim/rsem/paper_gb/sim_1/salmon/sample{5:8}"
)

sampleTable <- data.frame(
  samples = glue("sample{1:8}"),
  conditions = rep(c("s1", "s2"), each = 4),
  files_bam = c(bams_s1, bams_s2),
  files_fq = c(fqs_s1, fqs_s2),
  dirs_salmon = c(dirs_salmon_s1, dirs_salmon_s2),
  read_lengths = rep(101, 8),
  library_types = rep("paired-end", 8),
  strandedness = rep("no", 8)
)

tx2gene <- vroom(
  "~/doc/tx2gene/tx2gene_v32.txt", delim = "\t", col_types = cols(),
  col_names = c("chr", "gene_id", "gene_symbol", "transcript_id")
)

asp <- createASP(
  sampletable = sampleTable,
  tx2gene = tx2gene,
  gtf = "~/doc/reference/gtf/gencode.v32.annotation.gtf",
  gff = "~/doc/reference/gtf/gff/gencode.v32.annotation.gff3",
  fa = "~/doc/reference/fa/GRCh38.p13.genome.fa",
  fa_transcript = "~/doc/reference/fa/gencode.v32.transcripts.fa",
  genome_version = "hg38",
  dir_out = "~/projects/AS/analysis/paper_gb/sim_1/",
  nproc = 20,
  novel = F,
  write_log = T,
  parallel = T,
  np = 5,
  conda_path = "~/bin/anaconda3/"
)

asp <- run_asp(asp, tool = "rmats", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "cash", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "spladder", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "leafcutter", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "suppa", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "psichomics", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "majiq", run = F, read_results = T, parallel = F)
asp <- run_asp(asp, tool = "bandits", run = F, read_results = T, parallel = F)

asp <- sig_asp(
  asp, tool = c("rmats", "cash", "leafcutter", "suppa", "psichomics"),
  p_adj = 0.05, delta_psi = 0.05
)
asp <- sig_asp(
  asp, tool = "majiq", p_value = 0.05, delta_psi = 0.05
)
asp <- sig_asp(
  asp, tool = "spladder", p_adj = 0.05, effect_size_col = "log2FC_event_count", effect_size_threshod = 1
)
asp <- sig_asp(
  asp, tool = "bandits", p_adj = 0.05, effect_size_col = "DTU_measure", effect_size_threshod = 0.1
)

asp <- intersect_asp(asp)

asp <- plot_pre_asp(asp, org.Hs.eg.db)

asp <- plot_asp(asp, tool = "rmats", plot_samples = 4)
asp <- plot_asp(asp, tool = "cash", plot_samples = 4)
asp <- plot_asp(asp, tool = "spladder", plot_samples = 4)
asp <- plot_asp(asp, tool = "leafcutter", plot_samples = c("sample1", "sample2", "sample5", "sample6"))
asp <- plot_asp(asp, tool = "suppa", plot_samples = c("sample1", "sample5", "sample6"))
asp <- plot_asp(asp, tool = "psichomics", plot_samples = c("sample1"))
asp <- plot_asp(asp, tool = "majiq", plot_samples = 1)





