pkgs <- c("dplyr", "vroom", "tidyr", "tibble")
for (pkg in pkgs) suppressPackageStartupMessages(library(pkg, character.only = T))
setwd("inst/extdata/results/")
tx2gene <- read.delim("~/doc/tx2gene/tx2gene_v32.txt.gz", header = F)
genes_keep <- c("NBPF26", "HAT1", "INTS13", "NCOA3", "PIANP", "DROSHA", "PIGN", "NASP", "RPL35A", "TMUB2")
genes_id_keep <- unique(tx2gene[match(genes_keep, tx2gene[, 3]), 2])

##### BANDITS #####
dat <- read.delim("BANDITS/results_gene.txt")
dat2 <- read.delim("BANDITS/results_transcript.txt")

dat_keep <- dat[dat$Gene_id %in% genes_id_keep, ]
dat2_keep <- dat2[dat2$Gene_id %in% genes_id_keep, ]

write.table(
  dat_keep, "BANDITS/results_gene_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat2_keep, "BANDITS/results_transcript_keep.txt",
  sep = "\t", quote = F, row.names = F
)

##### CASH #####
dat <- read.delim("CASH/cash_G.s1vss2.alldiff.txt")

dat_keep <- dat[dat$AccID %in% genes_keep, ]

write.table(
  dat_keep, "CASH/cash_G.s1vss2.alldiff_keep.txt",
  sep = "\t", quote = F, row.names = F
)

##### LeafCutter #####
dat <- read.delim("LeafCutter/ds/leafcutter_cluster_significance.txt")
dat2 <- read.delim("LeafCutter/ds/leafcutter_effect_sizes.txt")
dat2 <- dat2 |>
  rowwise() |>
  mutate(tmp = rev(unlist(strsplit(intron, ":")))[1]) |>
  ungroup()

dat_keep <- dat[dat$genes %in% genes_keep, ]
dat_keep <- dat_keep |>
  rowwise() |>
  mutate(tmp = rev(unlist(strsplit(cluster, ":")))[1]) |>
  ungroup()
dat2_keep <- dat2[dat2$tmp %in% dat_keep$tmp, ]
dat_keep <- dat_keep |> select(!tmp)
dat2_keep <- dat2_keep |> select(!tmp)

write.table(
  dat_keep, "LeafCutter/ds/leafcutter_cluster_significance_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat2_keep, "LeafCutter/ds/leafcutter_effect_sizes_keep.txt",
  sep = "\t", quote = F, row.names = F
)

##### MAJIQ #####
dat <- read.delim("MAJIQ/voila/s1_s2_0.2.tsv", comment.char = "#")

dat_keep <- dat[dat$gene_name %in% genes_keep, ]

write.table(
  dat_keep, "MAJIQ/voila/s1_s2_0.2_keep.tsv",
  sep = "\t", quote = F, row.names = F
)

##### psichomics #####
dat <- vroom("psichomics/results.txt", delim = "\t")

dat_keep <- dat[dat$Gene %in% genes_keep, ]

write.table(
  dat_keep, "psichomics/results_keep.txt",
  sep = "\t", quote = F, row.names = F
)

##### rMATS #####
dat <- read.delim("rMATS/RI.MATS.JC.txt.gz")
dat2 <- read.delim("rMATS/A5SS.MATS.JC.txt.gz")
dat3 <- read.delim("rMATS/A3SS.MATS.JC.txt.gz")
dat4 <- read.delim("rMATS/MXE.MATS.JC.txt.gz")
dat5 <- read.delim("rMATS/SE.MATS.JC.txt.gz")

dat_keep <- dat[dat$geneSymbol %in% genes_keep, ]
dat2_keep <- dat2[dat2$geneSymbol %in% genes_keep, ]
dat3_keep <- dat3[dat3$geneSymbol %in% genes_keep, ]
dat4_keep <- dat4[dat4$geneSymbol %in% genes_keep, ]
dat5_keep <- dat5[dat5$geneSymbol %in% genes_keep, ]

write.table(
  dat_keep, "rMATS/RI.MATS.JC_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat2_keep, "rMATS/A5SS.MATS.JC_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat3_keep, "rMATS/A3SS.MATS.JC_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat4_keep, "rMATS/MXE.MATS.JC_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat5_keep, "rMATS/SE.MATS.JC_keep.txt",
  sep = "\t", quote = F, row.names = F
)

##### SplAdder #####
dat <- read.delim("SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_mult_exon_skip.tsv.gz")
dat2 <- read.delim("SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_mutex_exons.tsv.gz")
dat3 <- read.delim("SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_intron_retention.tsv.gz")
dat4 <- read.delim("SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_exon_skip.tsv.gz")
dat5 <- read.delim("SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_alt_5prime.tsv.gz")
dat6 <- read.delim("SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_alt_3prime.tsv.gz")

dat7 <- read.delim("SplAdder/graphs/merge_graphs_mult_exon_skip_C3.confirmed.txt.gz")
dat8 <- read.delim("SplAdder/graphs/merge_graphs_mutex_exons_C3.confirmed.txt.gz")
dat9 <- read.delim("SplAdder/graphs/merge_graphs_intron_retention_C3.confirmed.txt.gz")
dat10 <- read.delim("SplAdder/graphs/merge_graphs_exon_skip_C3.confirmed.txt.gz")
dat11 <- read.delim("SplAdder/graphs/merge_graphs_alt_5prime_C3.confirmed.txt.gz")
dat12 <- read.delim("SplAdder/graphs/merge_graphs_alt_3prime_C3.confirmed.txt.gz")

dat_keep <- dat[dat$gene %in% genes_id_keep, ]
dat2_keep <- dat2[dat2$gene %in% genes_id_keep, ]
dat3_keep <- dat3[dat3$gene %in% genes_id_keep, ]
dat4_keep <- dat4[dat4$gene %in% genes_id_keep, ]
dat5_keep <- dat5[dat5$gene %in% genes_id_keep, ]
dat6_keep <- dat6[dat6$gene %in% genes_id_keep, ]

dat7_keep <- dat7[dat7$gene_name %in% genes_id_keep, ]
dat8_keep <- dat8[dat8$gene_name %in% genes_id_keep, ]
dat9_keep <- dat9[dat9$gene_name %in% genes_id_keep, ]
dat10_keep <- dat10[dat10$gene_name %in% genes_id_keep, ]
dat11_keep <- dat11[dat11$gene_name %in% genes_id_keep, ]
dat12_keep <- dat12[dat12$gene_name %in% genes_id_keep, ]

write.table(
  dat_keep, "SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_mult_exon_skip_keep.tsv",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat2_keep, "SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_mutex_exons_keep.tsv",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat3_keep, "SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_intron_retention_keep.tsv",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat4_keep, "SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_exon_skip_keep.tsv",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat5_keep, "SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_alt_5prime_keep.tsv",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat6_keep, "SplAdder/graphs/testing_s1_vs_s2/test_results_extended_C3_alt_3prime_keep.tsv",
  sep = "\t", quote = F, row.names = F
)

write.table(
  dat7_keep, "SplAdder/graphs/merge_graphs_mult_exon_skip_C3.confirmed_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat8_keep, "SplAdder/graphs/merge_graphs_mutex_exons_C3.confirmed_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat9_keep, "SplAdder/graphs/merge_graphs_intron_retention_C3.confirmed_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat10_keep, "SplAdder/graphs/merge_graphs_exon_skip_C3.confirmed_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat11_keep, "SplAdder/graphs/merge_graphs_alt_5prime_C3.confirmed_keep.txt",
  sep = "\t", quote = F, row.names = F
)
write.table(
  dat12_keep, "SplAdder/graphs/merge_graphs_alt_3prime_C3.confirmed_keep.txt",
  sep = "\t", quote = F, row.names = F
)

##### SUPPA #####
dat <- read.delim("SUPPA/ds/s1_s2.A3.dpsi.gz")
dat2 <- read.delim("SUPPA/ds/s1_s2.A5.dpsi.gz")
dat3 <- read.delim("SUPPA/ds/s1_s2.AL.dpsi.gz")
dat4 <- read.delim("SUPPA/ds/s1_s2.AF.dpsi.gz")
dat5 <- read.delim("SUPPA/ds/s1_s2.SE.dpsi.gz")
dat6 <- read.delim("SUPPA/ds/s1_s2.RI.dpsi.gz")
dat7 <- read.delim("SUPPA/ds/s1_s2.MX.dpsi.gz")

dat <- dat |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")
dat2 <- dat2 |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")
dat3 <- dat3 |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")
dat4 <- dat4 |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")
dat5 <- dat5 |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")
dat6 <- dat6 |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")
dat7 <- dat7 |> rownames_to_column("event") |> separate(event, c("gene_id", "coord"), sep = ";")

dat_keep <- dat[dat$gene_id %in% genes_id_keep, ]
dat2_keep <- dat2[dat2$gene_id %in% genes_id_keep, ]
dat3_keep <- dat3[dat3$gene_id %in% genes_id_keep, ]
dat4_keep <- dat4[dat4$gene_id %in% genes_id_keep, ]
dat5_keep <- dat5[dat5$gene_id %in% genes_id_keep, ]
dat6_keep <- dat6[dat6$gene_id %in% genes_id_keep, ]
dat7_keep <- dat7[dat7$gene_id %in% genes_id_keep, ]

dat_keep <- dat_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())
dat2_keep <- dat2_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())
dat3_keep <- dat3_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())
dat4_keep <- dat4_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())
dat5_keep <- dat5_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())
dat6_keep <- dat6_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())
dat7_keep <- dat7_keep |> rowwise() |> mutate(ID = paste(gene_id, coord, sep = ";")) |> ungroup() |> select(!c(gene_id, coord)) |> select(ID, everything())

write.table(dat_keep, "SUPPA/ds/s1_s2.A3_keep.dpsi", sep = "\t", quote = F, row.names = F)
write.table(dat2_keep, "SUPPA/ds/s1_s2.A5_keep.dpsi", sep = "\t", quote = F, row.names = F)
write.table(dat3_keep, "SUPPA/ds/s1_s2.AL_keep.dpsi", sep = "\t", quote = F, row.names = F)
write.table(dat4_keep, "SUPPA/ds/s1_s2.AF_keep.dpsi", sep = "\t", quote = F, row.names = F)
write.table(dat5_keep, "SUPPA/ds/s1_s2.SE_keep.dpsi", sep = "\t", quote = F, row.names = F)
write.table(dat6_keep, "SUPPA/ds/s1_s2.RI_keep.dpsi", sep = "\t", quote = F, row.names = F)
write.table(dat7_keep, "SUPPA/ds/s1_s2.MX_keep.dpsi", sep = "\t", quote = F, row.names = F)




