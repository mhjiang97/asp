##### SANJUAN #####
myCreatedir(glue("{dir_out}/SANJUAN/"))

add_genome <- T
name_genome <- "v32_pa"
dir_sanjuan <- "~/bin/SANJUAN/"

library_type <- table(sampleTable$library_types) %>% .[. == max(.)] %>% names()
strandedness <- table(sampleTable$strandedness) %>% .[. == max(.)] %>% names()
l <- 1
if (library_type == "paired-end") l <- 2
if (strandedness %in% c("yes", "reverse")) {
  l <- paste0(l, "S")
} else {l <- paste0(l, "U")}

if (add_genome) {
  file.symlink(gtf, glue("{dir_sanjuan}/db/gtfs/{name_genome}.gtf"))
  cmds <- list()
  cmds[["building sanjuan annotation files"]] <- glue(
    "build_annotations_from_gtf.pl \\
    {dir_sanjuan}/db/gtfs/{name_genome}.gtf \\
    {name_genome} \\
    {dir_sanjuan}"
  )
  cmds[["samtools indexing fa"]] <- glue(
    "samtools faidx {fa}"
  )
  cmds[["generating genmoe file"]] <- glue(
    "awk -v OFS='\\t' '{print 'chr'$1,$2}' [fa].fai > [dir_sanjuan]/db/genomes/[name_genome].genome",
    .open = "[",
    .close = "]",
  )
  run_cmds(cmds)
}

cmds <- list()
cmds[["running sanjuan"]] <- glue(
  "sanjuan.pl \\
  -nprocs {nproc} \\
  -g {name_genome} \\
  -b S \\
  -g1 {basics$condition_1} \\
  -f1 {paste(basics$bams_1, collapse = ',')} \\
  -g2 {basics$condition_2} \\
  -f2 {paste(basics$bams_2, collapse = ',')} \\
  -o {dir_out}/SANJUAN/ \\
  -p phred33 \\
  -l {l} \\
  -c MC"
)
run_cmds(cmds)


##### JUM #####
myCreatedir(glue("{dir_out}/JUM/JUM_diff/"))

dir_jum <- "~/bin/JUM-JUMv2.0.2_short_user_argument_BSD/"
bam2sam <- F
JuncThreshold <- 5
IRthreshold <- 5
Condition1_fileNum_threshold <- 3
Condition2_fileNum_threshold <- 3
Cutoff <- 0.05
REF <- "~/doc/reference/genepred/gencode_v32.primary_assembly.protein_coding2.gpd"
dir_index_star <- "~/doc/reference/star_2.7.5a_pa"

Readlength <- table(sampleTable$read_lengths) %>% .[. == max(.)] %>% names() %>% as.numeric()
condition_1 <- unique(sampleTable$conditions)[1]
condition_2 <- unique(sampleTable$conditions)[2]
samples_1 <- sampleTable$samples[sampleTable$conditions == condition_1]
samples_2 <- sampleTable$samples[sampleTable$conditions == condition_2]
files_fq <- sapply(strsplit(sampleTable$files_fq, ';'), function(x) {paste(x, collapse = " ")})

cmds <- list()
cmds[["star mapping 1st pass"]] <- glue(
  "STAR --runThreadN {nproc} --genomeDir {dir_index_star} \\
  --outFileNamePrefix {dir_out}/JUM/{sampleTable$samples} \\
  --readFilesIn {files_fq} \\
  --readFilesCommand zcat \\
  --outSJfilterReads Unique 1>{dir_out}/JUM/_log.{sampleTable$samples}.star1 2>&1"
)
if (parallel) cmds[["star mapping 1st pass"]] <- myParallel2(cmds[["star mapping 1st pass"]])
run_cmds(cmds)
dir.create(glue("{dir_out}/JUM/1st_SJ/"))
system(glue("mv {dir_out}/JUM/*SJ.out.tab {dir_out}/JUM/*Log.final.out {dir_out}/JUM/*Log.progress.out {dir_out}/JUM/*Log.out {dir_out}/JUM/1st_SJ/"))
system(glue("rm {dir_out}/JUM/*Aligned.out.sam"))
files_sj <- list.files(glue("{dir_out}/JUM/1st_SJ/"), pattern = "SJ.out.tab", full.names = T) %>%
  paste(collapse = " ")
cmds <- list()
cmds[["star mapping 2nd pass"]] <- glue(
  "STAR --runThreadN {nproc} --genomeDir {dir_index_star} \\
  --outFileNamePrefix {dir_out}/JUM/{sampleTable$samples} \\
  --readFilesIn {files_fq} \\
  --readFilesCommand zcat \\
  --outSJfilterReads Unique \\
  --outSAMstrandField intronMotif \\
  --outFilterMultimapNmax 1 \\
  --sjdbFileChrStartEnd {files_sj} 1>{dir_out}/JUM/_log.{sampleTable$samples}.star2 2>&1"
)
cmds[["star mapping 2nd pass"]] <- myParallel2(cmds[["star mapping 2nd pass"]])
run_cmds(cmds)
cmds <- list()
cmds[["converting sams to bams"]] <- glue(
  "samtools view -@ {nproc} -bS {dir_out}/JUM/{sampleTable$samples}Aligned.out.sam >{dir_out}/JUM/{sampleTable$samples}Aligned.out.bam"
) %>% myParallel2()
cmds[["sorting bams"]] <- glue(
  "samtools sort -@ {nproc} -o {dir_out}/JUM/{sampleTable$samples}Aligned.out_sorted.bam -T {dir_out}/JUM/{sampleTable$samples}_temp {dir_out}/JUM/{sampleTable$samples}Aligned.out.bam"
) %>% myParallel2()
cmds[["indexing bams"]] <- glue(
  "samtools index -@ {nproc} {dir_out}/JUM/{sampleTable$samples}Aligned.out_sorted.bam"
) %>% myParallel2()
system(glue("rm {dir_out}/JUM/*Aligned.out.bam"))

if (bam2sam) {
  cmds <- list()
  cmds[["converting bams to sams"]] <- glue(
    'samtools view -h -@ {nproc} {sampleTable$files_bam} -O SAM \\
    -o {gsub(".bam",".sam",sampleTable$files_bam)}'
  )
  if (parallel) {
    cmds[["converting bams to sams"]] <- myParallel2(cmds[["converting bams to sams"]], np)
  }
  run_cmds(cmds, T)
}

for (s in sampleTable$samples) {
  file_bam <- sampleTable$files_bam[sampleTable$samples == s]
  file.symlink(file_bam, glue("{dir_out}/JUM/{s}Aligned.out_sorted.bam"))
  file.symlink(paste0(file_bam, ".bai"), glue("{dir_out}/JUM/{s}Aligned.out_sorted.bam.bai"))

  file_sam <- sampleTable$files_bam[sampleTable$samples == s] %>% gsub("\\.bam", "\\.sam", .)
  file.symlink(file_sam, glue("{dir_out}/JUM/{s}Aligned.out.sam"))

  file_junc <- glue("{dirname(file_bam)}/SJ.out.tab")
  file.symlink(file_junc, glue("{dir_out}/JUM/{s}SJ.out.tab"))
}

cmds <- list()
cmds[["step 1 running JUM_A"]] <- glue(
  "cd {dir_out}/JUM/ && bash {dir_jum}/JUM_A.sh \\
  --Folder {dir_jum} \\
  --JuncThreshold {JuncThreshold} \\
  --Condition1_fileNum_threshold {Condition1_fileNum_threshold} \\
  --Condition2_fileNum_threshold {Condition2_fileNum_threshold} \\
  --IRthreshold {IRthreshold} \\
  --Readlength {Readlength} \\
  --Thread {ceiling(nproc / nrow(sampleTable))} \\
  --Condition1SampleName {paste(samples_1, collapse = ',')} \\
  --Condition2SampleName {paste(samples_2, collapse = ',')} && \\
  rm -rf {dir_out}/JUM/temp_JUM_A_run/"
)
if (write_log) cmds[["running JUM_A"]] <- paste(cmds[["running JUM_A"]], glue("1>{dir_out}/JUM/_log.jum_a 2>&1"))
write.table(
  dplyr::select(sampleTable, samples, conditions) %>% dplyr::rename(condition = conditions),
  glue("{dir_out}/JUM/JUM_diff/experiment_design.txt"),
  quote = F, row.names = F, sep = "\t"
)
system(glue('sed -i "1s/samples//g" {dir_out}/JUM/JUM_diff/experiment_design.txt'))
cmds[["step 2 running R_script_JUM"]] <- glue(
  "cd {dir_out}/JUM/JUM_diff/ && Rscript {dir_jum}/R_script_JUM.R \\
  experiment_design.txt > outputFile.Rout 2> errorFile.Rout"
)
cmds[["step 3 running JUM_B"]] <- glue(
  "cd {dir_out}/JUM/JUM_diff/ && bash {dir_jum}/JUM_B.sh \\
  --Folder {dir_jum} \\
  --Test adjusted_pvalue \\
  --Cutoff {Cutoff} \\
  --TotalFileNum {nrow(sampleTable)} \\
  --Condition1_fileNum_threshold {Condition1_fileNum_threshold} \\
  --Condition2_fileNum_threshold {Condition2_fileNum_threshold} \\
  --Condition1SampleName {condition_1} \\
  --Condition2SampleName {condition_2}"
)
if (write_log) cmds[["step 3 running JUM_B"]] <- paste(cmds[["step 3 running JUM_B"]], glue("{dir_out}/JUM/JUM_diff/_log.jum_b 2>&1"))
cmds[["step4 running JUM_C"]] <- glue(
  "cd {dir_out}/JUM/JUM_diff/FINAL_JUM_OUTPUT_adjusted_pvalue_{Cutoff}/ && \\
  bash {dir_jum}/JUM_C.sh \\
  --Folder {dir_jum} \\
  --Test adjusted_pvalue \\
  --TotalConditoin1FileNum {length(samples_1)} \\
  --TotalConditoin2FileNum {length(samples_2)} \\
  --REF {REF}"
)
if (write_log) cmds[["step4 running JUM_C"]] <- paste(cmds[["step4 running JUM_C"]], glue("{dir_out}/JUM/JUM_diff/FINAL_JUM_OUTPUT_adjusted_pvalue_{Cutoff}/_log.jum_c 2>&1"))


##### jSplice #####
myCreatedir(glue("{dir_out}/jSplice/bed/"))

dir_jsplice <- "~/bin/jsplice/bin/"

cmds <- list()
cmds[["juncs to beds"]] <- glue(
  "{dir_jsplice}/starJxn2bed \\
  -f {c(basics$juncs_star_1, basics$juncs_star_2)} \\
  -o {dir_out}/jSplice/bed/{c(basics$samples_1, basics$samples_2)}.bed"
)

tbl_jsplice <- dplyr::select(sampleTable, conditions, files_bam)
tbl_jsplice$files_bed <- glue("{dir_out}/jSplice/bed/{c(basics$samples_1, basics$samples_2)}.bed")
tbl_jsplice$experiment <- "experiemnt"
tbl_jsplice <- dplyr::select(tbl_jsplice, experiment, conditions, files_bed, files_bam)
tbl_jsplice$files_bed <- path.expand(tbl_jsplice$files_bed)
tbl_jsplice$files_bam <- path.expand(tbl_jsplice$files_bam)
write.table(
  tbl_jsplice, glue("{dir_out}/jSplice/expdesign.txt"),
  sep = "\t", quote = F, row.names = F, col.names = F
)

cmds[["running jsplice"]] <- glue(
  "{dir_jsplice}/jsplice \\
  -o {dir_out}/jSplice/ \\
  -a {gtf} \\
  -d {dir_out}/jSplice/expdesign.txt \\
  -k 0 \\
  -n {ceiling(nproc/nrow(sampleTable))}"
)


##### Spanki #####
myCreatedir(glue("{dir_out}/Spanki/junctions_out/"))

astalavista <- T
file_astalavista <- NULL

cmds <- list()
cmds[["running astalavista"]] <- glue(
  "astalavista -t asta -i {gtf} -e [ASE,ASI] -o {dir_out}/Spanki/{basename(gtf)}_astalavista.gtf.gz \\
  && gunzip {dir_out}/Spanki/{basename(gtf)}_astalavista.gtf.gz"
)
cmds[["evaluating junction alignment"]] <- glue(
  "spankijunc -i {sampleTable$files_bam} -g {gtf} -f {fa} -o {dir_out}/Spanki/junctions_out/"
)

