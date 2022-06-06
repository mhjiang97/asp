cd ~/projects/AS/data/sim/rsem/paper_gb/sim_1/star/
nproc=0
for id in sample{1..8}
do
(samtools view \
-@ 20 \
-hb ${id}/${id}.SortedByCoord.bam \
-L ~/projects/pkgs/R/asp/inst/extdata/extract.bed \
> ~/projects/pkgs/R/asp/inst/extdata/${id}.SortedByCoord.bam && \
echo "${id} bam extraction was done" && \
samtools index ~/projects/pkgs/R/asp/inst/extdata/${id}.SortedByCoord.bam && \
samtools sort \
-@ 20 \
-O bam \
-T ${id} \
-n ~/projects/pkgs/R/asp/inst/extdata/${id}.SortedByCoord.bam \
> ~/projects/pkgs/R/asp/inst/extdata/${id}.SortedByName.bam && \
echo "${id} bam convertion was done" && \
bedtools bamtofastq \
-i ~/projects/pkgs/R/asp/inst/extdata/${id}.SortedByName.bam \
-fq >(gzip > ~/projects/pkgs/R/asp/inst/extdata/${id}_R1.fq.gz) \
-fq2 >(gzip > ~/projects/pkgs/R/asp/inst/extdata/${id}_R2.fq.gz) && \
echo "${id} fastq generation was done") \
1> ~/projects/pkgs/R/asp/inst/extdata/_log_${id} 2>&1 &

nproc=$((${nproc}+1))
if [ ${nproc} -ge 5 ]
then
wait
nproc=0
fi

done


cd ~/projects/pkgs/R/asp/inst/extdata/
nproc=0
for id in sample{1..8}
do
bash run.sh -c config.yaml -s ${id} 1>_log_rswp_${id} 2>&1 &

nproc=$((${nproc}+1))
if [ ${nproc} -ge 5 ]
then
wait
nproc=0
fi

done




