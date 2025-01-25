# Genome QC:
- CheckM: Assess the genome quality and summarise basic information of genomes, including genome size, completeness, contamination, N50, coding density, contig number.
- barrnap: Predict rRNA.
- tRNAscan-SE: Predict tRNA.
   
**Quality evaluation**  
Put all the genomes into one directory `$wd/03_genome/`. Calculate the quality score (QS) according to the formula `QS = Completeness - 5 * Contamination`, only genomes with QS > 50 are retained for downstream analysis.
```bash
# run checkm
time checkm lineage_wf \
   $wd/03_genome/ $wd/04_QC/checkM/ \
   -f $DIR/all_genomes/02_checkM/results.tsv --tab_table \
   --pplacer_threads 24 -t 24 -x fna >> $wd/log/checkM.log 2>&1
# extract import information
sed 's/,/\t/g' $wd/04_QC/checkM/storage/bin_stats_ext.tsv | awk -F "\t" '{print $1 "\t" $12 "\t" $13 "\t" $14 "\t" $16 "\t" $19 "\t" $23 "\t" $26}' | sed "s/'Completeness': //g" | sed "s/'Contamination': //g" | sed "s/'GC': //g" | sed "s/'Genome size': //g" | sed "s/'# contigs': //g" | sed "s/'N50 (contigs)': //g" | sed "s/'Coding density': //g" | awk -F "\t" -v OFS="\t" '{$9=$2-5*$3;print $0}' | awk 'BEGIN{print "GenomeID" "\t" "Completeness" "\t" "Contamination" "\t" "GC" "\t" "Genome_size" "\t" "Contig_num" "\t" "N50" "\t" "Coding_density" "\t" "Quality_score"}1' > $wd/04_QC/checkM_result.tsv
```

**rRNA prediction**  
Extract the potential rRNA sequences and summarise.  
```bash
time parallel -j 2 --xapply \
    "barrnap \
        --kingdom bac --threads 64 \
        --outseq $wd/04_QC/rRNA_prediction/{1}_rRNA.fna \
        $wd/03_genome/{1}.fna >> $wd/log/rna_prediction.log 2>&1;
    grep '>' $wd/04_QC/rRNA_prediction/{1}_rRNA.fna | sed 's/>//;s/::/\t/' | cut -f 1 | sort | uniq -c | sed 's/^/{1}\t/;s/ 5S/\t5S/;s/ 16S/\t16S/;s/ 23S/\t23S/' | sed 's/ //g' >> $wd/04_QC/rRNA_summary.txt" \
:::$(cat GenomeID)
```  

**tRNA prediction**  
Predict the tRNA, summary the tRNA type and number.
```bash
time parallel -j 20 --xapply \
    "tRNAscan-SE \
        -B -I -o $wd/04_QC/tRNA_prediction/{1}_output.txt \
        $wd/03_genome/{1}.fna --thread 4;
    cut -f 5 $wd/04_QC/tRNA_prediction/{1}_output.txt | grep -v -E 'tRNA|Type|---' | sort | uniq -c | sed 's/^[[:space:]]*//' | sed 's/ /\t/;s/^/{1}\t/'>> $wd/04_QC/tRNA_summary.txt" \
:::$(cat GenomeID)
```
