# Viral gene prediction and functional annotation:
- prodigal-gv was used to predict viral genes.
- VIBRANT was used to annotate
## Virus prediction:
**genomad**   
Use the `end-to-end` mode with the `--conservative` option to predict.   
The input files are fasta files of genomes or metagenome-assembled contig          
                        s.
```bas

# parallel run:
time parallel -j 3 --xapply \
    genomad end-to-end -t 32  \
        --conservative \
        --cleanup \
        --splits 32 $wd/03_genome/{1}.fna \
        $wd/12_viral_genome/genomad_output/{1} \
        --quiet \
        $db/genomad_db >> $wd/log/genomad.log 2>&1 ::: $(cat GenomeID)

# combine all the predicted virus and change header:
cat $wd/12_viral_genome/genomad_output/*fna > $wd/12_viral_genome/genomad_predicted_virus.fna
sed -i 's/>/>genomad_/' $wd/12_viral_genome/genomad_predicted_virus.fna
```

**virsorter**
use the second database, and only keep the genomes with level 1,2,4,5
```bash
# parallel run:
time parallel -j 3 --xapply \
wrapper_phage_contigs_sorter_iPlant.pl \
-f $wd/03_genome/{1}.fna \
--diamond \
--db 2 \
--wdir $wd/12_viral_genome/virSorter_output/{1} \
--ncpu 32 --data-dir $db/VirSorter/virsorter-data  >> $wd/log/virsorter.log 2>&1 ::: $(cat GenomeID)

# combine all the predicted virus and cange header：
cat $wd/12_viral_genome/virSorter_output/*/Predicted_viral_sequences/VIRSorter_cat-1.fasta > $wd/12_viral_genome/virsorter_predicted_cat-1.fasta
cat $wd/12_viral_genome/virSorter_output/*/Predicted_viral_sequences/VIRSorter_cat-2.fasta > $wd/12_viral_genome/virsorter_predicted_cat-2.fasta
cat $wd/12_viral_genome/virSorter_output/*/Predicted_viral_sequences/VIRSorter_cat-3.fasta > $wd/12_viral_genome/virsorter_predicted_cat-3.fasta
cat $wd/12_viral_genome/virSorter_output/*/Predicted_viral_sequences/VIRSorter_prophages_cat-4.fasta > $wd/12_viral_genome/virsorter_predicted_cat-4.fasta
cat $wd/12_viral_genome/virSorter_output/*/Predicted_viral_sequences/VIRSorter_prophages_cat-5.fasta > $wd/12_viral_genome/virsorter_predicted_cat-5.fasta
cat $wd/12_viral_genome/virSorter_output/*/Predicted_viral_sequences/VIRSorter_prophages_cat-6.fasta > $wd/12_viral_genome/virsorter_predicted_cat-6.fasta

sed -i 's/>/>virsorter_/' $wd/12_viral_genome/virsorter*.fasta

cat  $wd/12_viral_genome/virsorter_predicted_cat-1.fasta $wd/12_viral_genome/virsorter_predicted_cat-2.fasta \
 $wd/12_viral_genome/virsorter_predicted_cat-4.fasta  $wd/12_viral_genome/virsorter_predicted_cat-5.fasta > $wd/12_viral_genome/virsorter_predicted_virus.fna
```

**deepvirfinder**
The input files are fasta files of metagenome-assembled contigs.  
With the threshold p < 0.05 & score > 0.5 
```bash
# Sore the contigs.
time python $software/DeepVirFinder/dvf.py \
      -i $wd/02_bins/megahit/final.contigs.fa \
      -o $wd/12_viral_genome/deepvirfinder_output/ \
      -m $db/DeepVirFinder/models \
      -l 1000 \
      -c 24 >> $wd/log/dvf.log 2>&1

# the output file is a tab seperated table with contig name, length, score and pvalue.
# filter potential viral contigs with p < 0.05 and score > 0.5
awk -F "\t" '{if($3*1>0.5 && $4*1<0.05)print $0}' $wd/12_viral_genome/deepvirfinder_output/final.contigs.fa_gt1000bp_dvfpred.txt  > $wd/12_viral_genome/deepvirfinder_output/p0.05_s0.5_filtered_gt1000bp_dvfpred.txt

cut -f 1 $wd/12_viral_genome/deepvirfinder_output/p0.05_s0.5_filtered_gt1000bp_dvfpred.txt | tail -n+2 >  $wd/12_viral_genome/dvf_predicted_viral_contigID

# subset viral contigs and add label
seqtk subseq \
    $wd/02_bins/megahit/final.contigs.fa \
    $wd/12_viral_genome/dvf_predicted_viral_contigID > $wd/12_viral_genome/dvf_predicted_virus.fna

sed -i 's/^/dvf_/' $wd/12_viral_genome/dvf_predicted_viral_contigID
sed -i 's/>/dvf_/' $wd/12_viral_genome/dvf_predicted_virus.fna
```

## Remove host region and filter viral genomes:

**Host removal and quality evaluation:**
```bash
# cat all predicted virus together.
cat $wd/12_viral_genome/genomad_predicted_virus.fna $wd/12_viral_genome/virsorter_predicted_virus.fna $wd/12_viral_genome/dvf_predicted_virus.fna > $wd/12_viral_genome/all_predicted_virus.fna

# run checkV to remove host region
time checkv contamination \
        $wd/12_viral_genome/all_predicted_virus.fna \
        $wd/12_viral_genome/checkv_output/viral_host_remove/ \
        -t 32 -d /mnt/m2/dairui/software/checkv_data/checkv-db-v1.5/ >> $wd/log/checkv.log 2>&1

# run end-to-end mode to identify complete genomes and evaluate the genome quality.
cat $wd/12_viral_genome/checkv_output/viral_host_remove/*fna > $wd/12_viral_genome/checkv_output/all_predicted_virus_rmhost.fna

time checkv end_to_end \
    $wd/12_viral_genome/checkv_output/all_predicted_virus_rmhost.fna \
    $wd/12_viral_genome/checkv_output/end_to_end/ \
    -t 32 -d /mnt/m2/dairui/software/checkv_data/checkv-db-v1.5/ --remove_tmp >> $wd/log/checkv.log 2>&1
```
**Filter candidate viral genomes**
```bash
awk -F "\t" '{if($2*1>=1000 && $7/$5<0.3)print $0}' $wd/12_viral_genome/checkv_output/end_to_end/quality_summary.tsv | \
grep -E  "Complete|High-quality|Medium-quality" > $wd/12_viral_genome/viral_genomeID

seqtk subseq \
    $wd/12_viral_genome/checkv_output/all_predicted_virus_rmhost.fna \
    $wd/12_viral_genome/viral_genomeID > $wd/12_viral_genome/overMedium_predicted_virus_rmhost.fna
```

**Sequence dereplication**  
The script blastani.py and cluster.py from the MGV were used to perform sequence dereplication using a greedy, centroid-based algorithm based on 100% sequence similarity and 100% coverage of short sequences.
```bash
# blastdb build
time makeblastdb \
    -in $wd/12_viral_genome/overMedium_predicted_virus_rmhost.fna \
    -out $wd/12_viral_genome/blastdb/overMedium_predicted_virus_rmhost \
    -dbtype nucl

# all-vs-all blast
time blastn \
    -query $wd/12_viral_genome/overMedium_predicted_virus_rmhost.fna \
    -db $wd/12_viral_genome/blastdb/overMedium_predicted_virus_rmhost \
    -out $wd/12_viral_genome/overMedium_predicted_virus_rmhost_blast.tsv \
    -outfmt '6 std qlen slen' -max_target_seqs 30000 -perc_identity 90 -num_threads 96

# calculate ani
time python blastani.py \
    -i $wd/12_viral_genome/overMedium_predicted_virus_rmhost_blast.tsv \
    -o $wd/12_viral_genome/overMedium_predicted_virus_rmhost_ani.tsv

# 100% + 100% dereplicate
time python cluster.py \
    --fna $wd/12_viral_genome/overMedium_predicted_virus_rmhost.fna \
    --ani $wd/12_viral_genome/overMedium_predicted_virus_rmhost_ani.tsv \
    --out $wd/12_viral_genome/CRVC_genomes.tsv \
    --min_ani 100 --min_qcov 0 --min_tcov 100
```