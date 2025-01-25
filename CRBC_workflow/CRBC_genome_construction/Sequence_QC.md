# Sequence QC:
- For bacterial isolates, we use Trimmomatic to remove low quality reads
- For metagenomic sequensing data, the QC step including remove reads with low quality and reads derived from host plants, conducted by KneadData and Kraken2.

**Isolates**
```bash
# parallel run:
time parallel -j 4 --xapply \
    trimmomatic PE \
        -threads 16 -phred33 -summary $wd/01_trimmomatic/{1}_trimmomatic_summary \
        $wd/seq/{1}_1.fq.gz $wd/seq/{1}_2.fq.gz \
        $wd/01_trimmomatic/{1}_paired_1.fq.gz /dev/null \
        $wd/01_trimmomatic/{1}_paired_2.fq.gz /dev/null \
        SLIDINGWINDOW:4:20 MINLEN:100  >> $wd/log/trimmomatic.log 2>&1 \
::: $(cat Genome_ID)
```

**MAGs**
- *Step 1: QC and host remove*
```bash
# run KneadData:
time parallel -j 2 --xapply \
    kneaddata \
        -i $wd/seq/{1}_1.fq.gz -i $wd/seq/{1}_2.fq.gz \
        -o $wd/01_kneaddata/ \
        -db $host_db/bowtie_index/ \
        -t 36 -v \
        --trimmomatic $software/bin \
        --trimmomatic-options="SLIDINGWINDOW:4:20" \
        --trimmomatic-options="MINLEN:50" \
        --bowtie2  $software/bin \
        --bowtie2-options="--very-sensitive" \
        --reorder --remove-intermediate-output >> $wd/log/kneaddata.log 2>&1 \
::: $(cat $wd/SampleID)

# QC summary:
kneaddata_read_count_table --input $wd/01_kneaddata/ \
    --output $wd/01_kneaddata/kneaddata_statistics.txt
```

- *Step 2: Check host remove state*  
The sequences were classified using Kraken2 to check if the host sequence was removed completely in the first step, and if not, further host removement was conducted.  
```bash
# clean reads classification
time parallel -j 2 --xapply \
    kraken2 \
        --db $kraken2_db/ \
        --paired $wd/01_kneaddata/{1}_*_kneaddata_paired_*.fastq \
        --threads 32 --use-mpa-style --report-zero-counts \
        --report $wd/01_kneaddata/kraken2/{1}_report \
        --output $wd/01_kneaddata/kraken2/{1}_output >> $wd/log/further_kraken2.log 2>&1 \
::: $(cat $wd/SampleID)
```

- *Step 3: Further clean based on kraken2 results*  
Note the parameters `-t 4577 --exclude` means to find all reads not matching taxid 4577 (taxid of host *Zea mays*).
```bash
# run Kraken2 scripts to extact non-host reads
for file in $(cat $wd/SampleID); 
do echo $file; 
python3.8 krakentools/KrakenTools-master/extract_kraken_reads.py \
-k $wd/01_kneaddata/kraken2/"$file"_output \
-s $wd/01_kneaddata/"$file"_*_kneaddata_paired_1.fastq \
-s2 $wd/01_kneaddata/"$file"_*_kneaddata_paired_2.fastq \
-t 4577 --exclude --fastq-output \
-o $wd/01_kneaddata/further_clean/further_clean_"$file"_kneaddata_paired_1.fastq \
-o2 $wd/01_kneaddata/further_clean/further_clean_"$file"_kneaddata_paired_2.fastq; done >> $wd/log/further_clean.log 2>&1

# further clean summary
cd $wd/01_kneaddata/further_clean/
## seqkit summarise seq number and length
for file in $(cat $wd/SampleID)
do
echo $file
seqkit stat further_clean_"$file"_kneaddata_paired_1.fastq >> further_clean_reads_stats.txt
done
## summary
cat further_clean_reads_stats.txt | sed 's,further_clean_,,g' | sed 's,_kneaddata_paired_1.fastq,,g' | sed 's,  ,\t,g' | cut -f 1,4 | grep -v 'file' | awk 'BEGIN{print "SampleID" "\t" "further_clean_seq_sum"}1' > further_clean_sum.txt
```
