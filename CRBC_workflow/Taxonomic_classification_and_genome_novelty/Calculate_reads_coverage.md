# Calculate reads coverage of catalogs:
Evaluate the metagenomic reads use efficiency of each genomic catalogs.
   
**Metagenomic data normalize to 1G**  
The metagenomic sequencing data are 150-PE. 
```bash
for file in $(cat SampleID)
do
    echo $file
    seqtk sample -s 1013 $wd/01_kneaddata/further_clean/further_clean_"$file"_kneaddata_paired_1.fastq 3333333 > $wd/01_kneaddata/subset1G_reads/subset_1G_reads_further_clean_"$file"_paired_1.fastq
    seqtk sample -s 1013 $wd/01_kneaddata/further_clean/further_clean_"$file"_kneaddata_paired_2.fastq 3333333 > $wd/01_kneaddata/subset1G_reads/subset_1G_reads_further_clean_"$file"_paired_2.fastq
done
```

**Salmon index**   
```bash
time salmon index -p 96 \
    -t $db/ref_dataset.fasta \
    -i $wd/06_salmon/index/ref_dataset_index \
    --type puff  -k 31 >> $wd/log/salmon_index.log 2>&1
```  

**Salmon alignment**  
```bash
# salmon quant
time parallel -j 2 --xapply \
    salmon quant \
        -i $wd/06_salmon/index/ref_dataset_index \
        --libType A -p 32 --meta --validateMappings \
        -1 $wd/01_kneaddata/subset1G_reads/subset_1G_reads_further_clean_{1}_paired_1.fastq \
        -2 $wd/01_kneaddata/subset1G_reads/subset_1G_reads_further_clean_{1}_paired_2.fastq \
        -o $wd/06_salmon/alignment/{1} >> $wd/log/salmon_quant.log 2>&1 \
::: $(cat SampleID)

# summary
for i in $(cat SampleID)
do
echo -n $i"\t" >> $wd/06_salmon/mapping_ratio.txt
echo `grep "Mapping rate" $wd/06_salmon/alignment/$i/logs/salmon_quant.log | cut -f 2 -d '=' | sed 's/%//'`>> $wd/06_salmon/mapping_ratio.txt
done
```
