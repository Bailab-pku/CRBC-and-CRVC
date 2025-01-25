# Quantification of bacterial genomes in metagenomic samples:
- This script is used to measure the presence or absence of crop root bacterial genomes in macrogenomic samples, also  the abundance information, etc.
- Referring to the InStrain profile, we used representative genomes at the specie level as reference databases, to prevent mapping cofusion.

##  Mapping sequences to reference genome databases  

**Step 1: build bowtie2 index**
```bash
# all representative species were conbined.
cat $wd/05_drep/0.95_drep/dereplicated_genomes/*fna > $wd/05_drep/0.95_drep/representative_species.fna
# build index:
time bowtie2-build \
    -f --seed 218 --threads 96 --large-index \
    $wd/05_drep/0.95_drep/representative_species.fna \
    $wd/06_bowtie/index/representative_species_index >$wd/log/bowtie_index.log 2>&1
```

**Step 2: align the clean sequences to genomes**   
Use options `--non-deterministic` and `--very-sensitive`, with default end-to-end read alignment.
```bash
time parallel -j 2 --xapply \
    bowtie2 \
        --non-deterministic -q --threads 48 --very-sensitive \
        -x $wd/06_bowtie/index/representative_species_index \
        -1 $wd/01_kneaddata/further_clean/further_clean_{1}_kneaddata_paired_1.fastq \
        -2 $wd/01_kneaddata/further_clean/further_clean_{1}_kneaddata_paired_2.fastq  \
        -S $wd/06_bowtie/sam/{1}_end.sam >> $wd/log/bowtie_align.log 2>&1 \
    ::: $(cat SampleID)
```
**Step3: sort and convert to bam by reference**  
```bash
time cat SampleID | parallel -j 3 --xapply \
    'file={};wd="~"; \
    samtools sort -@ 36 -O bam $wd/06_bowtie/sam/"$file"_end.sam  -o $wd/06_bowtie/bam/"$file"_end_sorted.bam'
```

## Calculate the abundance of genomes
- **InStrain quick-profile**  
Use coverm profile, give the genome breadth, coverage and  
```bash
# generate a file to let inStrain know which scaffolds came from which genomes, from drep
time $wd/envs/drep/bin/parse_stb.py \
--reverse -f $wd/05_drep/0.95_drep/dereplicated_genomes/* \
-o $wd/05_drep/0.95_drep/representative_species.stb

# run quick-profile
time cat SampleID | parallel -j 4 --xapply \
    'wd="~";
    inStrain quick_profile \
        $wd/06_bowtie/bam/{1}_end_sorted.bam \
        $wd/05_drep/0.95_drep/representative_species.fna \
        -p 32 -o $wd/07_InStrain/quick_coverm/{1} \
        -s $wd/05_drep/representative_species.stb'

# for the results are seperated according to sample names,run in-house script to merge them together.
time python3 $scripts/merge_InStrain_output.py --input_dir $wd/07_InStrain/quick_coverm/ --output $wd/07_InStrain/InStrain_quick_profile

# output the breadth, coverage and count for each genome.
sed -i 's/\.fna//g;s/genome/GenomeID/' $wd/07_InStrain/InStrain_quick_profile*tsv
```
We considered a genome to be detected in a sample only if its breadth was greater than or equal to 50%, and reads count were used to calculate the relative abundance of each genome.


- **InStrain complete profile** 
Use the `--database_mode` to run InStrain profile, can get more information. The InStrain will automatically filter the ANI and low quality reads.
```bash
# prepare gene file for gene profile
tail -n+2 $wd/05_drep/0.95_drep/data_tables/de_replicated_genomeID.txt | parallel -j 24 --xapply "ln -s  $wd/07_prodigal/{1}.ffn $wd/temp/"
cat $wd/temp/*.ffn > $wd/07_InStrain/representative_species_gene.fna
rm $wd/temp/*.ffn

time cat SampleID | parallel -j 4 --xapply \
    'wd="~";
    inStrain quick_profile \
        $wd/06_bowtie/bam/{1}_end_sorted.bam \
        $wd/05_drep/0.95_drep/representative_species.fna \
        -p 32 -o $wd/07_InStrain/quick_coverm/{1} \
        -s $wd/05_drep/representative_species.stb'

# run profile:
time cat SampleID | parallel -j 4 --xapply \
    'wd="~";
        inStrain profile \
        $wd/06_bowtie/bam/{1}_end_sorted.bam \
        $wd/05_drep/0. 95_drep/representative_species.fna \
        -p 32 -o $wd/07_InStrain/complete_profile/{1} \
        -g $wd/07_InStrain/representative_species_gene.fna \
        --skip_plot_generation --database_mode \
        -s $wd/05_drep/representative_species.stb'
```
