# Genome assmeble:
- For bacterial isolates, SPAdes were used to *de novo* assemble genomes.
- For metagenomic sequensing data, Megahit were used to *de novo* assemble microbial contigs, and MetaWRAP pipeline were used for binning and refinement. 

**Isolates**  
Reads were first assemble into contigs using SPAdes, an seqkit was used to extract contigs length longer than 1 kb.
```bash
# parallel run:
time parallel -j 2 --xapply \
    spades.py \
        --pe1-1 $wd/01_trimmomatic/{1}_paired_1.fq.gz \
        --pe1-2 $wd/01_trimmomatic/{1}_paired_2.fq.gz \
        --isolate --careful --cov-cutoff auto \
        -o $wd/02_spades/temp/{1} >> $wd/log/{1}_SPAdes.log 2>&1; \
    seqkit seq -m 1000 $wd/02_spades/temp/{1}/contigs.fasta > $wd/02_spades/{1}.fna \
::: $(cat Genome_ID)
```

**MAGs**  
Construction of metagenome-assembled genomes including 3 steps, 1) assemble; 2) binning; 3) bin refinement.
- *Step 1: Contig assemble*
```bash
# co assemble, list all sample sequences, jeep the same order：
ls $wd/01_kneaddata/further_clean/further_clean_"$file"_kneaddata_paired_1.fastq > $wd/02_bins/root_R1_list
sed 's/_1.fastq/_2.fastq/' $wd/02_bins/root_R1_list > $wd/02_bins/root_R2_list
# run megahit:
time megahit \
    -t 64 --continue --min-contig-len 1000 --presets meta-large \
        -1 `cat root_R1_list | paste -s -d ','` \
        -2 `cat root_R2_list | paste -s -d ','` \
        -o $wd/02_bins/megahit/ >> $wd/log/megahit.log 2>&1
```

- *Step 2: Binning*  
Use both MetaBAT2 and MaxBin2 to binning the contigs.
```bash
# binning module:
time metawrap binning -t 48 \
    -o $wd/02_bins/binning \
    -a $wd/02_bins/megahit/final.contigs.fa \
    --metabat2 --maxbin2 \
    $wd/01_kneaddata/further_clean/*.fastq >> $wd/log/binning.log 2>&1
```

- *Step 3: Bin refinement*  
Improve the bin quality and remove redundant bins. The threshold of completeness and contamination are 50% and 10% respectively.
```bash
# refinement module
time metawrap bin_refinement -t 48 -c 50 -x 10 \
  -o $wd/02_bins/bin_refinement/ \
  -A $wd/02_bins/binning/metabat2_bins/ \
  -B $wd/02_bins/binning/maxbin2_bins/ >> $wd/log/bin_refinement.log 2>&1
```
