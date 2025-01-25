# Quantify the growth rates of genomes
Use the peak-to-trough ratio (PTR) to measure the growth dynamics of bacteria.
## iRep:  
The input file should be sorted sam files.
```bash
# convert bam to sam
time cat SampleID | parallel -j 3 --xapply \
    'wd="~";
    samtools view -h $wd/06_bowtie/bam/{}_end_sorted.bam > $wd/06_bowtie/sam/ref_sorted/{}_sorted.sam -@ 32'

# run iRep:
time iRep \
    -f $wd/05_drep/0.95_drep/dereplicated_genomes/*fna \
    -t 32 -s $wd/06_bowtie/sam/ref_sorted/*_sorted.sam \
    -o $wd/08_growth_rate/irep/iRep_result --no-plot
```


## DEMIC：
Similar to iRep.
 ```bash

time perl $software/DEMIC_v1.0.2/DEMIC.pl \
    -S  $wd/06_bowtie/sam/ref_sorted/ \
    -T 32 \
    -F $wd/05_drep/0.95_drep/dereplicated_genomes/ \
    -O $wd/08_growth_rate/demic/
 ```
