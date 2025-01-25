# Species and Genus level clustering:
Genomes from CRVC, IMG/VR, RefSeq, MGV and GOV2 were combined for species- and genus- level clustering.   
The clustering strategies and scripts refer to the methods in MGV.

## Species level clustering:
According to MIUViG standards, a threshold of 95% ANI and coverage of 85% of the shorter sequence were used for species-level (vOTU) clustering. 
1. All-vs-all local comparisons between viral sequences were performed using the blastn package from BLAST v2.5.0138. 
2. The script blastani.py and cluster.py from the MGV were used to perform species-level clustering using a greedy, centroid-based algorithm. A total of 7,653 species-level clusters were identified in the CRVC.   

**All-vs-all blast of all the genomes**  
The thresholds are 90% identity and 30,000 targeted sequences.
```bash
# blast database build:
time makeblastdb -in $wd/12_viral_genome/sequence.fna -out $wd/12_viral_genome/blastdb/blastdb -dbtype nucl

# all-vs-all blast:
time blastn -query $wd/12_viral_genome/sequence.fna \
        -db $wd/12_viral_genome/blastdb/blastdb \
        -num_threads 96 \
        -out $wd/12_viral_genome/blast.tsv \
        -outfmt '6 std qlen slen' -max_target_seqs 30000 -perc_identity 90
```

**ANI calculation and clustering**  
```bash
# calculate all-vs-all ani of whole genomes.
time python blastani.py -i $wd/12_viral_genome/blast.tsv -o $wd/12_viral_genome/ani.tsv

# species level sequence：
time python cluster.py \
    --fna $wd/12_viral_genome/sequence.fna \
    --ani $wd/12_viral_genome/ani.tsv \
    --out $wd/12_viral_genome/clusters_viral_species.tsv \
    --min_ani 95 --min_qcov 0 --min_tcov 85

# label the clusters:
cat -n $wd/12_viral_genome/clusters_viral_species.tsv | \
    sed 's/^/vOTU_/;s/ //g' | cut -f 1,3 | \
    awk -v OFS='\t' '{split($NF, arr, ","); for(i=1; i<=length(arr); i++) print $1, arr[i]}' \
> $wd/12_viral_genome/clusters_viral_species_split_withID.tsv
```

## Genus level clustering:
Clustering was performed based on protein similarity and the protein share network between genomes. Viruses from the CRVC and public databases were represented by representative genomes at the species level.
1. All-vs-all protein alignments were conducted using DIAMOND BLASTP v2.0.15.153116
2. Calculated the number of shared proteins between genomes, the proportion of proteins shared, and the average AAI of shared proteins between each genome pair. 
3. The score between genome pairs was calculated based on mincov× meanaai and used as the network edge between the genome pairs. 
4. Finally, genome pairs were filtered based on the set threshold, and mcl v14-137141 was used to perform clustering with the option ‘--abc -I’.   


   We benchmarked combinations of different filtering thresholds, including protein sharing ratios of 10%, 15%, 20%, 25%, and 30%, and note that when the viral genome was too large, at least 20 proteins needed to be shared between them. We also tested average AAI values of 20%, 30%, 40%, 50%, 60%, and 70% as well as MCL inflation factors of 1.1, 1.2, 1.4, 2.0, 4.0, and 6.0. We evaluated their combinations by comparing them against the taxonomy and clustering results derived from the RefSeq viruses, which served as our mock dataset, using Adjusted Mutual Information scores142. Finally, we selected a protein network sharing threshold at average AAI of 30% and an inflation factor of 4.0 for genus-level clustering, with precision and recall thresholds set to 0.85 and 0.88, respectively. Viral species- or genus-level clusters that could not be clustered with other viruses were defined as novel.
