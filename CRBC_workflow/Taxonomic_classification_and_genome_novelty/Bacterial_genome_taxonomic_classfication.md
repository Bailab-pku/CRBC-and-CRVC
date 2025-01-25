# Bacterial_genome_taxonomic_classfication:
Genome Taxonomy Database (GTDB) r207 was used to classify the genome taxonomy.
   
**GTDB annotation**  
Use the gtdbtk `classify_wf` to classify genomes. 
```bash
time gtdbtk classify_wf \
        --genome_dir $wd/03_genome/ \
        --out_dir $wd/04_gtdbtk_r207 \
        --cpus 32 --extension fna --full_tree --debug \
        --prefix CRBC >> $wd/log/gtdbtk_classify_wf.log 2>&1
```

**GTDB to NCBI**  
Use the `gtdb_to_ncbi_majority_vote.py` file to transform gtdb taxonomy into ncbi taxonomy. Download the corresponding GTDB r207 metadata first. 
```bash
time gtdb_to_ncbi_majority_vote.py \
 --gtdbtk_output_dir $wd/04_gtdbtk_r207/ \
 --ar53_metadata_file $db/GTDBr207/metadata/release207/ar53_metadata_r207.tsv \
 --bac120_metadata_file $db/GTDBr207/metadata/release207/bac120_metadata_r207.tsv \
 --output_file $wd/04_gtdbtk_r207/CRBC_gtdb_to_ncbi.txt \
 --gtdbtk_prefix CRBC
```  
**Generate genome taxonomy data**  
```bash
# a. gtdbtk-version
tail -n+2 $wd/04_gtdbtk_r207/CRBC_gtdb_to_ncbi.txt |cut -f 1-2|sed 's/;/\t/g'|sed '1 s/^/ID\tKingdom_gtdb\tPhylum_gtdb\tClass_gtdb\tOrder_gtdb\tFamily_gtdb\tGenus_gtdb\tSpecies_gtdb\n/' > $wd/04_gtdbtk_r207/CRBC_gtdb_tax.txt
# b. ncbi-version
tail -n+2 $wd/04_gtdbtk_r207/CRBC_gtdb_to_ncbi.txt |cut -f 1,3|sed 's/;/\t/g'|sed '1 s/^/ID\tKingdom_ncbi\tPhylum_ncbi\tClass_ncbi\tOrder_ncbi\tFamily_ncbi\tGenus_ncbi\tSpecies_ncbi\n/' > $wd/04_gtdbtk_r207/CRBC_ncbi_tax.txt
```
