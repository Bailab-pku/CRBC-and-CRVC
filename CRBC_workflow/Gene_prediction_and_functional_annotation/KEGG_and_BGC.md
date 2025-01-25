# Gene prediction and functional annotation:
Predict the potential genes, annotate the functional genes, predict the BGCs and anti-phage defense systems.
   
## Gene prediction:
Before functional annotation, predict the reliable protein-coding genes first. Output files including fasta files of genes and proteins, and gff file.
```bash
# run prodigal.
time parallel -j 24 --xapply \
    prodigal \
        -i $wd/03_genome/{1}.fna \
        -o $wd/07_prodigal/{1}.gff  \
        -a $wd/07_prodigal/{1}.faa \
        -d $wd/07_prodigal/{1}.ffn \
        -p single -f gff >> $wd/log/prodigal.log 2>&1 ::: $(cat GenomeID)

# summarise result, gene number, length and complete gene number:
echo "file num_seqs sum_len min_len avg_len max_len complete_num" > $wd/07_prodigal/CRBC_prodigal_summary.txt
cat GenomeID | parallel -j 24 --xapply \
    echo $(seqkit stat -b sep_genome/{1}.ffn | tail -1) $(grep -c "partial=00" sep_genome/{1}.ffn) | cut -f 1,4-9 -d ' ' >> $wd/07_prodigal/CRBC_prodigal_summary.txt
```
## KEGG annotation:
- **Step 1: KEGG protein annotation**
```bash
# cat all data together, speed up the annotation step.
cat $wd/07_prodigal/*faa > $wd/CRBC_prodigal_genes.faa
# build index:
time diamond makedb --threads 64 \
    --in $db/kegg/kegg_2021_April.fasta \
    --db $db/kegg/kegg_2021_April >> $wd/log/diamond_makedb.log 2>&1

# run diamond:
time diamond blastp \
    --query $wd/CRBC_prodigal_genes.faa \
    --db $db/kegg/kegg_2021_April.dmnd \
    --out $wd/08_diamond/CRBC_kegg_1.0e-5.fm6 \
    --outfmt 6 --threads 48 --log --max-target-seqs 1 \
    --evalue 1e-5 --sensitive --block-size 6 \
    --index-chunks 1 >> $wd/log/diamond.log 2>&1
```

- **Step 2: Generate genome annotation table**  
The input files were downloaded from the official KEGG website.
```bash
# Annotation results split by genome
cat GenomeID | parallel -j 48 --xapply \
    grep "{1}___" $wd/08_diamond/CRBC_kegg_1.0e-5.fm6 > $wd/08_diamond/temp/{1}_kegg_1.0e-5_id50.fm6

## keep the annotation with the thresholds of evalue 1.0e-5 and identity 50%
# summarise KO, pathway and superpathway annotation for each sample
cat GenomeID | parallel -j 48 --xapply \
    "python3 $script/fm6_to_KEGG_levels.py \
    --KEGG_blast_fm6 $wd/08_diamond/temp/{1}_kegg_1.0e-5_id50.fm6 \
    --Kegg_gene_to_KO_map $db/kegg/ko_gene_map.txt \
    --Kegg_KO_to_pathway $db/kegg/pathway_ko.txt \
    --Kegg_pathway_hierarchical_file $db/kegg/pathway_htext.txt \
    --evalue 1.0e-5 --identity 50 \
    --KO_output $wd/08_diamond/temp/{1}_KO_1.0e-5_id50plus.txt \
    --Pathway_output $wd/08_diamond/temp/{1}_pathway_1.0e-5_id50plus.txt \
    --SuperPath_output $wd/08_diamond/temp/{1}_superpath_1.0e-5_id50plus.txt \
    --Global_catagory_output $wd/08_diamond/temp/{1}_global_1.0e-5_id50plus.txt \
    --genomeID {1}"

# table integrate:
## KO
time python3 $script/merge_thousands_KO_files_into_one.py \
    --input_dir $wd/08_diamond/temp/ --pattern "KO_1.0e-5_id50plus.txt" \
    --how outer --key KO \
    --output $wd/08_diamond/CRBC_genome_KO_1.0e-5_id50plus.txt
## Pathway
time python3 $script/merge_thousands_KO_files_into_one.py \
    --input_dir $wd/08_diamond/temp/ --pattern "pathway_1.0e-5_id50plus.txt" \
    --how outer --key Pathway \
    --output $wd/08_diamond/CRBC_genome_pathway_1.0e-5_id50plus.txt
## Superpathway
time python3 $script/merge_thousands_KO_files_into_one.py \
    --input_dir $wd/08_diamond/temp/ --pattern "superpath_1.0e-5_id50plus.txt" \
    --how outer --key Superpathway \
    --output $wd/08_diamond/CRBC_genome_superpath_1.0e-5_id50plus.txt
## Global
time python3 $script/merge_thousands_KO_files_into_one.py \
    --input_dir $wd/08_diamond/temp/ --pattern "global_1.0e-5_id50plus.txt" \
    --how outer --key Global_catagory \
    --output $wd/08_diamond/CRBC_genome_global_1.0e-5_id50plus.txt

rm -r $wd/08_diamond/temp/
```

## Detect Biosynthetic Gene Clusters (BGCs)
- **Run antismash to predict potential BGCs**
```bash
## subset contigs with length >= 5kb
## run antismash
time parallel -j 2 --xapply \
    antismash -c 24 \
        --allow-long-headers --asf --minlength 5000 \
        --cb-general --cb-knownclusters --cb-subclusters --pfam2go \
        --output-dir $wd/09_BGCs/{1} \
        --genefinding-gff3 $wd/07_prodigal/{1}.gff \
        $wd/03_genome/{1}.fna" ::: $(cat GenomeID)

# summary results:
cd $wd/09_BGCs/
echo "GenomeID\tBGCID\tContig\tStart\tEnd\tLength" > $wd/CRBC_BGCs_basic_info.txt
for i in `cat GenomeID`; do
    for j in `ls $wd/09_BGCs/$i/*region*.gbk`;do
        genome_and_gbk=`ls $j | sed "s,/,\t,"`;contig=`ls $j | sed "s,/,\t,;s,\.region,\t," | cut -f 2`;start=`grep "Orig. start" $j | cut -f 3 -d ':' | sed 's/ //g'`;end=`grep "Orig. end" $j | cut -f 3 -d ':' | sed 's/ //g'`;
        printf "%s\t%s\t%s\t%s\t%s\n" "$genome_and_gbk" "$contig" "$start" "$end" "$((end-start))" >> $wd/CRBC_BGCs_basic_info.txt
    done
done
```
