# Predict defense systems and CRISPR array:
Predict the antiviral systems and crispr array in bacterial genomes.
   
## Antiviral systems
Predict the antiviral systems of each genome using DefenseFinder.
```bash
# run defense-finder
time parallel -j 4 --xapply \
    defense-finder run $wd/07_prodigal/{1}.faa \
        -w 24 -o $wd/10_defenseFinder/{1} \
        --models-dir $defenseFinder >> $wd/log/defenseFinder.log ::: $(cat GenomeID)

# summarise result, add genome ID:       
echo "GenomeID\tsys_id\ttype\tsubtype\tsys_beg\tsys_end\tprotein_in_syst\tgenes_count\tname_of_profiles_in_sys" > $wd/10_defenseFinder/CRBC_defense_finder_systems.tsv
time parallel -j 4 --xapply \
    sed "s/^/{1}\t/" $wd/10_defenseFinder/{1}/defense_finder_systems.tsv | tail -n+2 >> $wd/10_defenseFinder/CRBC_defense_finder_systems.tsv ::: $(cat SampleID)
```
## CRISPR array:  
**Step1: Predict crispr array using CRISPRCasFinder**  
```bash
cd $wd/11_crisprcasfinder/
ln -s $db/singularity/CrisprCasFinder.simg ./
# note that genome file should in the doc
mkdir -p $wd/11_crisprcasfinder/genome
for i in `cat GenomeID`
do
   cp $wd/03_genome/${i}.fna $wd/11_crisprcasfinder/genome/${i}.fna
    singularity exec -B $PWD CrisprCasFinder.simg \
       perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl \
       -so /usr/local/CRISPRCasFinder/sel392v2.so \
       -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 \
       -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
       -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv \
       -cpuM 48 -log -q -cas -def G -out $wd/11_crisprcasfinder/${i} \
       -in $wd/11_crisprcasfinder/genome/${i}.fna
   rm $wd/11_crisprcasfinder/genome/${i}.fna
 done
```
**Step2: Summarise the confidence level of crispr arrays**  
The array informations are in the document `TSV/Crisprs_REPORT.tsv`. Note that in the gff file, the CRISPR ID consists of the contigID, the start position of the array, so we need to generate a new ID accordingly.  
```bash
# summary the genome, contig and confidence level of crispr arrays.
echo "GenomeID\tContigID\tCRISPR_Id\tCRISPR_Start\tCRISPR_End\tCRISPR_Length\tConsensus_Repeat\tSpacers_Nb\tEvidence_Level" > $wd/11_crisprcasfinder/CRBC_Crisprs_REPORT.tsv

time cat GenomeID | parallel -j 24 --xapply \
"grep -v ^$ $wd/11_crisprcasfinder/{1}/TSV/Crisprs_REPORT.tsv | tail -n+2 | cut -f 2,5,6,7,8,11,15,27 | sed 's/^/{1}\t/' >> $wd/11_crisprcasfinder/CRBC_Crisprs_REPORT.tsv"

# add array name
awk -F "\t" -v OFS="\t" '{$10=$2"_"$4"_"$5;print $0}' $wd/11_crisprcasfinder/CRBC_Crisprs_REPORT.tsv | sed "s/Sequence_CRISPR_Start_CRISPR_End/ArrayID/" > $wd/11_crisprcasfinder/CRBC_Crisprs_REPORT_with_name.tsv

# filter arrays with a confidence level of 4
awk -F "\t" 'NR==1 || ($9==4)' $wd/11_crisprcasfinder/CRBC_Crisprs_REPORT_with_name.tsv > $wd/11_crisprcasfinder/CRBC_Crisprs_evidence4_REPORT_with_name.tsv
```
**Step3: summary spacer information**  
The spacer sequences are in the gff file.
```bash
# all predicted crispr array information.
cat $wd/11_crisprcasfinder/*/*_gff.txt > $wd/11_crisprcasfinder/CRBC_CRISPR_array.gff
## summarise spacers' information
grep CRISPRspacer $wd/11_crisprcasfinder/CRBC_CRISPR_array.gff | cut -f 1,4,5,9 |  sed 's/sequence=/\t/;s/;Name=/\t/;s/;Parent=/\t/;s/;ID=/\t/' | tr -s '\t'| cut -f 1-4,6 |  awk -F "\t" -v OFS="\t" '{$6=$5"___CRISPRspacer___"$2"___"$3;print $0}' | sed '1i ContigID\tCRISPR_Start\tCRISPR_End\tSequence\tArrayID\tSpacerID' > $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_metadata.txt
## generate fasta file of spacers:
paste -d "\n" <(cut -f 6 $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_metadata.txt | tail -n+2 | sed 's/^/>/' ) <(cut -f 4 $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_metadata.txt | tail -n+2) > $wd/11_crisprcasfinder/CRBC_CRISPR_spacer.fasta
```

**Step4: extract the spacers in level 4 arrays**  
```bash
# filter spacer ID:
grep -wf <(cut -f 10 $wd/11_crisprcasfinder/CRBC_Crisprs_evidence4_REPORT_with_name.tsv) $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_metadata.txt | cut -f 6 > $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_evidence4ID
# run seqtk to extract level 4 sequences.
time seqtk subseq $wd/11_crisprcasfinder/CRBC_CRISPR_spacer.fasta  $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_evidence4ID > $wd/11_crisprcasfinder/CRBC_CRISPR_spacer_evidence4.fasta
```
