# Species cluster and identify novel species:
- To evaluate the phylogeny of the CRBC, genomes were clustered at species level use dRep, with the threshold of ANI 95% and AF 30%
- To identify novel species and reduce computational costs, we compare our genomes with published genome databases use Mash and Mummer, with the same threshold as before.
   
## Cluster genomes at species level  
Note that since we want to compare the phylogenetic diversity of CRBC with published crop root bacteria, we include all crop root bacteria.
- **Step 1: prepare the files for clustering**
  1. Genome location
  2. Basic information of genomes, including contamination, completeness and N50, for calculating the drep score for genomes.
```bash
# genome information:
awk -F "\t" 'BEGIN{OFS=","}{sub(/$/,".fna", $1);if($9>1*50)print $1,$2,$3,$7}' $wd/04_QC/checkM_result.tsv | \
    sed 's/GenomeID.fna/genome/;s/Completeness/completeness/;s/Contamination/contamination/' > $wd/05_drep/QS50_genomeInfo.csv
# genome location
for i in `tail -n+2 QS50_genomeInfo.csv | cut -f 1 -d ','`
do
ls $wd/03_genome/${i} >> $wd/05_drep/QS50_genome_location.txt
done
```
- **Step 2: Genome de-replication**  
First, redundant genomes were removed with a threshold of 99.9% ANI.
```bash
time dRep dereplicate $wd/05_drep/0.999_drep/ \
    -g $wd/05_drep/QS50_genome_location.txt \
    --genomeInfo $wd/05_drep/QS50_genomeInfo.csv \
    -strW 0 -centW 0 -sa 0.999 -nc 0.30 \
    -p 64 -comp 50 -con 10 >> $wd/log/drep_0.999.log 2>&1
```
Representative genomes were selected based on genome scores (`data_tables/Sdb.csv`) and clustering information (`data_tables/Cdb.csv`), and isolate was surpassed over MAG.
```r
library(tidyverse)
library(reshape2)
wd="~"
setwd(paste(wd,"/05_drep/0.999_drep/",sep=""))
sdb<-read.csv("data_tables/Sdb.csv")
cdb<-read.csv("data_tables/Cdb.csv")
genome_type<-read.table(paste(wd,"/data/metadata.csv",sep=""),sep="\t",col.names=c("genome","type"))
genome_info<-merge(sdb,cdb, by="genome")%>%merge(genome_type,by="genome")
genome_info$type<-factor(genome_info$type, levels=c("Isolates","MAGs"))
# pick the representative genomes:
genome_rep<-genome_info%>%group_by(secondary_cluster)%>%arrange(type, -score)%>%do(head(.,n=1))%>%select(genome)
write.csv(genome_rep, "data_tables/de_replicated_genomeID.txt",quote=F,sep="\t",row.names=F)
```

- **Step 3: Species level cluster**    
Same as the step1 and 2, prepare genome data and run dRep, the only difference is to set the `-sa 0.95`.
```bash
# 1. location and genome info:
grep -wf $wd/05_drep/0.999_drep/data_tables/de_replicated_genomeID.txt $wd/05_drep/QS50_genome_location.txt > $wd/05_drep/QS50_de_repliated_genome_location.txt
grep -wf $wd/05_drep/0.999_drep/data_tables/de_repliated_genomeID.txt $wd/05_drep/QS50_genomeInfo.csv > $wd/05_drep/ QS50_de_repliated_genomeInfo.csv

# 2. run drep
time dRep dereplicate $wd/05_drep/0.95_drep/ \
    -g $wd/05_drep/QS50_de_repliated_genome_location.txt \
    --genomeInfo $wd/05_drep/QS50_de_repliated_genomeInfo.csv \
    -strW 0 -centW 0 -sa 0.95 -nc 0.30 \
    -p 64 -comp 50 -con 10 >> $wd/log/drep_0.95.log 2>&1

# 3. pick representative genomes as step2
# 4. soft link the representative genomes:
cd $wd/05_drep/0.95_drep/dereplicated_genomes
rm *
for i in `tail -n+2 $wd/05_drep/0.95_drep/data_tables/de_replicated_genomeID.txt`
do
ln -s $wd/03_genome/${i} ./
done
```
Species clusters with no pubished genomes were considered as CRBC novel species.

## Species level tree build: 
Construct the crop root bacterial phylogenetic tree.
```bash
# 1. GTDB-based classify these species-level genomes to generate the MSA data:
time gtdbtk classify_wf \
    --genome_dir $wd/05_drep/0.95_drep/dereplicated_genomes \
    --out_dir $wd/06_species_tree/ \
    --cpus 64 --extension fna --full_tree --debug \
    --prefix Crop_bac_species >>$wd/log/drep_gtdbtk.log 2>&1

# 2. Infer the tree:
time gtdbtk infer \
--msa_file $wd/06_species_tree/align/Crop_bac_species.bac120.user_msa.fasta.gz \
--out_dir $wd/06_species_tree/tree_build/ \
--cpu 64 --prefix Crop_bac_species >> $wd/log/tree_bac.log 2>&1
```
The phylogenetic tree was visualized using R package ggtree, phylogenetic diversity was calculated using R package picante.

## Identify novel species:
Use comparison of CRBC with GTDB as an example.
- **Step 1: Mash estimate genome distance**
```bash
# 1. build sketch of gtdb and CRBC
time mash sketch \
  -o $wd/05_novel_species/Mash/gtdb207_reference_sketch1000 \
  -l gtdb_r207_location.txt -p 64

time mash sketch \
  -o $wd/05_novel_species/Mash/CRBC_sketch1000 \
  -l CRBC_location.txt -p 64 &

# 2. Calculate mash distance
time mash dist -p 64 \
$wd/05_novel_species/Mash/CRBC_sketch1000.msh $wd/05_novel_species/Mash/gtdb207_reference_sketch1000.msh > $wd/05_novel_species/Mash/CRBC_GTDB_distances.tab

# 3. pick the genome pairs with lowest genome distance:
time for i in `cat CRBC_GenomeID`
do
grep -w $i $wd/05_novel_species/Mash/CRBC_GTDB_distances.tab | sort -gk3 | head -3 >> $wd/05_novel_species/Mash/CRBC_GTDB_closely_related.txt
done
```
- **Step 2: Calculate genome identity**
```bash
# 1. dnadiff from mummer
 time parallel -j 96 --xapply \
  dnadiff $wd/03_genome/{1} $gtdb_db/{2} -p $wd/05_novel_species/mummer/temp/{1}_{2} \
 ::: `cat $wd/05_novel_species/Mash/CRBC_GTDB_closely_related.txt | cut -f 1` \
 ::: `cat $wd/05_novel_species/Mash/CRBC_GTDB_closely_related.txt | cut -f 2`

# 2. summary identity and coverage:
cd $wd/05_novel_species/mummer/temp/
time parallel -j 96 --xapply \
echo -e '{1} {2} `grep "AlignedBases" {1}_{2}.report` `grep "AvgIdentity" {1}_{2}.report | tail -1`' >> $wd/05_novel_species/mummer/CRBC_GTDB_ANI_AF_temp.txt \
::: `cat $wd/05_novel_species/Mash/CRBC_GTDB_closely_related.txt | cut -f 1` \
::: `cat $wd/05_novel_species/Mash/CRBC_GTDB_closely_related.txt | cut -f 2`

# Summarise the result, pick the lowest coverage:
sed 's/(/ /g;s/)/ /g;s/%//g;s/.fna//g;s/  / /g' $wd/05_novel_species/mummer/CRBC_GTDB_ANI_AF_temp.txt | cat -A | cut -f 1-2,5,7,9 -d ' ' | awk '{if($3>=$4)$6=$4;else{$6=$3};print $0}' | sed '1i GenomeID Reference Ref_cov Qry_cov Identity Coverage' | sed 's/ /\t/g' > $wd/05_novel_species/CRBC_GTDB_id_cov_info.txt
```
