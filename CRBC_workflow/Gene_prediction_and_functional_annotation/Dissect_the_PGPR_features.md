## Dissect the PGPR potential among CRBC and published genomes

### 1. Conduct genome function annotation using diamond

#### a. Build diamond database

```
## first build a kegg protein database including all eukaryotes, prokaryotes and virus, please be noted the kegg database we used is the April 1st 2021 version
    cd kegg/genes/fasta
    gunzip eukaryotes.pep.gz
    gunzip prokaryotes.pep.gz
    gunzip T40000.pep.gz
    cat eukaryotes.pep prokaryotes.pep T40000.pep > KEGG_eukaryotes_prokaryotes_virus_protein.faa

## extract peptide that encode PGPR functions, for detailed PGPR function defination please refer to the supplenmental table

### prepare all PGPR peptide sequences

    grep -w -F -f PGPR_KO_list ~/db/ftp_download_KEGG_04_01_2021/ko_gene_map.txt  | cut -f 1 > All_PGPR_KeggGeneID_list #153498 KeggGeneID
    seqkit grep -f All_PGPR_KeggGeneID_list  KEGG_eukaryotes_prokaryotes_virus_protein.faa > All_PGPR_all_organisms_KeggGene_peptide.faa #

### add CPS and KS protein into above (note: cause lots of bacteria are able to synthesize GA, however, kegg database did not include as much as GA synthesis related proteins, so we manually download and add into the PGPR protein reference)
    > CPS accession ID:NP_768789, NP_106893, NP_443949, NP_659791, WP_003466962, AEQ94336, WP_020322919
    > KS accession ID: NP_768790,NP_106894,NP_443948,NP_659792,WP_003466963,AEQ94335
    conda activate entrez_direct
    esearch -db protein -query "$(cat KS_accession_list | tr '\n' ' ')" | efetch -format fasta > KS_protein_sequences.fasta
    esearch -db protein -query "$(cat CPS_accession_list | tr '\n' ' ')" | efetch -format fasta > CPS_protein_sequences.fasta

### Add CPS and KS as prefix in the header line for both faa file and then combine together to construct diamond db
    cat CPS_protein_sequences.fasta  KS_protein_sequences.fasta > CPS_KS_protein_seq.faa

## make diamond database
    db=/mnt/m1/liufang/db
    cat All_PGPR_all_organisms_KeggGene_peptide.faa CPS_KS_protein_seq.faa  > Kegg_PGPR_plus_GA_CPS_KS.faa
    time \
    /mnt/m2/dairui/anaconda3/envs/diamond2.0.15/bin/diamond makedb \
    --in Kegg_PGPR_plus_GA_CPS_KS.faa \
    --db $db/PGPR_DIAMOND_2023/All_Kegg_PGPR_plus_pub_CPS_KS_dmnd \
    --threads 96 >> log/diamond_makedb.log 2>&1
    grep '^>' CPS_KS_protein_seq.faa  | cut -d ' ' -f 1 | sed 's,^>,,g' > Bacteria_CPS_KS_KO_gene_map.txt
```

#### b. Conduct blast

```
mkdir log
mkdir diamond_out
time \
    diamond blastp \
    --query CRBD_all_9772_prodigal_combine.faa \
    --db $db/PGPR_DIAMOND_2023/All_Kegg_PGPR_plus_pub_CPS_KS_dmnd.dmnd \
    --threads 96 \
    --out diamond_out/CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot.fm6 \
    --outfmt 6 \
    --log \
    --max-target-seqs 1 \
    --evalue 1e-5 \
    --sensitive \
    --block-size 4 \
    --index-chunks 1 \
    >> log/CRBD_all_9772_genomes_faa_vs_PGPR_dmndblastp.log 2>&1
```

#### c. predict phosphatase subcellular location (SCL) (since not all acid or alkaline phosphatase will be secreted to extracellular, so we predicted their SCL)

```
### prepare the query protein sequence

## extract the geneID for alkaline phosphatase
    grep -w -F -f <(grep -w -F -f ALP_KO_list  ~/db/ftp_download_KEGG_04_01_2021/ko_gene_map.txt | cut -f 1) All_genomes_14242_kegg_1.0e-5_id50plus.fm6 | cut -f 1 > kegg_fm6_extracted_alkaline_phosphatase_gene_ID
    cat /mnt/m2/dairui/project/binning/MAG_finalization/all/04_prodigal/sep_genome_faa/*.faa > All_14242_genomes_prodigal.faa
    grep -w -F -f kegg_fm6_extracted_alkaline_phosphatase_gene_ID  All_14242_genomes_prodigal.faa | cut -d ' ' -f 1   | sed 's,^>,,g' > All_14242_genomes_prodigal_alkaline_phosphatase_geneID
    wc -l All_14242_genomes_prodigal_alkaline_phosphatase_geneID # 28930
    wc -l kegg_fm6_extracted_alkaline_phosphatase_gene_ID # 29936

#### NOTE!!! some of the geneID listed in the prodigal file could not be matched with that from kegg file

    paste -d "\t" <(grep -v -w -F -f All_14242_genomes_prodigal_alkaline_phosphatase_geneID kegg_fm6_extracted_alkaline_phosphatase_gene_ID | sed 's,___,\t,g' | cut -f 1) <(grep -v -w -F -f All_14242_genomes_prodigal_alkaline_phosphatase_geneID kegg_fm6_extracted_alkaline_phosphatase_gene_ID) | sed 's,\t,___,g'  > ALP_lack_geneID_in_prodigal_output_vs_kegg_fm6
    grep -w -F -f ALP_lack_geneID_in_prodigal_output_vs_kegg_fm6  All_14242_genomes_prodigal.faa | wc -l
    cat All_14242_genomes_prodigal_alkaline_phosphatase_geneID ALP_lack_geneID_in_prodigal_output_vs_kegg_fm6 > All_14242_genomes_prodigal_alkaline_phosphatase_geneID_update
    grep -w -F -f All_14242_genomes_prodigal_alkaline_phosphatase_geneID_update All_14242_genomes_prodigal.faa | sed 's,^>,,g' > All_14242_genomes_prodigal_alkaline_phosphatase_gene_header_list
    #seqkit grep -n -f All_14242_genomes_prodigal_alkaline_phosphatase_gene_header_list --max-mismatch 0 All_14242_genomes_prodigal_alkaline_phosphatase_protein.faa > All_14242_genomes_prodigal_alkaline_phosphatase_protein.faa
    seqtk subseq All_14242_genomes_prodigal.faa All_14242_genomes_prodigal_alkaline_phosphatase_gene_header_list > All_14242_genomes_prodigal_alkaline_phosphatase_protein.faa


##### extract acid phosphatase, di and tri esterase and phytase enzymes

    grep -w -F -f <(grep -w -F -f NACP_ditriesterase_phytase_KO_list  ~/db/ftp_download_KEGG_04_01_2021/ko_gene_map.txt | cut -f 1) All_genomes_14242_kegg_1.0e-5_id50plus.fm6 | cut -f 1 > kegg_fm6_extracted_other_P_enzymes_gene_ID
    grep -w -F -f kegg_fm6_extracted_other_P_enzymes_gene_ID  All_14242_genomes_prodigal.faa | cut -d ' ' -f 1   | sed 's,^>,,g' > All_14242_genomes_prodigal_other_P_enzymes_geneID
    wc -l All_14242_genomes_prodigal_other_P_enzymes_geneID # 37410
    wc -l kegg_fm6_extracted_other_P_enzymes_gene_ID #39205
    paste -d "\t" <(grep -v -w -F -f All_14242_genomes_prodigal_other_P_enzymes_geneID kegg_fm6_extracted_other_P_enzymes_gene_ID | sed 's,___,\t,g' | cut -f 1) <(grep -v -w -F -f All_14242_genomes_prodigal_other_P_enzymes_geneID kegg_fm6_extracted_other_P_enzymes_gene_ID) | sed 's,\t,___,g'  > lack_geneID_in_prodigal_output_vs_kegg_fm6
    grep -w -F -f lack_geneID_in_prodigal_output_vs_kegg_fm6  All_14242_genomes_prodigal.faa | wc -l
    cat All_14242_genomes_prodigal_other_P_enzymes_geneID lack_geneID_in_prodigal_output_vs_kegg_fm6 > All_14242_genomes_prodigal_other_P_enzymes_geneID_update
    wc -l All_14242_genomes_prodigal_other_P_enzymes_geneID_update #

    grep -w -F -f All_14242_genomes_prodigal_other_P_enzymes_geneID_update All_14242_genomes_prodigal.faa | sed 's,^>,,g' > All_14242_genomes_prodigal_other_P_enzymes_gene_header_list
    #seqkit grep -n -f All_14242_genomes_prodigal_other_P_enzymes_gene_header_list --max-mismatch 0 All_14242_genomes_prodigal_other_P_enzymes_protein.faa > All_14242_genomes_prodigal_other_P_enzymes_protein.faa
    seqtk subseq All_14242_genomes_prodigal.faa All_14242_genomes_prodigal_other_P_enzymes_gene_header_list > All_14242_genomes_prodigal_other_P_enzymes_protein.faa

#### predict the subcellular location for phosphatase

##### for alkaline phosphatase subset
    mdkir psort_output
    /mnt/m1/liufang/software/docker_psort/psortb -i All_14242_genomes_prodigal_alkaline_phosphatase_protein.faa -r /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/kegg/Alkaline_phosphatase/docker_psort/psort_output -p -n -a

    ## in total there are 29936 sequences were annotated to be alkaline phosphatase
    cd psort_output
    grep -c 'SeqID' 20231013055245_psortb_gramneg.txt #29936
    ## To extract the sequence ID for those predicted to be extracellular by psortb
    /bin/rm Subtract_ECSVM_and_SCL_BLAST.txt
    for i in $(seq 1 29936)
    do
    	grep -m $i "SeqID" 20231013055245_psortb_gramneg.txt| tail -n 1 | tr '\n' '\t' >> Subtract_ECSVM_and_SCL_BLAST.txt
    	grep -A 4 -m $i "SeqID" 20231013055245_psortb_gramneg.txt | tail -n 1 | tr '\n' '\t' >> Subtract_ECSVM_and_SCL_BLAST.txt
    	grep -A 11 -m $i "SeqID" 20231013055245_psortb_gramneg.txt | tail -n 1  >> Subtract_ECSVM_and_SCL_BLAST.txt
    done

    ## Generate the header list for exctracellular protein sequence
    grep 'Extracellular' Subtract_ECSVM_and_SCL_BLAST.txt | cut -f 1 -d '#' | sed 's,SeqID: ,,g' | sed 's, ,,g' > ALP_extracellular_prodigal_faa_gene_list.txt
    cat <(grep -v '___Pub' ALP_extracellular_prodigal_faa_gene_list.txt) <(grep '___Pub' ALP_extracellular_prodigal_faa_gene_list.txt | sed 's,___,\t,g' | cut -f 2,3 | sed 's,\t,___,g') > ALP_extracellular_prodigal_faa_gene_list_up.txt

    ## Summary of ALP extracellular gene list and matchness with that from kegg fm6
    wc -l ALP_extracellular_prodigal_faa_gene_list_up.txt #5720
    grep -w -F -f ALP_extracellular_prodigal_faa_gene_list_up.txt  All_genomes_14242_kegg_1.0e-5_id50plus.fm6 | wc -l #5720
    grep -v -w -F -f <(cut -f 2 /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/meta/Ray_May8_all_14242_genomes_full_metadata.txt)   <(cut -f 1 ALP_extracellular_all_14242_genomes_kegg_1.0e-5_id50plus.fm6 | sed 's,___,\t,g' | cut -f 1| sort | uniq)  # all contigs in the ALP fm6 file could be split into GenomeID and contig

    #!/bin/Rscript

    library(dplyr)
    library(reshape2)
    library(tidyr)
    library(stringr)

    ALP_extracellular_fm6<-read.table("ALP_extracellular_all_14242_genomes_kegg_1.0e-5_id50plus.fm6",sep="\t",header=FALSE)
    ALP_extracellular_fm6<-ALP_extracellular_fm6%>%select(c(V1,V2))
    colnames(ALP_extracellular_fm6)<-c("GeneID","KeggGeneID")
    ko_gene_map<-read.table("ko_gene_map.txt",sep="\t",header=TRUE)
    ALP_extracellular_fm6_addKO<-inner_join(ALP_extracellular_fm6,ko_gene_map)
    ALP_extracellular_fm6_addKO$GenomeID<-str_split(ALP_extracellular_fm6_addKO$GeneID,pattern="___",simplify=TRUE,n=2)[,1]
    ALP_extracellular_KO<-ALP_extracellular_fm6_addKO%>%select(c(KO,GenomeID))%>%unique()%>%mutate(Count=1)
    ALP_extracellular_KO_count<-ALP_extracellular_KO%>%spread(key=KO,value=Count)
    ALP_extracellular_KO_count[is.na(ALP_extracellular_KO_count)]=0
    write.table(ALP_extracellular_KO_count,"All_14242_genomes_ALP_extracellular_KO_count.txt",sep="\t",quote=FALSE,row.names=FALSE)

##### for other phosphatase subset
    mkdir psort_output_v2
    > NOTE < please be noted the term of "other_P_enzymes" includes the KOs belongs to non-specific acid phosphatse, Phytase, di- and tri- estrerase, such as K03788,K09474,K01078,K01083,K01093,K01126,K07048. ALP means alkaline phosphatase including K01077,K01113,K07093
    /mnt/m1/liufang/software/docker_psort/psortb -i All_14242_genomes_prodigal_other_P_enzymes_protein.faa -r /mnt/m3/liufang/NRgenome_CRBC_NCBI_IMG_ENA/kegg/Alkaline_phosphatase/docker_psort/psort_output_v2 -p -n -a
    cd psort_output_v2
    ## in total there are 29936 sequences were annotated to be other phosphatase
    grep -c 'SeqID' 20231013075650_psortb_gramneg.txt # 39205

    ## To extract the sequence ID for those predicted to be extracellular by psortb
    /bin/rm Subtract_ECSVM_and_SCL_BLAST.txt
    for i in $(seq 1 39205)
    do
    	grep -m $i "SeqID" 20231013075650_psortb_gramneg.txt| tail -n 1 | tr '\n' '\t' >> Subtract_ECSVM_and_SCL_BLAST.txt
    	grep -A 4 -m $i "SeqID" 20231013075650_psortb_gramneg.txt | tail -n 1 | tr '\n' '\t' >> Subtract_ECSVM_and_SCL_BLAST.txt
    	grep -A 11 -m $i "SeqID" 20231013075650_psortb_gramneg.txt | tail -n 1  >> Subtract_ECSVM_and_SCL_BLAST.txt
    done

    ## Generate the header list for exctracellular protein sequence
    grep 'Extracellular' Subtract_ECSVM_and_SCL_BLAST.txt | cut -f 1 -d '#' | sed 's,SeqID: ,,g' | sed 's, ,,g' > other_P_enzymes_extracellular_prodigal_faa_gene_list.txt
    cat <(grep -v '___Pub' other_P_enzymes_extracellular_prodigal_faa_gene_list.txt) <(grep '___Pub' other_P_enzymes_extracellular_prodigal_faa_gene_list.txt | sed 's,___,\t,g' | cut -f 2,3 | sed 's,\t,___,g') > other_P_enzymes_extracellular_prodigal_faa_gene_list_up.txt

    ## Summary of other_P_enzymes extracellular gene list and matchness with that from kegg fm6
    wc -l other_P_enzymes_extracellular_prodigal_faa_gene_list_up.txt # 2161
    grep -w -F -f other_P_enzymes_extracellular_prodigal_faa_gene_list_up.txt  All_genomes_14242_kegg_1.0e-5_id50plus.fm6 | wc -l # 2161
    grep -w -F -f other_P_enzymes_extracellular_prodigal_faa_gene_list_up.txt  All_genomes_14242_kegg_1.0e-5_id50plus.fm6 > other_P_enzymes_extracellular_all_14242_genomes_kegg_1.0e-5_id50plus.fm6

    #!/bin/Rscript

    library(dplyr)
    library(reshape2)
    library(tidyr)
    library(stringr)

    other_P_enzymes_extracellular_fm6<-read.table("other_P_enzymes_extracellular_all_14242_genomes_kegg_1.0e-5_id50plus.fm6",sep="\t",header=FALSE)
    other_P_enzymes_extracellular_fm6<-other_P_enzymes_extracellular_fm6%>%select(c(V1,V2))
    colnames(other_P_enzymes_extracellular_fm6)<-c("GeneID","KeggGeneID")
    ko_gene_map<-read.table("$db/ftp_download_KEGG_04_01_2021/ko_gene_map.txt",sep="\t",header=TRUE)
    other_P_enzymes_extracellular_fm6_addKO<-inner_join(other_P_enzymes_extracellular_fm6,ko_gene_map)
    other_P_enzymes_extracellular_fm6_addKO$GenomeID<-str_split(other_P_enzymes_extracellular_fm6_addKO$GeneID,pattern="___",simplify=TRUE,n=2)[,1]
    other_P_enzymes_extracellular_KO<-other_P_enzymes_extracellular_fm6_addKO%>%select(c(KO,GenomeID))%>%unique()%>%mutate(Count=1)
    other_P_enzymes_extracellular_KO_count<-other_P_enzymes_extracellular_KO%>%spread(key=KO,value=Count)
    other_P_enzymes_extracellular_KO_count[is.na(other_P_enzymes_extracellular_KO_count)]=0
    write.table(other_P_enzymes_extracellular_KO_count,"All_14242_genomes_other_P_enzymes_extracellular_KO_count.txt",sep="\t",quote=FALSE,row.names=FALSE)

##### Conbine ALP and Other_P_enzyme results together
    mkdir phosphorus_SCL_combine
    cd phosphorus_SCL_combine
    #!/bin/Rscript
    library(dplyr)
    library(reshape2)
    library(stringr)
    library(tidyr)

    ALP_KO_extracellular<-read.table("All_14242_genomes_ALP_extracellular_KO_count.txt",sep="\t",header=TRUE,quote="")
    other_P_enzymes_extracellular<-read.table("All_14242_genomes_other_P_enzymes_extracellular_KO_count.txt",sep="\t",header=TRUE,quote="")
    CRBD_9772_genomes_meta<-read.table("9772_CRBD_genomes_metadata.txt",sep="\t",header=TRUE,quote="")
    ALP_NSAP_phytase_esterase_extracellular_KO<-full_join(ALP_KO_extracellular,other_P_enzymes_extracellular)
    CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO<-left_join(CRBD_9772_genomes_meta%>%select(GenomeID),ALP_NSAP_phytase_esterase_extracellular_KO)
    CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO[is.na(CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO)]=0
    CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO%>%mutate(SUM=K01077+K01113+K07093+K01083+K01126)%>%filter(SUM>0)%>%dim() # 4164
    write.table(CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO,"CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO.txt",sep="\t",row.names=FALSE,quote=FALSE)
```

## d. phosphatase entries were updated based on subcellular location prediction results

```
    cat ko_gene_map.txt Bacteria_CPS_KS_KO_gene_map.txt > ko_gene_map_add_pub_CPS_KS.txt
    paste -d "\t" <(sed 's,.*___Pub,Pub,g' diamond_out/CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot.fm6 | cut -f 1 | sed 's,___,\t,g' | cut -f 1) <(sed 's,.*___Pub,Pub,g' diamond_out/CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot.fm6) > CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot_up.fm6

    #!/bin/Rscript

    library(dplyr)
    library(reshape2)
    library(stringr)

    PGPR<-read.table("CRBD_all_9772_genome_vs_PGPR_Dec20_dmndblastp_annot_up.fm6",sep="\t",header=FALSE,quote="")
    # ---------- PGPR KO count --- evalue<1.0e-5 ------------

    PGPR_up_evalue5<-PGPR%>%select(c(V1,V2,V3))%>%unique()
    dim(PGPR_up_evalue5) # 4205703 rows
    colnames(PGPR_up_evalue5)<-c("NRgenomeID","GeneID","KeggGeneID")
    ko_gene_map<-read.table("ko_gene_map_add_pub_CPS_KS.txt",sep="\t",header=TRUE,quote="") # cut -f 1 ko_gene_map.txt | sort | uniq -d | wc -l # 17060 of the KeggGeneId belongs to several different KOs, 16975039 rows in total
    ko_gene_map<-unique(ko_gene_map) # 16975039 unique rows
    summary(PGPR_up_evalue5$KeggGeneID%in%ko_gene_map$KeggGeneID)
    #   Mode    TRUE
    #logical 4205703
    PGPR_up_evalue5_addKO<-inner_join(PGPR_up_evalue5,ko_gene_map)
    dim(PGPR_up_evalue5_addKO) #4204713 rows and 4 columns, some of teh KeggGeneID lack of KO and some are assigned to more than one KO
    PGPR_up_evalue5_KO_count<-PGPR_up_evalue5_addKO%>%group_by(NRgenomeID,KO)%>%mutate(KO_count=n())%>%select(c(NRgenomeID,KO,KO_count))%>%unique()
    library(tidyr)
    PGPR_up_evalue5_KO_count_wide<-data.frame(PGPR_up_evalue5_KO_count%>%spread(key=NRgenomeID,value=KO_count),check.names=FALSE)
    PGPR_up_evalue5_KO_count_wide[is.na(PGPR_up_evalue5_KO_count_wide)]=0
    write.table(PGPR_up_evalue5_KO_count_wide,"CRBD_all_9772_genomes_PGPR_Dec20_evalue_5_KO_count.txt",sep="\t",row.names=FALSE,quote=FALSE)
    summary(data.frame(apply(PGPR_up_evalue5_KO_count_wide%>%select(!KO),2,function(x) sum(x>0))))
    CRBD_HQ_NRgenome_meta<-read.table("6109_CRBD_HQ_NRgenome_meta.txt",sep="\t",header=TRUE)
    CRBD_HQ_NRgenome_PGPR_evalue5_KO_count<-PGPR_up_evalue5_KO_count_wide%>%select(c(colnames(PGPR_up_evalue5_KO_count_wide)[colnames(PGPR_up_evalue5_KO_count_wide)%in%CRBD_HQ_NRgenome_meta$GenomeID]))
    rownames(CRBD_HQ_NRgenome_PGPR_evalue5_KO_count)<-as.character(PGPR_up_evalue5_KO_count_wide$KO)
    CRBD_HQ_NRgenome_PGPR_KO_SUM<-data.frame(PGPR_KO_SUM=apply(CRBD_HQ_NRgenome_PGPR_evalue5_KO_count,2,function(x) sum(x>0)))
    summary(CRBD_HQ_NRgenome_PGPR_KO_SUM)
    CRBD_HQ_NRgenome_PGPR_KO_SUM_meta<-inner_join(data.frame(CRBD_HQ_NRgenome_PGPR_KO_SUM,GenomeID=rownames(CRBD_HQ_NRgenome_PGPR_KO_SUM)),CRBD_HQ_NRgenome_meta%>%select(c(GenomeID,Kingdom..gtdb.,Phylum..gtdb.,Class..gtdb.,Order..gtdb.,Family..gtdb.,Genus..gtdb.,Species..gtdb.)))
    CRBD_HQ_NRgenome_PGPR_KO_SUM_meta%>%filter(PGPR_KO_SUM<=60)%>%group_by(Phylum..gtdb.)%>%mutate(temp=n())%>%select(c(Phylum..gtdb.,temp))%>%unique()
    # ---------- PGPR KO count --- evalue<1.0e-5 and identity>=50%-----------

    PGPR_up_evalue5_id50plus<-PGPR%>%filter(V4>=50)%>%select(c(V1,V2,V3))%>%unique()
    dim(PGPR_up_evalue5_id50plus) # 407810       3
    colnames(PGPR_up_evalue5_id50plus)<-c("NRgenomeID","GeneID","KeggGeneID")
    PGPR_up_evalue5_id50plus_addKO<-inner_join(PGPR_up_evalue5_id50plus,ko_gene_map)
    dim(PGPR_up_evalue5_id50plus_addKO) # 406677 rows and 4 columns,
    PGPR_up_evalue5_id50plus_KO_count<-PGPR_up_evalue5_id50plus_addKO%>%group_by(NRgenomeID,KO)%>%mutate(KO_count=n())%>%select(c(NRgenomeID,KO,KO_count))%>%unique()
    library(tidyr)
    PGPR_up_evalue5_id50plus_KO_count_wide<-data.frame(PGPR_up_evalue5_id50plus_KO_count%>%spread(key=NRgenomeID,value=KO_count))
    PGPR_up_evalue5_id50plus_KO_count_wide[is.na(PGPR_up_evalue5_id50plus_KO_count_wide)]=0
    unique(PGPR$V1)[!unique(PGPR$V1)%in%colnames(PGPR_up_evalue5_id50plus_KO_count_wide%>%select(!KO))]
    #NOTE!!! the results indicate that the above 12 genomes does not has any of the PGPR KO when filter based on evalue_5&id50plus
    write.table(PGPR_up_evalue5_id50plus_KO_count_wide,"CRBD_all_9772_genomes_PGPR_Dec20_evalue5_id50plus_KO_count.txt",sep="\t",row.names=FALSE,quote=FALSE)

    #!/bin/shell
    grep -w -v -F -f Phosphatase_SCL_results_link/KO_list_used_for_SCL_prediction CRBD_all_9772_genomes_PGPR_Dec20_evalue5_id50plus_KO_count.txt > CRBD_all_9772_genomes_PGPR_Dec20_exclude_phosphorus_SCL_evalue5_id50plus_KO_count.txt
    grep -w -F -f Phosphatase_SCL_results_link/KO_list_used_for_SCL_prediction CRBD_all_9772_genomes_PGPR_Dec20_evalue5_id50plus_KO_count.txt  | cut -f 1 # those below 9 KOs are removed since their SCL KO listed were updated based on psortb results

    ## now, we have both ``CRBD_all_9772_genomes_PGPR_Dec20_exclude_phosphorus_SCL_evalue5_id50plus_KO_count.txt`` and  ``CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO.txt`` ready, combine together we will analyze the PGPR distribution among CRBD
```
