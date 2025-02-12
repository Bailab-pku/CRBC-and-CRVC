## Rscript for main figure - The multifacet PGPR potential encoded within CRBD 
#----------------------------------
knitr::opts_chunk$set(echo = TRUE)

# load envs and setup main_theme

library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(vegan)
library(phyloseq)
library(ggpubr)
library(pheatmap)
library(DESeq2)
library(directlabels)
library(tidyr)
library(extrafont)
library(circlize)
#font_import(pattern = "Arial") # this only needed once and no need to import every time we run loadfonts(device='pdf')
loadfonts(device = "pdf")
library(RColorBrewer)
library(tidyheatmap)
library(cowplot)
library(stringr)
library(patchwork)
library(ggsci)
main_theme<-theme(panel.background = element_rect(fill = "transparent",colour = "black"),plot.background = element_rect(fill = "transparent", color = NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(family="Arial"), plot.title = element_text(family='Arial', hjust = .5,vjust = 3,size = 15,face="bold"), axis.text.x = element_text(family="Arial",size=10,angle = 90,hjust = 0,vjust = 0),axis.text.y = element_text(family="Arial",size=10),axis.title.x = element_text(family="Arial",size=12,face="bold"),axis.title.y = element_text(family="Arial",size=12,face="bold"), legend.title = element_text(family="Arial",size=12,face="bold"),legend.text =element_text(family="Arial",size=10),legend.position = "right",legend.background = element_rect(fill="transparent"),legend.box.background = element_rect(color="transparent", fill = "transparent"),legend.key = element_rect(fill = "transparent", colour = NA),strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),plot.margin=unit(c(1,1,1,1),"cm"))

taxa_color<-read.csv("/Users/fangliu/Documents/IGDB_Bai_lab/Database_and_color_pallete/Color_scheme/GTDB_phylum_color_Ray_define_top12_Sep_05_2023.csv",header = FALSE,col.names = c("PhyClass","Color"))
rownames(taxa_color)<-taxa_color$PhyClass

family_color_meta<-read.csv("/Users/fangliu/Documents/IGDB_Bai_lab/Database_and_color_pallete/Color_scheme/Family_color_meta.csv",header = TRUE)
rownames(family_color_meta)<-family_color_meta$Family


# Subset CRBD HQ NRgenome taxonomy


CRBD_genome_meta<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/9772_CRBD_genomes_metadata_0918.txt",sep="\t",header=TRUE,quote = "") #9772
summary(CRBD_genome_meta$GenomeID=="Pub_GCA_000468095.1_Pantoea_sp._AS-PWVM4_genomic") # TRUE=1
colnames(CRBD_genome_meta)<-gsub(pattern = "[.]",replacement = "_",colnames(CRBD_genome_meta))
colnames(CRBD_genome_meta)<-gsub(pattern = "__ncbi_","_ncbi",colnames(CRBD_genome_meta))
colnames(CRBD_genome_meta)<-gsub(pattern = "__gtdb_","_gtdb",colnames(CRBD_genome_meta))
CRBD_genome_meta_up<-inner_join(CRBD_genome_meta,CRBD_genome_meta%>%filter(Rep_NR_Genomes=="Yes")%>%select(GenomeID,Genome_ClusterID),by="Genome_ClusterID")%>%dplyr::rename(NRgenome_GenomeID=GenomeID.y)%>%dplyr::rename(GenomeID=GenomeID.x) # 9772
CRBD_genome_meta_up2<-inner_join(CRBD_genome_meta_up,CRBD_genome_meta_up%>%filter(Rep_NR_Species=="Yes")%>%select(GenomeID,Species_ClusterID),by="Species_ClusterID")%>%dplyr::rename(RepSpecies_GenomeID=GenomeID.y)%>%dplyr::rename(GenomeID=GenomeID.x) #9772
CRBD_NRgenome<-CRBD_genome_meta_up2%>%filter(Rep_NR_Genomes=="Yes") #7531 CRBD NRgenome
CRBD_HQ_NRgenome<-CRBD_NRgenome%>%filter(Quality_level=="High Quality") #6109 CRBD HQ NRgenome
# length(unique((CRBD_genome_meta_up2%>%filter(Quality_level=="High Quality"))$NRgenome_GenomeID)) #[1] 6111
ID1<-(CRBD_genome_meta_up2%>%filter(Quality_level=="High Quality")%>%filter(Rep_NR_Genomes=="Yes"))$NRgenome_GenomeID
ID2<-unique((CRBD_genome_meta_up2%>%filter(Quality_level=="High Quality"))$NRgenome_GenomeID)
ID2[!ID2%in%ID1] # this conflict is caused by the fact that the NRgenomeID is defined by the Quality score, but some times the Rep NRgenomeID has higher quality score but not met with the High Quality standard. So, defined always based on Rep_NR_Genomes and then filter to High Quality. 
rownames(CRBD_HQ_NRgenome)<-CRBD_HQ_NRgenome$NRgenome_GenomeID
CRBD_RepSpecies<-CRBD_genome_meta_up2%>%filter(Rep_NR_Species=="Yes") # 3044 CRBD RepSpecies
CRBD_HQ_RepSpecies<-CRBD_RepSpecies%>%filter(Quality_level=="High Quality") # 2659 CRBD HQ RepSpecies
row.names(CRBD_HQ_RepSpecies)<-CRBD_HQ_RepSpecies$RepSpecies_GenomeID

uniq_genus_CRBD_HQ_NRgenome<-unique(CRBD_HQ_NRgenome$Genus_gtdb)
uniq_genus_CRBD_HQ_RepSpecies<-unique(CRBD_HQ_RepSpecies$Genus_gtdb)
uniq_genus_CRBD_HQ_NRgenome[!uniq_genus_CRBD_HQ_NRgenome%in%uniq_genus_CRBD_HQ_RepSpecies]


CRBD_NRgenome_taxonomy<-CRBD_NRgenome%>%select(c(NRgenome_GenomeID,Phylum_gtdb,Class_gtdb,Order_gtdb,Family_gtdb,Genus_gtdb,Species_gtdb))
colnames(CRBD_NRgenome_taxonomy)<-c("NRgenome_GenomeID","Phylum","Class","Order","Family","Genus","Species")
CRBD_NRgenome_taxonomy<-CRBD_NRgenome_taxonomy%>%mutate(PhyClass=if_else(Phylum=="p__Proteobacteria",Class,Phylum))
CRBD_NRgenome_taxonomy$PhyClass<-gsub("p__","",CRBD_NRgenome_taxonomy$PhyClass)
CRBD_NRgenome_taxonomy$PhyClass<-gsub("c__","",CRBD_NRgenome_taxonomy$PhyClass)
CRBD_NRgenome_taxonomy$Genus<-gsub("g__","",CRBD_NRgenome_taxonomy$Genus)

CRBD_HQ_NRgenome_taxonomy<-CRBD_HQ_NRgenome%>%select(c(NRgenome_GenomeID,Phylum_gtdb,Class_gtdb,Order_gtdb,Family_gtdb,Genus_gtdb,Species_gtdb))
colnames(CRBD_HQ_NRgenome_taxonomy)<-c("NRgenome_GenomeID","Phylum","Class","Order","Family","Genus","Species")
CRBD_HQ_NRgenome_taxonomy<-CRBD_HQ_NRgenome_taxonomy%>%mutate(PhyClass=if_else(Phylum=="p__Proteobacteria",Class,Phylum))
CRBD_HQ_NRgenome_taxonomy$PhyClass<-gsub("p__","",CRBD_HQ_NRgenome_taxonomy$PhyClass)
CRBD_HQ_NRgenome_taxonomy$PhyClass<-gsub("c__","",CRBD_HQ_NRgenome_taxonomy$PhyClass)
CRBD_HQ_NRgenome_taxonomy$Genus<-gsub("g__","",CRBD_HQ_NRgenome_taxonomy$Genus)
CRBD_HQ_NRgenome_taxonomy_PhyClass_size<-data.frame(CRBD_HQ_NRgenome_taxonomy%>%group_by(PhyClass)%>%mutate(PhyClass_size=n())%>%select(c(PhyClass,PhyClass_size))%>%unique()%>%arrange(desc(PhyClass_size))) # the order is not consistent with Ray's definition, but in order to keep the figure legend consistent, will define the dominant PhyClass just the same as Ray's

PhyClass_dominant=c("Alphaproteobacteria","Gammaproteobacteria","Actinobacteriota","Firmicutes","Bacteroidota","Myxococcota","Patescibacteria","Spirochaetota","Fibrobacterota","Chloroflexota","Verrucomicrobiota","Acidobacteriota") # Patescibacteria does not have HQ NRgenome
CRBD_HQ_NRgenome_taxonomy<-CRBD_HQ_NRgenome_taxonomy%>%mutate(PhyClass_collapse=if_else(PhyClass%in%PhyClass_dominant,PhyClass,"Others"))
PhyClass_order<-c(PhyClass_dominant,"Others")
write.table(CRBD_HQ_NRgenome_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/Meta/CRBD_HQ_NRgenome_taxanomy_gtdb.txt",sep = "\t",quote=FALSE,row.names = FALSE)
rownames(CRBD_HQ_NRgenome_taxonomy)<-CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID # 6109 HQ RepSpecies
CRBD_HQ_NRgenome_taxonomy_PhyClass_size<-data.frame(CRBD_HQ_NRgenome_taxonomy%>%group_by(PhyClass_collapse)%>%summarise(PhyClass_size=n()))%>%arrange(desc(PhyClass_size))
HQ_NRgenome_taxonomy_up_sort<-inner_join(CRBD_HQ_NRgenome_taxonomy,CRBD_HQ_NRgenome_taxonomy_PhyClass_size)%>%arrange(desc(PhyClass_size))
rownames(HQ_NRgenome_taxonomy_up_sort)<-HQ_NRgenome_taxonomy_up_sort$NRgenome_GenomeID
HQ_NRgenome_taxonomy_up_sort$PhyClass_collapse<-factor(HQ_NRgenome_taxonomy_up_sort$PhyClass_collapse,levels = PhyClass_order[PhyClass_order!="Patescibacteria"]) #6109

PhyClass_col<-list(Color=taxa_color[c(1:12,15),]$Color)
names(PhyClass_col$Color)<-taxa_color[c(1:12,15),]$PhyClass

# CRBD HQ NRgenome KOs


## ----- combine input together -------
CRBD_all_9772_genomes_PGPR_minus_PSCL_KO_evalue5id50<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/PGPR_KO_count/CRBD_all_9772_genomes_PGPR_Dec20_exclude_phosphorus_SCL_evalue5_id50plus_KO_count.txt",sep="\t",header = TRUE,quote = "",row.names = 1,check.names = FALSE) # 126 PGPR KOs across 9761 genomes, PSCL means phosphorus related KOs that are predicted subcellular location
summary(rownames(CRBD_all_9772_genomes_PGPR_minus_PSCL_KO_evalue5id50)=="K13035") # none of them are 13035, which is nitralase but not contributed to IAA in out dataset
colnames(CRBD_all_9772_genomes_PGPR_minus_PSCL_KO_evalue5id50)<-gsub("Pub_GCA_000468095.1_Pantoea_sp._AS.PWVM4_genomic","Pub_GCA_000468095.1_Pantoea_sp._AS-PWVM4_genomic",colnames(CRBD_all_9772_genomes_PGPR_minus_PSCL_KO_evalue5id50))
summary(colnames(CRBD_all_9772_genomes_PGPR_minus_PSCL_KO_evalue5id50)%in%CRBD_genome_meta$GenomeID) # All is TRUE
#CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50<-CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50%>%filter(!KO%in%c("K01113","K07093","K03788","K09474","K01078","K01083","K01093","K01126","K07048"))%>%select(!KO) # 125 KOs, K01078 not exist, alkaline, acid, phytase, diesterase, triesterase, CPS(K0412) and KS (K04121) are absent.
CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50<-CRBD_all_9772_genomes_PGPR_minus_PSCL_KO_evalue5id50%>%select(c(CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID)) # 126 KOs across 6109 genomes
summary(colnames(CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50)%in%CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID) # 6109 all TRUE
CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50$KO<-row.names(CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50)

## Read into P related extracellular enzymes

P_extracellular_KO_count<-read.table("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/PGPR_KO_count/CRBD_9772_genomes_ALP_NSAP_phytase_esterase_extracellular_KO.txt",sep="\t",header=TRUE,row.names = 1)
P_extracellular_KO_count_t<-data.frame(t(P_extracellular_KO_count),KO=colnames(P_extracellular_KO_count),check.names = FALSE)[,colnames(CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50)] # 5 KOs
identical(colnames(P_extracellular_KO_count_t),colnames(CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50)) # TRUE
CRBC_HQ_NRgenome_PGPR_KO_evalue5id50<-data.frame(rbind(CRBC_HQ_NRgenome_PGPR_minus_PSCL_KO_evalue5id50,P_extracellular_KO_count_t),check.names = FALSE) # now 131 KOs 


## ----- annotation files ---

all_PGPR_KO_annot<-read.csv("/Users/fangliu/Documents/IGDB_Bai_lab/Database_and_color_pallete/KEGG/FL_version/PGPR_final_version_Jan28/updated_keep_core/update_0911/All_PGPR_combine_1007.csv",header = TRUE,sep = ',') #172 + phnCED = 175
all_PGPR_KO_annot<-all_PGPR_KO_annot%>%filter(KO!="")
all_PGPR_KO_annot$KO[duplicated(all_PGPR_KO_annot$KO)] # Six KOs involved in SA biosynthesis are duplicated since those genes are shared for SA and siderophore biosynthesis, "K01851" "K14759" "K02361" "K02552" "K04782" "K04781"
growth_PGPR_KO_annot<-all_PGPR_KO_annot%>%filter(PGP=="Growth")
growth_PGPR_KO_annot$Group<-factor(growth_PGPR_KO_annot$Group,levels = c("IAA","GA","CK"))
growth_PGPR_KO_annot$KO<-factor(growth_PGPR_KO_annot$KO,levels = as.character(growth_PGPR_KO_annot$KO))
rownames(growth_PGPR_KO_annot)<-growth_PGPR_KO_annot$KO
growth_PGPR_KO_annot<-growth_PGPR_KO_annot%>%select(!Details)

nutrient_PGPR_KO_annot<-all_PGPR_KO_annot%>%filter(PGP=="Nutrient")
nutrient_PGPR_KO_annot$Group<-factor(nutrient_PGPR_KO_annot$Group,levels = c("Nitrogen","Phosphorus","Iron"))
nutrient_PGPR_KO_annot$KO<-factor(nutrient_PGPR_KO_annot$KO,levels = as.character(nutrient_PGPR_KO_annot$KO))
rownames(nutrient_PGPR_KO_annot)<-nutrient_PGPR_KO_annot$KO
nutrient_PGPR_KO_annot<-nutrient_PGPR_KO_annot%>%select(!Details)

AdaptationR_KO_annot<-all_PGPR_KO_annot%>%filter(PGP=="Adaptation")
AdaptationR_KO_annot$Group<-factor(AdaptationR_KO_annot$Group,levels = c("ACC_deaminase","ET","SA","ABA","JA","Antibiotics","Cyclic_lipopeptides","VOC"))
AdaptationR_KO_annot$KO<-factor(AdaptationR_KO_annot$KO,levels = as.character(AdaptationR_KO_annot$KO))
rownames(AdaptationR_KO_annot)<-AdaptationR_KO_annot$KO
AdaptationR_KO_annot<-AdaptationR_KO_annot%>%select(!Details)

## ------- heatmap - cluster within phylum ---

CRBD_HQ_NRgenoem_all_PGPR_KO<-CRBC_HQ_NRgenome_PGPR_KO_evalue5id50[rownames(CRBC_HQ_NRgenome_PGPR_KO_evalue5id50)[rownames(CRBC_HQ_NRgenome_PGPR_KO_evalue5id50)%in%unique(all_PGPR_KO_annot$KO)],]%>%select(!KO) # 131 PGPR KOs across 6109 NRgenomes
CRBD_HQ_NRgenoem_all_PGPR_KO<-CRBD_HQ_NRgenoem_all_PGPR_KO[apply(CRBD_HQ_NRgenoem_all_PGPR_KO,1,sum)>0,] # 128 none zero after subset to CRBD HQ NRgenome, "K10760=IPT" "K12237=NRPS" "K24425=SDR9C7"

####------ KO annotation
all_PGPR_KO_annot_up<-all_PGPR_KO_annot%>%filter(KO%in%rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)) #134 including six duplicated KOs
all_PGPR_KO_annot_up<-all_PGPR_KO_annot_up%>%filter(Class!="SA_shared") # 128
rownames(all_PGPR_KO_annot_up)<-all_PGPR_KO_annot_up$KO
all_PGPR_KO_group<-all_PGPR_KO_annot_up[rownames(all_PGPR_KO_annot_up)[rownames(all_PGPR_KO_annot_up)%in%rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)],]%>%select(Group) # 128
all_PGPR_KO_group_list<-as.character(all_PGPR_KO_group$Group)
all_PGPR_col = list(Color=as.character(unique(all_PGPR_KO_annot_up$Color)))
names(all_PGPR_col$Color)<-as.character(unique(all_PGPR_KO_annot_up$Group))

####------ NRgenome annotation
GenomeID_annot<-HQ_NRgenome_taxonomy_up_sort[rownames(HQ_NRgenome_taxonomy_up_sort)[rownames(HQ_NRgenome_taxonomy_up_sort)%in%colnames(CRBD_HQ_NRgenoem_all_PGPR_KO)],]%>%select(PhyClass_collapse)
PhyClass = as.character(GenomeID_annot$PhyClass_collapse) #6109

library(ComplexHeatmap)
set.seed(1013)
CRBD_HQ_NRgenoem_all_PGPR_KO<-CRBD_HQ_NRgenoem_all_PGPR_KO[rownames(all_PGPR_KO_group),rownames(GenomeID_annot)] 
CRBD_HQ_NRgenoem_all_PGPR_KO[CRBD_HQ_NRgenoem_all_PGPR_KO>0]=1
write.csv(CRBD_HQ_NRgenoem_all_PGPR_KO,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_All_PGPR_CRBD_HQ_NRgenome_KO_binary.csv",sep = ',',row.names = TRUE,quote = FALSE)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/CRBD_HQ_NRgenome_all_PGPR_phylogenetic_distribution.pdf",width = 15,height = 10)
m = as.matrix(CRBD_HQ_NRgenoem_all_PGPR_KO)
dend1 = cluster_within_group(m,PhyClass)
Heatmap(m, cluster_columns = dend1, column_split = 12,col = colorspace::sequential_hcl(n=50,palette = "Dark Mint",rev = TRUE),
    row_title = "cluster_within_PhyClass",row_title_side = 'left',row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize=2),top_annotation = HeatmapAnnotation(PhyClass= PhyClass, col = list(PhyClass=PhyClass_col$Color)),left_annotation= rowAnnotation(all_PGPR_annot=all_PGPR_KO_group_list,col=list(all_PGPR_annot=all_PGPR_col$Color)),cluster_rows = FALSE,use_raster=FALSE)
dev.off()

## -------- nutrient only ---------

CRBD_HQ_NRgenome_nutrient_KO<-CRBD_HQ_NRgenoem_all_PGPR_KO[rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)[rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)%in%(all_PGPR_KO_annot%>%filter(PGP=="Nutrient"))$KO],] # 56 KOs across 6109 genomes
nutrient_PGPR_KO_annot_up<-data.frame(all_PGPR_KO_annot%>%filter(PGP=="Nutrient"))
rownames(nutrient_PGPR_KO_annot_up)<-nutrient_PGPR_KO_annot_up$KO
nutrient_PGPR_KO_annot_up<-nutrient_PGPR_KO_annot_up[rownames(CRBD_HQ_NRgenome_nutrient_KO),]
nutrient_PGPR_KO_group<-nutrient_PGPR_KO_annot_up[rownames(CRBD_HQ_NRgenome_nutrient_KO),]%>%select(Group) #56
nutrient_PGPR_KO_group_list<-as.character(nutrient_PGPR_KO_group$Group)

nutrient_PGPR_col = list(Color=unique(nutrient_PGPR_KO_annot_up$Color))
names(nutrient_PGPR_col$Color)<-as.character(unique(nutrient_PGPR_KO_annot_up$Group))

GenomeID_annot<-HQ_NRgenome_taxonomy_up_sort[rownames(HQ_NRgenome_taxonomy_up_sort)[rownames(HQ_NRgenome_taxonomy_up_sort)%in%colnames(CRBD_HQ_NRgenome_nutrient_KO)],]%>%select(PhyClass_collapse)
PhyClass = as.character(GenomeID_annot$PhyClass_collapse) #6109

library(ComplexHeatmap)
set.seed(1013)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/Nutrient_PGPR_KO_binary_phylogenetic_distribution.pdf",width = 15,height = 8)
m = as.matrix(CRBD_HQ_NRgenome_nutrient_KO)
dend1 = cluster_within_group(m,PhyClass)
Heatmap(m, cluster_columns = dend1, column_split = 12,col = colorspace::sequential_hcl(n=50,palette = "Dark Mint",rev = TRUE),
    row_title = "cluster_within_PhyClass",row_title_side = 'left',row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize=2),top_annotation = HeatmapAnnotation(PhyClass= PhyClass, col = list(PhyClass=PhyClass_col$Color)),left_annotation= rowAnnotation(nutrient_PGPR=nutrient_PGPR_KO_group_list,col=list(nutrient_PGPR=nutrient_PGPR_col$Color)),cluster_rows = FALSE,use_raster=FALSE)
dev.off()

## -------- Growth only ---------

CRBD_HQ_NRgenome_growth_KO<-CRBD_HQ_NRgenoem_all_PGPR_KO[rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)[rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)%in%(all_PGPR_KO_annot%>%filter(PGP=="Growth"))$KO],] # 34 KOs across 6109 genomes
growth_PGPR_KO_annot_up<-data.frame(all_PGPR_KO_annot%>%filter(PGP=="Growth"))
rownames(growth_PGPR_KO_annot_up)<-growth_PGPR_KO_annot_up$KO
growth_PGPR_KO_annot_up<-growth_PGPR_KO_annot_up[rownames(CRBD_HQ_NRgenome_growth_KO),]
growth_PGPR_KO_group<-growth_PGPR_KO_annot_up[rownames(CRBD_HQ_NRgenome_growth_KO),]%>%select(Group) # 34
growth_PGPR_KO_group_list<-as.character(growth_PGPR_KO_group$Group)

growth_PGPR_col = list(Color=unique(growth_PGPR_KO_annot_up$Color))
names(growth_PGPR_col$Color)<-as.character(unique(growth_PGPR_KO_annot_up$Group))

GenomeID_annot<-HQ_NRgenome_taxonomy_up_sort[rownames(HQ_NRgenome_taxonomy_up_sort)[rownames(HQ_NRgenome_taxonomy_up_sort)%in%colnames(CRBD_HQ_NRgenome_growth_KO)],]%>%select(PhyClass_collapse)
PhyClass = as.character(GenomeID_annot$PhyClass_collapse) #6109

library(ComplexHeatmap)
set.seed(1013)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Growth/Growth_PGPR_KO_binary_phylogenetic_distribution.pdf",width = 15,height = 8)
m = as.matrix(CRBD_HQ_NRgenome_growth_KO)
dend1 = cluster_within_group(m,PhyClass)
Heatmap(m, cluster_columns = dend1, column_split = 12,col = colorspace::sequential_hcl(n=50,palette = "Dark Mint",rev = TRUE),
    row_title = "cluster_within_PhyClass",row_title_side = 'left',row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize=2),top_annotation = HeatmapAnnotation(PhyClass= PhyClass, col = list(PhyClass=PhyClass_col$Color)),left_annotation= rowAnnotation(growth_PGPR=growth_PGPR_KO_group_list,col=list(growth_PGPR=growth_PGPR_col$Color)),cluster_rows = FALSE,use_raster=FALSE)
dev.off()

## -------- Adaptation ------

CRBD_HQ_NRgenome_adaptation_KO<-CRBD_HQ_NRgenoem_all_PGPR_KO[rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)[rownames(CRBD_HQ_NRgenoem_all_PGPR_KO)%in%(all_PGPR_KO_annot%>%filter(PGP=="Adaptation"))$KO],] # 44 KOs across 6109 genomes
adaptation_PGPR_KO_annot_up<-data.frame(all_PGPR_KO_annot%>%filter(PGP=="Adaptation"))
rownames(adaptation_PGPR_KO_annot_up)<-adaptation_PGPR_KO_annot_up$KO
adaptation_PGPR_KO_annot_up<-adaptation_PGPR_KO_annot_up[rownames(CRBD_HQ_NRgenome_adaptation_KO),]
adaptation_PGPR_KO_group<-adaptation_PGPR_KO_annot_up[rownames(CRBD_HQ_NRgenome_adaptation_KO),]%>%select(Group) # 38
adaptation_PGPR_KO_group_list<-as.character(adaptation_PGPR_KO_group$Group)

adaptation_PGPR_col = list(Color=unique(adaptation_PGPR_KO_annot_up$Color))
names(adaptation_PGPR_col$Color)<-as.character(unique(adaptation_PGPR_KO_annot_up$Group))

GenomeID_annot<-HQ_NRgenome_taxonomy_up_sort[rownames(HQ_NRgenome_taxonomy_up_sort)[rownames(HQ_NRgenome_taxonomy_up_sort)%in%colnames(CRBD_HQ_NRgenome_adaptation_KO)],]%>%select(PhyClass_collapse)
PhyClass = as.character(GenomeID_annot$PhyClass_collapse) #6109

library(ComplexHeatmap)
set.seed(1013)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Adaptation/Adaptation_PGPR_KO_binary_phylogenetic_distribution.pdf",width = 15,height = 8)
m = as.matrix(CRBD_HQ_NRgenome_adaptation_KO)
dend1 = cluster_within_group(m,PhyClass)
Heatmap(m, cluster_columns = dend1, column_split = 12,col = colorspace::sequential_hcl(n=50,palette = "Dark Mint",rev = TRUE),
    row_title = "cluster_within_PhyClass",row_title_side = 'left',row_names_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize=2),top_annotation = HeatmapAnnotation(PhyClass= PhyClass, col = list(PhyClass=PhyClass_col$Color)),left_annotation= rowAnnotation(adaptation_PGPR=adaptation_PGPR_KO_group_list,col=list(adaptation_PGPR=adaptation_PGPR_col$Color)),cluster_rows = FALSE,use_raster=FALSE)
dev.off()


# PGPR nutrient category - upset and Genus ratio
##  nutrient category - summarize

## -------- define module --------

Nif_KO<-c("K02586","K02591","K02588")
Anf_KO<-c("K02586","K02591","K02588","K00531")
Vnf_KO<-c("K22896","K22897","K22898","K22899")
P_dissolving_KO<-c("K00117")
C_P_lyase_KO<-c("K06164","K06166","K06165","K05780","K06162","K06163","K05781","K06167","K05774")
AEP_typeI_KO<-c("K05306","K03430")
AEP_typeII_KO<-c("K00206","K19670","K03430")
Phosphonoacetate_KO<-c("K19670")
Phytase_KO<-c("K01083","K01093")
phosphomonoesterase_KO<-c("K01077","K01113","K07093","K03788","K09474","K01078")
phospho_other_esterase_KO<-c("K07048","K01126")
# K01077,K01113,K07093,K03788,K09474,K01078,K01083,K01093,K01126,K07048
ICS_KO<-c("K01851","K02361","K02552","K14759")
IPL_KO<-c("K04782")
SAS_KO<-c("K04781")
Yesiniabactin_KO<-c("K04783","K04784","K04785","K04786") # 补充了K04786,结果有一个Pub_2908669403变化，有K04785但是没有K04786
Pyochelin_KO<-c("K12238","K12239","K12240","K12241")
Mycobactin_KO<-c("K04787","K04788","K04789","K04790","K04791","K04792","K04793")
siderophore_shared_type2_KO<-c("K01252","K00216","K02363")
Bacillibactin_KO<-c("K04780")
Vibriobactin_KO<-c("K04778","K12237")
Enterochelin_KO<-c("K02362","K02364","K24147")
Myxochelin_KO<-c("K15653")

## ----- N fixation---
CRBD_HQ_NRgenome_nutrient_Nif<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(Nif_KO)%>%mutate(Nif= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Nif$Nif[CRBD_HQ_NRgenome_nutrient_Nif$Nif<length(Nif_KO)]=0
CRBD_HQ_NRgenome_nutrient_Nif$Nif[CRBD_HQ_NRgenome_nutrient_Nif$Nif==length(Nif_KO)]=1
CRBD_HQ_NRgenome_nutrient_Nif<-CRBD_HQ_NRgenome_nutrient_Nif%>%select(Nif) # 1295 has
CRBD_HQ_NRgenome_nutrient_Anf<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(Anf_KO)%>%mutate(Anf= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Anf$Anf[CRBD_HQ_NRgenome_nutrient_Anf$Anf<length(Anf_KO)]=0
CRBD_HQ_NRgenome_nutrient_Anf$Anf[CRBD_HQ_NRgenome_nutrient_Anf$Anf==length(Anf_KO)]=1
CRBD_HQ_NRgenome_nutrient_Anf<-CRBD_HQ_NRgenome_nutrient_Anf%>%select(Anf) # 5 has Anf
CRBD_HQ_NRgenome_nutrient_Vnf<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Vnf_KO))%>%mutate(Vnf= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Vnf$Vnf[CRBD_HQ_NRgenome_nutrient_Vnf$Vnf<length(Vnf_KO)]=0
CRBD_HQ_NRgenome_nutrient_Vnf$Vnf[CRBD_HQ_NRgenome_nutrient_Vnf$Vnf==length(Vnf_KO)]=1
CRBD_HQ_NRgenome_nutrient_Vnf<-CRBD_HQ_NRgenome_nutrient_Vnf%>%select(Vnf) #none has Vnf
CRBD_HQ_NRgenome_nutrient_N_fixation<-data.frame(cbind(CRBD_HQ_NRgenome_nutrient_Nif,CRBD_HQ_NRgenome_nutrient_Anf,CRBD_HQ_NRgenome_nutrient_Vnf)%>%mutate(N_fixation=if_else(Nif+Anf+Vnf>0,1,0))%>%select(N_fixation)) # 1295 has


## ----- C_P_lyase ---
CRBD_HQ_NRgenome_nutrient_C_P_lyase<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(C_P_lyase_KO))%>%mutate(C_P_lyase= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_C_P_lyase$C_P_lyase[CRBD_HQ_NRgenome_nutrient_C_P_lyase$C_P_lyase<length(C_P_lyase_KO)]=0
CRBD_HQ_NRgenome_nutrient_C_P_lyase$C_P_lyase[CRBD_HQ_NRgenome_nutrient_C_P_lyase$C_P_lyase==length(C_P_lyase_KO)]=1
CRBD_HQ_NRgenome_nutrient_C_P_lyase<-CRBD_HQ_NRgenome_nutrient_C_P_lyase%>%select(C_P_lyase) # 1747 has

## ----- AEP_typeI ---
CRBD_HQ_NRgenome_nutrient_AEP_typeI<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(AEP_typeI_KO))%>%mutate(AEP_typeI= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_AEP_typeI$AEP_typeI[CRBD_HQ_NRgenome_nutrient_AEP_typeI$AEP_typeI<length(AEP_typeI_KO)]=0
CRBD_HQ_NRgenome_nutrient_AEP_typeI$AEP_typeI[CRBD_HQ_NRgenome_nutrient_AEP_typeI$AEP_typeI==length(AEP_typeI_KO)]=1
CRBD_HQ_NRgenome_nutrient_AEP_typeI<-CRBD_HQ_NRgenome_nutrient_AEP_typeI%>%select(AEP_typeI) # 455 has

## ----- AEP_typeII ---
CRBD_HQ_NRgenome_nutrient_AEP_typeII<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(AEP_typeII_KO))%>%mutate(AEP_typeII= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_AEP_typeII$AEP_typeII[CRBD_HQ_NRgenome_nutrient_AEP_typeII$AEP_typeII<length(AEP_typeII_KO)]=0
CRBD_HQ_NRgenome_nutrient_AEP_typeII$AEP_typeII[CRBD_HQ_NRgenome_nutrient_AEP_typeII$AEP_typeII==length(AEP_typeII_KO)]=1
CRBD_HQ_NRgenome_nutrient_AEP_typeII<-CRBD_HQ_NRgenome_nutrient_AEP_typeII%>%select(AEP_typeII) # 567 has

## ----- Phosphonoacetate ---
CRBD_HQ_NRgenome_nutrient_Phosphonoacetate<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Phosphonoacetate_KO))%>%mutate(Phosphonoacetate= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Phosphonoacetate$Phosphonoacetate[CRBD_HQ_NRgenome_nutrient_Phosphonoacetate$Phosphonoacetate<length(Phosphonoacetate_KO)]=0
CRBD_HQ_NRgenome_nutrient_Phosphonoacetate$Phosphonoacetate[CRBD_HQ_NRgenome_nutrient_Phosphonoacetate$Phosphonoacetate==length(Phosphonoacetate_KO)]=1
CRBD_HQ_NRgenome_nutrient_Phosphonoacetate<-CRBD_HQ_NRgenome_nutrient_Phosphonoacetate%>%select(Phosphonoacetate) # 1046 has

## ---- NSAP: non-specific acid Phosphatase ---
CRBD_HQ_NRgenome_nutrient_NSAP<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(c("K03788","K09474","K01078")))%>%mutate(NSAP= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_NSAP$NSAP[CRBD_HQ_NRgenome_nutrient_NSAP$NSAP>0]=1
CRBD_HQ_NRgenome_nutrient_NSAP<-CRBD_HQ_NRgenome_nutrient_NSAP%>%select(NSAP) # 929 has, after extracellular none NSAP

## ---- ALP: Alkaline Phosphatase ---

CRBD_HQ_NRgenome_nutrient_ALP<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(c("K01077","K01113","K07093")))%>%mutate(ALP= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_ALP$ALP[CRBD_HQ_NRgenome_nutrient_ALP$ALP>0]=1
CRBD_HQ_NRgenome_nutrient_ALP<-CRBD_HQ_NRgenome_nutrient_ALP%>%select(ALP) # 2176 has

## ---- Phosphodiesterase ---

CRBD_HQ_NRgenome_nutrient_Phosphodiesterase<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(c("K01126")))%>%mutate(Phosphodiesterase= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Phosphodiesterase$Phosphodiesterase[CRBD_HQ_NRgenome_nutrient_Phosphodiesterase$Phosphodiesterase>0]=1
CRBD_HQ_NRgenome_nutrient_Phosphodiesterase<-CRBD_HQ_NRgenome_nutrient_Phosphodiesterase%>%select(Phosphodiesterase) # 5740 has, update 16 

## ---- Phosphotriesterase ---

CRBD_HQ_NRgenome_nutrient_Phosphotriesterase<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(c("K07048")))%>%mutate(Phosphotriesterase= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Phosphotriesterase$Phosphotriesterase[CRBD_HQ_NRgenome_nutrient_Phosphotriesterase$Phosphotriesterase>0]=1
CRBD_HQ_NRgenome_nutrient_Phosphotriesterase<-CRBD_HQ_NRgenome_nutrient_Phosphotriesterase%>%select(Phosphotriesterase) # 1291 has, update 0 

## ---- Phytase ---

CRBD_HQ_NRgenome_nutrient_Phytase<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Phytase_KO))%>%mutate(Phytase= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Phytase$Phytase[CRBD_HQ_NRgenome_nutrient_Phytase$Phytase>0]=1
CRBD_HQ_NRgenome_nutrient_Phytase<-CRBD_HQ_NRgenome_nutrient_Phytase%>%select(Phytase) # 1695 has, update 984

## ---- P_dissolving ---
CRBD_HQ_NRgenome_nutrient_P_dissolving<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(P_dissolving_KO))%>%mutate(P_dissolving= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_P_dissolving$P_dissolving[CRBD_HQ_NRgenome_nutrient_P_dissolving$P_dissolving>0]=1
CRBD_HQ_NRgenome_nutrient_P_dissolving<-CRBD_HQ_NRgenome_nutrient_P_dissolving%>%select(P_dissolving) # 2899 has

## ---- this is version without alkaline phosphatase, diesterase and triesterase ---
CRBD_HQ_NRgenome_nutrient_phosphorus<-cbind(CRBD_HQ_NRgenome_nutrient_C_P_lyase,CRBD_HQ_NRgenome_nutrient_AEP_typeI,CRBD_HQ_NRgenome_nutrient_AEP_typeII,CRBD_HQ_NRgenome_nutrient_Phosphonoacetate,CRBD_HQ_NRgenome_nutrient_NSAP,CRBD_HQ_NRgenome_nutrient_Phytase,CRBD_HQ_NRgenome_nutrient_P_dissolving)
CRBD_HQ_NRgenome_nutrient_phosphorus<-CRBD_HQ_NRgenome_nutrient_phosphorus%>%mutate(Phosphorus=rowSums(across(where(is.numeric))))%>%select(Phosphorus)
CRBD_HQ_NRgenome_nutrient_phosphorus$Phosphorus[CRBD_HQ_NRgenome_nutrient_phosphorus$Phosphorus>0]=1 # 4695 has, update 4253 has

## -------- Siderophore ----------

## ---- ICS---
CRBD_HQ_NRgenome_nutrient_ICS<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(ICS_KO))%>%mutate(ICS= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_ICS$ICS[CRBD_HQ_NRgenome_nutrient_ICS$ICS>0]=1
CRBD_HQ_NRgenome_nutrient_ICS<-CRBD_HQ_NRgenome_nutrient_ICS%>%select(ICS) # has

## ---- IPL---
CRBD_HQ_NRgenome_nutrient_IPL<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(IPL_KO))%>%mutate(IPL= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_IPL$IPL[CRBD_HQ_NRgenome_nutrient_IPL$IPL>0]=1
CRBD_HQ_NRgenome_nutrient_IPL<-CRBD_HQ_NRgenome_nutrient_IPL%>%select(IPL) # has

## ---- SAS---

CRBD_HQ_NRgenome_nutrient_SAS<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(SAS_KO))%>%mutate(SAS= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_SAS$SAS[CRBD_HQ_NRgenome_nutrient_SAS$SAS>0]=1
CRBD_HQ_NRgenome_nutrient_SAS<-CRBD_HQ_NRgenome_nutrient_SAS%>%select(SAS) #  has

## ---- SA---
CRBD_HQ_NRgenome_nutrient_SA<-cbind(CRBD_HQ_NRgenome_nutrient_ICS,CRBD_HQ_NRgenome_nutrient_IPL,CRBD_HQ_NRgenome_nutrient_SAS)%>%mutate(SA=if_else(ICS+IPL==2|SAS==1,1,0))%>%select(SA) # has

## ---- Yesiniabactin---

CRBD_HQ_NRgenome_nutrient_Yesiniabactin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Yesiniabactin_KO))%>%mutate(Yesiniabactin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Yesiniabactin$Yesiniabactin[CRBD_HQ_NRgenome_nutrient_Yesiniabactin$Yesiniabactin<length(Yesiniabactin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Yesiniabactin$Yesiniabactin[CRBD_HQ_NRgenome_nutrient_Yesiniabactin$Yesiniabactin==length(Yesiniabactin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Yesiniabactin<-CRBD_HQ_NRgenome_nutrient_Yesiniabactin%>%select(Yesiniabactin) #  has

## ---- Pyochelin---

CRBD_HQ_NRgenome_nutrient_Pyochelin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Pyochelin_KO))%>%mutate(Pyochelin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Pyochelin$Pyochelin[CRBD_HQ_NRgenome_nutrient_Pyochelin$Pyochelin<length(Pyochelin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Pyochelin$Pyochelin[CRBD_HQ_NRgenome_nutrient_Pyochelin$Pyochelin==length(Pyochelin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Pyochelin<-CRBD_HQ_NRgenome_nutrient_Pyochelin%>%select(Pyochelin) #  has

## ---- Mycobactin---

CRBD_HQ_NRgenome_nutrient_Mycobactin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Mycobactin_KO))%>%mutate(Mycobactin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Mycobactin$Mycobactin[CRBD_HQ_NRgenome_nutrient_Mycobactin$Mycobactin<length(Mycobactin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Mycobactin$Mycobactin[CRBD_HQ_NRgenome_nutrient_Mycobactin$Mycobactin==length(Mycobactin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Mycobactin<-CRBD_HQ_NRgenome_nutrient_Mycobactin%>%select(Mycobactin) #  has

## ---- siderophore_shared_type2 ---

CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(siderophore_shared_type2_KO))%>%mutate(siderophore_shared_type2= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2$siderophore_shared_type2[CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2$siderophore_shared_type2<length(siderophore_shared_type2_KO)]=0
CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2$siderophore_shared_type2[CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2$siderophore_shared_type2==length(siderophore_shared_type2_KO)]=1
CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2<-CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2%>%select(siderophore_shared_type2) #  has

## ---- Bacillibactin---

CRBD_HQ_NRgenome_nutrient_Bacillibactin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Bacillibactin_KO))%>%mutate(Bacillibactin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Bacillibactin$Bacillibactin[CRBD_HQ_NRgenome_nutrient_Bacillibactin$Bacillibactin<length(Bacillibactin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Bacillibactin$Bacillibactin[CRBD_HQ_NRgenome_nutrient_Bacillibactin$Bacillibactin==length(Bacillibactin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Bacillibactin<-CRBD_HQ_NRgenome_nutrient_Bacillibactin%>%select(Bacillibactin) #  has

## ---- Vibriobactin---

CRBD_HQ_NRgenome_nutrient_Vibriobactin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Vibriobactin_KO))%>%mutate(Vibriobactin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Vibriobactin$Vibriobactin[CRBD_HQ_NRgenome_nutrient_Vibriobactin$Vibriobactin<length(Vibriobactin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Vibriobactin$Vibriobactin[CRBD_HQ_NRgenome_nutrient_Vibriobactin$Vibriobactin==length(Vibriobactin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Vibriobactin<-CRBD_HQ_NRgenome_nutrient_Vibriobactin%>%select(Vibriobactin) # none has

## ---- Enterochelin---

CRBD_HQ_NRgenome_nutrient_Enterochelin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Enterochelin_KO))%>%mutate(Enterochelin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Enterochelin$Enterochelin[CRBD_HQ_NRgenome_nutrient_Enterochelin$Enterochelin<length(Enterochelin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Enterochelin$Enterochelin[CRBD_HQ_NRgenome_nutrient_Enterochelin$Enterochelin==length(Enterochelin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Enterochelin<-CRBD_HQ_NRgenome_nutrient_Enterochelin%>%select(Enterochelin) #  has

## ---- Myxochelin---

CRBD_HQ_NRgenome_nutrient_Myxochelin<-data.frame(t(CRBD_HQ_NRgenome_nutrient_KO),check.names = FALSE)%>%select(any_of(Myxochelin_KO))%>%mutate(Myxochelin= rowSums(across(where(is.numeric))))
CRBD_HQ_NRgenome_nutrient_Myxochelin$Myxochelin[CRBD_HQ_NRgenome_nutrient_Myxochelin$Myxochelin<length(Myxochelin_KO)]=0
CRBD_HQ_NRgenome_nutrient_Myxochelin$Myxochelin[CRBD_HQ_NRgenome_nutrient_Myxochelin$Myxochelin==length(Myxochelin_KO)]=1
CRBD_HQ_NRgenome_nutrient_Myxochelin<-CRBD_HQ_NRgenome_nutrient_Myxochelin%>%select(Myxochelin) #  has

## ---- Iron ---

CRBD_HQ_NRgenome_nutrient_Iron<-cbind(CRBD_HQ_NRgenome_nutrient_SA,CRBD_HQ_NRgenome_nutrient_Yesiniabactin,CRBD_HQ_NRgenome_nutrient_Pyochelin,CRBD_HQ_NRgenome_nutrient_Mycobactin,CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2,CRBD_HQ_NRgenome_nutrient_Bacillibactin,CRBD_HQ_NRgenome_nutrient_Vibriobactin,CRBD_HQ_NRgenome_nutrient_Enterochelin,CRBD_HQ_NRgenome_nutrient_Myxochelin)%>%mutate(Iron=if_else((CRBD_HQ_NRgenome_nutrient_SA>0&(CRBD_HQ_NRgenome_nutrient_Yesiniabactin+CRBD_HQ_NRgenome_nutrient_Pyochelin+CRBD_HQ_NRgenome_nutrient_Mycobactin>0))|(CRBD_HQ_NRgenome_nutrient_siderophore_shared_type2&(CRBD_HQ_NRgenome_nutrient_Bacillibactin+CRBD_HQ_NRgenome_nutrient_Vibriobactin+CRBD_HQ_NRgenome_nutrient_Enterochelin+CRBD_HQ_NRgenome_nutrient_Myxochelin)),1,0))%>%select(Iron) # 327 has

## -------- All nutrient together -----------

CRBD_HQ_NRgenome_nutrient_combine<-cbind(CRBD_HQ_NRgenome_nutrient_N_fixation,CRBD_HQ_NRgenome_nutrient_phosphorus,CRBD_HQ_NRgenome_nutrient_Iron)
table(apply(CRBD_HQ_NRgenome_nutrient_combine,1,function(x) sum(x>0))) # Freq: 0    1    2    3 = 1621 3119 1351   18 
CRBD_HQ_NRgenome_nutrient_combine_taxonomy<-inner_join(data.frame(CRBD_HQ_NRgenome_nutrient_combine,NRgenome_GenomeID=row.names(CRBD_HQ_NRgenome_nutrient_combine)),CRBD_HQ_NRgenome_taxonomy)%>%mutate(Freq=N_fixation+Phosphorus+Iron)
data.frame(CRBD_HQ_NRgenome_nutrient_combine_taxonomy%>%filter(Freq==3))%>%group_by(Genus)%>%mutate(SUM=n())%>%select(c(Genus,SUM))%>%unique()%>%as.data.frame() # "Bradyrhizobium" "Allorhizobium"  "Klebsiella"     "Enterobacter"   "Kosakonia"      "Phytobacter"    "Pantoea"        "Paenibacillus" 
#           Genus SUM
# 1 Bradyrhizobium   1
# 2  Allorhizobium   1
# 3     Klebsiella   8
# 4   Enterobacter   2
# 5      Kosakonia   3
# 6    Phytobacter   1
# 7        Pantoea   1
# 8  Paenibacillus   1
write.table(CRBD_HQ_NRgenome_nutrient_combine_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/R_write_CRBD_HQ_NRgenome_nutrient_Freq_taxonomy.csv",quote = FALSE,row.names = FALSE)


##  nutrient category - upset

## ------- upset plot --------
combination_matrix = make_comb_mat(CRBD_HQ_NRgenome_nutrient_combine,mode = 'distinct')
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/Nutrient_PGPR_upset_summary_v3.pdf",width = 6,height = 3)
UpSet(combination_matrix,set_order = c("N_fixation","Phosphorus","Iron"),comb_col = rev(c("#425712", "#6e8a2f", "#99b55b"))[comb_degree(combination_matrix)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

## ------- upset plot --------
combination_matrix_v2<-combination_matrix[comb_degree(combination_matrix)>=1]
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/Nutrient_PGPR_upset_summary_v3.pdf",width = 6,height = 3)
UpSet(combination_matrix_v2,set_order = c("N_fixation","Phosphorus","Iron"),comb_col = rev(c("#425712", "#6e8a2f", "#99b55b"))[comb_degree(combination_matrix)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

## ------- Family composition ------

## extract combination labels and its genomeID

Genome_list<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Genome_list[[i]]=extract_comb(combination_matrix_v2,comb)
}

Comb_label<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Comb_label[[i]]=rep(comb,each=length(extract_comb(combination_matrix_v2, comb)))
}

combination_genome_composition<-data.frame(GenomeID=unlist(Genome_list),Comb_label=unlist(Comb_label))
combination_genome_composition$Comb_label<-factor(combination_genome_composition$Comb_label,levels = unique(combination_genome_composition$Comb_label))
combination_genome_composition_taxa<-left_join(combination_genome_composition,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
combination_genome_composition_taxa$PhyClass_collapse<-factor(combination_genome_composition_taxa$PhyClass_collapse,levels = PhyClass_order)
combination_genome_composition_taxa$Count=1

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/Upset_phyClass_composition_nutrient.pdf",height = 3.5,width = 6)
ggplot(combination_genome_composition_taxa,aes(x=Comb_label,y=Count,fill=PhyClass_collapse))+geom_bar(stat='identity',position = 'stack',width = 0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Genus -----

combination_genome_Genus_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Genus)%>%mutate(GenusSum=n())%>%select(c(Comb_label,Genus,GenusSum))%>%unique())
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_A","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_B","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_E","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_F","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Bacillus_A","Bacillus",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Arthrobacter_I","Arthrobacter",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("^$","Others",combination_genome_Genus_count$Genus)
combination_genome_Genus_count<-combination_genome_Genus_count%>%mutate(Genus_up=if_else(GenusSum>=40,Genus,"Others"))
combination_genome_Genus_count$Genus_up<-factor(combination_genome_Genus_count$Genus_up,levels = c(unique(combination_genome_Genus_count$Genus_up)[unique(combination_genome_Genus_count$Genus_up)!="Others"],"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/Upset_Genus_composition_Nutrient.pdf",height = 5,width = 8)
ggplot(combination_genome_Genus_count,aes(x=Comb_label,y=GenusSum,fill=Genus_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#b2798d','#235931','#CD9BCD','#00CDCD','#ff1ce0','#CBD588','#f075f4','#f7cdcf','#85c1b8','#ccff60','#0f9032','#DA5724','#009897','#757f77','#aa0012','#121457','#03ffe6','#f7ff03','#ff8400','#c2edc0','#fcfbc0','#fce3d7','#d7dbf7','#ffb444','#4cb7cf','#ff9082','black','blue','#bfba9b','#d3dbd3'))+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=2))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Family -----

combination_genome_Family_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Family)%>%mutate(FamilySum=n())%>%select(c(Comb_label,Family,FamilySum))%>%unique())
combination_genome_Family_count$Family<-gsub("^$","Others",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("f__","",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_G","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_H","Bacillaceae",combination_genome_Family_count$Family)

combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(Family%in%family_color_meta$Family,Family,"Others"))
combination_genome_Family_count$Family_up<-factor(combination_genome_Family_count$Family_up,levels = as.character(family_color_meta$Family)[as.character(family_color_meta$Family)%in%unique(combination_genome_Family_count$Family_up)])

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/Upset_Family_composition_Nutrient.pdf",height = 6,width = 6)
ggplot(combination_genome_Family_count,aes(x=Comb_label,y=FamilySum,fill=Family_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values = family_color_meta[levels(combination_genome_Family_count$Family_up),]$Color)+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))+labs(y="NRgenome number")
#scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#d3dbd3'))
dev.off()


# PGPR growth category - upset and Genus ratio
##  growth category - summarize

## -------- define module --------

IAA_TAM_KO<-c(c("K01593","K22433"),c("K00274","K11182"),c("K00128","K14085","K00149","K11817","K22417"))
IAA_IAM_KO<-c("K00466",c("K01426","K21801"))
IAA_IPyA_IPDC_KO<-c(c("K00832","K00838","K05821","K14265","K16903","K03334"),"K04103",c("K00128","K14085","K00149","K11817","K22417"))
IAA_IPyA_YUCCA_KO<-c(c("K00832","K00838","K05821","K14265","K16903","K03334"),"K11816")
IAA_IAN_typeI_KO<-c(c("K11812","K11813"),"K11868",c("K01501")) # "K13035" was removed, not change anything since K13035 does not present in our CRBD genomes
IAA_IAN_typeII_KO<-c(c("K11812","K11813"),"K11868",c("K01721","K20807"),c("K01426","K21801"))
CK_typeI_KO<-c("K00791","K06168","K06169","K22522")
CK_typeII_KO<-c(c("K10760","K22871"),"K10717","K22522")
GA_KO<-c(c("K13789","K13787","K00804"),"K04120","K04121",c("K20657"),"K21118","K21117","K21116") # K20657 has CPS_KS activity for fungi

## ----- IAA_TAM ---
CRBD_HQ_NRgenome_growth_IAA_TAM<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(IAA_TAM_KO))%>%mutate(IAA_TAM=if_else(((K01593)>0)&(K00274>0)&(K00128>0),1,0))%>%select(IAA_TAM) # 1086 has 

## ----- IAA_IAM ---

CRBD_HQ_NRgenome_growth_IAA_IAM<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(IAA_IAM_KO))%>%mutate(IAA_IAM=if_else(((K00466)>0)&(K01426+K21801>0),1,0))%>%select(IAA_IAM) # 180 has

## ----- IAA_IPyA_IPDC ---

CRBD_HQ_NRgenome_growth_IAA_IPyA_IPDC<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(IAA_IPyA_IPDC_KO))%>%mutate(IAA_IPyA_IPDC=if_else((K00832>0)&(K04103>0)&(K00128>0),1,0))%>%select(IAA_IPyA_IPDC) # 380 has

## ----- IAA_IPyA_YUCCA ---

CRBD_HQ_NRgenome_growth_IAA_IPyA_YUCCA<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(IAA_IPyA_YUCCA_KO))%>%mutate(IAA_IPyA_YUCCA=if_else((K00832>0)&(K11816>0),1,0))%>%select(IAA_IPyA_YUCCA) # 291 has

## ----- IAA_IAN_typeI ---

CRBD_HQ_NRgenome_growth_IAA_IAN_typeI<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(IAA_IAN_typeI_KO)) # only K01501 left, none has full path, only can synthesized from IAN

## ----- IAA_IAN_typeII ---

CRBD_HQ_NRgenome_growth_IAA_IAN_typeII<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(IAA_IAN_typeII_KO))# only K01721 K20807 K01426 K21801 left, none has full path, only can synthesized from IAN


## ----- CK_typeI ---
CRBD_HQ_NRgenome_growth_CK_typeI<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(CK_typeI_KO))%>%mutate(CK_typeI =if_else(((K00791+K06169+K22522)==3),1,0))%>%select(CK_typeI) # 222

## ----- CK_typeII ---
CRBD_HQ_NRgenome_growth_CK_typeII<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(CK_typeII_KO))%>%mutate(CK_typeII =if_else(((K22871)>0&(K10717+K22522)==2),1,0))%>%select(CK_typeII) # K10760 is absent, so, here the formula change from (K10760+K22871)>0  to (K22871)>0, none of the path is full

## ----- GA ---
CRBD_HQ_NRgenome_growth_GA<-data.frame(t(CRBD_HQ_NRgenome_growth_KO),check.names = FALSE)%>%select(any_of(GA_KO))%>%mutate(GA=if_else(((K13789+K13787)>=1)&((K04120+K04121+K21117+K21116+K21118)==5),1,0))%>%select(GA) # 105 has


## -------- All growth together -----------

CRBD_HQ_NRgenome_growth_IAA<-cbind(CRBD_HQ_NRgenome_growth_IAA_IAM,CRBD_HQ_NRgenome_growth_IAA_TAM,CRBD_HQ_NRgenome_growth_IAA_IPyA_IPDC)%>%mutate(IAA=rowSums(across(where(is.numeric))))%>%select(IAA)
CRBD_HQ_NRgenome_growth_IAA$IAA[CRBD_HQ_NRgenome_growth_IAA$IAA>0]=1
table(CRBD_HQ_NRgenome_growth_IAA$IAA) # 1609 HQ NRgenome has IAA biosynthesis potential
CRBD_HQ_NRgenome_growth_combine<-cbind(CRBD_HQ_NRgenome_growth_IAA,CRBD_HQ_NRgenome_growth_CK_typeI,CRBD_HQ_NRgenome_growth_GA)%>%dplyr::rename(CK=CK_typeI)
table(apply(CRBD_HQ_NRgenome_growth_combine,1,function(x) sum(x>0))) # Freq: 0    1    2  = 4314 1654  141 
CRBD_HQ_NRgenome_growth_combine_taxonomy<-inner_join(data.frame(CRBD_HQ_NRgenome_growth_combine,NRgenome_GenomeID=row.names(CRBD_HQ_NRgenome_growth_combine)),CRBD_HQ_NRgenome_taxonomy)%>%mutate(Freq=IAA+CK+GA)
data.frame(CRBD_HQ_NRgenome_growth_combine_taxonomy%>%filter(Freq==2))%>%group_by(Genus)%>%mutate(SUM=n())%>%select(c(Genus,SUM))%>%unique()%>%as.data.frame() 
#            Genus SUM
# 1 Bradyrhizobium  25
# 2  Mesorhizobium  35
# 3      Rhizobium  21
# 4  Pseudomonas_E  59
# 5  Acinetobacter   1
write.table(CRBD_HQ_NRgenome_growth_combine_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/R_write_CRBD_HQ_NRgenome_nutrient_Freq_taxonomy.csv",quote = FALSE,row.names = FALSE)

temp<-inner_join(data.frame(CRBD_HQ_NRgenome_growth_combine,NRgenomeID=row.names(CRBD_HQ_NRgenome_growth_combine)),CRBD_HQ_NRgenome_taxonomy,by=c("NRgenomeID"="NRgenome_GenomeID"))%>%filter(GA>0)%>%select(c(Genus))%>%unique()


##  growth category - upset


## ------- upset plot --------
combination_matrix = make_comb_mat(CRBD_HQ_NRgenome_growth_combine,mode = 'distinct')

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Growth/growth_PGPR_upset_summary.pdf",width = 4,height = 2.5)
UpSet(combination_matrix,set_order = c("IAA","CK","GA"),comb_col = rev(c("#425712", "#6e8a2f", "#99b55b"))[comb_degree(combination_matrix)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

combination_matrix_v2<-combination_matrix[comb_degree(combination_matrix)>=1]
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Growth/growth_PGPR_upset_summary_v2.pdf",width = 4,height = 2.5)
UpSet(combination_matrix_v2,set_order = c("IAA","CK","GA"),comb_col = rev(c("#425712", "#6e8a2f", "#99b55b"))[comb_degree(combination_matrix_v2)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

## ------- Family composition ------

## extract combination labels and its genomeID

Genome_list<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Genome_list[[i]]=extract_comb(combination_matrix_v2,comb)
}

Comb_label<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Comb_label[[i]]=rep(comb,each=length(extract_comb(combination_matrix_v2, comb)))
}

combination_genome_composition<-data.frame(GenomeID=unlist(Genome_list),Comb_label=unlist(Comb_label))
combination_genome_composition$Comb_label<-factor(combination_genome_composition$Comb_label,levels = unique(combination_genome_composition$Comb_label))
combination_genome_composition_taxa<-left_join(combination_genome_composition,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
combination_genome_composition_taxa$PhyClass_collapse<-factor(combination_genome_composition_taxa$PhyClass_collapse,levels = PhyClass_order)
combination_genome_composition_taxa$Count=1

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Growth/Upset_phyClass_composition_growth.pdf",height = 3,width = 6)
ggplot(combination_genome_composition_taxa,aes(x=Comb_label,y=Count,fill=PhyClass_collapse))+geom_bar(stat='identity',position = 'stack',width = 0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Genus -----

combination_genome_Genus_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Genus)%>%mutate(GenusSum=n())%>%select(c(Comb_label,Genus,GenusSum))%>%unique())
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_A","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_B","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_E","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_F","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Bacillus_A","Bacillus",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Arthrobacter_I","Arthrobacter",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("^$","Others",combination_genome_Genus_count$Genus)
combination_genome_Genus_count<-combination_genome_Genus_count%>%mutate(Genus_up=if_else(GenusSum>=10,Genus,"Others"))
combination_genome_Genus_count$Genus_up<-factor(combination_genome_Genus_count$Genus_up,levels = c(unique(combination_genome_Genus_count$Genus_up)[unique(combination_genome_Genus_count$Genus_up)!="Others"],"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Growth/Upset_Genus_composition_Growth.pdf",height = 5,width = 8)
ggplot(combination_genome_Genus_count,aes(x=Comb_label,y=GenusSum,fill=Genus_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#b2798d','#235931','#CD9BCD','#00CDCD','#ff1ce0','#CBD588','#f075f4','#f7cdcf','#85c1b8','#ccff60','#0f9032','#DA5724','#009897','#757f77','#aa0012','#121457','#03ffe6','#f7ff03','#ff8400','#c2edc0','#fcfbc0','#fce3d7','#d7dbf7','#ffb444','#4cb7cf','#ff9082','black','blue','#bfba9b','#d3dbd3'))+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=2))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Family -----

combination_genome_Family_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Family)%>%mutate(FamilySum=n())%>%select(c(Comb_label,Family,FamilySum))%>%unique())
combination_genome_Family_count$Family<-gsub("^$","Others",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("f__","",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_G","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_H","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(Family%in%family_color_meta$Family,Family,"Others"))
combination_genome_Family_count$Family_up<-factor(combination_genome_Family_count$Family_up,levels = as.character(family_color_meta$Family)[as.character(family_color_meta$Family)%in%unique(combination_genome_Family_count$Family_up)])

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Growth/Upset_Family_composition_Growth.pdf",height = 6,width = 6)
ggplot(combination_genome_Family_count,aes(x=Comb_label,y=FamilySum,fill=Family_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values = family_color_meta[levels(combination_genome_Family_count$Family_up),]$Color)+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))+labs(y="NRgenome number")
#scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#d3dbd3'))
dev.off()


# PGPR adaptation category - upset and Genus ratio
##  adaptation category - summarize

## -------- define module --------

ACC_deaminase_KO<-c("K01505")
EFE_microbe_KO<-c("K21815")
EFE_plant_KO<-c("K00789",c("K01762","K20772"),"K05933")
#SA has summarized in siderophore sections
JA_KO<-c("K19246","K26082","K10525","K05894","K00232","K10527","K07513") # ß-oxidation genes are not clearly listed in review article
ABA_KO<-c("K09838","K14594","K09840","K09841","K09842","K15631")
ISR_pyocynin_KO<-c("K13063","K20261","K06998","K20260","K20262","K21103","K20940")
ISR_DAPG_KO<-c("K22838","K22839","K22840")
ISR_surfactin_KO<-c("K15654","K15655","K15656")
ISR_fengycin_KO<-c("K15664","K15665","K15666","K15667","K15668")
ISR_syringomycin_KO<-c("K16125")
ISR_VOC_KO<-c("K01652","K01653","K11258","K01575")

## ----- ACC_deaminase ---
CRBD_HQ_NRgenome_adaptation_ACC_deaminase<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(K01505)%>%mutate(ACC_deaminase=if_else(K01505>0,1,0))%>%select(ACC_deaminase) # 2486 has vs 2502

## ----- EFE_microbe ---
CRBD_HQ_NRgenome_adaptation_EFE_microbe<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(EFE_microbe_KO))%>%mutate(EFE_microbe=if_else(K21815>0,1,0))%>%select(EFE_microbe) # 22 has vs 31

## ----- EFE_plant ---
CRBD_HQ_NRgenome_adaptation_EFE_plant<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(EFE_plant_KO)) # lack of ACC oxidase K05933, none has full path

## ----- JA ---
CRBD_HQ_NRgenome_adaptation_JA<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(JA_KO)) # lack of "K26082" and K00232 none has full path

## ----- ABA ---
CRBD_HQ_NRgenome_adaptation_ABA<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ABA_KO)) # lack of K14594, none has full path

## ----- ISR_pyocynin  ---
CRBD_HQ_NRgenome_adaptation_ISR_pyocynin<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ISR_pyocynin_KO))%>%mutate(ISR_pyocynin=if_else((K13063+K20261+K06998+K20260+K20262+K21103+K20940)==length(ISR_pyocynin_KO),1,0))%>%select(ISR_pyocynin) # 3 has full path

## ----- ISR_DAPG  ---
CRBD_HQ_NRgenome_adaptation_ISR_DAPG<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ISR_DAPG_KO))%>%mutate(ISR_DAPG=if_else((K22838+K22839+K22840)==length(ISR_DAPG_KO),1,0))%>%select(ISR_DAPG) # 27 has full path

## ----- ISR_surfactin  ---
CRBD_HQ_NRgenome_adaptation_ISR_surfactin<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ISR_surfactin_KO))%>%mutate(ISR_surfactin=if_else((K15654+K15655+K15656)==3,1,0))%>%select(ISR_surfactin) # 41 has

## ----- ISR_fengycin  ---
CRBD_HQ_NRgenome_adaptation_ISR_fengycin<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ISR_fengycin_KO))%>%mutate(ISR_fengycin=if_else((K15664+K15665+K15666+K15667+K15668)==length(ISR_fengycin_KO),1,0))%>%select(ISR_fengycin) # 27 has

## ----- ISR_syringomycin  ---
CRBD_HQ_NRgenome_adaptation_ISR_syringomycin<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ISR_syringomycin_KO)) # K15664 is absent

## ----- ISR_VOC  ---
CRBD_HQ_NRgenome_adaptation_ISR_VOC<-data.frame(t(CRBD_HQ_NRgenome_adaptation_KO),check.names = FALSE)%>%select(any_of(ISR_VOC_KO))%>%mutate(ISR_VOC=if_else((K01652+K01653+K11258+K01575)==length(ISR_VOC_KO),1,0))%>%select(ISR_VOC) # 84 has

## ----- ISR  ---

CRBD_HQ_NRgenome_adaptation_ISR<-cbind(CRBD_HQ_NRgenome_adaptation_ISR_pyocynin,CRBD_HQ_NRgenome_adaptation_ISR_DAPG,CRBD_HQ_NRgenome_adaptation_ISR_fengycin,CRBD_HQ_NRgenome_adaptation_ISR_VOC)%>%mutate(ISR=rowSums(across(where(is.numeric))))%>%select(ISR)
CRBD_HQ_NRgenome_adaptation_ISR$ISR[CRBD_HQ_NRgenome_adaptation_ISR$ISR>0]=1 # 139 has full ISR

## ----- adaptation combine ---

CRBD_HQ_NRgenome_adaptation_combine<-cbind(CRBD_HQ_NRgenome_adaptation_ACC_deaminase,CRBD_HQ_NRgenome_adaptation_EFE_microbe,CRBD_HQ_NRgenome_nutrient_SA,CRBD_HQ_NRgenome_adaptation_ISR)%>%dplyr::rename(EFE=EFE_microbe)%>%select(!ISR) #JA，ABA and Ethylene plant does not exist
table(apply(CRBD_HQ_NRgenome_adaptation_combine,1,function(x) sum(x>0))) # Freq:0    1    2    3  = 3339 2635  131    4
CRBD_HQ_NRgenome_adaptation_combine_taxonomy<-inner_join(data.frame(CRBD_HQ_NRgenome_adaptation_combine,NRgenome_GenomeID=row.names(CRBD_HQ_NRgenome_adaptation_combine)),CRBD_HQ_NRgenome_taxonomy)%>%mutate(Freq=ACC_deaminase+SA+EFE)
data.frame(CRBD_HQ_NRgenome_adaptation_combine_taxonomy%>%filter(Freq==3))%>%group_by(Genus)%>%mutate(SUM=n())%>%select(c(Genus,SUM))%>%unique()%>%as.data.frame() 
#          Genus SUM
# 1      Lentzea   1
# 2 Streptomyces   2
# 3                1
write.table(CRBD_HQ_NRgenome_adaptation_combine_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Nutrient/R_write_CRBD_HQ_NRgenome_nutrient_Freq_taxonomy.csv",quote = FALSE,row.names = FALSE)


##  adaptation category - upset

## ------- upset plot --------
combination_matrix = make_comb_mat(CRBD_HQ_NRgenome_adaptation_combine%>%select(c(ACC_deaminase,EFE,SA)),mode = 'distinct')

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Adaptation/adaptation_PGPR_upset_summary.pdf",width = 6,height = 3)
UpSet(combination_matrix,set_order = c("ACC_deaminase","EFE","SA"),comb_col = rev(c("#425712", "#6e8a2f", "#99b55b"))[comb_degree(combination_matrix)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

## ------- upset plot --------
combination_matrix_v2<-combination_matrix[comb_degree(combination_matrix)>=1]
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Adaptation/adaptation_PGPR_upset_summary_v2.pdf",width = 6,height = 3)
UpSet(combination_matrix_v2,set_order = c("ACC_deaminase","EFE","SA"),comb_col = rev(c("#425712", "#6e8a2f", "#99b55b"))[comb_degree(combination_matrix)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

## ------- Family composition ------

## extract combination labels and its genomeID

Genome_list<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Genome_list[[i]]=extract_comb(combination_matrix_v2,comb)
}

Comb_label<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Comb_label[[i]]=rep(comb,each=length(extract_comb(combination_matrix_v2, comb)))
}

combination_genome_composition<-data.frame(GenomeID=unlist(Genome_list),Comb_label=unlist(Comb_label))
combination_genome_composition$Comb_label<-factor(combination_genome_composition$Comb_label,levels = unique(combination_genome_composition$Comb_label))
combination_genome_composition_taxa<-left_join(combination_genome_composition,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
combination_genome_composition_taxa$PhyClass_collapse<-factor(combination_genome_composition_taxa$PhyClass_collapse,levels = PhyClass_order)
combination_genome_composition_taxa$Count=1

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Adaptation/Upset_phyClass_composition_adaptation.pdf",height = 3,width = 6)
ggplot(combination_genome_composition_taxa,aes(x=Comb_label,y=Count,fill=PhyClass_collapse))+geom_bar(stat='identity',position = 'stack',width = 0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Genus -----

combination_genome_Genus_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Genus)%>%mutate(GenusSum=n())%>%select(c(Comb_label,Genus,GenusSum))%>%unique())
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_A","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_B","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_E","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_F","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Bacillus_A","Bacillus",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Arthrobacter_I","Arthrobacter",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("^$","Others",combination_genome_Genus_count$Genus)
combination_genome_Genus_count<-combination_genome_Genus_count%>%mutate(Genus_up=if_else(GenusSum>=10,Genus,"Others"))
combination_genome_Genus_count$Genus_up<-factor(combination_genome_Genus_count$Genus_up,levels = c(unique(combination_genome_Genus_count$Genus_up)[unique(combination_genome_Genus_count$Genus_up)!="Others"],"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Adaptation/Upset_Genus_composition_Adaptation.pdf",height = 6,width = 10)
ggplot(combination_genome_Genus_count,aes(x=Comb_label,y=GenusSum,fill=Genus_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#b2798d','#235931','#CD9BCD','#00CDCD','#ff1ce0','#CBD588','#f075f4','#f7cdcf','#85c1b8','#ccff60','#0f9032','#DA5724','#009897','#757f77','#aa0012','#121457','#03ffe6','#f7ff03','#ff8400','#c2edc0','#fcfbc0','#fce3d7','#d7dbf7','#ffb444','#4cb7cf','#ff9082','black','blue','#bfba9b','#d3dbd3'))+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=2))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Family -----

combination_genome_Family_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Family)%>%mutate(FamilySum=n())%>%select(c(Comb_label,Family,FamilySum))%>%unique())
combination_genome_Family_count$Family<-gsub("^$","Others",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("f__","",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_G","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_H","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(Family%in%family_color_meta$Family,Family,"Others"))
combination_genome_Family_count$Family_up<-factor(combination_genome_Family_count$Family_up,levels = as.character(family_color_meta$Family)[as.character(family_color_meta$Family)%in%unique(combination_genome_Family_count$Family_up)])

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Adaptation/Upset_Family_composition_Adaptation.pdf",height = 6,width = 6)
ggplot(combination_genome_Family_count,aes(x=Comb_label,y=FamilySum,fill=Family_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values = family_color_meta[levels(combination_genome_Family_count$Family_up),]$Color)+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))+labs(y="NRgenome number")
#scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#d3dbd3'))
dev.off()


# All PGPR category combine -- upset and genome level heatmap
##  all PGPR category - upset

## Summarize the number of genomes harbor individual functions
CRBD_HQ_NRgenome_PGPR_combine<-cbind(CRBD_HQ_NRgenome_nutrient_combine,CRBD_HQ_NRgenome_growth_combine,CRBD_HQ_NRgenome_adaptation_combine)%>%dplyr::rename(Nitrogen=N_fixation)%>%dplyr::rename(Ethylene=EFE)
row.names(CRBD_HQ_NRgenome_PGPR_combine)<-row.names(CRBD_HQ_NRgenome_nutrient_combine)
write.csv(CRBD_HQ_NRgenome_PGPR_combine,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_All_PGPR_9_category_combine_table.csv",sep=',',quote = FALSE,row.names = TRUE)
PGPR_Freq<-data.frame(Freq=apply(CRBD_HQ_NRgenome_PGPR_combine,1,function(x) sum(x>0)))
Freq_summary<-data.frame(PGPR_Freq%>%group_by(Freq)%>%mutate(SUM=n())%>%unique())
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/PGPR_Freq_histogram_at_individual_genome_level.pdf",width = 3,height = 2.5)
p<-ggplot(PGPR_Freq,aes(x=Freq))+geom_histogram(bins = 20,color="grey",fill="#1e5955",size=0.1)+theme_bw()
data_counts <- ggplot_build(p)$data[[1]]%>%filter(y>0)
p + geom_text(data = data_counts, aes(x = x, y = y, label = y),
              vjust = -0.5, color = "#249c94", size = 2)+labs(x="PGP richness",y="Number of genomes")+ylim(0,2500)+theme(axis.text = element_text(size=6))
dev.off()

## summarize CRBC versus published PGP contribution

CRBD_NRgenome_sourceSum_per_NRgenome<-data.frame(CRBD_genome_meta_up2%>%group_by(NRgenome_GenomeID,Source)%>%mutate(Source_per_NRgenome=n())%>%select(c(NRgenome_GenomeID,Source,Source_per_NRgenome))%>%unique())%>%filter(NRgenome_GenomeID%in%CRBD_HQ_NRgenome$GenomeID) # 6111, since two of the NRgenome: Pub_2849560528" "Iso_Md_B527b2" has both CRBC and Published genomes, when group_by they will be summed twice
CRBD_NRgenome_Sum_per_NRgenome<-data.frame(CRBD_genome_meta_up2%>%group_by(NRgenome_GenomeID)%>%mutate(Sum_per_NRgenome=n())%>%select(c(NRgenome_GenomeID,Sum_per_NRgenome))%>%unique())%>%filter(NRgenome_GenomeID%in%CRBD_HQ_NRgenome$GenomeID) # 6109
CRBD_NRgenome_CRBC_ratio<-left_join(CRBD_NRgenome_Sum_per_NRgenome,CRBD_NRgenome_sourceSum_per_NRgenome)%>%filter(Source=="CRBC")%>%mutate(CRBC_ratio=Source_per_NRgenome/Sum_per_NRgenome) # 4052 NRgenomes has CRBC genomes
CRBD_NRgenome_Published_ratio<-left_join(CRBD_NRgenome_Sum_per_NRgenome,CRBD_NRgenome_sourceSum_per_NRgenome)%>%filter(Source=="Published")%>%mutate(Published_ratio=Source_per_NRgenome/Sum_per_NRgenome) # 2061 NRgenomes has Published genomes
CRBD_NRgenome_Published_ratio$NRgenome_GenomeID[CRBD_NRgenome_Published_ratio$NRgenome_GenomeID%in%CRBD_NRgenome_CRBC_ratio$NRgenome_GenomeID] 
CRBD_NRgenome_CRBC_vs_Published_ratio<-full_join(CRBD_NRgenome_CRBC_ratio%>%select(c(NRgenome_GenomeID,CRBC_ratio)),CRBD_NRgenome_Published_ratio%>%select(NRgenome_GenomeID,Published_ratio)) ## only two genomes, including NRgenome of "Pub_2849560528" "Iso_Md_B527b2" has both CRBC and Published
CRBD_NRgenome_CRBC_vs_Published_ratio%>%filter(CRBC_ratio==0.5&Published_ratio==0.5) 
#     NRgenome_GenomeID CRBC_ratio Published_ratio
# 1    Pub_2849560528        0.5             0.5
# 2     Iso_Md_B527b2        0.5             0.5

CRBD_NRgenome_CRBC_vs_Published_ratio[is.na(CRBD_NRgenome_CRBC_vs_Published_ratio)]=0 # 6109
CRBD_NRgenome_CRBC_vs_Published_ratio$NRgenome_GenomeID[!CRBD_NRgenome_CRBC_vs_Published_ratio$NRgenome_GenomeID%in%CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID] # "Iso_Wt_W2680b1"         "MAG_Rc_HNminiLNg4Mb035"
CRBD_HQ_NRgenome_CRBC_vs_Published_ratio<-CRBD_NRgenome_CRBC_vs_Published_ratio%>%filter(NRgenome_GenomeID%in%CRBD_HQ_NRgenome_taxonomy$NRgenome_GenomeID) #6109
CRBD_HQ_NRgenome_source_category<-CRBD_HQ_NRgenome_CRBC_vs_Published_ratio%>%mutate(Source_cat=case_when(CRBC_ratio>0&Published_ratio>0~"Both",CRBC_ratio>0&Published_ratio==0~"CRBC",CRBC_ratio==0&Published_ratio>0~"Published"))
table(CRBD_HQ_NRgenome_source_category$Source_cat)
     # Both      CRBC Published 
     #    2      4048      2059

CRBD_HQ_NRgenome_PGPR_combine_add_source_category<-inner_join(data.frame(NRgenome_GenomeID=row.names(CRBD_HQ_NRgenome_PGPR_combine),CRBD_HQ_NRgenome_PGPR_combine),CRBD_HQ_NRgenome_source_category%>%select(c(NRgenome_GenomeID,Source_cat)))
rownames(CRBD_HQ_NRgenome_PGPR_combine_add_source_category)<-CRBD_HQ_NRgenome_PGPR_combine_add_source_category$NRgenome_GenomeID
CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource<-CRBD_HQ_NRgenome_PGPR_combine_add_source_category%>%select(!NRgenome_GenomeID)%>%group_by(Source_cat)%>%mutate_at(vars(-Source_cat),function(x) sum(x>0))%>%unique()%>%as.data.frame()
rownames(CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource)<-CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource$Source_cat
CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2<-CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource%>%select(!Source_cat)%>%t()%>%as.data.frame()%>%mutate(Classes=(c(rep("Nutrient",3),rep("Growth",3),rep("Stress",3))))
CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2$PGP=row.names(CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2)
CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt<-melt(CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2,id.vars = c("Classes","PGP"),variable.name = "Source",value.name = "Presence_SUM")
CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt$Classes<-factor(CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt$Classes,levels = c("Nutrient","Growth","Stress"))
CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt$PGP<-factor(CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt$PGP,levels =c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene") )

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/CRBC_vs_Published_PGP_presence_genome_number.pdf",height = 4,width = 5)
ggplot(data=CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt,aes(x=PGP,y=Presence_SUM,fill=Source,color=Source,group=Source))+geom_bar(stat = 'identity',position = "dodge")+ggforce::facet_row(~Classes,space = 'free',scales = 'free_x')+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))+scale_color_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))+labs(x="PGP",y="HQ NRgenome number")
dev.off()
write.csv(CRBD_HQ_NRgenome_PGPR_combine_summarize_presence_Freq_perSource_v2_melt,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBC_vs_Published_PGP_presence_genome_number.csv",row.names = FALSE,quote = FALSE)


## ------ N,P,Fe,IAA,GA,CK,ACC_deaminase,Ethylene and SA (ISR removed)  ----
CRBD_HQ_NRgenome_PGPR_combine<-cbind(CRBD_HQ_NRgenome_nutrient_combine,CRBD_HQ_NRgenome_growth_combine,CRBD_HQ_NRgenome_adaptation_combine)%>%dplyr::rename(Nitrogen=N_fixation)%>%dplyr::rename(Ethylene=EFE)
summary(apply(CRBD_HQ_NRgenome_PGPR_combine,1,function(x) sum(x>0))>=2)
#    Mode   FALSE    TRUE 
# logical    3112    2997 
summary(apply(CRBD_HQ_NRgenome_PGPR_combine,1,function(x) sum(x>0))==0)
#    Mode   FALSE    TRUE 
# logical    5231     878 
combination_matrix = make_comb_mat(CRBD_HQ_NRgenome_PGPR_combine,mode = 'distinct')
write.table(CRBD_HQ_NRgenome_PGPR_combine,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_HQ_NRgenome_PGPR_combine.txt",quote = FALSE,row.names = TRUE)

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/All_PGPR_upset_summary_all.pdf",width = 6,height = 3)
UpSet(combination_matrix,set_order = c("Nitrogen","Phosphorus","Iron","IAA","GA","CK","ACC_deaminase","Ethylene","SA"),comb_col = rev(c("#c9d1b6","#e6edd5","#252829","#2f7782","#86b1b8", "#42b38f","#6da392","#8fa35b","#abb88c"))[comb_degree(combination_matrix)],pt_size = unit(1, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

combination_matrix_v2<-combination_matrix[comb_degree(combination_matrix)>=1]
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/All_PGPR_upset_summary_degree_1plus.pdf",width = 6,height = 3)
UpSet(combination_matrix_v2,set_order = c("Nitrogen","Phosphorus","Iron","IAA","GA","CK","ACC_deaminase","Ethylene","SA"),comb_col = rev(c("#c9d1b6","#e6edd5","#252829","#2f7782","#86b1b8", "#42b38f","#6da392","#8fa35b","#abb88c"))[comb_degree(combination_matrix_v2)],pt_size = unit(1, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

## extract combination labels and its genomeID

Genome_list<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Genome_list[[i]]=extract_comb(combination_matrix_v2,comb)
}

Comb_label<-list()
comb_list<-labels(comb_degree(combination_matrix_v2))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Comb_label[[i]]=rep(comb,each=length(extract_comb(combination_matrix_v2, comb)))
}

combination_genome_composition<-data.frame(GenomeID=unlist(Genome_list),Comb_label=unlist(Comb_label))
combination_genome_composition$Comb_label<-factor(combination_genome_composition$Comb_label,levels = unique(combination_genome_composition$Comb_label))
combination_genome_composition_taxa<-left_join(combination_genome_composition,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
combination_genome_composition_taxa$PhyClass_collapse<-factor(combination_genome_composition_taxa$PhyClass_collapse,levels = PhyClass_order)
combination_genome_composition_taxa$Count=1

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_phyClass_composition_All_PGPR_degree_1plus.pdf",height = 5,width = 10)
ggplot(combination_genome_composition_taxa,aes(x=Comb_label,y=Count,fill=PhyClass_collapse))+geom_bar(stat='identity',position = 'stack',width = 0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

##------Family level ---

combination_genome_Family_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Family)%>%mutate(FamilySum=n())%>%select(c(Comb_label,Family,FamilySum))%>%unique())
combination_genome_Family_count$Family<-gsub("^$","Others",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("f__","",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_G","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_H","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(FamilySum>=40,Family,"Others"))
combination_genome_Family_count$Family_up<-factor(combination_genome_Family_count$Family_up,levels = c(unique(combination_genome_Family_count$Family_up)[unique(combination_genome_Family_count$Family_up)!="Others"],"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_Family_composition_All_PGPR_degree_1plus.pdf",height = 6,width = 10)
ggplot(combination_genome_Family_count,aes(x=Comb_label,y=FamilySum,fill=Family_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values = c(pal_simpsons()(16),pal_jama()(3),"#8ae9eb","grey"))+theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_line(size = 0.1, color = "#e3e3e3"))+guides(fill=guide_legend(ncol=1))+labs(y="NRgenome number")
#scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#d3dbd3'))
dev.off()

## ------ N,P,Fe,IAA,GA,CK,ACC_deaminase,Ethylene and SA at least two PGP----

combination_matrix_v3<-combination_matrix_v2[comb_degree(combination_matrix_v2)>=2]

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/All_PGPR_upset_summary_degress_2plus.pdf",width = 6,height = 3)
UpSet(combination_matrix_v3,set_order = c("Nitrogen","Phosphorus","Iron","IAA","GA","CK","ACC_deaminase","Ethylene","SA"),comb_col = rev(c("#c9d1b6","#e6edd5","#252829","#2f7782","#86b1b8", "#42b38f","#6da392","#8fa35b","#abb88c"))[comb_degree(combination_matrix_v3)],pt_size = unit(1, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

Genome_list<-list()
comb_list<-labels(comb_degree(combination_matrix_v3))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Genome_list[[i]]=extract_comb(combination_matrix_v3,comb)
}

Comb_label<-list()
comb_list<-labels(comb_degree(combination_matrix_v3))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Comb_label[[i]]=rep(comb,each=length(extract_comb(combination_matrix_v3, comb)))
}

combination_genome_composition<-data.frame(GenomeID=unlist(Genome_list),Comb_label=unlist(Comb_label))
combination_genome_composition$Comb_label<-factor(combination_genome_composition$Comb_label,levels = unique(combination_genome_composition$Comb_label))
combination_genome_composition_taxa<-left_join(combination_genome_composition,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
combination_genome_composition_taxa$PhyClass_collapse<-factor(combination_genome_composition_taxa$PhyClass_collapse,levels = PhyClass_order)
combination_genome_composition_taxa$Count=1

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_phyClass_composition_All_PGPR_degree_2plus.pdf",height = 5,width = 10)
ggplot(combination_genome_composition_taxa,aes(x=Comb_label,y=Count,fill=PhyClass_collapse))+geom_bar(stat='identity',position = 'stack',width = 0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

##------Family level ---

combination_genome_Family_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Family)%>%mutate(FamilySum=n())%>%select(c(Comb_label,Family,FamilySum))%>%unique())
combination_genome_Family_count$Family<-gsub("^$","Others",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("f__","",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_G","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_H","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(FamilySum>=20,Family,"Others"))
combination_genome_Family_count$Family_up<-factor(combination_genome_Family_count$Family_up,levels = c(unique(combination_genome_Family_count$Family_up)[unique(combination_genome_Family_count$Family_up)!="Others"],"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_Family_composition_All_PGPR_degree_2plus.pdf",height = 6,width = 10)
ggplot(combination_genome_Family_count,aes(x=Comb_label,y=FamilySum,fill=Family_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values = c(pal_simpsons()(16)[c(2,1,3:16)],pal_jama()(3),"#8ae9eb","grey"))+theme(axis.text.x = element_text(angle = 90),panel.grid.major = element_line(size = 0.1, color = "#e3e3e3"))+guides(fill=guide_legend(ncol=1))+labs(y="NRgenome number")
#scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#d3dbd3'))
dev.off()


# All PGPR category combine -- genus mean RA level

dim(CRBD_HQ_NRgenome_PGPR_combine) # 6109 HQ across 9 Category
CRBD_HQ_NRgenome_PGPR_combine_taxonomy<-left_join(data.frame(CRBD_HQ_NRgenome_PGPR_combine,NRgenomeID=rownames(CRBD_HQ_NRgenome_PGPR_combine)),CRBD_HQ_NRgenome_taxonomy,by=c("NRgenomeID"="NRgenome_GenomeID"))
rownames(CRBD_HQ_NRgenome_PGPR_combine_taxonomy)<-CRBD_HQ_NRgenome_PGPR_combine_taxonomy$NRgenomeID

## ---- CRBD NRgenome source
CRBD_HQ_NRgenome_source_category [1:5,]

## calculate PGPR percentage per Genus

CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up<-inner_join(CRBD_HQ_NRgenome_PGPR_combine_taxonomy,CRBD_HQ_NRgenome_CRBC_vs_Published_ratio,by=c("NRgenomeID"="NRgenome_GenomeID"))
summary(CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up$Genus=="") # 155 CRBD HQ NRgenome does not has genus name
CRBD_HQ_NRgenome_Sum_per_Genus<-CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up%>%filter(Genus!="")%>%select(c(Nitrogen,Phosphorus,Iron,IAA,GA,CK,ACC_deaminase,Ethylene,SA,Genus))%>%group_by(Genus)%>%mutate(Sum_per_Genus=n())%>%select(c(Genus,Sum_per_Genus))%>%unique #443 genus
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus<-data.frame(CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up%>%filter(Genus!="")%>%select(c(Nitrogen,Phosphorus,Iron,IAA,GA,CK,ACC_deaminase,Ethylene,SA,Genus))%>%group_by(Genus)%>%mutate_at(vars(-c(Genus)),function(x) sum(x>0)/sum(x>=0)*100))
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus[is.na(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus)]=0


## extract the best RepSpecies as the representative genome for each Genus

CRBD_HQ_RepSpecies_sub<-data.frame(CRBD_HQ_RepSpecies%>%select(RepSpecies_GenomeID,Genus_gtdb,Score)%>%group_by(Genus_gtdb)%>%mutate(MaxScore=max(Score))%>%filter(Score==MaxScore))  #442 genus with g__ unclassified
summary(CRBD_HQ_RepSpecies_sub$Genus_gtdb=="g__") # 441 genus and 1 g__ unclassified
CRBD_HQ_RepSpecies_sub$Genus_gtdb<-gsub("g__","",CRBD_HQ_RepSpecies_sub$Genus_gtdb)
CRBD_HQ_RepSpecies_sub<-CRBD_HQ_RepSpecies_sub%>%filter(Genus_gtdb!="")%>%dplyr::rename(Genus=Genus_gtdb)
Genus_RepSpecies_map<-CRBD_HQ_RepSpecies_sub%>%select(RepSpecies_GenomeID,Genus)

## subset RepSpecies tree to genus representative ones

## ---- IRBC rep species tree ---

library(ggtree)
library(ape)
CRBD_RepSpecies_tree<-read.tree("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PhyloTree/CRBD_3044.unrooted.tree")
tree_order<-fortify(CRBD_RepSpecies_tree)%>%filter(isTip=="TRUE")%>%arrange(y) # 3044 RepSpecies
RepSpGenome_order_intree<-as.character(tree_order$label) # 3044
RepSpecies_order<-RepSpGenome_order_intree
nodes_to_trim<-RepSpecies_order[!RepSpecies_order%in%Genus_RepSpecies_map$RepSpecies_GenomeID] # 2603
CRBD_HQ_RepGenus_tree<-drop.tip(CRBD_RepSpecies_tree, nodes_to_trim)
length((fortify(CRBD_HQ_RepGenus_tree)%>%filter(isTip=="TRUE")%>%arrange(y))$label) # 441 CRBD HQ RepGenus
HQ_RepGenus_tree_order<-as.character((fortify(CRBD_HQ_RepGenus_tree)%>%filter(isTip=="TRUE")%>%arrange(y))$label) # 441
write.tree(CRBD_HQ_RepGenus_tree,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PhyloTree/CRBD_HQ_441_RepGenus_RepSpecies_unrooted.tree")
write.tree(CRBD_HQ_RepGenus_tree,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot/CRBD_HQ_441_RepGenus_RepSpecies_unrooted.tree")


# Generate the annotation file for RepGenus tree

source("/Users/fangliu/Documents/IGDB_Bai_lab/Script_backup/table2itol/table2itol-master/table2itol.R")
#Apparently this script is running in interactive mode. You could now generate iTOL files by setting some 'infiles' variable to a vector of file names and then calling:create_itol_files(infiles)
setwd("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot")

create_itol_files(infiles = "CRBD_HQ_NRgenome_PGPR_percentage_per_Genus.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")
create_itol_files(infiles = "CRBD_HQ_RepGenus_CRBC_and_Published_present_or_absent.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")
create_itol_files(infiles = "CRBD_HQ_RepGenus_taxonomy.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")

## PGPR percentage
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up<-data.frame(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus%>%filter(Genus%in%Genus_RepSpecies_map$Genus)%>%unique()) #since 443 genera were kept in HQ NRgenome (after remove the g__) and 441 left for CRBD HQ RepSpecies, reason is that RepSpecies were chosed by quality score, and the highest quality score is MQ but no HQ, so, there is some conflict. Be noted, here genus_unclassifed are removed since they does not have information, mixed from different unassigned levels(Phylum, Class all the way to Genus)
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up<-inner_join(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up,Genus_RepSpecies_map)
rownames(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up)<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up$RepSpecies_GenomeID
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up[HQ_RepGenus_tree_order,]
write.csv(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot/CRBD_HQ_NRgenome_PGPR_percentage_per_Genus.csv",row.names = FALSE,quote = FALSE)

## CRBC vs Published has present or absent
#CRBD_HQ_RepGenus_CRBC_vs_Published_presence<-data.frame(inner_join(CRBD_HQ_NRgenome_CRBC_vs_Published_ratio,CRBD_HQ_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,Genus)))%>%group_by(Genus)%>%mutate_at(vars(-c(NRgenome_GenomeID,Genus)),function(x) sum(x))%>%select(!c(NRgenome_GenomeID))%>%unique()) #444 this is not very correctly defined,since source of genomes should not defined based only on high quality genomes, all MQ genomes were used to define the source of each genus
CRBD_RepGenus_CRBC_vs_Published_presence<-data.frame(inner_join(CRBD_NRgenome_CRBC_vs_Published_ratio,CRBD_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,Genus)))%>%group_by(Genus)%>%mutate_at(vars(-c(NRgenome_GenomeID,Genus)),function(x) sum(x))%>%select(!c(NRgenome_GenomeID))%>%unique())  # 444 genera # 6109 HQ NRgenomes

CRBD_HQ_RepGenus_CRBC_vs_Published_presence<-CRBD_RepGenus_CRBC_vs_Published_presence%>%filter(Genus%in%Genus_RepSpecies_map$Genus) # to filter the genus to only those left in Genus_RepSpecies_map
CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up<-inner_join(Genus_RepSpecies_map,CRBD_HQ_RepGenus_CRBC_vs_Published_presence)%>%mutate(CRBC_presence=if_else(CRBC_ratio>0,"Present","Absent"))%>%mutate(Published=if_else(Published_ratio>0,"Present","Absent"))%>%select(!c(CRBC_ratio,Published_ratio))
rownames(CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up)<-CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up$RepSpecies_GenomeID
CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up<-CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up[HQ_RepGenus_tree_order,]
write.csv(CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot/CRBD_HQ_RepGenus_CRBC_and_Published_present_or_absent.csv",row.names = FALSE,quote = FALSE)
CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up<-CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up%>%mutate(Genus_source_cat=case_when(CRBC_presence=="Present"&Published=="Present"~"Both",CRBC_presence=="Present"&Published=="Absent"~"CRBC",CRBC_presence=="Absent"&Published=="Present"~"Published"))
table(CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up$Genus_source_cat) # HQ Genus_source_cat
     # Both      CRBC Published 
     #  129       246        66

## CRBD HQ RepGenus PhyClass
CRBD_HQ_RepGenus_taxonomy<-inner_join(Genus_RepSpecies_map,data.frame(CRBD_HQ_NRgenome_taxonomy%>%select(c(Genus,PhyClass_collapse))%>%unique()))
rownames(CRBD_HQ_RepGenus_taxonomy)<-CRBD_HQ_RepGenus_taxonomy$RepSpecies_GenomeID
CRBD_HQ_RepGenus_taxonomy<-CRBD_HQ_RepGenus_taxonomy[HQ_RepGenus_tree_order,]
write.csv(CRBD_HQ_RepGenus_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot/CRBD_HQ_RepGenus_taxonomy.csv",row.names = FALSE,quote = FALSE)

########### ------- subset to 441 genus -------- ############

create_itol_files(infiles = "CRBD_HQ_NRgenome_PGPR_percentage_per_Genus.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")
create_itol_files(infiles = "CRBD_HQ_RepGenus_CRBC_and_Published_present_or_absent.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")
create_itol_files(infiles = "CRBD_HQ_RepGenus_taxonomy.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")

# All PGPR 9 categories - HQ_NRgenome PGPR prevalence

head(CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up) # 6109 HQ NRgenome
write.csv(CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up_for_SuppFig_3b_v1.csv",row.names = FALSE,quote = FALSE)
CRBD_HQ_NRgenome_PhyClassSUM<-CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up%>%select(c(PhyClass,NRgenomeID))%>%group_by(PhyClass)%>%mutate(PhyClass_SUM=n())%>%select(c(PhyClass,PhyClass_SUM))%>%unique()%>%as.data.frame()
CRBD_HQ_NRgenome_PhyClass_PGP_SUM<-CRBD_HQ_NRgenome_PGPR_combine_taxonomy_up%>%select(Nitrogen:SA,PhyClass)%>%group_by(PhyClass)%>%mutate_at(vars(-PhyClass),function(x) sum(x>0))%>%unique()%>%as.data.frame()
CRBD_HQ_NRgenome_PGP_SUM_PhyClassSUM<-inner_join(CRBD_HQ_NRgenome_PhyClassSUM,CRBD_HQ_NRgenome_PhyClass_PGP_SUM)%>%as.data.frame()
row.names(CRBD_HQ_NRgenome_PGP_SUM_PhyClassSUM)<-CRBD_HQ_NRgenome_PGP_SUM_PhyClassSUM$PhyClass
CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio<-CRBD_HQ_NRgenome_PGP_SUM_PhyClassSUM%>%mutate_at(vars(-c(PhyClass,PhyClass_SUM)),~./PhyClass_SUM*100)%>%arrange(desc(PhyClass_SUM))
CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio$PhyClass<-factor(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio$PhyClass,levels = rev(c(unique(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio$PhyClass)[!unique(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio$PhyClass)%in%rev(c(taxa_color$PhyClass[c(1:6,8:12)],"Deinococcota","Campylobacterota","Nitrospirota","Cyanobacteria","Planctomycetota"))],taxa_color$PhyClass[c(1:6,8:12)],"Deinococcota","Campylobacterota","Nitrospirota","Cyanobacteria","Planctomycetota")))

HQ_NRgenome_PhyClass_SUM<-ggplot(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio,aes(y=PhyClass,x=PhyClass_SUM))+geom_bar(color="grey",fill="#1e5955",stat = 'identity',width = 0.5,size=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 0,hjust = 1))
HQ_NRgenome_PhyClass_SUM_top6<-ggplot(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio%>%filter(PhyClass%in%c(taxa_color$PhyClass[c(1:6)])),aes(y=PhyClass,x=PhyClass_SUM))+geom_bar(color="grey",fill="#1e5955",stat = 'identity',width = 0.5,size=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 0,hjust = 1))

CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt<-melt(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio%>%select(!PhyClass_SUM),id.vars = c("PhyClass"),variable.name = "PGP",value.name = "Percentage")
CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt$PhyClass<-factor(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt$PhyClass,levels = c(taxa_color$PhyClass[c(1:6,8:12)],"Deinococcota","Campylobacterota","Nitrospirota","Cyanobacteria","Planctomycetota"))
PGP_vs_Classes<-data.frame(PGP=c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene"),Classes=c(rep("Nutrient",3),rep("Growth",3),rep("Stress",3)))
CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt<-inner_join(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt,PGP_vs_Classes)
CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt$PGP<-factor(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt$PGP,levels = c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene"))
CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt$Classes<-factor(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt$Classes,levels = c("Nutrient","Growth","Stress"))
#HQ_NRgenome_PGP_PhyClassRatio<-ggplot(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt,aes(x=PhyClass,y=Percentage,fill=PGP,color=PGP))+geom_bar(stat = 'identity',position = "dodge",width = 0.8)+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("#8a8803","#a8a628","#bfbe69","#eb7d00","#f7a548","#ffc98c","#1f8056","#4ba37d","#7bbda0"))+scale_color_manual(values = c("#8a8803","#a8a628","#bfbe69","#eb7d00","#f7a548","#ffc98c","#1f8056","#4ba37d","#7bbda0"))+ggforce::facet_row(vars(Classes))
HQ_NRgenome_PGP_PhyClassRatio<-ggplot(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt,aes(x=PGP,y=Percentage,fill=PhyClass,color=PhyClass))+geom_bar(stat = 'identity',position = "dodge",width = 0.6)+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+scale_color_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+ggforce::facet_row(vars(Classes),space = "free",scales = 'free_x')
HQ_NRgenome_PGP_PhyClassRatio_top6<-ggplot(CRBD_HQ_NRgenome_PGP_SUM_PhyClassRatio_melt%>%filter(PhyClass%in%c(taxa_color$PhyClass[c(1:6)])),aes(x=PGP,y=Percentage,fill=PhyClass,color=PhyClass))+geom_bar(stat = 'identity',position = "dodge",width = 0.6)+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+scale_color_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+ggforce::facet_row(vars(Classes),space = "free",scales = 'free_x')+ylim(0,100)

library(patchwork)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/HQ_NRgenome_PGP_PhyClassRatio.pdf",width = 10,height = 6)
HQ_NRgenome_PGP_PhyClassRatio+HQ_NRgenome_PhyClass_SUM
dev.off()
library(patchwork)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/HQ_NRgenome_PGP_PhyClassRatio_top6.pdf",width = 10,height = 6)
HQ_NRgenome_PGP_PhyClassRatio_top6+HQ_NRgenome_PhyClass_SUM_top6
dev.off()

# All PGPR 9 categories - GenusSUM source dissection


## PGP encoding Freq
CRBD_HQ_NRgenome_PGPR_binary_per_Genus<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up%>%mutate_at(vars(-c(Genus,RepSpecies_GenomeID)),~ ifelse(. != 0, 1, 0))
CRBD_HQ_NRgenome_PGPR_binary_per_Genus_taxonomy<-inner_join(CRBD_HQ_NRgenome_PGPR_binary_per_Genus,CRBD_HQ_NRgenome_taxonomy%>%select(c(Genus,PhyClass))%>%unique()) #441 HQ genus except g__
write.csv(CRBD_HQ_NRgenome_PGPR_binary_per_Genus_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_HQ_NRgenome_PGPR_binary_per_Genus_taxonomy_for_SuppFig_3b_v2.csv",quote = FALSE,row.names = FALSE)
CRBD_HQ_RepGenus_PhyClassSUM<-CRBD_HQ_NRgenome_PGPR_binary_per_Genus_taxonomy%>%select(c(PhyClass,Genus))%>%group_by(PhyClass)%>%mutate(PhyClass_SUM=n())%>%select(c(PhyClass,PhyClass_SUM))%>%unique()%>%as.data.frame()
CRBD_HQ_RepGenus_PhyClass_PGP_SUM<-CRBD_HQ_NRgenome_PGPR_binary_per_Genus_taxonomy%>%select(!Genus)%>%group_by(PhyClass)%>%mutate_at(vars(-PhyClass),function(x) sum(x>0))%>%unique()%>%as.data.frame()
CRBD_HQ_RepGenus_PGP_SUM_PhyClassSUM<-inner_join(CRBD_HQ_RepGenus_PhyClassSUM,CRBD_HQ_RepGenus_PhyClass_PGP_SUM)%>%as.data.frame()
row.names(CRBD_HQ_RepGenus_PGP_SUM_PhyClassSUM)<-CRBD_HQ_RepGenus_PGP_SUM_PhyClassSUM$PhyClass
CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio<-CRBD_HQ_RepGenus_PGP_SUM_PhyClassSUM%>%mutate_at(vars(-c(PhyClass,PhyClass_SUM)),~./PhyClass_SUM*100)%>%arrange(desc(PhyClass_SUM))
CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio$PhyClass<-factor(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio$PhyClass,levels = rev(c(taxa_color$PhyClass[c(1:6,8:12)],"Deinococcota","Campylobacterota","Nitrospirota","Cyanobacteria","Planctomycetota")))

HQ_RepGenus_PhyClass_SUM<-ggplot(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio,aes(y=PhyClass,x=PhyClass_SUM))+geom_bar(color="grey",fill="#1e5955",stat = 'identity',width = 0.5,size=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 0,hjust = 1))
HQ_RepGenus_PhyClass_SUM_top6<-ggplot(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio%>%filter(PhyClass%in%c(taxa_color$PhyClass[c(1:6)])),aes(y=PhyClass,x=PhyClass_SUM))+geom_bar(color="grey",fill="#1e5955",stat = 'identity',width = 0.5,size=0.1)+theme_bw()+theme(axis.text.x = element_text(angle = 0,hjust = 1))

CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt<-melt(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio%>%select(!PhyClass_SUM),id.vars = c("PhyClass"),variable.name = "PGP",value.name = "Percentage")
CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt$PhyClass<-factor(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt$PhyClass,levels = c(taxa_color$PhyClass[c(1:6,8:12)],"Deinococcota","Campylobacterota","Nitrospirota","Cyanobacteria","Planctomycetota"))
PGP_vs_Classes<-data.frame(PGP=c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene"),Classes=c(rep("Nutrient",3),rep("Growth",3),rep("Stress",3)))
CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt<-inner_join(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt,PGP_vs_Classes)
CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt$PGP<-factor(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt$PGP,levels = c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene"))
CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt$Classes<-factor(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt$Classes,levels = c("Nutrient","Growth","Stress"))
#HQ_RepGenus_PGP_PhyClassRatio<-ggplot(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt,aes(x=PhyClass,y=Percentage,fill=PGP,color=PGP))+geom_bar(stat = 'identity',position = "dodge",width = 0.8)+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("#8a8803","#a8a628","#bfbe69","#eb7d00","#f7a548","#ffc98c","#1f8056","#4ba37d","#7bbda0"))+scale_color_manual(values = c("#8a8803","#a8a628","#bfbe69","#eb7d00","#f7a548","#ffc98c","#1f8056","#4ba37d","#7bbda0"))+ggforce::facet_row(vars(Classes))
HQ_RepGenus_PGP_PhyClassRatio<-ggplot(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt,aes(x=PGP,y=Percentage,fill=PhyClass,color=PhyClass))+geom_bar(stat = 'identity',position = "dodge",width = 0.6)+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+scale_color_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+ggforce::facet_row(vars(Classes),space = "free",scales = 'free_x')
HQ_RepGenus_PGP_PhyClassRatio_top6<-ggplot(CRBD_HQ_RepGenus_PGP_SUM_PhyClassRatio_melt%>%filter(PhyClass%in%c(taxa_color$PhyClass[c(1:6)])),aes(x=PGP,y=Percentage,fill=PhyClass,color=PhyClass))+geom_bar(stat = 'identity',position = "dodge",width = 0.6)+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+scale_color_manual(values = c(taxa_color$Color[1:12],"#aba2fa","#43d154","#b6b8b6","#e6e6e6"))+ggforce::facet_row(vars(Classes),space = "free",scales = 'free_x')+ylim(0,100)

library(patchwork)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/HQ_RepGenus_PGP_PhyClassRatio.pdf",width = 10,height = 6)
HQ_RepGenus_PGP_PhyClassRatio+HQ_RepGenus_PhyClass_SUM
dev.off()
library(patchwork)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/HQ_RepGenus_PGP_PhyClassRatio_top6.pdf",width = 10,height = 6)
HQ_RepGenus_PGP_PhyClassRatio_top6+HQ_RepGenus_PhyClass_SUM_top6
dev.off()

## PGP encoding Freq distribution
rownames(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up)<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up$Genus
PGP_encoding_Genus_Freq<-data.frame(Freq=apply(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up%>%select(!c(Genus,RepSpecies_GenomeID)),1,function(x) sum(x>0)))
PGP_encoding_Genus_Freq_SUM<-PGP_encoding_Genus_Freq%>%group_by(Freq)%>%mutate(Freq_SUM=n())%>%unique()%>%as.data.frame()
colnames(PGP_encoding_Genus_Freq_SUM)<-c("PGP_number","Frequency")
PGP_encoding_Genus_Freq%>%filter(Freq>=6)
#                Freq
# Bradyrhizobium    7
# Pseudomonas_E     6
# Burkholderia      6
# Streptomyces      6

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/PGPR_Freq_histogram_at_genus_level.pdf",width = 3,height = 2.5)
p<-ggplot(PGP_encoding_Genus_Freq_SUM,aes(x=PGP_number,y=Frequency))+geom_bar(color="grey",fill="#1e5955",stat = 'identity',width = 0.5,size=0.1)+theme_bw()
data_counts <- ggplot_build(p)$data[[1]]%>%filter(y>0)
p + geom_text(data = data_counts, aes(x = x, y = y, label = y),
              vjust = -0.5, color = "#249c94", size = 2)+labs(x="PGP richness",y="Number of genus")+ylim(0,150)+theme(axis.text = element_text(size=6))
dev.off()

## PGP genus source
head(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up) # 441
head(CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up) #441, CRBC 245, Published 63 and Both 133
table(CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up$Genus_source_cat)
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source<-inner_join(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up%>%select(!RepSpecies_GenomeID),CRBD_HQ_RepGenus_CRBC_vs_Published_presence_up%>%select(c(Genus,Genus_source_cat)))
row.names(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source)<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source$Genus
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy<-inner_join(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source,CRBD_HQ_NRgenome_taxonomy%>%select(!c(NRgenome_GenomeID,Species))%>%unique)
write.csv(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_PGPR_HQ_perGenus_Source_infor.csv",row.names = FALSE,quote = FALSE)
CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source%>%select(!Genus)%>%group_by(Genus_source_cat)%>%mutate_at(vars(-Genus_source_cat),function(x) sum(x>0))%>%unique()%>%as.data.frame()
row.names(CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source)<-CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source$Genus_source_cat
CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up<-CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source%>%select(!Genus_source_cat)%>%t()%>%as.data.frame()%>%mutate(PGP=colnames(CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_up)[1:9])%>%mutate(Classes=c(rep("Nutrient",3),rep("Growth",3),rep("Stress",3)))%>%melt(id.vars=c("PGP","Classes"),value.name = "Genus_SUM",variable.name="Source")
CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up$Classes<-factor(CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up$Classes,levels = c("Nutrient","Growth","Stress"))
CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up$Source<-factor(CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up$Source,levels = c("CRBC","Both","Published"))
CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up$PGP<-factor(CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up$PGP,levels =c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene") )
write.csv(CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBC_vs_Published_PGP_presence_Genus_number.csv",row.names = FALSE,quote = FALSE)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/CRBC_vs_Published_PGP_presence_Genus_number.pdf",width=5,height = 3.5)
ggplot(data=CRBD_HQ_NRgenome_PGPR_Genus_SUM_by_Source_up,aes(x=PGP,y=Genus_SUM,fill=Source,color=Source,group=Source))+geom_bar(stat = 'identity',position = "dodge",width = 0.8)+ggforce::facet_row(~Classes,space = 'free',scales = 'free_x')+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))+scale_color_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))
dev.off()


# All PGPR 9 categories - FamilySUM source dissection

CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy[1:5,] # 441 
CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy$Family<-gsub("f__","",CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy$Family)
CRBD_HQ_RepFamily_list<-(CRBD_genome_meta_up2%>%filter(Genus_gtdb%in%paste("g__",Genus_RepSpecies_map$Genus,sep = ""))%>%filter(Quality_level=="High Quality")%>%select(Family_gtdb)%>%unique()%>%as.data.frame())$Family_gtdb
# 130 families, double checked if direct filtering the HQ family from NRgenome, it also yields of 130 familys, using HQ_repGenus also yields 130 families. 
CRBD_RepFamily_CRBC_vs_Published_presence<-data.frame(inner_join(CRBD_NRgenome_CRBC_vs_Published_ratio,CRBD_NRgenome_taxonomy%>%select(c(NRgenome_GenomeID,Family)))%>%group_by(Family)%>%mutate_at(vars(-c(NRgenome_GenomeID,Family)),function(x) sum(x))%>%select(!c(NRgenome_GenomeID))%>%unique())  # 205-1 (f__)
CRBD_RepFamily_CRBC_vs_Published_presence_up<-CRBD_RepFamily_CRBC_vs_Published_presence%>%mutate(CRBC_presence=if_else(CRBC_ratio>0,"Present","Absent"))%>%mutate(Published=if_else(Published_ratio>0,"Present","Absent"))%>%select(!c(CRBC_ratio,Published_ratio))
rownames(CRBD_RepFamily_CRBC_vs_Published_presence_up)<-CRBD_RepFamily_CRBC_vs_Published_presence_up$Family
CRBD_HQ_RepFamily_CRBC_vs_Published_presence_up<-CRBD_RepFamily_CRBC_vs_Published_presence_up%>%filter(Family%in%CRBD_HQ_RepFamily_list)
CRBD_HQ_Family_source_cat<-CRBD_HQ_RepFamily_CRBC_vs_Published_presence_up%>%mutate(Family_source_cat=case_when(CRBC_presence=="Present"&Published=="Present"~"Both",CRBC_presence=="Present"&Published=="Absent"~"CRBC",CRBC_presence=="Absent"&Published=="Present"~"Published"))
CRBD_HQ_Family_source_cat$Family_up<-gsub("f__","",CRBD_HQ_Family_source_cat$Family)
rownames(CRBD_HQ_Family_source_cat)<-CRBD_HQ_Family_source_cat$Family_up
CRBD_HQ_Family_source_cat<-CRBD_HQ_Family_source_cat%>%select(!Family)%>%dplyr::rename(Family=Family_up)
table(CRBD_HQ_Family_source_cat$Family_source_cat)
#        Both   CRBC    Published 
#         51     70         9 
write.csv(CRBD_HQ_Family_source_cat,"Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_HQ_RepFamily_CRBC_and_Published_present_or_absent.csv",row.names = FALSE,quote = FALSE)
#CRBD_HQ_Family_source_cat<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy%>%group_by(Family,Genus_source_cat)%>%mutate(Family_SUM=n())%>%select(c(Family,Family_SUM))%>%unique()%>%spread(key=Genus_source_cat,value = Family_SUM)%>%mutate_all(~replace(., is.na(.), 0))%>%mutate(Family_source_cat=case_when(Both>0~"Both",Both==0&CRBC>0&Published>0~"Both",Both==0&CRBC>0&Published==0~"CRBC",Both==0&CRBC==0&Published>0~"Published"))%>%select(c(Family,Family_source_cat))%>%as.data.frame() # 130 families, this is not correctly defined
CRBD_HQ_NRgenome_PGPR_Family_SUM<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy%>%select(c(Nitrogen,Phosphorus,Iron,IAA,GA,CK,ACC_deaminase,Ethylene,SA,Family))%>%group_by(Family)%>%mutate_at(vars(-Family),function(x) sum(x>0))%>%unique()%>%as.data.frame()
rownames(CRBD_HQ_NRgenome_PGPR_Family_SUM)<-CRBD_HQ_NRgenome_PGPR_Family_SUM$Family
CRBD_HQ_NRgenome_PGPR_Family_SUM<-CRBD_HQ_NRgenome_PGPR_Family_SUM%>%select(!Family)
CRBD_HQ_NRgenome_PGPR_Family_SUM[CRBD_HQ_NRgenome_PGPR_Family_SUM>0]=1
CRBD_HQ_NRgenome_PGPR_Family_SUM_add_source<-inner_join(data.frame(CRBD_HQ_NRgenome_PGPR_Family_SUM,Family=row.names(CRBD_HQ_NRgenome_PGPR_Family_SUM)),CRBD_HQ_Family_source_cat)
rownames(CRBD_HQ_NRgenome_PGPR_Family_SUM_add_source)<-CRBD_HQ_NRgenome_PGPR_Family_SUM_add_source$Family
write.csv(CRBD_HQ_NRgenome_PGPR_Family_SUM_add_source,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_PGPR_HQ_perFamily_Source_info.csv",row.names = FALSE,quote = FALSE)

CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM<-CRBD_HQ_NRgenome_PGPR_Family_SUM_add_source%>%select(!c(Family,CRBC_presence,Published))%>%group_by(Family_source_cat)%>%mutate_at(vars(-Family_source_cat),function(x) sum(x>0))%>%unique()%>%as.data.frame()
rownames(CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM)<-CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM$Family_source_cat
CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up<-CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM%>%select(!Family_source_cat)%>%t()%>%as.data.frame()%>%mutate(PGP=colnames(CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM)[1:9])%>%mutate(Classes=c(rep("Nutrient",3),rep("Growth",3),rep("Stress",3)))%>%melt(id.vars=c("PGP","Classes"),value.name = "Family_SUM",variable.name="Source")
CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up$Classes<-factor(CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up$Classes,levels = c("Nutrient","Growth","Stress"))
CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up$Source<-factor(CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up$Source,levels = c("CRBC","Both","Published"))
CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up$PGP<-factor(CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up$PGP,levels =c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene") )
write.csv(CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBC_vs_Published_PGP_presence_Family_number.csv",row.names = FALSE,quote = FALSE)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/CRBC_vs_Published_PGP_presence_Family_number.pdf",width=5,height = 3.5)
ggplot(data=CRBD_HQ_NRgenome_PGPR_percentage_per_Family_Source_SUM_up,aes(x=PGP,y=Family_SUM,fill=Source,color=Source,group=Source))+geom_bar(stat = 'identity',position = "dodge",width = 0.8)+ggforce::facet_row(~Classes,space = 'free',scales = 'free_x')+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))+scale_color_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))
dev.off()


# All PGPR 9 categories - PhyClassSUM source dissection

CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy[1:5,] # 441 
CRBD_HQ_PhyClass_source_cat<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy%>%group_by(PhyClass,Genus_source_cat)%>%mutate(PhyClass_SUM=n())%>%select(c(PhyClass,PhyClass_SUM))%>%unique()%>%spread(key=Genus_source_cat,value = PhyClass_SUM)%>%mutate_all(~replace(., is.na(.), 0))%>%mutate(PhyClass_source_cat=case_when(Both>0~"Both",Both==0&CRBC>0&Published>0~"Both",Both==0&CRBC>0&Published==0~"CRBC",Both==0&CRBC==0&Published>0~"Published"))%>%select(c(PhyClass,PhyClass_source_cat))%>%as.data.frame() # 16 PhyClasses

CRBD_HQ_NRgenome_PGPR_PhyClass_SUM<-CRBD_HQ_NRgenome_PGPR_percentage_per_Genus_add_source_taxonomy%>%select(c(Nitrogen,Phosphorus,Iron,IAA,GA,CK,ACC_deaminase,Ethylene,SA,PhyClass))%>%group_by(PhyClass)%>%mutate_at(vars(-PhyClass),function(x) sum(x>0))%>%unique()%>%as.data.frame()
rownames(CRBD_HQ_NRgenome_PGPR_PhyClass_SUM)<-CRBD_HQ_NRgenome_PGPR_PhyClass_SUM$PhyClass
CRBD_HQ_NRgenome_PGPR_PhyClass_SUM<-CRBD_HQ_NRgenome_PGPR_PhyClass_SUM%>%select(!PhyClass)
CRBD_HQ_NRgenome_PGPR_PhyClass_SUM[CRBD_HQ_NRgenome_PGPR_PhyClass_SUM>0]=1
CRBD_HQ_NRgenome_PGPR_PhyClass_SUM_add_source<-inner_join(data.frame(CRBD_HQ_NRgenome_PGPR_PhyClass_SUM,PhyClass=row.names(CRBD_HQ_NRgenome_PGPR_PhyClass_SUM)),CRBD_HQ_PhyClass_source_cat)
rownames(CRBD_HQ_NRgenome_PGPR_PhyClass_SUM_add_source)<-CRBD_HQ_NRgenome_PGPR_PhyClass_SUM_add_source$PhyClass
write.csv(CRBD_HQ_NRgenome_PGPR_PhyClass_SUM_add_source,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_PGPR_HQ_perPhyClass_Source_info.csv",row.names = FALSE,quote = FALSE)
CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM<-CRBD_HQ_NRgenome_PGPR_PhyClass_SUM_add_source%>%select(!PhyClass)%>%group_by(PhyClass_source_cat)%>%mutate_at(vars(-PhyClass_source_cat),function(x) sum(x>0))%>%unique()%>%as.data.frame()
rownames(CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM)<-CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM$PhyClass_source_cat
CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up<-CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM%>%select(!PhyClass_source_cat)%>%t()%>%as.data.frame()%>%mutate(PGP=colnames(CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM)[1:9])%>%mutate(Classes=c(rep("Nutrient",3),rep("Growth",3),rep("Stress",3)))%>%melt(id.vars=c("PGP","Classes"),value.name = "PhyClass_SUM",variable.name="Source")
CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up$Classes<-factor(CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up$Classes,levels = c("Nutrient","Growth","Stress"))
CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up$Source<-factor(CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up$Source,levels = c("CRBC","Both","Published"))
CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up$PGP<-factor(CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up$PGP,levels =c("Phosphorus","Nitrogen","Iron","IAA","GA","CK","ACC_deaminase","SA","Ethylene") )
write.csv(CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBC_vs_Published_PGP_presence_PhyClass_number.csv",row.names = FALSE,quote = FALSE)
pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/CRBC_vs_Published_PGP_presence_PhyClass_number.pdf",width=5,height = 3.5)
ggplot(data=CRBD_HQ_NRgenome_PGPR_percentage_per_PhyClass_Source_SUM_up,aes(x=PGP,y=PhyClass_SUM,fill=Source,color=Source,group=Source))+geom_bar(stat = 'identity',position = "dodge",width = 0.8)+ggforce::facet_row(~Classes,space = 'free',scales = 'free_x')+theme(strip.text = element_text(size=10,colour = "white"),strip.background = element_rect(fill="#557982", colour="grey"),axis.text.x =element_text(angle = 90,hjust = 1))+scale_fill_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))+scale_color_manual(values = c("#61bfbe","#A9C7AF","#F0CEA0"))
dev.off()


# All PGPR three classes -- upset

## ----- upset plot -----
CRBD_HQ_NRgenome_PGPR_combine_3classes<-CRBD_HQ_NRgenome_PGPR_combine%>%mutate(Nutrient=if_else((Nitrogen+Phosphorus+Iron)>0,1,0))%>%mutate(Growth=if_else((IAA+GA+CK)>0,1,0))%>%mutate(Stress=if_else((ACC_deaminase+Ethylene+SA)>0,1,0))%>%select(c(Nutrient,Growth,Stress)) # 4488 nutrient,1795 growth and 2770 Stress
summary(apply(CRBD_HQ_NRgenome_PGPR_combine_3classes,1,function(x) sum(x>0))>=2) # logical    3430    2679 (TRUE)
summary(apply(CRBD_HQ_NRgenome_PGPR_combine_3classes,1,function(x) sum(x>0))==1) # logical    3557 (FALSE)    2552 (TRUE)
summary(apply(CRBD_HQ_NRgenome_PGPR_combine_3classes,1,function(x) sum(x>0))==2) # logical    4573    1536 (TRUE)
summary(apply(CRBD_HQ_NRgenome_PGPR_combine_3classes,1,function(x) sum(x>0))==3) # logical    4966 (FALSE)    1143 (TRUE)
summary(apply(CRBD_HQ_NRgenome_PGPR_combine_3classes,1,function(x) sum(x>0))==0) # logical     5231     878 (TRUE)

write.table(CRBD_HQ_NRgenome_PGPR_combine_3classes,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/R_write_CRBD_HQ_NRgenome_PGPR_combine_3classes.txt",quote = FALSE)
combination_matrix = make_comb_mat(CRBD_HQ_NRgenome_PGPR_combine_3classes,mode = 'distinct')

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/All_PGPR_upset_summary_3classes.pdf",width = 4,height = 3)
UpSet(combination_matrix,set_order = c("Nutrient","Growth","Stress"),comb_col = rev(c("#19635a","#45857e","#6ca69f"))[comb_degree(combination_matrix)],pt_size = unit(2, "mm"),lwd=0.5,bg_pt_col = 'white',bg_col = '#f5f7f7')
dev.off()

##  ---- stacked taxa barplot for PGPR combination - PhyClass -----

## extract combination labels and its genomeID

Genome_list<-list()
comb_list<-labels(comb_degree(combination_matrix))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Genome_list[[i]]=extract_comb(combination_matrix,comb)
}

Comb_label<-list()
comb_list<-labels(comb_degree(combination_matrix))
for (i in 1:length(comb_list)){
  comb<-comb_list[i]
  Comb_label[[i]]=rep(comb,each=length(extract_comb(combination_matrix, comb)))
}

combination_genome_composition<-data.frame(GenomeID=unlist(Genome_list),Comb_label=unlist(Comb_label))
combination_genome_composition$Comb_label<-factor(combination_genome_composition$Comb_label,levels = unique(combination_genome_composition$Comb_label))
combination_genome_composition_taxa<-left_join(combination_genome_composition,CRBD_HQ_NRgenome_taxonomy,by=c("GenomeID"="NRgenome_GenomeID"))
combination_genome_composition_taxa$PhyClass_collapse<-factor(combination_genome_composition_taxa$PhyClass_collapse,levels = PhyClass_order)
combination_genome_composition_taxa$Count=1

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_phyClass_composition_PGPR_3classes.pdf",height = 5,width = 6)
ggplot(combination_genome_composition_taxa,aes(x=Comb_label,y=Count,fill=PhyClass_collapse))+geom_bar(stat='identity',position = 'stack',width = 0.7)+scale_fill_manual(values = taxa_color[PhyClass_order,]$Color)+theme_bw()+theme(axis.text.x = element_text(angle = 90))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Genus -----

combination_genome_Genus_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Genus)%>%mutate(GenusSum=n())%>%select(c(Comb_label,Genus,GenusSum))%>%unique())
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_A","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_B","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_E","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Pseudomonas_F","Pseudomonas",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Bacillus_A","Bacillus",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("Arthrobacter_I","Arthrobacter",combination_genome_Genus_count$Genus)
combination_genome_Genus_count$Genus<-gsub("^$","Others",combination_genome_Genus_count$Genus)
combination_genome_Genus_count<-combination_genome_Genus_count%>%mutate(Genus_up=if_else(GenusSum>=20,Genus,"Others"))
combination_genome_Genus_count$Genus_up<-factor(combination_genome_Genus_count$Genus_up,levels = c(unique(combination_genome_Genus_count$Genus_up)[unique(combination_genome_Genus_count$Genus_up)!="Others"],"Others"))

top_10_Genus<-data.frame()
 for (comb in as.character(unique(combination_genome_Genus_count$Comb_label))){
   i=which(comb==unique(combination_genome_Genus_count$Comb_label))
   sub=combination_genome_Genus_count%>%filter(Comb_label==comb)
   length=if_else(dim(sub)[1]>=5,5,as.numeric(dim(sub)[1]))
   top_10_Genus[1:length,i]=data.frame(sub%>%arrange(desc(GenusSum)))$Genus[1:length]
 }
 top_10_Genus_uniq<-unique(as.character(as.matrix(top_10_Genus)))
combination_genome_Genus_count<-combination_genome_Genus_count%>%mutate(Genus_up=if_else(Genus%in%top_10_Genus_uniq,Genus,"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_Genus_composition_PGPR_3classes.pdf",height = 6,width = 6)
ggplot(combination_genome_Genus_count,aes(x=Comb_label,y=GenusSum,fill=Genus_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values =c('#ff9082','#45b1f9','#c45cc0','#f9f500','#98d66d','#91749e','#f76fc3','#85c1b8','#a584ff','#ffb444','#7ebfe5','#cec0c9','#467584','#005ff9','#8569D5','#ff4c05','#673770','#D14285','#bc8c38','#bcba6d','#b2798d','#235931','#CD9BCD','#00CDCD','#ff1ce0','#CBD588','#f075f4','#f7cdcf','#85c1b8','#ccff60','#0f9032','#DA5724','#009897','#757f77','#aa0012','#121457','#03ffe6','#f7ff03','#ff8400','#c2edc0','#fcfbc0','#fce3d7','#d7dbf7','#ffb444','#4cb7cf','#ff9082','black','blue','#bfba9b','#d3dbd3'))+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=2))
dev.off()

##  ---- stacked taxa barplot for PGPR combination - Family -----

combination_genome_Family_count<-data.frame(combination_genome_composition_taxa%>%group_by(Comb_label,Family)%>%mutate(FamilySum=n())%>%select(c(Comb_label,Family,FamilySum))%>%unique())
combination_genome_Family_count$Family<-gsub("^$","Others",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("f__","",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_G","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count$Family<-gsub("Bacillaceae_H","Bacillaceae",combination_genome_Family_count$Family)
combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(FamilySum>=40,Family,"Others"))
combination_genome_Family_count$Family_up<-factor(combination_genome_Family_count$Family_up,levels = c(unique(combination_genome_Family_count$Family_up)[unique(combination_genome_Family_count$Family_up)!="Others"],"Others"))


# top_10_Family<-data.frame()
# for (comb in as.character(unique(combination_genome_Family_count$Comb_label))){
#   i=which(comb==unique(combination_genome_Family_count$Comb_label))
#   sub=combination_genome_Family_count%>%filter(Comb_label==comb)
#   length=if_else(dim(sub)[1]>=5,5,as.numeric(dim(sub)[1]))
#   top_10_Family[1:length,i]=data.frame(sub%>%arrange(desc(FamilySum)))$Family[1:length]
# }
# top_10_Family_uniq<-unique(as.character(as.matrix(top_10_Family)))
#combination_genome_Family_count<-combination_genome_Family_count%>%mutate(Family_up=if_else(Family%in%top_10_Family_uniq,Family,"Others"))

pdf("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/All_PGPR/Upset_Family_composition_PGPR_3classes_upcolor.pdf",height = 6,width = 8)
ggplot(combination_genome_Family_count,aes(x=Comb_label,y=FamilySum,fill=Family_up))+geom_bar(stat='identity',position = 'stack',width = 0.7)+theme_bw()+scale_fill_manual(values =family_color_meta$Color)+theme(axis.text.x = element_text(angle = 90))+guides(fill=guide_legend(ncol=1))+labs(y="NRgenome number")
#scale_fill_manual(values = c(pal_simpsons()(16),pal_jama()(3),"#8ae9eb","grey"))
Family_color<-data.frame(Family=levels(combination_genome_Family_count$Family_up),Color=c(pal_simpsons()(16),pal_jama()(3),"#8ae9eb","grey"))
write.csv(Family_color,"/Users/fangliu/Documents/IGDB_Bai_lab/Database_and_color_pallete/Color_scheme/Family_color_meta.csv",row.names = FALSE,quote = FALSE)
dev.off()


# All PGPR three classes -- Generate the annotation file

source("/Users/fangliu/Documents/IGDB_Bai_lab/Script_backup/table2itol/table2itol-master/table2itol.R")
#Apparently this script is running in interactive mode. You could now generate iTOL files by setting some 'infiles' variable to a vector of file names and then calling:create_itol_files(infiles)
setwd("/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot")
CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy<-inner_join(data.frame(NRgenomeID=rownames(CRBD_HQ_NRgenome_PGPR_combine_3classes),CRBD_HQ_NRgenome_PGPR_combine_3classes),CRBD_HQ_NRgenome_taxonomy,by=c("NRgenomeID"="NRgenome_GenomeID"))%>%filter(Genus!="")%>%select(Nutrient,Growth,Stress,Genus,NRgenomeID)
rownames(CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy)<-CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy$NRgenomeID
CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy<-CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy%>%select(!NRgenomeID) ## 155 CRBD HQ NRgenome does not has genus name
CRBD_HQ_NRgenome_Sum_per_Genus<-CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy%>%group_by(Genus)%>%mutate(Sum_per_Genus=n())%>%select(c(Genus,Sum_per_Genus))%>%unique #443 genus
CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus<-data.frame(CRBD_HQ_NRgenome_PGPR_combine_3classes_taxonomy%>%group_by(Genus)%>%mutate_at(vars(-c(Genus)),function(x) sum(x>0)/sum(x>=0)*100)%>%unique())
CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus[is.na(CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus)]=0

## PGPR_3classes percentage

CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up<-data.frame(CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus%>%filter(Genus%in%Genus_RepSpecies_map$Genus)%>%unique()) #since 443 genera were kept in HQ NRgenome and 441 left for CRBD HQ RepSpecies
CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up<-inner_join(CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up,Genus_RepSpecies_map)
rownames(CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up)<-CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up$RepSpecies_GenomeID
CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up<-CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up[HQ_RepGenus_tree_order,]
write.csv(CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus_up,"/Users/fangliu/Documents/IGDB_Bai_lab/NRgenome/NRgenome_finalize_06_12_2023/PGPR/Dec_20_CRBD_HQ_NRgenome_PGPR_evalue5id50/Tree_annot/CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus.csv",row.names = FALSE,quote = FALSE)
create_itol_files(infiles = "CRBD_HQ_NRgenome_PGPR_3classes_percentage_per_Genus.csv",identifier = "RepSpecies_GenomeID",label = "RepSpecies_GenomeID",separator = ",")
#----------------------------------

