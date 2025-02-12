# figure2：
rm(list=ls())

library(patchwork)
library(ggtree)
library(ggtreeExtra)
library(tidyverse)
library(picante)
library(phytools)
library(multcompView)
library(PMCMRplus)

setwd("~")
# input data:
#-----------------------
CRBD_metadata<-read.table("doc/9772_CRBD_genomes_metadata.txt",header=T,sep="\t",quote="",check.names = F)
CRBD_tree<-read.tree("doc/CRBD_3044.unrooted.tree")
# novel species info
novelty_CRBC<-read.table("doc/CRBC_novelty.txt",header=T,sep="\t",check.names = F)
# color code:
color_code<-read.csv("doc/color_code.csv",header=T)
#----------------------

# organize data:
#----------------------
CRBD_Species_metadata<-CRBD_metadata%>%filter(Rep_NR_Species=="Yes")%>%
  merge(novelty_CRBC%>%select(GenomeID, `Novel to GTDB`),by="GenomeID")%>%rename(label=GenomeID)%>%mutate(taxonomy=`Phylum (gtdb)`)
CRBD_Species_metadata[CRBD_Species_metadata$taxonomy=="p__Proteobacteria",]$taxonomy=CRBD_Species_metadata[CRBD_Species_metadata$taxonomy=="p__Proteobacteria",]$`Class (gtdb)`

CRBD_Species_metadata$taxonomy=gsub("p__","",CRBD_Species_metadata$taxonomy)
CRBD_Species_metadata$taxonomy=gsub("c__","",CRBD_Species_metadata$taxonomy)

# Keep only top 11 abundant taxa:
Tax_freq<-CRBD_Species_metadata%>%select(taxonomy)%>%table()%>%as.data.frame()%>%arrange(-Freq)

# 2992 species
idx=CRBD_Species_metadata$taxonomy%in%c(head(Tax_freq$taxonomy,n=12)) & CRBD_Species_metadata$taxonomy!=""
CRBD_Species_metadata[!idx,]$taxonomy="Others"

Tax_freq<-CRBD_Species_metadata%>%select(taxonomy)%>%table()%>%as.data.frame()%>%arrange(-Freq)
Tax_level<-c(as.character(Tax_freq[Tax_freq$taxonomy!="Others",]$taxonomy), "Others")

# order of taxa
CRBD_Species_metadata$taxonomy<-factor(CRBD_Species_metadata$taxonomy, levels = Tax_level)

# other information
Species_size<-CRBD_metadata%>%filter(Rep_NR_Genomes=="Yes")%>%select(Species_ClusterID)%>%table()%>%as.data.frame()%>%rename(species_size=Freq)
CRBD_Species_metadata<-CRBD_Species_metadata%>%merge(Species_size, by="Species_ClusterID")%>%mutate(Genome_group=novel_to_CRBD_Pub)

CRBD_Species_metadata[CRBD_Species_metadata$Genome_group=="",]$Genome_group="published isolates and MAGs"
CRBD_Species_metadata[CRBD_Species_metadata$Genome_group=="Yes" & CRBD_Species_metadata$`Genome type`=="Isolate",]$Genome_group="CRBC novel isolates"
CRBD_Species_metadata[CRBD_Species_metadata$Genome_group=="Yes" & CRBD_Species_metadata$`Genome type`=="MAG",]$Genome_group="CRBC novel MAGs"

CRBD_Species_metadata[CRBD_Species_metadata==""]="no"
# 2318 CRBC species.
CRBC_species<-CRBD_metadata%>%filter(Source=="CRBC")%>%select(Species_ClusterID)%>%unique()
CRBD_Species_annotation<-CRBD_Species_metadata%>%select(label, Genome_group, `Quality level`,Species_ClusterID,taxonomy, species_size, novel_to_Plant_root, `Novel to GTDB`)%>%
  mutate(CRBC="no")

CRBD_Species_annotation[CRBD_Species_annotation$Species_ClusterID%in%CRBC_species$Species_ClusterID,]$CRBC="Yes"
rownames(CRBD_Species_annotation)=CRBD_Species_annotation$label
#----------------------

# 1. tree plot:
#------------------
CRBD_tree_meta<-full_join(CRBD_tree,CRBD_Species_annotation%>%select(label,Genome_group),by="label")

tree_plot<-ggtree(CRBD_tree_meta,layout = "fan",open.angle = 15)+
  geom_aline(aes(color=Genome_group),linetype = "solid",size=0.2)+
  scale_color_manual(values=c("#61bfbe", "#96ceb4","#F0CEA0"))

description_plot<-CRBD_Species_annotation%>%dplyr::select(label,taxonomy,`Quality level`)%>%
  mutate(count=1)

species_size_plot<-CRBD_Species_annotation%>%select(label, species_size)%>%mutate(log10size=log10(species_size))

novelty_plot<-CRBD_Species_annotation%>%select(label,`Novel to GTDB`)%>%mutate(count=1)
novelty_plot[novelty_plot$`Novel to GTDB`=="Yes",]$`Novel to GTDB`="novel_GTDB"

tax_color=color_code$Color
names(tax_color)=color_code$Label
  
p<-rotate_tree(tree_plot,angle=20)+
  geom_treescale(linesize = 3,fontsize = 20,width=0.4)+
  geom_fruit(data=description_plot,geom=geom_bar,
             stat="identity",width=0.7,
             aes(y=label,x=count,fill=taxonomy),pwidth=0.05)+
  geom_fruit(data=species_size_plot,geom=geom_bar,
             stat="identity",width=0.7,
             aes(y=label,x=log10size),pwidth=0.25)+
  geom_fruit(data=novelty_plot,geom=geom_bar,
             stat="identity",width=0.7,
             aes(y=label,x=count,fill=`Novel to GTDB`),pwidth=0.05)+
  scale_fill_manual(values=tax_color)+theme(legend.position = 'none')

ggsave("figure2/Fig2A_Crop_root_bac_3044_tree.pdf", p, width=100,height = 100,limitsize = F)
#------------------

# 2. tree phylogenetic diversity calculate:
#------------------
## trans to rooted tree
CRBD_tree_rooted<-multi2di(CRBD_tree)
### generate richness table：
mimic_table<-reshape2::dcast(CRBD_metadata%>%mutate(Group=paste(Source, `Genome type`,sep="_")), Species_ClusterID~Group)%>%
  mutate(CRBD_Pub=Published_Isolate+Published_MAG,all=CRBC_Isolate+CRBC_MAG+CRBD_Pub)%>%
  merge(CRBD_Species_annotation%>%select(label, Species_ClusterID),by="Species_ClusterID")%>%
  select(-c("Published_Isolate", "Published_MAG","Species_ClusterID" ))

rownames(mimic_table)=mimic_table$label
mimic_table=mimic_table%>%select(-label)

# a. calculate total diversity
mimic_table_final<-mimic_table%>%mutate_all(as.numeric)%>%
  mutate(all_noMAG=all-CRBC_MAG)%>%as.data.frame()

mimic_table_t=t(mimic_table_final)

PD_table=pd(mimic_table_t,CRBD_tree_rooted,include.root = T)
PD_table$PD=round(PD_table$PD,digits = 2)

PD_table_t<-data.frame(t(PD_table),check.names = F)

PD_table_t<-PD_table_t%>%mutate(CRBC_gain=all-CRBD_Pub,
                                CRBC_Isolate_gain=all_noMAG-CRBD_Pub,
                                CRBC_MAG_gain=all-all_noMAG)

PD_table_2<-data.frame(t(PD_table_t))
PD_table_2$PD_total=PD_table_2[rownames(PD_table_2)=="all",]$PD
PD_table_2$proportion<-100*PD_table_2$PD/PD_table_2$PD_total

PD_table_2$group=rownames(PD_table_2)
PD_table_2

PD_table_2$group=gsub("CRBC_Isolate","CRBC isolates",PD_table_2$group)
PD_table_2$group=gsub("CRBC_MAG","CRBC MAGs",PD_table_2$group)
PD_table_2$group=gsub("CRBD_Pub","Published genomes",PD_table_2$group)
PD_table_2$group=gsub("_"," ",PD_table_2$group)

# b. for venn  phylo:
mimic_table<-reshape2::dcast(CRBD_metadata%>%mutate(Group=paste(Source, `Genome type`,sep="_")), Species_ClusterID~Group)%>%
  mutate(CRBD_Pub=Published_Isolate+Published_MAG,all=CRBC_Isolate+CRBC_MAG+CRBD_Pub,
         Iso_MAG=CRBC_Isolate+CRBC_MAG,Iso_Pub=CRBC_Isolate+CRBD_Pub,MAG_Pub=CRBC_MAG+CRBD_Pub
  )%>%
  merge(CRBD_Species_annotation%>%select(label, Species_ClusterID),by="Species_ClusterID")%>%
  select(-c("Published_Isolate", "Published_MAG","Species_ClusterID" ))

rownames(mimic_table)=mimic_table$label
mimic_table=mimic_table%>%select(-label)

mimic_table_t=t(mimic_table)

PD_table=pd(mimic_table_t,CRBD_tree_rooted,include.root = T)
PD_table$PD=round(PD_table$PD,digits = 2)

PD_table_t<-data.frame(t(PD_table),check.names = F)
PD_table_t

PD_table_cal<-PD_table_t%>%mutate(Pub_Uni=all-Iso_MAG,Iso_Uni=all-MAG_Pub,MAG_Uni=all-Iso_Pub,
                                  Pub_Iso_Uni=MAG_Pub-CRBC_MAG-Pub_Uni, Pub_MAG_Uni=Iso_Pub-CRBC_Isolate-Pub_Uni,Iso_MAG_Uni=Iso_Pub-CRBD_Pub-Iso_Uni,
                                  All_Shared=CRBD_Pub-Pub_Uni-Pub_Iso_Uni-Pub_MAG_Uni)

PD_table_2<-data.frame(t(PD_table_cal))
PD_table_2$PD_total=PD_table_2[rownames(PD_table_2)=="all",]$PD
PD_table_2$proportion<-100*PD_table_2$PD/PD_table_2$PD_total

PD_table_2$group=rownames(PD_table_2)
PD_table_2

PD_table_2_final<-data.frame(t(PD_table_2%>%select(!group)))


venn.plot <- VennDiagram::draw.triple.venn(
  area1 = round(PD_table_2_final[4,]$CRBD_Pub,1),
  area2 = round(PD_table_2_final[4,]$CRBC_Isolate,1),
  area3 = round(PD_table_2_final[4,]$CRBC_MAG,1),
  n12 = round(PD_table_2_final[4,]$Pub_Iso_Uni+PD_table_2_final[4,]$All_Shared,1),
  n23 = round(PD_table_2_final[4,]$Iso_MAG_Uni+PD_table_2_final[4,]$All_Shared,1),
  n13 = round(PD_table_2_final[4,]$Pub_MAG_Uni+PD_table_2_final[4,]$All_Shared,1),
  n123 =round(PD_table_2_final[4,]$All_Shared,1),
  category = c("Published", "CRBC isolates", "CRBC MAGs"),
  cex = 1.5,
  cat.cex = 1.5,
  cat.fontfamily = "sans", 
  fontfamily = "sans"
)

ggsave("figure2/Fig2B_venn_plot_of_CRBD_PD.pdf",venn.plot,width = 4,height = 4)

## c. taxonomy Grouped bar plot
PC_taxonomy<-CRBD_Species_metadata
PC_taxonomy$taxonomy=PC_taxonomy$`Phylum (gtdb)`
PC_taxonomy[PC_taxonomy$taxonomy=="p__Proteobacteria",]$taxonomy=PC_taxonomy[PC_taxonomy$taxonomy=="p__Proteobacteria",]$`Class (gtdb)`

PC_taxonomy$taxonomy=gsub("p__","",PC_taxonomy$taxonomy)
PC_taxonomy$taxonomy=gsub("c__","",PC_taxonomy$taxonomy)
rownames(PC_taxonomy)=PC_taxonomy$label

mimic_table_tax<-mimic_table_final%>%merge(PC_taxonomy%>%select(taxonomy),by="row.names")
mimic_table_tax[mimic_table_tax$CRBC_Isolate!=0,]$CRBC_Isolate=1
mimic_table_tax[mimic_table_tax$CRBC_MAG!=0,]$CRBC_MAG=1
mimic_table_tax[mimic_table_tax$CRBD_Pub!=0,]$CRBD_Pub=1
mimic_table_tax[mimic_table_tax$all!=0,]$all=1
mimic_table_tax[mimic_table_tax$all_noMAG!=0,]$all_noMAG=1

## PC
mimic_table_PC=mimic_table_tax
rownames(mimic_table_PC)=mimic_table_PC$Row.names
i="Alphaproteobacteria"
data=mimic_table_PC%>%filter(taxonomy==i)%>%
  select(all,CRBD_Pub,CRBC_Isolate,CRBC_MAG,all_noMAG)%>%
  t()%>%as.data.frame()

PD_table=pd(data,CRBD_tree_rooted,include.root = T)%>%
  mutate(PD=round(PD,digits = 2), taxonomy=c(i))
PD_table$group=rownames(PD_table)

for (i in unique(mimic_table_PC$taxonomy)) {
  data=mimic_table_PC%>%filter(taxonomy==i)%>%
    select(all,CRBD_Pub,CRBC_Isolate,CRBC_MAG,all_noMAG)%>%
    t()%>%as.data.frame()
  
  PD_table_temp=pd(data,CRBD_tree_rooted,include.root = T)%>%
    mutate(PD=round(PD,digits = 2), taxonomy=c(i))
  PD_table_temp$group=rownames(PD_table_temp)
  
  PD_table=rbind(PD_table,PD_table_temp)
  
}

PD_table<-PD_table%>%unique()%>%as.data.frame()

PD_table$group=gsub("CRBC_Isolate","CRBC isolates",PD_table$group)
PD_table$group=gsub("CRBC_MAG","CRBC MAGs",PD_table$group)
PD_table$group=gsub("CRBD_Pub","Published genomes",PD_table$group)
PD_table$group=gsub("all_noMAG","Published and CRBC isolates",PD_table$group)
PD_table$group=gsub("_"," ",PD_table$group)

PD_table_dcast<-reshape2::dcast(PD_table, taxonomy~group, value.var = "PD")
PD_table_dcast<-PD_table_dcast%>%mutate(`CRBC gain`=all- `Published genomes`,
                                        `CRBC isolates gain`=`Published and CRBC isolates`-`Published genomes`,
                                        `CRBC MAGs gain`=all-`Published and CRBC isolates`)

Species_PC_sum<-mimic_table_PC%>%
  group_by(taxonomy)%>%summarise_if(is.numeric, sum)%>%mutate(taxonomy_info=paste(taxonomy," (",all,")",sep=""))%>%
  arrange(-all)%>%select(taxonomy,taxonomy_info)

plot_PG<-PD_table_dcast%>%select(taxonomy,`Published genomes`, `CRBC isolates gain`, `CRBC MAGs gain`)%>%
  reshape2::melt(id.vars="taxonomy")%>%
  merge(Species_PC_sum, by="taxonomy")
colnames(plot_PG)=c("taxonomy","group", "PG","taxonomy_info" )
plot_PG<-plot_PG%>%mutate(PD_proportion=100*PG/max(PD_table_2$PD))

plot_PG$taxonomy_info=factor(plot_PG$taxonomy_info, levels=c(Species_PC_sum$taxonomy_info))

plot_PG$group=factor(plot_PG$group, levels=c("CRBC MAGs gain","CRBC isolates gain", "Published genomes"))

PG_label<-PD_table%>%filter(group=="all")%>%
  mutate(PD_proportion=100*PD/max(PD_table_2$PD),taxonomy_info=paste(taxonomy,"(",SR,")",sep=""))

### only top 20
taxonomy_list<-Species_PC_sum%>%filter(taxonomy!="")%>%head(n=20)%>%select(taxonomy_info)

plot_PG2<-plot_PG%>%filter(taxonomy_info%in%taxonomy_list$taxonomy_info)
plot_PG2$taxonomy_info<-factor(plot_PG2$taxonomy_info, levels = c(taxonomy_list$taxonomy_info))

PG_label2<-PG_label%>%filter(taxonomy_info%in%taxonomy_list$taxonomy_info)
PG_label2$taxonomy_info<-factor(PG_label2$taxonomy_info, levels = c(taxonomy_list$taxonomy_info))


p<-ggplot(plot_PG2,aes(taxonomy_info, PD_proportion))+
  scale_fill_manual(values = c( "#96ceb4", "#61bfbe","#F0CEA0"))+
  geom_bar(aes(fill=group),alpha=0.8,stat="identity",position = "stack",color="black")+
  labs(y="phylogenetic diversity contribution at each phylum (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),axis.text=element_text(size=7,color = "black"),
        axis.title.x = element_blank(),legend.position = "bottom")

p

ggsave("figure2/Fig2B_PC_PG_distribution_top20_left_gtdb.pdf",p,width = 3.5,height=6)
#------------------

# 3. genome coverage:
#------------------
catalog_coverage<-read.table("doc/genome_coverage.txt",header=T,sep="\t")

## data_out_of_genome_catdalog:
catalog_42_coverage<-catalog_coverage%>%filter(Note=="Extra 42")
catalog_42_coverage_melt<-reshape2::melt(catalog_42_coverage[,1:9], id.var='SampleID_ori')
colnames(catalog_42_coverage_melt)=c("SampleID", "group", "mapping_ratio" )

## sort the group
catalog_42_coverage_melt_group<-catalog_42_coverage_melt%>%group_by(group)%>%summarise_if(is.numeric,mean)%>%arrange(-mapping_ratio)
catalog_42_coverage_melt$group<-factor(catalog_42_coverage_melt$group, levels=catalog_42_coverage_melt_group$group)

## kruskal test & dunn test
kruskal<-kruskal.test(mapping_ratio~group, data=catalog_42_coverage_melt)
kruskal

multi_test<-kwAllPairsDunnTest(mapping_ratio~group, 
                               data=catalog_42_coverage_melt,p.adjust.method="fdr")
multi_test

multi_test_pvalue<-as.data.frame(multi_test$p.value)%>%mutate(group1=rownames(.))%>%
  reshape2::melt(id.var="group1")%>%mutate(group=paste(group1,variable,sep="-"))%>%
  dplyr::select(group, value)%>%filter(value>0)

dif<-c(multi_test_pvalue$value)
names(dif)<-c(multi_test_pvalue$group)

# label the difference
cld <- multcompLetters(dif)

catalog_sig<-cld$Letters%>%as.data.frame()%>%rename(Letters=".")%>%mutate(group=rownames(.))
catalog_max<-catalog_42_coverage_melt%>%group_by(group)%>%summarise_if(is.numeric, max)

catalog_sig<-catalog_sig%>%merge(catalog_max, by="group")

catalog_42_coverage_melt$group<-factor(catalog_42_coverage_melt$group, levels=c("NCBI","GTDB","GEM","UHGG","OMD",
                                                                                "CRBD_Pub","CRBC","CRBD"))

p<-ggplot(catalog_42_coverage_melt, aes(group,mapping_ratio, fill=group, color=group))+
  geom_boxplot(alpha=1, color="black", outlier.colour = NA,linewidth=0.5)+
  labs(y="Percentage of root microbiomes coverd by genome reference (%)")+
  geom_jitter(height=0,width = 0.2,alpha=0.6,size=1)+
  geom_text(data=catalog_sig,aes(label=Letters), color="black")+
  scale_fill_manual(values = c("#666633", "#CC9900", "#3AA646","#f17d80", "#737495", "#F0CEA0", "#61bfbe","#428584"))+
  scale_color_manual(values = c("#666633","#CC9900","#3AA646","#f17d80", "#737495", "#F0CEA0", "#61bfbe","#428584"))+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text = element_text(size=7,color="black"),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),legend.position = "none")

p

ggsave("figure2/FigC_mapping_ratio_of_datasets_outCatalog.pdf",p, width=4,height = 5.5)
#------------------

