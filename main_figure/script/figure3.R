# Figure3：
rm(list=ls())

library(patchwork)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(tidyverse)
library(picante)
library(phytools)
library(reshape2)
library(RColorBrewer)
library(ComplexUpset)

setwd("~")

# input data:
#-----------------------
CRBD_metadata<-read.table("doc/9772_CRBD_genomes_metadata.txt",header=T,sep="\t",quote="",check.names = F)
CRBD_name_changed<-read.table("doc/CRBD_GenomeID.txt",header=T,sep="\t")
CRBD_metadata$Genome_ClusterID<-gsub("Genome","Gen", CRBD_metadata$Genome_ClusterID)
CRBD_metadata$Species_ClusterID<-gsub("Species","bOTU", CRBD_metadata$Species_ClusterID)
CRBD_metadata<-CRBD_metadata%>%merge(CRBD_name_changed,by="GenomeID")

PGP_of_genomes<-read.csv("doc/HQ_NRgenome_PGPR_combine_taxonomy.csv",check.names = F)

BGC_table<-read.csv("doc/BGC_table.txt",sep="\t",check.names = F)
GCF_table<-read.csv("doc/GCF_table.txt",header=T,sep="\t",check.names = F)

# genus level phylogenetic tree：
CRBD_genus_tree<-read.tree("doc/CRBD_HQ_441_RepGenus_RepSpecies_unrooted.tree")
g<-ggtree(CRBD_genus_tree)

tree_data <- get_taxa_name(g)
color_code<-read.csv("doc/color_code.csv",header = T)
#-----------------------

# organize data：
#-----------------------
# filter HQ NR genomes:
HQ_NR_CRBD_meta<-CRBD_metadata%>%filter(GenomeID%in%PGP_of_genomes$NRgenomeID)%>%mutate(Genus=gsub("g__","",`Genus (gtdb)`))
Genus_Rep_genome<-HQ_NR_CRBD_meta%>%filter(GenomeID%in%CRBD_genus_tree$tip.label)
Genus_Rep_genome$GenomeID<-factor(Genus_Rep_genome$GenomeID, levels=c(tree_data))
Genus_Rep_genome<-Genus_Rep_genome%>%arrange(GenomeID)

Genus_size=HQ_NR_CRBD_meta%>%select(Genus)%>%table()%>%as.data.frame()%>%filter(Genus%in%Genus_Rep_genome$Genus)

PGP_of_genomes$Family=gsub("f__","",PGP_of_genomes$Family)
#-----------------------

#  1.plot tree:
#-----------------------
Genus_PGP<-PGP_of_genomes%>%mutate(Genus_size=CRBC_ratio+Published_ratio)%>%group_by(PhyClass,Genus)%>%
  summarise_if(is.numeric,sum)%>%filter(Genus%in%Genus_Rep_genome$Genus)%>%
  mutate(across(where(is.numeric),~./Genus_size))%>%as.data.frame()%>%
  merge(Genus_Rep_genome%>%select(GenomeID, Genus),by="Genus")

Genus_PGP$PhyClass=ifelse(Genus_PGP$PhyClass=="Deinococcota","Others",
                          ifelse(Genus_PGP$PhyClass%in%color_code$Label,Genus_PGP$PhyClass, "Others"))

Genus_PGP$Genus<-factor(Genus_PGP$Genus, levels=Genus_Rep_genome$Genus)
Genus_PGP<-Genus_PGP%>%arrange(Genus)
rownames(Genus_PGP)=Genus_PGP$GenomeID

# Genus group:
Genus_PGP$Group=ifelse(Genus_PGP$CRBC_ratio!=0,"CRBC genomes","Published genomes")

Genus_gh=Genus_PGP[,3:11]
Genus_gh_melt<-Genus_PGP[,c(3:11,15)]%>%reshape2::melt(id.var="GenomeID",variable.name="Type",value.name = "Prop")
Genus_gh_melt$Type<-factor(Genus_gh_melt$Type,levels=c("Phosphorus","Nitrogen","Iron",
                                                       "IAA","GA","CK",
                                                       "ACC_deaminase","SA","Ethylene"))

## BGC comp:
BGC_table_each_genome<-BGC_table%>%mutate(SUM=1)%>%group_by(GenomeID)%>%summarise_if(is.numeric,sum)%>%select(GenomeID,SUM)

Genus_BGC<-BGC_table%>%merge(HQ_NR_CRBD_meta%>%select(GenomeID, Genus),by="GenomeID")%>%
  merge(BGC_table_each_genome,by="GenomeID")%>%mutate(Count=1,Comp=Count/SUM)%>%
  group_by(GenomeID,Genus,`BiG SCAPE class`)%>%summarise_if(is.numeric,sum)%>%
  group_by(Genus,`BiG SCAPE class`)%>%summarise_if(is.numeric,mean)%>%select(Genus, `BiG SCAPE class`,Comp)%>%
  merge(Genus_Rep_genome%>%select(GenomeID,Genus),by="Genus")%>%as.data.frame()

Genus_BGC_norm<-Genus_BGC%>%merge(Genus_BGC%>%group_by(GenomeID)%>%summarise_if(is.numeric,sum)%>%rename(SUM=Comp),by="GenomeID")%>%
  mutate(Comp_mean=Comp/SUM)

Genus_BGC_norm$`BiG SCAPE class`=factor(Genus_BGC_norm$`BiG SCAPE class`, levels=c("NRPS","PKS-NRP_Hybrids","PKSI","PKSother",
                                                                               "RiPPs","Saccharides","Terpene", "Others"))

# GCF number median
GCF_median<-BGC_table%>%select(GenomeID, `GCF id`)%>%unique()%>%select(GenomeID)%>%table()%>%
  as.data.frame()%>%rename(SUM=Freq)%>%merge(HQ_NR_CRBD_meta%>%select(GenomeID,Genus),by="GenomeID")%>%
  group_by(Genus)%>%summarise_if(is.numeric,median)%>%
  merge(Genus_Rep_genome%>%select(GenomeID,Genus),by="Genus")%>%as.data.frame()

CRBD_genus_tree_meta<-full_join(CRBD_genus_tree,Genus_PGP%>%rename(label=GenomeID)%>%select(label,PhyClass),by="label")

p<-ggtree(CRBD_genus_tree_meta)+layout_fan(angle=180)+
  geom_aline(aes(color=PhyClass),linetype = "solid",size=1)+
  geom_fruit(data=Genus_PGP%>%select(GenomeID,Group),geom=geom_bar,
             stat="identity",aes(y=GenomeID,x=1,fill=Group),pwidth = 0.1,offset = 0.01)+
  scale_color_manual(values=c("Alphaproteobacteria"="#FFA500",
                              "Actinobacteriota"="#8bc24c",
                              "Gammaproteobacteria"= "#FFD700",
                              "Bacteroidota"= "#2196F3",
                              "Firmicutes"= "#8A4B08",
                              "Myxococcota"= "#d9534f",
                              "Patescibacteria"="#86519D",
                              "Spirochaetota"= "#EE4C97",
                              "Chloroflexota"= "#3B3B98",
                              "Verrucomicrobiota"= "#2B5A41",
                              "Fibrobacterota"="#b4bb72",
                              "Acidobacteriota"="#b23256",
                              "Others"= "#cfcecc"))+
  scale_fill_manual(values=c(`CRBC genomes` = "#61bfbe",
                             `Published genomes`="#F0CEA0"), guide = "none")+
  new_scale_fill()+
  geom_fruit(data=Genus_gh_melt%>%rename(label=GenomeID)%>%filter(Type%in%c("Phosphorus","Nitrogen","Iron")),geom = geom_tile,
             aes(y=label,x=Type, fill=Prop),axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"),pwidth = 0.2,offset = 0.05)+
  geom_fruit(data=Genus_gh_melt%>%rename(label=GenomeID)%>%filter(Type%in%c("IAA","GA","CK")),geom = geom_tile,
             aes(y=label,x=Type, fill=Prop),axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"),pwidth = 0.2,offset = 0.05)+
  geom_fruit(data=Genus_gh_melt%>%rename(label=GenomeID)%>%filter(Type%in%c("ACC_deaminase","SA","Ethylene")),geom = geom_tile,
             aes(y=label,x=Type, fill=Prop),axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"),pwidth = 0.2,offset = 0.05)+
  scale_fill_gradientn(colours=colorRampPalette(c("white","#edaa53", "#388E3C"))(100), guide = "none")+
  new_scale_fill()+
  geom_fruit(data=Genus_BGC_norm%>%rename(label=GenomeID),geom=geom_bar,aes(y=label,x=Comp_mean,fill=`BiG SCAPE class`),
             stat="identity",width=0.8,pwidth = 1,offset = 0.1,axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"))+
  geom_fruit(data=GCF_median%>%rename(label=GenomeID),geom=geom_bar,aes(y=label,x=SUM),
             stat="identity",width=0.8, fill="#B44916", pwidth = 0.5,offset = 0.02)+
  scale_fill_manual(values=c("NRPS"="#E6996C",
                             "PKS-NRP-Hybrids"="#DCADA9",
                             "PKSI"="#C3817D",
                             "PKSother"="#AF605D",
                             "RiPPs"="#D8CF71",
                             "Saccharides"="#D4AED0",
                             "Terpene"="#599F69",
                             "Others"="#68A7B8"))+
  geom_treescale(linesize = 3,fontsize = 20,width=1)

p

ggsave("figure3/figure3A_Genus_level_functional_tree.pdf", p, width=10,height = 10,limitsize = F)
#-----------------------



















