# figure5：
rm(list=ls())

library(tidyverse)
library(ggVennDiagram)
library(reshape2)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(multcompView)
library(PMCMRplus)

setwd("~")

# input data:
#-----------------------
CRBD_metadata<-read.table("doc/9772_CRBD_genomes_metadata.txt",header=T,sep="\t",quote="",check.names = F)
CRBD_name<-read.table("doc/CRBD_GenomeID.txt",header=T,sep="\t")
CRBD_name$GenomeID_standard<-gsub("Pub_","",CRBD_name$GenomeID_standard)
CRBD_metadata<-CRBD_metadata%>%merge(CRBD_name,by="GenomeID")

CRVC_metadata<-read.table("doc/CRVC_metadata.txt",header=T,sep="\t",check.names = F)
connect_origin<-read.table("doc/CRVC_predicted_phage_blast_result_SpeciesLevelSummary.txt",sep="\t",header=T)
CRVC_host_virus_connection_species<-read.table("doc/CRVC_host_virus_connection_species.txt",header=T,sep="\t",check.names = F)%>%
  merge(connect_origin%>%select(contigID, Viral_GenusID)%>%unique(),by="contigID")

PC_cluster_similarity<-read.table("doc/rep_Genus_viral_protein_id50_cov80_genome_PC_jaccard_dcast.tsv",
                                                                  header=T, row.names = 1,sep="\t")

CRVC_novelty<-read.table("doc/CRVC_Genus_novelty.txt",header=T,sep="\t",check.names = F)

virus_pi_table<-read.table("doc/virus_pi_table.txt",header=T,sep="\t",check.names = F)

Sample_metadata<-read.csv("doc/metagenome_metadata.txt",sep="\t",check.names = F)%>%filter(Dataset!=0)
Sample_metadata$Crop=Sample_metadata$Host
Sample_metadata[Sample_metadata$Host=="Medicago truncatula",]$Crop="Medicago"
Sample_metadata[Sample_metadata$Host=="Oryza sativa L.",]$Crop="Rice"
Sample_metadata[Sample_metadata$Host=="Zea mays L.",]$Crop="Maize"
Sample_metadata[Sample_metadata$Host=="Triticum aestivum L.",]$Crop="Wheat"

Sample_metadata$Crop<-factor(Sample_metadata$Crop, levels=c("Wheat","Rice","Maize","Medicago"))
rownames(Sample_metadata)=Sample_metadata$SampleID_ori

Sample_metadata<-Sample_metadata%>%arrange(Crop, Location)
Root_metadata<-Sample_metadata%>%filter(Habitat=="Root")

Color_code<-read.csv("doc/color_code.csv",header = T)

color_code=Color_code$Color
names(color_code)=Color_code$Label
#-----------------------

# Figure 5A:
#-----------------------
# 1. host range:
Genus_rep_genome<-CRVC_metadata%>%filter(`Viral genome`%in%rownames(PC_cluster_similarity))%>%select(GenusID, `Viral genome`)%>%rename(Viral_GenusID=GenusID)

CRBD_metadata$Taxonomy=ifelse(CRBD_metadata$`Phylum (gtdb)`=="p__Proteobacteria",CRBD_metadata$`Class (gtdb)`,CRBD_metadata$`Phylum (gtdb)`)
CRBD_metadata$Taxonomy=gsub("p__","", CRBD_metadata$Taxonomy)
CRBD_metadata$Taxonomy=gsub("c__","", CRBD_metadata$Taxonomy)

CRVC_host_virus_connection_tax<-CRVC_host_virus_connection_species%>%merge(CRBD_metadata%>%select(GenomeID_standard,Taxonomy)%>%rename(`Host bacteria`=GenomeID_standard),by="Host bacteria")

CRVC_host_virus_connection_size<-CRVC_host_virus_connection_tax%>%select(`Viral genome matched`, Viral_GenusID)%>%unique()%>%select(Viral_GenusID)%>%table()%>%as.data.frame()

CRVC_host_range<-CRVC_host_virus_connection_tax%>%select(`Viral genome matched`, Viral_GenusID, Taxonomy)%>%unique()%>%mutate(count=1)%>%
  group_by(Viral_GenusID, Taxonomy)%>%summarise_if(is.numeric,sum)%>%merge(CRVC_host_virus_connection_size,by="Viral_GenusID")%>%
  mutate(consistency=count/Freq)%>%filter(consistency>=0.7)%>%merge(Genus_rep_genome,by="Viral_GenusID")

# 2. temperate rate:
CRVC_genus_temperate<-CRVC_metadata%>%filter(Lifestyle=="lysogenic")%>%select(GenusID)%>%table()%>%as.data.frame()%>%
  merge(CRVC_novelty%>%select(GenusID, `Genus size`),by="GenusID",all.y=T)%>%mutate(temperate_pro=Freq/`Genus size`)%>%
  rename(Viral_GenusID=GenusID)%>%merge(Genus_rep_genome,by="Viral_GenusID")
CRVC_genus_temperate[is.na(CRVC_genus_temperate)]=0

# 3. Genus size
CRVC_Genus_size<-CRVC_novelty%>%select(GenusID, `Genus size`)%>%rename(Viral_GenusID=GenusID)%>%
  mutate(log_Size=log10(`Genus size`))%>%merge(Genus_rep_genome,by="Viral_GenusID")

# 4. novelty：
CRVC_novelty_summary<-CRVC_novelty%>%select(!`Genus size`)%>%melt(id.vars="GenusID")%>%filter(value=="Yes")%>%
  select(GenusID)%>%table()%>%as.data.frame()%>%filter(Freq==5)%>%rename(Viral_GenusID=GenusID)%>%merge(Genus_rep_genome,by="Viral_GenusID",all.y=T)
CRVC_novelty_summary[is.na(CRVC_novelty_summary)]=0

CRVC_novelty_summary$novelty=ifelse(CRVC_novelty_summary$Freq==5,"Yes","no")
CRVC_novelty_summary$novelty<-factor(CRVC_novelty_summary$novelty, levels=c('Yes',"no"))
CRVC_novelty_summary<-CRVC_novelty_summary%>%rename(novelty_code=Freq)

top_host=names(table(CRVC_host_range$Taxonomy)%>%sort()%>%tail(n=7))
CRVC_host_range$Taxonomy=ifelse(CRVC_host_range$Taxonomy%in%top_host,CRVC_host_range$Taxonomy, "Others")
CRVC_host_range$Taxonomy<-factor(CRVC_host_range$Taxonomy,levels = c(rev(top_host), "Others"))

# 2917 tree：
idx=colSums(PC_cluster_similarity)==dim(PC_cluster_similarity)[1]-1

PC_cluster<-hclust(as.dist(PC_cluster_similarity[!idx,!idx]),method = "average")

PC_tree<-as.phylo(PC_cluster)
PC_tree<-midpoint.root(PC_tree)

p<-ggtree(PC_tree,size=0.1, layout = "fan",open.angle = 180,branch.length = 'none')+
  geom_fruit(data=CRVC_host_range%>%select(`Viral genome`, Taxonomy),axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"),
             geom=geom_bar,pwidth=0.1, mapping=aes(y=`Viral genome`, x=0.1, fill=Taxonomy),stat="identity")+
  scale_fill_manual(values=color_code,guide='none')+
  new_scale_fill()+
  geom_fruit(data=CRVC_genus_temperate%>%select(`Viral genome`, temperate_pro),axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"),
             geom=geom_bar, pwidth = 0.1, mapping=aes(y=`Viral genome`, x=0.1, fill=temperate_pro),stat="identity")+
  scale_fill_gradient2(low="white",mid = "#D4C9A8", high = "#527C5A",midpoint = 0.5,guide='none')+
  new_scale_fill()+
  geom_fruit(data=CRVC_Genus_size%>%select(`Viral genome`,log_Size),axis.params = list(axis="xy",hjust=1, vjust=1,text.size=2,line.color="black"),grid.params = list(vline=F,size=0.1,linetype=2,nb),
             geom=geom_bar,pwidth=0.3, mapping=aes(y=`Viral genome`, x=log_Size),stat="identity")+
  new_scale_fill()+
  geom_fruit(data=CRVC_novelty_summary%>%select(`Viral genome`,novelty),offset = 0,
             axis.params = list(axis="y",hjust=1, vjust=1,line.color="black"),
             geom=geom_bar,pwidth=0.1, mapping=aes(y=`Viral genome`, x=1, fill=novelty),stat="identity")+
  scale_fill_manual(values = c("Yes"="#E09220","no"="white"))+
  theme(legend.position = "bottom")

p

ggsave("figure5/Fig5A_PC_phage_genus_tree_jaccard_2917.pdf",p,width=20,height=10)
#-----------------------

# Figure 5B:
#-----------------------
CRVC_genus<-CRVC_metadata%>%select(GenusID)%>%unique()%>%rename(ClusterID=GenusID)
CRVC_genus$ClusterID<-gsub("Genus","Cluster",CRVC_genus$ClusterID)
Refseq<-read.table("doc/Refseq_rep_cluster.tsv",header=T)
MGV<-read.table("doc/MGV_rep_cluster.tsv",header=T)
GOV2<-read.table("doc/GOV2_rep_cluster.tsv",header=T)
IMG_Root<-read.table("doc/IMG_Root_rep_cluster.tsv",header=T)
IMG_Soil<-read.table("doc/IMG_Soil_rep_cluster.tsv",header=T)

# compare environment:
env_genus<-list(CRVC_3097=CRVC_genus$ClusterID,
                Refseq_3615=Refseq$ClusterID,
                MGV_5697=MGV$ClusterID,
                GOV2_8441=GOV2$ClusterID,
                IMG_Root_and_Soil_29375=union(IMG_Soil$ClusterID,IMG_Root$ClusterID))

p<-ggVennDiagram(env_genus, label_size = 4, edge_size = 0.5, label = "count")+
  scale_fill_gradient(low="white",high = "white",name = "species overlap")+
  scale_color_manual(values=c("red","black","black","black","black"))+
  labs(title="Genus level clusters")+
  theme(panel.background = element_blank(),panel.grid = element_blank(),
        strip.background = element_rect(fill=NA),legend.position = 'none',
        legend.text = element_text(size=7),legend.title= element_text(size=8),
        plot.title = element_text(hjust=0.5,size=8))

p

ggsave("figure5/Fig5B_venn_plot_genus_level.pdf",p,width=4,height=4)
#-----------------------

# Figure 5C:
#-----------------------
Viral_pi_melt<-virus_pi_table%>%
  melt(id.vars=c("contigID","SpeciesID","Host domain","Lifestyle" , "Group","Representatives"),variable.name = "SampleID_ori",value.name = "pi")

Viral_pi_melt$Group=factor(Viral_pi_melt$Group, levels=c("Multi_crops_multizonal","Single_crop_multizonal","Single_crop_regional","rare"))
# filter exist pi
Viral_pi_melt<-na.omit(Viral_pi_melt)

CRVC_pi<-Viral_pi_melt%>%merge(Root_metadata%>%select(SampleID,Crop,Dataset)%>%rename(SampleID_ori=SampleID),by="SampleID_ori")

CRVC_pi_grouped<-CRVC_pi%>%group_by(SampleID_ori,Group, Crop, Dataset)%>%summarise_if(is.numeric,mean)

all_compare_result<-data.frame(stringsAsFactors = F)

for (i in unique(CRVC_pi_grouped$Crop)) {
  sub_data<-CRVC_pi_grouped%>%filter(Crop==i)
  multi_test<-kwAllPairsDunnTest(pi~Group, 
                                 data=sub_data,p.adjust.method="fdr")
  
  multi_test_pvalue<-as.data.frame(multi_test$p.value)%>%mutate(group1=rownames(.))%>%
    reshape2::melt(id.var="group1",value.name = "p")%>%mutate(group=paste(group1,variable,sep="-"))%>%
    dplyr::select(group, p)
  
  multi_test_pvalue<-na.omit(multi_test_pvalue)
  
  dif<-c(multi_test_pvalue$p)
  names(dif)<-c(multi_test_pvalue$group)
  # label
  cld <- multcompLetters(dif,reversed = T)
  catalog_sig<-cld$Letters%>%as.data.frame()%>%rename(Letters=".")%>%mutate(group=rownames(.))%>%mutate(Crop=i)
  catalog_max<-sub_data%>%group_by(Group)%>%summarise_if(is.numeric, max)%>%select(Group, pi)
  
  catalog_sig<-catalog_sig%>%rename(Group=group)%>%merge(catalog_max, by="Group")
  
  
  all_compare_result<-rbind(all_compare_result,catalog_sig)
}

CRVC_pi_grouped$Crop<-factor(CRVC_pi_grouped$Crop, levels=c("Wheat","Rice","Maize","Medicago"))
all_compare_result$Crop<-factor(all_compare_result$Crop, levels=c("Wheat","Rice","Maize","Medicago"))
p<-ggplot(CRVC_pi_grouped,aes(Group, pi,fill=Group))+
  geom_jitter(width=0.1,size=0.5,aes(color=Group),alpha=0.5)+
  geom_boxplot(outlier.colour = NA,linewidth=0.1,width=0.7)+
  geom_text(data=all_compare_result,aes(Group, pi,label=Letters),size=3)+
  labs(y="Microdiversity (average pi)",x="all root virus")+
  scale_fill_brewer(palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
  facet_wrap(~Crop,nrow=1)+theme_bw()+
  theme(text = element_text(color="black",size=7),panel.background = element_blank(),
        panel.grid = element_blank(),strip.background  = element_rect(fill=NA),
        legend.position = "bottom",legend.text = element_text(size=7),legend.title= element_text(size=8),axis.text.x = element_blank(), 
        plot.title = element_text(hjust=0.5,size=8))

p

ggsave("figure5/Fig5C_Samples_average_pi.pdf", p, width = 5, height = 4)
#-----------------------