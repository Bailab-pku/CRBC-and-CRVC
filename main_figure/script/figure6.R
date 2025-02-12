# figure6：
rm(list=ls())

library(tidyverse)
library(reshape2)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(rstatix)
library(ggpubr)
library(ComplexHeatmap)

setwd("~")

# input data:
#-----------------------
# virus RA table
Viral_tpmean_table<-read.table("doc/Viral_tpmean_table.txt",sep="\t",header=T,check.names = F)

Sample_metadata<-read.csv("doc/metagenome_metadata.txt",sep="\t",check.names = F)%>%filter(Dataset!=0)
Sample_metadata$Crop=Sample_metadata$Host
Sample_metadata[Sample_metadata$Host=="Medicago truncatula",]$Crop="Medicago"
Sample_metadata[Sample_metadata$Host=="Oryza sativa L.",]$Crop="Rice"
Sample_metadata[Sample_metadata$Host=="Zea mays L.",]$Crop="Maize"
Sample_metadata[Sample_metadata$Host=="Triticum aestivum L.",]$Crop="Wheat"

Sample_metadata$Crop<-factor(Sample_metadata$Crop, levels=c("Wheat","Rice","Maize","Medicago"))
rownames(Sample_metadata)=Sample_metadata$SampleID_ori

Sample_metadata<-Sample_metadata%>%arrange(Crop, Location)
# defense system RA table
Defensefinder_TPM<-read.table("doc/Defensefinder_TPM.txt",sep="\t",header=T,check.names = F)
Defense_machanism<-read.csv("doc/defense_system_machanism.csv")%>%filter(Mechanism!="Unknown")

# bacteria related tables
CRBD_metadata<-read.table("doc/9772_CRBD_genomes_metadata.txt",header=T,sep="\t",quote="",check.names = F)
CRBD_metadata$Species_ClusterID=gsub("Species","bOTU",CRBD_metadata$Species_ClusterID)
# novel species  info。
novelty_CRBC<-read.table("doc/CRBC_novelty.txt",header=T,sep="\t",check.names = F)

CRBD_abun<-read.table("doc/CRBD_abun.txt",header=T,sep="\t")
CRBD_top_tree<-read.tree("doc/CRBD_330.unrooted.tree")
CRBD_top_tree$tip.label<-CRBD_metadata[match(CRBD_top_tree$tip.label,CRBD_metadata$GenomeID),]$Species_ClusterID

CRBD_CRVC_connection<-read.table("doc/CRBD_CRVC_connection.txt",header=T,sep="\t")

# color:
Color_code<-read.csv("doc/color_code.csv",header = T)

color_code=Color_code$Color
names(color_code)=Color_code$Label
#-----------------------

# 统计prevalence等是按地点而不是dataset合并的，所以需要额外统计一个： 
#-----------------------
# 同一个地点合并：
Sample_metadata$group=Sample_metadata$Dataset
Sample_metadata[Sample_metadata$group%in%c("dataset_7","dataset_8"),]$group="dataset_7&8"
Sample_metadata[Sample_metadata$group%in%c("dataset_10","dataset_11"),]$group="dataset_10&11"
Sample_metadata[Sample_metadata$group%in%c("dataset_1","dataset_2","dataset_3"),]$group="dataset_1&2&3"

Root_metadata<-Sample_metadata%>%filter(Habitat=="Root")
#-------------------------


# figure 6A 
#-----------------------
Phage_tpmean_table<-Viral_tpmean_table%>%filter(`Host domain`=="Bacteria")
CRVC_abun_meta<-colSums(Phage_tpmean_table[,-(1:6)])%>%as.data.frame()%>%rename(tpmean=".")%>%
  mutate(SampleID=rownames(.))%>%merge(Sample_metadata, by="SampleID")

i="tpmean"  
sub_data<-CRVC_abun_meta
stat.test <- sub_data%>%
  group_by(Crop) %>%
  wilcox_test(tpmean~Habitat) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")

stat.test <- stat.test %>%
  add_xy_position(x = "Crop", dodge = 0.8)
stat.test$y.position=1.05*max(sub_data$tpmean)

set.seed(1)
p<-ggplot(sub_data,aes(x=Crop,y=tpmean,color=Habitat))+
  geom_boxplot(aes(fill=Habitat),color='black',alpha=0.7,outlier.colour = NA,lwd=0.1,width=0.7)+
  geom_point(position=position_jitterdodge(jitter.width = 0.2),aes(color=Habitat),size=0.5,alpha=0.5)+
  scale_fill_manual(values=c("#04916e","#754415"))+scale_color_manual(values = c("#04916e","#754415"))+
  theme_bw()+
  #facet_wrap(~Crop,nrow = 1)+
  labs(x="",y="Relative abundance (%) \n", title=paste("all",i,sep=" "))+
  stat_pvalue_manual(stat.test,  label = "p.adj", tip.width=0.1, size = 1.5,color = 'black')+
  theme(text = element_text(color="black",size=7),panel.background = element_blank(),panel.grid = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),strip.background = element_blank(),
        legend.position = "bottom",legend.text = element_text(size=7),legend.title= element_text(size=8),
        plot.title = element_text(hjust=0.5,size=8))

p  

ggsave(paste("figure6/Fig6A_RA_boxplot_",i,"_phage.pdf",sep=""),p,width=3,height = 3)
#-----------------------

# figure 6B
#-----------------------
Defensefinder_TPM_meta<-Defensefinder_TPM%>%melt(id.vars="system",variable.name = "SampleID",value.name = "tpm")%>%
  merge(Defense_machanism%>%select(system, Mechanism)%>%unique(),by="system")%>%
  merge(Sample_metadata%>%select(SampleID, Crop, Dataset,Habitat),by="SampleID")

Defensefinder_system_tpm_median<-Defensefinder_TPM_meta%>%
  group_by(system, Mechanism,Dataset,Crop, Habitat)%>%summarise_if(is.numeric, median)%>%
  group_by(Mechanism,system)%>%summarise_if(is.numeric,max)%>%
  filter(tpm>=0.0001)

dataset_abundance<-Defensefinder_TPM_meta%>%
  filter(system%in%Defensefinder_system_tpm_median$system)%>%
  group_by(system,Habitat,Dataset, Crop)%>%summarise_if(is.numeric,mean)%>%as.data.frame()

min(dataset_abundance[dataset_abundance$tpm>0,]$tpm)

log2_FC_dataset<-dataset_abundance%>%filter(Habitat=="Root")%>%select(system,Crop, Dataset,tpm)%>%
  rename(Root=tpm)%>%merge(dataset_abundance%>%filter(Habitat=="Soil")%>%select(system,Crop, Dataset,tpm)%>%
                                            rename(Soil=tpm),by=c("system", "Dataset","Crop"))%>%
  mutate(LogFC=log10((Root+0.00000001)/(Soil+0.00000001)))

none_0_system<-log2_FC_dataset%>%group_by(system)%>%summarise_if(is.numeric,mean)%>%filter(LogFC!=0)%>%select(system)

log2_FC_dataset$Dataset<-factor(log2_FC_dataset$Dataset,levels = unique(Sample_metadata$Dataset))

# normalize：
log2_FC_dataset_dcast<-dcast(log2_FC_dataset%>%filter(system%in%none_0_system$system),system~Dataset,value.var = "LogFC")
rownames(log2_FC_dataset_dcast)=log2_FC_dataset_dcast$system

log2_FC_dataset_dcast_scale<-log2_FC_dataset_dcast%>%select(!system)%>%t()%>%as.data.frame()%>%
  scale(center=F)%>%t()%>%as.data.frame()

# add annotation：
anno_col<-Sample_metadata%>%select(Dataset,Crop)%>%unique()
rownames(anno_col)=anno_col$Dataset
anno_row<-Defense_machanism%>%filter(system%in%log2_FC_dataset$system)%>%select(system,Mechanism)%>%unique()
rownames(anno_row)=anno_row$system  
anno_row$system<-factor(anno_row$system, levels=rownames(log2_FC_dataset_dcast_scale))
anno_row<-anno_row%>%arrange(system)

bk<-seq(-max(abs(log2_FC_dataset_dcast_scale)),max(abs(log2_FC_dataset_dcast_scale)),length.out=100)

left_anno = HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=c("#F64E60","#61C0BF","#FFC300")),
                                                 labels=c("Abortive Infection","Degrading nucleic acids","DNA synthesis inhibition"),
                                                 labels_gp = gpar(col = "black", fontsize = 7)),
                              show_legend = T, which = "row")

row_split_df<-anno_row%>%select(Mechanism)

top_anno = HeatmapAnnotation(df=anno_col%>%select(Crop),
                             col=list(Crop=c(Medicago="#DF8F44FF",Maize="#6A6599FF",Rice="#00A1D5FF",Wheat="#B24745FF")))

bk2<-seq(min(log2_FC_dataset_dcast_scale),max(log2_FC_dataset_dcast_scale),length.out=100)

pdf("figure6/Fig6B_heatmap_of_defense_system_TPM.pdf",width=10,height = 6)
Heatmap(as.matrix(log2_FC_dataset_dcast_scale),
        row_split = row_split_df, 
        left_annotation = left_anno,
        top_annotation = top_anno,
        border = "black",
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7),
        rect_gp = gpar(col="white",lwd=1),
        column_title = NULL,
        row_title = NULL,
        col = c(colorRampPalette(c("#754415", "#B3A492", "#D6C7AE", "white"))(length(bk2[bk2<=0])),
                colorRampPalette(c("white","#99B080", "#748E63","#183D3D"))(length(bk2[bk2>=0]))),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2"
)

dev.off()
#-----------------------

# figure 6C
#-----------------------
CRBD_top_abun_meta<-CRBD_abun%>%filter(Label=="High abundant")%>%melt(id.vars = c("Representatives","SpeciesID","Group","Label"),variable.name = "SampleID",value.name = "RA")%>%
  merge(Root_metadata, by="SampleID")

Crop_prevalence<-Root_metadata%>%select(group)%>%table()%>%as.data.frame()%>%rename(size=Freq)

# 1. grouped RA
CRBD_top_RA_sep_host<-CRBD_top_abun_meta%>%filter(RA!=0)%>%
  group_by(Representatives, SpeciesID, Crop, group)%>%summarise_if(is.numeric,mean)%>%
  group_by(Representatives, SpeciesID, Crop)%>%summarise_if(is.numeric,mean)%>%
  select(Representatives, SpeciesID, Crop,RA)

# 2. grouped prevalence
CRBD_top_prevalence<-CRBD_top_abun_meta%>%filter(RA>0)%>%mutate(prevalence_count=1)%>%
  group_by(SpeciesID, Group,group)%>%summarise_if(is.numeric,sum)%>%
  merge(Crop_prevalence,by="group")%>%mutate(prevalence=100*prevalence_count/size)%>%
  group_by(SpeciesID, Group)%>%summarise_if(is.numeric,mean)%>%
  select(SpeciesID,prevalence,Group)

## 3. phage species infect these bacteria,89 species：
CRBD_top_phage_species<-CRBD_CRVC_connection%>%filter(Bacterial_SpeciesID%in%CRBD_top_RA_sep_host$SpeciesID)%>%
  select(Bacterial_SpeciesID, Lifestyle)%>%table()%>%as.data.frame()%>%filter(Freq>0)%>%rename(phage_Species_num=Freq)

## CRBD taxonomy and novelty:
CRBD_top_RA<-novelty_CRBC%>%filter(SpeciesID%in%CRBD_top_RA_sep_host$SpeciesID)
CRBD_top_RA$Taxonomy<-ifelse(CRBD_top_RA$`Phylum (gtdb)`=="p__Proteobacteria",CRBD_top_RA$`Class (gtdb)`,CRBD_top_RA$`Phylum (gtdb)`)
CRBD_top_RA$Taxonomy=gsub("p__","", CRBD_top_RA$Taxonomy)
CRBD_top_RA$Taxonomy=gsub("c__","", CRBD_top_RA$Taxonomy)
CRBD_top_RA$`Family (gtdb)`=gsub("f__","", CRBD_top_RA$`Family (gtdb)`)

# top 9 tax:
top_9_tax<-table(CRBD_top_RA$Taxonomy)%>%as.data.frame()%>%arrange(-Freq)%>%head(n=9)
CRBD_top_RA$Taxonomy_adj<-ifelse(CRBD_top_RA$Taxonomy%in%top_9_tax$Var1, CRBD_top_RA$Taxonomy, "Others")

Family_list<-c("Burkholderiaceae","Xanthomonadaceae","Aeromonadaceae","Enterobacteriaceae",
               "Pseudomonadaceae","Rhizobiaceae","Rhizobiaceae_A")

CRBD_top_RA$Family<-CRBD_top_RA$`Family (gtdb)`
idx=CRBD_top_RA$Family%in%Family_list
CRBD_top_RA[!idx,]$Family="not highlight"

# highlight nodes,Enterobacterales order:
highlight_node<-data.frame(nodes=c(388,586,569,563,575,354),
                           Family=c("Rhizobiaceae","Xanthomonadaceae","Aeromonadaceae",
                                    "Enterobacteriaceae", "Pseudomonadaceae","Burkholderiaceae"))

# color family
CRBD_top_tree_meta<-full_join(midpoint.root(CRBD_top_tree),CRBD_top_RA%>%rename(label=SpeciesID)%>%
                                select(label, Family),by="label")

p<-ggtree(midpoint.root(CRBD_top_tree),layout = "fan",open.angle = 15, branch.length = 'none')

tree_plot<-rotate_tree(p,angle=135)+
  geom_highlight(data=highlight_node,aes(node=nodes,fill=Family))+
  geom_fruit(data=CRBD_top_RA%>%select(!Family),geom=geom_point, pwidth = 0.1, 
             mapping=aes(y=SpeciesID, x=0, color=Taxonomy_adj,shape=`Crop root novel species`))+
  geom_fruit(data=CRBD_top_RA_sep_host%>%filter(Crop=="Wheat"),offset = 0.06, color="white",
             geom=geom_tile, pwidth = 0,fill="#980100", mapping=aes(y=SpeciesID, x=1, fill=Crop,alpha=log10(RA)))+
  geom_fruit(data=CRBD_top_RA_sep_host%>%filter(Crop=="Rice"),offset = 0.04, color="white",
             geom=geom_tile, pwidth = 0,fill="#006C90", mapping=aes(y=SpeciesID, x=1, fill=Crop,alpha=log10(RA)))+
  geom_fruit(data=CRBD_top_RA_sep_host%>%filter(Crop=="Maize"),offset = 0.04, color="white",
             geom=geom_tile, pwidth = 0,fill="#3D3976", mapping=aes(y=SpeciesID, x=1, fill=Crop,alpha=log10(RA)))+
  geom_fruit(data=CRBD_top_RA_sep_host%>%filter(Crop=="Medicago"),offset = 0.04, color="white",
             geom=geom_tile, pwidth = 0,fill="#D16200", mapping=aes(y=SpeciesID, x=1, fill=Crop,alpha=log10(RA)))+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5))+
  guides(fill=guide_legend(ncol=2))+
  # prevalence：
  geom_fruit(data=CRBD_top_prevalence,axis.params = list(axis="xy",hjust=1, vjust=1),offset = 0.06,
             geom=geom_bar, pwidth = 0.2,color="white", mapping=aes(y=SpeciesID, x=prevalence,fill=Group),stat="identity")+
  # phage species：
  geom_fruit(data=CRBD_top_phage_species,axis.params = list(axis="y",hjust=1, vjust=1),offset = 0.06,
             geom=geom_point, alpha=0.6,pwidth = 0.1, mapping=aes(y=Bacterial_SpeciesID, x=Lifestyle,size=log10(phage_Species_num),color=Lifestyle))+
  scale_shape_manual(values=c(21,16))+
  guides(color=guide_legend(ncol=2))+
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
                              "Others"= "#cfcecc",
                              "lysogenic"="#527C5A",
                              "lytic"="#E76F51",
                              "unknown"="#cfcecc"))+
  scale_fill_manual(values=c("Rhizobiaceae"='#FADFE8',
                             "Pseudomonadaceae"='#F4F1E2',
                             "Aeromonadaceae"='#E2F4F4',
                             "Enterobacteriaceae"='#EEE2F4',
                             "Xanthomonadaceae"='#E3F4E2',
                             "Burkholderiaceae"='#FEFFBB94',
                             "Multi-crops multizonal"="#7FC97F",
                             "Single-crop multizonal"="#BEAED4",
                             "Single-crop regional"="#FDC086",
                             "Rare"="#FFFF99"))+
  scale_size_continuous(range=c(0.5,3))

tree_plot

ggsave("figure6/Fig6C_high_abun_bac_phage_connection_sep_host.pdf",tree_plot,width=20,height=10)
#-----------------------



