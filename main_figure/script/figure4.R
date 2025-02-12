# figure4：
rm(list=ls())

library(tidyverse)
library(ggVennDiagram)
library(patchwork)
library(ggforce)
setwd("~")

# input data:
#-----------------------
Sample_metadata<-read.csv("doc/metagenome_metadata.txt",sep="\t",check.names = F)%>%filter(Dataset!=0)

# PC and Superpathway: 
Phylum_df<-read.csv("doc/Phylum_table.csv",row.names = 1)%>%mutate(Taxonomy=rownames(.))
Class_df<-read.csv("doc/Class_table.csv",row.names = 1)%>%mutate(Taxonomy=rownames(.))
Superpath_df<-read.csv("doc/SuperPath_table.csv",row.names = 1)%>%mutate(Superpath=rownames(.))

# microbiome related pathways:
Pathway_htext<-read.csv("doc/Total_vs_Bac_KOnum_within_pathway_htext_keep_set.csv",header = T)%>%filter(Keep_or_not=="Keep")
Path_information<-Pathway_htext%>%select(!Pathway)%>%select(!`Keep_or_not`)%>%unique()
Pathway_info<-read.table("doc/pathway_ko.txt",sep='\t',header=T)
Color_code<-read.csv("doc/color_code.csv",header = T)

# pathway conditions:
host_consistent_pathway<-read.csv("doc/four_host_diff_pathway_summary.txt",sep="\t")
host_consistent_ko<-read.table("doc/four_host_diff_KO_summary.txt",header=T,sep="\t")
diff_KO_level<-read.csv("doc/diff_KO_level.txt",sep="\t")

# compare result:
pathway_compare_result<-read.table("doc/compare_Pathway.txt",header=T,sep="\t")
KO_compare_result<-read.table("doc/compare_Pathway.txt",header=T,sep="\t")
#-----------------------

# organize data
#-----------------------
# 1. Crop name change:
Sample_metadata$Crop=Sample_metadata$Host
Sample_metadata[Sample_metadata$Host=="Medicago truncatula",]$Crop="Medicago"
Sample_metadata[Sample_metadata$Host=="Oryza sativa L.",]$Crop="Rice"
Sample_metadata[Sample_metadata$Host=="Zea mays L.",]$Crop="Maize"
Sample_metadata[Sample_metadata$Host=="Triticum aestivum L.",]$Crop="Wheat"

Sample_metadata$Crop<-factor(Sample_metadata$Crop, levels=c("Wheat","Rice","Maize","Medicago"))
rownames(Sample_metadata)=Sample_metadata$SampleID_ori

Sample_metadata<-Sample_metadata%>%arrange(Crop, Location)
Root_metadata<-Sample_metadata%>%filter(Habitat=="Root")

# 2. generate PC table:
proteo_list<-c("Alphaproteobacteria","Gammaproteobacteria")
proteo_sub_class_table<-Class_df%>%filter(Taxonomy%in%proteo_list)

proteo_pc_table<-rbind(proteo_sub_class_table,Phylum_df%>%filter(Taxonomy=="Proteobacteria"))
rownames(proteo_pc_table)=proteo_pc_table$Taxonomy
proteo_pc_table_t<-proteo_pc_table%>%select(!Taxonomy)%>%t()%>%as.data.frame()%>%
  mutate(Proteobacteria_Others=Proteobacteria-(Alphaproteobacteria+Gammaproteobacteria))%>%
  select(!Proteobacteria)%>%t()%>%as.data.frame()%>%mutate(Taxonomy=rownames(.))

PC_df<-rbind(proteo_pc_table_t,Phylum_df%>%filter(Taxonomy!="Proteobacteria"))

# count top 11 abundant tax:
PC_sum_top11<-PC_df%>%mutate(Taxonomy=rownames(.))%>%reshape2::melt(id.vars="Taxonomy")%>%group_by(Taxonomy)%>%
  summarise_if(is.numeric, mean)%>%arrange(-value)%>%filter(Taxonomy!="unassigned" & Taxonomy!="Proteobacteria_Others")%>%head(n=11)

PC_table_top11=PC_df%>%mutate(Taxonomy=rownames(.))
idx=PC_table_top11$Taxonomy%in%PC_sum_top11$Taxonomy

PC_table_top11[!idx,]$Taxonomy="Others"

PC_table_top11<-PC_table_top11%>%group_by(Taxonomy)%>%summarise_if(is.numeric,sum)
PC_table_top11$Taxonomy<-factor(PC_table_top11$Taxonomy, levels=c(PC_sum_top11$Taxonomy,"Others"))

PC_table_top11_meta<-PC_table_top11%>%reshape2::melt(id.vars="Taxonomy")%>%
  rename(SampleID=variable, abundance=value)%>%merge(Root_metadata,by="SampleID")%>%
  group_by(Dataset,Crop,Taxonomy)%>%summarise_if(is.numeric,mean)

# 3. generate superpathway table:
rownames(Superpath_df)=Superpath_df$Superpath
Superpath_df_sub<-Superpath_df%>%filter(Superpath%in%Path_information$Superpathway)
Superpath_df_norm<-1000000*t(t(Superpath_df_sub[,-408])/colSums(Superpath_df_sub[,-408]))%>%as.data.frame()

Superpathway_sum_top11<-Superpath_df_norm%>%mutate(Superpathway=rownames(.))%>%reshape2::melt(id.vars="Superpathway")%>%group_by(Superpathway)%>%
  summarise_if(is.numeric, mean)%>%arrange(-value)%>%head(n=17)

Superpathway_table_top11<-Superpath_df_norm%>%mutate(Superpathway=rownames(.))
idx=Superpathway_table_top11$Superpathway%in%Superpathway_sum_top11$Superpathway

Superpathway_table_top11[!idx,]$Superpathway="Others"

Superpathway_table_top11<-Superpathway_table_top11%>%group_by(Superpathway)%>%summarise_if(is.numeric,sum)
Superpathway_table_top11$Superpathway<-factor(Superpathway_table_top11$Superpathway, levels=c(Superpathway_sum_top11$Superpathway,"Others"))

Superpathway_table_top11_meta<-Superpathway_table_top11%>%reshape2::melt(id.vars="Superpathway")%>%
  rename(SampleID=variable, abundance=value)%>%merge(Root_metadata,by="SampleID")%>%
  group_by(Dataset,Crop,Superpathway)%>%summarise_if(is.numeric,mean)

color_code=Color_code$Color
names(color_code)=Color_code$Label
#-----------------------

# figure 4A & 4B:
#-----------------------
# 4A
p<-ggplot(PC_table_top11_meta,aes(Dataset, abundance, fill=Taxonomy))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = color_code)+
  facet_grid(~Crop,scales = "free_x",)+
  labs(y="Relative abundance (%)",x="")+
  guides(fill=guide_legend(nrow=6))+
  theme_bw()+theme(text = element_text(color="black",size=7),panel.background = element_blank(),panel.grid = element_blank(),
                   strip.background = element_rect(fill=NA),
                   #axis.ticks.x = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1,vjust=1),
                   legend.position = "bottom",legend.text = element_text(size=6),legend.title= element_text(size=8),legend.key.size = unit(2,'mm'),
                   plot.title = element_text(hjust=0.5,size=8))
p


ggsave("figure4/Fig4A_taxonomic_composition_of_crop_root.pdf",p,width=4,height = 3)

# 4B
p<-ggplot(Superpathway_table_top11_meta,aes(Dataset, abundance, fill=Superpathway))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=color_code)+
  facet_grid(~Crop,scales = "free_x")+
  labs(y="Relative abundance (%)",x="")+
  guides(fill=guide_legend(nrow=6))+
  theme_bw()+theme(text = element_text(color="black",size=7),panel.background = element_blank(),panel.grid = element_blank(),
                   strip.background = element_rect(fill=NA),
                   #axis.ticks.x = element_blank(),
                   axis.text.x = element_text(angle = 45,hjust = 1,vjust=1),
                   legend.position = "bottom",legend.text = element_text(size=6),legend.title= element_text(size=8),legend.key.size = unit(2,'mm'),
                   plot.title = element_text(hjust=0.5,size=8))
p

ggsave("figure4/Fig4B_functional_composition_of_crop_root.pdf",p,width=4,height = 3)
#-----------------------

# figure 4C:
#-----------------------
host_consistent_pathway$Crop<-factor(host_consistent_pathway$Crop, levels=c("Wheat","Rice","Maize","Medicago"))

# 4 host co-enriched:
enriched_data<-host_consistent_pathway%>%filter(Level=="Enriched")%>%reshape2::dcast(Pathway~Crop,value.var="Pathway")
rownames(enriched_data)=enriched_data$Pathway
enriched_list<-enriched_data%>%select(!Pathway)%>%as.list()
enriched_list<-lapply(enriched_list, na.omit)

p<-ggVennDiagram(enriched_list, label_size = 4, edge_size = 0.5, edge_lty = "solid",label_alpha = 0)+
  scale_fill_gradient2(low="white", midpoint = 6,high = "#01665E",name = "Pathway_description count")+
  scale_color_manual(values=c(rep("black",4)))+
  ggtitle(paste("Four host root enriched Pathway",sep=" "))+
  theme(legend.position = 'none',plot.title=element_text(hjust=0.5))

p

ggsave(paste("figure4/Fig4C_four_host_root_enriched_Path.pdf",sep=""),p, width=4,height=4)
#-----------------------

# figure 4D & E:
#-----------------------
# 4D
count_of_enriched_path<-host_consistent_pathway%>%filter(Level=="Enriched")%>%select(Pathway)%>%table()%>%as.data.frame()
# boxplot data:
pathway_boxplot_data<-pathway_compare_result%>%select(Pathway,Crop,Dataset,Root,Soil)%>%reshape2::melt(id.vars=c("Pathway","Crop", "Dataset"))%>%rename(Habitat=variable, relative_abundance=value)%>%
  mutate(relative_abundance=relative_abundance/10^4)%>%filter(Pathway%in%count_of_enriched_path[count_of_enriched_path$Freq==4,]$Pathway)

# sort by relative abundance：
pathway_sort<-pathway_boxplot_data%>%filter(Habitat=="Root")%>%group_by(Pathway,Crop)%>%
  summarise_if(is.numeric, mean)%>%group_by(Pathway)%>%summarise_if(is.numeric, mean)%>%arrange(relative_abundance)

pathway_boxplot_data$Pathway<-factor(pathway_boxplot_data$Pathway, levels=pathway_sort$Pathway)
pathway_boxplot_data$Habitat<-factor(pathway_boxplot_data$Habitat,levels = c("Soil","Root"))

p<-ggplot(pathway_boxplot_data,aes(x=Pathway,y=relative_abundance,fill=Habitat,color=Habitat))+
  geom_boxplot(color="black",linewidth=0.1,outlier.color = NA)+
  geom_jitter(size=0.1,width=0.1,height=0)+
  coord_flip()+
  scale_fill_manual(values=c("#754415","#35978F"))+
  scale_color_manual(values=c("#754415","#35978F"))+
  #facet_grid(Level~.,scales="free",space="free")+
  scale_y_continuous(position="right")+
  labs(y="relative abundance (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=6,angle=45,color="black",hjust=0.1),
        axis.title.y = element_blank(),axis.title.x=element_text(size=7),strip.text = element_blank(),
        panel.background = element_blank(),panel.grid = element_blank(),legend.position = "bottom",
        legend.title = element_text(size=7),legend.text = element_text(size=5),legend.key.size = unit(4,'mm'))

p
ggsave(paste("figure4/Fig4D_box_of_Path.pdf",sep=""),p, width=5,height=3)

# 4E:
path_information<-inner_join(Pathway_htext,Pathway_info,by="Pathway")%>%select(!Pathway)%>%select(!`Keep_or_not`)%>%unique()

diff_path_KO<-diff_KO_level%>%merge(path_information,by="KO")%>%rename(Pathway=Pathway_description)%>%
  merge(host_consistent_pathway,by="Pathway")%>%filter(Pathway%in%pathway_boxplot_data$Pathway)%>%
  select(KO,Enriched_level, Pathway,Level)%>%mutate(count=1)%>%group_by(Enriched_level,Level, Pathway)%>%
  summarise_if(is.numeric,sum)%>%as.data.frame()

diff_path_KO$Pathway<-factor(diff_path_KO$Pathway, levels=pathway_sort$Pathway)
diff_path_KO$Enriched_level<-as.factor(diff_path_KO$Enriched_level)

p<-ggplot(diff_path_KO,aes(x=Pathway,y=count,fill=Enriched_level))+
  geom_bar(color="black",linewidth=0.1,stat="identity",position="fill")+
  coord_flip()+
  scale_fill_manual(values = c( "4"="#01665E", "3"="#35978F", "2"="#80CDC1", "1"="#C7EAE5", "0"="white",
                                "-1"="#F6E8C3", "-2"="#DFC27D", "-3"="#BF812D", "-4"="#8C550A"))+
  facet_grid(Level~.,scales="free",space="free")+
  scale_y_continuous(position="right")+
  labs(y="distribution of KOs (%)",fill="Enriched KO frequency in 4 crops")+
  theme_bw()+
  guides(fill=guide_legend(nrow=2,byrow = T))+
  theme(axis.text.x = element_text(size=6,angle=45,color="black",hjust=0.1),axis.text.y = element_blank(),
        axis.title.y = element_blank(),axis.title.x=element_text(size=7),strip.text = element_blank(),
        panel.background = element_blank(),panel.grid = element_blank(),legend.position = "bottom",
        legend.title = element_text(size=7),legend.text = element_text(size=5),legend.key.size = unit(4,'mm'))

p

ggsave(paste("figure4/Fig4E_bar_of_Path.pdf",sep=""),p, width=3,height=3)
#-----------------------

# figure 4F & G:
#-----------------------
conserved_pattern<-diff_KO_level%>%filter(Enriched_level%in%c(-4,4))%>%merge(path_information,by="KO",all.x = T)%>%
  select(KO, Enriched_level, Superpathway)%>%unique()

conserved_pattern[is.na(conserved_pattern$Superpathway),]$Superpathway="Unknown"

# top 10：
conserved_pattern_top8<-conserved_pattern%>%group_by(Enriched_level,Superpathway)%>%mutate(count=1)%>%summarise_if(is.numeric,sum)%>%
  as.data.frame()%>%arrange(-count)%>%filter(Superpathway!="Unknown")

conserved_root_enriched<-conserved_pattern_top8%>%filter(Enriched_level==4)
conserved_root_enriched[9:dim(conserved_root_enriched)[1],]$Superpathway="Others"

conserved_root_enriched<-conserved_root_enriched%>%group_by(Superpathway)%>%summarise_if(is.numeric,sum)%>%arrange(-count)
conserved_root_enriched$Superpathway<-factor(conserved_root_enriched$Superpathway, levels=c(conserved_root_enriched[conserved_root_enriched$Superpathway!="Others",]$Superpathway,"Others"))
conserved_root_enriched<-conserved_root_enriched%>%arrange(Superpathway)

num_enrich=length(unique(conserved_pattern[conserved_pattern$Enriched_level==4 & conserved_pattern$Superpathway!="Unknown",]$KO))

p<-ggplot(conserved_root_enriched)+
  geom_arc_bar(stat="pie",aes(x0=0,y0=0,r0=1,r=2,amount=count,fill=Superpathway))+
  annotate("text", x = 0, y = 0, label = num_enrich, size = 6)+
  labs(title="root consistently enriched KOs")+
  scale_fill_manual(values = color_code)+
  theme_minimal()+
  theme(axis.title = element_blank(),axis.text = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),legend.text = element_text(size=7),
        legend.key.size = unit(5,'mm'))

p

ggsave("figure4/Fig4F_enrich_KO_superpathway_circle.pdf",p,width=6,height = 4)

conserved_root_depleted<-conserved_pattern_top8%>%filter(Enriched_level==-4)
conserved_root_depleted[9:dim(conserved_root_depleted)[1],]$Superpathway="Others"

conserved_root_depleted<-conserved_root_depleted%>%group_by(Superpathway)%>%summarise_if(is.numeric,sum)%>%arrange(-count)
conserved_root_depleted$Superpathway<-factor(conserved_root_depleted$Superpathway, levels=c(conserved_root_depleted[conserved_root_depleted$Superpathway!="Others",]$Superpathway,"Others"))
conserved_root_depleted<-conserved_root_depleted%>%arrange(Superpathway)

num_deplete=length(unique(conserved_pattern[conserved_pattern$Enriched_level==-4 & conserved_pattern$Superpathway!="Unknown",]$KO))

p<-ggplot(conserved_root_depleted)+
  geom_arc_bar(stat="pie",aes(x0=0,y0=0,r0=1,r=2,amount=count,fill=Superpathway))+
  annotate("text", x = 0, y = 0, label = num_deplete, size = 6)+
  labs(title="root consistently depleted KOs")+
  scale_fill_manual(values = color_code)+
  theme_minimal()+
  #coord_polar(theta = "y")+
  theme(axis.title = element_blank(),axis.text = element_blank(),panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5,size=12),legend.text = element_text(size=7),
        legend.key.size = unit(5,'mm'))

p

ggsave("figure4/Fig4G_deplete_KO_superpathway_circle.pdf",p,width=6,height = 4)
#-----------------------
