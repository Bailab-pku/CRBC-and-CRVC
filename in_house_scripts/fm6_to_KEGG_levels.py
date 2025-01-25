#!/usr/bin/env python3
# @Author : Fang Liu
# @Email: fliu21@genetics.ac.cn
# @File : fm6_to_KEGG_levels.py

import pandas as pd
#import dask.dataframe as dd
import sys
from optparse import OptionParser

parser=OptionParser()
parser.add_option("--KEGG_blast_fm6",dest="KEGG_blast_fm6",help="This is the blast output file in format 6",metavar="input file")
parser.add_option("--genomeID",dest="genomeID",help="This genomeID will be used as the header after data transformation",metavar="genomeID")
parser.add_option("--Kegg_gene_to_KO_map",dest="Gene_to_KO",help="This is the map file between KEGG gene and KEGG KO",metavar="inpu0t file")
parser.add_option("--Kegg_KO_to_pathway",dest="Kegg_KO_to_Pathway",help="This is the hierarchical file of KEGG KO",metavar="input file")
parser.add_option("--Kegg_pathway_hierarchical_file",dest="Pathway_htext",help="This is the hierarchical file of KEGG pathways",metavar="input file")
parser.add_option("--KO_output",dest="KO_count",help="This is the KO count output after join count table with Gene_to_KO map file",metavar="output file")
parser.add_option("--Pathway_output",dest="Pathway_count",help="This is the Pathway count output after join KO count table with KEGG hierarchical file",metavar="output file")
parser.add_option("--SuperPath_output",dest="SuperPathway_count",help="This is the count table on top of KEGG pathway",metavar="output file")
parser.add_option("--Global_catagory_output",dest="Global_catagory_count",help="This is the count table on top of the top catagory of KEGG pathway",metavar="output file")
parser.add_option("--evalue",dest="evalue",help="This is the evalue cut off set to extract NRgenes that belong to specific pathway,such as [1.0e-30]",type='float')
parser.add_option("--identity",dest="identity",help="This is the identity threshold that is set to exclude annotation iterms that below this cutoff, such as [95.00]",type='float')

(options, args) = parser.parse_args()
KEGG_blast_fm6=options.KEGG_blast_fm6
genomeID=options.genomeID
gene_KO_map=options.Gene_to_KO
KO_path_map=options.Kegg_KO_to_Pathway
path_htext=options.Pathway_htext
evalue=options.evalue
identity=options.identity

KO_out_filename=options.KO_count
Pathway_out_filename=options.Pathway_count
SuperPathway_out_filename=options.SuperPathway_count
Global_catagory_out_filename=options.Global_catagory_count

# read in kegg fm6 file
fm6=pd.read_csv(KEGG_blast_fm6,sep="\t",header=None,names=["NRgeneID","KEGG_gene","Identity","Align_len","Mismatches","Gaps","q_start","q_end","s_start","s_end","evalue","Bit_score"])
fm6_up1=fm6.loc[(fm6["Identity"]>=identity)&(fm6["evalue"]<=evalue),["NRgeneID","KEGG_gene"]]
fm6_up1["Count"]=1
fm6_up2=fm6_up1.iloc[:,[1,2]]
gene_KO_map=pd.read_csv(gene_KO_map,sep="\t",header=0)
KO_path_map=pd.read_csv(KO_path_map,sep="\t",header=0)
path_htext=pd.read_csv(path_htext,sep="\t",header=0)

# inner merge KEGG mapping file with genome fm6
KO_df=pd.merge(fm6_up2,gene_KO_map,left_on="KEGG_gene",right_on="KeggGeneID",how="inner").groupby(["KO"]).sum()
KO_df.columns.values[0]=str(genomeID)
KO_df.to_csv(str(KO_out_filename),sep="\t",index=True)
#KO_df["KO"]=KO_df.index
ko_df=pd.merge(KO_df,KO_path_map,on="KO").groupby(["Pathway"]).sum()
Hierarchy_df=pd.merge(ko_df,path_htext,on="Pathway",how="inner")
Pathway_df=Hierarchy_df.loc[:,["Pathway_description",str(genomeID)]]
Pathway_df.columns=["Pathway",str(genomeID)]
Pathway_df_up=Pathway_df.groupby(["Pathway"]).sum()# updated as three biofilm need to be further merged and sum 
Pathway_df_up.to_csv(str(Pathway_out_filename),sep="\t",index=True) #updated
Superpath_df=Hierarchy_df.groupby(["Superpathway"]).sum()
Superpath_df.columns.values[0]=str(genomeID)
Superpath_df.to_csv(str(SuperPathway_out_filename),sep="\t",index=True)
Global_df=Hierarchy_df.groupby(["Global_catagory"]).sum()
Global_df.columns.values[0]=str(genomeID)
Global_df.to_csv(str(Global_catagory_out_filename),index=True,sep="\t")

