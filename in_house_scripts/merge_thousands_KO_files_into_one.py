#!/usr/bin/env python3
# @Author : Fang Liu
# @Email: fliu21@genetics.ac.cn
# @File : merge_thousands_KO_files_into_one.py

import pandas as pd
import os
import functools as ft
import sys
import fnmatch as fn

from optparse import OptionParser

parser=OptionParser()


parser.add_option("--input_dir",dest="input_dir",help="This is the input directory for files that you want to merge into one",metavar="input file")
parser.add_option("--pattern",dest="pattern",help="The pattern of input files that you want to merge together ",metavar="input file pattern")
parser.add_option("--how",dest="how",help="This used to define how you want you files merged,such as ['left','right','outer','inner','cross']",metavar="merge type")
parser.add_option("--key",dest="key",help="This is used to define the key, on which column the files were merged together",metavar="parameter")
parser.add_option("--output",dest="output",help="This is the file name for the merged output file",metavar="output file")

(options, args) = parser.parse_args()

input_dir=options.input_dir
pattern=options.pattern
how=options.how
key=options.key
output=options.output

all_files = os.listdir(str(input_dir))
input_file=fn.filter(all_files,"*"+str(pattern))
bulk_df_list=[]
for i in range(len( input_file)):
	temp_df=pd.read_csv(str(input_dir)+"/"+input_file[i],sep="\t",header=0)
	bulk_df_list.append(temp_df)

bulk_df_list
merged_KO=ft.reduce(lambda left,right: pd.merge(left,right,how=how,on=key),bulk_df_list)
merged_KO.to_csv(str(output),sep="\t",index=False,na_rep="NA")
