# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:19:51 2019

@author: liushang
"""
print('>===starting FPKM treating===<')
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')
print('***import the basic packages successfully***')
parser=argparse.ArgumentParser(description='packages for the fpkm treating in transcriptome')
parser.add_argument('-df','--dataframe',type=str,help='the path of dataframe file')
parser.add_argument('-rf','--result_file',type=str,help='the path of result file')
parser.add_argument('-op','--option',type=str,help='the option for expression files',
                    choices=['htseq','cufflink','cuffdiff'])
parser.add_argument('-p','--plot',type=str,choices=['vocano','density','length_plot'])
parser.add_argument('-g','--gff',type=str)
args=parser.parse_args()
result=args.result_file
gff_file=args.gff
if args.option=='htseq':
    dataframe=pd.read_csv(args.dataframe,sep='\t',header=None)
else:
    dataframe=pd.read_csv(args.dataframe,sep='\t')
print('the dataframe file has been loaded successfully')
if args.option=='htseq':
    def get_gff(input_gff,mid_gff,item):
        with open(input_gff,'r') as file:
            gff_list=[]
            for line in file:
                line=line.strip()
                if line[0]=='#':
                    continue
                else:
                    gff_list.append(line)
        with open(mid_gff,'w') as file:
            for i in gff_list:
                file.write(i+'\n')
        del gff_list
        table=pd.read_csv(mid_gff,sep='\t',header=None)
        exon_table=table[table[2]==item]
        gene_exon_list=[]
        for i in exon_table[8]:
            i=i.split(';')[0][3:]
            i=i.split('.')[0]
            gene_exon_list.append(i)
        exon_table[8]=gene_exon_list
        gene_exon_list_new=list(set(gene_exon_list))
        del gene_exon_list
        exon_table[9]=exon_table[4]-exon_table[3]+1
        gene_exon_length=[]
        append_number=0
        for i in gene_exon_list_new:
            i=i.strip()
            temp=exon_table[exon_table[8]==i]
            sum_length=0
            verbose_length=0
            start,end=[],[]
            for i,j in zip(temp[3],temp[4]):
                start.append(i)
                end.append(j)
                sum_length+=j-i+1
            if len(start)==1:
                gene_exon_length.append(sum_length)
            else:          
                for i,j in zip(start[1:],end[:-1]):
                    if j<i:
                        continue
                    else:
                        verbose_length+=j-i+1
                gene_exon_length.append(sum_length-verbose_length)
            append_number+=1
            if append_number%1000==0:
                print('the records has been extracted:%d'%append_number)
        del exon_table
        del table
        return gene_exon_list_new,gene_exon_length
    print('starting the gff_file loading, this will cost some time')
    gene_list,gene_exon_length=get_gff(gff_file,result+'/new_mid.gff3','exon')
    print('the gff file has been loaded successfully')
if args.option=='htseq':   
    def read_count_to_fpkm(read_count_df,target_fpkm_file):
        global gene_list,gene_exon_length
        pre_fpkm=[]
        for r_n,rc in zip(read_count_df[0],read_count_df[1]):
            for gn,gl in zip(gene_list,gene_exon_length):
                if r_n==gn:
                    pre_fpkm.append((r_n,rc,gl))
        rc_sum=0
        for i in pre_fpkm:
            rc_sum+=i[1]
        fpkm_dic={}
        for i in pre_fpkm:
            fpkm_dic[i[0]]=(i[1]/i[2])*1e9/rc_sum
        fpkm=pd.DataFrame.from_dict(fpkm_dic,orient='index')
        fpkm.to_csv(target_fpkm_file,header=None,sep='\t')
        del pre_fpkm,read_count_df,fpkm_dic,fpkm
    read_count_to_fpkm(dataframe,result+'/read_fpkm.csv')
    print('read counts has beeen translated into fpkm')
elif args.option=='cufflink': 
    def convert_cufflink_fpkm_to_tpkm(cufflink_df,tpkm_file):
        fpkm=cufflink_df['FPKM'].sum()
        cufflink_df['TPKM']=cufflink_df['FPKM']*1e6/fpkm
        cufflink_df.to_csv(tpkm_file,index=False)
        del cufflink_df
    convert_cufflink_fpkm_to_tpkm(dataframe,result+'/cufflink_tpkm.csv')
    print('fpkm has been transformed into TPKM')
elif args.option=='cuffdiff':
    def cuff_diff_fpkm_to_tpkm(fpkm_diff,tpkm_file):
        fpkm1=fpkm_diff['value_1'].sum()
        fpkm2=fpkm_diff['value_2'].sum()
        fpkm_diff['TPKM1']=fpkm_diff['value_1']*1e6/fpkm1
        fpkm_diff['TPKM2']=fpkm_diff['value_2']*1e6/fpkm2
        fpkm_diff.to_csv(tpkm_file,index=False)
        del fpkm_diff
    cuff_diff_fpkm_to_tpkm(dataframe,result+'/tpkm_diff.csv') 
    print('fpkm has been transformed into TPKM')
    def extract_target_gene(df,result):
        result_df=df[(df['q_value']<=0.05)&(df['log2(fold_change)']>=1)]
        result_df.to_csv(result+'/target_gene_diff.csv',index=False)
    extract_target_gene(dataframe,result)
    print('have extract genes with significant difference')
else:
    exit('the choices must be in htseq,cuffdiff and cufflink')
if (args.plot=='vocano')&(args.option=='cuffdiff'):    
    def vocano_plot(dataframe,result):
            import seaborn as sns
            import matplotlib.pyplot as plt
            dataframe['treated_q']=-np.log10(dataframe['q_value'])
            dataframe['class']='normal'
            dataframe.loc[(dataframe['log2(fold_change)']>1)&(dataframe['q_value']<0.05),'class']='up'
            dataframe.loc[(dataframe['log2(fold_change)']<-1)&(dataframe['q_value']<0.05),'class']='down'
            plt.figure(figsize=(15,10))
            ax=sns.scatterplot(x='log2(fold_change)',y='treated_q',
                               hue='class',
                               hue_order=('down','normal','up'),
                               palette=('blue','grey','red'),
                               data=dataframe)
            ax.set_ylabel('-log10(q_value)',fontweight='bold')
            ax.set_xlabel('foldchange',fontweight='bold')
            plt.legend(loc='best')
            plt.savefig(result+'/vocano.png')
            plt.savefig(result+'/vocano.svg')
            plt.clf()
    vocano_plot(dataframe,result)
    print('vocano plot has been performed successfully')
elif (args.plot=='length_plot'):
    def exon_length_plot(number_list,result):
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,15))
        plt.hist(number_list)
        plt.savefig(result+'/exon_length.svg')
        plt.savefig(result+'/exon_length.png')
        plt.clf()
    exon_length_plot(gene_exon_length,result)
    print('length distribution of genes\' exons has been ploted')
elif (args.plot=='density')&(args.option=='cuffdiff'):
    def density_for_cuffdiff(df,result):
        import matplotlib.pyplot as plt
        import numpy as np
        np.log2(df['value_1']+1).plot(kind='kde')
        np.log2(df['value_2']+1).plot(kind='kde')
        plt.savefig(result+'/density.svg')
        plt.savefig(result+'/density.png')
        plt.clf()
    density_for_cuffdiff(dataframe,result)
    print('density of pair fpkm values has been displayed')
elif (args.plot=='density')&(args.option=='cufflink'):
    def density_for_cufflink(df,result):
        import matplotlib.pyplot as plt
        import numpy as np
        np.log2(df['FPKM']+1).plot(kind='kde')
        plt.savefig(result+'/cufflink_density.svg')
        plt.savefig(result+'/cufflink_density.png')
        plt.clf()
    density_for_cufflink(dataframe,result)
    print('density of single fpkm values has been displayed')
else:
    exit('the plot should be consist with dataframe type')
print('>==the analysis has been finished===<')