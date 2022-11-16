import os,sys,re
from multiprocessing import Pool
import copy
import numpy
def help():
    print("python3 xx.py gff cnv_file_from_FREEC label> output_file")
#       "cnv_file" is tab delimited output file of control-freec ending in "CNVs"
#       "cnv_file" consists of 5 columns: scaffold/chromosome   start   stop    copy_number     gain or loss
#
#       Example "cnv_file" (no header line necessary):
#       scaffold79      60000   77999   4       gain
#       scaffold79      492000  509999  5       gain
try:
    gff_file=sys.argv[1]
    cnv_raw=sys.argv[2]
    sample_name=sys.argv[3]
except:
    help()
    sys.exit()
gene_dict={}
with open(gff_file) as fp:
    for line in fp:
        line=line.strip()
        if not line or re.search("^#",line):
            continue
        line=line.split()
        if line[2] != "mRNA":
            continue
        gene=re.search("ID=([^;]*)",line[8]).group(1)
        gene_dict[gene]=[line[0],int(line[3]),int(line[4])]
cnv_dict={}
with open(cnv_raw) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        line=line.split()
        cnv_dict.setdefault(line[0],[]).append([int(line[1])+1,int(line[2])+1,int(line[3])])
for chr_ in cnv_dict:
    cnv_dict[chr_].sort(key=lambda x:x[0])
def gene_cn_count(gene_pos,cnvs):
    CV_list=[]
    gene_len=gene_pos[1]-gene_pos[0]+1
    while cnvs:
        if gene_pos[1] < cnvs[0][0]:
            break
        elif gene_pos[0] > cnvs[0][1]:
            cnvs.pop(0)
        elif gene_pos[0] <= cnvs[0][1] and gene_pos[1] >= cnvs[0][0]:
            cnv_len=cnvs[0][1]-cnvs[0][0]+1
            if gene_pos[1] <= cnvs[0][1] and gene_pos[0] >= cnvs[0][0]:
                CV_list.append(cnvs.pop(0)[2])
            elif gene_pos[1] >= cnvs[0][1] and gene_pos[0] <= cnvs[0][0]:
                CV_list.append(cnvs.pop(0)[2])
            elif gene_pos[1] >= cnvs[0][1] and gene_pos[0] >= cnvs[0][0]:
                over_len=cnvs[0][1]-gene_pos[0]+1
                if over_len/gene_len >= 0.5 or over_len/cnv_len >0.5:
                    CV_list.append(cnvs.pop(0)[2])
                else:
                    cnvs.pop(0)
            elif gene_pos[1] <= cnvs[0][1] and gene_pos[0] <= cnvs[0][0]:
                over_len=gene_pos[1]-cnvs[0][0]+1
                if over_len/gene_len >= 0.5 or over_len/cnv_len >0.5:
                    CV_list.append(cnvs.pop(0)[2])
                else:
                    cnvs.pop(0)
    return CV_list
def gene_CN(gene):
    gene_cn_dict={}
    gene_pos=gene_dict[gene][1:]
    chr_=gene_dict[gene][0]
    if chr_ not in cnv_dict:
        return gene_cn_dict
    cnvs=copy.deepcopy(cnv_dict[chr_])
    CV_list=gene_cn_count(gene_pos,cnvs)
    if not CV_list:
        return gene_cn_dict
    CN=numpy.mean(CV_list)
    if CN%1==0:
        CN=str(int(CN))
    else:
        CN="%.2f"%(CN)
    CV_COUNT=str(len(CV_list))
    PER_CV=','.join([str(x) for x in CV_list])
    gene_cn_dict[gene]=[gene,chr_,str(gene_pos[0]),str(gene_pos[1]),CN,CV_COUNT,PER_CV]
    return gene_cn_dict
po1=Pool(50)
result=po1.map(gene_CN,sorted(gene_dict.keys()))
po1.close()
head=['GENE_ID','CHR','GENE_START','GENE_STOP','CN','CV_COUNT','PER_CV','LABEL']
print('\t'.join(head))
for i in result:
    if not i:
        continue
    for gene in i:
        print('\t'.join(i[gene])+'\t'+sample_name)
