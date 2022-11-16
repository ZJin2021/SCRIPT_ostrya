#!/bin/python3
#***************
#zhangjin lzu
#***************
import os,sys,re
import gzip
import copy
from scipy import stats
from multiprocessing import Pool
try:
    fst_file=sys.argv[1]
    #CHROM   POS     WEIR_AND_COCKERHAM_FST
    bed_file=sys.argv[2]
    vcf_file=sys.argv[3]
    pop_file=sys.argv[4]#only for selet pop
    theads=int(sys.argv[5])
except:
    print("python3 xx.py fst_file(vcftools-single) bed_file vcf_file pop_file(select_pop_samples) theads ")
    sys.exit()
try:
    sys.argv[6]=='less' 
    run_t='more' 
	
except:
    run_t='less' 
    #the genome-wide average A/B, which was computed as the sum of A and B values across all genes analyzed.

global snp_fsts
snp_fsts={}
with open(fst_file) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        if re.search('^CHROM',line):
            continue
        line=line.split()
        if line[2]=="-nan":
            continue

        
        if float(float(line[2]))>=0.95:
            snp_fsts.setdefault(line[0],{}).setdefault(int(line[1]),float(line[2]))
samples=[]

with open(pop_file) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        samples.append(line)
        
opener=gzip.open if vcf_file.endswith(".gz") else open
def read_vcf(file_):
    with opener(file_,'rt') as fp:
        for line in fp:
            line=line.strip()
            if re.search("^##",line):
                continue
            yield line
def sample_loci(line_,samples):
    line=line_.split()
    global all_loci
    all_loci=[line.index(x) for x in samples]

def find_polymorphic(line_):
    line=line_.split()
    all_ref=0
    all_alt=0
    for i in all_loci:
        type_=line[i][0:3]
        if '0' in type_ and '1' not in type_:
            all_ref+=2
        elif '1' in type_ and '0' not in type_:
            all_alt+=2
        elif '1' in type_ and '0' in type_:
             all_ref+=1
             all_alt+=1
    if (all_ref+all_alt)==0:
        return
    rato=all_ref/(all_ref+all_alt)
    if rato >=0.1 and rato <=0.9:
        snp_fsts.setdefault(line[0],{}).setdefault(int(line[1]),0.5)

def dict2list(dict_):
    result_dict={}
    for chr_ in dict_:
        result_dict[chr_]=[]
        for start in dict_[chr_]:
            result_dict[chr_].append([start,dict_[chr_][start]])
    return result_dict

lines=read_vcf(vcf_file)
sample_info=next(lines)
sample_loci(sample_info,samples)
for line in lines:
    find_polymorphic(line)
snp_fsts=dict2list(snp_fsts)


gene_pos={}
with open(bed_file) as fp:
    #gene_name\tchr\tpostart-end;posstat-end
    for line in fp:
        line=line.strip()
        if not line:
            continue
        line=line.split()
        pos=[x.split("-") for x in  line[2].split(';')]
        for i in range(len(pos)):
            for j in range(len(pos[i])):
                pos[i][j]=int(pos[i][j])
        #{chr:{gene:[[1,2][5,6]]}}
        gene_pos.setdefault(line[1],{}).setdefault(line[0],pos)

def fixed_loci(contig):
    congtig_dict={}
    contig_snp_fsts=snp_fsts[contig]
    contig_snp_fsts.sort(key=lambda x : x[0])
    for gene in gene_pos[contig]:
        pos=copy.deepcopy(gene_pos[contig][gene])
        pos.sort(key=lambda x: x[0])
        gene_fst=[]
        #[[snp_pos,fst]]
        for snp in contig_snp_fsts:
            if not pos:
                break
            m=0
            while m==0:
                if not pos:
                    break
                if snp[0] < pos[0][0]:
                    m=1
                elif snp[0] >= pos[0][0] and snp[0] <= pos[0][1]:
                    snp=snp+[contig]
                    gene_fst.append(snp)
                    m=1
                elif snp[0] > pos[0][1]:
                    pos.pop(0)
        if len(gene_fst)==0:
            continue
        congtig_dict[gene]=gene_fst
    return congtig_dict

congtigs=list(set(snp_fsts.keys())&set(gene_pos.keys()))
gene_snp_fst_list=[]
po1=Pool(theads)
gene_snp_fst_list+=po1.map(fixed_loci,congtigs)
po1.close()
gene_snp_fst_dict={}
#{gene:[[snp_pos,fst]]}
for i in gene_snp_fst_list:
    gene_snp_fst_dict.update(i)

all_loci=0
all_fixed=0
file_fix=open('gene_fixed_loci.txt','w')
print('#gene\tchr\tpos\tfst',file=file_fix)
gene_fixed_dict={}
#{gene:[gene_all_fixed,gene_all_loci]}
for gene in gene_snp_fst_dict:
    gene_all_loci=len(gene_snp_fst_dict[gene])
    gene_all_fixed=0
    all_loci+=len(gene_snp_fst_dict[gene])
    for snp in gene_snp_fst_dict[gene]:
        if snp[1]>=0.95:
            gene_all_fixed+=1
            all_fixed+=1
            print(gene+'\t'+snp[2]+'\t'+str(snp[0])+'\t'+str(snp[1]),file=file_fix)
    gene_fixed_dict[gene]=[gene_all_fixed,gene_all_loci]
file_fix.close()

if run_t=='less':
    all_nofix=str(all_loci-all_fixed)
    all_fixed=str(all_fixed)
elif run_t=='more':
    all_gene_num=len(gene_fixed_dict)
    all_nofix=str(int(((all_loci-all_fixed)/all_gene_num)//1))
    all_fixed=str(int((all_fixed/all_gene_num)//1))

with open("gene_fixed_and_polymorphic","w") as fp:
    print("#gene\tfixed_count\tpolymorphic_count\tall_fixed\tall_polymorphic",file=fp)
    for gene in gene_fixed_dict:
        loci_info=gene_fixed_dict[gene]
        gene_fixed=str(loci_info[0])
        gene_nofix=str(loci_info[1]-loci_info[0])
        print("\t".join([gene,gene_fixed,gene_nofix,all_fixed,all_nofix]),file=fp)

os.system("Rscript /data/00/user/user153/script/07.differentiation/HKA/chisq.test.R gene_fixed_and_polymorphic")
