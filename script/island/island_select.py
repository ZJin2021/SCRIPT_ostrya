import os,sys,re
from decimal import *
import numpy as np
import random
import scipy
import scipy.stats
import math
from multiprocessing import Pool

def help():
    print("python3 xx.py window_fst snp_fst  window step [--top=? top_num default=0.01] [--min=?  min_snp_num default=10] [--per=? permute_num default=1000000] [--T=? theads default=1]")
    print("window_fst and snp_fst calculated by vcftools")
    sys.exit()
def read_window_fst(window_fst_file,window,min_snp_num):
    #CHROM   BIN_START       BIN_END N_VARIANTS      WEIGHTED_FST    MEAN_FST
    windows_fst=[]
    #[chr,start,SNP_num,Mean_fst]
    with open(window_fst_file) as fp:
        for line in fp:
            line=line.strip()
            if not line or re.search("^CHROM",line):
                continue
            line=line.split()
            if int(line[2])-int(line[1])+1 !=window:
                continue
            if int(line[3])<min_snp_num:
                continue
            if float(line[-1])<=0:
                continue
            windows_fst.append([line[0],line[1],int(line[3]),float(line[-1])])
    return windows_fst
def read_snp_fst(snp_fst_file):
    snp_fst=[]
    with open(snp_fst_file):
        with open(snp_fst_file) as fp:
            for line in fp:
                line=line.strip()
                if not line or re.search("^CHROM",line):
                    continue
                line=line.split()
                if line[2]=="-nan":
                    continue
                snp_fst.append(float(line[2]))
    return snp_fst
def top_window_fst(windows_fst,top_num):
    windows_fst.sort(key=lambda x:-x[3])
    Tnum=int(len(windows_fst)*top_num//1)
    m=0
    for i in range(len(windows_fst)):
        if m<Tnum:
            windows_fst[i].append("TOP")
        else:
            windows_fst[i].append("BACKGROUND")
        m+=1
    #[chr,start,SNP_num,Mean_fst,type]
    return windows_fst
def genomic_sample_snp(num_):
    back_snp_fst=random.sample(snp_fst, num_)
    mean_=np.mean(back_snp_fst)
    return mean_
def print_result(windows_fst):
    with open("permute_window_snp_fst.result","w") as fp:
        #[chr,start,SNP_num,Mean_fst,type,all_mean,all_std,R_std,kstest_p,normaltest_p,P_value]
        print("CHR\tSTART\tSNP_num\tFST\ttype\tMEAN_FST\tstd\tR_std\tkstest_p\tnormaltest_p\tP_value",file=fp)
        for window_ in windows_fst:
            print("\t".join([str(x) for x in window_]),file=fp)
def main(window_fst_file,snp_fst_file,window,step,top_num,min_snp_num,permute_num,theads):
    windows_fst=read_window_fst(window_fst_file,window,min_snp_num)
    global snp_fst
    snp_fst=read_snp_fst(snp_fst_file)
    windows_fst=top_window_fst(windows_fst,top_num)
    snp_num_window_fst={}
    for i in range(len(windows_fst)):
        window_=windows_fst[i]
        window_snp_num=window_[2]
        if window_snp_num in snp_num_window_fst:
            windows_fst[i]=windows_fst[i]+snp_num_window_fst[window_snp_num]
        else:
            po1=Pool(theads)
            result_=po1.map(genomic_sample_snp,[window_snp_num]*permute_num)
            po1.close()
            all_mean=np.mean(result_)
            all_var=np.var(result_)
            all_std=math.sqrt(all_var)
            R_std=math.sqrt(all_var*permute_num/(permute_num-1))
            P_value=len([x for x in result_ if x >=window_[3]])/permute_num
            kstest_p=scipy.stats.kstest(result_,'norm')[1]
            normaltest_p=scipy.stats.normaltest(result_, axis=None)[1]
            windows_fst[i]=windows_fst[i]+[all_mean,all_std,R_std,kstest_p,normaltest_p,P_value]
            snp_num_window_fst[window_snp_num]=[all_mean,all_std,R_std,kstest_p,normaltest_p,P_value]
    print_result(windows_fst)
if __name__=="__main__":
    try:
        window_fst_file=sys.argv[1]
        snp_fst_file=sys.argv[2]
        window=int(sys.argv[3])
        step=int(sys.argv[4])
    except:
        help()
    top_num=0.01
    min_snp_num=10
    permute_num=1000000
    theads=1
    if len(sys.argv)>5:
        options=sys.argv[5:]
        for option in options:
            if "--top" in option:
                top_num=float(re.search("--top=([^ ]*)$",option).group(1))
            elif "--min" in option:
                min_snp_num=int(re.search("--min=([^ ]*)$",option).group(1))
            elif "--per" in option:
                permute_num=int(re.search("--per=([^ ]*)$",option).group(1))
            elif "--T" in option:
                theads=int(re.search("--T=([^ ]*)$",option).group(1))
    main(window_fst_file,snp_fst_file,window,step,top_num,min_snp_num,permute_num,theads)
