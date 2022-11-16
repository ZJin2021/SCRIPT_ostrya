import os,sys,re
from decimal import *

def help():
    print("python3 xx.py genomic_FDR.txt out_file [--fdr=? [FDR/std_norm_FDR/R_std_norm_FDR default=R_std_norm_FDR]] [--p=? p_value default=0.01")
    sys.exit()

try:
    genomic_FDR_txt=sys.argv[1]
    out_file=sys.argv[2]
except:
    help()

fdr=-1
p_value=0.01
type_loci={'FDR':-3,'std_norm_FDR':-2,'R_std_norm_FDR':-1}
if len(sys.argv) >3:
    for option in sys.argv[3:]:
        if '--fdr' in option:
            fdr_type=re.search("^--fdr=(.*)$",option).group(1)
            fdr=type_loci[fdr_type]
        elif '--p' in option:
            p_value=float(re.search("^--p=([\w\W]*)$",option).group(1))
island_window={}

with open(genomic_FDR_txt) as fp:
    for line in fp:
        line=line.strip()
        if not line or re.search('^CHR',line):
            continue
        line=line.split()
        if line[5] !='TOP':
            continue
        if Decimal(line[fdr])>p_value:
            continue
        island_window.setdefault(line[0],[]).append([int(line[1]),int(line[2])])

def merge_window(windows):
    merge_windows=[]
    windows.sort(key=lambda x:x[0])
    last_start=windows[0][0]
    last_end=windows[0][1]
    if len(windows)>1:
        for window in windows[1:]:
            if window[0]<=last_end+1:
                last_end=window[1]
            else:
                window_len='%dK'%(Decimal(last_end-last_start+1)/1000//1)
                merge_windows.append([str(last_start),str(last_end),window_len])
                last_start=window[0]
                last_end=window[1]
    window_len='%dK'%(Decimal(last_end-last_start+1)/1000//1)
    merge_windows.append([str(last_start),str(last_end),window_len])
    return merge_windows

island_merge_window={}
for chr_ in island_window:
    island_merge_window[chr_]=merge_window(island_window[chr_])
with open(out_file,"w") as fp:
    print("CHROM\tSTART\tEND\tLENGTH",file=fp)
    for chr_ in island_merge_window:
        for window in island_merge_window[chr_]:
            print(chr_+'\t'+'\t'.join(window),file=fp)
