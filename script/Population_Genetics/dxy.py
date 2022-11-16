import os,sys,re
import slide_window
import vcf_operating
import read_group
import read_dict
from multiprocessing import Pool
#https://github.com/mfumagalli/ngsPopGen/blob/9ee3a6d5e733c1e248e81bfc21514b0527da967b/scripts/getDxy.pl
def help():
    print("python3 xx.py vcf two_group_list fasta_fai outfile window step")
    print("two_group_list: classA\tHH1\nclassB\tHH2\nclassB\tHH3")
    sys.exit()
def maf_1_0(list_):
    count_0=sum([x.count('0') for x in list_])
    count_1=sum([x.count('1') for x in list_])
    all_count=count_0+count_1
    if all_count==0:
        return [0,0]
    maf_1=count_1/all_count
    maf_0=count_0/all_count
    return [maf_0,maf_1]

def each_maf(line):
    types1=[line[x][:3] for x in loci1]
    types2=[line[x][:3] for x in loci2]
    all_maf=maf_1_0(types1)+maf_1_0(types2)
    return all_maf

def each_dxy(line_):
    line=line_.split()
    chr_,pos=line[0],int(line[1])
    all_maf=each_maf(line)
    dxy=all_maf[0]*all_maf[3]+all_maf[1]*all_maf[2]
    return [chr_,pos,dxy]

def every_dxy(vcfs):
    head=next(vcfs).split()
    global loci1,loci2
    loci1=[head.index(x) for x in group_dict[group1]]
    loci2=[head.index(x) for x in group_dict[group2]]
    po1=Pool(theads)
    all_dxy=po1.map(each_dxy,vcfs)
    po1.close()
    dxy_dict={}
    for list_ in all_dxy:
        if list_[2]==0:
            continue
        dxy_dict.setdefault(list_[0],{}).setdefault(list_[1],list_[2])
    return dxy_dict

def window_dxy(cal_window):
    chr_=cal_window[0]
    if chr_ not in dxy_snp_pos_dict:
        return
    chr_snp_pos=dxy_snp_pos_dict[chr_]
    N_SNP,ALL_DXY=0,0
    for pos in chr_snp_pos:
        if pos<cal_window[1]:
            continue
        elif cal_window[1]<=pos<=cal_window[2]:
            N_SNP+=1
            ALL_DXY+=dxy_dict[chr_][pos]
        elif pos >cal_window[2]:
            break
    if ALL_DXY==0 or N_SNP==0:
        return
    ALL_DXY=ALL_DXY/(cal_window[2]-cal_window[1]+1)
    return [cal_window[0],str(cal_window[1]),str(cal_window[2]),str(N_SNP),str(ALL_DXY)]

def write_result(all_window_dxy,out_file):
    all_window_dxy.sort(key=lambda x:(x[0],int(x[1])))
    with open(out_file,'w') as fp:
        print('CHROM\tSTART\tEND\tN_SNP\tDXY',file=fp)
        for win in all_window_dxy:
            print('\t'.join(win),file=fp)

def main(vcf_file,group_file,fasta_fai,out_file,window,step):
    global group_dict,group1,group2
    group_dict=read_group.Read_group1(group_file)
    [group1,group2]=list(group_dict.keys())
    vcfs=For_vcf.vcf_read2(vcf_file)
    global dxy_dict,dxy_snp_pos_dict
    dxy_dict=every_dxy(vcfs)
    dxy_snp_pos_dict={}
    for chr_ in dxy_dict:
        dxy_snp_pos_dict[chr_]=sorted(dxy_dict[chr_].keys())
    chr_len_dict=read_dict.read_fai(fasta_fai)
    all_windows=[]
    #[[chr,start,end]]
    for chr_ in chr_len_dict:
        all_windows+=slide_window.slide_window2(chr_,window,step,1,chr_len_dict[chr_])
    po2=Pool(theads)
    all_window_dxy=po2.map(window_dxy,all_windows)
    po2.close()
    #all_window_dxy=list(filter(None,all_window_dxy))
    all_window_dxy=[x for x in all_window_dxy if x !=None]
    write_result(all_window_dxy,out_file)


if __name__=="__main__":
    try:
        [vcf_file,group_file,fasta_fai,out_file]=sys.argv[1:5]
        window=int(sys.argv[5])
        step=int(sys.argv[6])
    except:
        help()
    theads=30
    For_vcf=vcf_operating.for_vcf()
    main(vcf_file,group_file,fasta_fai,out_file,window,step)
