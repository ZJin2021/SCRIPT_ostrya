import gzip,sys,re
vcf_file=sys.argv[1]
pop_file=sys.argv[2]
miss=float(sys.argv[3]) if len(sys.argv) > 3 else 0.2

opener=gzip.open if vcf_file.endswith('.gz') else open
out_pop=[]
groups={}
with open(pop_file) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        line=line.split()
        if line[1]=='Out_group':
            out_pop.append(line[0])
        else:
            groups.setdefault(line[1],[]).append(line[0])
no_list=['.','./.','.|.']

with opener(vcf_file,'rt') as fp:
    for line in fp:
        line=line.strip()
        if re.search('^#',line):
            print(line)
            if re.search('^#CHROM',line):
                line=line.split()
                sample_loci={}
                for i in range(9,len(line)):
                    sample_loci[line[i]]=i
            continue
        line_=line.split()
        if out_pop:
            N_out=0
            for sample in out_pop:
                if line_[sample_loci[sample]].split(':')[0] not in no_list:
                    N_out=1
                    break
            if N_out==0:
                continue
        N_group=0
        for pop in groups:
            pop_N=0
            for sample in groups[pop]:
                if line_[sample_loci[sample]].split(':')[0] in no_list:
                    pop_N+=1
            if pop_N/len(groups[pop]) > miss:
                N_group=1
                break
        if N_group==1:
            continue
        print(line)
