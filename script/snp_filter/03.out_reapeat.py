import os,sys,re
import gzip
sv_intraval={}
with open(sys.argv[1]) as fp:
    for line in fp:
        line=line.strip()
        if not line or re.search('^#',line):
            continue
        line=line.split()
        start=int(line[1])
        end=int(line[2])
        sv_intraval.setdefault(line[0],[]).append([start,end])

for chr_ in sv_intraval:
    sv_intraval[chr_].sort(key=lambda x:x[0])

with gzip.open(sys.argv[2],'rt') as fp:
    last_chr=0
    last_loci=0
    for line_ in fp:
        line_=line_.strip()
        if not line_ or re.search('^#',line_):
            print(line_)
            continue
        line=line_.split()
        if line[0]!=last_chr:
            last_chr=line[0]
            last_loci=0
        if last_chr not in sv_intraval:
            print(line_)
            continue
        while last_loci<len(sv_intraval[last_chr]):
            if int(line[1]) < sv_intraval[last_chr][last_loci][0]:
                print(line_)
                break
            elif int(line[1])>=sv_intraval[last_chr][last_loci][0] and int(line[1])<= sv_intraval[last_chr][last_loci][1]:
                break
            elif int(line[1])>sv_intraval[last_chr][last_loci][1]:
                last_loci+=1
                continue
        if last_loci>=len(sv_intraval[last_chr]):
            print(line_)
