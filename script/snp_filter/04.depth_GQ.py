import os,sys,re,gzip
depth_file=sys.argv[1]
snp_file=sys.argv[2]
GQ_num=float(sys.argv[3])
depth={}
with open(depth_file) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        line=line.split()
        depth[line[0]]=[float(line[1])/3,float(line[1])*3]
with gzip.open(snp_file,"rt") as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        if re.search("^##",line):
            print(line)
            continue
        if re.search("^#CHROM",line):
            print(line)
            smaple_dep={}
            line=line.split()
            for i in range(9,len(line)):
                smaple_dep[i]=depth[line[i]]
            continue
        line=line.split()
        for i in range(9,len(line)):
            if line[i][0:3] in ["./.",".|.",'.']:
                continue
            sample_snp=line[i].split(":")
            if sample_snp[3]=='.' :
                if float(sample_snp[2])>=smaple_dep[i][0] and float(sample_snp[2])<=smaple_dep[i][1]:
                    continue
                else:
                    sample_snp[0]='./.'
                    line[i]=':'.join(sample_snp)
            else :
                if float(sample_snp[2])>=smaple_dep[i][0] and float(sample_snp[2])<=smaple_dep[i][1] and float(sample_snp[3]) >=GQ_num:
                    continue
                else:
                    sample_snp[0]='./.'
                    line[i]=':'.join(sample_snp)
        print('\t'.join(line))
