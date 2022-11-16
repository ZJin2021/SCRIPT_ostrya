import os,sys,re
import gzip
with gzip.open(sys.argv[1],'rt') as fp:
    for line in fp:
        line=line.strip()
        if re.search("^#",line):
            print(line)
            continue
        line=line.split()
        if len(line[4].split(','))>1:
            continue
        print('\t'.join(line))
