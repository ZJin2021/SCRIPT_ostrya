import sys,os,re
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)
def help():
    print("python3 xx.py merge_gene_cn sample_files type(int) > out_file")
try:
    merge_gene_cn=sys.argv[1]
    sample_sort_file=sys.argv[2]
    type_=int(sys.argv[3])
#GENE_ID              CHR     GENE_START      GENE_STOP       CN      CV_COUNT        PER_CV    LABEL
#OreG0000090.1   scaffold1021    1232461         1232733      6.50       2          9,4        Oja01
# 0                   1          2               3             4        5             6         7
except:
    help()
    sys.exit()
sort_samples=[]
with open(sample_sort_file) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        sort_samples.append(line)
gene_cns={}
with open(merge_gene_cn) as fp:
    for line in fp:
        line=line.strip()
        if not line or re.search("^GENE_ID",line):
            continue
        line=line.split()
        gene_cns.setdefault(line[7],{}).setdefault(line[0],line[4])
no_sample=set(sort_samples)-set(gene_cns.keys())
if no_sample:
    logging.warning('merge_gene_cn is not include %s,the all genes will be  %s'%(' '.join(list(no_sample)),type_))
all_genes=set()
for sample in gene_cns:
    all_genes=all_genes|set(gene_cns[sample].keys())
print("GENE_ID\t"+'\t'.join(sort_samples))
for gene in all_genes:
    print_line=[gene]
    for sample in sort_samples:
        if sample in gene_cns:
            if gene in gene_cns[sample]:
                print_line.append(gene_cns[sample][gene])
            else:
                print_line.append(str(type_))
        else:
            print_line.append(str(type_))
    print("\t".join(print_line))
logging.info("all done")
