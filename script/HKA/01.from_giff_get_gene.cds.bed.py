import os,sys,re
import gzip
gff=sys.argv[1]
opener=gzip.open if gff.endswith('.gz') else open
#chr exonerate gene/mRNA/exon/CDS pos_start pos_end . +/- ./0/1/2
with opener(gff,'rt') as fp:
    gene_cds_dict={}
    for line in fp:
        line=line.strip()
        if not line:
            continue
        if re.search('^#',line):
            continue
        line=line.split()
        if line[2]=='CDS':
            gene_id=re.search('Parent=([^;]*)',line[-1]).group(1)
            cds_pos=[line[3],line[4]]
            gene_cds_dict.setdefault(gene_id,[line[0]]).append(cds_pos)
for gene_id in gene_cds_dict:
    contig=gene_cds_dict[gene_id][0]
    cds_pos=gene_cds_dict[gene_id][1:]
    cds_pos.sort(key=lambda x:int(x[0]))
    cds_pos2=['-'.join(x) for x in cds_pos]
    print(gene_id+'\t'+contig+'\t'+';'.join(cds_pos2))
