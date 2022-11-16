#<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>
#this is for group_reorhanization_rate
#use : R3.6.1 FastEPRR 1.0 beagle 4.1
#conda avtivate R3
#<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>
import os,sys,re
import gzip
from pathlib import Path
from multiprocessing import Pool
def help():
    print('use:python xx.py O=utdir V=vcf_file R=fa.fai W=window-size (K) S=step-size (K)  C=control-snp  T=process_number Phase=T[F]')
    print('the process_number at most 40 and shold not much than the number of your chmomosome')
    print('control-snp :winDXThreshold , Minimum threshold of sum of folded doubletons and folded xtons for each window. The window will be excluded if the sum of folded doubletons and folded xtons is smaller than winDXThreshold. Default is 10 and the minimum is 2. \n. It is the minimum number of ξ2′+ ξí′ for one window. The windows are excluded if the number of ξ2′ + ξí′ of these windows is less than winSNPThreshold. If missing, the default value is 30')

out_dir='.'
vcf_file,fa_fai='',''
window_size,step_size=10,10
control_snp=10
num=10
Phase=True
for opt in sys.argv[1:]:
    if opt.startswith('O='):
        out_dir=opt.split('=')[1]
    elif opt.startswith('V='):
        vcf_file=opt.split('=')[1]
    elif opt.startswith('R='):
        fa_fai=opt.split('=')[1]
    elif opt.startswith('W='):
        window_size=int(opt.split('=')[1])
    elif opt.startswith('S='):
        step_size=int(opt.split('=')[1])
    elif opt.startswith('C='):
        control_snp=int(opt.split('=')[1])
    elif opt.startswith('T='):
        num=int(opt.split('=')[1])
    elif opt=="Phase=F":
        Phase=False

if not vcf_file or not fa_fai:
    help()
    sys.exit()

chr_len={}
with open(fa_fai) as fp:
    for line in fp:
        line=line.strip()
        if not line:
            continue
        line=line.split()
        chr_len[line[0]]=int(line[1])

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
vcf_name = re.search('([^/]*$)',vcf_file).group(1)
vcf_path = os.path.abspath(vcf_file)
os.chdir(out_dir)
dir_dict = ['01.vcf','02.split','step1','step2','step3']
for dir_n in dir_dict:
    if not Path(dir_n).is_dir():
        os.mkdir(dir_n)
if not Phase:
    os.system('ln -s %s 01.vcf/gt.impute.vcf.gz'%(vcf_path))
else:
    os.system('ln -s %s 01.vcf'%(vcf_path))
    if not Path('01.vcf/gt.impute.vcf.gz').is_file():
        print('now Phase')
        os.system('java -jar -Xmx20G /data/00/software/beagle/beagle.18May20.d20.jar gt=01.vcf/%s out=01.vcf/impute nthreads=20'%(vcf_name))
        print('phase down')

new_vcf_name = 'gt.impute.vcf.gz'

if not Path('02.split/IO').is_file():
    os.chdir('01.vcf')
    print('start split')
    vcf_head = []
    vcf_dict = {}
    if not Phase:
        with gzip.open(new_vcf_name,'rt') as fp:
            for line in fp:
                line=line.strip()
                if not line:
                    continue
                fis_str = re.search('^#',line)
                if fis_str:
                    vcf_head.append(line)
                    continue
                line=line.split('\t')
                line[5]='.'
                info=re.search('(AC=\d+)[\w\W]*(AN=\d+)',line[7])
                line[7]=info.group(1)+';'+info.group(2)
                line[8]='GT'
                for i in range(9,len(line)):
                    line[i]=line[i].split(':')[0]
                    if line[i] in ['.','./.']:
                         line[i]='.|.'
                    else:
                        line[i]=re.sub('/','|',line[i])
                vcf_dict.setdefault(line[0],[]).append('\t'.join(line))
    else:
        with gzip.open(new_vcf_name,'rt') as fp:
            for line in fp:
                line=line.strip()
                if not line:
                    continue
                if re.search('^#',line):
                    vcf_head.append(line)
                    continue
                line_=line.split('\t')
                vcf_dict.setdefault(line_[0],[]).append(line)
    vcf_keys = sorted(vcf_dict.keys())
    for key_ in vcf_keys:
        if key_ not in chr_len:
            raise ValueError('Some chr in vcf but not in fa.fai, please check the two file')
    os.chdir('..')
os.chdir('02.split')
def split_vcf(vcf_lines,vcf_key):
    global vcf_head
    vcf_w = open('%s.vcf'%(vcf_key),'w')
    for i in range(len(vcf_head)):
        print(vcf_head[i],file = vcf_w)
    while vcf_lines:
        vcf_line = vcf_lines.pop(0)
        print(vcf_line,file = vcf_w)
    vcf_w.close()
    os.system('bgzip %s.vcf'%(vcf_key))

if not Path('IO').is_file():
    po1 = Pool(num)
    for vcf_key in vcf_keys:
        vcf_lines = vcf_dict.pop(vcf_key)
        po1.apply_async(split_vcf,(vcf_lines,vcf_key,))
    po1.close()
    po1.join()
    os.system('touch IO')
os.chdir('..')
print('split done')

vcf_keys=[x.split('.')[0] for x in os.listdir("02.split") if x.endswith('.vcf.gz')]
job_num = len(vcf_keys)

def step1_sys(contig):
    vcfFilePath = '02.split/%s.vcf.gz'%(contig)
    step1 = 'library(FastEPRR); FastEPRR_VCF_step1(vcfFilePath="%s", winLength="%d", stepLength="%d", winDXThreshold=%d, erStart = "0.001", erEnd = "%s", srcOutputFilePath="step1/%s.txt")'%(vcfFilePath,window_size,step_size,control_snp,chr_len[contig]/1000,contig)
    #print(step1)
    #qualThreshold,default=20
    os.system("echo '%s' | Rscript -"%(step1))

def step2_sys(i):
    global job_num
    step2 = 'library(FastEPRR); FastEPRR_VCF_step2(srcFolderPath="step1", jobNumber=%d, currJob=%d, DXOutputFolderPath="step2" )'%(job_num,i+1)
    os.system("echo '%s' | Rscript -"%(step2))

print("start step1")
po2 = Pool(num)
for vcf_key in vcf_keys:
    po2.apply_async(step1_sys,(vcf_key,))
#po2.map(step1_sys,vcf_keys)
po2.close()
po2.join()

print("step1 done")

print("start step2")
po3 = Pool(num)
for i in range(job_num):
    po3.apply_async(step2_sys,(i,))
po3.close()
po3.join()
print("step2 done")

print("start step3")
step3 = 'library(FastEPRR); FastEPRR_VCF_step3(srcFolderPath="step1",DXFolderPath= "step2", finalOutputFolderPath="step3")'
os.system("echo '%s' | Rscript -"%(step3))
print("step3 done")
print("all done")

