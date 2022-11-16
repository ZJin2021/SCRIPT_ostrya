####
## Script to preform permutation test on the specified matrices of gene CN counts
####
setwd("C:\\Users\\16640\\Desktop\\R\\oju\\cnv\\SAMPLE")
library(devtools)
library(Hmisc)

#columns: BK_SAMN01057691	BK_SRR7813602	BB_SRR7758718	BB_SAMN02256313	BB_SAMN02256315	BB_SAMN02256316	BB_SAMN02256317	BB_SAMN02256318	BB_SAMN02256319	BB_SAMN02256320	BB_SAMN02256321	BB_SAMN02256322	PB_SAMN02261811	PB_SAMN02261819	PB_SAMN02261821	PB_SAMN02261826	PB_SAMN02261840	PB_SAMN02261845	PB_SAMN02261851	PB_SAMN02261853	PB_SAMN02261854	PB_SAMN02261856	PB_SAMN02261858	PB_SAMN02261865	PB_SAMN02261868	PB_SAMN02261870	PB_SAMN02261871	PB_SAMN02261878	PB_SAMN02261880

cnvFC = read.csv(file="gene_cn_matrix.txt", sep = '\t', row.names = 1, header=TRUE)

#define sample lists

samp=c("Oja01","Oja02","Oja03","Oja04","Oja05","Oja06","Oja07","Oja08","Oja09","Oja10","Oja11","Oja12","Oja13","Oja14","Oja15","Oja16","Oja17","Oja18","Oja19","Oja20","Oja21","Oja23","Oja24","Oja25","Oja26","Oja27","Omu01","Omu02","Omu03","Omu04","Omu05","Omu06","Omu07","Omu08","Omu09","Omu10","Omu11","Omu12","Omu13","Omu14")
spec=c("Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Oja","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu","Omu")

#spec=c("Black Bear","Black Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Brown Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear","Polar Bear")


##General format for Vst calculation for each row
####
##Vp*n==[var(BB)*pBB+var(PB)*pPB]
##[var(BBPB) - (Vp*n/pPBBB)]/var(BBPB)
####
##[var(BBPB) - ([var(BB)*pBB+var(PB)*pPB]/pPBBB)]/var(BBPB)
###
##[var(BBPB) - ([var(BB)*9+var(PB)*17]/26)]/var(BBPB)
####

#calculate per-gene CN variance in BB, PB and BBPB
V1all = t(apply(cnvFC, 1, function(y) c(var(y[1:26])*25/26, var(y[27:40])*13/14, var(y[1:40])*39/40))) ##equal var.p in Excel

#calculate per-gene Vst
VstFCall = (V1all[,3] - ((V1all[,1]*26 + V1all[,2]*14)/40))/V1all[,3]

###
#Remove all genes with variance == 0 to create a reduced matric containing on CN variable genes
###

NonzeroVarFC = names(which(V1all[,3]>0))

cnvFCreduced = cnvFC[(rownames(cnvFC) %in% NonzeroVarFC),]
####
#calculate Vst (this is redundant from above but it only takes a few seconds)

V1 = t(apply(cnvFCreduced, 1, function(y) c(var(y[1:26])*25/26, var(y[27:40])*13/14, var(y[1:40])*39/40)))

VstFCreduced = (V1[,3] - ((V1[,1]*26 + V1[,2]*14)/40))/V1[,3]

#shuffle reduced CN matrix
cnvFCrdm = cnvFCreduced[, c(sample(1:ncol(cnvFCreduced)))]

#calculate Vst on shuffled matrix

V1rdm = t(apply(cnvFCrdm , 1, function(y) c(var(y[1:26])*25/26, var(y[27:40])*13/14, var(y[1:40])*39/40)))

VstFCrdm = (V1rdm[,3] - ((V1rdm[,1]*26 + V1rdm[,2]*14)/40))/V1rdm[,3]

VstPerFC = VstFCrdm

#repeat x1000
for (i in 1:9999){

cnvFCrdm = cnvFCreduced[, c(sample(1:ncol(cnvFCreduced)))]

V1rdm = t(apply(cnvFCrdm , 1, function(y) c(var(y[1:26])*25/26, var(y[27:40])*13/14, var(y[1:40])*39/40)))

VstFCrdm = (V1rdm[,3] - ((V1rdm[,1]*26 + V1rdm[,2]*14)/40))/V1rdm[,3]

#add vector to permutation test matrix
VstPerFC = cbind(VstFCrdm, VstPerFC)
}

############
##Extract cutoffs for each gene from permutations
############
VstFCcutoff.50 = apply(VstPerFC, 1, function(x) quantile(x, probs=.50))
VstFCcutoff.95 = apply(VstPerFC, 1, function(x) quantile(x, probs=.95))
VstFCcutoff.99 = apply(VstPerFC, 1, function(x) quantile(x, probs=.99))

VstCutoffsFC = cbind(VstFCcutoff.50, VstFCcutoff.95, VstFCcutoff.99, VstFCreduced)
write.table(VstCutoffsFC,file="gene_vcf.count",sep = "\t")
