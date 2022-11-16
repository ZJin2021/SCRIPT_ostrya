#取消科学计数法
options(scipen = 200)
args <- commandArgs(TRUE)
data1=read.table(args[1],header = T)
#CHR START SNP_num FST type MEAN_FST std R_std kstest_p normaltest_p P_value
data2=data1[data1$type=="TOP",]
p_fdr <- function(data_){
    std_norm_P=1-pnorm(data_$FST,mean=data_$MEAN_FST,sd=data_$std)
    R_std_norm_P=1-pnorm(data_$FST,mean=data_$MEAN_FST,sd=data_$R_std)
    FDR=p.adjust(data_$P_value,method = "fdr")
    std_norm_FDR=p.adjust(std_norm_P,method = "fdr")
    R_std_norm_FDR=p.adjust(R_std_norm_P,method = "fdr")
    end=data_$START+9999
    data_<-cbind.data.frame(data_[,1:2],END=end,data_[,3:length(data_[1,])])
    data_<-cbind.data.frame(data_,std_norm_P=std_norm_P,R_std_norm_P=R_std_norm_P,FDR=FDR,std_norm_FDR=std_norm_FDR,R_std_norm_FDR=R_std_norm_FDR)
    return(data_)
}
data1<-p_fdr(data1)
data2<-p_fdr(data2)
write.table(data1,file="whole_genomic_FDR.txt",quote = FALSE,row.names = FALSE,
              col.names = T,sep="\t")
write.table(data2,file="TOP_genomic_FDR.txt",quote = FALSE,row.names = FALSE,
	    col.names = T,sep="\t")
