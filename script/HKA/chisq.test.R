args <- commandArgs(TRUE)
data <- read.table(args[1],header = F)
person_TEST <- function(list_){
  m1 <- matrix(as.numeric(list_[c(2:5)]),ncol = 2)
  test <- chisq.test(m1,correct = T)
  p1 <- test$p.value
  return(p1)
}
data$p_value <- apply(data,1,person_TEST)
colnames(data) <- c("#gene","fixed_count","polymorphic_count","all_fixed","all_polymorphic","p_value")
write.table(data,file="select_gene_HKA.txt",row.names = F,col.names = T,quote = F,sep="\t")
