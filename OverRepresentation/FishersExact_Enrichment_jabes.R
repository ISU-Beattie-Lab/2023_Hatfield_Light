# Purpose: Identify the functional categories that are over-represented among
# BphP1-dependent photoresponsive genes in Pss B728a.

#=============================================================================

counts<-read.csv("./BphP1DepDEGs_function_counts.csv",header=T)

counts<-counts[!is.na(counts[,1]),]
counts<-counts[,!is.na(counts[1,])]
colnames(counts)[2]<-"Func.Cat"

head(counts)


pval<-NULL
for(row in 2:nrow(counts)){
  p.val<-NULL
  for(col in 3:ncol(counts)){
    
    out.DE<-counts[1,col]-counts[row,col]; out.EE<-counts[1,1]-counts[row,1]-out.DE
    in.DE<-counts[row,col]; in.EE<-counts[row,1]-in.DE
    tab<-matrix(c(out.DE,out.EE,in.DE,in.EE),2,2)
    colnames(tab)<-c("Out","In"); rownames(tab)<-c("DE","EE")
    dir<-sign(in.DE/in.EE-out.DE/out.EE)
    p.val<-c(p.val,dir*fisher.test(tab)$p.value)
  }
  pval<-rbind(pval,p.val)
}

colnames(pval)<-colnames(counts)[3:ncol(counts)]
rownames(pval)<-counts[-1,2]


options(scipen = 6)

pval

write.csv(pval, file="./BphP1DepDEGs_function_pval.csv")

source("./jabes.txt")


qvalues=apply(abs(pval),2,jabes.q)*sign(pval)


qvalues

write.csv(qvalues, file="./BphP1DepDEGs_function_qval.csv")
