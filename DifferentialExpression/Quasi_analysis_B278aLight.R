# Purpose: Identify genes that are differentially expressed in Pss B728a and isogenic mutants upon 
# exposure to blue, red, or far-red light, or maintenance in the dark.

# Usage: specify a LABEL DESIGNATION FOR OUTPUT FILES (line 74)
#=============================================================================

################ Load Required Libraries ################
library(xtable)
library(QuasiSeq)


################ Data Import ################
## Import data from csv file and save as .rds
sampleName <- paste0(rep(c("1A" ,"1B" ,"1C" ,"1D" ,"1E" ,"2A" ,"2B" ,"2C" ,"2D" ,"3A" ,"3B" ,"3C" ,"3D" ,"4B" ,"4D","4H" ,"4J" ,"5B" ,"5D" ,"6C" ,"7C" ,"8C"),each=3),'_',rep(paste0('rep',1:3),22))
dat0 <- read.csv('./read_counts.csv',header=FALSE)
colnames(dat0) <- c('GeneID',sampleName)
rownames(dat0) <- dat0[,1]
dat0 <- dat0[,-1]
saveRDS(dat0,file='count_data.rds')

## Define how the data are organized using treatment codes
dat <- readRDS('count_data.rds')
info <- data.frame(geno=rep(c("1A" ,"1B" ,"1C" ,"1D" ,"1E" ,"2A" ,"2B" ,"2C" ,"2D" ,"3A" ,"3B" ,"3C" ,"3D" ,"4B" ,"4D","4H" ,"4J" ,"5B" ,"5D" ,"6C" ,"7C" ,"8C"),each=3),reps=rep(paste0('rep',1:3),22))


################ Defining Pairwise Comparisons ################
## Define strain x treatment pairs to be compared (g1 = strain x treatment in the numerator; g2 = strain x treatment in the denominator)
comp <- data.frame(g1=c("1A","1B" ,"1C" ,"1E" ,rep("1A" ,2) ,"1B", rep(c("2A" ,"2B" ,"2C") ,2) ,"2D" ,rep(c("3A" ,"3B" ,"3C") ,2) ,"3D" ,rep("4B" ,2) ,"4D" ,rep("5B" ,2) ,"5D" ,rep("6C" ,2) ,"7C" ,"8C"),g2=c(rep("1D" ,4) ,"1B" ,rep("1C" ,2) ,rep("2D" ,3) ,"1A" ,"1B" ,"1C" ,"1D" ,rep("3D" ,3) ,"1A" ,"1B" ,"1C" ,"1D" ,"4D" ,"1B" ,"1D" ,"5D" ,"1B" ,"1D" ,"1C" ,"2C" ,rep("1C" ,2)))

nc <- nrow(comp)
Ngeno <- 31


################ Gene selection ################
## Filter genes with average read count of <1 across all samples (exclude from all subsequent analysis)
flag <- rowMeans(dat)>1 # True: 5181
write.table(rownames(dat)[!flag],file='removed_genes.txt',quote=FALSE,sep='\n',row.names=FALSE,col.names=FALSE)


################ Count normalization ################
## Normalize read counts using upper-quartile (UQ) method
counts <- dat[flag,]
size <- apply(counts,2,quantile,0.75)


################ QuasiDE Analysis ################
## Define experimental design factors
trt <- as.factor(info$geno)
Rep <- as.factor(info$rep)

## Define model designs
design.list <- vector("list",nc+1)
design.list[[1]] <- model.matrix(~Rep+trt)
for(k in 1:nc) design.list[[k+1]] <- model.matrix(~Rep+as.factor(gsub(comp[k,1],comp[k,2],trt)))
test.mat <- matrix(c(rep(1,nc),2:(nc+1)),ncol=2)

## Fit the model
fit <- QL.fit(counts,design.list,test.mat,log.offset=log(size),Model='NegBin',print.progress=TRUE,bias.fold.tolerance=1)

##Save model results
saveRDS(fit,'Quasi_fit.rds')

## Visualize model fitting results
res <- QL.results(fit,Plot=TRUE)


############### Examine how many genes have q-values less than 0.05 ############
Qsig <- sapply(res$Q.values,function(qval) colSums(qval<0.05))
rownames(Qsig) <- paste0(comp[,1],'-',comp[,2])
xtable(Qsig,digits=0)


############### Specify label designation for output files ############
append <- 'Date_or_other_identifier'


############ Generate histogram of p-values for each pairwise comparison ############
pdf(paste0('./figures/hist_pval_',append,'.pdf'))
for(i in 1:nrow(comp))
  hist(res$P.values$QLSpline[,i],xlab='P values - QLSpline',main=paste(comp[i,1],'vs',comp[i,2]))
dev.off()


############ Generate volcano plots (natural-log[FC] vs -log10[pval]) for each pairwise comparison ############
## Calculate natural-log[FC] values for each comparison
g1 <- levels(factor(info$geno))
g2 <- g1[order(tolower(g1))]
o1 <- o2 <- rep(NA,nc)
for(i in 1:nc) {
  o1[i] <- which(g2==comp[i,1])+2
  o2[i] <- which(g2==comp[i,2])+2
}
bmat <- fit$coef
part2 <- NULL
for(i in 1:nc){
  if(o2[i]>3) tmp <- bmat[,o1[i]] - bmat[,o2[i]] else tmp <- bmat[,o1[i]]
  part2 <- cbind(part2,tmp)
}
colnames(part2) <- paste0('logfc_',comp[,1],'_vs_',comp[,2])
bmat <- fit$coefficients
logfc <- part2
pmat <- res$P.values$QLSpline
qmat <- res$Q.values$QLSpline
logpval <- -log(pmat,10)

i <- 20

## Generate PDFs of volcano plots
pdf(paste0('./figures/volcano_',append,'.pdf'))
for(i in 1:nc){
  plot(logfc[,i],logpval[,i],xlab='log(FC)',ylab='-log10(Pval)',main=paste(comp[i,1],'vs',comp[i,2]))
  sq <- max(pmat[which(qmat[, i] <= 0.05), i],na.rm=TRUE)
  print(sum(pmat[,i]<=sq))
  abline(h=-log(sq,10),col='red')
}
dev.off()
##



############ Generate output file containing natural-log(FC) and q-values ############
## q-value matrix
part1 <- res$Q.val$QLSpline
colnames(part1) <- paste0('qval_',comp[,1],'_vs_',comp[,2])

## natural-log(FC) matrix
g1 <- levels(factor(info$geno))
g2 <- g1[order(tolower(g1))]
o1 <- o2 <- rep(NA,nc)
for(i in 1:nc) {
  o1[i] <- which(g2==comp[i,1])+2
  o2[i] <- which(g2==comp[i,2])+2
}

bmat <- fit$coef
part2 <- NULL
for(i in 1:nc){
  if(o2[i]>3) tmp <- bmat[,o1[i]] - bmat[,o2[i]] else tmp <- bmat[,o1[i]]
  part2 <- cbind(part2,tmp)
}
colnames(part2) <- paste0('logfc_',comp[,1],'_vs_',comp[,2])

write.table(cbind(part1,part2),file=paste0('./QuasiSeq_results_',append,'.txt'),quote=FALSE,sep='\t')


############ Generate output file containing normalized read counts (log-transformed UQ-normalized) ############
dat <- readRDS('count_data.rds')
flag <- rowMeans(dat)>1 # True: 5181
lc <- log(dat[flag,]+0.05) # Add 0.05 to all counts to avoid errors while taking log transformation
size <- apply(lc,2,quantile,0.75)
out <- t(t(lc)-size)
write.table(out,file=paste0('./Normalized_data_add05.txt'),quote=FALSE,sep='\t')



##### Visual check that output files were saved correctly #####
head(out[,1:3])
lc[1:3,1]-size[1]



