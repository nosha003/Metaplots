# qsub -I -l walltime=5:00:00,nodes=1:ppn=6,mem=60G

setwd("/home/springer/nosha003/database/B73v4/te")
#b_tir_fam <- read.delim("B73_tir_fam_count20_nonnested_20Dec2018.txt", sep="\t", header=F)
#b_ltr_fam <- read.delim("B73_ltr_fam_count20_nonnested_20Dec2018.txt", sep="\t", header=F)
b_tir_fam <- read.delim("B73.structuralTEv2.2018.12.20.filteredTE.TIRnonnested.gff3", sep="\t", header=F)
b_ltr_fam <- read.delim("B73.structuralTEv2.2018.12.20.filteredTE.LTRnonnested.gff3", sep="\t", header=T)
b_tir_fam$TE <- substr(b_tir_fam$V9, 4, 24)
b_ltr_fam$TE <- substr(b_ltr_fam$V9, 4, 24)

#b_ltr_fam <- read.delim("B73_ltr_fam_count20_nonnested_16Aug2018.txt", header=F, sep="\t")
#b_tir_fam <- read.delim("B73_tir_fam_count20_nonnested_16Aug2018.txt", header=F, sep="\t")

setwd("/home/springer/nosha003/methylation_spreading")
ltr_orientation <- read.delim("B73_ltr_closestTEtoBin_orientation_20Dec2018.txt", sep="\t", header=T)
tir_orientation <- read.delim("B73_tir_closestTEtoBin_orientation_20Dec2018.txt", sep="\t", header=T)
ltr_orientation$TE <- substr(ltr_orientation$TE, 4, 24)
ltr_orientation2 <- unique(ltr_orientation[,c(8,10)])
tir_orientation$TE <- substr(tir_orientation$TE, 4, 24)
tir_orientation2 <- unique(tir_orientation[,c(8,10)])

# CG
library(fields)
library(dplyr)
setwd("/scratch.global/nosha003/methylation_spreading")
data2 <- read.table(file="b73_v2_CG_TIR_metaplot.txt", header=T, sep="\t")
data3 <- read.table(file="b73_v2_CG_LTR_metaplot.txt", header=T, sep="\t")
data2$TE <- substr(data2$gene, 4, 24)
data3$TE <- substr(data3$gene, 4, 24)
data2$fam <- substr(data2$TE, 1, 8)
data3$fam <- substr(data3$TE, 1, 8)

data2_orient <- left_join(data2, tir_orientation2, by="TE")
data3_orient <- left_join(data3, ltr_orientation2, by="TE")

#data2.1 <- subset(data2_orient, data2$fam %in% b_tir_fam$V2)
#data3.1 <- subset(data3_orient, data3$fam %in% b_ltr_fam$V2)
data2.1 <- subset(data2_orient, data2$TE %in% b_tir_fam$TE)
data3.1 <- subset(data3_orient, data3$TE %in% b_ltr_fam$TE)
data2.1 <- data2.1 %>% mutate(dist_orient = ifelse(correction == -1, relative_distance*-1, relative_distance))
data3.1 <- data3.1 %>% mutate(dist_orient = ifelse(correction == -1, relative_distance*-1, relative_distance))

#look just at flanking 1kb of FGS
data.sub2=subset(data2.1,data2.1$dist_orient < 2000 & data2.1$dist_orient > -1000)
data.sub3=subset(data3.1,data3.1$dist_orient < 2000 & data3.1$dist_orient > -1000)

#get 100 bins across the relative gene distance and get stats on your data column
stats.test2=stats.bin(data.sub2$relative_distance,data.sub2$count,N=100)
stats.test3=stats.bin(data.sub3$relative_distance,data.sub3$count,N=100)

p.2=cbind(matrix(stats.test2$centers,ncol=1),stats.test2$stats["mean",])
p.3=cbind(matrix(stats.test3$centers,ncol=1),stats.test3$stats["mean",])

pdf("b73_nonnested_methylation_CG_metaplot_n100_oriented_20Dec.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='CG methylation across TE')
lines(p.2,col="black",lwd=1)
lines(p.3,col="blue",lwd=1)
xline(0,lty=2,col='grey')
xline(1000,lty=2,col='grey')
legend("topright",c('TIR', 'LTR'),col=c("black","blue"),lty=c(1,1,1),lwd=2,cex=0.7)
dev.off()


# CHG
library(dplyr)
library(fields)
setwd("/scratch.global/nosha003/methylation_spreading")
data2 <- read.table(file="b73_v2_CHG_TIR_metaplot.txt", header=T, sep="\t")
data3 <- read.table(file="b73_v2_CHG_LTR_metaplot.txt", header=T, sep="\t")
data2$TE <- substr(data2$gene, 4, 24)
data3$TE <- substr(data3$gene, 4, 24)
data2$fam <- substr(data2$TE, 1, 8)
data3$fam <- substr(data3$TE, 1, 8)

data2_orient <- left_join(data2, tir_orientation2, by="TE")
data3_orient <- left_join(data3, ltr_orientation2, by="TE")

#data2.1 <- subset(data2_orient, data2$fam %in% b_tir_fam$V2)
#data3.1 <- subset(data3_orient, data3$fam %in% b_ltr_fam$V2)
data2.1 <- subset(data2_orient, data2$TE %in% b_tir_fam$TE)
data3.1 <- subset(data3_orient, data3$TE %in% b_ltr_fam$TE)
data2.1$dist_orient <- data2.1$relative_distance
data3.1$dist_orient <- data3.1$relative_distance 

#look just at flanking 1kb of FGS
data.sub2=subset(data2.1,data2.1$dist_orient < 2000 & data2.1$dist_orient > -1000)
data.sub3=subset(data3.1,data3.1$dist_orient < 2000 & data3.1$dist_orient > -1000)

#get 100 bins across the relative gene distance and get stats on your data column
stats.test2=stats.bin(data.sub2$relative_distance,data.sub2$count,N=100)
stats.test3=stats.bin(data.sub3$relative_distance,data.sub3$count,N=100)

p.2=cbind(matrix(stats.test2$centers,ncol=1),stats.test2$stats["mean",])
p.3=cbind(matrix(stats.test3$centers,ncol=1),stats.test3$stats["mean",])

pdf("b73_nonnested_methylation_CHG_metaplot_n100_oriented_20Dec.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='CHG methylation across TE')
lines(p.2,col="black",lwd=1)
lines(p.3,col="blue",lwd=1)
xline(0,lty=2,col='grey')
xline(1000,lty=2,col='grey')
legend("topright",c('TIR', 'LTR'),col=c("black","blue"),lty=c(1,1,1),lwd=2,cex=0.7)
dev.off()



# CHH
library(dplyr)
library(fields)
setwd("/scratch.global/nosha003/methylation_spreading")
data2 <- read.table(file="b73_v2_CHH_TIR_metaplot.txt", header=T, sep="\t")
data3 <- read.table(file="b73_v2_CHH_LTR_metaplot.txt", header=T, sep="\t")
data2$TE <- substr(data2$gene, 4, 24)
data3$TE <- substr(data3$gene, 4, 24)
data2$fam <- substr(data2$TE, 1, 8)
data3$fam <- substr(data3$TE, 1, 8)

data2_orient <- left_join(data2, tir_orientation2, by="TE")
data3_orient <- left_join(data3, ltr_orientation2, by="TE")

data2.1 <- subset(data2_orient, data2$TE %in% b_tir_fam$TE)
data3.1 <- subset(data3_orient, data3$TE %in% b_ltr_fam$TE)
data2.1$dist_orient <- data2.1$relative_distance 
data3.1 <- data3.1 %>% mutate(dist_orient = ifelse(correction == -1, relative_distance*-1, relative_distance))

#look just at flanking 1kb of FGS
data.sub2=subset(data2.1,data2.1$dist_orient < 2000 & data2.1$dist_orient > -1000)
data.sub3=subset(data3.1,data3.1$dist_orient < 2000 & data3.1$dist_orient > -1000)

#get 100 bins across the relative gene distance and get stats on your data column
stats.test2=stats.bin(data.sub2$relative_distance,data.sub2$count,N=100)
stats.test3=stats.bin(data.sub3$relative_distance,data.sub3$count,N=100)

p.2=cbind(matrix(stats.test2$centers,ncol=1),stats.test2$stats["mean",])
p.3=cbind(matrix(stats.test3$centers,ncol=1),stats.test3$stats["mean",])

pdf("b73_nonnested_methylation_CHH_metaplot_n100_oriented_20Dec.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='CHH methylation across TE')
lines(p.2,col="black",lwd=1)
lines(p.3,col="blue",lwd=1)
xline(0,lty=2,col='grey')
xline(1000,lty=2,col='grey')
legend("topright",c('TIR', 'LTR'),col=c("black","blue"),lty=c(1,1,1),lwd=2,cex=0.7)
dev.off()

pdf("b73_nonnested_methylation_CHH_metaplot_n100_oriented_20Dec_scale.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,0.2),xlab="",ylab='read count',main='CHH methylation across TE')
lines(p.2,col="black",lwd=1)
lines(p.3,col="blue",lwd=1)
xline(0,lty=2,col='grey')
xline(1000,lty=2,col='grey')
legend("topright",c('TIR', 'LTR'),col=c("black","blue"),lty=c(1,1,1),lwd=2,cex=0.7)
dev.off()


# /home/springer/nosha003/scripts/bw_mC_metaplot_CGCHG_B73_20Dec.R
