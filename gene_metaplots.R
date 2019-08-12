# get gene coordinates
setwd("/home/springer/nosha003/database/B73v4")
gene <- read.delim("Zea_mays.AGPv4.32.gene.sort2.gff3", header=F, sep="\t", stringsAsFactors = F)
gene$geneID <- substr(gene$V9, 9, 22)

setwd("/home/springer/nosha003/sarah/TEpolymorphism")
category <- read.delim("genes_lists_for_metaplots.txt", header=T, sep="\t", stringsAsFactors = F)

library(dplyr)
library(tidygenomics)
gene_category <- left_join(category, gene, by="geneID")
gene_gff <- gene_category[,c(3:11, 1:2)]
gene_gff_na <- na.omit(gene_gff)
setwd("/home/springer/nosha003/sarah/TEpolymorphism")
write.table(gene_gff_na, "genes_lists_for_metaplots_coordinates.gff3", quote=F, row.names=F, sep="\t")


# unix --> closest gene to each 100bp bin
sed '1d' /home/springer/nosha003/sarah/TEpolymorphism/genes_lists_for_metaplots_coordinates.gff3 | sort -k 1,1 -k 4,4n > /home/springer/nosha003/sarah/TEpolymorphism/genes_lists_for_metaplots_coordinates_sort.gff3

module load bedtools
bedtools closest -D b -t first -k 1 -b /home/springer/nosha003/sarah/TEpolymorphism/genes_lists_for_metaplots_coordinates_sort.gff3 -a /home/springer/nosha003/wgbs_schmitz/tiles_output/B73_V2_100bp_ratio.txt > /home/springer/nosha003/sarah/TEpolymorphism/genes_lists_for_metaplots_bin.txt

# dataframe
perl /home/springer/nosha003/methylation_spreading/metaplot_gene.pl -i /home/springer/nosha003/sarah/TEpolymorphism/genes_lists_for_metaplots_bin.txt -o /scratch.global/nosha003/genes_lists_for_metaplots_bin_df.txt


# metaplots
awk '{if ($15 == "syntenic") print $0}' /scratch.global/nosha003/genes_lists_for_metaplots_bin_df.txt > /scratch.global/nosha003/genes_lists_for_metaplots_bin_df_syntenic.txt
awk '{if ($15 == "non.syntenic") print $0}' /scratch.global/nosha003/genes_lists_for_metaplots_bin_df.txt > /scratch.global/nosha003/genes_lists_for_metaplots_bin_df_nonsyntenic.txt
awk '{if ($15 == "in.shared.TE") print $0}' /scratch.global/nosha003/genes_lists_for_metaplots_bin_df.txt > /scratch.global/nosha003/genes_lists_for_metaplots_bin_df_sharedTE.txt
awk '{if ($15 == "in.variable.TE") print $0}' /scratch.global/nosha003/genes_lists_for_metaplots_bin_df.txt > /scratch.global/nosha003/genes_lists_for_metaplots_bin_df_varialeTE.txt

# generate metaplot of data 
# qsub -I -l walltime=5:00:00,nodes=1:ppn=6,mem=60G

library(fields)
setwd("/scratch.global/nosha003/")
#setwd("/Volumes/CBS/Groups/LAB-springer/Jaclyn/Graduate_projects/methylation/spreading/")
q1 <- read.table(file="genes_lists_for_metaplots_bin_df_syntenic.txt", header=F, sep="\t")
q2 <- read.table(file="genes_lists_for_metaplots_bin_df_nonsyntenic.txt", header=F, sep="\t")
q3 <- read.table(file="genes_lists_for_metaplots_bin_df_sharedTE.txt", header=F, sep="\t")
q4 <- read.table(file="genes_lists_for_metaplots_bin_df_varialeTE.txt", header=F, sep="\t")
colnames(q1) <- c("chr", "binstart", "binstop", "binid", "cg", "chg", "chh", "TE", "TEstart", "TEend", "TEsize", "distance", "relative_distance", "real_distance", "category") 
colnames(q2) <- c("chr", "binstart", "binstop", "binid", "cg", "chg", "chh", "TE", "TEstart", "TEend", "TEsize", "distance", "relative_distance", "real_distance", "category") 
colnames(q3) <- c("chr", "binstart", "binstop", "binid", "cg", "chg", "chh", "TE", "TEstart", "TEend", "TEsize", "distance", "relative_distance", "real_distance", "category") 
colnames(q4) <- c("chr", "binstart", "binstop", "binid", "cg", "chg", "chh", "TE", "TEstart", "TEend", "TEsize", "distance", "relative_distance", "real_distance", "category") 

look just at flanking 1kb of FGS
data.sub1=subset(q1,q1$relative_distance < 2000 & q1$relative_distance > -1000)
data.sub2=subset(q2,q2$relative_distance < 2000 & q2$relative_distance > -1000)
data.sub3=subset(q3,q3$relative_distance < 2000 & q3$relative_distance > -1000)
data.sub4=subset(q4,q4$relative_distance < 2000 & q4$relative_distance > -1000)

get 100 bins across the relative gene distance and get stats on your data column
stats.test1.cg=stats.bin(data.sub1$relative_distance,data.sub1$cg,N=600)
stats.test1.chg=stats.bin(data.sub1$relative_distance,data.sub1$chg,N=600)
stats.test1.chh=stats.bin(data.sub1$relative_distance,data.sub1$chh,N=600)

stats.test2.cg=stats.bin(data.sub2$relative_distance,data.sub2$cg,N=600)
stats.test2.chg=stats.bin(data.sub2$relative_distance,data.sub2$chg,N=600)
stats.test2.chh=stats.bin(data.sub2$relative_distance,data.sub2$chh,N=600)

stats.test3.cg=stats.bin(data.sub3$relative_distance,data.sub3$cg,N=600)
stats.test3.chg=stats.bin(data.sub3$relative_distance,data.sub3$chg,N=600)
stats.test3.chh=stats.bin(data.sub3$relative_distance,data.sub3$chh,N=600)

stats.test4.cg=stats.bin(data.sub4$relative_distance,data.sub4$cg,N=600)
stats.test4.chg=stats.bin(data.sub4$relative_distance,data.sub4$chg,N=600)
stats.test4.chh=stats.bin(data.sub4$relative_distance,data.sub4$chh,N=600)

p.1.cg=cbind(matrix(stats.test1.cg$centers,ncol=1),stats.test1.cg$stats["mean",])
p.1.chg=cbind(matrix(stats.test1.chg$centers,ncol=1),stats.test1.chg$stats["mean",])
p.1.chh=cbind(matrix(stats.test1.chh$centers,ncol=1),stats.test1.chh$stats["mean",])

p.2.cg=cbind(matrix(stats.test2.cg$centers,ncol=1),stats.test2.cg$stats["mean",])
p.2.chg=cbind(matrix(stats.test2.chg$centers,ncol=1),stats.test2.chg$stats["mean",])
p.2.chh=cbind(matrix(stats.test2.chh$centers,ncol=1),stats.test2.chh$stats["mean",])

p.3.cg=cbind(matrix(stats.test3.cg$centers,ncol=1),stats.test3.cg$stats["mean",])
p.3.chg=cbind(matrix(stats.test3.chg$centers,ncol=1),stats.test3.chg$stats["mean",])
p.3.chh=cbind(matrix(stats.test3.chh$centers,ncol=1),stats.test3.chh$stats["mean",])

p.4.cg=cbind(matrix(stats.test4.cg$centers,ncol=1),stats.test4.cg$stats["mean",])
p.4.chg=cbind(matrix(stats.test4.chg$centers,ncol=1),stats.test4.chg$stats["mean",])
p.4.chh=cbind(matrix(stats.test4.chh$centers,ncol=1),stats.test4.chh$stats["mean",])

setwd("/scratch.global/nosha003/")

# make plots for each gene set with CG, CHG, and CHH on the same plot

pdf("sarah_gene_metaplot_syntenic.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='syntenic')
lines(p.1.cg,col=1,lwd=1)
lines(p.1.chg,col=2,lwd=1)
lines(p.1.chh,col=3,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('CG', 'CHG', 'CHH'),col=c(1,2,3),lty=1,lwd=2,cex=0.7)
dev.off()

pdf("sarah_gene_metaplot_nonsyntenic.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='non.syntenic')
lines(p.2.cg,col=1,lwd=1)
lines(p.2.chg,col=2,lwd=1)
lines(p.2.chh,col=3,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('CG', 'CHG', 'CHH'),col=c(1,2,3),lty=1,lwd=2,cex=0.7)
dev.off()

pdf("sarah_gene_metaplot_shared.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='in.shared.TE')
lines(p.3.cg,col=1,lwd=1)
lines(p.3.chg,col=2,lwd=1)
lines(p.3.chh,col=3,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('CG', 'CHG', 'CHH'),col=c(1,2,3),lty=1,lwd=2,cex=0.7)
dev.off()

pdf("sarah_gene_metaplot_variable.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='in.variable.TE')
lines(p.4.cg,col=1,lwd=1)
lines(p.4.chg,col=2,lwd=1)
lines(p.4.chh,col=3,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('CG', 'CHG', 'CHH'),col=c(1,2,3),lty=1,lwd=2,cex=0.7)
dev.off()


# make plots with all gene sets on the same plot (separate plots for CG, CHG, and CHH)
pdf("sarah_gene_metaplot_cg.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='CG')
lines(p.1.cg,col=1,lwd=1)
lines(p.2.cg,col=2,lwd=1)
lines(p.3.cg,col=3,lwd=1)
lines(p.4.cg,col=4,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('syntenic', 'non.syntenic', 'in.shared.TE', 'in.variable.TE'),col=c(1,2,3,4),lty=1,lwd=2,cex=0.7)
dev.off()

pdf("sarah_gene_metaplot_chg.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,1),xlab="",ylab='read count',main='CHG')
lines(p.1.chg,col=1,lwd=1)
lines(p.2.chg,col=2,lwd=1)
lines(p.3.chg,col=3,lwd=1)
lines(p.4.chg,col=4,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('syntenic', 'non.syntenic', 'in.shared.TE', 'in.variable.TE'),col=c(1,2,3,4),lty=1,lwd=2,cex=0.7)
dev.off()

pdf("sarah_gene_metaplot_chh.pdf")
par(mfrow=c(2,1))
par(mar = rep(2,4))
plot(x=NULL, y=NULL,xlim=c(-1000,2000),ylim=c(0,0.04),xlab="",ylab='read count',main='CHH')
lines(p.1.chh,col=1,lwd=1)
lines(p.2.chh,col=2,lwd=1)
lines(p.3.chh,col=3,lwd=1)
lines(p.4.chh,col=4,lwd=1)
xline(0,lty=2,col='black')
xline(1000,lty=2,col='black')
legend("topright",c('syntenic', 'non.syntenic', 'in.shared.TE', 'in.variable.TE'),col=c(1,2,3,4),lty=1,lwd=2,cex=0.7)
dev.off()
