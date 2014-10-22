#source("ggplot-heatmap.R")
source("heatmap.3.R")
#library(gplots)
library(vegan)
library(RColorBrewer)
#library(gplots)
#library(DESeq2)
#library(fdrtool)
#library(MASS)
library(ggplot2)
#library(reshape)

columns_to_use <- 29

pdf("./13.figures.heatmap.pdf")

d <- read.table("12.persistent_otus.xls", header=T, row.names=1, sep="\t")

raw_otus <- d[,1:columns_to_use]
head(raw_otus)

sample_read_sums <- read.table("13.figures.read_sums.xls", header=T, row.names=1, sep="\t")
sample_read_sums

# normalize the raw_otu matrix to % of reads
otus_pcnt <- scale(raw_otus, center=FALSE, scale=t(sample_read_sums))
head(otus_pcnt)

write.table(otus_pcnt, file="13.figures.otus_percent.tab", quote=F, sep="\t")

#heatmap.3(as.matrix(otus_pcnt), as.matrix(otus_pcnt), key=T, distfun=vegdist, col=brewer.pal(9, "YlGnBu"), hclustfun = function(x) hclust(x,method = 'average'), scale='none', na.rm=T, na.color="white", density.info='none', trace='none', KeyValueName="Percent of reads")

# normalize the otu maxtrix by the maximum number of reads
max_seqs <- max(sample_read_sums)
otus <- scale(raw_otus, center=FALSE, scale=t(sample_read_sums))
otus <- otus * max_seqs
colSums(otus)
head(otus)

log_otus <- log(otus)
# artificial baseline of the 0 entries
log_otus[log_otus == -Inf] <- NA
head(log_otus)

write.table(log_otus, file="13.figures.log_otus.tab", quote=F, sep="\t")

heatmap.3(as.matrix(log_otus), as.matrix(otus_pcnt), key=T, distfun=vegdist, col=brewer.pal(9, "YlGnBu"), hclustfun = function(x) hclust(x,method = 'average'), scale='none', na.rm=T, na.color="white", density.info='none', trace='none', KeyValueName="ln Reads", symkey=F, symbreaks=F, sepcolor="lightgrey", rowsep=c(1:23), colsep=c(1:35))

log_otus_high <- log_otus[,c("sediment.rep.1", "H1.1", "H1.2", "H1.3", "H1.4", "H1.5", "H1.6", "H1.7", "H1.8", "H1.9", "H1.20", "H1.21", "H1.23")]
heatmap.3(as.matrix(log_otus_high), as.matrix(log_otus_high), key=T, Rowv=F, Colv=F, col=brewer.pal(9, "YlGnBu"), scale='none', na.rm=T, na.color="white", density.info='none', trace='none', KeyValueName="Log10 Reads", symkey=F, symbreaks=F, sepcolor="lightgrey", rowsep=c(1:23), colsep=c(1:35))

log_otus_low <- log_otus[,c("sediment.rep.1", "L1.1", "L1.2", "L1.3", "L1.4", "L1.5", "L1.6", "L1.9", "L1.10", "L1.11", "L1.20", "L1.21", "L1.23")]
heatmap.3(as.matrix(log_otus_low), as.matrix(log_otus_low), key=T, Rowv=F, Colv=F, col=brewer.pal(9, "YlGnBu"), scale='none', na.rm=T, na.color="white", density.info='none', trace='none', KeyValueName="Log10 Reads", symkey=F, symbreaks=F, sepcolor="lightgrey", rowsep=c(1:23), colsep=c(1:35))

# stacked histogram
library(latticeExtra)

otus_reshaped <- read.table("13.figures.persistent_otus.reshaped.tab", sep='\t', quote='', comment.char='', header=F)
names(otus_reshaped) <- c("expt", "week", "sample", "otu", "reads", "log10_reads")
otus_reshaped <- cbind(otus_reshaped, otus_reshaped$reads)
colnames(otus_reshaped)[7]<-"percent"
head(otus_reshaped)

for(sample in levels(otus_reshaped$sample)) {
	otus_reshaped$percent[otus_reshaped$sample == sample] <- otus_reshaped$reads[otus_reshaped$sample == sample] / sample_read_sums[[sample]]
}

#pl <- cloud(reads ~ week * otu | expt, data=otus_reshaped, screen=list(z=-40, x=-25), zoom=1.1, col.facet="grey", xbase=0.6, ybase=0.6, par.settings=list(box.3d=list(col="transparent")), aspect=c(1.5, 0.75), panel.aspect=0.75, panel.3d.cloud=panel.3dbars)
#print(pl)

otus_reshaped_high <- otus_reshaped[otus_reshaped$expt == "H1", ]
head(otus_reshaped_high)
#ggplot(otus_reshaped_high, aes(x = factor(week), y=percent, fill=otu)) + geom_bar(stat = "identity")

otus_reshaped_low <- otus_reshaped[otus_reshaped$expt == "L1", ]
head(otus_reshaped_low)
#ggplot(otus_reshaped_low, aes(x = factor(week), y=percent, fill=otu)) + geom_bar(stat = "identity")

dev.off()

# diversity estimates
# load the data from variableO2_illumina_analyses.xls
d <- read.table("variableO2_illumina_analyses.xls", sep='\t', quote='', comment.char='', header=T, row.names=1);
otus <- d[,1:columns_to_use]
head(otus)

otus.t <- t(otus)

pdf("13.figures.change.pdf", width=3.5, height=3.5);

# change over time
inlp <- read.table("data/otus.pop.dat.combined_replicates.reshape", sep='\t', quote='', comment.char='', header=F);
names(inlp) <- c("Sample", "Weeks", "OTUs", "stdev", "replicates");
levels(inlp$Sample) <- c("HO", "LO")
ggplot(inlp, aes(x=Weeks, y=OTUs, colour=Sample)) +
	geom_errorbar(aes(ymin=OTUs-stdev, ymax=OTUs+stdev), width=1) +
	theme_bw() + theme(legend.position = c(.75, .75)) +
    geom_point(shape=1) + geom_line()

insh <- read.table("data/otus.shannon.dat.combined_replicates.reshape", sep='\t', quote='', comment.char='', header=F);
names(insh) <- c("Sample", "Weeks", "Shannon", "stdev", "replicates");
levels(insh$Sample) <- c("HO", "LO")
ggplot(insh, aes(x=Weeks, y=Shannon, colour=Sample)) +
	geom_errorbar(aes(ymin=Shannon-stdev, ymax=Shannon+stdev), width=.1) +
	theme_bw() + theme(legend.position = c(.75, .75)) +
    geom_point(shape=1) + geom_line()


dev.off()

