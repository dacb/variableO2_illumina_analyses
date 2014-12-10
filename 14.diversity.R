library(vegan)
#source("http://www.jennajacobs.org/R/rarefaction.txt")
source("rarefaction.txt")

columns_to_use <- 29

d <- read.table("variableO2_illumina_analyses.xls", header=T, row.names=1, sep="\t")

raw_otus <- d[,1:columns_to_use]

dmat <- t(raw_otus)

H <- diversity(dmat)

write.table(H, file="14.Shannon_Weaver.alpha.tab", quote=F, sep="\t")

z <- betadiver(dmat, "z")

print(z)

#write.table(z, file="14.beta.tab", quote=F, sep="\t")

rf <- rarefy(dmat, min(rowSums(dmat)))

print(rf)

pdf("14.rarefaction.pdf")
rfc <- rarefaction(dmat)

rfc$richness
rfc$SE
rfc$subsample

write.table(rfc$richness, file="14.richness.tab", quote=F, sep="\t")
