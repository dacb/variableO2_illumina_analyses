library(vegan)
library(ggplot2)

columns_to_use <- 30

d <- read.table("variableO2_illumina_analyses.xls", header=T, row.names=1, sep="\t")

raw_otus <- d[,1:columns_to_use]

dmat <- t(raw_otus)

H <- diversity(dmat)

# to generate the next lines, I used:
# printf "5\n15\n25\n50\n75\n" | awk '{ for (i = 10; i <= 15; i+=5) printf("O2_%d_%d <- c(H[\"X%dA.%d\"], H[\"X%dB.%d\"], H[\"X%dC.%d\"])\n", $1, i, $1, i, $1, i, $1, i) }'

O2_5_10 <- c(H["X5A.10"], H["X5B.10"], H["X5C.10"])
O2_5_15 <- c(H["X5A.15"], H["X5B.15"], H["X5C.15"])
O2_15_10 <- c(H["X15A.10"], H["X15B.10"], H["X15C.10"])
O2_15_15 <- c(H["X15A.15"], H["X15B.15"], H["X15C.15"])
O2_25_10 <- c(H["X25A.10"], H["X25B.10"], H["X25C.10"])
O2_25_15 <- c(H["X25A.15"], H["X25B.15"], H["X25C.15"])
O2_50_10 <- c(H["X50A.10"], H["X50B.10"], H["X50C.10"])
O2_50_15 <- c(H["X50A.15"], H["X50B.15"], H["X50C.15"])
O2_75_10 <- c(H["X75A.10"], H["X75B.10"], H["X75C.10"])
O2_75_15 <- c(H["X75A.15"], H["X75B.15"], H["X75C.15"])

df <- data.frame(
	O2 <- factor(c(5, 5, 15, 15, 25, 25, 50, 50, 75, 75)),
	alpha <- c(mean(O2_5_10), mean (O2_5_15), mean(O2_15_10), mean(O2_15_15), mean(O2_25_10), mean(O2_25_15),
			mean(O2_50_10), mean(O2_50_15), mean(O2_75_10), mean(O2_75_15)),
	Week <- factor(c(10, 16, 10, 16, 10, 16, 10, 16, 10, 16)),
	sd <- c(sd(O2_5_10), sd (O2_5_15), sd(O2_15_10), sd(O2_15_15), sd(O2_25_10), sd(O2_25_15),
			sd(O2_50_10), sd(O2_50_15), sd(O2_75_10), sd(O2_75_15))
)

# Define the top and bottom of the errorbars
limits <- aes(ymax = alpha + sd, ymin=alpha - sd)

pdf("15.plot_diversity.pdf")

p <- ggplot(df, aes(colour=Week, y=alpha, x=O2))
p + geom_point() + geom_errorbar(limits, width=0.2) + ylab(expression(paste("Shannon-Weaver ", alpha, " diversity")))

dev.off()
