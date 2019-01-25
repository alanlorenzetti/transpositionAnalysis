#!/usr/bin/Rscript

library(ggplot2) ; theme_set(theme_bw())
library(ggridges)
library(gridExtra)
library(dplyr)
library(scales)

df = read.delim("stats.txt", header=T, stringsAsFactors=F)
k = 2500000
palette=hue_pal()(6)
palette=c(palette,"grey")
palette=rev(palette)

df$barcode[df$barcode == "barcode01"] = "NRC1_0p"
df$barcode[df$barcode == "barcode02"] = "Mutant1_0p"
df$barcode[df$barcode == "barcode03"] = "Mutant2_0p"
df$barcode[df$barcode == "barcode04"] = "Mutant3_0p"
df$barcode[df$barcode == "barcode05"] = "Mutant1_20p"
df$barcode[df$barcode == "barcode06"] = "Mutant2_20p"
df$barcode[df$barcode == "barcode06"] = "unclassified"

df$barcode = factor(df$barcode, levels = rev(c("NRC1_0p", "Mutant1_0p", "Mutant2_0p", "Mutant3_0p", "Mutant1_20p", "Mutant2_20p", "unclassified")))

dfsum = group_by(df, barcode)
dfsum = summarise(dfsum, throughput=sum(length))
dfsum$barcode = factor(dfsum$barcode, levels = rev(c("NRC1_0p", "Mutant1_0p", "Mutant2_0p", "Mutant3_0p", "Mutant1_20p", "Mutant2_20p", "unclassified")))

p1 = ggplot(df, aes(x=as.factor(barcode), fill=barcode)) + geom_bar(color="black") + ylim(c(0,80000)) + coord_flip() + theme(text=element_text(size=12), axis.title.y=element_blank()) + scale_fill_manual(values = palette, guide=F)
p2 = ggplot(df, aes(y=as.factor(barcode), x=log10(length),stat="identity", fill=barcode)) + geom_density_ridges2() + xlim(c(1,5)) + theme(text=element_text(size=12), axis.text.y=element_blank(), axis.title.y=element_blank()) + scale_fill_manual(values = palette, guide=F)
p3 = ggplot(df, aes(y=as.factor(barcode), x=meanQuality, stat="identity", fill=barcode)) + geom_density_ridges2() + xlim(c(20,100)) + theme(text=element_text(size=12), axis.title.y=element_blank()) + scale_fill_manual(values = palette, guide=F)
p4 = ggplot(dfsum, aes(x=barcode, y=throughput/k, fill=barcode)) + geom_col(color="black") + ylim(c(0,200)) + ylab(label="throughput/genomeSize") + coord_flip() + theme(text=element_text(size=12), axis.text.y=element_blank(), axis.title.y=element_blank()) + scale_fill_manual(values = palette, guide=F)

grid.arrange(p1, p2, p3, p4)
