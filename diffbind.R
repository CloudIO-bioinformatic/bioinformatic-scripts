library(DiffBind)
#samples <- read.csv("sea_urchin.csv",header=TRUE,sep="\t")
#urchin <- dba(sampleSheet=samples)
#urchin <- dba.count(urchin, bUseSummarizeOverlaps=FALSE)
#png(filename = "diffbind1.png", width = 900, height = 1200)
#plot(urchin)
#dev.off()

#info <- dba.show(urchin)
#libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
#rownames(libsizes) <- info$ID

#png(filename = "diffbind2.png", width = 900, height = 1200)
#plot(urchin)
#dev.off()


#ahora uno por uno
vs_blastula <- read.csv("vs_blastula.csv",header=TRUE,sep="\t")
vs_gastrula <- read.csv("vs_gastrula.csv",header=TRUE,sep="\t")
vs_prismatica <- read.csv("vs_prismatica.csv",header=TRUE,sep="\t")
vs_pluteus <- read.csv("vs_pluteus.csv",header=TRUE,sep="\t")


urchin_vs_blastula <- dba(sampleSheet=vs_blastula)
urchin_vs_gastrula <- dba(sampleSheet=vs_gastrula)
urchin_vs_prismatica <- dba(sampleSheet=vs_prismatica)
urchin_vs_pluteus <- dba(sampleSheet=vs_pluteus)

urchin_vs_blastula <- dba.count(urchin_vs_blastula, bUseSummarizeOverlaps=FALSE)
urchin_vs_gastrula <- dba.count(urchin_vs_gastrula, bUseSummarizeOverlaps=FALSE)
urchin_vs_prismatica <- dba.count(urchin_vs_prismatica, bUseSummarizeOverlaps=FALSE)
urchin_vs_pluteus <- dba.count(urchin_vs_pluteus, bUseSummarizeOverlaps=FALSE)

#
info <- dba.show(urchin_vs_blastula)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

png(filename = "diffbind_vs_blastula1.png", width = 900, height = 1200)
dba.plotHeatmap(urchin_vs_blastula)
dev.off()

png(filename = "diffbind_vs_blastula2.png", width = 900, height = 1200)
plot(libsizes)
dev.off()
#
info <- dba.show(urchin_vs_gastrula)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

png(filename = "diffbind_vs_gastrula1.png", width = 900, height = 1200)
dba.plotHeatmap(urchin_vs_gastrula)
dev.off()

png(filename = "diffbind_vs_gastrula2.png", width = 900, height = 1200)
plot(libsizes)
#
info <- dba.show(urchin_vs_prismatica)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

png(filename = "diffbind_vs_prismatica1.png", width = 900, height = 1200)
dba.plotHeatmap(urchin_vs_prismatica)
dev.off()

png(filename = "diffbind_vs_prismatica2.png", width = 900, height = 1200)
plot(libsizes)
dev.off()
#
info <- dba.show(urchin_vs_pluteus)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID

png(filename = "diffbind_vs_pluteus1.png", width = 900, height = 1200)
dba.plotHeatmap(urchin_vs_pluteus)
dev.off()

png(filename = "diffbind_vs_pluteus2.png", width = 900, height = 1200)
plot(libsizes)
dev.off()

