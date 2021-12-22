library(ChIPseeker)
#library(clusterProfiler)
library('GenomicFeatures')
library(ggupset)
library(ggplotify)
library(ggimage)

#creando txdb de S. purpuratus
txdb_spur <- makeTxDbFromGFF("/kraken/cquevedo/memoria/data-atacseq/ngmerged/spurpuratus_v5_GCF.gff3", dataSource="Echinobase", organism="Strongylocentrotus purpuratus")

#leyendo peaks
file <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.01/genome-wide_peaks/idr-blastula"

peak <- readPeakFile(file)


#png(filename = paste("chip1",".png"), width = 1920/2, height = 1080/2)
#covplot(peak, weightCol="X1000")
#dev.off()


promoter <- getPromoters(TxDb=txdb_spur, upstream=3500, downstream=3500)

tagMatrix <- getTagMatrix(peak, windows=promoter)

png(filename = paste("chip2",".png"), width = 500, height = 1080)
tagHeatmap(tagMatrix, xlim=c(-3500, 3500), color="red")
dev.off()

png(filename = paste("chip3",".png"), width = 1920/2, height = 1080/2)
plotAvgProf(tagMatrix, xlim=c(-3500, 3500),xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()
##estos
png(filename = paste("chip4",".png"), width = 400, height = 800)
plotAvgProf(tagMatrix, xlim=c(-2000, 2000), conf = 0.95, resample = 1000)
dev.off()

peakAnno <- annotatePeak(file, tssRegion=c(-3500,3500),TxDb=txdb_spur)

png(filename = paste("chip5",".png"), width = 1920/2, height = 1080/2)
plotAnnoPie(peakAnno)
dev.off()

png(filename = paste("chip6",".png"), width = 1920/2, height = 1080/2)
plotAnnoBar(peakAnno)
dev.off()

png(filename = paste("chip7",".png"), width = 1920/2, height = 1080/2)
upsetplot(peakAnno)
dev.off()

png(filename = paste("chip8",".png"), width = 1920/2, height = 1080/2)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()




#########################
#IDR 0.05

library(ChIPseeker)
library('GenomicFeatures')
library(ggupset)
library(ggplotify)
library(ggimage)

#creando txdb de S. purpuratus
txdb_spur <- makeTxDbFromGFF("/kraken/cquevedo/memoria/data-atacseq/ngmerged/spurpuratus_v5_GCF.gff3", dataSource="Echinobase", organism="Strongylocentrotus purpuratus")

#rutas de peaks
file1 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.05/filtered_peaks_idr/idr-blastula-filtered.narrowPeak"
file2 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.05/filtered_peaks_idr/idr-gastrula-filtered.narrowPeak"
file3 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.05/filtered_peaks_idr/idr-prismatica-filtered.narrowPeak"
file4 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.05/filtered_peaks_idr/idr-pluteus-filtered.narrowPeak"

files <- NULL

files$blastula <- file1
files$gastrula <- file2
files$prismatic <- file3
files$pluteus <- file4

#peak <- readPeakFile(files[[4]])
promoter <- getPromoters(TxDb=txdb_spur, upstream=2000, downstream=2000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb_spur,tssRegion=c(-2000, 2000), verbose=FALSE)

#comparativo
png(filename = "idr-0.05-comparative1.png", width = 1600, height = 800)
plotAvgProf(tagMatrixList, xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = "idr-0.05-comparative4.png", width = 900, height = 600)
plotAnnoBar(peakAnnoList)
dev.off()

png(filename = "idr-0.05-comparative5.png", width = 900, height = 600)
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci relative to TSS")
dev.off()

#individual linea
png(filename = paste("idr-0.05-blastula1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[1]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = paste("idr-0.05-gastrula1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[2]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = paste("idr-0.05-prismatica1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[3]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = paste("idr-0.05-pluteus1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[4]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

#individual pie

png(filename = "idr-0.05-blastula2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[1]])
dev.off()

png(filename = "idr-0.05-gastrula2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[2]])
dev.off()

png(filename = "idr-0.05-prismatica2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[3]])
dev.off()

png(filename = "idr-0.05-pluteus2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[4]])
dev.off()



######################################################
### IDR 0.01

library(ChIPseeker)
library('GenomicFeatures')
library(ggupset)
library(ggplotify)
library(ggimage)

#creando txdb de S. purpuratus
txdb_spur <- makeTxDbFromGFF("/kraken/cquevedo/memoria/data-atacseq/ngmerged/spurpuratus_v5_GCF.gff3", dataSource="Echinobase", organism="Strongylocentrotus purpuratus")

#rutas de peaks
file1 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.01/peaks/filtered_peak_idr/idr-blastula-filtered.narrowPeak"
file2 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.01/peaks/filtered_peak_idr/idr-gastrula-filtered.narrowPeak"
file3 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.01/peaks/filtered_peak_idr/idr-prismatica-filtered.narrowPeak"
file4 <- "/kraken/cquevedo/memoria/data-atacseq/ngmerged/result/fseq2_10000_p0.01/peaks/filtered_peak_idr/idr-pluteus-filtered.narrowPeak"

files <- NULL

files$blastula <- file1
files$gastrula <- file2
files$prismatic <- file3
files$pluteus <- file4

#peak <- readPeakFile(files[[4]])
promoter <- getPromoters(TxDb=txdb_spur, upstream=2000, downstream=2000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb_spur,tssRegion=c(-2000, 2000), verbose=FALSE)

#comparativo
png(filename = "idr-0.01-comparative1.png", width = 1600, height = 800)
plotAvgProf(tagMatrixList, xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()


png(filename = "idr-0.01-comparative2.png", width = 1600, height = 800)
plotAvgProf2(files, upstream=2000,downstream=2000,TxDb=txdb_spur,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = "idr-0.01-comparative4.png", width = 900, height = 600)
plotAnnoBar(peakAnnoList)
dev.off()

png(filename = "idr-0.01-comparative5.png", width = 900, height = 600)
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci relative to TSS")
dev.off()

#individual linea
png(filename = paste("idr-0.01-blastula1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[1]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = paste("idr-0.01-gastrula1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[2]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = paste("idr-0.01-prismatica1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[3]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

png(filename = paste("idr-0.01-pluteus1.png",".png"), width = 400, height = 800)
plotAvgProf(tagMatrixList[[4]], xlim=c(-2000, 2000), conf = 0.95, resample = 500,xlab = "Genomic Region (5'->3')",ylab = "Peak Count Frequency",facet = "column")
dev.off()

#individual pie

png(filename = "idr-0.01-blastula2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[1]])
dev.off()

png(filename = "idr-0.01-gastrula2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[2]])
dev.off()

png(filename = "idr-0.01-prismatica2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[3]])
dev.off()

png(filename = "idr-0.01-pluteus2.png", width = 900, height = 1200)
plotAnnoPie(peakAnnoList[[4]])
dev.off()


#guardar anotaciones

write.table(as.data.frame(peakAnnoList[[1]]@anno), file="blastula_peakAnno.txt", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(peakAnnoList[[2]]@anno), file="gastrula_peakAnno.txt", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(peakAnnoList[[3]]@anno), file="prismatica_peakAnno.txt", quote = F, row.names = F, sep = "\t")
write.table(as.data.frame(peakAnnoList[[4]]@anno), file="pluteus_peakAnno.txt", quote = F, row.names = F, sep = "\t")






