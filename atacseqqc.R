library(ATACseqQC)
library(Rsamtools)

bamfile1 <- BamFile("SRR5310895_nodups_noMT.bam")
bamfile2 <- BamFile("SRR5310901_nodups_noMT.bam")
bamfile3 <- BamFile("SRR5310907_nodups_noMT.bam")
bamfile4 <- BamFile("SRR5310913_nodups_noMT.bam")
bamfile5 <- BamFile("SRR5310904_nodups_noMT.bam")
bamfile6 <- BamFile("SRR5310910_nodups_noMT.bam")

bamfile1.labels <- gsub("_nodups_noMT.bam", "", basename(bamfile1))
bamfile2.labels <- gsub("_nodups_noMT.bam", "", basename(bamfile2))
bamfile3.labels <- gsub("_nodups_noMT.bam", "", basename(bamfile3))
bamfile4.labels <- gsub("_nodups_noMT.bam", "", basename(bamfile4))
bamfile5.labels <- gsub("_nodups_noMT.bam", "", basename(bamfile5))
bamfile6.labels <- gsub("_nodups_noMT.bam", "", basename(bamfile6))


png(filename = paste(bamfile1.labels,".png"), width = 1920/2, height = 1080/2)
fragSize <- fragSizeDist(bamfile1$path, bamfile1.labels)
dev.off()
png(filename = paste(bamfile2.labels,".png"), width = 1920/2, height = 1080/2)
fragSize <- fragSizeDist(bamfile2$path, bamfile1.labels)
dev.off()
png(filename = paste(bamfile3.labels,".png"), width = 1920/2, height = 1080/2)
fragSize <- fragSizeDist(bamfile3$path, bamfile1.labels)
dev.off()
png(filename = paste(bamfile4.labels,".png"), width = 1920/2, height = 1080/2)
fragSize <- fragSizeDist(bamfile4$path, bamfile1.labels)
dev.off()
png(filename = paste(bamfile5.labels,".png"), width = 1920/2, height = 1080/2)
fragSize <- fragSizeDist(bamfile5$path, bamfile1.labels)
dev.off()
png(filename = paste(bamfile6.labels,".png"), width = 1920/2, height = 1080/2)
fragSize <- fragSizeDist(bamfile6$path, bamfile1.labels)
dev.off()




