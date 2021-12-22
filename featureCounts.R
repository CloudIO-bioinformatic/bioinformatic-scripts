library(Rsubread)
filename1 <- "SRR653788.bam"
filename2 <- "SRR653789.bam"
filename3 <- "SRR653790.bam"
filename4 <- "SRR653791.bam"
filename5 <- "SRR653792.bam"
filename6 <- "SRR653793.bam"
filename7 <- "SRR653794.bam"
filename8 <- "SRR653795.bam"


#cargando los 8 .bam
featurecounts_output1 <- featureCounts(filename1,annot.inbuilt="mm10")
featurecounts_output2 <- featureCounts(filename2,annot.inbuilt="mm10")
featurecounts_output3 <- featureCounts(filename3,annot.inbuilt="mm10")
featurecounts_output4 <- featureCounts(filename4,annot.inbuilt="mm10")
featurecounts_output5 <- featureCounts(filename5,annot.inbuilt="mm10")
featurecounts_output6 <- featureCounts(filename6,annot.inbuilt="mm10")
featurecounts_output7 <- featureCounts(filename7,annot.inbuilt="mm10")
featurecounts_output8 <- featureCounts(filename8,annot.inbuilt="mm10")


#extrayendo los % de reads asignados a genes
#también se puede hacer usando la info que entrega stat solamente.
#reads_percentn <- c(featurecounts_outputn$stat$filenamen[1]*100/sum(featurecounts_outputn$stat$filenamen))
reads_percent1 <- sum(featurecounts_output1$counts)*100/sum(featurecounts_output1$stat$SRR653788.bam)
reads_percent2 <- sum(featurecounts_output2$counts)*100/sum(featurecounts_output2$stat$SRR653789.bam)
reads_percent3 <- sum(featurecounts_output3$counts)*100/sum(featurecounts_output3$stat$SRR653790.bam)
reads_percent4 <- sum(featurecounts_output4$counts)*100/sum(featurecounts_output4$stat$SRR653791.bam)
reads_percent5 <- sum(featurecounts_output5$counts)*100/sum(featurecounts_output5$stat$SRR653792.bam)
reads_percent6 <- sum(featurecounts_output6$counts)*100/sum(featurecounts_output6$stat$SRR653793.bam)
reads_percent7 <- sum(featurecounts_output7$counts)*100/sum(featurecounts_output7$stat$SRR653794.bam)
reads_percent8 <- sum(featurecounts_output8$counts)*100/sum(featurecounts_output8$stat$SRR653795.bam)


#generando los archivos de counts
write.table(featurecounts_output1$counts, file = paste(filename1,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output2$counts, file = paste(filename2,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output3$counts, file = paste(filename3,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output4$counts, file = paste(filename4,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output5$counts, file = paste(filename5,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output6$counts, file = paste(filename6,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output7$counts, file = paste(filename7,"counts", sep = "."),quote=FALSE,sep="\t")
write.table(featurecounts_output8$counts, file = paste(filename8,"counts", sep = "."),quote=FALSE,sep="\t")
library(ggplot2)
all_reads_percent <- data.frame(c(reads_percent1,reads_percent2,reads_percent3,reads_percent4,reads_percent5,reads_percent6,reads_percent7,reads_percent8))
all_bam <- data.frame(c("SRR653788","SRR653789","SRR653790","SRR653791","SRR653792","SRR653793","SRR653794","SRR653795"))
graph_data <- cbind(all_reads_percent,all_bam)
colnames(graph_data) <- c("genes_reads_percent","bam")


png("pregunta1.png",width=1024,height=720,units="px",bg="white")

ggplot(data=graph_data, aes(y=genes_reads_percent,x=bam,fill=genes_reads_percent)) + geom_bar(stat = "identity")+labs(title='Gráfico de porcentaje de reads asignados a genes, por cada archivo BAM',fill= "",x='Nombres archivos .bam', y='% de reads asignados a genes.')+
theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))+
guides(shape = guide_legend(override.aes = list(size = 5)))
dev.off()

























