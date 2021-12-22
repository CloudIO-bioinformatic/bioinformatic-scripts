library(ggplot2)
library(edgeR)
run1 <- read.delim("SRR653788.bam.counts", header=FALSE, sep="\t")
run2 <- read.delim("SRR653789.bam.counts", header=FALSE, sep="\t")
run3 <- read.delim("SRR653790.bam.counts", header=FALSE, sep="\t")
run4 <- read.delim("SRR653791.bam.counts", header=FALSE, sep="\t")
run5 <- read.delim("SRR653792.bam.counts", header=FALSE, sep="\t")
run6 <- read.delim("SRR653793.bam.counts", header=FALSE, sep="\t")
run7 <- read.delim("SRR653794.bam.counts", header=FALSE, sep="\t")
run8 <- read.delim("SRR653795.bam.counts", header=FALSE, sep="\t")


allcounts <- merge(run1,run2,by="V1")
allcounts <- merge(allcounts,run3,by="V1")
allcounts <- merge(allcounts,run4,by="V1")
allcounts <- merge(allcounts,run5,by="V1")
allcounts <- merge(allcounts,run6,by="V1")
allcounts <- merge(allcounts,run7,by="V1")
allcounts <- merge(allcounts,run8,by="V1")

colnames(allcounts) <- c("Gene","SRR653788","SRR653789","SRR653790","SRR653791","SRR653792","SRR653793","SRR653794","SRR653795")

#identificadores deben ser los nombres de fila
row.names(allcounts) <- allcounts$Gene
allcounts$Gene <- NULL

################################################################################
#2.1 seria entre wt de 4wk y wt de 17wk, mutant de 4wk y mutant de 4wk

#creo una nueva variable para quitarle los que sobran, por que no estoy seguro como funciona DGEList
pregunta2.1wt <- allcounts
pregunta2.1wt$SRR653788 <- NULL
pregunta2.1wt$SRR653789 <- NULL
pregunta2.1wt$SRR653792 <- NULL
pregunta2.1wt$SRR653793 <- NULL

pregunta2.1mut <- allcounts
pregunta2.1mut$SRR653790 <- NULL
pregunta2.1mut$SRR653791 <- NULL
pregunta2.1mut$SRR653794 <- NULL
pregunta2.1mut$SRR653795 <- NULL

DGE_tiempo_wt <- DGEList(counts=pregunta2.1wt,group=c(1,1,2,2))
DGE_tiempo_mut <- DGEList(counts=pregunta2.1mut,group=c(1,1,2,2))

DGE_tiempo_wt <- calcNormFactors(DGE_tiempo_wt)
DGE_tiempo_mut <- calcNormFactors(DGE_tiempo_mut)

DGE_tiempo_wt <- estimateDisp(DGE_tiempo_wt)
DGE_tiempo_mut <- estimateDisp(DGE_tiempo_mut)

tiempo_exact_wt <- exactTest(DGE_tiempo_wt)
tiempo_exact_mut <- exactTest(DGE_tiempo_mut)

diffexpresion_tiempo_wt <- topTags(tiempo_exact_wt,n=nrow(pregunta2.1wt))
diffexpresion_tiempo_mut <- topTags(tiempo_exact_mut,n=nrow(pregunta2.1mut))

tabla_tiempo_wt <- diffexpresion_tiempo_wt$table
tabla_tiempo_mut <- diffexpresion_tiempo_mut$table

tabla_tiempo_wt$Category <- "Not significant"
tabla_tiempo_mut$Category <- "Not significant"

tabla_tiempo_wt[which(abs(tabla_tiempo_wt$logFC)>=2 & tabla_tiempo_wt$FDR<=0.05),]$Category <- "Significant"
tabla_tiempo_mut[which(abs(tabla_tiempo_mut$logFC)>=2 & tabla_tiempo_mut$FDR<=0.05),]$Category <- "Significant"

#> length(tabla_tiempo_wt[which(abs(tabla_tiempo_wt$logFC)>=2 & tabla_tiempo_wt$FDR<=0.05),]$logFC)
#[1] 83
#> length(tabla_tiempo_mut[which(abs(tabla_tiempo_mut$logFC)>=2 & tabla_tiempo_mut$FDR<=0.05),]$logFC)
#[1] 662

#graficos de volcano
png("pregunta2_1_wt.png",width=1024,height=720,units="px",bg="white")
ggplot(tabla_tiempo_wt)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))
dev.off()
png("pregunta2_1_mut.png",width=1024,height=720,units="px",bg="white")
ggplot(tabla_tiempo_mut)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))
dev.off()
################################################################################
#2.2 seria entre wt de 4wk y mutant de 4wk, wt de 17wk y mutant de 17wk
pregunta2.2_17wk <- allcounts
pregunta2.2_17wk$SRR653792 <- NULL
pregunta2.2_17wk$SRR653793 <- NULL
pregunta2.2_17wk$SRR653794 <- NULL
pregunta2.2_17wk$SRR653795 <- NULL

pregunta2.2_4wk <- allcounts
pregunta2.2_4wk$SRR653788 <- NULL
pregunta2.2_4wk$SRR653789 <- NULL
pregunta2.2_4wk$SRR653790 <- NULL
pregunta2.2_4wk$SRR653791 <- NULL

DGE_genotipo_17wk <- DGEList(counts=pregunta2.2_17wk,group=c(1,1,2,2))
DGE_genotipo_4wk <- DGEList(counts=pregunta2.2_4wk,group=c(1,1,2,2))

DGE_genotipo_17wk <- calcNormFactors(DGE_genotipo_17wk)
DGE_genotipo_4wk <- calcNormFactors(DGE_genotipo_4wk)

DGE_genotipo_17wk <- estimateDisp(DGE_genotipo_17wk)
DGE_genotipo_4wk <- estimateDisp(DGE_genotipo_4wk)

genotipo_exact_17wk <- exactTest(DGE_genotipo_17wk)
genotipo_exact_4wk <- exactTest(DGE_genotipo_4wk)

diffexpresion_genotipo_17wk <- topTags(genotipo_exact_17wk,nrow(pregunta2.2_17wk))
diffexpresion_genotipo_4wk <- topTags(genotipo_exact_4wk,nrow(pregunta2.2_4wk))

tabla_genotipo_17wk <- diffexpresion_genotipo_17wk$table
tabla_genotipo_4wk <- diffexpresion_genotipo_4wk$table

tabla_genotipo_17wk$Category <- "Not significant"
tabla_genotipo_4wk$Category <- "Not significant"

tabla_genotipo_17wk[which(abs(tabla_genotipo_17wk$logFC)>=2 & tabla_genotipo_17wk$FDR<=0.05),]$Category <- "Significant"
tabla_genotipo_4wk[which(abs(tabla_genotipo_4wk$logFC)>=2 & tabla_genotipo_4wk$FDR<=0.05),]$Category <- "Significant"
#length(tabla_genotipo_17wk[which(abs(tabla_genotipo_17wk$logFC)>=2 & tabla_genotipo_17wk$FDR<=0.05),]$logFC)
#[1] 412
#length(tabla_genotipo_4wk[which(abs(tabla_genotipo_4wk$logFC)>=2 & tabla_genotipo_4wk$FDR<=0.05),]$logFC)
#[1] 16

#graficos de volcano
png("pregunta2_2_17wk.png",width=1024,height=720,units="px",bg="white")
ggplot(tabla_genotipo_17wk)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))
dev.off()
png("pregunta2_2_4wk.png",width=1024,height=720,units="px",bg="white")
ggplot(tabla_genotipo_4wk)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))
dev.off()
################################################################################
#2.2 seria entre mut de 17wk  y wt de 4wk

newallcounts <- allcounts

newallcounts$SRR653790 <- NULL
newallcounts$SRR653791 <- NULL
newallcounts$SRR653792 <- NULL
newallcounts$SRR653793 <- NULL

DGE_pregunta2.3 <- DGEList(counts=newallcounts,group=c(1,1,2,2))
DGE_pregunta2.3 <- calcNormFactors(DGE_pregunta2.3)

DGE_pregunta2.3 <- estimateDisp(DGE_pregunta2.3)
pregunta2.3_exact <- exactTest(DGE_pregunta2.3)
diffexpresion_pregunta2.3 <- topTags(pregunta2.3_exact,n=nrow(newallcounts))
tabla_pregunta2.3 <- diffexpresion_pregunta2.3$table
tabla_pregunta2.3$Category <- "Not significant"
tabla_pregunta2.3[which(tabla_pregunta2.3$logFC>=2 & tabla_pregunta2.3$FDR<=0.05),]$Category <- "Up-regulated"
#> length(tabla_pregunta2.3[which(tabla_pregunta2.3$logFC>=2 & tabla_pregunta2.3$FDR<=0.05),]$logFC)
#[1] 59

#graficos de volcano
png("pregunta2_3.png",width=1024,height=720,units="px",bg="white")
ggplot(tabla_pregunta2.3)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))
dev.off()
####################################################################################
#generando archivos pregunta 2.3
preg2.3wt <- tabla_tiempo_wt[which(tabla_tiempo_wt$logFC>=2 & tabla_tiempo_wt$FDR<=0.05),]
preg2.3mut <- tabla_tiempo_mut[which(tabla_tiempo_mut$logFC>=2 & tabla_tiempo_mut$FDR<=0.05),]
preg2.3sobexpr <- tabla_pregunta2.3[which(tabla_pregunta2.3$logFC>=2 & tabla_pregunta2.3$FDR<=0.05),]
write.table(preg2.3wt, file = "preg2.3wt", quote=FALSE,sep="\t")
write.table(preg2.3mut , file = "preg2.3mut", quote=FALSE,sep="\t")
write.table(preg2.3sobexpr , file = "preg2.3sobexpr", quote=FALSE,sep="\t")

#############################################################################################################################
#DGE_respuesta_fenotipo <- DGEList(counts=allcounts,group=c(1,1,2,2,1,1,2,2))
#calcular la normalización
#DGE_respuesta_tiempo <- calcNormFactors(DGE_respuesta_tiempo)
#DGE_respuesta_fenotipo <- calcNormFactors(DGE_respuesta_fenotipo)

#Estimo la disperción
#DGE_respuesta_tiempo <- estimateDisp(DGE_respuesta_tiempo)
#DGE_respuesta_fenotipo <- estimateDisp(DGE_respuesta_fenotipo)

#ejecuto función exactTest
#tiempo_exact <- exactTest(DGE_respuesta_tiempo)
#fenotipo_exact <- exactTest(DGE_respuesta_fenotipo)

#obtener resultados de expresión diferencial
#diffexpresion_tiempo <- topTags(tiempo_exact,n=nrow(allcounts))
#diffexpresion_fenotipo <- topTags(fenotipo_exact,n=nrow(allcounts))

#traspaso las tablas de información a nuevas variables
#tabla_tiempo <- diffexpresion_tiempo$table
#tabla_fenotipo <- diffexpresion_fenotipo$table

#creo una columna llamada Category, la que inicialmente parte como no significativa.
#tabla_tiempo$Category <- "Not significant"
#tabla_fenotipo$Category <- "Not significant"

#filtro en los datos con los criterios de cambio de expresión
#tabla_tiempo[which(abs(tabla_tiempo$logFC)>=2 & tabla_tiempo$FDR<=0.05),]$Category <- "Significant"
#tabla_fenotipo[which(abs(tabla_fenotipo$logFC)>=2 & tabla_fenotipo$FDR<=0.05),]$Category <- "Significant"

#grafíco
#ggplot(tabla_tiempo)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))

#ggplot(tabla_fenotipo)+geom_point(aes(x=logFC,y=-log10(FDR),fill=Category),shape=21,size=3)+theme_classic()+geom_vline(xintercept=c(-2,2),linetype="dashed")+geom_hline(yintercept=-log10(0.05),linetype="dashed")+scale_y_continuous(expand=c(0,0))+scale_fill_manual(values=c("grey","blue"))


















































