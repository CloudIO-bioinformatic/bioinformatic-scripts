#HEATMAP
library(ggplot2)
#set.seed(1)
data <- read.table("reverse_proteinas_graficar_en_R.csv",header=FALSE,sep="\t")
data2 <- read.table("proteinas_graficar_en_R.csv",header=FALSE,sep="\t")
data$V1 <- as.character(data$V1)
data$V1 <- factor(data$V1, levels=unique(data$V1))
data$V3 <- factor(data2$V3, levels=unique(data2$V3))

#agregando info de estado
data$V4 <- ""
data[which(data$V1=="otx"),]$V4 <- "Maternal"
data[which(data$V1=="beta-catenin"),]$V4 <- "Maternal"
data[which(data$V1=="wnt6"),]$V4 <- "Maternal"
data[which(data$V1=="hnf6"),]$V4 <- "PMC"
data[which(data$V1=="pmar1"),]$V4 <- "PMC"
data[which(data$V1=="nBeta-TCF"),]$V4 <- "PMC"
data[which(data$V1=="blimp1/krox"),]$V4 <- "PMC"
data[which(data$V1=="hesc"),]$V4 <- "PMC"
data[which(data$V1=="alx1"),]$V4 <- "PMC"
data[which(data$V1=="ets1"),]$V4 <- "PMC"
data[which(data$V1=="tbr"),]$V4 <- "PMC"
data[which(data$V1=="tgif"),]$V4 <- "PMC"
data[which(data$V1=="erg"),]$V4 <- "PMC"
data[which(data$V1=="hex"),]$V4 <- "PMC"
data[which(data$V1=="tel"),]$V4 <- "PMC"
data[which(data$V1=="dri"),]$V4 <- "PMC"
data[which(data$V1=="ficolin"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="cyp"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r1"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r2"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r3"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r4"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r5"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r6"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="msp130r7"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm32"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm50"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm29"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm30a"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm30b"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm30c"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm30d"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="sm30e"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="c-lect13"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="c-lect14"),]$V4 <- "Skeletogenesis"
data[which(data$V1=="pm27"),]$V4 <- "Skeletogenesis"

png(filename = "heatmap.png", width = 900, height = 1020)
ggplot(data, aes(x=V3, y=V1, fill=log(V2+1)))+
geom_tile(color = "white",lwd = 0.7,linetype = 1)+
theme_classic()+
scale_fill_gradient(low = "lightblue", high = "red")+
labs(y = "", x="",title="",fill="Log(TPM+1)")+
theme(axis.text.x = element_text(angle = 35, hjust = 1, size=12, color="black"))+
theme(axis.text.y = element_text(hjust = 1, size=12, color="black"))+
theme(legend.text = element_text(size=12))+
guides(fill = guide_colourbar(barwidth = 1,barheight = 15))+
facet_grid(V4 ~ ., scales = "free", space = "free")+
theme(strip.text.y = element_text(angle=90,size = 12))
dev.off()

#heatmap con tpm promedio
blastula <- NULL
blastula$gene <- data[which(data$V3 == "Blastula_1"),]$V1
blastula$stage <- data[which(data$V3 == "Blastula_1"),]$V4
blastula$tpm <- (data[which(data$V3 == "Blastula_1"),]$V2 + data[which(data$V3 == "Blastula_2"),]$V2 / 2)
blastula$name <- "Blastula"
blastula <- data.frame(blastula)

gastrula <- NULL
gastrula$gene <- data[which(data$V3 == "Gastrula_1"),]$V1
gastrula$stage <- data[which(data$V3 == "Gastrula_1"),]$V4
gastrula$tpm <- (data[which(data$V3 == "Gastrula_1"),]$V2 + data[which(data$V3 == "Gastrula_2"),]$V2 / 2)
gastrula$name <- "Gastrula"
gastrula <- data.frame(gastrula)

prismatic <- NULL
prismatic$gene <- data[which(data$V3 == "Prismatic__larva_1"),]$V1
prismatic$stage <- data[which(data$V3 == "Prismatic__larva_1"),]$V4
prismatic$tpm <- (data[which(data$V3 == "Prismatic__larva_1"),]$V2 + data[which(data$V3 == "Prismatic__larva_2"),]$V2 / 2)
prismatic$name <- "Prismatic"
prismatic <- data.frame(prismatic)

pluteus <- NULL
pluteus$gene <- data[which(data$V3 == "Pluteus__larva_1"),]$V1
pluteus$stage <- data[which(data$V3 == "Pluteus__larva_1"),]$V4
pluteus$tpm <- (data[which(data$V3 == "Pluteus__larva_1"),]$V2 + data[which(data$V3 == "Pluteus__larva_2"),]$V2 / 2)
pluteus$name <- "Pluteus"
pluteus <- data.frame(pluteus)

newdata2 <- rbind(blastula, gastrula, prismatic, pluteus)
newdata <- rbind(pluteus, prismatic, gastrula, blastula)
newdata$gene <- as.character(newdata$gene)
newdata$gene <- factor(newdata$gene, levels=unique(newdata$gene))
newdata$name <- factor(newdata2$name, levels=unique(newdata2$name))

png(filename = "heatmap_mean.png", width = 900, height = 1020)
ggplot(newdata, aes(x=name, y=gene, fill=log(tpm+1)))+
geom_tile(color = "white",lwd = 0.7,linetype = 1)+
theme_classic()+
scale_fill_gradient(low = "lightblue", high = "red")+
labs(y = "", x="",title="",fill="Log(TPM+1)")+
theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=20, color="black"))+
theme(axis.text.y = element_text(hjust = 1, size=15, color="black"))+
theme(legend.text = element_text(size=15))+
theme(legend.title = element_text(size=20))+
guides(fill = guide_colourbar(barwidth = 1,barheight = 15))+
facet_grid(stage ~ ., scales = "free", space = "free")+
theme(strip.text.y = element_text(angle=90,size = 15))
dev.off()
#GEOM_LINE
q#layout(matrix(c(1,1,2,3,4), 2, 2, byrow = TRUE),
 #  widths=c(3,1), heights=c(1,2))
data2 <- read.table("proteinas_graficar_en_R.csv",header=FALSE,sep="\t")
data2$V3 <- factor(data2$V3, levels=unique(data2$V3))
data2$V1 <- factor(data$V1, levels=unique(data$V1))
otx <- data2[which(data2$V1 == "gene-otx"),]
betacatenin <- data2[which(data2$V1 == "gene-beta-catenin"),]
wnt6 <- data2[which(data2$V1 == "gene-wnt6"),]
png(filename = "geomline-genes-maternal.png", width = 900, height = 250)
ggplot(otx, aes(x=V3,y=log(V2+1),group=V1,color =V1))+geom_line()+geom_line(data=betacatenin,aes(x=V3,y=log(V2+1),group=V1,color =V1))+geom_line(data=wnt6,aes(x=V3,y=log(V2+1),group=V1,color =V1))+labs(y = "Log(TPM+1)", size=15, x="",title="Maternal stage",colour="Gene")+theme_classic()
dev.off()
















#RADARCHART

library(ggpubr)
library(fmsb)
library(ggplot2)
midata <- read.table("proteinas_biomineralizacion_Spurv5_TPM_names.csv",header=TRUE, sep="\t")
png(filename = "radarchart-genes-maternal.png", width = 900, height = 450)
op <- par(mfrow=c(1, 3))
rownames(midata) <- midata$gene
midata$gene <- NULL
gen1 <- midata[1:1,]
max <- apply(gen1,1,max)
gen1 <- rbind(rep(max,5) , rep(0,5) , gen1)
gen_title <- as.character(rownames(midata[1:1,]))
radarchart( gen1[,c(1, 8:2)],pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.4),cglcol="black",axistype=1, plwd=4 , plty=1, cglty=1, axislabcol="black", caxislabels=seq(0,max,max/4),cglwd=2,vlcex=2,title=gen_title,cex.main = 3)


gen2 <- midata[2:2,]
max <- apply(gen2,1,max)
gen2 <- rbind(rep(max,5) , rep(0,5) , gen2)
gen_title <- as.character(rownames(midata[2:2,]))
radarchart( gen2[,c(1, 8:2)],pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.4),cglcol="black",axistype=1, plwd=4 , plty=1, cglty=1, axislabcol="black", caxislabels=seq(0,max,max/4),cglwd=2,vlcex=2,title=gen_title,cex.main = 3)

gen3 <- midata[3:3,]
max <- apply(gen3,1,max)
gen3 <- rbind(rep(max,5) , rep(0,5) , gen3)
gen_title <- as.character(rownames(midata[3:3,]))

radarchart( gen3[,c(1, 8:2)],pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.4),cglcol="black",axistype=1, plwd=4 , plty=1, cglty=1, axislabcol="black", caxislabels=seq(0,max,max/4),cglwd=2,vlcex=2,title=gen_title,cex.main = 3)
par(op)
dev.off()

