library(ggplot2)
library(ggpubr)

# comando para contar identificadores únicos:
#awk 'BEGIN{OFS=FS="\t"}{print $1}' bd.blast |sort|uniq -c|awk 'BEGIN{OFS=FS=" "}($1==1){print $0}'|wc -l
#a <- read.table("a",header=FALSE,sep="\t")
#b <- read.table("b",header=FALSE,sep="\t")
#c <- read.table("c",header=FALSE,sep="\t")
a <- read.table("lns2.graphdata",header=FALSE,sep="\t")
b <- read.table("lns3.graphdata",header=FALSE,sep="\t")
c <- read.table("lns4.graphdata",header=FALSE,sep="\t")
c$V3 <- sample(colors(),12)
b$V3 <- sample(colors(),12)
a$V3 <- sample(colors(),12)
grap_c <- ggplot(c,aes(x=V1,y=V2))+geom_bar(stat="identity",color="black",fill=c$V3)+geom_text(aes(label = V2), vjust = -0.3, size = 3.5)+theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(y = "# of hits", x="",title="LNS4.")+ylim(0, 1250)
grap_b <- ggplot(b,aes(x=V1,y=V2))+geom_bar(stat="identity",color="black",fill=b$V3)+geom_text(aes(label = V2), vjust = -0.3, size = 3.5)+theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(y = "# of hits",x="Databases",title="LNS3.")+ylim(0, 1250)
grap_a <- ggplot(a,aes(x=V1,y=V2))+geom_bar(stat="identity",color="black",fill=a$V3)+geom_text(aes(label = V2), vjust = -0.3, size = 3.5)+theme_classic()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(y = "# of hits",x="",title="LNS2.")+ylim(0, 1250)
figure <- ggarrange(grap_a, grap_b, grap_c,labels = c("A", "B", "C"),ncol = 3, nrow = 1)
#png(filename = "databasegraph.png", width = 900, height = 650)
png(filename = "databasegraphlirima.png", width = 900, height = 650)
figure
dev.off()

#mayor a 4
#sample	crispr	spacers	quality
#1_S1_L001	55	44	4
#1_S1_L001	88	70	3
#1_S1_L001	97	77	2
#H3_S7_L002	0	0	4
#H3_S7_L002	0	0	3
#H3_S7_L002	11	10	2
#H0_S6_L002	0	0	4
#H0_S6_L002	3	2	3
#H0_S6_L002	14	11	2

a <- read.table("crisprdetectgraph.R",header=TRUE, sep="\t")
a$Quality <- paste("Quality: ",a$quality)
crisprgraph <- ggplot(a, aes(sample,  crispr, shape=Quality)) + geom_point(size=4)+theme_classic()+labs(y = "# of direct repeats",x="",title="CRISPR.")+ylim(0, 100)+geom_text(aes(label = crispr), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")
spacer <- ggplot(a, aes(sample,  spacers, shape = Quality)) + geom_point(size=4)+theme_classic()+labs(y = "# of spacers",x="",title="Spacers.")+ylim(0, 100)+geom_text(aes(label = spacers), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")
figure2 <- ggarrange(crisprgraph,spacer,labels = c("A", "B"),ncol = 2, nrow = 1)
png(filename = "crisprdetectgraph.png", width = 800, height = 600)
figure2
dev.off()


#quality4+
#sample	hits	database	type
#1_S1_L001	1	pads	direct_repeat
#1_S1_L001	20	pads	binding_site
#1_S1_L001	51	pads	region_repeat
#H3_S7_L002	0	pads	direct_repeat
#H3_S7_L002	0	pads	binding_site
#H3_S7_L002	0	pads	region_repeat
#H0_S6_L002	0	pads	direct_repeat
#H0_S6_L002	0	pads	binding_site
#H0_S6_L002	0	pads	region_repeat
#1_S1_L001	24	crisprcasdb	direct_repeat
#1_S1_L001	2	crisprcasdb	binding_site
#1_S1_L001	0	crisprcasdb	region_repeat
#H3_S7_L002	0	crisprcasdb	direct_repeat
#H3_S7_L002	0	crisprcasdb	binding_site
#H3_S7_L002	0	crisprcasdb	region_repeat
#H0_S6_L002	0	crisprcasdb	direct_repeat
#H0_S6_L002	0	crisprcasdb	binding_site
#H0_S6_L002	0	crisprcasdb	region_repeat

a <- read.table("padscasdb.R",header=TRUE, sep="\t")


png(filename = "padscasdbgraph.png", width = 800, height = 600)
ggplot(a, aes(sample,  hits, shape=type,colour=database)) + geom_point(size=4)+theme_classic()+labs(y = "# of hits",x="",title="PADS & CRISPR-casdb (with CRISPRDetectQuality>=4).")+ylim(0, 60)+geom_text(aes(label = hits), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")
dev.off()


library(ggplot2)
library(ggpubr)
a <- read.table("crisprlirima.graphdata",header=TRUE,sep="\t")

crisprgraph <- ggplot(a, aes(sample,  crispr)) + geom_bar(stat="identity",,width=0.2,color="black",fill="blue",position = position_dodge(width=0.2))+theme_classic()+labs(y = "# of direct repeats",x="",title="CRISPR (CRISPRDetectQuality >= 4).")+ylim(0, 180)+geom_text(aes(label = crispr), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))
spacer <- ggplot(a, aes(sample,  spacers)) + geom_bar(stat="identity",width=0.2,color="black",fill="red",position = position_dodge(width=0.2))+theme_classic()+labs(y = "# of spacers",x="",title="Spacers (CRISPRDetectQuality >= 4).")+ylim(0, 180)+geom_text(aes(label = spacers), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))
figure3 <- ggarrange(crisprgraph,spacer,labels = c("A)", "B)"),ncol = 2, nrow = 1)
png(filename = "crisprdetectgraphlirima.png", width = 800, height = 600)
figure3
dev.off()


a <- read.table("crisprsalar.graphdata",header=TRUE,sep="\t")

crisprgraph <- ggplot(a, aes(sample,  crispr)) + geom_bar(stat="identity",,width=0.2,color="black",fill="blue",position = position_dodge(width=0.2))+theme_classic()+labs(y = "# of direct repeats",x="",title="CRISPR (CRISPRDetectQuality >= 4).")+ylim(0, 180)+geom_text(aes(label = crispr), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))
spacer <- ggplot(a, aes(sample,  spacers)) + geom_bar(stat="identity",width=0.2,color="black",fill="red",position = position_dodge(width=0.2))+theme_classic()+labs(y = "# of spacers",x="",title="Spacers (CRISPRDetectQuality >= 4).")+ylim(0, 180)+geom_text(aes(label = spacers), vjust = -1.1, size = 3)+geom_hline(yintercept=0,color="red", linetype="dotted")+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))
figure4 <- ggarrange(crisprgraph,spacer,labels = c("A)", "B)"),ncol = 2, nrow = 1)
png(filename = "crisprdetectgraphsalar.png", width = 800, height = 600)
figure4
dev.off()

#figurejuntas <- ggarrange(figure3,figure4,labels = c("Lirima", "Salar"),ncol = 2, nrow = 1)


#tipos de crisp-cas en lns2
#comandos:
#awk 'BEGIN{OFS=FS="\t"}{print $2}' dr_pads.blastx |awk 'BEGIN{OFS=FS=","}{print $3}'|sort|uniq -c > x.graphdata
#library(ggplot2)
#library(ggpubr)
#library(ggradar)
#library(scales)
#library(dplyr)
#a <- data.frame(read.table("tipoCASlns2.graphdata",header=FALSE,sep="\t"))
#a_radar <- a %>% 
 # as_tibble(rownames = "X") %>% 
  #mutate_at(vars(-X), rescale) %>% 
  #tail(1) %>% 
  #select(1:7)

#install.packages("fmsb")
library(fmsb)
a <- data.frame(read.table("tipoCASlns2.graphdata",header=TRUE,sep="\t"))
max_min <- data.frame(group=c(41,0),Cas1=c(41,0), Cas3=c(41,0), Cas4=c(41,0),Cas6=c(41,0),Cas10=c(41,0),Cas10d=c(41,0), Cmr4=c(41,0),Csa3=c(41,0),Csc1=c(41,0), Csc2=c(41,0), Cse1=c(41,0), Csm3=c(41,0))
rownames(a)<-c("LNS2","LNS3","LNS4")
rownames(max_min) <- c("Max", "Min")
df <- rbind(max_min, a)

#crear funcion radarbonito

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 1,
                                        caxislabels = NULL, title = NULL){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "black", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title
  )
}

# Define colors and titles
colors <- c("#00AFBB", "#E7B800", "#FC4E07")
titles <- c("LNS2", "LNS3", "LNS4")
# Reduce plot margin using par()
# Split the screen in 3 parts
op <- par(mar = c(1, 1, 1, 1))
#op <- par(mar = c(2, 2, 2, 2))
par(mfrow = c(1,3))
# Create the radar chart

# Create the radar chart
#for(i in 1:3){
#  create_beautiful_radarchart(
#    data = df[c(1, 2, i+2), ][2:12], caxislabels = c(0, 11, 22, 33, 44),
#    color = colors[i]
#    )
#}
# Add an horizontal legend
#legend(
#  x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
#  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
#  text.col = "black", cex = 1, pt.cex = 1.5
#  )
#par(op)

png(filename = "radarchartlirima.png", width = 800, height = 600)
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart(
  data = df[2:12], caxislabels = c(0, 11, 22, 33, 44),
  color = c("#00AFBB", "#E7B800", "#FC4E07"),title="Types of CAS genes found in Lirima."
)
# Add an horizontal legend
legend(
  x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
  )
par(op)

dev.off()







# para huasco el radar chart :D:D:D:D:D

library(fmsb)
a <- data.frame(read.table("tipoCASsalar.graphdata",header=TRUE,sep="\t"))
max_min <- data.frame(group=c(1,0),Cas1=c(1,0), Cas3=c(1,0), Cas4=c(1,0),Cas6=c(1,0),Cas10=c(1,0),Cas10d=c(1,0), Cmr4=c(1,0),Csa3=c(1,0),Csc1=c(1,0), Csc2=c(1,0), Cse1=c(1,0), Csm3=c(1,0))
rownames(a)<-c("1_S1_L001", "H0_S6_L002", "H3_S7_L002")
rownames(max_min) <- c("Max", "Min")
df <- rbind(max_min, a)

#crear funcion radarbonito

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 1,
                                        caxislabels = NULL, title = NULL){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "black", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title
  )
}

# Define colors and titles
colors <- c("#00AFBB", "#E7B800", "#FC4E07")
titles <- c("1_S1_L001", "H0_S6_L002", "H3_S7_L002")
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(1,3))

png(filename = "radarcharthuasco.png", width = 800, height = 600)
op <- par(mar = c(1, 2, 2, 2))
# Create the radar charts
create_beautiful_radarchart(
  data = df[2:12], caxislabels = c(0, 0.1, 0.4, 0.7, 1),
  color = c("#00AFBB", "#E7B800", "#FC4E07"),title="Types of CAS genes found in Salar de Huasco."
)
# Add an horizontal legend
legend(
  x = "bottom", legend = rownames(df[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
  )
par(op)

dev.off()



#tam assembly

library(ggplot2)
library(ggpubr)
data <- read.table("assembly_size.graphdata",header=TRUE,sep="\t")
graph1 <- ggplot(data, aes(Assembly,  Size.MB.))+ geom_bar(stat="identity",width=0.2,color="black",fill="blue",position = position_dodge(width=0.2))+theme_classic()+labs(y = "Size (MB)",x="",title="Assembly size.")+ylim(0, 150)+geom_text(aes(label = Size.MB.), vjust = -1.1, size = 3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))

graph2 <- ggplot(data, aes(Assembly,  N.contigs))+ geom_bar(stat="identity",width=0.2,color="black",fill="green",position = position_dodge(width=0.2))+theme_classic()+labs(y = "# of contigs",x="",title="Number of contigs per assembly.")+ylim(0, 83000)+geom_text(aes(label = N.contigs), vjust = -1.1, size = 3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))

graph3 <- ggplot(data, aes(Assembly,  Max.size.contig))+ geom_bar(stat="identity",width=0.2,color="black",fill="yellow",position = position_dodge(width=0.2))+theme_classic()+labs(y = "Size of contigs (bp)",x="",title="Maximum contig size per assembly.")+ylim(0, 155000)+geom_text(aes(label = Max.size.contig), vjust = -1.1, size = 3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))

graph4 <- ggplot(data, aes(Assembly,  Contig.classification))+ geom_bar(stat="identity",width=0.2,color="black",fill="red",position = position_dodge(width=0.2))+theme_classic()+labs(y = "% of contig classification",x="",title="Percentage of taxonomic classification of contigs using CAT.")+ylim(0, 100)+geom_text(aes(label = Contig.classification), vjust = -1.1, size = 3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))

graph5 <- ggplot(data, aes(Assembly,  Nucleotide.representation))+ geom_bar(stat="identity",width=0.2,color="black",fill="purple",position = position_dodge(width=0.2))+theme_classic()+labs(y = "% of nucleotide representation",x="",title="Percentage of nucleotide representation in each assembly.")+ylim(0, 100)+geom_text(aes(label = Nucleotide.representation), vjust = -1.1, size = 3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))

graph6 <- ggplot(data, aes(Assembly,  CRISPR))+ geom_bar(stat="identity",width=0.2,color="black",fill="skyblue",position = position_dodge(width=0.2))+theme_classic()+labs(y = "# of CRISPR predict",x="",title="Number of CRISPR predict per assembly.")+ylim(0, 170)+geom_text(aes(label = CRISPR), vjust = -1.1, size = 3)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+theme(axis.title.x = element_text(vjust=-0.5, colour="black", size=rel(1.5))) +
theme(axis.title.y = element_text(vjust=1.5, colour="black", size=rel(1.5)))+
theme (axis.text.x = element_text(colour="black", size=rel(1.5)),axis.text.y = element_text(colour="black", size=rel(1.5), hjust=0.5))


png(filename = "comparisongraph.png", width = 1250, height = 870)
ggarrange(graph1, graph2, graph3,graph4,graph5,graph6,labels = c("A)", "B)", "C)","D)","E)","F)"),ncol = 3, nrow = 2)
dev.off()




































