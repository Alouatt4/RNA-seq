  ############## RNA-seq ############## 
 #              samples              #
#####################################


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("pasilla")
BiocManager::install("pcaExplorer")
BiocManager::install("airway")
install.packages("rmarkdown")
install.packages("markdown")


library(DESeq2)
library(pasilla)
library(pcaExplorer)
library(rmarkdown)
library(markdown)

setwd("")
dir <- ""
files <- grep("*htseq.txt",list.files(),value=TRUE)

# informando arquivo de metadados e coluna das condições
metadata_file <- read.csv(file.path("meta-data-samples-fvar.txt"), sep = "\t", row.names = 1)

phase = metadata_file$phase

sampleTable <- data.frame(sample = rownames(metadata_file),
                          arquivo = files,
                          condicao = phase)

# transformar a tabela em uma matriz (conhecida como objeto)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       dir = dir,
                                       design = ~ condicao)
ddsHTSeq

#################################################################


#pre-filtragem
#cortar linhas que tenham menos de 10 reads
corte <- rowSums(counts(ddsHTSeq)) >= 10
dds_filtrado <- ddsHTSeq[corte,]

#analise de expressao
DEG <- DESeq(dds_filtrado)
res_DEG <- results(DEG)
res_DEG
resultsNames(DEG)

#gerar tabela de resultados
table_DEG <- results(DEG, contrast = c("condicao", "nurse", "forager"))

#exportar tabela DEG
write.table(table_DEG, ".../resultado-deseq2.csv", 
            row.names = TRUE, sep = ",")

#tabela de expressão
vsd <- vst(DEG, blind=FALSE)
table_vsd <- assay(vsd)
head(assay(vsd),3)
write.table(table_vsd, ".../cntgs-vsd-nurse-vs-forager.csv", 
            row.names = TRUE, sep = ",")
#
ntd <- normTransform(DEG)
table_ntd <- assay(ntd)
write.table(table_ntd, "C.../cntgs-ntd-nurse-vs-forager.csv", 
            row.names = TRUE, sep = ",")
#
rld <- rlog(DEG, blind=FALSE)
table_rld <- assay(rld)
write.table(table_rld, ".../cntgs-rld-nurse-vs-forager.csv", 
            row.names = TRUE, sep = ",")

#pca

pcaExplorer(dds = DEG, dst = vsd)



######################## end #######################

#plot manual pca
hclust(dist(t(table_vsd))) -> htree
plot(htree)

pca <- prcomp(t(table_vsd))

pca.dat <- as.data.frame(pca$x)
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)

########## color-test ##########
#pca
display.brewer.pal(n=9, name="RdBu")
colorsPal <- brewer.pal(n=9, name="RdBu")
colors = c(colorsPal[1],colorsPal[2],colorsPal[8],colorsPal[9])

pca.dat.plot <- pca.dat
pca.dat.plot

pca.dat.plot$stages <- c(rep("nurse",3),rep("forager",3))

#plot pca ggplot2
library(ggplot2)
# library(ggpubr)
pca.dat.plot
pca.plot <- ggplot(pca.dat.plot, aes(PC1, PC2,color=stages))+
  geom_point(size=5) +
  geom_text(label = pca.dat.plot$stages, col="black",
            nudge_x = 1, nudge_y = 4) +
  labs(x = paste0("PC1: ",pca.var.percent[1], "%"),
       y= paste0("PC2: ",pca.var.percent[2], "%"))+
  scale_color_manual(values = c("nurse" = colors[1] , "forager" = colors[2])) +
  stat_ellipse(geom = "polygon", aes(fill = stages), type ="t", alpha = 0.2) +  # Adiciona elipses
  theme_bw()
pca.plot


#heatmap | use "t"
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <-  vsd$condicao
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheat <- pheatmap(sampleDistMatrix,
                  clustering_distance_rows=sampleDists,
                  clustering_distance_cols=sampleDists,
                  col=colors)

#volcano-plot
resF1F3 <- results(DEG)
resF1F3
library(EnhancedVolcano)
p1 <- EnhancedVolcano(resF1F3 ,
                      lab = rownames(resF1F3),
                      x = 'log2FoldChange',
                      y = 'pvalue',
                      title = 'Nurse versus Forager',
                      pCutoff = 10e-1, #ver p-value para selecionar os DEGs OU ver na tabela deseq2
                      FCcutoff = 2,
                      pointSize = 3.0,
                      labSize = 6.0)
p1
#VER A SUBSTITUIÇÃO DE "LOC" PELO NOME DOS GENES

# filtrar DEGs
head(resF1F3)
db <- as.data.frame(resF1F3)
head(db)

#omitir os NA da tabela
db2 <- na.omit(db)

summary(db2)

#ajustar o p-value de corte
db_005 <- db2[db2$padj < 0.005 ,]
summary(db_005)

#ajustar o fold-change de corte
db_FC <- db2[db2$log2FoldChange > 2 | db2$log2FoldChange < -2 ,]
db_FC_005 <- db_005[db_005$log2FoldChange > 2 | db_005$log2FoldChange < -2 ,]

#ver info da tabela
dim(db_FC_005)
summary(db_FC_005)
dim(db)

db_FC_005
dim(db2)
dim(db_FC)
summary(db_FC)
dim(db_005)



head(db_005)
db_005
#
#

#coletar genes especificos para analise
exp <- as.data.frame(assay(vsd))
df <- as.data.frame(assay(ddsHTSeq))
listaIDs = c("LOC122528793","LOC122534943","LOC122536424")
head(df)
df[rownames(df) == listaIDs[3],]


exp$gene <- row.names(exp)
rownames(exp) <- NULL
expDF <- data.frame()

##################
#coletar genes especificos para analise
##################
for (i in 1:length(listaIDs)){
  print(i)
  ID = listaIDs[i]
  info =  exp[exp$gene == ID,]
  expDF <- rbind(expDF,info)
}
expDF

##############################################################
#genes selecionados com base no corte "p-value e/ou FoldChange
exp <- as.data.frame(assay(vsd))
head(exp)
dim(exp)
#inserir nome dos genes no gráfico
selected_genes <- rownames(db_FC_005)
selected_genes

exp_selected <- exp[exp$gene %in% selected_genes,]
dim(exp_selected)
exp_selected


m <- t(exp_selected[,1:6])
pheatmap(m)
#nomear os genes no heatmap
colnames(m) <- selected_genes
#plot heatmap
pheatmap(m)



exp[exp$gene == listaIDs[3],]

##################
