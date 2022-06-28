#par(mfrow=c(3,1))
#sessionInfo()
#install.packages("https://cran.r-project.org/src/contrib/Archive/GetoptLong/GetoptLong_0.1.8.tar.gz", repo = NULL)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#library("ComplexHeatmap")

####Normal echantillons k=3
seq1=read.table("normal/3mer/TCG-38-4625-.fa")
seq2=read.table("normal/3mer/TCG-38-4626-.fa")
seq3=read.table("normal/3mer/TCG-38-4627-.fa")
seq4=read.table("normal/3mer/TCG-38-4632-.fa")
seq5=read.table("normal/3mer/TCG-44-2655-.fa")
seq6=read.table("normal/3mer/TCG-44-2657-.fa")
seq7=read.table("normal/3mer/TCG-44-266-.fa")
seq8=read.table("normal/3mer/TCG-44-2662-.fa")
seq9=read.table("normal/3mer/TCG-44-2665-.fa")
seq10=read.table("normal/3mer/TCG-44-2668-.fa")
seq1$V1= as.character(seq1$V1)

matricen3=cbind(seq1$V2,seq2$V2,seq3$V2,seq4$V2,seq5$V2,seq6$V2,seq7$V2,seq8$V2,seq9$V2,seq10$V2)
matricen3 = t(matricen3)
rownames(matricen3)=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10")
colnames(matricen3)=c(seq1$V1)


##tumor echantillon k=3

seqt1=read.table("tumor/3mer/TCG-05-4244-0.fa")
seqt2=read.table("tumor/3mer/TCG-05-4249-0.fa")
seqt3=read.table("tumor/3mer/TCG-05-4250-0.fa")
seqt4=read.table("tumor/3mer/TCG-05-4382-0.fa")
seqt5=read.table("tumor/3mer/TCG-05-4384-0.fa")
seqt6=read.table("tumor/3mer/TCG-05-4389-0.fa")
seqt7=read.table("tumor/3mer/TCG-05-4390-0.fa")
seqt8=read.table("tumor/3mer/TCG-05-4395-0.fa")
seqt9=read.table("tumor/3mer/TCG-05-4396-0.fa")
seqt10=read.table("tumor/3mer/TCG-05-4397-0.fa")

seqt1$V1= as.character(seqt1$V1)

matricet3=cbind(seqt1$V2,seqt2$V2,seqt3$V2,seqt4$V2,seqt5$V2,seqt6$V2,seqt7$V2,seqt8$V2,seqt9$V2,seqt10$V2)
matricet3 = t(matricet3)
rownames(matricet3)=c("seqt1","seqt2","seqt3","seqt4","seqt5","seqt6","seqt7","seqt8","seqt9","seqt10")
colnames(matricet3)=c(seqt1$V1)

#
##normalisation 
#

totaln=rowSums(matricen3)
totalt=rowSums(matricet3)
for (i in 1:nrow(matricen3)){
  for (j in 1:ncol(matricen3)){
    matricen3[i,j]=matricen3[i,j]/totaln[i]  
    matricet3[i,j]=matricet3[i,j]/totalt[i]
  }
}
freq1=matricen3[1,]
freq2=matricen3[2,]
freq3=matricen3[3,]
freq4=matricen3[4,]
freq5=matricen3[5,]
freq6=matricen3[6,]
freq7=matricen3[7,]
freq8=matricen3[8,]
freq9=matricen3[9,]
freq10=matricen3[10,]

protein3=read.table("normal/3mer/gencode/protein_coding.fa")
transcri3=read.table("normal/3mer/gencode/transcri.fa")

names = protein3$V1
protein3 = cbind(protein3$V2)
rownames(protein3) = names
tot = sum(protein3)
for (i in 1:nrow(protein3)){
  protein3[i]=protein3[i]/tot
}
protein3=protein3[order(-protein3),]

names = transcri3$V1
transcri3 = cbind(transcri3$V2)
rownames(transcri3) = names
tot = sum(transcri3)
for (i in 1:nrow(transcri3)){
  transcri3[i]=transcri3[i]/tot
}
transcri3
protein3
echant = rbind(matricen3,matricet3)
freq1=freq1[order(-freq1)]
freq2=freq2[order(-freq2)]
freq3=freq3[order(-freq3)]
freq4=freq4[order(-freq4)]
freq5=freq5[order(-freq5)]
freq6=freq6[order(-freq6)]
freq7=freq7[order(-freq7)]
freq8=freq8[order(-freq8)]
freq9=freq9[order(-freq9)]
freq10=freq10[order(-freq10)]
transcri3=transcri3[order(-transcri3),]




protein3=protein3[order(-protein3),]
protein3
protein3[1:10]
transcri3[1:10]
head(transcri3)
freq1
write.table(freq10, file='test.tsv', quote=FALSE, sep='\t')
head(freq1)
head(freq2)
head(freq3)
head(freq4)
head(freq5)
head(freq6)
head(freq7)
head(freq8)
head(freq9)
head(freq10)

#
##Contrôle fréquence
#
seq1[order(-seq1$V2),]
seq1[order(seq1$V2),]
#ne semble pas correct mais vu que c'est canonique résultat biaisé
#
##passage au log
#
for (i in 1:nrow(matricen3)){
  for (j in 1:ncol(matricen3)){
    matricen3[i,j]=log(matricen3[i,j])
    matricet3[i,j]=log(matricet3[i,j])
  }
}
#
### Analyse statistique
#

moyenne1= c()
moyenne2= c()
pvalue2= c()


for (i in 1:ncol(matricen3)){
  stat=t.test(matricen3[1:5,i],matricen3[6:10,i])
  moyenne1=c(moyenne1,as.numeric(stat$estimate[1]))
  moyenne2=c(moyenne2,as.numeric(stat$estimate[2]))
  pvalue2= c(pvalue2,stat$p.value)
}

padj2=p.adjust(pvalue2, method = "BH")
statnn3 = cbind(moyenne1,moyenne2,pvalue2,padj2)
rownames(statnn3)=c(seq1$V1)
statnn3 <- statnn3[order(statnn3[,4]),] # tri en fonction de pvalue adjuster

moyenne1= c()
moyenne2= c()
pvalue2= c()


for (i in 1:ncol(matricet3)){
  stat=t.test(matricet3[1:5,i],matricet3[6:10,i])
  moyenne1=c(moyenne1,as.numeric(stat$estimate[1]))
  moyenne2=c(moyenne2,as.numeric(stat$estimate[2]))
  pvalue2= c(pvalue2,stat$p.value)
}

padj2=p.adjust(pvalue2, method = "BH")
stattt3 = cbind(moyenne1,moyenne2,pvalue2,padj2)
rownames(stattt3)=c(seq1$V1)
stattt3 <- stattt3[order(stattt3[,4]),] # tri en fonction de pvalue adjuster



meany= c()
meanx= c()
pvalue= c()

for (i in 1:ncol(matricen3)){
  stat=t.test(matricen3[,i],matricet3[,i])
  meanx=c(meanx,as.numeric(stat$estimate[1]))
  meany=c(meany,as.numeric(stat$estimate[2]))
  pvalue= c(pvalue,stat$p.value)
}

padj=p.adjust(pvalue, method = "BH")
mstat3 = cbind(meanx,meany,pvalue,padj)
rownames(mstat3)=c(seq1$V1)
test3 =mstat3
test3 = cbind(test3,t(echant))
test3 <- test3[order(test3[,4]),]
mstat3 <- mstat3[order(mstat3[,4]),] # tri en fonction de pvalue adjuster

###tri en fonction des pvalue adjusté < 0.05
i=1
while (test3[i,4]< 0.05){
  j=i
  i=i+1
}
heatmap(test3[1:j,5:24])
#
##Top 3 kmer 
#

  #
  ##### Analyse 4_mer
  ####Normal echantillons k=4
  
  seq41=read.table("normal/4mer/TCG-38-4625-.fa")
  seq42=read.table("normal/4mer/TCG-38-4626-.fa")
  seq43=read.table("normal/4mer/TCG-38-4627-.fa")
  seq44=read.table("normal/4mer/TCG-38-4632-.fa")
  seq45=read.table("normal/4mer/TCG-44-2655-.fa")
  seq46=read.table("normal/4mer/TCG-44-2657-.fa")
  seq47=read.table("normal/4mer/TCG-44-266-.fa")
  seq48=read.table("normal/4mer/TCG-44-2662-.fa")
  seq49=read.table("normal/4mer/TCG-44-2665-.fa")
  seq410=read.table("normal/4mer/TCG-44-2668-.fa")
  
  seq41$V1= as.character(seq41$V1)
  
  matricen4=cbind(seq41$V2,seq42$V2,seq43$V2,seq44$V2,seq45$V2,seq46$V2,seq47$V2,seq48$V2,seq49$V2,seq410$V2)
  matricen4=t(matricen4)
  rownames(matricen4)=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10")
  colnames(matricen4)=c(seq41$V1)
  
  
  ##tumor echantillon k=4
  
  seqt41=read.table("tumor/4mer/TCG-05-4244-0.fa")
  seqt42=read.table("tumor/4mer/TCG-05-4249-0.fa")
  seqt43=read.table("tumor/4mer/TCG-05-4250-0.fa")
  seqt44=read.table("tumor/4mer/TCG-05-4382-0.fa")
  seqt45=read.table("tumor/4mer/TCG-05-4384-0.fa")
  seqt46=read.table("tumor/4mer/TCG-05-4389-0.fa")
  seqt47=read.table("tumor/4mer/TCG-05-4390-0.fa")
  seqt48=read.table("tumor/4mer/TCG-05-4395-0.fa")
  seqt49=read.table("tumor/4mer/TCG-05-4396-0.fa")
  seqt410=read.table("tumor/4mer/TCG-05-4397-0.fa")
  
  seqt41$V1= as.character(seqt41$V1)
  
  matricet4=cbind(seqt41$V2,seqt42$V2,seqt43$V2,seqt44$V2,seqt45$V2,seqt46$V2,seqt47$V2,seqt48$V2,seqt49$V2,seqt410$V2)
  matricet4=t(matricet4)
  rownames(matricet4)=c("seqt1","seqt2","seqt3","seqt4","seqt5","seqt6","seqt7","seqt8","seqt9","seqt10")
  colnames(matricet4)=c(seqt41$V1)
  
  
  
  
  #
  ##normalisation 
  #
  
  totaln=rowSums(matricen4)
  totalt=rowSums(matricet4)
  for (i in 1:nrow(matricen4)){
    for (j in 1:ncol(matricen4)){
      matricen4[i,j]=matricen4[i,j]/totaln[i]  
      matricet4[i,j]=matricet4[i,j]/totalt[i]
    }
  }
  
  freq1=matricen4[1,]
  freq2=matricen4[2,]
  freq3=matricen4[3,]
  freq4=matricen4[4,]
  freq5=matricen4[5,]
  freq6=matricen4[6,]
  freq7=matricen4[7,]
  freq8=matricen4[8,]
  freq9=matricen4[9,]
  freq10=matricen4[10,]
  freq1=freq1[order(-freq1)]
  freq2=freq2[order(-freq2)]
  freq3=freq3[order(-freq3)]
  freq4=freq4[order(-freq4)]
  freq5=freq5[order(-freq5)]
  freq6=freq6[order(-freq6)]
  freq7=freq7[order(-freq7)]
  freq8=freq8[order(-freq8)]
  freq9=freq9[order(-freq9)]
  freq10=freq10[order(-freq10)]
  
  protein4=read.table("normal/4mer/gencode/pc.fa")
  transcri4=read.table("normal/4mer/gencode/transcri.fa")
  names = protein4$V1
  protein4 = cbind(protein4$V2)
  rownames(protein4) = names
  tot = sum(protein4)
  for (i in 1:nrow(protein4)){
    protein4[i]=protein4[i]/tot
  }
  protein4=protein4[order(-protein4),]
  
  names = transcri4$V1
  transcri4 = cbind(transcri4$V2)
  rownames(transcri4) = names
  tot = sum(transcri4)
  for (i in 1:nrow(transcri4)){
    transcri4[i]=transcri4[i]/tot
  }
  transcri4=transcri4[order(-transcri4),]
  
  
  
  head(protein4)
  head(transcri4)
  head(freq1)
  head(freq2)
  head(freq3)
  head(freq4)
  head(freq5)
  head(freq6)
  head(freq7)
  head(freq8)
  head(freq9)
  head(freq10)
  
  #
  ##Contrôle fréquence
  #
  #head(seq41[order(-seq41$V2),])
  #head(seq41[order(seq41$V2),])
  #ne semble pas correct mais vu que c'est canonique résultat biaisé
  #
  ##passage au log
  #
  for (i in 1:nrow(matricen4)){
    for (j in 1:ncol(matricen4)){
      matricen4[i,j]=log(matricen4[i,j])
      matricet4[i,j]=log(matricet4[i,j])
    }
  }
  #
  ### Analyse statistique
  #
  
  moyenne1= c()
  moyenne2= c()
  pvalue2= c()
  
  
  for (i in 1:ncol(matricen4)){
    stat=t.test(matricen4[1:5,i],matricen4[6:10,i])
    moyenne1=c(moyenne1,as.numeric(stat$estimate[1]))
    moyenne2=c(moyenne2,as.numeric(stat$estimate[2]))
    pvalue2= c(pvalue2,stat$p.value)
  }
  
  padj2=p.adjust(pvalue2, method = "BH")
  statnn4 = cbind(moyenne1,moyenne2,pvalue2,padj2)
  rownames(statnn4)=c(seq41$V1)
  statnn4 <- statnn4[order(statnn4[,4]),] # tri en fonction de pvalue adjuster
  
  moyenne1= c()
  moyenne2= c()
  pvalue2= c()
  
  
  for (i in 1:ncol(matricet4)){
    stat=t.test(matricet4[1:5,i],matricet4[6:10,i])
    moyenne1=c(moyenne1,as.numeric(stat$estimate[1]))
    moyenne2=c(moyenne2,as.numeric(stat$estimate[2]))
    pvalue2= c(pvalue2,stat$p.value)
  }
  
  padj2=p.adjust(pvalue2, method = "BH")
  stattt4 = cbind(moyenne1,moyenne2,pvalue2,padj2)
  rownames(stattt4)=c(seq41$V1)
  stattt4 <- stattt4[order(stattt4[,4]),]
  
  
  
  meany= c()
  meanx= c()
  pvalue= c()
  
  for (i in 1:ncol(matricen4)){
    stat=t.test(matricen4[,i],matricet4[,i])
    meanx=c(meanx,as.numeric(stat$estimate[1]))
    meany=c(meany,as.numeric(stat$estimate[2]))
    pvalue= c(pvalue,stat$p.value)
  }
  
  padj=p.adjust(pvalue, method = "BH")
  mstat4= cbind(meanx,meany,pvalue,padj)
  rownames(mstat4)=c(seq41$V1)
  mstat4 <- mstat4[order(mstat4[,4]),] # tri en fonction de pvalue adjuster


#
##### Analyse 6_mer
#

###Normal echantillons k=6

seq61=read.table('normal/6mer/TCG-38-4625-.fa')
seq62=read.table('normal/6mer/TCG-38-4626-.fa')
seq63=read.table('normal/6mer/TCG-38-4627-.fa')
seq64=read.table('normal/6mer/TCG-38-4632-.fa')
seq65=read.table('normal/6mer/TCG-44-2655-.fa')
seq66=read.table('normal/6mer/TCG-44-2657-.fa')
seq67=read.table('normal/6mer/TCG-44-266-.fa')
seq68=read.table('normal/6mer/TCG-44-2662-.fa')
seq69=read.table('normal/6mer/TCG-44-2665-.fa')
seq610=read.table('normal/6mer/TCG-44-2668-.fa')

seq61$V1= as.character(seq61$V1)

matricen6=cbind(seq61$V2,seq62$V2,seq63$V2,seq64$V2,seq65$V2,seq66$V2,seq67$V2,seq68$V2,seq69$V2,seq610$V2)
matricen6= t(matricen6)
rownames(matricen6)=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10")
colnames(matricen6)=c(seq61$V1)

##tumor echantillon k=6

seqt61=read.table("tumor/6mer/TCG-05-4244-0.fa")
seqt62=read.table("tumor/6mer/TCG-05-4249-0.fa")
seqt63=read.table("tumor/6mer/TCG-05-4250-0.fa")
seqt64=read.table("tumor/6mer/TCG-05-4382-0.fa")
seqt65=read.table("tumor/6mer/TCG-05-4384-0.fa")
seqt66=read.table("tumor/6mer/TCG-05-4389-0.fa")
seqt67=read.table("tumor/6mer/TCG-05-4390-0.fa")
seqt68=read.table("tumor/6mer/TCG-05-4395-0.fa")
seqt69=read.table("tumor/6mer/TCG-05-4396-0.fa")
seqt610=read.table("tumor/6mer/TCG-05-4397-0.fa")

seqt61$V1= as.character(seqt61$V1)

matricet6=cbind(seqt61$V2,seqt62$V2,seqt63$V2,seqt64$V2,seqt65$V2,seqt66$V2,seqt67$V2,seqt68$V2,seqt69$V2,seqt610$V2)
matricet6 = t(matricet6)
rownames(matricet6)=c("seqt1","seqt2","seqt3","seqt4","seqt5","seqt6","seqt7","seqt8","seqt9","seqt10")
colnames(matricet6)=c(seqt61$V1)
#
##Normalisation
#
totaln=rowSums(matricen6)
totalt=rowSums(matricet6)
for (i in 1:nrow(matricen6)){
  for (j in 1:ncol(matricen6)){
    matricen6[i,j]=matricen6[i,j]/totaln[i]  
    matricet6[i,j]=matricet6[i,j]/totalt[i]
  }
}

freq1=matricen6[1,]
freq2=matricen6[2,]
freq3=matricen6[3,]
freq4=matricen6[4,]
freq5=matricen6[5,]
freq6=matricen6[6,]
freq7=matricen6[7,]
freq8=matricen6[8,]
freq9=matricen6[9,]
freq10=matricen6[10,]
freq1=freq1[order(-freq1)]
freq2=freq2[order(-freq2)]
freq3=freq3[order(-freq3)]
freq4=freq4[order(-freq4)]
freq5=freq5[order(-freq5)]
freq6=freq6[order(-freq6)]
freq7=freq7[order(-freq7)]
freq8=freq8[order(-freq8)]
freq9=freq9[order(-freq9)]
freq10=freq10[order(-freq10)]

protein6=read.table("normal/6mer/gencode/protin_c.fa")
transcri6=read.table("normal/6mer/gencode/transcri.fa")
names = protein6$V1
protein6 = cbind(protein6$V2)
rownames(protein6) = names
tot = sum(protein6)
for (i in 1:nrow(protein6)){
  protein6[i]=protein6[i]/tot
}
protein6=protein6[order(-protein6),]

names = transcri6$V1
transcri6 = cbind(transcri6$V2)
rownames(transcri6) = names
tot = sum(transcri6)
for (i in 1:nrow(transcri6)){
  transcri6[i]=transcri6[i]/tot
}
transcri6=transcri6[order(-transcri6),]


write.table(freq10, file='test.tsv', quote=FALSE, sep='\t')
head(protein6)
head(transcri6)
head(freq1)
head(freq2)
head(freq3)
head(freq4)
head(freq5)
head(freq6)
head(freq7)
head(freq8)
head(freq9)
head(freq10)

#
##passage au log
#
for (i in 1:nrow(matricen6)){
  for (j in 1:ncol(matricen6)){
    matricen6[i,j]=log(matricen6[i,j])
    matricet6[i,j]=log(matricet6[i,j])
  }
}

#
##Analyse stat avec t.test
#

moyenne1= c()
moyenne2= c()
pvalue2= c()


for (i in 1:ncol(matricen6)){
  stat=t.test(matricen6[1:5,i],matricen6[6:10,i])
  moyenne1=c(moyenne1,as.numeric(stat$estimate[1]))
  moyenne2=c(moyenne2,as.numeric(stat$estimate[2]))
  pvalue2= c(pvalue2,stat$p.value)
}

padj2=p.adjust(pvalue2, method = "BH")
statnn6 = cbind(moyenne1,moyenne2,pvalue2,padj2)
rownames(statnn6)=c(seq61$V1)
statnn6 <- statnn6[order(statnn6[,4]),] # tri en fonction de pvalue adjuster

moyenne1= c()
moyenne2= c()
pvalue2= c()


for (i in 1:ncol(matricet6)){
  stat=t.test(matricet6[1:5,i],matricet6[6:10,i])
  moyenne1=c(moyenne1,as.numeric(stat$estimate[1]))
  moyenne2=c(moyenne2,as.numeric(stat$estimate[2]))
  pvalue2= c(pvalue2,stat$p.value)
}

padj2=p.adjust(pvalue2, method = "BH")
stattt6 = cbind(moyenne1,moyenne2,pvalue2,padj2)
rownames(stattt6)=c(seq61$V1)
stattt6 <- stattt6[order(stattt6[,4]),]

meany= c()
meanx= c()
pvalue= c()

for (i in 1:ncol(matricen6)){
  stat=t.test(matricen6[,i],matricet6[,i])
  meanx=c(meanx,as.numeric(stat$estimate[1]))
  meany=c(meany,as.numeric(stat$estimate[2]))
  pvalue= c(pvalue,stat$p.value)
}

padj=p.adjust(pvalue, method = "BH")
mstat6 = cbind(meanx,meany,pvalue,padj)
rownames(mstat6)=c(seq61$V1)
mstat6 <- mstat6[order(mstat6[,4]),] # tri en fonction de pvalue adjuster


top6mer = mstat6[1:50,]
top4mer = mstat4[1:44,]
top3mer = mstat3[1:13,]

heatmap(t(top3mer))
heatmap(t(top4mer))
heatmap(t(top6mer))
?heatmap




protein4=read.table("normal/4mer/gencode/pc.fa")
transcri4=read.table("normal/4mer/gencode/transcri.fa")
head(protein4[order(-protein4$V2),])
head(seq41[order(-seq41$V2),])
head(transcri4[order(-transcri4$V2),])
head(seq42[order(-seq42$V2),])


protein6=read.table("normal/6mer/gencode/protin_c.fa")
transcri6=read.table("normal/6mer/gencode/transcri.fa")
head(protein6[order(-protein6$V2),])
head(seq61[order(-seq61$V2),])
head(transcri6[order(-transcri6$V2),])
head(seq62[order(-seq62$V2),])

head(protein3)


column_ha = HeatmapAnnotation(foo1 = runif(32), bar1 = anno_barplot(runif(32)))
row_ha = rowAnnotation(foo2 = runif(32), bar2 = anno_barplot(runif(32)))
Heatmap(mstat3, name = "noraml/tumeur", top_annotation = column_ha, right_annotation = row_ha)



