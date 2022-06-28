#install.packages("BiocManager")
library(RColorBrewer)
library(ComplexHeatmap)

seq1=read.table("3mer/normal/normal_1.fa")
seq2=read.table("3mer/normal/normal_10.fa")
seq3=read.table("3mer/normal/normal_11.fa")
seq4=read.table("3mer/normal/normal_12.fa")
seq5=read.table("3mer/normal/normal_14.fa")
seq6=read.table("3mer/normal/normal_16.fa")
seq7=read.table("3mer/normal/normal_17.fa")
seq8=read.table("3mer/normal/normal_18.fa")
seq9=read.table("3mer/normal/normal_19.fa")
seq10=read.table("3mer/normal/normal_20.fa")
seq11=read.table("3mer/normal/normal_21.fa")
seq12=read.table("3mer/normal/normal_22.fa")
seq13=read.table("3mer/normal/normal_25.fa")
seq14=read.table("3mer/normal/normal_26.fa")
seq15=read.table("3mer/normal/normal_27.fa")
seq16=read.table("3mer/normal/normal_32.fa")
seq17=read.table("3mer/normal/normal_33.fa")
seq18=read.table("3mer/normal/normal_34.fa")
seq19=read.table("3mer/normal/normal_35.fa")
seq20=read.table("3mer/normal/normal_36.fa")
seq21=read.table("3mer/normal/normal_37.fa")
seq22=read.table("3mer/normal/normal_38.fa")
seq23=read.table("3mer/normal/normal_39.fa")
seq24=read.table("3mer/normal/normal_40.fa")
seq25=read.table("3mer/normal/normal_41.fa")
seq26=read.table("3mer/normal/normal_42.fa")
seq27=read.table("3mer/normal/normal_43.fa")
seq28=read.table("3mer/normal/normal_44.fa")
seq29=read.table("3mer/normal/normal_45.fa")
seq30=read.table("3mer/normal/normal_46.fa")

seq1$V1= as.character(seq1$V1)

matricen3=cbind(seq1$V2,seq2$V2,seq3$V2,seq4$V2,seq5$V2,seq6$V2,seq7$V2,seq8$V2,seq9$V2,seq10$V2,seq11$V2,
                seq12$V2,seq13$V2,seq14$V2,seq15$V2,seq16$V2,seq17$V2,seq18$V2,seq19$V2,seq20$V2,
                seq21$V2,seq22$V2,seq23$V2,seq24$V2,seq25$V2,seq26$V2,seq27$V2,seq28$V2,seq29$V2,seq30$V2)
matricen3 = t(matricen3)
rownames(matricen3)=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10","seq11","seq12","seq13","seq14","seq15","seq16","seq17","seq18","seq19","seq20","seq21","seq22","seq23","seq24","seq25","seq26","seq27","seq28","seq29","seq30")
colnames(matricen3)=c(seq1$V1)


##tumor echantillon k=3

seqt1=read.table("3mer/tumor/tumor_1.fa")
seqt2=read.table("3mer/tumor/tumor_10.fa")
seqt3=read.table("3mer/tumor/tumor_11.fa")
seqt4=read.table("3mer/tumor/tumor_14.fa")
seqt5=read.table("3mer/tumor/tumor_16.fa")
seqt6=read.table("3mer/tumor/tumor_17.fa")
seqt7=read.table("3mer/tumor/tumor_18.fa")
seqt8=read.table("3mer/tumor/tumor_19.fa")
seqt9=read.table("3mer/tumor/tumor_20.fa")
seqt10=read.table("3mer/tumor/tumor_21.fa")
seqt11=read.table("3mer/tumor/tumor_23.fa")
seqt12=read.table("3mer/tumor/tumor_24.fa")
seqt13=read.table("3mer/tumor/tumor_25.fa")
seqt14=read.table("3mer/tumor/tumor_26.fa")
seqt15=read.table("3mer/tumor/tumor_27.fa")
seqt16=read.table("3mer/tumor/tumor_46.fa")
seqt17=read.table("3mer/tumor/tumor_47.fa")
seqt18=read.table("3mer/tumor/tumor_48.fa")
seqt19=read.table("3mer/tumor/tumor_67.fa")
seqt20=read.table("3mer/tumor/tumor_68.fa")
seqt21=read.table("3mer/tumor/tumor_69.fa")
seqt22=read.table("3mer/tumor/tumor_70.fa")
seqt23=read.table("3mer/tumor/tumor_71.fa")
seqt24=read.table("3mer/tumor/tumor_72.fa")
seqt25=read.table("3mer/tumor/tumor_73.fa")
seqt26=read.table("3mer/tumor/tumor_74.fa")
seqt27=read.table("3mer/tumor/tumor_75.fa")
seqt28=read.table("3mer/tumor/tumor_76.fa")
seqt29=read.table("3mer/tumor/tumor_77.fa")
seqt30=read.table("3mer/tumor/tumor_78.fa")

seqt1$V1= as.character(seqt1$V1)

matricet3=cbind(seqt1$V2,seqt2$V2,seqt3$V2,seqt4$V2,seqt5$V2,seqt6$V2,seqt7$V2,seqt8$V2,seqt9$V2,seqt10$V2,
                seqt11$V2,seqt12$V2,seqt13$V2,seqt14$V2,seqt15$V2,seqt16$V2,seqt17$V2,seqt18$V2,seqt19$V2,seqt20$V2,
                seqt21$V2,seqt22$V2,seqt23$V2,seqt24$V2,seqt25$V2,seqt26$V2,seqt27$V2,seqt28$V2,seqt29$V2,seqt30$V2)
matricet3 = t(matricet3)
rownames(matricet3)=c("seqt1","seqt2","seqt3","seqt4","seqt5","seqt6","seqt7","seqt8","seqt9","seqt10",
                      "seqt11","seqt12","seqt13","seqt14","seqt15","seqt16","seqt17","seqt18","seqt19","seqt20",
                      "seqt21","seqt22","seqt23","seqt24","seqt25","seqt26","seqt27","seqt28","seqt29","seqt30")
colnames(matricet3)=c(seqt1$V1)

totaln=rowSums(matricen3)
totalt=rowSums(matricet3)
for (i in 1:nrow(matricen3)){
  for (j in 1:ncol(matricen3)){
    matricen3[i,j]=matricen3[i,j]/totaln[i]  
    matricet3[i,j]=matricet3[i,j]/totalt[i]
  }
}
echant = rbind(matricen3,matricet3)



rownames(echant)=c("normal_1","normal_10","normal_11","normal_12","normal_14","normal_16","normal_17","normal_18","normal_19","normal_20","normal_21","normal_22","normal_25","normal_26","normal_27","normal_32","normal_33","normal_34","normal_35","normal_36","normal_37","normal_38","normal_39","normal_40","normal_41","normal_42","normal_43","normal_44","normal_45","normal_46",
                   "tumor_1","tumor_10","tumor_11","tumor_14","tumor_16","tumor_17","tumor_18","tumor_19","tumor_20","tumor_21","tumor_23","tumor_24","tumor_25","tumor_26","tumor_27","tumor_46","tumor_47","tumor_48","tumor_67","tumor_68","tumor_69","tumor_70","tumor_71","tumor_72","tumor_73","tumor_74","tumor_75","tumor_76","tumor_77","tumor_78")

for (i in 1:nrow(matricen3)){
  for (j in 1:ncol(matricen3)){
    matricen3[i,j]=log(matricen3[i,j])
    matricet3[i,j]=log(matricet3[i,j])
  }
}
meanx=c()
meany=c()
pvalue=c()
padj=c()
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
while (test3[i,4]< 0.5){
  j=i
  i=i+1
}
j
my_group=cbind(c("normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal",
                 "tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor"))
colnames(my_group)=c("Type")
my_group=data.frame(my_group)

colors <- c("red","green")
names(colors)<-c("tumor","normal")
colors<-list('Type' = c('tumor' = 'cyan', 'normal' = 'red'))

colAnn <- HeatmapAnnotation(df = my_group ,col = colors,which = 'col')

Heatmap(t(scale(t(test3[1:j,5:64]))),column_title ="Echantillon Normal/Tumor",row_title ="3_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25),column_labels = colnames(test3[,5:64]),row_km = 2,column_km = 2,top_annotation = colAnn)


#
##### Analyse 4_mer
####Normal echantillons k=4

seq1=read.table("4mer/normal/normal_1.fa")
seq2=read.table("4mer/normal/normal_10.fa")
seq3=read.table("4mer/normal/normal_11.fa")
seq4=read.table("4mer/normal/normal_12.fa")
seq5=read.table("4mer/normal/normal_14.fa")
seq6=read.table("4mer/normal/normal_16.fa")
seq7=read.table("4mer/normal/normal_17.fa")
seq8=read.table("4mer/normal/normal_18.fa")
seq9=read.table("4mer/normal/normal_19.fa")
seq10=read.table("4mer/normal/normal_20.fa")
seq11=read.table("4mer/normal/normal_21.fa")
seq12=read.table("4mer/normal/normal_22.fa")
seq13=read.table("4mer/normal/normal_25.fa")
seq14=read.table("4mer/normal/normal_26.fa")
seq15=read.table("4mer/normal/normal_27.fa")
seq16=read.table("4mer/normal/normal_32.fa")
seq17=read.table("4mer/normal/normal_33.fa")
seq18=read.table("4mer/normal/normal_34.fa")
seq19=read.table("4mer/normal/normal_35.fa")
seq20=read.table("4mer/normal/normal_36.fa")
seq21=read.table("4mer/normal/normal_37.fa")
seq22=read.table("4mer/normal/normal_38.fa")
seq23=read.table("4mer/normal/normal_39.fa")
seq24=read.table("4mer/normal/normal_40.fa")
seq25=read.table("4mer/normal/normal_41.fa")
seq26=read.table("4mer/normal/normal_42.fa")
seq27=read.table("4mer/normal/normal_43.fa")
seq28=read.table("4mer/normal/normal_44.fa")
seq29=read.table("4mer/normal/normal_45.fa")
seq30=read.table("4mer/normal/normal_46.fa")

seq1$V1= as.character(seq1$V1)

matricen4=cbind(seq1$V2,seq2$V2,seq3$V2,seq4$V2,seq5$V2,seq6$V2,seq7$V2,seq8$V2,seq9$V2,seq10$V2,seq11$V2,
                seq12$V2,seq13$V2,seq14$V2,seq15$V2,seq16$V2,seq17$V2,seq18$V2,seq19$V2,seq20$V2,
                seq21$V2,seq22$V2,seq23$V2,seq24$V2,seq25$V2,seq26$V2,seq27$V2,seq28$V2,seq29$V2,seq30$V2)
matricen4 = t(matricen4)
rownames(matricen4)=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10","seq11","seq12","seq13","seq14","seq15","seq16","seq17","seq18","seq19","seq20","seq21","seq22","seq23","seq24","seq25","seq26","seq27","seq28","seq29","seq30")
colnames(matricen4)=c(seq1$V1)


##tumor echantillon k=4

seqt1=read.table("4mer/tumor/tumor_1.fa")
seqt2=read.table("4mer/tumor/tumor_10.fa")
seqt3=read.table("4mer/tumor/tumor_11.fa")
seqt4=read.table("4mer/tumor/tumor_14.fa")
seqt5=read.table("4mer/tumor/tumor_16.fa")
seqt6=read.table("4mer/tumor/tumor_17.fa")
seqt7=read.table("4mer/tumor/tumor_18.fa")
seqt8=read.table("4mer/tumor/tumor_19.fa")
seqt9=read.table("4mer/tumor/tumor_20.fa")
seqt10=read.table("4mer/tumor/tumor_21.fa")
seqt11=read.table("4mer/tumor/tumor_23.fa")
seqt12=read.table("4mer/tumor/tumor_24.fa")
seqt13=read.table("4mer/tumor/tumor_25.fa")
seqt14=read.table("4mer/tumor/tumor_26.fa")
seqt15=read.table("4mer/tumor/tumor_27.fa")
seqt16=read.table("4mer/tumor/tumor_46.fa")
seqt17=read.table("4mer/tumor/tumor_47.fa")
seqt18=read.table("4mer/tumor/tumor_48.fa")
seqt19=read.table("4mer/tumor/tumor_67.fa")
seqt20=read.table("4mer/tumor/tumor_68.fa")
seqt21=read.table("4mer/tumor/tumor_69.fa")
seqt22=read.table("4mer/tumor/tumor_70.fa")
seqt23=read.table("4mer/tumor/tumor_71.fa")
seqt24=read.table("4mer/tumor/tumor_72.fa")
seqt25=read.table("4mer/tumor/tumor_73.fa")
seqt26=read.table("4mer/tumor/tumor_74.fa")
seqt27=read.table("4mer/tumor/tumor_75.fa")
seqt28=read.table("4mer/tumor/tumor_76.fa")
seqt29=read.table("4mer/tumor/tumor_77.fa")
seqt30=read.table("4mer/tumor/tumor_78.fa")

seqt1$V1= as.character(seqt1$V1)

matricet4=cbind(seqt1$V2,seqt2$V2,seqt3$V2,seqt4$V2,seqt5$V2,seqt6$V2,seqt7$V2,seqt8$V2,seqt9$V2,seqt10$V2,
                seqt11$V2,seqt12$V2,seqt13$V2,seqt14$V2,seqt15$V2,seqt16$V2,seqt17$V2,seqt18$V2,seqt19$V2,seqt20$V2,
                seqt21$V2,seqt22$V2,seqt23$V2,seqt24$V2,seqt25$V2,seqt26$V2,seqt27$V2,seqt28$V2,seqt29$V2,seqt30$V2)
matricet4 = t(matricet4)
rownames(matricet4)=c("seqt1","seqt2","seqt3","seqt4","seqt5","seqt6","seqt7","seqt8","seqt9","seqt10",
                      "seqt11","seqt12","seqt13","seqt14","seqt15","seqt16","seqt17","seqt18","seqt19","seqt20",
                      "seqt21","seqt22","seqt23","seqt24","seqt25","seqt26","seqt27","seqt28","seqt29","seqt30")
colnames(matricet4)=c(seqt1$V1)
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

echant4 = rbind(matricen4,matricet4)
rownames(echant4)=c("normal_1","normal_10","normal_11","normal_12","normal_14","normal_16","normal_17","normal_18","normal_19","normal_20","normal_21","normal_22","normal_25","normal_26","normal_27","normal_32","normal_33","normal_34","normal_35","normal_36","normal_37","normal_38","normal_39","normal_40","normal_41","normal_42","normal_43","normal_44","normal_45","normal_46",
                   "tumor_1","tumor_10","tumor_11","tumor_14","tumor_16","tumor_17","tumor_18","tumor_19","tumor_20","tumor_21","tumor_23","tumor_24","tumor_25","tumor_26","tumor_27","tumor_46","tumor_47","tumor_48","tumor_67","tumor_68","tumor_69","tumor_70","tumor_71","tumor_72","tumor_73","tumor_74","tumor_75","tumor_76","tumor_77","tumor_78")

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
rownames(mstat4)=c(seq1$V1)
test4 =mstat4
test4 = cbind(test4,t(echant4))
test4 <- test4[order(test4[,4]),]
mstat4 <- mstat4[order(mstat4[,4]),] # tri en fonction de pvalue adjuster

###tri en fonction des pvalue adjusté < 0.05
i=1
while (test4[i,4]< 0.05){
  j=i
  i=i+1
}
j
#my_group=cbind(c("normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor"))
my_group=cbind(c("normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal",
                 "tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor"))

colnames(my_group)=c("Type")
my_group=data.frame(my_group)

colors <- c("red","green")
names(colors)<-c("tumor","normal")
colors<-list('Type' = c('tumor' = 'cyan', 'normal' = 'red'))

colAnn <- HeatmapAnnotation(df = my_group ,col = colors,which = 'col')

Heatmap(t(scale(t(test4[1:j,5:64]))),column_title ="Echantillon Normal/Tumor",row_title ="4_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25),column_labels = colnames(test4[,5:64]),row_km = 2,column_km = 2,top_annotation = colAnn)


#Heatmap(t(scale(t(test4[1:j,5:64]))),column_title ="Echantillon Normal/Tumor",row_title ="4_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25),column_labels = colnames(test4[,5:64]),row_km = 2,column_km = 2,top_annotation = colAnn)

#Heatmap(t(scale(test4[1:j,5:64])))

#
##### Analyse 6_mer
####Normal echantillons k=6

seq1=read.table("6mer/normal/normal_1.fa")
seq2=read.table("6mer/normal/normal_10.fa")
seq3=read.table("6mer/normal/normal_11.fa")
seq4=read.table("6mer/normal/normal_12.fa")
seq5=read.table("6mer/normal/normal_14.fa")
seq6=read.table("6mer/normal/normal_16.fa")
seq7=read.table("6mer/normal/normal_17.fa")
seq8=read.table("6mer/normal/normal_18.fa")
seq9=read.table("6mer/normal/normal_19.fa")
seq10=read.table("6mer/normal/normal_20.fa")
seq11=read.table("6mer/normal/normal_21.fa")
seq12=read.table("6mer/normal/normal_22.fa")
seq13=read.table("6mer/normal/normal_25.fa")
seq14=read.table("6mer/normal/normal_26.fa")
seq15=read.table("6mer/normal/normal_27.fa")
seq16=read.table("6mer/normal/normal_32.fa")
seq17=read.table("6mer/normal/normal_33.fa")
seq18=read.table("6mer/normal/normal_34.fa")
seq19=read.table("6mer/normal/normal_35.fa")
seq20=read.table("6mer/normal/normal_36.fa")
seq21=read.table("6mer/normal/normal_37.fa")
seq22=read.table("6mer/normal/normal_38.fa")
seq23=read.table("6mer/normal/normal_39.fa")
seq24=read.table("6mer/normal/normal_40.fa")
seq25=read.table("6mer/normal/normal_41.fa")
seq26=read.table("6mer/normal/normal_42.fa")
seq27=read.table("6mer/normal/normal_43.fa")
seq28=read.table("6mer/normal/normal_44.fa")
seq29=read.table("6mer/normal/normal_45.fa")
seq30=read.table("6mer/normal/normal_46.fa")

seq1$V1= as.character(seq1$V1)

matricen6=cbind(seq1$V2,seq2$V2,seq3$V2,seq4$V2,seq5$V2,seq6$V2,seq7$V2,seq8$V2,seq9$V2,seq10$V2,seq11$V2,
                seq12$V2,seq13$V2,seq14$V2,seq15$V2,seq16$V2,seq17$V2,seq18$V2,seq19$V2,seq20$V2,
                seq21$V2,seq22$V2,seq23$V2,seq24$V2,seq25$V2,seq26$V2,seq27$V2,seq28$V2,seq29$V2,seq30$V2)
matricen6 = t(matricen6)
rownames(matricen6)=c("seq1","seq2","seq3","seq4","seq5","seq6","seq7","seq8","seq9","seq10","seq11","seq12","seq13","seq14","seq15","seq16","seq17","seq18","seq19","seq20","seq21","seq22","seq23","seq24","seq25","seq26","seq27","seq28","seq29","seq30")
colnames(matricen6)=c(seq1$V1)


##tumor echantillon k=6

seqt1=read.table("6mer/tumor/tumor_1.fa")
seqt2=read.table("6mer/tumor/tumor_10.fa")
seqt3=read.table("6mer/tumor/tumor_11.fa")
seqt4=read.table("6mer/tumor/tumor_14.fa")
seqt5=read.table("6mer/tumor/tumor_16.fa")
seqt6=read.table("6mer/tumor/tumor_17.fa")
seqt7=read.table("6mer/tumor/tumor_18.fa")
seqt8=read.table("6mer/tumor/tumor_19.fa")
seqt9=read.table("6mer/tumor/tumor_20.fa")
seqt10=read.table("6mer/tumor/tumor_21.fa")
seqt11=read.table("6mer/tumor/tumor_23.fa")
seqt12=read.table("6mer/tumor/tumor_24.fa")
seqt13=read.table("6mer/tumor/tumor_25.fa")
seqt14=read.table("6mer/tumor/tumor_26.fa")
seqt15=read.table("6mer/tumor/tumor_27.fa")
seqt16=read.table("6mer/tumor/tumor_46.fa")
seqt17=read.table("6mer/tumor/tumor_47.fa")
seqt18=read.table("6mer/tumor/tumor_48.fa")
seqt19=read.table("6mer/tumor/tumor_67.fa")
seqt20=read.table("6mer/tumor/tumor_68.fa")
seqt21=read.table("6mer/tumor/tumor_69.fa")
seqt22=read.table("6mer/tumor/tumor_70.fa")
seqt23=read.table("6mer/tumor/tumor_71.fa")
seqt24=read.table("6mer/tumor/tumor_72.fa")
seqt25=read.table("6mer/tumor/tumor_73.fa")
seqt26=read.table("6mer/tumor/tumor_74.fa")
seqt27=read.table("6mer/tumor/tumor_75.fa")
seqt28=read.table("6mer/tumor/tumor_76.fa")
seqt29=read.table("6mer/tumor/tumor_77.fa")
seqt30=read.table("6mer/tumor/tumor_78.fa")

seqt1$V1= as.character(seqt1$V1)

matricet6=cbind(seqt1$V2,seqt2$V2,seqt3$V2,seqt4$V2,seqt5$V2,seqt6$V2,seqt7$V2,seqt8$V2,seqt9$V2,seqt10$V2,
                seqt11$V2,seqt12$V2,seqt13$V2,seqt14$V2,seqt15$V2,seqt16$V2,seqt17$V2,seqt18$V2,seqt19$V2,seqt20$V2,
                seqt21$V2,seqt22$V2,seqt23$V2,seqt24$V2,seqt25$V2,seqt26$V2,seqt27$V2,seqt28$V2,seqt29$V2,seqt30$V2)
matricet6 = t(matricet6)
rownames(matricet6)=c("seqt1","seqt2","seqt3","seqt4","seqt5","seqt6","seqt7","seqt8","seqt9","seqt10",
                      "seqt11","seqt12","seqt13","seqt14","seqt15","seqt16","seqt17","seqt18","seqt19","seqt20",
                      "seqt21","seqt22","seqt23","seqt24","seqt25","seqt26","seqt27","seqt28","seqt29","seqt30")
colnames(matricet6)=c(seqt1$V1)

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


echant6 = rbind(matricen6,matricet6)
#rownames(echant6)= c("normal1","normal2","normal3","normal4","normal5","normal6","normal7","normal8","normal9","normal10","tumor1","tumor2","tumor3","tumor4","tumor5","tumor6","tumor7","tumor8","tumor9","tumor10")
rownames(echant6)=c("normal_1","normal_10","normal_11","normal_12","normal_14","normal_16","normal_17","normal_18","normal_19","normal_20","normal_21","normal_22","normal_25","normal_26","normal_27","normal_32","normal_33","normal_34","normal_35","normal_36","normal_37","normal_38","normal_39","normal_40","normal_41","normal_42","normal_43","normal_44","normal_45","normal_46",
                    "tumor_1","tumor_10","tumor_11","tumor_14","tumor_16","tumor_17","tumor_18","tumor_19","tumor_20","tumor_21","tumor_23","tumor_24","tumor_25","tumor_26","tumor_27","tumor_46","tumor_47","tumor_48","tumor_67","tumor_68","tumor_69","tumor_70","tumor_71","tumor_72","tumor_73","tumor_74","tumor_75","tumor_76","tumor_77","tumor_78")


#
##passage au log
#
for (i in 1:nrow(matricen6)){
  for (j in 1:ncol(matricen6)){
    matricen6[i,j]=log(matricen6[i,j])
    matricet6[i,j]=log(matricet6[i,j])
  }
}

meany= c()
meanx= c()
pvalue= c()

for (i in 1:ncol(matricen6)){
  stat6=t.test(matricen6[,i],matricet6[,i])
  meanx=c(meanx,as.numeric(stat6$estimate[1]))
  meany=c(meany,as.numeric(stat6$estimate[2]))
  pvalue= c(pvalue,stat6$p.value)
}

padj=p.adjust(pvalue, method = "BH")
mstat6 = cbind(meanx,meany,pvalue,padj)
rownames(mstat6)=seq6$V1

test6 =mstat6
test6 = cbind(test6,t(echant6))
test6 <- test6[order(test6[,4]),]

mstat6 <- mstat6[order(mstat6[,4]),]



##figure

i=1
while (test6[i,4]< 0.00005){
  j=i
  i=i+1
}
j=2080

#my_group=cbind(c("normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor"))
my_group=cbind(c("normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal",
                 "tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor","tumor"))

colnames(my_group)=c("Type")
my_group=data.frame(my_group)

colors <- c("red","green")
names(colors)<-c("tumor","normal")
colors<-list('Type' = c('tumor' = 'cyan', 'normal' = 'red'))

colAnn <- HeatmapAnnotation(df = my_group ,col = colors,which = 'col')

#pour alpha = 0.001
j=236
Heatmap(t(scale(t(test6[1:j,5:64]))),column_title ="Echantillon Normal/Tumor",row_title ="6_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25),column_labels = colnames(test6[,5:64]),row_km = 2,column_km = 2,top_annotation = colAnn)
#Heatmap(test6[1:j,5:24],column_title ="Echantillon Normal/Tumor",row_title ="6_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = terrain.colors(256),column_labels = colnames(test6[,5:24]),row_km = 4,column_km = 2)
#Heatmap(test6[1:j,5:24])
#pour alpha = 0.0001
j=75
Heatmap(t(scale(t(test6[1:j,5:24]))),column_title ="Echantillon Normal/Tumor",row_title ="6_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25),column_labels = colnames(test6[,5:24]),row_km = 2,column_km = 2,top_annotation = colAnn)
#Heatmap(test6[1:j,5:24],column_title ="Echantillon Normal/Tumor",row_title ="6_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = terrain.colors(256),column_labels = colnames(test6[,5:24]),row_km = 3,column_km = 2)
#heatmap(test6[1:j,5:24])
#pour alpha = 0.00005
j=44
#col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))


Heatmap(t(scale(t(test6[1:j,5:64]))),column_title ="Echantillon Normal/Tumor",row_title ="6_mer",name="frequence k_mer",cluster_rows = TRUE,cluster_columns = TRUE, cluster_row_slices = TRUE,cluster_column_slices =TRUE,  col = colorRampPalette(brewer.pal(8, "RdYlBu"))(25),column_labels = colnames(test6[,5:64]),row_km = 2,column_km = 2,top_annotation = colAnn)

library(questionr)

test=shapiro.test(matricen6[1,])
test$p.value
p_val = c()
for (i in 1:nrow(matricet6)){
  test=shapiro.test(matricet6[i,])
  p_val=c(p_val,test$p.value)
}
p_val
shapiro.test(matricen6[1,])



