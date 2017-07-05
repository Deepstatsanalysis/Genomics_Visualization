# This is a code for Transcript Visualization. The input is Gene Symbol and 
# The output is a plot which shows Transcripts, Exons, Introns,
# Coding sites and other informations

  annot_plot<-function(Gene_symb){
  
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(ggplot2)
  require(org.Hs.eg.db)
  
  Gene_name<-Gene_symb
  txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  ids<-select(org.Hs.eg.db, keys=Gene_name, keytype = "SYMBOL", columns = "ENTREZID")
  Q<-select(txdb, keys=ids$ENTREZID, keytype = "GENEID",
            columns = c("TXCHROM" ,"TXID","TXSTRAND" ,"TXSTART","TXEND", "EXONID","EXONSTART","EXONEND"))
  Q1<-select(txdb, keys=ids$ENTREZID, keytype = "GENEID",columns = c("CDSSTART","CDSEND","TXID"))
    
  num_line<-unique(Q$TXID)
  num_level<-seq(1:length(num_line))
  re<-re1<-0
  for (i in 1:length(num_line)){
    re[i]<- sum(Q$TXID==num_line[i])
    re1[i]<-sum(Q1$TXID==num_line[i])
  }
  
  num_levels1<-rep(num_level, re)
  num_levels2<-rep(num_level, re1)
  group1<-rep(seq(1:nrow(Q)),2)
  
  num_desc1<-rep(num_levels1,2)
  group2<-rep(seq(1:nrow(Q1)),2)
  num_desc2<-rep(num_levels2,2)
  
  Ex=c(Q$EXONSTART,Q$EXONEND)
  Tr=c(Q$TXSTART,Q$TXEND)
  Tr1=c(Q1$TXSTART,Q1$TXEND)
  CD=c(Q1$CDSSTART,Q1$CDSEND)
  s<-data.frame(Ex,  Tr,  group1, num_desc1, num_Chr=rep(Q$TXCHROM,2),  Chr_strand=rep(Q$TXSTRAND,2))
  s1<-data.frame(CD, group2, num_desc2)
  
  ggplot(s)+geom_line(aes(x=Tr, y=num_desc1, group=group1))+
    geom_line(aes(x=Ex, y=num_desc1, group=group1),color="blue", size=12, alpha=.3)+
    geom_line(data=s1,aes(x=CD, y=num_desc2, group=group2),color="red", size=6, alpha=.5)+
    labs(x=Q$TXCHROM[1], y="Transcripts")+ 
    theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
    ylim(0, length(num_line)+1)+
    annotate("text", x = min(s$Tr)+(max(s$Tr)-min(s$Tr))/10, y = seq(length(num_line)+0.4,length(num_line)+1,len=3)
             , label = c(paste0("Gene :",Gene_name),paste0("Gene ID :",Q$GENEID[1]),paste0("Tx Strand: ",Q$TXSTRAND[1])),size=3.5)+
    annotate("text", x = max(s$Tr)-(max(s$Tr)-min(s$Tr))/50, y = num_level+0.25,label = paste0("TX_ID :",num_line),size=3.5)

  
  }
