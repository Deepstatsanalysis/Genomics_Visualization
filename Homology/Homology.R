#Problem 4 and Problem 5

# First:Solving Problem 4
require(BSgenome)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(ggplot2)
require(BSgenome.Hsapiens.UCSC.hg19)
require(BSgenome.Mmusculus.UCSC.mm10)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)

#Section 1
#Human Genome Length
HG<-seqinfo(Hsapiens)
HG1<-as.data.frame(HG)
HGL<-sum(as.numeric(HG1[,1]))
#Mouse Genome Length
MG<-seqinfo(Mmusculus)
MG1<-as.data.frame(MG)
MGL<-sum(as.numeric(MG1[,1]))

HAF<-matrix(0, length(Hsapiens), 5); GC_content_Hs<-0
for (i in 1:length(Hsapiens)){
  HAF[i,]<-as.matrix(alphabetFrequency(Hsapiens[[i]],   as.prob=T, baseOnly=T))
  GC_content_Hs[i]<-letterFrequency(Hsapiens[[i]],  letters="CG",as.prob=T)
}

MAF<-matrix(0, length(Mmusculus), 5); GC_content_Mm<-0
for (i in 1:length(Mmusculus)){
  MAF[i,]<-as.matrix(alphabetFrequency(Mmusculus[[i]],   as.prob=T, baseOnly=T))
  GC_content_Mm[i]<-letterFrequency(Mmusculus[[i]],  letters="CG",as.prob=T)
}
colnames(HAF)<-c("A","C","G","T","Others")
colnames(MAF)<-c("A","C","G","T","Others")

Hs_summary<-data.frame(chr_names = row.names(HG1),chr_length=HG1[,1], HAF, GC_content=GC_content_Hs)
Mm_summary<-data.frame(chr_names = row.names(MG1),chr_length=MG1[,1], MAF, GC_content=GC_content_Mm)

X<-Hs_summary[1:25,]
X$chr_name <- factor(X$chr_name, levels = X$chr_name)
ggplot(X, aes(x=chr_name,y=GC_content))+geom_bar(stat="identity",fill="white", colour="black")+
  xlab("Chromosomes Name")+ylab("GC Content")+ggtitle("GC Content of Human Chromosomes")

Y<-Mm_summary[1:22,]
Y$chr_name <- factor(Y$chr_name, levels = Y$chr_name)
ggplot(Y, aes(x=chr_name,y=GC_content))+geom_bar(stat="identity",fill="white", colour="black")+
  xlab("Chromosomes Name")+ylab("GC Content")+ggtitle("GC Content of Mouse Chromosomes")

dev.copy(pdf,"GC_Content_Human.pdf", width=8, height=6)
dev.off()

#GC Content of Human CHR4 and Mouse CHR3 CHR5 CHR8
GC_content_H4 <- letterFrequency(Hsapiens$chr4,  letters="CG",as.prob=T)
GC_content_M3 <- letterFrequency(Mmusculus$chr3, letters="CG",as.prob=T)
GC_content_M5 <- letterFrequency(Mmusculus$chr5, letters="CG",as.prob=T)
GC_content_M8 <- letterFrequency(Mmusculus$chr8, letters="CG",as.prob=T)

##..........................................................................................
#Section 2
Txdb_h<-TxDb.Hsapiens.UCSC.hg19.knownGene
Txdb_m<-TxDb.Mmusculus.UCSC.mm10.knownGene
#Chromosomes in Human Genome Database
rownames(HG1)
seqlevels(Txdb_h)
#Chromosomes in Mouse Genome Database
rownames(MG1)
seqlevels(Txdb_m)

seqlevels(Txdb_h)<- paste0("chr",c(1:22,"X","Y","M"))
seqlevels(Txdb_m)<- paste0("chr",c(1:19,"X","Y","M"))
##...........................................................................................
#Section 3
CDH<-cds(Txdb_h)
CDM<-cds(Txdb_m)

CDH_pos<-cds(Txdb_h, vals<-list(tx_strand="+"))
CDH_neg<-cds(Txdb_h, vals<-list(tx_strand="-"))

CDM_pos<-cds(Txdb_m, vals<-list(tx_strand="+"))
CDM_neg<-cds(Txdb_m, vals<-list(tx_strand="-"))

#Number of Human and Mouse Coding Regions
HCds<-length(CDH)
MCds<-length(CDM)

#Number of Overlaps in Human and Mouse Coding Regions for both Positive and Negative strands
#cdh_seqs1 <- getSeq(Hsapiens, Cd_H )
FH<-findOverlaps(CDH,CDH,ignore.strand=FALSE)
FH_pos<-findOverlaps(CDH_pos,CDH_pos)
FH_neg<-findOverlaps(CDH_neg,CDH_neg)

(length(FH_pos)-length(CDH_pos))/2
(length(FH_neg)-length(CDH_neg))/2
(length(FH)-length(CDH))/2


FM<-findOverlaps(CDM,CDM,ignore.strand=FALSE)
FM_pos<-findOverlaps(CDM_pos,CDM_pos)
FM_neg<-findOverlaps(CDM_neg,CDM_neg)

(length(FM_pos)-length(CDM_pos))/2
(length(FM_neg)-length(CDM_neg))/2
(length(FM)-length(CDM))/2

CDH_non <- reduce(CDH, with.revmap=TRUE, ignore.strand=FALSE)
revmap <- mcols(CDH_non)$revmap     # an IntegerList

CDM_non <- reduce(CDM, with.revmap=TRUE, ignore.strand=FALSE)
revmap1 <- mcols(CDM_non)$revmap    # an IntegerList

#Find Overlaps between Two Strands Without Considering Each Strand Seperately.
CDH_nonboth<-reduce(CDH, with.revmap=TRUE, ignore.strand=TRUE)
CDM_nonboth<-reduce(CDM, with.revmap=TRUE, ignore.strand=TRUE)

getSeq(Hsapiens,  CDH_non)
getSeq(Mmusculus, CDM_non)

##.......................................................................................
#Section 4
H4<-CDH_non[seqnames(CDH_non) == 'chr4']
M3<-CDM_non[seqnames(CDM_non) == 'chr3']
M5<-CDM_non[seqnames(CDM_non) == 'chr5']
M8<-CDM_non[seqnames(CDM_non) == 'chr8']

H4_seq<-getSeq(Hsapiens, H4)
M3_seq<-getSeq(Mmusculus, M3)
M5_seq<-getSeq(Mmusculus, M5)
M8_seq<-getSeq(Mmusculus, M8)

H4_seq1<-DNAStringSet(unlist(H4_seq))
names(H4_seq1)<-"Human_Chr4_CDS"
M3_seq1<-DNAStringSet(unlist(M3_seq))
names(M3_seq1)<-"Mouse_Chr3_CDS"
M5_seq1<-DNAStringSet(unlist(M5_seq))
names(M5_seq1)<-"Mouse_Chr5_CDS"
M8_seq1<-DNAStringSet(unlist(M8_seq))
names(M8_seq1)<-"Mouse_Chr8_CDS"
Mouse_gen_CDS<-DNAStringSet(c(M3_seq1, M5_seq1, M8_seq1))


#Comparing CDS Length with Chr Length
length_comp1<-length(H4_seq1[[1]])/seqlengths(Hsapiens)[4]
length_comp2<-length(M3_seq1[[1]])/seqlengths(Mmusculus)[3]
length_comp3<-length(M5_seq1[[1]])/seqlengths(Mmusculus)[5]
length_comp4<-length(M8_seq1[[1]])/seqlengths(Mmusculus)[8]


GC_content_CDH4 <- letterFrequency(H4_seq1[[1]],  letters="CG",as.prob=T)
GC_content_CDM3 <- letterFrequency(M3_seq1[[1]],  letters="CG",as.prob=T)
GC_content_CDM5 <- letterFrequency(M5_seq1[[1]],  letters="CG",as.prob=T)
GC_content_CDM8 <- letterFrequency(M8_seq1[[1]],  letters="CG",as.prob=T)


Ch_name<-rep(c("Human_chr4","Mouse_chr3", "Mouse_chr5", "Mouse_chr8"),2)
GC_compare<-c(GC_content_H4,GC_content_M3,GC_content_M5,GC_content_M8,
              GC_content_CDH4,GC_content_CDM3,GC_content_CDM5,GC_content_CDM8)
groups<-rep(c("Whole_Chr", "CDS"),each=4)
T<-data.frame(Ch_name,"GC_Content"=GC_compare ,groups)
ggplot(T, aes(x=Ch_name, y=GC_Content, fill=groups)) +
  geom_bar(stat="identity",position="dodge",color="black")+scale_fill_brewer(palette="Pastel1")+
  geom_text(aes(label=round(GC_Content,2)), vjust=1.5,position=position_dodge(.9), size=4)

S<-data.frame(Ch_name[1:4], length_compare=c(length_comp1,length_comp2,length_comp3,length_comp4)
              ,Chr_length=c(seqlengths(Hsapiens)[4],seqlengths(Mmusculus)[3],seqlengths(Mmusculus)[5],seqlengths(Mmusculus)[8]))
ggplot(S, aes(x=Ch_name.1.4., y=Chr_length, size=length_compare))+geom_point()
dev.copy(pdf,"GC_compare.pdf", width=8, height=6)
dev.off()
##........................................................................................
# Problem 5
require(genoPlotR)
H_Chr4<-DNAStringSet(Hsapiens$chr4 )
names(H_Chr4)<-"Human_Chr4"
M_Chr3<-DNAStringSet(Mmusculus$chr3)
names(M_Chr3)<-"Mouse_Chr3"
M_Chr5<-DNAStringSet(Mmusculus$chr5 )
names(M_Chr5)<-"Mouse_Chr5"
M_Chr8<-DNAStringSet(Mmusculus$chr8 )
names(M_Chr8)<-"Mouse_Chr8"
Mouse_gen<-DNAStringSet(c(M_Chr3, M_Chr5, M_Chr8))

writeXStringSet(H_Chr4, "H_Chr4.fa")
writeXStringSet(M_Chr3, "M_Chr3.fa")
writeXStringSet(M_Chr5, "M_Chr5.fa")
writeXStringSet(M_Chr8, "M_Chr8.fa")
writeXStringSet(Mouse_gen, "Mouse_gen.fa")


cmd<-"makeblastdb -in H_Chr4.fa -dbtype nucl -title Chr4_Human -out Chr4_Human -parse_seqids"
cmd1<-"blastn -query Mouse_gen.fa -db Chr4_Human -out Human_Mouse.txt -outfmt 6 -evalue 1e-50" 
system(cmd)
system(cmd1)

# Designing filter. you can choose Alignment length and which_mouse Chromosom 
# as input to this filter
filter_align=function(L,which_chrom){
  # L Can be chosen 100, 1000 or Other Numbers
  # which_chrom can be chosen as c("Mouse_Chr3","Mouse_Chr5","Mouse_Chr8") and determine Mouse Chr  
  U<- read.table("Human_Mouse.txt")
  colnames(U)<- c("q.id","s.id","identity","align_L","mismacth","gap open","q.start",
                  "q.end","s.start","s.end", "evalue","bitscore")
  U<-U[U$align_L>L,]            
  U<-U[U$q.id==which_chrom,]
  return(U)
}

# Designing Filter For Improving Representation of genoplot (Removing Repeat Regions in Query and Subject)
filter_rep<-function(N)
{
  #N can be The Table (Output of blast or output of blast after filter_align)
  a<-d<-list()
  v<-rep(1,nrow(N))
    for (i in 1:nrow(N)){
      A1=N$q.start[i]
      B1=N$q.end[i]
      A2=N$s.start[i]
      B2=N$s.end[i]
      z<-(N$q.start==A1 & N$q.end==B1)
      z1<-(N$s.start==A2 & N$s.end==B2)
      if (sum(z)>1){
        a[[i]]<-which(z==TRUE)
      }
      if (sum(z1)>1){
        d[[i]]<-which(z1==TRUE)
    }
  }
  b<-unique(c(unlist(a),unlist(d)))
  v[b]=0
  U2<-N[v==1,]
  return(U2)
}
# if you are dealing with Mouse Chr3, Choose "M_Chr3.fa" and so on.
whi_chrom1<-"M_Chr3.fa"     #It can be chosen by User from c("M_Chr3.fa","M_Chr5.fa","M_Chr8.fa")
U<-filter_align(1000, "Mouse_Chr3")
U<-filter_rep(U)

write.table(U, "Human_Mouse2.txt",row.names = FALSE,col.names = FALSE, sep="\t")
E<-read_comparison_from_blast("Human_Mouse2.txt")
s2<-read_dna_seg_from_fasta(whi_chrom1)
s3<-read_dna_seg_from_fasta("H_Chr4.fa")

plot_gene_map(dna_segs=list(s2,s3),comparisons=list(E))

dev.copy(pdf,"genoplot_Human_Mouse.pdf", width=8, height=6)
dev.off()
##.......................................................................................
## Ploting Genoplot For CDS Sites
writeXStringSet(H4_seq1, "H_Chr4_CDS.fa")
writeXStringSet(M3_seq1, "M_Chr3_CDS.fa")
writeXStringSet(M5_seq1, "M_Chr5_CDS.fa")
writeXStringSet(M8_seq1, "M_Chr8_CDS.fa")
writeXStringSet(Mouse_gen_CDS, "Mouse_gen_CDS.fa")


cmd2<-"makeblastdb -in H_Chr4_CDS.fa -dbtype nucl -title Chr4_Human_CDS -out Chr4_Human_CDS -parse_seqids"
cmd3<-"blastn -query Mouse_gen_CDS.fa -db Chr4_Human_CDS -out HM_CDS.txt -outfmt 6 -evalue 1e-50" 
system(cmd2)
system(cmd3)

filter_align_CDS=function(L,which_chrom){
  # L Can be chosen 100, 1000 or Other Numbers
  #which_chrom can be chosen by User from c("Mouse_Chr3_CDS","Mouse_Chr5_CDS","Mouse_Chr8_CDS")
  U<- read.table("HM_CDS.txt")
  colnames(U)<- c("q.id","s.id","identity","align_L","mismacth","gap open","q.start",
                  "q.end","s.start","s.end", "evalue","bitscore")
  U<-U[U$align_L>L,]            
  U<-U[U$q.id==which_chrom,]
  return(U)
}

whi_chr_CDS<-"M_Chr8_CDS.fa"   #It can be chosen by User c("M_Chr3_CDS.fa","M_Chr5_CDS.fa","M_Chr8_CDS.fa")
U1<-filter_align_CDS(1000, "Mouse_Chr8_CDS")
U1<-filter_rep(U1)

write.table(U1, "HM_CDS1.txt",row.names = FALSE,col.names = FALSE, sep="\t")
E<-read_comparison_from_blast("HM_CDS1.txt")
first_Chr<-read_dna_seg_from_fasta(whi_chr_CDS)
sec_Chr<-read_dna_seg_from_fasta("H_Chr4_CDS.fa")
plot_gene_map(dna_segs=list(first_Chr,sec_Chr),comparisons=list(E))

##........................................................................................
#Section 5 (Problem 5)
# For Whole Chromosoms
V1<-filter_align(1000, "Mouse_Chr3")
V1<-filter_rep(V1)
V2<-filter_align(1000, "Mouse_Chr5")
V2<-filter_rep(V2)
V3<-filter_align(1000, "Mouse_Chr8")
V3<-filter_rep(V3)

# if you want to use scatter plot for each Chromosoms:
#Human_start_point<-V1$s.start;  Mouse_start_point<- V1$q.start

Human_start_point<-c(V1$s.start,V2$s.start,V3$s.start) 
Mouse_start_point<- c(V1$q.start,V2$q.start+seqlengths(Mmusculus)[3]
                      , V3$q.start+seqlengths(Mmusculus)[3]+seqlengths(Mmusculus)[5])
align_L<-c(V1$align_L,V2$align_L,V3$align_L) 
V4<- data.frame(Human_start_point,Mouse_start_point, align_L)
ggplot(V4, aes(x=Human_start_point, y=Mouse_start_point, size=align_L)) + 
  geom_point()

#For Coding Sites
V5<-filter_align_CDS(100, "Mouse_Chr3_CDS")
V5<-filter_rep(V5)
V6<-filter_align_CDS(100, "Mouse_Chr5_CDS")
V6<-filter_rep(V6)
V7<-filter_align_CDS(100, "Mouse_Chr8_CDS")
V7<-filter_rep(V7)

Human_start_point1<-c(V5$s.start,V6$s.start,V7$s.start) 
Mouse_start_point1<- c(V5$q.start,V6$q.start+length(M3_seq[[1]])
                      , V7$q.start+length(M3_seq[[1]])+length(M5_seq[[1]]))
align_L1<-c(V5$align_L,V6$align_L,V7$align_L) 
V8<- data.frame(Human_start_point1,Mouse_start_point1, align_L1)
ggplot(V8, aes(x=Human_start_point1, y=Mouse_start_point1, size=align_L1)) + 
  geom_point()


