library('DESeq2')
library('tximport')
files=c("C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\scer1\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\scer2\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\spar1\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\spar2\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\scer_ribo1.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\scer_ribo2.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\Spar_ribo1.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\Spar_ribo2.tsv")

names(files)=c("scer_rna1","scer_rna2","Spar_rna1","Spar_rna2","scer_ribo1","scer_ribo2","Spar_ribo1","Spar_ribo2")

txdat = tximport(files, type ="kallisto", txOut=TRUE)

coldata=data.frame(condition=c("scer_rna_1","scer_rna_2","spar_rna_1","spar_rna_2","scer_ribo_1","scer_ribo_2","spar_ribo_1","spar_ribo_2"))
rownames(coldata) =names(files)
coldata
dds2 = DESeqDataSetFromTximport(txdat, colData=coldata, design=~condition)
#Get the expressions
des=DESeq(dds2)
res = results(des)
plotMA(res, ylim=c(-2,2))



write.table(res,file="te_result.txt", sep="\t", quote=FALSE)
