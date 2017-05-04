library('DESeq2')
library('tximport')
directory="C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna"

#Require a named vector of files to import
files=c("C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\scer1\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\scer2\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\spar1\\abundance.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_rna\\spar2\\abundance.tsv")

names(files)=c("scer_rna1","scer_rna2","Spar_rna1","Spar_rna2")
files

#Reading in files
txdat = tximport(files, type ="kallisto", txOut=TRUE)

#Conditions are the rep names for the yeasts: e.g, S.cerevisiae_1 
coldata=data.frame(condition=c("scer_rna_1","scer_rna_2","spar_rna_1","spar_rna_2"))

rownames(coldata) =names(files)
coldata

#Grouping samples by the condition column of coldata
dds = DESeqDataSetFromTximport(txdat, colData=coldata, design=~ condition)

#Get the expressions
des=DESeq(dds)
res = results(des)

#write it to a file
write.table(res,file="mRna_seq_result.txt", sep="\t", quote=FALSE)
