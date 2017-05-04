library('DESeq2')
library('tximport')
#C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo
#Require a named vector of files to import
files=c("C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\scer_ribo1.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\scer_ribo2.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\Spar_ribo1.tsv","C:\\Users\\tayaz\\Desktop\\Seminar\\final\\R_ribo\\Spar_ribo2.tsv")

names(files)=c("scer_ribo1","scer_ribo2","Spar_ribo1","Spar_ribo2")
files

#Reading in files
txdat = tximport(files, type ="kallisto", txOut=TRUE)

#Conditions are the rep names for the yeasts: e.g, S.cerevisiae_1 
coldata=data.frame(condition=c("scer_ribo_1","scer_ribo_2","spar_ribo_1","spar_ribo_2"))

rownames(coldata) =names(files)
coldata

#Grouping samples by the condition column of coldata
dds = DESeqDataSetFromTximport(txdat, colData=coldata, design=~ condition)

#Get the expressions
des=DESeq(dds)
res = results(des)

#write it to a file
write.table(res,file="ribo_seq_result.txt", sep="\t", quote=FALSE)

