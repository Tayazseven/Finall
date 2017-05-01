import roman
from Bio import SeqIO
import csv
import re

chr_ = []
chr_1 = []
chr_2=[]
chr_f=[]
init=[]
end_=[]
gene=[]
strand=[]

data_file = open("S_cerevisiae_genes.bed", "r") #reading the tsv file
data_reader =csv.reader(data_file, delimiter = '\t')
for row in data_reader:
    chr_.append(row[0])
    init.append(row[1])
    end_.append(row[2])
    gene.append(row[3])
    strand.append(row[5])
for x in chr_:
    line = re.sub('chr', '', x)
    chr_1.append(line)
for x in chr_1:
    chr_2.append(roman.fromRoman(x))
for x in chr_2:
    chr_f.append("Scer_" + str(x))

base = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
def reverse_complement(seq):
    return ''.join(map(base.get, seq.upper()[::-1]))

gene_f=[]
for x in gene:
    gene_f.append(">" + str(x))

transcriptom = []
record = SeqIO.to_dict(SeqIO.parse('S_cerevisiae.fa', 'fasta'))
for i in range(len(init)):
    fwseq = record[chr_f[i]].seq[int(init[i]):int(end_[i])]
    if strand[i] == '+':
        #print "forward sequence:  " + fwseq
        transcriptom.append(fwseq)
    else:
        transcriptom.append(reverse_complement(fwseq))
    #print seq
with open('genes_cerevisiae.fa', 'wb') as f:  # writing to tsv file
    a = csv.writer(f, delimiter='\n')
    a.writerows(zip(gene_f, transcriptom))
