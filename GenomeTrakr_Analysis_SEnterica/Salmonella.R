#Database was downloaded on 6-22-20
#----------------------
GenomeTrakrData<-read.csv(file = "pathogens_SE.csv",sep = ",",header=FALSE,fill=TRUE,stringsAsFactors = FALSE)
colnames(GenomeTrakrData)=c("Organism Group","Strain","Serovar","Isolate","Create Date","Assembly","N50","Level","TaxID","Length","SRA Release Date")
GenomeTrakrData<-GenomeTrakrData[-1,]
GenomeTrakrData$N50<-as.numeric(GenomeTrakrData$N50)
GenomeTrakrData$Length<-as.numeric(GenomeTrakrData$Length)
#Determining quality of this dataset and what is annotated versus what is not
sum(GenomeTrakrData$Serovar == "")
# 133027/281550 = 47% do not have a serovar assigned
length(table(GenomeTrakrData$Serovar))
#Out of 133027 serovar labels, 1548 different labels exist
sum(GenomeTrakrData$Assembly == "")
# 22082/281550 = 8% do not have an assembly
# 259468 have an assembly

#For figuring out the distribution of levels, Scaffold = Half-filled circle, Chromsome = 3/4-filled circle
#Contig = 1/4 filled circle, Complete Genome = Full filled circle

GTD_Scaffold<-GenomeTrakrData[which(GenomeTrakrData$Level=="Scaffold"),] #3481 Assemblies
GTD_Chromosome<-GenomeTrakrData[which(GenomeTrakrData$Level=="Chromosome"),] #179 Assemblies
GTD_Contig<-GenomeTrakrData[which(GenomeTrakrData$Level=="Contig"),] #254871 assemblies
GTD_Complete_Genome<-GenomeTrakrData[which(GenomeTrakrData$Level=="Complete Genome"),] #936 Assemblies
save(GTD_Chromosome,GTD_Complete_Genome,GTD_Contig,GTD_Scaffold,file = "GenomeTrakrData_Raw.Rdata")

setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/16s_Files")
which(Assemblies_For_Download=="GCA_012530335.1")
y<-364
for(y in Assemblies_For_Download[366])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  RefSeq_DataFile<-entrez_fetch(db = "nuccore",id = refseq_data[[1]][1][[1]][1],rettype = "gb")
  write(RefSeq_DataFile,file=paste(y,".txt",sep = ""))
  Sys.sleep(2)
}

Files_To_Parse<-list.files(path = "~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/16s_Files")
source('~/Desktop/GenomeTrakr_Analysis/Scripts/Parse_16s_Numbers.R')
Raw_DataTable16s<-map_dfr(Files_To_Parse, Parse_16s_Numbers)
DataTable_16s<-as.data.frame(Raw_DataTable16s[,2:4])
Numbers_16s_GoodAssemblies<-as.numeric(as.vector(DataTable_16s[,2]))
hist(Numbers_16s_GoodAssemblies)
#Out of 1115 files, 1090 were readable andcontain rRNA data. 
#Table of resulting data
table(Numbers_16s_GoodAssemblies)
#1059 Assemblies have 7 16s sequences
#-----

#First filter out the files which do not have 7 6s copy numbers
Files_With_7Copies<-Files_To_Parse[which(Numbers_16s_GoodAssemblies==7)]
Files_With_7Copies_Fixed<-gsub(x = Files_With_7Copies,pattern = ".txt",replacement = "")

#Need to fetch all of the FASTAs, this loop will pull the largest complete genome FASTA files based 
#on sequence length to make sure the biggest assembly is downloaded
y<-Files_With_7Copies_Fixed[1]

setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Assembly_FASTAs")
source('~/Desktop/GenomeTrakr_Analysis/Scripts/Sequence_Length_Check.R')
for(y in Files_With_7Copies_Fixed)
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  CDS_FASTA<-entrez_fetch(db = "nuccore",id = names(Sequence_Lengths_Sorted[1]),rettype = "fasta")
  write(CDS_FASTA,file=paste(y,".txt",sep = ""))
}

which(Files_With_7Copies_Fixed=="GCA_002777175.1")
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/CDS_FASTAs")
for(y in Files_With_7Copies_Fixed[509:1115])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  CDS_FASTA<-entrez_fetch(db = "nuccore",id = names(Sequence_Lengths_Sorted[1]),rettype = "fasta_cds_na")
  write(CDS_FASTA,file=paste(y,".txt",sep = ""))
}


#-------- working with core genes now
Potential_CG_List<-as.list(NULL)
for(i in 1:1059)
{
  fastaFile <- readDNAStringSet(list.files("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Assembly_FASTAs")[i],format = "fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  SplitFASTA <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
  PotentialCGs<-str_extract(string = SplitFASTA[,1],pattern = "(?<=\\[gene=)[:alpha:]*")
  Potential_CG_List[[i]]<-PotentialCGs
}

Potential_CG_List_Collapsed<-map(.x = Potential_CG_List,.f = ~paste(.,collapse=""))

Potential_Core_Genes<-Potential_CG_List[[6]]
Potential_Core_Genes2<-na.omit(Potential_Core_Genes)

TestList<-as.list(NULL)

for(i in Potential_Core_Genes2)
{
  True_Value<-table(map_chr(.x = Potential_CG_List_Collapsed,.f = ~str_detect(string = .,pattern = i)))[2]
  TestList[[i]]<-True_Value
}
TestList2<-unlist(TestList)
table(TestList2)
TestList3<-TestList2[which(TestList2>1056)]
SE_Core_Genes_Pre<-names(TestList3)
SE_Core_Genes<-gsub(pattern = ".TRUE",replacement = "",x = SE_Core_Genes_Pre)
#core genes list is set now
###########################

library("Biostrings")
SE_Core_Genes_Adj<-paste("gene=",SE_Core_Genes,sep="")
zz<-1:1059

Core_Genes_List<-as.list(NULL)
source('~/Desktop/GenomeTrakr_Analysis/Scripts/Core_Gene_Extractor.R')

for(i in zz)
{
  fastaFile <- readDNAStringSet(list.files("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/CDS_FASTAs")[i],format = "fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  SplitFASTA <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
  Core_Genes_Table<-map_dfr(SE_Core_Genes_Adj,Core_Gene_Extractor)
  Core_Genes_Table_Final<-Core_Genes_Table[!duplicated(Core_Genes_Table$seq_name),]
  FixedGeneNames<-str_extract(string = Core_Genes_Table_Final[,1],pattern = "(?<=\\[gene=)[:alpha:]*")
  AssemblyName<-paste(list.files("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/CDS_FASTAs")[i])
  AssemblyNameFixed<-gsub(pattern = ".txt",replacement = "",x = AssemblyName)
  Core_Genes_Table_Final[,1]<-paste(AssemblyNameFixed,FixedGeneNames,sep = " ")
  Core_Genes_Table_Final2<-Core_Genes_Table_Final[!(duplicated(Core_Genes_Table_Final[,1]) | duplicated(Core_Genes_Table_Final[,1], fromLast = TRUE)), ]
  Core_Genes_List[[i]]<-Core_Genes_Table_Final2
}

save(Core_Genes_List,file = "Core_Genes_List_SE.Rdata")

Core_Genes_List_Cut<-Core_Genes_List[Assemblies_With_CGs_Index]

#we want to apply over data frames first - then apply over possible gene names
source('~/Desktop/GenomeTrakr_Analysis/Scripts/CoreGeneBuilder.R')
Core_Gene_Assemblies<-map(.x = SE_Core_Genes, .f = ~ Core_Gene_Builder(.))
save(Core_Gene_Assemblies,file = "Core_Gene_Assemblies.Rdata")
Core_Gene_Assemblies_NoZero<-map_int(.x = Core_Gene_Assemblies,.f = ~ length(rownames(.)))
#we are going to focus on 89 core genes, leaving us with 1057 total assemblies

source('~/Desktop/GenomeTrakr_Analysis/Scripts/Assemblies_Missing_Genes.R')
Assemblies_To_Keep<-map_int(.x = Files_With_7Copies_Fixed,.f = Assemblies_Missing_Genes)
names(Assemblies_To_Keep)<-Files_With_7Copies_Fixed
Assemblies_With_CGs_Index<-which(Assemblies_To_Keep==89)
#1045 total assemblies remain

#second version to focus on final stuff 
Core_Gene_Assemblies_Cut<-Core_Gene_Assemblies[which(Core_Gene_Assemblies_NoZero==1045)]
source('~/Desktop/GenomeTrakr_Analysis/Scripts/FixColNames.R')
Core_Gene_Assemblies_Cut2<-map(.x = Core_Gene_Assemblies_Cut,.f = FixColNames)

#final is 1045 assemblies used, 109 core genes included


setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Alignment_FASTAs")
for(i in 1:109)
{
  Data_Table<-Core_Gene_Assemblies_Cut2[[i]]
  dat2fasta(dat = Data_Table,outfile = paste(i,".fasta",sep = ""))
}

#generating mafft alignments using command line
for(i in 11:109)
{
  system(command = paste("/usr/local/bin/mafft --retree 2 --inputorder /users/meech/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Alignment_FASTAs/",i,".fasta > /users/meech/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/MAFFT_Outputs/",i,"_aligned.fasta",sep = ""))
}

MAFFT_List<-as.list(NULL)

for (i in 1:109)
{
  fastaFile <- readDNAStringSet(paste("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/MAFFT_Outputs/",i,"_aligned.fasta",sep=""),format = "fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  Aligned_FASTA <- data.frame(seq_name, sequence,stringsAsFactors = FALSE)
  MAFFT_List[[i]]<-Aligned_FASTA
}

StringMaker<-function(xx)
{
  FullString<-map_chr(.x = MAFFT_List,.f = ~paste(.[xx,2]))
  FullString2<-paste(FullString,sep="",collapse="")
  return(data.frame(FullString2))
}

Test2<-Core_Gene_Assemblies_Cut2[[1]][,1]
Test<-1:1045
Maybe<-map_dfr(.x = Test,.f = StringMaker)
MaybeFinal<-cbind(Test2,Maybe)
colnames(MaybeFinal)=c("seq.name","seq.text")
dat2fasta(dat = MaybeFinal,outfile = "RAXML_Input.fasta")


source('~/Desktop/GenomeTrakr_Analysis/Scripts/Sequence_Length_Check.R')
for(y in Files_With_7Copies_Fixed[11:1059])
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  Nuc_Input<-names(Sequence_Lengths_Sorted[1])
  system(command = paste('wget -O GFF3s/',y,'.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=',Nuc_Input,'"',sep = ""))
}

#after running raxml, reinput the data
RAXML_Output<-readChar("RAxML_result.test.out2",nchars = 2e7)
vert.tree<-read.tree(text=RAXML_Output)

#just in case
save(Core_Genes_List_Cut,Potential_CG_List,file = "Backup.Rdata")


#Getting into constructing the 16s parsing table
#Use Assembly_Names_Subset fo downloads
#Use Limited_FASTAs to test

source('~/Desktop/GenomeTrakr_Analysis/Scripts/Sequence_Length_Check.R')
for(y in Assembly_Names_Subset)
{
  assembly_search<-entrez_search(db = "assembly",term=y)
  refseq_data<-entrez_link(dbfrom = "assembly",id = assembly_search$ids,db = "nucleotide")
  Sequences_To_Scan<-refseq_data[[1]][1][[1]]
  Sequence_Length_Results<-map_int(Sequences_To_Scan, Sequence_Length_Check)
  names(Sequence_Length_Results)<-Sequences_To_Scan
  Sequence_Lengths_Sorted<-sort(Sequence_Length_Results,decreasing = TRUE)
  Nuc_Input<-names(Sequence_Lengths_Sorted[1])
  system(command = paste('wget -O GFF3s/',y,'.gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=',Nuc_Input,'"',sep = ""))
}

#first need to rename headers of raw fasta files
library("seqinr", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files")
for(y in Files_With_7Copies_Fixed)
{
  FASTA_To_Rename<-readDNAStringSet(filepath = paste("Raw_FASTAs/",y,".txt",sep=""),format="fasta")
  names(FASTA_To_Rename)<-"chr"
  write.fasta(sequences = paste(FASTA_To_Rename),names = names(FASTA_To_Rename),file.out = paste("Raw_FASTAs_CHR/",y,".txt",sep=""))
}

#removing files that wouldnt download left us with 1088 / 1109
#now looping to generate BED files
library("stringr", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files")
GFFs_To_Process_Raw<-list.files(path = "GFF3s/")
GFFs_To_Process<-gsub(pattern = ".gff3",replacement = "",x = GFFs_To_Process_Raw)

for(y in GFFs_To_Process[515:1059])
{
  Test<-read.csv(file = paste("GFF3s/",y,".gff3",sep=""),sep = "\r",header = FALSE)
  Test2<-as.data.frame(Test[which(str_detect(string = Test$V1,pattern = "16S ribosomal RNA")),])
  Test3<-as.data.frame(Test2[c(1,3,5,7,9,11,13),])
  Position_Table<-map_dfc(.x = Test3[,1],.f = ~str_split(string = .,pattern = "\t"))
  Position_Table_Cut<-t(Position_Table[4:5,1:7])
  First_Column<-rep("chr",7)
  Last_Column<-c("16s_1","16s_2","16s_3","16s_4","16s_5","16s_6","16s_7")
  Final_Position_Table<-as.data.frame(cbind(First_Column,Position_Table_Cut,Last_Column))
  write.table(x = Final_Position_Table,file = paste("GFF3s_Positions/",y,".bed",sep=""),sep="\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
#bedtools section to make 16s data table

GFFs_To_Process_Raw2<-list.files(path = "GFF3s_Positions/")
GFFs_To_Process2<-gsub(pattern = ".bed",replacement = "",x = GFFs_To_Process_Raw2)

setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files")
for(y in GFFs_To_Process2)
{
  system(command = paste("bedtools getfasta -fi Raw_FASTAs_CHR/",y,".txt -bed GFF3s_Positions/",y,".bed -fo 16s_Sequences_BT/",y,".txt",sep=""))
  system(command = paste("rm Raw_FASTAs_CHR/",y,".txt.fai",sep=""))
}

#script to merge 16s files into one data frame
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/16s_Sequences_BT")
source('~/Desktop/GenomeTrakr_Analysis/Scripts/Merger_16s.R')
Full_16s_Data<-map(GFFs_To_Process2, Merger_16s)

#need to fix bad nucleotides
source('~/Desktop/GenomeTrakr_Analysis/Scripts/GSubber.R')
Full_16s_Data_Fixed<-map(.x = Full_16s_Data,.f = ~map_chr(.x = .,.f = GSubber))

#Finding index values for strands with - that need to be RCed
source('~/Desktop/GenomeTrakr_Analysis/Scripts/Strand_Searcher.R')
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files")
Strand_Indices<-map(.x = GFFs_To_Process2,.f = Strand_Searcher)
#actual functions that do the reverse complement fixing - input range equal to number of files
source('~/Desktop/GenomeTrakr_Analysis/Scripts/ReverseComplementFixer.R')
source('~/Desktop/GenomeTrakr_Analysis/Scripts/ReverseComplementer.R')

Full_16s_Data_Fixed2<-Full_16s_Data_Fixed[c(1:21,23:1057)]
Full_16s_Data_Fixed<-Full_16s_Data_Fixed2
Strand_Indices2<-Strand_Indices[c(1:21,23:1057)]
Strand_Indices<-Strand_Indices2

Full_16s_Data_Fixed2<-Full_16s_Data_Fixed[c(1:939,941:1056)]
Full_16s_Data_Fixed<-Full_16s_Data_Fixed2
Strand_Indices2<-Strand_Indices[c(1:939,941:1056)]
Strand_Indices<-Strand_Indices2
GFFs_To_Process3<-GFFs_To_Process2[c(1:21,23:940,942:1057)]

Full_16s_Data_RC<-map(.x = 1:1055,.f = ReverseComplementFixer)
save(Full_16s_Data_RC,file = "Full_16s_Data_RC.Rdata")

library("dada2", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

## New scoring system idea
source('~/Desktop/GenomeTrakr_Analysis/Scripts/Advanced_HAM_Calculator.R')

#strings for testing 
Sequence_Of_Interest<-Full_16s_Data_RC[[which(list.files("Data_Files/16s_Sequences_BT")=="GCA_000803705.1.txt")]]
#Meat sample sequences
Seq1<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGATTAGGTCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)
Seq2<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGATTAGGTCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)
Seq3<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAGCACAGAGAGCTTGCTCTCGGGTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTTGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAAGAATCCAGAGATGGATTTGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)
Seq4<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGCGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAATTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)
Seq5<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)
Seq6<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCACAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAATTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)
Seq7<-c("TCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTGGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCAGATGTGCCCAGATGGGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGTGTTGTGGTTAATAACCGCAGCAATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTTGAGTCTTGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACAAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCTACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTAGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAAGAATCCAGAGATGGATTTGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTTAGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTG",1)

MultiSequence_To_Test<-as.data.frame(c(Sequence_Of_Interest[2:8]))

Random<-c(1,1,1,1,1,1,1)
MultiSequence_To_Test<-cbind(Test_Sequence,Random)
#actual function to form table of sums - worse but faster version
Multi_Sum_Table<-map_dfc(.x = 1:7,.f = ~Advanced_HAM_Summer(TestSequence = MultiSequence_To_Test[.,1],Multi = MultiSequence_To_Test[.,2]))
Multi_Sum_Final<-map_int(.x = 1:1055,.f = ~sum(Multi_Sum_Table[.,]))

Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
Extra_Column<-rep(x = c(1),times = 1055)
Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")

Test_Sequence<-MultiSequence_To_Test[,1]
#fancier version
StartTime<-Sys.time()
#This line executes the above functions and applies them over the database of reference assembly 16-23s data, setting the range of reference assemblies to test with .x = 
Multi_Sum_Final2<-map_int(.x = 1:1055,.f = ~Super_HAM_TableMaker(TestSequence = Test_Sequence,ReferenceSequence = Full_16s_Data_RC[[.]]))
EndTime<-Sys.time()
EndTime-StartTime

#testing other assemblies
Sequence_Of_Interest<-Full_16s_Data_RC[[which(list.files("Data_Files/16s_Sequences_BT")=="GCA_000756465.1.txt")-1]]
Test_Sequence<-Sequence_Of_Interest[2:8]

#making taxid table
BigTable<-rbind(GTD_Chromosome,GTD_Complete_Genome)
TaxID_Table<-BigTable[,c(6,9)]
Interesting_Indices<-c(which(TaxID_Table$TaxID=="108619"),which(TaxID_Table$TaxID=="143221"),which(TaxID_Table$TaxID=="2583588"),which(TaxID_Table$TaxID=="28150"),which(TaxID_Table$TaxID=="28901"),which(TaxID_Table$TaxID=="58095"),which(TaxID_Table$TaxID=="58712"),which(TaxID_Table$TaxID=="59201"),which(TaxID_Table$TaxID=="595"),which(TaxID_Table$TaxID=="611"),which(TaxID_Table$TaxID=="90370"))

which(TaxID_Table$TaxID=="108619")
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="149539")]<-"Orange"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="143221")]<-"Purple1"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="2583588")]<-"Purple2"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="28150")]<-"Purple3"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="28901")]<-"Purple4"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="58095")]<-"Purple5"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="58712")]<-"Purple6"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="59201")]<-"Purple7"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="595")]<-"Purple8"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="611")]<-"Purple9"
TaxID_Table$TaxID[which(TaxID_Table$TaxID=="90370")]<-"Purple10"

TaxID_Table$TaxID[!Interesting_Indices]<-"Blank"

TaxID_Table_Cut<-TaxID_Table[Interesting_Indices,]


##########making comparisons of different distance metrics
#these lines make the scatter plot
CoP_Phylo<-cophenetic.phylo(x = vert.tree)
which(vert.tree$tip.label=="GCA_001890445.1")
Phylo_Distances<-CoP_Phylo[699,]
PD_Labels<-as.data.frame(cbind(vert.tree$tip.label,Phylo_Distances))
PD_Labels$Phylo_Distances<-as.numeric(PD_Labels$Phylo_Distances)

Multi_Mapping_Table_Plot_Test<-Multi_Mapping_Table_Plot
rownames(Multi_Mapping_Table_Plot_Test)=Multi_Mapping_Table_Plot_Test[,1]
rownames(PD_Labels)=PD_Labels[,1]
Scatter_Table<-merge(Multi_Mapping_Table_Plot_Test,PD_Labels,by=0,all=TRUE)


range<-Multi_Mapping_Table_Plot$value
range2<-range[which(range<150)]

#generate distances:
#minhash first
i<-list.files()[1]
for(i in list.files()[208:1060])
{
  system(command = paste("~/Desktop/mash-OSX64-v2.2/mash dist GCA_001890445.1.txt ",i," > ~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/MinHASH_Calcs/",i,sep = ""))
}
MASH_List<-as.list(NULL)
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/MinHASH_Calcs")
for(j in 1:1059)
{
  MASH_Table<-read.table(file=list.files()[j])
  MASH_List[[j]]<-MASH_Table[,3]
}
MASH_Names<-list.files()
MASH_Names2<-gsub(pattern = ".txt",replacement = "",x = MASH_Names)
MASH_Scores<-unlist(MASH_List)
MASH_Final_Table<-cbind(MASH_Names2,MASH_Scores)
rownames(MASH_Final_Table)<-MASH_Final_Table[,1]
rownames(Scatter_Table)<-Scatter_Table[,1]
Scatter_Table2<-merge(Scatter_Table,MASH_Final_Table,by=0,all=TRUE)

#now try to generate blast scores
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Raw_FASTAs")

for(i in list.files()[22:1060])
{
  system(command = paste('/usr/local/ncbi/blast/bin/blastn -query GCA_001890445.1.txt -subject ',i,' -task blastn -outfmt "6 pident nident" > ~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/BLAST_Calcs/',i,sep = ''))
}

library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cores=detectCores()
registerDoParallel(cores)

foreach(i = list.files()[101:1060]) %dopar% {
  system(command = paste('/usr/local/ncbi/blast/bin/blastn -query GCA_001890445.1.txt -subject ',i,' -task blastn -outfmt "6 pident nident" > ~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/BLAST_Calcs/',i,sep = ''))
}

BLAST_List<-as.list(NULL)
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/BLAST_Calcs")
j<-1
for(j in 22:729)
{
  BLAST_Table<-read.table(file=list.files()[j])
  Test2<-BLAST_Table$V1*BLAST_Table$V2
  Test3<-sum(Test2)
  Test4<-Test3/sum(BLAST_Table$V2)
  BLAST_List[[j]]<-Test4
}
BLAST_Names<-list.files()
BLAST_Names2<-gsub(pattern = ".txt",replacement = "",x = BLAST_Names)
BLAST_Scores<-unlist(BLAST_List)
BLAST_Final_Table<-cbind(BLAST_Names2,BLAST_Scores)
rownames(BLAST_Final_Table)<-BLAST_Final_Table[,1]
rownames(Scatter_Table2)<-Scatter_Table2[,1]
Scatter_Table3<-merge(Scatter_Table2,BLAST_Final_Table,by=0,all=TRUE)

Scatter_Table_Trim<-Scatter_Table3[,c(4,5,8,10,12)]
colnames(Scatter_Table_Trim)<-c("label","Score_16s","Phylo_Distance","MASH_Distance","BLAST_Distance")
Scatter_Table_Trim$MASH_Distance<-as.numeric(Scatter_Table_Trim$MASH_Distance)
Scatter_Table_Trim$BLAST_Distance<-as.numeric(Scatter_Table_Trim$BLAST_Distance)
Scatter_Table_Trim$BLAST_Distance<-100-Scatter_Table_Trim$BLAST_Distance

#testing section
ggplot(Scatter_Table_Trim, aes(x=MASH_Distance, y=BLAST_Distance)) +
  geom_point(size=2, shape=23) 


#phylo distance: 0 -> 0.7
#score_16s: 0 -> 100

#Trying to color some serovars
TaxID_Table2<-High_Quality_Assemblies[,c(6,3)]
Good_Serovars<-c(which(TaxID_Table2$Serovar=="Enteritidis"),which(TaxID_Table2$Serovar=="Typhi"),which(TaxID_Table2$Serovar=="Typhimurium"),which(TaxID_Table2$Serovar=="Infantis"),which(TaxID_Table2$Serovar=="Heidelberg"),which(TaxID_Table2$Serovar=="Newport"),which(TaxID_Table2$Serovar=="1,4,[5],12:i:-"),which(TaxID_Table2$Serovar=="Anatum"),which(TaxID_Table2$Serovar=="Bareilly"),which(TaxID_Table2$Serovar=="Montevideo"),which(TaxID_Table2$Serovar=="Senftenberg"),which(TaxID_Table2$Serovar=="Indiana"),which(TaxID_Table2$Serovar=="Agona"),which(TaxID_Table2$Serovar=="Saintpaul"),which(TaxID_Table2$Serovar=="Dublin"),which(TaxID_Table2$Serovar=="H58"),which(TaxID_Table2$Serovar=="Thompson"),which(TaxID_Table2$Serovar=="Paratyphi A"),which(TaxID_Table2$Serovar=="Kentucky"),which(TaxID_Table2$Serovar=="4,[5],12:i:-"))
TaxID_Table3<-TaxID_Table2[Good_Serovars,]

Testing_List<-as.list(NULL)
Assemblies_Of_Interest<-c("GCA_006088735.1","GCA_005144925.1","GCA_009756555.1","GCA_009756455.1","GCA_005160385.1","GCA_003177235.1","GCA_006697045.2","GCA_001953675.1","GCA_001690075.1","GCA_007765995.2","GCA_001623625.1","GCA_011465945.1","GCA_001305835.1","GCA_003052785.2","GCA_001890445.1","GCA_004847885.1","GCA_004848885.1","GCA_004848885.1","GCA_000272815.2","GCA_000625495.2","GCA_000750435.1","GCA_002220345.1","GCA_003184425.1","GCA_003589785.1","GCA_000192085.1","GCA_008504985.2","GCA_000973665.2","GCA_000188955.5","GCA_002953175.1","GCA_900478065.1","GCA_006165225.1","GCA_000503845.1","GCA_004768585.1","GCA_003073535.1","GCA_006517075.1","GCA_012052445.1","GCA_000011885.1","GCA_001048035.2","GCA_900205255.1")
Assemblies_Of_Interest2<-c("GCA_000385905.1","GCA_003718495.1","GCA_001165785.2","GCA_003718615.1","GCA_001104885.3","GCA_000026565.1","GCA_000486165.2","GCA_012054045.1","GCA_001448785.2","GCA_008505445.2","GCA_011057955.1","GCA_000018385.1","GCA_003031875.1","GCA_004358925.1","GCA_001409175.1","GCA_009730015.1","GCA_002386325.1","GCA_003325255.1","GCA_002234515.1","GCA_006165245.1","GCA_003710145.1","GCA_002313105.1","GCA_008727615.1","GCA_003325235.1","GCA_004136095.1","GCA_007972245.1","GCA_002094915.1","GCA_003324815.1","GCA_009648815.1","GCA_009667745.1","GCA_003324935.1","GCA_000486365.2","GCA_003986635.1","GCA_000493295.2","GCA_003325035.1","GCA_005576575.1","GCA_000341425.1","GCA_000756465.1","GCA_006113225.2","GCA_000439255.1")
Full_AssembliesOI<-c(Assemblies_Of_Interest,Assemblies_Of_Interest2)
AssemblyMappingData<-c(Testing_List_FirstSet,Testing_List)

for(x in 1:40)
{
  Name<-Assemblies_Of_Interest2[x]
  TESTSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3==Name)]][2:8]
  Multi_Sum_Final<-map_int(.x = 1:1055,.f = ~Super_HAM_TableMaker(TestSequence = TESTSequence,ReferenceSequence = Full_16s_Data_RC[[.]]))
  Testing_List[[x]]<-Multi_Sum_Final
}
#Making value vs phylo distance plot for test data set

Test_List2<-as.list(NULL)

for(i in 1:79)
{
  Multi_Sum_Final<-AssemblyMappingData[[i]]
  Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
  Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
  Extra_Column<-rep(x = c(1),times = 1055)
  Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
  colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")
  
  TestName<-Full_AssembliesOI[i]
  TIP<-which(vert.tree$tip.label==TestName)
  Phylo_Distances<-CoP_Phylo[TIP,]
  PD_Labels<-as.data.frame(cbind(vert.tree$tip.label,Phylo_Distances))
  PD_Labels$Phylo_Distances<-as.numeric(PD_Labels$Phylo_Distances)
  
  Multi_Mapping_Table_Plot_Test<-Multi_Mapping_Table_Plot
  rownames(Multi_Mapping_Table_Plot_Test)=Multi_Mapping_Table_Plot_Test[,1]
  rownames(PD_Labels)=PD_Labels[,1]
  Scatter_Table<-merge(Multi_Mapping_Table_Plot_Test,PD_Labels,by=0,all=TRUE)
  Scatter_Table_Cut<-Scatter_Table[,c(3,6)]
  Test_List2[[i]]<-Scatter_Table_Cut
}

#Making value vs phylo distance plot for test data set (median version)

Test_List_Median<-as.list(NULL)

for(i in 1:79)
{
  Multi_Sum_Final<-Median_List[[i]]
  Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
  Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
  Extra_Column<-rep(x = c(1),times = 1055)
  Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
  colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")
  
  TestName<-Full_AssembliesOI[i]
  TIP<-which(vert.tree$tip.label==TestName)
  Phylo_Distances<-CoP_Phylo[TIP,]
  PD_Labels<-as.data.frame(cbind(vert.tree$tip.label,Phylo_Distances))
  PD_Labels$Phylo_Distances<-as.numeric(PD_Labels$Phylo_Distances)
  
  Multi_Mapping_Table_Plot_Test<-Multi_Mapping_Table_Plot
  rownames(Multi_Mapping_Table_Plot_Test)=Multi_Mapping_Table_Plot_Test[,1]
  rownames(PD_Labels)=PD_Labels[,1]
  Scatter_Table<-merge(Multi_Mapping_Table_Plot_Test,PD_Labels,by=0,all=TRUE)
  Scatter_Table_Cut<-Scatter_Table[,c(3,6)]
  Test_List_Median[[i]]<-Scatter_Table_Cut
}

Test_List3<-map(.x = Test_List2,.f = ~.[Subset_Indices,])
Big_Frame2<-rbindlist(Test_List3)

Test_List4<-map(.x = Test_List_Median,.f = ~.[Subset_Indices,])
Big_Frame_Median<-rbindlist(Test_List4)

Big_Frame_Median_Scaled<-Big_Frame_Median
Big_Frame_Median_Scaled$value<-Big_Frame_Median_Scaled$value*7

ggplot(Big_Frame2, aes(x=value, y=Phylo_Distances)) +
  geom_point(size=2, shape=23)+
  xlim(0,50)+
  ylim(0,0.2)+
  ggtitle("Heatmap of Phylogenetic Distance vs. 16s Mismatch Scores (Sum)")+
  xlab("16s Mismatch Score (Sum Calculation)")+
  ylab("Phylogenetic Distance")+
  theme(plot.title = element_text(hjust = 0.5))
  
  
ggplot(Big_Frame2, aes(x=value, y=Phylo_Distances) ) +
  geom_hex(bins = 150) +
  scale_fill_continuous(type = "viridis",limits=c(0,10)) +
  theme_bw()+
  xlim(0,250)+
  ylim(0,0.75)+
  ggtitle("Heatmap of Phylogenetic Distance vs. 16s Mismatch Scores (Sum)")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("16s Mismatch Score (Sum Calculation)")+
  ylab("Phylogenetic Distance")

ggplot(Big_Frame_Median, aes(x=value, y=Phylo_Distances)) +
  geom_point(size=2, shape=23)+
  xlim(0,250)+
  ylim(0,1)

ggplot(Big_Frame_Median_Scaled, aes(x=value, y=Phylo_Distances) ) +
  geom_hex(bins = 200) +
  scale_fill_continuous(type = "viridis",limits=c(0,10)) +
  theme_bw()+
  xlim(0,250)+
  ylim(0,0.75)+
  ggtitle("Heatmap of Phylogenetic Distance vs. 16s Mismatch Scores (Median)")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("16s Mismatch Score (Median Calculation)")+
  ylab("Phylogenetic Distance")

#playing around with changes from median points
Big_Frame_Interesting<-Big_Frame2[which(Big_Frame2$value>175),]

GFFs_To_Process2[558]
#GCA_012052445.1
TESTSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_005144925.1")]][2:8]
TESTSequence2<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_003324815.1")]][2:8]

#Examples with weird 16s: GCA_005144985.1, GCA_004919305.1, GCA_005222465.1


ggplot(Big_Frame, aes(x=value, y=Phylo_Distances) ) +
  geom_bin2d(bins = 500) +
  scale_fill_continuous(type = "viridis",limits=c(0,100)) +
  theme_bw()+
  xlim(0,50)+
  ylim(0,0.2)
#loop for recalculating test data set with median values
Median_List<-as.list(NULL)

for(x in 47:79)
{
  Name<-Full_AssembliesOI[x]
  TESTSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3==Name)]][2:8]
  Multi_Sum_Final<-map_int(.x = 1:1055,.f = ~Super_Ham_TableMaker_Median(TestSequence = TESTSequence,ReferenceSequence = Full_16s_Data_RC[[.]]))
  Median_List[[x]]<-Multi_Sum_Final
}

length(Full_AssembliesOI)
TESTSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_005144925.1")]][2:8]
TESTSequence2<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_003324815.1")]][2:8]

Outlier_Points_Median<-which(Big_Frame_Median_Scaled$value>130 & Big_Frame_Median_Scaled$Phylo_Distances<0.3)
Outlier_Points_Sum<-which(Big_Frame2$value>175 & Big_Frame2$Phylo_Distances<0.2)

Outlier_Points_Median<-Outlier_Points_Median[c(1,74:144)]
Outlier_Points_Sum<-Outlier_Points_Sum[c(1:2,75:286)]

Full_AssembliesOI[1]
TestSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_009730015.1")]][2:8]
ReferenceSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_006088735.1")]][2:8]

Row1<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[2], vec=TRUE, band=1000))
Row2<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[3], vec=TRUE, band=1000))
Row3<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[4], vec=TRUE, band=1000))
Row4<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[5], vec=TRUE, band=1000))
Row5<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[6], vec=TRUE, band=1000))
Row6<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[7], vec=TRUE, band=1000))
Row7<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[8], vec=TRUE, band=1000))
Hamming_Table_Result<-rbind(Row1,Row2,Row3,Row4,Row5,Row6,Row7)


