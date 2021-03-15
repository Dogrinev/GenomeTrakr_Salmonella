#GCA_001975245.1
#GCA_000280315.2
#GCA_000487775.2
#GCA_900635725.1
#GCA_901457615.1

setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica")
#Finding String
Sequence_Of_Interest<-Full_16s_Data_RC[[which(list.files("Data_Files/16s_Sequences_BT")=="GCA_000818075.1.txt")]]
Test_Sequence<-Sequence_Of_Interest[2:8]

#fancier version
library("dada2", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
library("tidyverse", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

StartTime<-Sys.time()
#This line executes the above functions and applies them over the database of reference assembly 16-23s data, setting the range of reference assemblies to test with .x = 
Multi_Sum_Final<-map_int(.x = 1:1055,.f = ~Super_HAM_TableMaker(TestSequence = Test_Sequence,ReferenceSequence = Full_16s_Data_RC[[.]]))
EndTime<-Sys.time()
EndTime-StartTime

Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
Extra_Column<-rep(x = c(1),times = 1055)
Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")

#Generate the plot here

CoP_Phylo<-cophenetic.phylo(x = vert.tree)
which(vert.tree$tip.label=="GCA_000818075.1")
Phylo_Distances<-CoP_Phylo[998,]
PD_Labels<-as.data.frame(cbind(vert.tree$tip.label,Phylo_Distances))
PD_Labels$Phylo_Distances<-as.numeric(PD_Labels$Phylo_Distances)

Multi_Mapping_Table_Plot_Test<-Multi_Mapping_Table_Plot
rownames(Multi_Mapping_Table_Plot_Test)=Multi_Mapping_Table_Plot_Test[,1]
rownames(PD_Labels)=PD_Labels[,1]
Scatter_Table<-merge(Multi_Mapping_Table_Plot_Test,PD_Labels,by=0,all=TRUE)


#minhashing scoring section
#directory should be Raw_FASTAs
#put in assembly file from the top to calculate
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Raw_FASTAs")

#before running this, change first assembly and 
for(i in list.files(path = "~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Raw_FASTAs")[1:1059])
{
  system(command = paste("~/Desktop/mash-OSX64-v2.2/mash dist GCA_901457615.1.txt ",i," > ~/Desktop/Examples_For_Meeting/5_MinHASH/",i,sep = ""))
}

MASH_List<-as.list(NULL)
setwd("~/Desktop/Examples_For_Meeting/5_MinHASH/")
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

Scatter_Table_Trim<-Scatter_Table2[,c(3,4,7,9)]
colnames(Scatter_Table_Trim)<-c("label","Score_16s","Phylo_Distance","MASH_Distance")
Scatter_Table_Trim$MASH_Distance<-as.numeric(Scatter_Table_Trim$MASH_Distance)
Scatter_Table_setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica")
#Finding String
Sequence_Of_Interest<-Full_16s_Data_RC[[which(list.files("Data_Files/16s_Sequences_BT")=="GCA_901457615.1.txt")-1]]
Test_Sequence<-Sequence_Of_Interest[2:8]

#fancier version
library("dada2", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
library("tidyverse", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

StartTime<-Sys.time()
#This line executes the above functions and applies them over the database of reference assembly 16-23s data, setting the range of reference assemblies to test with .x = 
Multi_Sum_Final<-map_int(.x = 1:1055,.f = ~Super_HAM_TableMaker(TestSequence = Test_Sequence,ReferenceSequence = Full_16s_Data_RC[[.]]))
EndTime<-Sys.time()
EndTime-StartTime

Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
Extra_Column<-rep(x = c(1),times = 1055)
Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")

#Generate the plot here

CoP_Phylo<-cophenetic.phylo(x = vert.tree)
which(vert.tree$tip.label=="GCA_901457615.1")
Phylo_Distances<-CoP_Phylo[290,]
PD_Labels<-as.data.frame(cbind(vert.tree$tip.label,Phylo_Distances))
PD_Labels$Phylo_Distances<-as.numeric(PD_Labels$Phylo_Distances)

Multi_Mapping_Table_Plot_Test<-Multi_Mapping_Table_Plot
rownames(Multi_Mapping_Table_Plot_Test)=Multi_Mapping_Table_Plot_Test[,1]
rownames(PD_Labels)=PD_Labels[,1]
Scatter_Table<-merge(Multi_Mapping_Table_Plot_Test,PD_Labels,by=0,all=TRUE)


#minhashing scoring section
#directory should be Raw_FASTAs
#put in assembly file from the top to calculate
setwd("~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Raw_FASTAs")

#before running this, change first assembly and 
for(i in list.files(path = "~/Desktop/GenomeTrakr_Analysis_SEnterica/Data_Files/Raw_FASTAs")[1:1059])
{
  system(command = paste("~/Desktop/mash-OSX64-v2.2/mash dist GCA_901457615.1.txt ",i," > ~/Desktop/Examples_For_Meeting/5_MinHASH/",i,sep = ""))
}

MASH_List<-as.list(NULL)
setwd("~/Desktop/Examples_For_Meeting/5_MinHASH/")
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

Scatter_Table_Trim<-Scatter_Table2[,c(3,4,7,9)]
colnames(Scatter_Table_Trim)<-c("label","Score_16s","Phylo_Distance","MASH_Distance")
Scatter_Table_Trim$MASH_Distance<-as.numeric(Scatter_Table_Trim$MASH_Distance)
Scatter_Table_GCA_901457615.1<-Scatter_Table_Trim<-Scatter_Table_Trim

##################################################################
######### Trying analysis again with e.coli 16-23s data##########
##################################################################

setwd("~/Desktop/GenomeTrakr_Analysis")
#Finding String
Sequence_Of_Interest<-Full_16_23s_Data_RC[[81]]
Test_Sequence<-Sequence_Of_Interest[2:8]

#fancier version
library("dada2", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")
library("tidyverse", lib.loc="/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

StartTime<-Sys.time()
#This line executes the above functions and applies them over the database of reference assembly 16-23s data, setting the range of reference assemblies to test with .x = 
Multi_Sum_Final<-map_int(.x = 1:1077,.f = ~Super_HAM_TableMaker(TestSequence = Test_Sequence,ReferenceSequence = Full_16_23s_Data_RC[[.]]))
EndTime<-Sys.time()
EndTime-StartTime

Names_List<-as.list(NULL)
for(i in 1:1077)
{
  Names_List[[i]]<-Full_16_23s_Data_RC[[i]][1]
}
Names_16_23s<-unlist(Names_List)

Multi_Mapping_Table<-as.data.frame(cbind(Names_16_23s,Multi_Sum_Final))
Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
Extra_Column<-rep(x = c(1),times = 1077)
Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")

#Generate the plot here
setwd("~/Desktop")
RAXML_Output<-readChar("RAxML_result.test.out",nchars = 2e7)
vert.tree_1623<-read.tree(text=RAXML_Output)

CoP_Phylo_1623<-cophenetic.phylo(x = vert.tree_1623)
which(vert.tree_1623$tip.label=="GCA_000803705.1")
Phylo_Distances<-CoP_Phylo_1623[1069,]
PD_Labels<-as.data.frame(cbind(vert.tree_1623$tip.label,Phylo_Distances))
PD_Labels$Phylo_Distances<-as.numeric(PD_Labels$Phylo_Distances)

Multi_Mapping_Table_Plot_Test<-Multi_Mapping_Table_Plot
rownames(Multi_Mapping_Table_Plot_Test)=Multi_Mapping_Table_Plot_Test[,1]
rownames(PD_Labels)=PD_Labels[,1]
Scatter_Table<-merge(Multi_Mapping_Table_Plot_Test,PD_Labels,by=0,all=TRUE)


#minhashing scoring section
#directory should be Raw_FASTAs
#put in assembly file from the top to calculate
setwd("~/Desktop/GenomeTrakr_Analysis/Data_Files/Raw_FASTAs")

#before running this, change first assembly and 
for(i in list.files(path = "~/Desktop/GenomeTrakr_Analysis/Data_Files/Raw_FASTAs")[1:1173])
{
  system(command = paste("~/Desktop/mash-OSX64-v2.2/mash dist GCA_000803705.1.txt ",i," > ~/Desktop/Examples_For_Meeting/7_MinHASH/",i,sep = ""))
}

MASH_List<-as.list(NULL)
setwd("~/Desktop/Examples_For_Meeting/7_MinHASH/")
for(j in 1077:1173)
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

Scatter_Table_Trim<-Scatter_Table2[,c(3,4,7,9)]
colnames(Scatter_Table_Trim)<-c("label","Score_16s","Phylo_Distance","MASH_Distance")
Scatter_Table_Trim$MASH_Distance<-as.numeric(Scatter_Table_Trim$MASH_Distance)
Scatter_Table_GCA_000803705.1<-Scatter_Table_Trim<-Scatter_Table_Trim

ggplot(Scatter_Table_GCA_000280315.2, aes(x=Score_16s, y=Phylo_Distance)) +
  geom_point(size=2, shape=23)+
  xlim(0,500)

ggplot(Scatter_Table_Trim, aes(x=MASH_Distance, y=Phylo_Distance)) +
  geom_point(size=2, shape=23)

########################################################
########################################################
########################################################
#Median Experiment

StartTime<-Sys.time()
#This line executes the above functions and applies them over the database of reference assembly 16-23s data, setting the range of reference assemblies to test with .x = 
Multi_Median_Final<-map_int(.x = 1:1055,.f = ~Super_HAM_TableMaker_Testing(TestSequence = Test_Sequence,ReferenceSequence = Full_16s_Data_RC[[.]]))
EndTime<-Sys.time()
EndTime-StartTime

Min_Median_Comparison<-as.data.frame(cbind(Multi_Sum_Final,Multi_Median_Final))
Min_Median_Comparison$Multi_Sum_Final<-as.numeric(Min_Median_Comparison$Multi_Sum_Final)
Min_Median_Comparison$Multi_Median_Final<-as.numeric(Min_Median_Comparison$Multi_Median_Final)
rownames(Min_Median_Comparison)<-Names_16_23s

ggplot(Min_Median_Comparison, aes(x=Multi_Sum_Final, y=Multi_Median_Final)) +
  geom_point(size=2, shape=23)+
  xlim(0,400)

#GCA_006514375.1 vs. GCA_000803705.1  (845 vs. 81)
Full_16_23s_Data_RC[[845]][1]
nchar(Full_16_23s_Data_RC[[81]][1])

#for final mapping
## Getting the nodes of interest
which(vert.tree$tip.label=="GCA_009650355.1")
which(vert.tree$tip.label=="GCA_009650395.1")
which(vert.tree$tip.label=="GCA_009650375.1")
which(vert.tree$tip.label=="GCA_002504125.1")

which(vert.tree$tip.label=="GCA_001975645.1")
which(vert.tree$tip.label=="GCA_002234795.1")

nodes_of_interest <- c(702)

## Getting the edges connecting to the nodes
edge_parent <- rownames(edge_table)[edge_table[,1] %in% nodes_of_interest]
edge_child  <- rownames(edge_table)[edge_table[,2] %in% nodes_of_interest]

which(vert.tree$tip.label=="GCA_009728795.1") #449
which(vert.tree$tip.label=="GCA_003030125.1") #448

which(vert.tree$tip.label=="GCA_008807355.1") #447

which(vert.tree$tip.label=="GCA_009650355.1") #702

Phylo_Distances<-CoP_Phylo[1045,]
Phylo_Distances[1002]


which(vert.tree$tip.label=="GCA_001104165.2") #1002
which(vert.tree$tip.label=="GCA_003718395.1") #46
which(vert.tree$tip.label=="GCA_000007545.1") #1045


ggplot(Scatter_Table, aes(x=value, y=Phylo_Distances)) +
  geom_point(size=2, shape=23)+
  xlim(0,250)+
  ylim(0,0.6)
