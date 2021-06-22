#Testing Section

D1<-unlist(Ancestors(vert.tree,node = 325))
D2<-unlist(Ancestors(vert.tree,node = 326))
D3<-unlist(Ancestors(vert.tree,node = 327))
D4<-unlist(Ancestors(vert.tree,node = 328))

Ancestor_List<-as.list(NULL)
Ancestor_List[[1]]<-D1
Ancestor_List[[2]]<-D2
Ancestor_List[[3]]<-D3
Ancestor_List[[4]]<-D4

Sum_List<-as.list(NULL)
for(i in 1:length(D2))
{
  Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(D2[i],.))
  Sum_List[[i]]<-Test_Vec
}

Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
Target_Value<-min(Result_Vec_NoNA)
Closest_Parent_Node<-D2[which(Result_Vec==Target_Value)]

#Testing on real data, need to make more automated
Testing_List
Multi_Sum_Final<-Testing_List[[31]]
Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
Extra_Column<-rep(x = c(1),times = 1055)
Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")

TableOI<-Multi_Mapping_Table_Plot[which(Multi_Mapping_Table_Plot$value<10),]
TipsOI<-TableOI[,1]
TipNumbersOI<-map_int(.x = TipsOI,.f = ~which(vert.tree$tip.label==.))
#Check here that the cluster makes sense
TipNumbersOI_Filtered<-TipNumbersOI
Ancestor_List<-map(.x = TipNumbersOI_Filtered,.f = ~unlist(Ancestors(vert.tree,node = .)))
Sum_List<-as.list(NULL)
for(i in 1:length(Ancestor_List[[1]]))
{
  Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(Ancestor_List[[1]][i],.))
  Sum_List[[i]]<-Test_Vec
}

Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
Target_Value<-min(Result_Vec_NoNA)
Closest_Parent_Node<-Ancestor_List[[1]][which(Result_Vec==Target_Value)]
Descendant_Values<-unlist(Descendants(vert.tree,Closest_Parent_Node))
Descendant_Assemblies<-vert.tree$tip.label[Descendant_Values]
NumberVector<-rep(1,length(TipsOI))
Target_Table1<-as.data.frame(cbind(TipsOI,NumberVector))
Target_Indices<-which(Target_Table1$TipsOI %in% Descendant_Assemblies)
Target_Table1[Target_Indices,2]<-2
colnames(Target_Table1)<-c("Assembly","Target")

#Trying to make it work with multiple groups
#First make table of edge lengths
Edge_Length_Table<-as.data.frame(cbind(vert.tree$edge,vert.tree$edge.length))
Test<-unlist(Ancestor_List)
Unique_Test<-unique(Test)
Test2<-Edge_Length_Table[which(Edge_Length_Table$V2 %in% Unique_Test),]

#trying to find edge length cutoff from tree
Big_Descendant_List<-map(.x = 1:2087,.f = ~Descendants(x = vert.tree,node = .))
Big_Descendant_EdgeLengths<-map(.x = Big_Descendant_List,.f = ~Edge_Length_Table[which(Edge_Length_Table$V2 %in% .[[1]]),])
Descendant_Sums<-unlist(map(.x = Big_Descendant_EdgeLengths,.f = ~sum(.[,3])))
hist(Descendant_Sums,breaks = 100,ylim = c(0,50),xlim = c(0,0.6))

#plotting: 1254 + 1269
Descendant_Values1<-unlist(Descendants(vert.tree,2043))
Descendant_Values2<-unlist(Descendants(vert.tree,1269))
Descendant_Values<-c(Descendant_Values1,Descendant_Values2)
Descendant_Assemblies<-vert.tree$tip.label[Descendant_Values]
NumberVector<-rep(1,length(Descendant_Assemblies))
Target_Table1<-as.data.frame(cbind(Descendant_Assemblies,NumberVector))
Target_Indices<-which(Target_Table1$Descendant_Assemblies %in% TipsOI)
Target_Table1[Target_Indices,2]<-2
colnames(Target_Table1)<-c("Assembly","Target")


#tree clustering section
Tree_Clusters<-read.csv(file = "treecluster2.txt",sep = "\t")

TableOI<-Tree_Clusters[which(Tree_Clusters$ClusterNumber==1),]
TipsOI<-TableOI[,1]
TipNumbersOI<-map_int(.x = TipsOI,.f = ~which(vert.tree$tip.label==.))
Ancestor_List<-map(.x = TipNumbersOI,.f = ~unlist(Ancestors(vert.tree,node = .)))
Sum_List<-as.list(NULL)
for(i in 1:length(Ancestor_List[[1]]))
{
  Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(Ancestor_List[[1]][i],.))
  Sum_List[[i]]<-Test_Vec
}

Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
Target_Value<-min(Result_Vec_NoNA)
Closest_Parent_Node<-Ancestor_List[[1]][which(Result_Vec==Target_Value)]
Closest_Parent_Node

#rooted tree?
TableOI<-Tree_Clusters[which(Tree_Clusters$ClusterNumber==23),]
TipsOI<-TableOI[,1]
TipNumbersOI<-map_int(.x = TipsOI,.f = ~which(rooted.vert.tree$tip.label==.))
Ancestor_List<-map(.x = TipNumbersOI,.f = ~unlist(Ancestors(rooted.vert.tree,node = .)))
Sum_List<-as.list(NULL)
for(i in 1:length(Ancestor_List[[1]]))
{
  Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(Ancestor_List[[1]][i],.))
  Sum_List[[i]]<-Test_Vec
}

Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
Target_Value<-min(Result_Vec_NoNA)
Closest_Parent_Node<-Ancestor_List[[1]][which(Result_Vec==Target_Value)]
Closest_Parent_Node

#plotting target family colors (salmonella)
p <- ggtree(rooted.vert.tree) + 
  xlim(.0000001, 5) +
  no_legend() +
  geom_point2(aes(subset=node==1048), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1052), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2068), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2051), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2050), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2043), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1123), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1126), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2042), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1377), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1155), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1240), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1255), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1334), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1345), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1351), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1354), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1366), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1046), color='darkgreen', size=5) 

p %<+% Target_Table1 + 
  geom_tiplab(aes(fill = factor(Target)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.10, "lines"), # amount of padding around the labels
              label.size = 0.01) 

#plotting target family colors (salmonella) - second part
p <- ggtree(rooted.vert.tree) + 
  xlim(.0000001, 5) +
  no_legend() +
  geom_point2(aes(subset=node==1048), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1052), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2068), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2051), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2050), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2043), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1123), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1126), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==2042), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1155), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1378), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1240), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1424), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1762), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1256), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1432), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1435), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1500), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1336), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1345), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1354), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1366), color='darkgreen', size=5) +
  geom_point2(aes(subset=node==1046), color='darkgreen', size=5) 

p %<+% Target_Table1 + 
  geom_tiplab(aes(fill = factor(Target)),
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.10, "lines"), # amount of padding around the labels
              label.size = 0.01) 


#want to look into test set to see results as a whole
#Cutoffs will be: 9, 14, 30
Big_Results_List<-as.list(NULL)
#start for loop here
Node_Distances<-dist.nodes(vert.tree)
i<-12

for(i in Scale_29)
{
  TestSet<-AssemblyMappingData[[i]]
  TestSet_UnderCutoff<-TestSet[which(TestSet<9)]
  Multi_Sum_Final<-AssemblyMappingData[[i]]
  Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
  Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
  Extra_Column<-rep(x = c(1),times = 1055)
  Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
  colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")
  TableOI<-Multi_Mapping_Table_Plot[which(Multi_Mapping_Table_Plot$value<9),]
  TipsOI<-TableOI[,1]
  TipNumbersOI<-map(.x = TipsOI,.f = ~which(vert.tree$tip.label==.))
  TipNumbersOI_Filtered<-TipNumbersOI[which(TipNumbersOI>=0)]
  Ancestor_List<-map(.x = TipNumbersOI_Filtered,.f = ~unlist(Ancestors(vert.tree,node = .)))
  Sum_List<-as.list(NULL)
  for(j in 1:length(Ancestor_List[[1]]))
  {
    Test_Vec<-map_int(.x = Ancestor_List,.f = ~match(Ancestor_List[[1]][j],.))
    Sum_List[[j]]<-Test_Vec
  }
  Result_Vec<-map_int(.x = Sum_List,.f = ~sum(.))
  Result_Vec_NoNA<-Result_Vec[!is.na(Result_Vec)]
  Target_Value<-min(Result_Vec_NoNA)
  Closest_Parent_Node<-Ancestor_List[[1]][which(Result_Vec==Target_Value)]
  Descendant_Values<-unlist(Descendants(vert.tree,Closest_Parent_Node))
  #Final_Data<-c(length(TestSet_UnderCutoff),length(Descendant_Values))
  Depth_Range<-Node_Distances[Closest_Parent_Node,Descendant_Values]
  Depth<-max(Depth_Range)
  Sum_Branch_Length<-sum(Depth_Range)
  Original_Query<-which(vert.tree$tip.label==Full_AssembliesOI[i])
  Original_Distance<-Node_Distances[Original_Query,unlist(TipNumbersOI)]
  Closest_OD<-sort(Original_Distance)[2]
  Final_Data<-c(length(TestSet_UnderCutoff),Depth,Sum_Branch_Length,Closest_OD)
  Big_Results_List[[i]]<-Final_Data
}

Cutoff9 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
colnames(Cutoff9) <- c("Hit_Number","Max_Depth","Sum_Branch_Lengths","Min_Distance_Original")
Cutoff9_NoHits<-Cutoff9[which(Cutoff9$Hit_Number==1),]
Cutoff9_Hits<-Cutoff9[which(Cutoff9$Hit_Number>1),]
length(which(Cutoff9_Hits$Min_Distance_Original<1e-4))/length(Cutoff9_Hits$Min_Distance_Original)

Cutoff14 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
colnames(Cutoff14) <- c("Hit_Number","Max_Depth","Sum_Branch_Lengths","Min_Distance_Original")
Cutoff14_NoHits<-Cutoff14[which(Cutoff14$Hit_Number==1),]
Cutoff14_Hits<-Cutoff14[which(Cutoff14$Hit_Number>1),]
length(which(Cutoff14_Hits$Min_Distance_Original<1e-4))/length(Cutoff14_Hits$Min_Distance_Original)

Cutoff29 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
colnames(Cutoff29) <- c("Hit_Number","Max_Depth","Sum_Branch_Lengths","Min_Distance_Original")
Cutoff29_NoHits<-Cutoff29[which(Cutoff29$Hit_Number==1),]
Cutoff29_Hits<-Cutoff29[which(Cutoff29$Hit_Number>1),]
length(which(Cutoff29_Hits$Min_Distance_Original<1e-4))/length(Cutoff29_Hits$Min_Distance_Original)

CutoffScale <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
colnames(CutoffScale) <- c("Hit_Number","Max_Depth","Sum_Branch_Lengths","Min_Distance_Original")
CutoffScale_NoHits<-CutoffScale[which(CutoffScale$Hit_Number==1),]
CutoffScale_Hits<-CutoffScale[which(CutoffScale$Hit_Number>1),]
length(which(CutoffScale_Hits$Min_Distance_Original<1e-4))/length(CutoffScale_Hits$Min_Distance_Original)


hist(Cutoff9_Hits$Max_Depth,breaks=100)
FSR_Hits$Hit_Number<-as.numeric(FSR_Hits$Hit_Number)
hist(FSR_Hits$Max_Depth,breaks=100)

#sliding version
Rescale_9 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
Scale_14<-which(Rescale_9$V1==1)
Rescale_14 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
Scale_19<-which(Rescale_14$V1==1)
Rescale_19 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
Scale_24<-which(Rescale_19$V1==1)
Rescale_24 <- as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
Scale_29<-which(Rescale_24$V1==1)
Final_Scaling_Result<-as.data.frame(matrix(unlist(Big_Results_List), ncol = 4, byrow = TRUE))
colnames(Final_Scaling_Result) <- c("Hit_Number","Max_Depth","Sum_Branch_Lengths","Min_Distance_Original")
FSR_NoHits<-Final_Scaling_Result[which(Final_Scaling_Result$Hit_Number==1),]
FSR_Hits<-Final_Scaling_Result[which(Final_Scaling_Result$Hit_Number>1),]
length(which(FSR_Hits$Min_Distance_Original<1e-4))/length(FSR_Hits$Min_Distance_Original)


#is the assembly finding itself
Percentage_List<-as.list(NULL)
i<-20
for(i in Scale_29)
{
  TargetAssembly<-Full_AssembliesOI[i]
  TestSet<-AssemblyMappingData[[i]]
  TestSet_UnderCutoff<-TestSet[which(TestSet<29)]
  Multi_Sum_Final<-AssemblyMappingData[[i]]
  Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
  Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
  Extra_Column<-rep(x = c(1),times = 1055)
  Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
  colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")
  TableOI<-Multi_Mapping_Table_Plot[which(Multi_Mapping_Table_Plot$value<29),]
  TipsOI<-TableOI[,1]
  Target_Serovars<-TaxID_Table3[which(TaxID_Table3$Assembly %in% TipsOI),]
  Original_Serovar<-TaxID_Table3[which(TaxID_Table3$Assembly %in% TargetAssembly),]
  SerovarName<-Original_Serovar[1,2]
  Percentage_List[[i]]<-length(which(Target_Serovars$Serovar==SerovarName))/length(Target_Serovars$Serovar)
}

Percentage_List9<-Percentage_List
Percentage_List9_Trim<-Percentage_List9[which(Cutoff9$Hit_Number>1)]
Percent_Serovars_Matching_Cutoff9<-unlist(Percentage_List9_Trim)
hist(Percent_Serovars_Matching_Cutoff9,breaks = 40)

Percentage_List14<-Percentage_List
Percentage_List14_Trim<-Percentage_List14[which(Cutoff14$Hit_Number>1)]
Percent_Serovars_Matching_Cutoff14<-unlist(Percentage_List14_Trim)
hist(Percent_Serovars_Matching_Cutoff14,breaks = 20)

Percentage_List30<-Percentage_List
Percentage_List30_Trim<-Percentage_List30[which(Cutoff30$Hit_Number>1)]
Percent_Serovars_Matching_Cutoff30<-unlist(Percentage_List30_Trim)
hist(Percent_Serovars_Matching_Cutoff30,breaks = 20)

Percentage_ListSc<-Percentage_List
Percentage_ListSc_Trim<-Percentage_ListSc[which(CutoffScale$Hit_Number>1)]
Percent_Serovars_Matching_CutoffSc<-unlist(Percentage_ListSc)
hist(Percent_Serovars_Matching_CutoffSc,breaks = 40)
