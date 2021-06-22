#want to look into test set to see results as a whole
#Cutoffs will be: 9, 14, 29
Big_Results_List<-as.list(NULL)
#start for loop here
Node_Distances<-dist.nodes(vert.tree)

for(i in Scale_19)
{
  TestSet<-Median_List[[i]]
  TestSet_UnderCutoff<-TestSet[which(TestSet<19)]
  Multi_Sum_Final<-Median_List[[i]]
  Multi_Mapping_Table<-as.data.frame(cbind(GFFs_To_Process3,Multi_Sum_Final))
  Multi_Mapping_Table$Multi_Sum_Final<-as.numeric(Multi_Mapping_Table$Multi_Sum_Final)
  Extra_Column<-rep(x = c(1),times = 1055)
  Multi_Mapping_Table_Plot<-cbind(Multi_Mapping_Table,Extra_Column)
  colnames(Multi_Mapping_Table_Plot)<-c("label","value","category")
  TableOI<-Multi_Mapping_Table_Plot[which(Multi_Mapping_Table_Plot$value<19),]
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
  Original_Distance<-Node_Distances[Original_Query,Descendant_Values]
  Closest_OD<-sort(Original_Distance)[2]
  Final_Data<-c(length(TestSet_UnderCutoff),Depth,Sum_Branch_Length,Closest_OD)
  Big_Results_List[[i]]<-Final_Data
}

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
FSR_Median_NoHits<-Final_Scaling_Result[which(Final_Scaling_Result$Hit_Number==1),]
FSR_Median_Hits<-Final_Scaling_Result[which(Final_Scaling_Result$Hit_Number>1),]
length(which(FSR_Hits$Min_Distance_Original<1e-4))/length(FSR_Hits$Min_Distance_Original)

FSR_Hits$Hit_Number<-as.numeric(FSR_Hits$Hit_Number)
hist(FSR_Hits$Max_Depth,breaks=100)

#is the assembly finding itself
Percentage_List<-as.list(NULL)
i<-20
for(i in Scale_29)
{
  TargetAssembly<-Full_AssembliesOI[i]
  TestSet<-Median_List[[i]]
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

Percentage_ListSc<-Percentage_List
Percentage_ListSc_Trim<-Percentage_ListSc[which(CutoffScale$Hit_Number>1)]
Percent_Serovars_Matching_CutoffSc<-unlist(Percentage_ListSc_Trim)
hist(Percent_Serovars_Matching_CutoffSc,breaks = 40)












#these scrips are the extras for doing calculations using medians instead of sums
Super_Ham_TableMaker_Median<-function(TestSequence,ReferenceSequence)
{
  Row1<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[2], vec=TRUE, band=1000))
  Row2<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[3], vec=TRUE, band=1000))
  Row3<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[4], vec=TRUE, band=1000))
  Row4<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[5], vec=TRUE, band=1000))
  Row5<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[6], vec=TRUE, band=1000))
  Row6<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[7], vec=TRUE, band=1000))
  Row7<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[8], vec=TRUE, band=1000))
  Hamming_Table_Result<-rbind(Row1,Row2,Row3,Row4,Row5,Row6,Row7)
  BestMatch<-min(map_int(.x = 1:5040,.f = ~Super_HAM_Combinator_Median(Hamming_Table = Hamming_Table_Result,Combination = as.integer(CombinationTable[.,]))))
  return(BestMatch)
}

Super_HAM_Combinator_Median<-function(Hamming_Table,Combination)
{
  Sum_Output<-median(c(Hamming_Table[1,Combination[1]],Hamming_Table[2,Combination[2]],Hamming_Table[3,Combination[3]],Hamming_Table[4,Combination[4]],Hamming_Table[5,Combination[5]],Hamming_Table[6,Combination[6]],Hamming_Table[7,Combination[7]]))
  return(Sum_Output)
}

#nicer plots for lab meeting
CutoffScale[which(CutoffScale$Hit_Number==(-5)),2]<-(-0.03)

ggplot(CutoffScale, aes(x=Hit_Number)) + 
  geom_histogram(binwidth=1,color="black",fill="white")+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14,16),limits = c(0,16),expand = c(0,0))+
  scale_x_continuous(breaks = c(0,25,50,75,100,125,150,175,200,225,250))+
  ylab("Frequency")+
  xlab("Hit Number")+
  ggtitle("Number of Hits Using Scaling Threshold")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(CutoffScale, aes(x=Max_Depth)) + 
  geom_histogram(binwidth=0.005,color="black",fill="white")+
  scale_y_continuous(expand = c(0,0),limits = c(0,40))+
  scale_x_continuous()+
  ylab("Frequency")+
  xlab("Depth of Farthest Hit from MRCA (Edge Length)")+
  ggtitle("Maximum Depth Using Scaling Threshold")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(CutoffScale, aes(x=Max_Depth)) + 
  geom_histogram(binwidth=0.005,color="black",fill="white")+
  scale_y_continuous(expand = c(0,0),limits = c(0,40))+
  scale_x_continuous()+
  ylab("Frequency")+
  xlab("Depth of Farthest Hit from MRCA (Edge Length)")+
  ggtitle("Maximum Depth Using Scaling Threshold")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

CutoffScale<-cbind(CutoffScale,Percent_Serovars_Matching_CutoffSc)

ggplot(CutoffScale, aes(x=Percent_Serovars_Matching_CutoffSc)) + 
  geom_histogram(binwidth=0.01,color="black",fill="white")+
  scale_y_continuous(expand = c(0,0),limits = c(0,40),breaks=c(0,5,10,15,20,25,30,35,40))+
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
  ylab("Frequency")+
  xlab("Percentage Matching")+
  ggtitle("Percentage of Serovars Matching Original Query")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

ggplot(data=Cumulative_Plot, aes(x=V1, y=V2, group=1)) +
  geom_line(linetype = "dashed")+
  geom_point()+
  scale_y_continuous(limits = c(0,1),breaks=c(0,0.2,0.4,0.6,0.8,1))+
  ggtitle("Percent of Serovars Matching Position of Original Query")+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("-log10(Edge Length to Closest Hit)")+
  ylab("Percentage")

TESTSequence<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_012052445.1")]][2:8]
TESTSequence2<-Full_16s_Data_RC[[which(GFFs_To_Process3=="GCA_004919305.1")]][2:8]



