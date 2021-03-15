Super_HAM_Placer<-function(TestSequence,ReferenceSequence)
{
  Row1<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[2], vec=TRUE, band=1000))
  Row2<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[3], vec=TRUE, band=1000))
  Row3<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[4], vec=TRUE, band=1000))
  Row4<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[5], vec=TRUE, band=1000))
  Row5<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[6], vec=TRUE, band=1000))
  Row6<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[7], vec=TRUE, band=1000))
  Row7<-map_int(.x = 1:7,.f = ~nwhamming(s1 = TestSequence[.],s2 = ReferenceSequence[8], vec=TRUE, band=1000))
  Hamming_Table_Result<-rbind(Row1,Row2,Row3,Row4,Row5,Row6,Row7)
  Combination_Integers<-map_int(.x = 1:5040,.f = ~Super_HAM_Combinator(Hamming_Table = Hamming_Table_Result,Combination = as.integer(CombinationTable[.,])))
  Best_Order<-as.integer(CombinationTable[which(Combination_Integers==min(Combination_Integers))[1],])
  return(Best_Order)
}

Concatenate_By_Order<-function(Order,RefSequence)
{
  Sequence<-RefSequence[2:8]
  Concatenated_Output<-paste(Sequence[Order[1]],Sequence[Order[2]],Sequence[Order[3]],Sequence[Order[4]],Sequence[Order[5]],Sequence[Order[6]],Sequence[Order[7]],sep = "")
  return(Concatenated_Output)
}