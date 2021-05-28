Sliding_Hit_Finder_Hits<-function(MappingData,Original)
{
  if(sort(MappingData)[2] <= 9)
  {
    TestSet_UnderCutoff<-which(MappingData<=9)
    TipsOI<-GFFs_Reference[TestSet_UnderCutoff]
    Remove_Index<-which(TipsOI == Original)
    FinalTipsOI<-TipsOI[-Remove_Index]
    return(FinalTipsOI)
  }
  
  if(sort(MappingData)[2] <= 14)
  {
    TestSet_UnderCutoff<-which(MappingData<=14)
    TipsOI<-GFFs_Reference[TestSet_UnderCutoff]
    Remove_Index<-which(TipsOI == Original)
    FinalTipsOI<-TipsOI[-Remove_Index]
    return(FinalTipsOI)
  }
  
  if(sort(MappingData)[2] <= 19)
  {
    TestSet_UnderCutoff<-which(MappingData<=19)
    TipsOI<-GFFs_Reference[TestSet_UnderCutoff]
    Remove_Index<-which(TipsOI == Original)
    FinalTipsOI<-TipsOI[-Remove_Index]
    return(FinalTipsOI)
  }
  
  if(sort(MappingData)[2] <= 24)
  {
    TestSet_UnderCutoff<-which(MappingData<=24)
    TipsOI<-GFFs_Reference[TestSet_UnderCutoff]
    Remove_Index<-which(TipsOI == Original)
    FinalTipsOI<-TipsOI[-Remove_Index]
    return(FinalTipsOI)
  }
  
  if(sort(MappingData)[2] <= 29)
  {
    TestSet_UnderCutoff<-which(MappingData<=29)
    TipsOI<-GFFs_Reference[TestSet_UnderCutoff]
    Remove_Index<-which(TipsOI == Original)
    FinalTipsOI<-TipsOI[-Remove_Index]
    return(FinalTipsOI)
  }
}