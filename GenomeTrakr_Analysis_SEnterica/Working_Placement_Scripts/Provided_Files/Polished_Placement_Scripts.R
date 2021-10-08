#It is reccomended to make a new directory for each analysis in order to organize the files for epa-ng 
#and the resulting outputs. All of the code can be run in the same working directory for consistency 
#and consolidation of all data. 

#Required package loading
library(tidyverse)
library(phangorn)
library(ggtree)
library(dada2)
library(phylotools)
library(treeio)
library(Biostrings)

source("Placement_Functions.R")
load("Placement_Required_Files.Rdata")

#The first step is to prepare an aligned version of the input sequences. These will be done
#with the MAFFT tool. Installation instructions for Mac OS X and Windows can be found at
#https://mafft.cbrc.jp/alignment/software/source.html 

#After installation, the following MAFFT command will be used in order to align the set of query sequences to 
#our pre-aligned 16s reference sequence database.

#Raw_Query.fasta should be a fasta file containing exactly 7 Salmonella 16s ASVs, in a non-aligned format. These
#sequences should contain only nucleotide characters (A,G,C,T) and no dashes (-). 

system(command = paste("mafft --add Raw_Query.fasta --keeplength Full_Alignment.txt > Full_Alignment_With_Query.fasta"))

#After this alignment command is executed, the aligned query will be the last 7 sequences in the completed file 
#(Full_Alignment_With_Query.fasta). These 7 sequences need to be copied into a new text file separate from the larger
#alignment file and renamed to "Aligned_Query.fasta". This file should contain only the 7 aligned sequences from the 
#previous MAFFT alignment step and no other sequences for the following analysis steps. 

#This first command will prepare the data for you given paths to the raw query FASTA and the aligned query FASTA. This
#command will output two new files into your working directory which are needed for the EPA-NG command step. In 
#this function, set the PathToRawQuery and PathToAlignedQuery to the appropriate paths for the query data file and
#the aligned query file prepared with the previous command. 

Data_Preparation(PathToRawQuery = "Raw_Query.fasta",PathToAlignedQuery = "Aligned_Query.fasta",ReferenceFile = Full_16s_Data_Trim,Aligned_ReferenceFile = Full_16s_Data_Aligned)

#EPA-NG section: this command outside of R will require that all of the input files are present in the same working
#directory. If consistently working in one directory, the Data_Preparation function will have already loaded the 
#necessary files into the correct location. 

#EPA-NG (https://github.com/Pbdas/epa-ng) can be installed by two convenient methods:

#With Conda:
#conda install -c bioconda epa-ng

#With Homebrew:
#brew install brewsci/bio/epa-ng

#After installation of EPA-NG, run the following command in the same working directory to generate the placement results.
#data file. 
system(command = paste("epa-ng --ref-msa concatref.fasta --tree full.tree --query query.fasta --model info.raxml.bestModel --filter-max 100"))

#This function requires that the working directory is set to the location of the EPA-NG jplace output file which
#it will read and calculate an MRCA from the placement results. The MRCA reported here is used in further calculations
#to predict the serovar and plot the actual phylogeny for visualization. 
Clade_MRCA<-Clade_Hit_Finder_Pendant_Final(Pendant_Multi = 1.5,Tree = vert.tree.correct)

#Sero_Table_Final will contain a shortened list of serovar results, and Serovar_Result will report
#the predicted serovar if there is a dominant serovar representing >30% of the results observed. This function
#will output a table containing the following information for each resulting serovar: 1. The number of matches 
#in the final resulting clade 2. The fraction of total matches which that serovar represents 3. The maximum depth,
#or the distance from the MRCA to the farthest hit in the clade (calculated for each serovar) 4. The sum of branch 
#lengths, or the total sum of branch lengths from each hit to the MRCA (calculated for each serovar)
Sero_Table_Final<-Placement_Results_Output(MRCA = Clade_MRCA,Tree = vert.tree.correct)
Serovar_Result<-Sero_Result(Sero_Table_Final)
#Alternatively, a full serovar report can be obtained by running the alternative function below. Placement_Results_Output
#will trim any results with below 5% representation in the results table, in order to see all results including
#low percentage results, the Placement_Results_Output_Full function will need to be used. 
Sero_Table_Full<-Placement_Results_Output_Full(MRCA = Clade_MRCA,Tree = vert.tree.correct)

#The following function will generate a circular phylogeny plot for easier visualization and color the path from 
#the clade MRCA to the resulting hits in red. After running the phylogeny plotting function below, plot the 
#Phylogeny_Plot object to generate the actual plot. Because of the size of the actual phylogeny, it is reccomended 
#to save the resulting plot as a PNG in dimensions of at least 15,000 x 15,000 pixels. An outside viewer should be
#used to open the PNG and zoom in to analyze the resulting phylogeny figure. Larger dimensions will increase clarity 
#and allow better zooming at the cost of larger image size. 
Phylogeny_Plot<-Phylogeny_Plotting(MRCA = Clade_MRCA,InitialTable = ColoringTable,Tree = vert.tree.correct)

#The plotting function from ggtree sometimes throws an error which reads "Error in UseMethod("depth") : no applicable
#method for 'depth' applied to an object of class "NULL"". Re-running the exact same code does not cause an error
#the second time. 

