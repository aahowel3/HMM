args = commandArgs(trailingOnly=TRUE)


####CHECK IN PARAMETERS.TXT WHICH VITERBI STATE ASSIGNMENT (1,2,3) WAS FOR THE IES GROUP - SHOULD HAVE LOW BUT NOT ZERO MEAN
iesvit=2

library(tidyverse)
#does this work with the pre-split IES file 
data3_chromosomes=read.csv(args[1],sep="\t")
  
#import viterbi data file 
datav=read.csv(args[2],sep="\t",header=FALSE)

test_subset=data3_chromosomes

#for each row (x) in the IES list dataset, subset the viterbi dataset between x row column[11] (IES in chr start) and column[12] (IES in chr end) of the IES list dataset
#Gives you a list of dataframes, named by the row index of the test_subset and the V4 column is the state position 
f = function(x){  
  subset(datav, ((as.numeric(x[11]) < datav$V3) & 
                   (datav$V3 < as.numeric(x[12]))))  
  
}

list_pulledout_viterbis=apply(test_subset, 1, f)

#rename list of dataframes to IES name 
#IES name alone is NOT a unique identifier bc supercontigs map to multiple locations in the chromosome assembly 
#IES-chr is ALSO not a unique identifer bc some contigs are used to assemble the same chromosme multiple times 
#IES-chr-IES_in_chr_startposition should be unique 
#create new column within the IES data list and rename list of dataframes w that 
#FUCKMYLIFE the IES identifers arent even unique in data1 they use the same goddamn identifier for 2 IESs in the same supercontig 
#but theres only 2 and they do start in diff positions in the same chr so the identifer IES_chr_startpos is still unique
#need to rename >names in the IES fasta 
test_subset$uniqueidentifier = paste0(test_subset$IES_ID, "-", test_subset$chromosome_name, "-IES_in_chr_startpos_", test_subset$IES_in_chromosome_start)
rename=as.vector(test_subset$uniqueidentifier)
names(list_pulledout_viterbis) = rename

#now calcualte % of 2s in each pulled out region
percentages <- lapply(list_pulledout_viterbis, function(x){
  twos=sum(x$V5==iesvit)
  totallength=length(x$V5)
  twos/totallength*100
})


x=as.data.frame(percentages)
y=gather(x)

print(y)


