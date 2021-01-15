#BiocManager::install("STAN")
library(STAN)

#CHANGE TO TEST DIFFERENT STATE NUMBERS, DISTRIBUTIONS, BW ITERATIONS 
#variables parts of script 
number_of_states=3
distribution="NegativeBinomial"   
bw_iterations=10
experiment="_nb_10_3state"



#run script 
#read in data from csv
data=read.csv("AncestorsGE_tomic_coverage.txt", sep="\t",header=FALSE)
colnames(data)=c("chromosome", "bp_pos", "coverage")

#data structure this package takes is a list of matrices 
#force dataframe to matric and then convert integers to doubles
coverage_matrix=matrix(data$coverage, nrow=length(data$coverage),ncol=1)
mode(coverage_matrix) <- "double"
#create list of matrices 
coverage_list=list(coverage_matrix)

#Negative Binomial 
#intialize model 
hmm_ex = initHMM(coverage_list, nStates=number_of_states, method=distribution)
#E-M algo to find best parameters 
hmm_fitted = fitHMM(coverage_list, hmm_ex, maxIters = bw_iterations)
#state annotation per bp 
nbViterbi=getViterbi(hmm_fitted, coverage_list)


#get transistion matrix
trans=Transitions(hmm_fitted)

#get start probabilities 
inits=InitProb(hmm_fitted)

#get emission parameters for each state 
ems=EmissionParams(hmm_fitted)

#likelihood of model 
loglik=getLogLik(hmm_fitted, coverage_list)

parameterfilename=paste("parameter", experiment, ".txt",sep="")
sink(parameterfilename)
print("Emission Parameters")
print(ems)
print("Transistion Probabilities")
print(trans)
print("Initial Probabilities")
print(inits)
print("Logliklihood")
print(loglik)
sink()

#convert vitberi sequence to a column in R
#nbViterbiList <- unlist(unlist(nbViterbi))
#data$viterbi = nbViterbiList
#output it and then download it back to R later to create graph - illustrate state differences between models\
#pick random 25,000 region in middle of chromosome 2 with at least 2-3 IESs 
#middle of chromosome will likely look best - less assembly junk
#outputfile_name=paste("viterbidata", experiment, ".txt",sep="") 
#write.table(data, file=outputfile_name,sep="\t",row.names=TRUE, col.names=TRUE,quote=FALSE)





