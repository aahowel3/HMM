#BiocManager::install("STAN")
library(STAN)
#these numbers are drawn from /work/aahowel3/HMM/parameter_nb_25_2state.txt file - which is generated from the stan_hmm.R 
States=c("MAC","IES")

x=8.35468e-146
startProbs = c(1-x,x)

transitions=c(0.9995179399, 0.0004820599, 0.0009551858, 0.9990448144)
transitions = matrix(transitions,nrow=2, dimnames =list(c("MAC","IES"), c("MAC","IES"))) 

#####################################################################
#this parts not terribly important 
####length of simulations estimations come from average scaffold length from Hamilton 2016 
#setwd("~/Documents/flowsortdata/retention_scores")
#test STAN hmm parameter estimation on scaffold length simualtion 
#estimage average MIC scaffold lenght from Hamilton 2016 
contigs=read.csv("contig_to_chromosome.csv")
#remove assigned scaffolds 
contigs <-contigs[!(contigs$Chromosome.superassembly=="unassigned"),]
contigs$length=abs(contigs$Start.coord. - contigs$End.coord.)
mean(contigs$length)
#average scaffold length = 192994.6bp 
#####################################################################

#length of simuation - start with chromosome length go to whole genome 
length=1000

  states   = c()
  emission = c()
  states   = c(states, sample(States,1,prob=startProbs))
#it starts at pos 2 instead of pos 1 bc state 1 is det by startProbs not the transProbs   
  for(i in 2:length)
  {
    #sample 1 thing from the States box based on the transistion probabilities of those states 
    state  = sample(States, 1, prob=transitions[states[i-1],])
      states = c(states, state)
  }
  #
  for(i in states)
  {
    #arguments in nbinom are # of experiments (outputs), mean, shape parameter
    #have to put them INSIDE the loop so it will generate a different number each time
    MAC_nbinom=rnbinom(1,m=115.323730,s=1.6469341)
    IES_nbinom=rnbinom(1,m=3.819932,s=0.7509493)
    emi=ifelse(i=="MAC",MAC_nbinom,IES_nbinom)
    emission = c(emission, emi)
  }

  states
  emission

#convert     
sim_emissions=data.frame(emission)  
sim_emissions$bp_pos=seq.int(nrow(sim_emissions))
sim_emissions$chromosome="scaffold"

sim_em=sim_emissions[,c("chromosome","bp_pos","emission")]
colnames(sim_em)=c("chromosome", "bp_pos", "coverage")

#####
#run your simulated emissions through the STAN parameter estimators - should get out what you put in god willing 
#####
#data structure this package takes is a list of matrices 
#force dataframe to matric and then convert integers to doubles
coverage_matrix=matrix(sim_em$coverage, nrow=length(sim_em$coverage),ncol=1)

mode(coverage_matrix) <- "double"
#create list of matrices 
coverage_list=list(coverage_matrix)

#Negative Binomial 
#intialize model 
hmm_ex = initHMM(coverage_list, nStates=2, method="NegativeBinomial")
#E-M algo to find best parameters 
hmm_fitted = fitHMM(coverage_list, hmm_ex, maxIters = 50)
#state annotation per bp 
nbViterbi=getViterbi(hmm_fitted, coverage_list)

#convert vitberi sequence to a column in R
nbViterbiList <- unlist(unlist(nbViterbi))
sim_em$viterbi = nbViterbiList

#get transistion matrix
Transitions(hmm_fitted)

#get start probabilities 
InitProb(hmm_fitted)

#get emission parameters for each state 
EmissionParams(hmm_fitted)

#likelihood of model 
getLogLik(hmm_fitted, coverage_list)




