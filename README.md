# HMM for identifying excision boundaries and estimating ploidy

sequencing data from /work/aahowel3/Ancestor GES 1-20 A and B 
combined L001-L004 Ancestor files generated using Joelle's script
https://github.com/joellejohnson/Cartwrightlabprojects/blob/master/READme%20Anc%20.pdf 
all scripts run in hines /work/aahowel3/HMM

#Identifying IESs in pooled GE sequencing data 
HMM_dataprocessing.sh - aligns Ancestor R1 and R2 GEs one at a time to the reference rather than concatenating them and aligning all at once
the samtools rmdup processing you did probably is fine but in reality I dont think you needed to do this 
uses samtools depth to create coverage file - this is the emissions data fed into the HMM
#idiot the output bam, sam, and coverage folders output to the main /work/aahowel3/Ancestor folder not into the HMM folder, why did you do that 

spits out AncestorsGE_tomic_coverage.txt in /work/aahowel3/HMM - coverage depth every position in all chromosomes 
######AncestorGE_tomic_coverage.txt is now for coverage data INCLUDING dup reads - cut out rmdup line of code in samtools file
######AncestorGE_tomic_coverage_2.txt is now for coverage data generated using the samtools rmdup command 
#same with the folders coverage and coverage2
R script stan_hmm_nb.R use the STAN package to take AncestorsGE_tomic_coverage.txt input and estimate parameters and output the viterbi assignment for each position. 
model type (NB v. Poisson), number of states and number of BW iterations manually changed within the script 
must be run in tmux window and within conda enviornment where you have STAN installed 
install STAN using 
source("https://bioconductor.org/biocLite.R")
biocLite("STAN")

the output of stan_hmm_nb.R is viterbidata_nb_25_2state.txt and viterbidata_nb_25_3state.txt (you can change parameters - model, bw iterations, # of states by changing variable header section of stan_hmm_nb.R) 

modified stan_hmm_nb.R so that it also prints the parameter output to a file instead of the command line - parameter_experimentvariables.txt


#Check if the viterbi algo predits KNOWN IESs 
stan_hmm_nb.R gives you model parameters - you can check if theyre any good by seeing how they assign states within known IES regions 

validate_IESpositions.sh uses the output of stan_hmm.R to create folders of the output broken up by chromosome (ex. /viterbidata_nb_25_2state or /viterbidata_nb_25_3state)

validate_IESpositions.sh then loops the chrX_viterbi file simultaneously with the chrX_IES_in_mic.tsv file in flowsortdata/retention_scores using the script validate_IESpositions.R 
##IMPORTANT## to check which viterbi state (1,2,3 ect) stan_HMM assigned to each group (IES v. Mac v. zero coverage) - check the parameters.txt file - modify this appropriately in the validate_IESpositions.R file. 

command line is run as bash validate_IESpositions.sh (no input) in the main HMM folder where the stan_hmm.R output is 

#Visualizing results of validate_IESpositions.sh 
validate_IESpositions_graphing.R see how much of each IES the viterbi assignment captured (in % of correctly state assigned bps) - can compare to different model runs (2 v. 3 states, etc) 

#Identifying IESs in individual GE sequencing data 
in the /work/aahowel3/HMM/indiviudal_GEs/coverage_files folder the script coverage.sh generates indiviudal coverage.txt files from the bams generated in the /work/aahowel3/Ancestor folder 

#Working script to simulate observations in order to test the Viterbi and BW algorithims of the STAN package is simulate_hmm_onemorego.R on your local. 
