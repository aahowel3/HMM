# HMM for identifying excision boundaries and estimating ploidy

sequencing data from /work/aahowel3/Ancestor 1-20 A and B 
all scripts run in hines /work/aahowel3/HMM

#Identifying IESs in pooled GE sequencing data 
HMM_dataprocessing.sh - aligns Ancestor R1 and R2 GEs one at a time to the reference rather than concatenating them and aligning all at once
uses samtools depth to create coverage file - this is the emissions data fed into the HMM
#idiot the output bam, sam, and coverage folders output to the main /work/aahowel3/Ancestor folder not into the HMM folder, why did you do that 

spits out AncestorsGE_tomic_coverage.txt in /work/aahowel3/HMM - coverage depth every position in all chromosomes 
R script stan_hmm_nb.R use the STAN package to take AncestorsGE_tomic_coverage.txt input and estimate parameters and output the viterbi assignment for each position. 
model type (NB v. Poisson), number of states and number of BW iterations manually changed within the script 
must be run in tmux window and within conda enviornment where you have STAN installed 
install STAN using 
source("https://bioconductor.org/biocLite.R")
biocLite("STAN")

the output of stan_hmm_nb.R is viterbidata_nb_25_2state.txt and viterbidata_nb_25_3state.txt 
the folders /viterbidata_nb_25_2state and /viterbidata_nb_25_3state contain the respective main files broken up by chromosome 

modified stan_hmm_nb.R so that it also prints the parameter output to a file instead of the command line - parameter_experimentvariables.txt

#Check if the viterbi algo predits KNOWN IESs 
using the viterbi data per chromosome in either /viterbidata_nb_25_2state or /viterbidata_nb_25_3state, loop the chrX_viterbi file simultaneously with the chrX_IES_in_mic.tsv file in flowsortdata/retention_scores using the scripts validate_IESpositions.R and validate_IESpositions.sh 

command line is run as bash validate_IESpositions.sh /viterbidata_nb_25_2state/*viterbidata* inside the viterbidata_nb_25_2state_validated folder which will loops chrx_viterbi and chrX_IES files through the R script validate_IESpositions.R
3 state version is run as bash validate_IESpositions.sh /viterbidata_nb_25_3state/*viterbidata* inside the viterbidata_nb_25_3state_validated folder

#Identifying IESs in individual GE sequencing data 

in the /work/aahowel3/HMM/indiviudal_GEs/coverage_files folder the script coverage.sh generates indiviudal coverage.txt files from the bams generated in the /work/aahowel3/Ancestor folder 
