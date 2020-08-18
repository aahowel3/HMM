# HMM for identifying excision boundaries and estimating ploidy

sequencing data from /work/aahowel3/Ancestor 1-20 A and B 
all scripts run in hines /work/aahowel3/HMM

#Identifying IESs in pooled GE sequencing data 
HMM_dataprocessing.sh - aligns Ancestor R1 and R2 GEs one at a time to the reference rather than concatenating them and aligning all at once
uses samtools depth to create coverage file - this is the emissions data fed into the HMM
#idiot the output bam, sam, and coverage folders output to the main Ancestor folder not into the HMM folder, why did you do that 

spits out AncestorsGE_tomic_coverage.txt - coverage depth every position in all chromosomes 
R script stan_hmm.R and stan_hmm_nb.R use the STAN package to take AncestorsGE_tomic_coverage.txt input and estimate parameters and output the viterbi assignment for each position. 




