# HMM for identifying excision boundaries and estimating ploidy
sequencing data from /work/aahowel3/Ancestor GES 1-20 A and B combined L001-L004 Ancestor files generated using Joelle's script
https://github.com/joellejohnson/Cartwrightlabprojects/blob/master/READme%20Anc%20.pdf 
all scripts run in hines /work/aahowel3/HMM

# READ FIRST - Issue with alignment files in /Ancestors
# All updated alignments are in /work/aahowel3/Alignments 
The ancestor GE mic and mac alignment files in /Ancestors are incomplete - 11_A in particular has an empty sam file. This is due to the reads being incorrectly ordered in the R1s and R2s when concatenating the Lane 1-4s. The order can be fixed with the script repairfqs.sh that uses the bbmap package.
In /Alignments there are 3 folders - fastqs, mic, and mac 
Script ancestors_aligntoX.sh in the mic and mac folder generate a sam and bam folder - in bam folder are final indexed and sorted bam files 

##############Scripts below here need to be updated to pull from new bamfiles in /work/aahowel3/Alignments 
# Identifying IESs in pooled GE sequencing data 
HMM_dataprocessing.sh - aligns Ancestor R1 and R2 GEs one at a time to the reference rather than concatenating them and aligning all at once
uses samtools depth to create coverage file - this is the emissions data fed into the HMM
The output bam, sam, and coverage folders output to the main /work/aahowel3/Ancestors folder 

Final output of HMM_dataprocessing.sh is AncestorsGE_tomic_coverage.txt in /work/aahowel3/HMM - coverage depth every position in all chromosomes 
AncestorGE_tomic_coverage.txt and coverage folder is now for coverage data INCLUDING dup reads - removed the rmdup line of code in samtools file
AncestorGE_tomic_coverage_2.txt and coverage2 is now for coverage data generated using the samtools rmdup command 
Did not physically remove or comment out the samtools rmdup line in the HMM_dataprocessing script - but by using the intermediate sorted.bam files instead of the rmdup_sorted.bam files to create AncestorGE_tomic_coverage.txt the effect is the same 

R script stan_hmm_nb.R use the STAN package to take AncestorsGE_tomic_coverage.txt input and estimate parameters and output the viterbi assignment for each position. 
model type (Negative Binomial v. Poisson), number of states and number of Baum-Welch iterations can changed within the script in the variable header section 
must be run in tmux window and within conda enviornment where you have STAN installed 
install STAN using 
source("https://bioconductor.org/biocLite.R")
biocLite("STAN")

the output of stan_hmm_nb.R is viterbidata_exerpimentalvariables_state.txt (you can change parameters - model, bw iterations, # of states by changing variable header section of stan_hmm_nb.R) 
modified stan_hmm_nb.R so that it also prints the parameter output to a file instead of the command line - parameter_experimentvariables.txt


# Check if the viterbi algorithim predits known (Hamilton et. al, 2016) IESs 
stan_hmm_nb.R gives you model parameters - you can check if the are accurate predictions by seeing how they assign states within known IES regions 

validate_IESpositions.sh uses the output of stan_hmm.R to create folders of the output broken up by chromosome (ex. /viterbidata_nb_25_2state or /viterbidata_nb_25_3state)

validate_IESpositions.sh then loops the chrX_viterbi file simultaneously with the chrX_IES_in_mic.tsv file in flowsortdata/retention_scores using the script validate_IESpositions.R 
IMPORTANT to check which viterbi state (1,2,3 ect) stan_HMM assigned to each group (IES v. Mac v. zero coverage) - check the parameters.txt file - modify this appropriately in the validate_IESpositions.R file. 

command line is run as bash validate_IESpositions.sh (no input argument) in the main HMM folder where the stan_hmm.R output is 

# Visualizing results of validate_IESpositions.sh 
Concatenate all the chrX tsv files together per folder for the input into script Validate_IESpositions_graphing.R
Validate_IESpositions_graphing.R plots how much of each IES the viterbi assignment captured (in % of correctly state assigned bps) - can compare to different model runs (2 v. 3 states, etc) 

# Check if the viterbi algorithim predits Novel IESs 
novelIEScalls.R and novelIEScalls.sh work similarly to validate_IESpositions.R/.sh - checks proportions of viterbi IES state assignments in regions outside Hamilition 2016 defined IESs. 
If there is a small proportion of high % regions they may be novel IESs - if there is a large proportion of high % regions the algorithim may be overreaching (as with a 3 state algorithim where both state 1 and 3 are called as IESs) 
Modified so that between IES intervals - i.e MDSs - cover nested/overlapping IESs 


# Identifying IESs in individual GE sequencing data 
in the /work/aahowel3/HMM/indiviudal_GEs/coverage_files folder the script coverage.sh generates indiviudal coverage.txt files from the bams generated in the /work/aahowel3/Ancestor folder - modified to now pull from sorted.bam files instead of rmdup.sorted.bam files

# Simulate observations in order to test the Viterbi and BW algorithims of the STAN package
simulate_hmm.R uses a Monte Carlo simulation to generate emissions data from the estimated parameters of the STAN package 
Input parameters of the simulations should be the output of the Viterbi and BW algorithims of the STAN package

# Check how much of an issue dips in coverage NEAR an IES will be for single samples
Updated using /work/Alignments
In HMM/coverage_IESboundaries script coverage_atmacexcisions_2.sh uses macexcisionsite.bed file (concatenated from the 5 files in flowsortdata/bam_IRS2 and first column removed to make it bed format) to pull reads from samtools merged ancestors GE aligned to mac reference from /work/Alignments/mac folder that cross the boundary into the file reads_atexcisionsites.bam 

Next step is finding where those reads end up in an ancestor GE to mic reference alignment - in searchreadsmic.sh breaks up reads_atexcisionsites.bam by first and second in pair, uses those names as text files to pull out those reads in the ancestor to mic alignments - sort by first and second in pair and count of those, how many are unmapped?

results - around just 0.14% (rough) of reads are dropped from these regions from MAC to MIC alignment - should not be causing artifical dips 

##############Scripts above here need to be updated to pull from new bamfiles in /work/aahowel3/Alignments 

# A different question - per sample, what is the coverage at each documented MAC excision site? 
Updated using /work/Alignments
If there is minimal coverage at a MAC excision site - could indicate an excision shift 
in coverage_IESboundariers/reads_atmacex_perancestor - reads_macex_peranc.sh - creates PERANC_reads_atexcisionssites.bam from the POOLED_reads_atexcisionssites.bam
calculate_macexcov.sh - for each PERANC_reads_atexcisionssites.bam file, for each MAC excision site - runs it through macex.R to calculate average depth per macex site - appends to an output file per GE ancestor ANC_macexcov.txt - file output per ancestor is a list of MACEX sites and their depths 

results - inconclusive metric - coverage is too low at sites and variable between A/B replicates - need to do single sequencing or more complex HMM to quantify excision shfiting

# What is the coverage per ancestor GE in a MIC alignment and a MAC alignment - this will tell us which to pull for deep sequencing  
Updated using /work/Alignments
In HMM/overallcoverage avgcov.sh and avgcov2.sh (one for running MIC and MAC covs simultaneously) - calculates samtools depth for each ancestor and averages it using avgcov.R - the final output is a single list of 1) anc 2) overall coverage - 2 final files one for MIC cov and MAC cov 

Validating these coverage calculations with another coverage tool in gatk - run using gatkcovestimates.sh 
