# HMM for identifying excision boundaries and estimating ploidy

Studies indicate (Abyzov et al., [2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3106330/); Feng et al., [2017](https://academic.oup.com/nar/article/45/16/9481/4037355)) that the retention of IESs and their excision boundaries exhibit variability during conjugation, making mutation detection at these regions unreliable. To identify these complex regions in each MA sample I have implemented a Hidden Markov Model (HMM) to characterize genomic regions based on sequencing depth. Preliminary results indicate a 2-state, Poisson distributed HMM with a minimum of 50 Baum-Welch iterations is the best fit model. By incorporating information from the HMM into Denovogear the false positive rate of mutation discovery can be further decreased by including the likelihood of multiple genotypes at a variably retained/excised IES boundary.

![Alt text](https://github.com/aahowel3/HMM/blob/master/HMM_model_resize.png)

### Identifying IESs in pooled GE sequencing data 
`HMM_dataprocessing.sh` - aligns Ancestor R1 and R2 GEs one at a time to the reference rather than concatenating them and aligning all at once
uses samtools depth to create coverage file - this is the emissions data fed into the HMM
The output bam, sam, and coverage folders output to the main `/Ancestors` folder 

Final output of `HMM_dataprocessing.sh` is AncestorsGE_tomic_coverage.txt in `/HMM` - coverage depth every position in all chromosomes 
`AncestorGE_tomic_coverage.txt` and coverage folder is now for coverage data INCLUDING dup reads - removed the rmdup line of code in samtools file
`AncestorGE_tomic_coverage_2.txt` and coverage2 is now for coverage data generated using the samtools rmdup command 
Did not physically remove or comment out the samtools rmdup line in the HMM_dataprocessing script - but by using the intermediate sorted.bam files instead of the rmdup_sorted.bam files to create AncestorGE_tomic_coverage.txt the effect is the same 

R script `stan_hmm_nb.R` use the STAN package to take `ancestors_mic_aligned_coverage.txt` input and estimate parameters and output the viterbi assignment for each position. 
model type (Negative Binomial v. Poisson), number of states and number of Baum-Welch iterations can changed within the script in the variable header section 

The output of `stan_hmm_nb.R` is `viterbidata_exerpimentalvariables_state.txt` (you can change parameters - model, bw iterations, # of states by changing variable header section of stan_hmm_nb.R) 
modified `stan_hmm_nb.R` so that it also prints the parameter output to a file instead of the command line - `parameter_experimentvariables.txt`

![Alt text](https://github.com/aahowel3/HMM/blob/master/HMM_modelsample.png)

### Validate if the viterbi algorithim predits known (Hamilton et. al, 2016) IESs 
`stan_hmm_nb.R` gives you model parameters - you can check if the are accurate predictions by seeing how they assign states within known IES regions 

`validate_IESpositions.sh` uses the output of `stan_hmm.R` to create folders of the output broken up by chromosome (ex. /viterbidata_nb_25_2state or /viterbidata_nb_25_3state)

`validate_IESpositions.sh` then loops the chrX_viterbi file simultaneously with the chrX_IES_in_mic.tsv file in flowsortdata/retention_scores using the script `validate_IESpositions.R` 
IMPORTANT to check which viterbi state (1,2,3 ect) stan_HMM assigned to each group (IES v. Mac v. zero coverage) - check the parameters.txt file - modify this appropriately in the validate_IESpositions.R file. 

### Visualizing results of validate_IESpositions.sh 
Concatenate all the chrX tsv files together per folder for the input into script `Validate_IESpositions_graphing.R`
`Validate_IESpositions_graphing.R` plots how much of each IES the viterbi assignment captured (in % of correctly state assigned bps) - can compare to different model runs (2 v. 3 states, etc) 

![Alt text](https://github.com/aahowel3/HMM/blob/master/detectable_IES_resize.png)

### Check if the viterbi algorithim predits Novel IESs 
`novelIEScalls.R` and `novelIEScalls.sh` work similarly to `validate_IESpositions.R/.sh` - checks proportions of viterbi IES state assignments in regions outside Hamilition 2016 defined IESs. 
If there is a small proportion of high % regions they may be novel IESs - if there is a large proportion of high % regions the algorithim may be overreaching (as with a 3 state algorithim where both state 1 and 3 are called as IESs) 
Modified so that between IES intervals - i.e MDSs - cover nested/overlapping IESs - using IRanges function 

![Alt text](https://github.com/aahowel3/HMM/blob/master/novel_IES_resize.png)

### Simulate observations in order to test the Viterbi and BW algorithims of the STAN package
`simulate_hmm.R` uses a Monte Carlo simulation to generate emissions data from the estimated parameters of the STAN package 
Input parameters of the simulations should be the output of the Viterbi and BW algorithims of the STAN package

### Check how much of an issue dips in coverage NEAR an IES will be for single samples
In HMM/coverage_IESboundaries script `coverage_atmacexcisions_2.sh` uses macexcisionsite.bed file (concatenated from the 5 files in flowsortdata/bam_IRS2 and first column removed to make it bed format) to pull reads from samtools merged ancestors GE aligned to mac reference from /work/Alignments/mac folder that cross the boundary into the file reads_atexcisionsites.bam 

Next step is finding where those reads end up in an ancestor GE to mic reference alignment - in `searchreadsmic.sh` breaks up reads_atexcisionsites.bam by first and second in pair, uses those names as text files to pull out those reads in the ancestor to mic alignments - sort by first and second in pair and count of those, how many are unmapped?

Results - around just 0.14% (rough) of reads are dropped from these regions from MAC to MIC alignment - should not be causing artifical dips 

### Testing new validation metrics in HMM3
Anc 4 is the best coverage sample combined from replicate A and B (88x) - test original IES viterbi calling metric, Reed's junction likelihood metric using the forward algo, and try an additional validation with Pairties 

R script `forwardalgo_junlik.R` in local HMM - all you need from here is the STAN getPosterior function - otherwise same as `stanhmm.R`
`pulledoutviterbis.R` - all continuous stretches of vit2 with RLE function 
`HMMM_smoothing.R` - condensing those continuous stretches with IRanges
`EM_3comp_scratch.R` - removes large coverage levels from data 

### All updated alignments are in `/Alignments` 
The ancestor GE mic and mac alignment files in `/Ancestors` are incomplete - 11_A in particular has an empty sam file. This is due to the reads being incorrectly ordered in the R1s and R2s when concatenating the Lane 1-4s. The order can be fixed with the script `repairfqs.sh` that uses the bbmap package.
In `/Alignments` there are 3 folders - fastqs, mic, and mac 
Script `ancestors_aligntoX.sh` in the mic and mac folder generate a sam and bam folder - in bam folder are final indexed and sorted bam files 
