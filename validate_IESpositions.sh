#USAGE: command line input should look like $bash validate_IESpositions.sh 
#no input - just change experiment_run variable to do the IES validation for each type of run 

#change this each type of viterbi run you do
experiment_run=_nb_25_2state




#actual script
#split the viterbi file output of stan_hmm_nb.R into the 5 chromosome files
mkdir viterbi${experiment_run}
cd viterbi${experiment_run}
awk -v extension=_viterbidata${experiment_run}.txt '{print >> ($2 extension)}' ../viterbidata${experiment_run}.txt
rm *bp*


#you are in the experiment run folder with 5 files - each of the 1-5 viterbi files correspond to 1-5 tsv files with info on each chromosomes IES
#use the R script to parse each pair of viterbi and IES files 
#check if in the regions there are IES - the viterbi state assignment is correct
for arg in *viterbidata*
do
        file=$(basename "$arg" _viterbidata${experiment_run}.txt)
	Rscript /work/aahowel3/HMM/validate_IESpositions.R /work/aahowel3/flowsortdata/2931489_Howell/retention_scores/"${file}_IESs_inmic.tsv" /work/aahowel3/HMM/viterbi${experiment_run}/"${file}_viterbidata${experiment_run}.txt" > "${file}_validatedIESs.tsv"
done 
