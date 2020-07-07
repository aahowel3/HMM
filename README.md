# HMM for identifying excision boundaries and estimating ploidy

sequencing data from /work/aahowel3/Ancestor 1-20 A and B 
all scripts run in hines /work/aahowel3/HMM

#Identifying IESs in pooled GE sequencing data 
HMM_dataprocessing.sh - aligns Ancestor R1 and R2 GEs one at a time to the reference rather than concatenating them and aligning all at once
uses samtools depth to create coverage file - this is the emissions data fed into the HMM
