# **description of the Mutation equilibrium project** 

### **dependancies**
- I was using a conda env to perform these analyses. 
- I created a .yml file using conda env export > mutEq_env.yml

### **there are  afew notebooks that should be run first**
1. download the exon data 
    - using the notebook in the Human_exon_data/ directory

2. extract the count of each triplet in the genom 
    - HumanTripletCounts

3. create the triplet mutability model 
    - Human_mutability_model folder 

## **description of the remaining folders** 
- NonCoding_Coding_sims
    - the simulations involving coding and non coding seqeunce (ie. figure 2: sliding window, figure S2)
    - also some supplemental analyses: 
        - the analysis of the effect of adjacent triplets 
        - we explored a 5-mer model in earlier reviews but did not include in the final MS
- Intronic_simredo 
    - simulations involving only noncoding regions (ie figure 1: scatter plot of triplet counts and change)
    - also compared the equilibirum triplet freqeuncies to those found in intergenic regions in he genome (figure S4)
- comparing_equilibrium_point
    - calculating equilibrium point (@ 2x mutational coverage)
    - comparing the equilibrium of different seqeunce classes (ie figure S1). 
- codon_usage_bias 
    - the simulations addresing codon usage bias (ie figure S3)
- WholeGenomeMutabilityExonIntron
    - calculating the mutability ofgenomic coding and noncoding regions (ie figure S5)
- nonoverlapSimulation 
    - creating a simplified simulation that highlights the complexity in our original work (figure S6)

