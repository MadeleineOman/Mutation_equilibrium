# **description of the folder where we performed standalone noncoding simulations**

1. generate and simulate random seqeunce to equilibrium 
- run the script_redoSim_2021_04_18_outputDNA_counts_newParallel.py script 
    - first argument = name of sim to print to file 
    - second argument = dna length 
    - third argument = list of mutation coverages to print (ie 0 = initial conditions, 0.1 = 10% mutation coverage, 2 = 200%)
    - bash command: "seq 10 | parallel python script_redoSim_2021_04_18_outputDNA_counts_newParallel.py cov2x_2021_04_18_Trial 100000 [0,0.1,2]"
    
2. wrangle output seqs into a dataframe 
- 1.2 notebook 

3. plot the outputed dataframe 
- 1.3 notebook 
    - plot available in plots/
- 1.31 notebook egenrates the plots we used to make the )alightly) adjusted figure for SMBE 2021

4. peform a regression for the final triplet frequencies 
- done for all triplets, and also without the CpGs
- 1.4 notebook 

5. compare the equilibirum triplet freqeuncies to those found in intergenic regions in he genome 
- 2 & 2.1 notebooks 