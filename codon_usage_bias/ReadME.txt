# **description of the codon usage bias directory**

step 0: extract the RSCU of the genome 
- this is dont in the generating_codon_usage_bias_data folder

step 1: need to iterate through the exons and create the modified dna strings 
- this is done in the 3_nonSynCodingSim_non100kbpStart.ipynb notebook 

step 2: need to simulate the modified exons 
- this is done using the Script_reuse_DNAoutputFinal_countsMuts_multExons_coveageInput.py script and the following command 
- while read name; do parallel -j 10 -i sh -c "python3.5 Script_reuse_DNAoutputFinal_countMuts_multExons_coverageInput.py CodingUsage_cds"$name"_Trial{} 6 data/redo_multDif_exons_2021_02_07/CodingUsage_cds"$name"_DNA_gen0.txt  data/redo_multDif_exons_2021_02_07/CodingUsage_cds"$name"_MAP.txt" -- {1..10}; done < data/redo_multDif_exons_2021_02_07/cdsNames_list.txt\
- # while bash loop for iterating through the file "data/redo_multDif_exons_2021_02_07/cdsNames_list.txt" that contaisn the exon names
- #sys_argv[1] = CodingUsage_cds"$name"_Trial{}
- #sys_argv[2] = its 6x coverage 
- # the file paths for dna_gen0 and map to use must be written but have the "$name" in them 
- #<file.txt is the cdsNmes file that the loop is iterating through 

step 3 : need to wranlge the output of the simulation to compare human codon usage and codon chnage during the simulation
- this is done in the 3.1 notebook nad plotted in 3.2 (plots not included here as not inlcluded in our final supplemental)

step 4: wranlge the output and compare human codon usage to the RSCU found in the late stages of the simulation
- this is done in the 4.1 & 4.2 notebook and plotted in the plots/