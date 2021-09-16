parallel -j 12 -i sh -c "python3.5 Script_reset_DNA.py DifferentDNAstring_Trial{} 800000 " -- {6..10}
parallel -j 12 -i sh -c "python3.5 Script_reset_DNA.py DifDNA800kgen_Trial{} 800000 " -- {1..5}

!parallel -j 12 -i sh -c "python3.5 Script_reuse_DNA.py Parallel_Trial{} 600000 Scripting_data/Files_to_use_for_parallel_DNA_gen0_cds.txt Scripting_data/Files_to_use_for_parallel_DNA_map_cds_invariant50.txt" -- {1..12}

!python3.5 Script_reuse_DNA.py Parallel_Trial{} 600000 Scripting_data/Files_to_use_for_parallel_DNA_gen0_cds.txt Scripting_data/Files_to_use_for_parallel_DNA_map_cds_invariant50.txt


parallel -j 8 -i sh -c "python3.5 Script_reuse_DNA.py Parallel_Trial{} 600000 Scripting_data/Files_to_use_for_parallel_DNA_gen0_cds.txt Scripting_data/Files_to_use_for_parallel_DNA_map_cds_invariant50.txt" -- {5..10}

-j 8 -i sh -c "python3.5 Script_reuseDNA.py Parallel_Trial{} 600000 Scripting_data/Files_to_use_for_parallel_gen0_cds.txt Scripting_data/Files_to_use_for_parallel_DNA_map_cds_invariant50.txt" -- {5..10}


