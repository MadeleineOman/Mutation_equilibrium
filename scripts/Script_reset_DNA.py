# IMPORT STATEMENTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import random
from random import randint
from typing import List, Dict
from numpy.random import choice
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import numpy
from Bio import SeqIO
from tqdm import tqdm
import copy 
import math 
import json
import sys 

#MAKING STATIC OBJECTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

aa_conversion_dict = {"TTT": "F", "TTC": "F", \
"TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", \
"ATT": "I", "ATC": "I", "ATA": "I", \
"ATG": "M", \
"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", \
"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S", \
"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", \
"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", \
"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", \
"TAT": "Y", "TAC": "Y", \
"CAT": "H", "CAC": "H", \
"CAA": "Q", "CAG": "Q", \
"AAT": "N", "AAC": "N", \
"AAA": "K", "AAG": "K", \
"GAT": "D", "GAC": "D", \
"GAA": "E", "GAG": "E", \
"TGT": "C", "TGC": "C", \
"TGG": "W", \
"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", \
"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G", \
"TAA": "*", "TAG": "*", "TGA": "*"} \

#FUNCTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def reverse_complement(dna):
    """
    Str --> str 
    Note that this function will revComp everything: need appropruate if conditional to make sure you are 64-->32 not 64--> complementary 64 
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "fuckup: used the revcomp"

def base_to_mutate(DNA, current_weights): 
    """
    (DNA: str, weights: Dict) -> int
    
    Will return the index of a random choice of a base that will be mutated based on probabilities given by the
    weights dictionary (i.e. the model)
    ex. "ATCGTA" --> index 3 ("G") will mutate
    """
    #DNA STR --> WEIGHTS LIST INDEXING DOUBLE CHECK 
    indices = []
    indices = list(i for i in range(1, len(DNA)-1))
    if (len(indices)) != (len(current_weights)): 
        print("ERROR, FUCKED UP")
        print("length of indices = " + str(len(indices)))
        print("length of weights = " + str(len(current_weights)))
    
    # NORMALIZE THE POPULATION OF WEIGHTS    
    total_freq = sum(current_weights) 
    normalized_weights = copy.copy(current_weights)
    for index, value in enumerate(current_weights):
        normalized_weights[index] = value/total_freq 

    # DRAW THE INDEX OF THE BASE THAT WILL BE MUTATED
    base_index = choice(indices, p=normalized_weights) 
    
    # RETURN THE INDEX
    return base_index

def current_frequencies(DNA, weights, index, current_weights): #and by weights, she means the "model" dictionary 
    """
    (DNA: str, weights list, index = chosen site) -> updated weights List 
    only chnages the 3 bases around the index chosen 
    """
    #CREATE THE 3 SOURROUNDING TRIPLETS
    triplet_preceding = DNA[index-2:index+1] #dont worry this indexing is double checked 
    triplet_center = DNA[index-1:index+2]
    triplet_following = DNA[index:index+3]
    triplet_list = [triplet_preceding, triplet_center, triplet_following]
    
    #64-->32 
    for triplet in triplet_list: 
        if triplet not in weights: 
            i = triplet_list.index(triplet)
            triplet_list[i] = reverse_complement(triplet)
    
    #GET THE WEIGHTS FOR EACH TRIPLET AND UPDATE WEIGHTS (LIST) OBJECT 
    current_weights[index-2]= weights[triplet_list[0]][0]    #note that it's the expectedindex -1 as the index = position in string (from 1 to n-1) but here the items are the codons ...? 
    current_weights[index-1]= weights[triplet_list[1]][0]
    current_weights[index] = weights[triplet_list[2]][0]
    
    return current_weights

def base_to_turn_into(i, model, DNA):
    """
    (i: int, weights: Dict, DNA: str) -> str (the base A/T/C/G)
    """
      
    # GRAB THE TRIPLET AND COMPLEMENT IT IF NEED BE.
    triplet = DNA[i - 1:i+2] 
    if triplet not in model:
            triplet = reverse_complement(triplet)
            
    # GENERATE A LIST OF THE TRIPLETS AND WEIGHTS OF THE THREE TRIPLETS THE CHOSEN TRIPLET CAN TURN INTO 
    triplets = list(key for key in model[triplet][1].keys())# ex. "ATC" --> ["AAC", "AGC", "ACC"] (triplets)
    current_weights_to_turn_into = list(value for value in model[triplet][1].values()) # ex. "ATCGTA" --> [0.3, 0.5, 0.1, 0.1] (weights of triplet base changes) (these are already normalized!)
    
    # MAKE AND RETURN THE CHOICE OF BASE TO TURN INTO.
    return choice(triplets, p = current_weights_to_turn_into)[1] # <-- the [1] will select the central base (middle N base out of strign 3)

def codingYesNo(ab_index, mb, DNA, DNA_map):
    """
    (ab_index: int, mb = base A/T/C/G, DNA: str, DNA_map: List) -> List or None
    
    Will check if the ancestral base is part of the coding/ non-coding/invariant region.
    If part of a noncoding region, it will return "None".
    If part of a coding region, it will return the ancestral and mutated codons in a list with two entries [ancestral_codon, mutated_codon] 
    If part of an invariant region, returns int(3)
    """
    if DNA_map[ab_index] == "N": # if it's not non-coding (then it's coding)...
        return None
    elif DNA_map[ab_index] == "3": #if invariant 
        return 3
    else:
        if DNA_map[ab_index] == "0": # if we're at the first position...
            codon = DNA[ab_index] + DNA[ab_index + 1] + DNA[ab_index + 2] # take the 012 codon
            m_codon = mb + DNA[ab_index + 1] + DNA[ab_index + 2] # make the mutated codon
        elif DNA_map[ab_index] == "1": # if we're at the second position...
            codon = DNA[ab_index - 1] + DNA[ab_index] + DNA[ab_index + 1] # take the 012 codon
            m_codon = DNA[ab_index - 1] + mb + DNA[ab_index + 1] # make the mutated codon
        elif DNA_map[ab_index] == "2": # if we're at the third position...
            codon = DNA[ab_index - 2] + DNA[ab_index - 1] + DNA[ab_index] # take the 012 codon
            m_codon = DNA[ab_index - 2] + DNA[ab_index - 1] + mb # make the mutated codon
        
        return [codon, m_codon] # if it's coding, return the ancestral and mutated codons

def blosum(a_m_codons, codon_aa_conversionDict, blosum_dict):
    """
    (a_m_codons: List [ancestral_codon, mutated_codon], codon_aa_conversionDict = triplet-->aa dict, blosum_dict = blosum matrix ) -> bool (True/False     
    Will return True if the change is acceptable and False if the change is not acceptable.  
    ex. ["CGT", "CAT"] --> R to H is accepted at 0.6, or 60%. If the function accepts it, it will return True.
    """
    
    a_aa = codon_aa_conversionDict[a_m_codons[0]] #ancestral a.a 
    m_aa = codon_aa_conversionDict[a_m_codons[1]]  #mutated a.a 
    
    return choice([True, False], p = [blosum_dict[a_aa][m_aa], (1-(blosum_dict[a_aa][m_aa]))]) 

def mutate_chromosome(DNA, ab_index, mb):
    """
    (DNA:str, ab_index:number, site chosen, mb:str) -> str (mutated DNA)
    Will return a new DNA string with a random base mutated to a new base as given by the previous functions.
    """
    return DNA[:ab_index] + mb + DNA[ab_index + 1:]

def average_mutability(current_weights):
    """
    (weights: list) -> float
    needs to come before the current_frequencies function 
    Will return a single mutability score which is the average of the triplet mutation frequencies in the string.
    ex. "ATCGTA" --> [0.02, 0.31, 0.15, 0.07] --> (0.02+0.31+0.15+0.07)/4 --> 0.1375
    """
    mean = np.mean(current_weights)
    return mean # return the average (number) 

def exon_mutability(DNA_map,current_weights):
    """
    (DNA_map = str, current_weights = list) -> float
    Same function as average_mutability, but for CODING REGIONs ONLY (why you need the map) 
    """
    #SUM MUTABILITY ACROSS THE ENTIRE DNA MAP STRING USING DNA_MAP IF CONDITIONAL 
    mutability = 0 
    count = 0 
    for i in range(1, len(DNA_map)-1): # this will exclude the first base and the last base
        if DNA_map[i] != "N": # if we are in the coding region
            mutability += current_weights[i-1]# add the mutability of the coding triplet to the total --> and YES ive chedked this, dont freak out the indices match up... 
            count += 1             
    return mutability/count 

def intron_mutability(DNA_map, current_weights):
    """
    (DNA_map = str, current_weights = list) -> float
    Same function as average_mutability, but for nonCODING REGIONs ONLY (why you need the map)
    """
    #SUM MUTABILITY ACROSS THE ENTIRE DNA MAP STRING USING DNA_MAP IF CONDITIONAL 
    mutability = 0 
    count = 0
    for i in range(1, len(DNA_map)-1): # this will exclude the first base and the last base
        if DNA_map[i] == "N": # if we are in the coding region
            mutability += current_weights[i-1]# add the mutability of the coding triplet to the total --> and YES ive chedked this, dont freak out the indices match up... 
            count += 1
    return mutability/count 

#MAKING VARIABLE OBJECTS 

#dictionary of CDS greater than a certain size (here >4k) 
cds_file_obj = open("Human_exon_data/Homo_sapiens.GRCh38.cds.all.fa")
total_cds_dict = {}
for gene in SeqIO.parse(cds_file_obj, "fasta"): 
    if len(gene.seq) >= 4000: 
        cds_map = "012"*int(len(gene.seq)/3)
        if len(cds_map) == len(gene.seq) and "N" not in gene.seq: 
            total_cds_dict[gene.id] = [str(gene.seq), cds_map]
print("total number of CDS where len >= 4k is "+str(len(total_cds_dict)))                                #! print statement? 

#blosum matrix 
blosum_dict = {'T': {'T': 1.0, 'Q': 0.5, 'C': 0.4, 'M': 0.5, 'A': 0.6, 'B': 0.5, 'K': 0.5, 'E': 0.5, 'P': 0.4, 'R': 0.4, 'G': 0.3, 'J': 0.4, 'F': 0.3, 'X': 0.5, 'Z': 0.5, '*': 0.0, 'N': 0.6, 'V': 0.5, 'H': 0.4, 'W': 0.2, 'L': 0.4, 'S': 0.7, 'D': 0.4, 'I': 0.5, 'Y': 0.4}, \
 'Q': {'T': 0.5, 'Q': 1.0, 'C': 0.2, 'M': 0.6, 'A': 0.5, 'B': 0.5, 'K': 0.7, 'E': 0.8, 'P': 0.4, 'R': 0.7, 'G': 0.3, 'J': 0.3, 'F': 0.2, 'X': 0.5, 'Z': 1.0, '*': 0.0, 'N': 0.6, 'V': 0.3, 'H': 0.7, 'W': 0.3, 'L': 0.3, 'S': 0.5, 'D': 0.5, 'I': 0.2, 'Y': 0.3}, \
 'C': {'T': 0.4, 'Q': 0.2, 'C': 1.0, 'M': 0.4, 'A': 0.5, 'B': 0.2, 'K': 0.2, 'E': 0.0, 'P': 0.2, 'R': 0.1, 'G': 0.2, 'J': 0.4, 'F': 0.3, 'X': 0.5, 'Z': 0.1, '*': 0.0, 'N': 0.2, 'V': 0.4, 'H': 0.1, 'W': 0.2, 'L': 0.4, 'S': 0.4, 'D': 0.1, 'I': 0.4, 'Y': 0.2}, \
 'M': {'T': 0.5, 'Q': 0.6, 'C': 0.4, 'M': 1.0, 'A': 0.4, 'B': 0.2, 'K': 0.4, 'E': 0.3, 'P': 0.3, 'R': 0.4, 'G': 0.2, 'J': 0.8, 'F': 0.5, 'X': 0.5, 'Z': 0.4, '*': 0.0, 'N': 0.3, 'V': 0.6, 'H': 0.3, 'W': 0.4, 'L': 0.8, 'S': 0.4, 'D': 0.2, 'I': 0.7, 'Y': 0.4}, \
 'A': {'T': 0.6, 'Q': 0.5, 'C': 0.5, 'M': 0.4, 'A': 1.0, 'B': 0.4, 'K': 0.5, 'E': 0.5, 'P': 0.5, 'R': 0.4, 'G': 0.6, 'J': 0.4, 'F': 0.3, 'X': 0.5, 'Z': 0.5, '*': 0.0, 'N': 0.4, 'V': 0.5, 'H': 0.4, 'W': 0.2, 'L': 0.4, 'S': 0.7, 'D': 0.3, 'I': 0.4, 'Y': 0.3}, \
 'B': {'T': 0.5, 'Q': 0.5, 'C': 0.2, 'M': 0.2, 'A': 0.4, 'B': 1.0, 'K': 0.5, 'E': 0.7, 'P': 0.3, 'R': 0.4, 'G': 0.4, 'J': 0.1, 'F': 0.2, 'X': 0.5, 'Z': 0.6, '*': 0.0, 'N': 1.0, 'V': 0.2, 'H': 0.5, 'W': 0.0, 'L': 0.1, 'S': 0.6, 'D': 1.0, 'I': 0.1, 'Y': 0.2}, \
 'K': {'T': 0.5, 'Q': 0.7, 'C': 0.2, 'M': 0.4, 'A': 0.5, 'B': 0.5, 'K': 1.0, 'E': 0.6, 'P': 0.4, 'R': 0.8, 'G': 0.4, 'J': 0.3, 'F': 0.2, 'X': 0.5, 'Z': 0.7, '*': 0.0, 'N': 0.6, 'V': 0.3, 'H': 0.5, 'W': 0.1, 'L': 0.3, 'S': 0.5, 'D': 0.5, 'I': 0.2, 'Y': 0.3}, \
 'E': {'T': 0.5, 'Q': 0.8, 'C': 0.0, 'M': 0.3, 'A': 0.5, 'B': 0.7, 'K': 0.6, 'E': 1.0, 'P': 0.4, 'R': 0.5, 'G': 0.3, 'J': 0.2, 'F': 0.1, 'X': 0.5, 'Z': 1.0, '*': 0.0, 'N': 0.5, 'V': 0.3, 'H': 0.5, 'W': 0.1, 'L': 0.2, 'S': 0.5, 'D': 0.7, 'I': 0.2, 'Y': 0.2}, \
 'P': {'T': 0.4, 'Q': 0.4, 'C': 0.2, 'M': 0.3, 'A': 0.5, 'B': 0.3, 'K': 0.4, 'E': 0.4, 'P': 1.0, 'R': 0.3, 'G': 0.3, 'J': 0.2, 'F': 0.2, 'X': 0.5, 'Z': 0.4, '*': 0.0, 'N': 0.3, 'V': 0.3, 'H': 0.3, 'W': 0.1, 'L': 0.2, 'S': 0.4, 'D': 0.3, 'I': 0.2, 'Y': 0.2}, \
 'R': {'T': 0.4, 'Q': 0.7, 'C': 0.1, 'M': 0.4, 'A': 0.4, 'B': 0.4, 'K': 0.8, 'E': 0.5, 'P': 0.3, 'R': 1.0, 'G': 0.3, 'J': 0.3, 'F': 0.2, 'X': 0.5, 'Z': 0.6, '*': 0.0, 'N': 0.5, 'V': 0.3, 'H': 0.6, 'W': 0.2, 'L': 0.3, 'S': 0.5, 'D': 0.3, 'I': 0.2, 'Y': 0.3}, \
 'G': {'T': 0.3, 'Q': 0.3, 'C': 0.2, 'M': 0.2, 'A': 0.6, 'B': 0.4, 'K': 0.4, 'E': 0.3, 'P': 0.3, 'R': 0.3, 'G': 1.0, 'J': 0.1, 'F': 0.1, 'X': 0.5, 'Z': 0.3, '*': 0.0, 'N': 0.5, 'V': 0.1, 'H': 0.3, 'W': 0.2, 'L': 0.1, 'S': 0.5, 'D': 0.4, 'I': 0.1, 'Y': 0.1}, \
 'J': {'T': 0.4, 'Q': 0.3, 'C': 0.4, 'M': 0.8, 'A': 0.4, 'B': 0.1, 'K': 0.3, 'E': 0.2, 'P': 0.2, 'R': 0.3, 'G': 0.1, 'J': 1.0, 'F': 0.6, 'X': 0.5, 'Z': 0.2, '*': 0.0, 'N': 0.2, 'V': 0.7, 'H': 0.2, 'W': 0.3, 'L': 1.0, 'S': 0.3, 'D': 0.1, 'I': 0.9, 'Y': 0.4}, \
 'F': {'T': 0.3, 'Q': 0.2, 'C': 0.3, 'M': 0.5, 'A': 0.3, 'B': 0.2, 'K': 0.2, 'E': 0.1, 'P': 0.2, 'R': 0.2, 'G': 0.1, 'J': 0.6, 'F': 1.0, 'X': 0.5, 'Z': 0.2, '*': 0.0, 'N': 0.2, 'V': 0.4, 'H': 0.4, 'W': 0.6, 'L': 0.6, 'S': 0.3, 'D': 0.1, 'I': 0.5, 'Y': 0.9}, \
 'X': {'T': 0.5, 'Q': 0.5, 'C': 0.5, 'M': 0.5, 'A': 0.5, 'B': 0.5, 'K': 0.5, 'E': 0.5, 'P': 0.5, 'R': 0.5, 'G': 0.5, 'J': 0.5, 'F': 0.5, 'X': 0.5, 'Z': 0.5, '*': 0.0, 'N': 0.5, 'V': 0.5, 'H': 0.5, 'W': 0.5, 'L': 0.5, 'S': 0.5, 'D': 0.5, 'I': 0.5, 'Y': 0.5}, \
 'Z': {'T': 0.5, 'Q': 1.0, 'C': 0.1, 'M': 0.4, 'A': 0.5, 'B': 0.6, 'K': 0.7, 'E': 1.0, 'P': 0.4, 'R': 0.6, 'G': 0.3, 'J': 0.2, 'F': 0.2, 'X': 0.5, 'Z': 1.0, '*': 0.0, 'N': 0.5, 'V': 0.3, 'H': 0.6, 'W': 0.2, 'L': 0.2, 'S': 0.5, 'D': 0.7, 'I': 0.2, 'Y': 0.3}, \
 '*': {'T': 0.0, 'Q': 0.0, 'C': 0.0, 'M': 0.0, 'A': 0.0, 'B': 0.0, 'K': 0.0, 'E': 0.0, 'P': 0.0, 'R': 0.0, 'G': 0.0, 'J': 0.0, 'F': 0.0, 'X': 0.0, 'Z': 0.0, '*': 0.7, 'N': 0.0, 'V': 0.0, 'H': 0.0, 'W': 0.0, 'L': 0.0, 'S': 0.0, 'D': 0.0, 'I': 0.0, 'Y': 0.0}, \
 'N': {'T': 0.6, 'Q': 0.6, 'C': 0.2, 'M': 0.3, 'A': 0.4, 'B': 1.0, 'K': 0.6, 'E': 0.5, 'P': 0.3, 'R': 0.5, 'G': 0.5, 'J': 0.2, 'F': 0.2, 'X': 0.5, 'Z': 0.5, '*': 0.0, 'N': 1.0, 'V': 0.2, 'H': 0.6, 'W': 0.1, 'L': 0.2, 'S': 0.6, 'D': 0.7, 'I': 0.2, 'Y': 0.3}, \
 'V': {'T': 0.5, 'Q': 0.3, 'C': 0.4, 'M': 0.6, 'A': 0.5, 'B': 0.2, 'K': 0.3, 'E': 0.3, 'P': 0.3, 'R': 0.3, 'G': 0.1, 'J': 0.7, 'F': 0.4, 'X': 0.5, 'Z': 0.3, '*': 0.0, 'N': 0.2, 'V': 1.0, 'H': 0.2, 'W': 0.3, 'L': 0.6, 'S': 0.4, 'D': 0.1, 'I': 0.9, 'Y': 0.3}, \
 'H': {'T': 0.4, 'Q': 0.7, 'C': 0.1, 'M': 0.3, 'A': 0.4, 'B': 0.5, 'K': 0.5, 'E': 0.5, 'P': 0.3, 'R': 0.6, 'G': 0.3, 'J': 0.2, 'F': 0.4, 'X': 0.5, 'Z': 0.6, '*': 0.0, 'N': 0.6, 'V': 0.2, 'H': 1.0, 'W': 0.3, 'L': 0.2, 'S': 0.4, 'D': 0.4, 'I': 0.2, 'Y': 0.7}, \
 'W': {'T': 0.2, 'Q': 0.3, 'C': 0.2, 'M': 0.4, 'A': 0.2, 'B': 0.0, 'K': 0.1, 'E': 0.1, 'P': 0.1, 'R': 0.2, 'G': 0.2, 'J': 0.3, 'F': 0.6, 'X': 0.5, 'Z': 0.2, '*': 0.0, 'N': 0.1, 'V': 0.3, 'H': 0.3, 'W': 1.0, 'L': 0.3, 'S': 0.2, 'D': 0.0, 'I': 0.2, 'Y': 0.8}, \
 'L': {'T': 0.4, 'Q': 0.3, 'C': 0.4, 'M': 0.8, 'A': 0.4, 'B': 0.1, 'K': 0.3, 'E': 0.2, 'P': 0.2, 'R': 0.3, 'G': 0.1, 'J': 1.0, 'F': 0.6, 'X': 0.5, 'Z': 0.2, '*': 0.0, 'N': 0.2, 'V': 0.6, 'H': 0.2, 'W': 0.3, 'L': 1.0, 'S': 0.3, 'D': 0.1, 'I': 0.7, 'Y': 0.4}, \
 'S': {'T': 0.7, 'Q': 0.5, 'C': 0.4, 'M': 0.4, 'A': 0.7, 'B': 0.6, 'K': 0.5, 'E': 0.5, 'P': 0.4, 'R': 0.5, 'G': 0.5, 'J': 0.3, 'F': 0.3, 'X': 0.5, 'Z': 0.5, '*': 0.0, 'N': 0.6, 'V': 0.4, 'H': 0.4, 'W': 0.2, 'L': 0.3, 'S': 1.0, 'D': 0.5, 'I': 0.3, 'Y': 0.3}, \
 'D': {'T': 0.4, 'Q': 0.5, 'C': 0.1, 'M': 0.2, 'A': 0.3, 'B': 1.0, 'K': 0.5, 'E': 0.7, 'P': 0.3, 'R': 0.3, 'G': 0.4, 'J': 0.1, 'F': 0.1, 'X': 0.5, 'Z': 0.7, '*': 0.0, 'N': 0.7, 'V': 0.1, 'H': 0.4, 'W': 0.0, 'L': 0.1, 'S': 0.5, 'D': 1.0, 'I': 0.1, 'Y': 0.2}, \
 'I': {'T': 0.5, 'Q': 0.2, 'C': 0.4, 'M': 0.7, 'A': 0.4, 'B': 0.1, 'K': 0.2, 'E': 0.2, 'P': 0.2, 'R': 0.2, 'G': 0.1, 'J': 0.9, 'F': 0.5, 'X': 0.5, 'Z': 0.2, '*': 0.0, 'N': 0.2, 'V': 0.9, 'H': 0.2, 'W': 0.2, 'L': 0.7, 'S': 0.3, 'D': 0.1, 'I': 1.0, 'Y': 0.4}, \
 'Y': {'T': 0.4, 'Q': 0.3, 'C': 0.2, 'M': 0.4, 'A': 0.3, 'B': 0.2, 'K': 0.3, 'E': 0.2, 'P': 0.2, 'R': 0.3, 'G': 0.1, 'J': 0.4, 'F': 0.9, 'X': 0.5, 'Z': 0.3, '*': 0.0, 'N': 0.3, 'V': 0.3, 'H': 0.7, 'W': 0.8, 'L': 0.4, 'S': 0.3, 'D': 0.2, 'I': 0.4, 'Y': 1.0}}

#mutability model 
model = {'GTC': [0.38170347003154576, {'GCC': 0.581267217630854, 'GGC': 0.17447199265381083, 'GAC': 0.24426078971533516}], 'TGA': [0.44338655339094774, {'TAA': 0.6247109349190618, 'TTA': 0.12685827552031714, 'TCA': 0.2484307895606211}], 'TAT': [0.544423228125351, {'TTT': 0.12131215184650299, 'TGT': 0.8021456571074892, 'TCT': 0.07654219104600785}], 'CGC': [0.9428571428571428, {'CTC': 0.06703397612488522, 'CCC': 0.027548209366391185, 'CAC': 0.9054178145087236}], 'ATT': [0.48758198043221157, {'AGT': 0.1133406835722161, 'ACT': 0.7592061742006616, 'AAT': 0.1274531422271224}], 'GCA': [0.49819102749638206, {'GTA': 0.6016702977487292, 'GAA': 0.24473493100944083, 'GGA': 0.15359477124183007}], 'CGT': [0.9583888149134487, {'CTT': 0.04341785342132685, 'CCT': 0.03403959708232025, 'CAT': 0.922542549496353}], 'CCA': [0.4589957500393515, {'CTA': 0.6903292181069959, 'CGA': 0.20781893004115226, 'CAA': 0.10185185185185185}], 'GGA': [0.5196784458214705, {'GTA': 0.17015791169835642, 'GAA': 0.602964872703835, 'GCA': 0.22687721559780857}], 'AAA': [0.3198252625708709, {'ACA': 0.2961348445219413, 'ATA': 0.17988956698634118, 'AGA': 0.5239755884917175}], 'AAC': [0.3978541712283775, {'ACC': 0.17776554760594387, 'AGC': 0.6692350027517887, 'ATC': 0.15299944964226747}], 'CTC': [0.3202682875707399, {'CGC': 0.21662303664921467, 'CCC': 0.5530104712041884, 'CAC': 0.23036649214659685}], 'AGG': [0.5123493090267568, {'ACG': 0.22955523672883785, 'AAG': 0.6662840746054519, 'ATG': 0.10416068866571017}], 'AGA': [0.4717741935483871, {'ACA': 0.34816490698843644, 'AAA': 0.47611865258924085, 'ATA': 0.1757164404223228}], 'AGC': [0.5081906865451868, {'ACC': 0.20101412531691415, 'AAC': 0.6356392611372691, 'ATC': 0.16334661354581673}], 'CAA': [0.36173285198555954, {'CTA': 0.13822355289421157, 'CGA': 0.5768463073852296, 'CCA': 0.28493013972055886}], 'CTT': [0.35645079041305455, {'CGT': 0.2341440152598951, 'CCT': 0.6132570338578922, 'CAT': 0.15259895088221268}], 'CAC': [0.3995351785336996, {'CCC': 0.20359598096245374, 'CTC': 0.18244315177154943, 'CGC': 0.6139608672659969}], 'TAG': [0.38808618504435993, {'TCG': 0.16525146962769433, 'TGG': 0.6962769431743958, 'TTG': 0.13847158719790986}], 'ACA': [0.5516478655164787, {'ATA': 0.6555733761026463, 'AGA': 0.17802726543704891, 'AAA': 0.16639935846030474}], 'CTG': [0.4073735527117611, {'CGG': 0.20905011219147346, 'CAG': 0.1480927449513837, 'CCG': 0.6428571428571429}], 'ACT': [0.5261127596439169, {'AGT': 0.24844895657078397, 'ATT': 0.5964467005076142, 'AAT': 0.1551043429216018}], 'TTA': [0.3860215053763441, {'TGA': 0.20533227218463987, 'TAA': 0.22682053322721846, 'TCA': 0.5678471945881417}], 'CGA': [0.9197916666666667, {'CTA': 0.026896942242355604, 'CAA': 0.9210079275198187, 'CCA': 0.052095130237825596}], 'GTA': [0.4143262045864468, {'GGA': 0.12686567164179105, 'GAA': 0.1455223880597015, 'GCA': 0.7276119402985075}], 'GGC': [0.5527913809990206, {'GCC': 0.1630049610205528, 'GTC': 0.23848334514528702, 'GAC': 0.5985116938341601}], 'AGT': [0.5261127596439169, {'ATT': 0.1551043429216018, 'ACT': 0.24844895657078397, 'AAT': 0.5964467005076142}], 'AAT': [0.48758198043221157, {'AGT': 0.7592061742006616, 'ATT': 0.1274531422271224, 'ACT': 0.1133406835722161}], 'CCG': [0.943577893317928, {'CTG': 0.943265306122449, 'CGG': 0.04673469387755102, 'CAG': 0.01}], 'ACC': [0.6345660930062248, {'AAC': 0.3147720715522216, 'ATC': 0.532602423542989, 'AGC': 0.15262550490478938}], 'TCC': [0.5196784458214705, {'TTC': 0.602964872703835, 'TGC': 0.22687721559780857, 'TAC': 0.17015791169835642}], 'TCT': [0.4717741935483871, {'TTT': 0.47611865258924085, 'TGT': 0.34816490698843644, 'TAT': 0.1757164404223228}], 'CAG': [0.4073735527117611, {'CTG': 0.1480927449513837, 'CGG': 0.6428571428571429, 'CCG': 0.20905011219147346}], 'TTG': [0.36173285198555954, {'TCG': 0.5768463073852296, 'TGG': 0.28493013972055886, 'TAG': 0.13822355289421157}], 'ATC': [0.4466903598400711, {'ACC': 0.5862754848334162, 'AAC': 0.2973644952759821, 'AGC': 0.11636001989060167}], 'CGG': [0.943577893317928, {'CTG': 0.01, 'CAG': 0.943265306122449, 'CCG': 0.04673469387755102}], 'TTC': [0.30323054331864907, {'TGC': 0.2391041162227603, 'TCC': 0.5617433414043583, 'TAC': 0.19915254237288135}], 'TAA': [0.3860215053763441, {'TGA': 0.5678471945881417, 'TTA': 0.22682053322721846, 'TCA': 0.20533227218463987}], 'ATG': [0.5741935483870968, {'ACG': 0.716724286949006, 'AAG': 0.1497407087294728, 'AGG': 0.13353500432152118}], 'GGT': [0.6345660930062248, {'GTT': 0.3147720715522216, 'GCT': 0.15262550490478938, 'GAT': 0.532602423542989}], 'TGC': [0.49819102749638206, {'TTC': 0.24473493100944083, 'TCC': 0.15359477124183007, 'TAC': 0.6016702977487292}], 'GTG': [0.3995351785336996, {'GAG': 0.18244315177154943, 'GCG': 0.6139608672659969, 'GGG': 0.20359598096245374}], 'GGG': [0.5453669813138123, {'GAG': 0.6450017661603674, 'GTG': 0.11974567290709998, 'GCG': 0.2352525609325327}], 'ACG': [0.9583888149134487, {'ATG': 0.922542549496353, 'AAG': 0.04341785342132685, 'AGG': 0.03403959708232025}], 'TGG': [0.4589957500393515, {'TCG': 0.20781893004115226, 'TTG': 0.10185185185185185, 'TAG': 0.6903292181069959}], 'GAA': [0.30323054331864907, {'GTA': 0.19915254237288135, 'GGA': 0.5617433414043583, 'GCA': 0.2391041162227603}], 'GAC': [0.38170347003154576, {'GCC': 0.17447199265381083, 'GTC': 0.24426078971533516, 'GGC': 0.581267217630854}], 'TAC': [0.4143262045864468, {'TTC': 0.1455223880597015, 'TGC': 0.7276119402985075, 'TCC': 0.12686567164179105}], 'GAT': [0.4466903598400711, {'GTT': 0.2973644952759821, 'GCT': 0.11636001989060167, 'GGT': 0.5862754848334162}], 'ATA': [0.544423228125351, {'ACA': 0.8021456571074892, 'AAA': 0.12131215184650299, 'AGA': 0.07654219104600785}], 'GCT': [0.5081906865451868, {'GTT': 0.6356392611372691, 'GGT': 0.20101412531691415, 'GAT': 0.16334661354581673}], 'CAT': [0.5741935483870968, {'CTT': 0.1497407087294728, 'CGT': 0.716724286949006, 'CCT': 0.13353500432152118}], 'CTA': [0.38808618504435993, {'CGA': 0.16525146962769433, 'CAA': 0.13847158719790986, 'CCA': 0.6962769431743958}], 'CCT': [0.5123493090267568, {'CTT': 0.6662840746054519, 'CGT': 0.22955523672883785, 'CAT': 0.10416068866571017}], 'TTT': [0.3198252625708709, {'TCT': 0.5239755884917175, 'TGT': 0.2961348445219413, 'TAT': 0.17988956698634118}], 'TCG': [0.9197916666666667, {'TGG': 0.052095130237825596, 'TTG': 0.9210079275198187, 'TAG': 0.026896942242355604}], 'AAG': [0.35645079041305455, {'ATG': 0.15259895088221268, 'ACG': 0.2341440152598951, 'AGG': 0.6132570338578922}], 'GTT': [0.3978541712283775, {'GCT': 0.6692350027517887, 'GGT': 0.17776554760594387, 'GAT': 0.15299944964226747}], 'TCA': [0.44338655339094774, {'TAA': 0.12685827552031714, 'TGA': 0.2484307895606211, 'TTA': 0.6247109349190618}], 'CCC': [0.5453669813138123, {'CTC': 0.6450017661603674, 'CGC': 0.2352525609325327, 'CAC': 0.11974567290709998}], 'GAG': [0.3202682875707399, {'GTG': 0.23036649214659685, 'GCG': 0.21662303664921467, 'GGG': 0.5530104712041884}], 'GCG': [0.9428571428571428, {'GAG': 0.06703397612488522, 'GTG': 0.9054178145087236, 'GGG': 0.027548209366391185}], 'GCC': [0.5527913809990206, {'GTC': 0.5985116938341601, 'GGC': 0.1630049610205528, 'GAC': 0.23848334514528702}], 'TGT': [0.5516478655164787, {'TCT': 0.17802726543704891, 'TAT': 0.6555733761026463, 'TTT': 0.16639935846030474}]}


#setting DNA length. also note that exon_length_desired is the minumum 
length_DNA = 100000
exon_length_desired = int(length_DNA/2)
intron_total_length = int(length_DNA/2)

#choose the exons to use in this simulation (total length should be about exon_length_desired: but minumum)
exons_to_use_dict = {}  
exon_length_realized = 0 
while exon_length_realized < exon_length_desired: 
    exon_key = random.choice(list(total_cds_dict))
    exon_string = total_cds_dict[exon_key][0]
    exon_map = total_cds_dict[exon_key][1]
    
    exons_to_use_dict[exon_key] = [exon_string, exon_map]
    exon_length_realized += len(exon_string)

#creating and printing key attributes of the DNA string we're gonna make 
number_of_exons = len(exons_to_use_dict)
print("number_of_exons is " + str(number_of_exons))
intron_gap = int(intron_total_length/number_of_exons)
print("intron_gap is " + str(intron_gap))

#making intron string from random choice of A/T/C/G. note that every intron is the same. 
intron_string = ""
for pos in range(0,intron_gap): 
    intron_string += choice(["A", "T", "C", "G"])
intron_map_string = "N"*len(intron_string)
print("note that every intron is the same")

#creating the DNA and DNA_map string while inserting cds 
DNA = ""
DNA_map = ""
exons_insertion_position_dict = {}
for exon_key, exon_item_list in exons_to_use_dict.items(): 
    DNA += intron_string
    
    exons_insertion_position_dict[exon_key] = [len(DNA)-1]
    DNA += exon_item_list[0]
    exons_insertion_position_dict[exon_key].append(len(DNA)-1)
    print("length of one of the exons is "+str(len(exon_item_list[0])))
    
    DNA_map += intron_map_string
    DNA_map += exon_item_list[1]

#randomly choosing half of CDS to be invariant sites (3 in the DNA_map)
number_codingSites_turnedInvariant = 0 
while number_codingSites_turnedInvariant <= exon_length_realized/2: 
    site_chosen = random.randrange(len(DNA))
    if (DNA_map[site_chosen]) != "N" and (DNA_map[site_chosen]) != "3" : 
        DNA_map = DNA_map[:site_chosen] + "3" + DNA_map[site_chosen+1:]
        number_codingSites_turnedInvariant += 1
        
#creating the list of initial weights 
weights = []
for i in range(1,len(DNA)-1): # dont want it to run all the way to the end: need to -1 so last centre index is second last in DNA string 
    if (DNA[i-1:i+2]) in model: 
        weights.append(model[(DNA[i-1:i+2])][0])
    else: 
        weights.append(model[reverse_complement(DNA[i-1:i+2])][0])

#writing initial conditions to file 
NameForThisSim = sys.argv[1]
text_file = open("Scripting_data_trial2/{t}_DNA_gen0_cds.txt".format(t=NameForThisSim), "w")                                               #! need to change file creation: use the {}.format... file structure? 
n = text_file.write(DNA)
text_file.close()
text_file = open("Scripting_data_trial2/{t}_DNA_map_cds_invariant50.txt".format(t=NameForThisSim), "w")
n = text_file.write(DNA_map)
text_file.close()
text_file = open("Scripting_data_trial2/{t}_exon_insertion_dict.txt".format(t=NameForThisSim), "w")
n = text_file.write(json.dumps(exons_insertion_position_dict))
text_file.close()
text_file = open("Scripting_data_trial2/{t}_exons_toUseSeqeunce_dict.txt".format(t=NameForThisSim), "w")
n = text_file.write(json.dumps(exons_to_use_dict))
text_file.close()

#setting the number of generations 
number_gens = int(sys.argv[2])

#SIMULATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#creating the black plotting lists 
av_mutability = []
intron_mut = []
exon_mut = []

for x in tqdm(range(number_gens)):
    
    #selecting the base to mutate 
    ab_index = base_to_mutate(DNA, weights) # ex. "mutate position 25 ("G")"   
    
    if len(DNA)-ab_index >= 3 and ab_index >= 3 : #making sure dont pick the ends 
        
        mb = base_to_turn_into(ab_index, model, DNA) # ex. "position 25 ("G") will become "A"    

        #determining if alloweed to mutate. if so, then mutating 
        if codingYesNo(ab_index, mb, DNA, DNA_map) == None:
            mDNA = mutate_chromosome(DNA, ab_index, mb) # always mutate 
            DNA = mDNA # the mutated DNA is the new DNA 
        elif codingYesNo(ab_index, mb, DNA, DNA_map) != None and codingYesNo(ab_index, mb, DNA, DNA_map) != 3: # if the mutation is in a CODING REGION...
            a_m_codons = codingYesNo(ab_index, mb, DNA, DNA_map) # check to see if the change is acceptable
            if blosum(a_m_codons, aa_conversion_dict, blosum_dict) == True: # if the change is acceptable...
                mDNA = mutate_chromosome(DNA, ab_index, mb) # mutate it (otherwise do nothing)
                DNA = mDNA # the mutated DNA is the new DNA    

        #updating the weights lits 
        weights = current_frequencies(DNA, model, ab_index, weights)

        #updating the mutability values 
        av_mutability.append(average_mutability(weights)) # AVERAGE MUTABILITY OF WHOLE STRING
        exon_mut.append(exon_mutability(DNA_map, weights)) # AVEARAGE MUTABILITY OF EXON (DNA_map, index, current_weights, total_intron_length)
        intron_mut.append(intron_mutability(DNA_map, weights)) # AVERAGE MUTABILITY OF INTRON  

# WRITING DATA TO FILES 

file = open("Scripting_data_trial2/{t}_avMutability_cds_invariant50.txt".format(t=NameForThisSim), "w+")
file.write(str(av_mutability))
file.close()

file = open("Scripting_data_trial2/{t}_intronMut_cds_invariant50.txt".format(t=NameForThisSim), "w+")
file.write(str(intron_mut))
file.close()

file = open("Scripting_data_trial2/{t}_exonMut_cds_invariant50.txt".format(t=NameForThisSim), "w+")
file.write(str(exon_mut))
file.close()

file = open("Scripting_data_trial2/{t}_DNA_final_cds_invariant50.txt".format(t=NameForThisSim), "w+")
file.write(DNA)
file.close()