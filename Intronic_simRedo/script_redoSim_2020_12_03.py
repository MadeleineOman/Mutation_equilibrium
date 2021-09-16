#command settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# trialN(ame) = (sys.argv[1])
# dna_length = int(sys.argv[2])
# prop_muts =int(sys.argv[3])



#imports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy.random import choice
import copy 
from tqdm import tqdm
import sys 
from datetime import datetime
import time
import timeit
import numpy as np 
import json


#defining funtions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def base_to_mutate(DNA, current_weights, indices): 
    """
    (DNA: str, weights: Dict) -> int
    
    Will return the index of a random choice of a base that will be mutated based on probabilities given by the
    weights dictionary (i.e. the model)
    ex. "ATCGTA" --> index 3 ("G") will mutate
    """

    # NORMALIZE THE POPULATION OF WEIGHTS    
    total_freq = sum(current_weights) - current_weights[0] - current_weights[-1] # remove the start and end weight
    normalized_weights = ["error"]*len(current_weights)
    for index, value in enumerate(current_weights):
        normalized_weights[index] = value/total_freq 
    
    normalized_weights[0] = 0
    normalized_weights[-1] = 0
    # DRAW THE INDEX OF THE BASE THAT WILL BE MUTATED
    base_index = choice(indices, p=normalized_weights) 
    

    # RETURN THE INDEX
    return base_index 


#command settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
trialN = (sys.argv[1])
dna_length = int(sys.argv[2])
prop_muts =int(sys.argv[3])

#informing the graph 
sim_details = ""

#making the dna ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DNA  = choice(["A", "T", "C", "G"], size = dna_length)
DNA = "".join(DNA)
DNA_initial = copy.copy(DNA)
sim_details += str(dna_length/1000000) + "MB random seq \n "

#triplets
triplets = []
for i_1 in ["A", "T", "G", "C"]: 
    for i_2 in ["A", "T", "G", "C"]: 
        for i_3 in ["A", "T", "G", "C"]: 
            triplets.append(i_1+i_2+i_3)


#mutability model 
model = json.load(open("../Human_mutability_model/Model_2020_12_02_genomeWide.txt"))
#triplet-count dict 
triplet_chosen_count_dict = {}
triplet_into_count_dict = {}
for triplet in triplets: 
    triplet_chosen_count_dict[triplet] = [0,0,0]  
    triplet_into_count_dict[triplet] = [0,0,0]

#making the weights 
current_weights = []
for i in range(1,len(DNA)-1): 
    triplet = DNA[i-1:i+2]
    current_weights.append(model[triplet][0])

indices = [i for i in range(1, len(DNA)-1)]
mut_indices = []

#choosing number of muts 
sim_details += str(prop_muts)+"X prop muts \n"


#runnign the sim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for i in tqdm(range(int(len(DNA)*prop_muts))):
#     s0 = time.time()
    base_index = base_to_mutate(DNA, current_weights, indices)
    
#     s1 = time.time()
#     t1 = time.time() - s0
    
#     curr_index = indices.index(base_index)
#     current_weights = current_weights[0:curr_index -2]+current_weights[curr_index+3:]
    
    mut_indices.append(base_index)
#     s2 = time.time()
#     t2 = time.time() - s1
    #adding the count for "chosen to mutate" in coutns dict 
    c_triplet = DNA[base_index-1: base_index+2]
    c_triplet_left = DNA[base_index-2: base_index+1]
    c_triplet_right = DNA[base_index: base_index+3]
#     s3 = time.time()
#     t3 = time.time() - s2
    
    triplet_chosen_count_dict[c_triplet_left][0] += 1
    triplet_chosen_count_dict[c_triplet][1] += 1
    triplet_chosen_count_dict[c_triplet_right][2] += 1
    
    
    #accurate model INTO probability ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~OR ^ 
    curr_into_bases = []
    curr_into_bases_probs = []
    for into_base, prob in model[c_triplet][1].items(): 
        curr_into_bases.append(into_base)
        curr_into_bases_probs.append(prob)
    mb =  choice(curr_into_bases, p = curr_into_bases_probs)[1]

    #adding the count for "mutated into thiss" in the triplet dict 
    m_triplet = c_triplet[0]+mb+c_triplet[2]
    m_triplet_left = c_triplet_left[0:2]+mb
    m_triplet_right = mb+c_triplet_right[1:3]
    
    triplet_into_count_dict[m_triplet_left][0] += 1 
    triplet_into_count_dict[m_triplet][1] += 1 
    triplet_into_count_dict[m_triplet_right][2] += 1 

    #chanaging the dna 
    DNA = DNA[:base_index]+mb+DNA[base_index+1:]

    #updating the weights 
    current_weights[base_index-2]= model[m_triplet_left][0]    
    current_weights[base_index-1]= model[m_triplet][0]
    current_weights[base_index] = model[m_triplet_right][0]
    



text_file = open("data/{t}_triplet_into_count_dict.txt".format(t=trialN), "w")
n = text_file.write(json.dumps(triplet_into_count_dict))
text_file.close()

text_file = open("data/{t}_triplet_chosen_count_dict.txt".format(t=trialN), "w")
n = text_file.write(json.dumps(triplet_chosen_count_dict))
text_file.close()

file = open("data/{t}_DNA.txt".format(t=trialN), "w")
file.write(str(DNA))
file.close()

