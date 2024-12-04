import os
import math

# Read .dssp files and append ss and aa-sequence to string
def init_dataset(folder_path, type):
    aa = ""
    ss = ""
    files = os.listdir(folder_path)
    num_training = 1990 # 1990 training files, 10 validation files
    if type == 1:
        sel_files = files[:len(files)-10]
        print("Number of training files: " + str(len(sel_files)))
    elif type == 2:
        sel_files = files[len(files)-10:]
        print("Number of validation files: " + str(len(sel_files)))
    for file in sel_files:
        file_path = os.path.join(folder_path, file)
        with open(file_path, "r") as file:
            lines = file.readlines()
            for i in range(len(lines)):
                line = str(lines[i])
                line = line.strip()
                if (i%2 == 0):
                    ss += line
                else:
                    aa += line
    filtered_aa = ""
    filtered_ss = ""
    for char1, char2 in zip(aa, ss):
        if char1 != 'X':
            filtered_aa += char1
            filtered_ss += char2
    N = len(filtered_aa)
    return filtered_aa, filtered_ss, N

# Initialize empty AA dictionary and string of AA for indices
def aa_init():
    aa_dictionary = {
    'A': 0, 'R': 0, 'N': 0, 'D': 0, 'C': 0,
    'E': 0, 'Q': 0, 'G': 0, 'H': 0, 'I': 0,
    'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0,
    'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
    }
    aa_string = "ARNDCEQGHILKMFPSTWYV"
    return aa_dictionary, aa_string

def ss_init():
    ss_dictionary = {
        'H':0, 'E':0, '-':0
    }
    ss_string = "HE-"
    return ss_dictionary, ss_string

# Obtain total number of amino acids + their frequency in the training dataset
def aa_frequency_train(aa, aa_dict):
    frequency_per_aa = aa_dict
    for char in aa:
        if char != 'X':
            frequency_per_aa[char]+=1
    return frequency_per_aa

# Obtain total number of secondary sequence predictions + their frequency in the training dataset
def ss_frequency_train(ss, ss_dict):
    frequency_per_ss = ss_dict
    for char in ss:
        frequency_per_ss[char]+=1
    return frequency_per_ss

# Count the total frequency of aa and ss pairs together
def ss_aa_pairs(aa,ss,N):
    pairs = {}
    for i in range(N):
        structure = ss[i]
        acid = aa[i]
        if structure == "-":
            structure = "O"
        pair = (acid + "_" + structure)
        if acid == "X":
            break
        if (pair not in pairs):
            pairs[pair] = 0
        pairs[pair] += 1
    return pairs
    
def window(aa, ss, N, aa_string, ss_string, aa_freq, ss_freq):
    # Obtain total occurrences of each AA at each position in the window
    H_table = [[0 for i in range(17)] for j in range(20)] # Alpha
    E_table = [[0 for i in range(17)] for j in range(20)] # Beta
    O_table = [[0 for i in range(17)] for j in range(20)] # Coil

    for i in range(8, N-8): # Begin at 8th index
        for j in range(-8,9):
            z = j+8
            aa_j = aa[i+j]
            ss_j = ss[i+j]
            if aa_j != "X":
                aa_idx = aa_string.index(aa_j)
                if ss_j == "H":
                    H_table[aa_idx][z] +=1
                elif ss_j == "E":
                    E_table[aa_idx][z] += 1
                elif ss_j == "-":
                    O_table[aa_idx][z] += 1
            else:
                break

    tables = {
        "H": H_table,
        "E": E_table,
        "O": O_table
    }

    # Storing final tables
    results = {}
    for i in range(20):
        a = aa_freq[aa_string[i]]
        for j in range(3):
            struct = ss_string[j]
            if struct == "-":
                struct = "O"
            table = tables[struct]
            final_table = struct + "_table_final"
            
            # Determine values for each variable to calculate log odds:
                # log(CN/AS)
                # C = Count number of AA with SS
                # N = Total number of AA
                # A = Count number of specific AA in db
                # S = Count number of specific structure in db

            s = ss_freq[ss_string[j]]
            final = [[(math.log((x*N)/(a*s))) for x in row] for row in table]
            results[final_table] = final
    return results

# Print formatted tables
def print_table(table):
    for row in table:
        print(row)
