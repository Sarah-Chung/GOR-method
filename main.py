import source as src

DSSP_PATH = "./M8_train_dataset"
SS_STRING = "HE-"
AA_STRING = "ARNDCEQGHILKMFPSTWYV"

# General workflow to obtain 20x17 tables with information scores
def workflow(folder_path, type):
    aa, ss, N = src.init_dataset(folder_path, type)
    aa_dictionary, aa_string = src.aa_init()
    ss_dictionary, ss_string = src.ss_init()
    frequency_per_aa = src.aa_frequency_train(aa, aa_dictionary)
    frequency_per_ss = src.ss_frequency_train(ss, ss_dictionary)
    results = src.window(aa, ss, N, aa_string, ss_string, frequency_per_aa, frequency_per_ss)
    return aa, ss, results

# Validation process
def validate(validation_results, aa, ss, ss_string, aa_string):
    results = {}
    predicted_ss = ""
    i = 0
    for table in validation_results.values():
        struct = ss_string[i]
        if struct == "-":
            struct = "O"
        summed_table_id = struct + "_sum"
        summed_tables = [sum(x) for x in table]
        results[summed_table_id] = summed_tables
        i+=1
    
    H_table = results["H_sum"]
    E_table = results["E_sum"]
    O_table = results["O_sum"]

    for acid in aa:
        acid_idx = aa_string.index(acid)
        max_score = max(H_table[acid_idx], E_table[acid_idx], O_table[acid_idx])
        if max_score == H_table[acid_idx]:
            predicted_ss += "H"
        elif max_score == E_table[acid_idx]:
            predicted_ss += "E"
        elif max_score == O_table[acid_idx]:
            predicted_ss += "-"

    # Prepare accuracy calculation
    total_ss = len(ss)
    match = 0
    for char1, char2 in zip(ss, predicted_ss):
        if char1 == char2:
            match+=1
    
    # Final accuracy calculation
    accuracy_raw = match/total_ss
    accuracy = round(accuracy_raw*100,2)
    return match, accuracy


if __name__ == "__main__":
    # TRAINING
    training_results = workflow(DSSP_PATH,1)

    # VALIDATION
    aa, ss, validation_results = workflow(DSSP_PATH,2)
    matches, accuracy = validate(validation_results, aa, ss, SS_STRING, AA_STRING)

    # VALIDATION SCORE
    print("Accuracy: " + str(accuracy) + "%")
