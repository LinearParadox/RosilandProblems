import Basic_Functions

str1 = input()
str2 = input()
consensus_matrix = Basic_Functions.consensus_dict("rosalind_cons (2).txt")
consensus_sequence = Basic_Functions.consensus_string(consensus_matrix)
write_file = open("consens_seq.txt", "w")
write_file.write(consensus_sequence + "\n")
for keys in consensus_matrix.items():
    write_file.write(keys[0]+": ")
    for nums in keys[1]:
        write_file.write(str(nums) + " ")
    write_file.write("\n")
write_file.close()

