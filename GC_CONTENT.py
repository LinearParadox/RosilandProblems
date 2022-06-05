import Basic_Functions

def find_gc(file):
    gc_dict = {}
    fasta_list = Basic_Functions.fasta_parse(file)
    for items in fasta_list:
        gc = 0
        for nucleotides in str(items):
            if nucleotides == "G" or nucleotides == "C":
                gc += 1
        gc_dict[items.get_id()] = gc / len(str(items))
    max_ID = max(gc_dict, key=gc_dict.get)
    max_value =max(gc_dict.values())
    return max_ID, max_value

key_val = find_gc("data/rosalind_gc (1).txt")
print(key_val[0])
print(key_val[1]*100)