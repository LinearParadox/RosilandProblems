

def reverse_comp(sequence):
    nucleotide_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    sequence = sequence[::-1]
    new_sequence = ""
    for nucleotides in sequence:
        new_sequence += nucleotide_dict.get(nucleotides)
    return new_sequence

def count_nucleotides(sequence):
    nucleotide_dict = {"A": 0, "T": 0, "G": 0, "C": 0}
    for nucleotides in sequence:
        nucleotide_dict[nucleotides] = nucleotide_dict[nucleotides]+1
    return nucleotide_dict

def to_rna(sequence):
    return str(sequence).replace("T", "U")

def parse_fasta(fastas):
    id = ""
    sequence = ""
    fastas_list = []
    for lines in fastas:
        if lines[0] == '>':
            if sequence != "" or id != "":
                fastas_list.append(FASTASequence(id, sequence))
                sequence = ""
            id = lines[1:]
        else:
            sequence = sequence + lines



class FASTASequence:
    def __init__(self, sequence: str, id: str):
        self.fasta_sequence = sequence
        self.id = id



