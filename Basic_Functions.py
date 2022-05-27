

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


def fasta_parse(file_path):
    seq = ""
    fasta_list = []
    with open(file_path) as sequences:
        ident = sequences.readline()[1:].strip()
        line = sequences.readline()
        while line != "":
            if line[0] == ">":
                fasta_list.append(FASTASequence(seq, ident))
                ident = line[1:].strip()
                seq = ""
                line = sequences.readline()
            else:
                seq = seq + line.strip()
                line = sequences.readline()
    return fasta_list



class FASTASequence:
    def __init__(self, fasta_sequence: str, fasta_id: str):
        self.fasta_sequence = fasta_sequence
        self.fasta_id = fasta_id

    def __str__(self):
        return self.fasta_sequence

    def __reversed__(self):
        return reverse_comp(self.fasta_sequence)





