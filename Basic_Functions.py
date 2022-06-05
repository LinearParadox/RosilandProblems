amino_acid_dict = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
                   "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
                   "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
                   "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
                   "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
                   "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                   "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
                   "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
                   "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
                   "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
                   "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
                   "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
                   "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                   "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
                   "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
                   "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }


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
        nucleotide_dict[nucleotides] = nucleotide_dict[nucleotides] + 1
    return nucleotide_dict


def hamming_distance(string1, string2):
    dist = 0
    for letters in range(0, len(string1)):
        if string1[letters] != string2[letters]:
            dist += 1
    return dist


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
        seq_set = set(seq)
        if seq_set == {"A", "G", "C", "U"} == set(seq) or seq_set == set(amino_acid_dict.values()) or \
                seq_set == {"A", "G", "C", "T"}:
            fasta_list.append(FASTASequence(seq, ident))

    return fasta_list


def find_substrings(string, substring, start=0):
    location_list = []
    last_pos = string.find(substring, start)
    if last_pos != -1:
        location_list.append(last_pos + 1)
        location_list = location_list + find_substrings(string, substring, last_pos + 1)
    return location_list


def consensus_dict(fasta_list):
    sequences = fasta_parse(fasta_list)
    sequence_matrix = []
    for fastas in sequences:
        sequence_matrix.append(list(str(fastas)))
    consensus_dict = {"A": ([0] * len(sequence_matrix[0])),
                      "C": [0] * len(sequence_matrix[0]),
                      "G": [0] * len(sequence_matrix[0]),
                      "T": [0] * len(sequence_matrix[0])
                      }
    for pos in range(0, len(consensus_dict["C"])):
        for elements in sequence_matrix:
            consensus_dict[elements[pos]][pos] += 1
    return consensus_dict


def consensus_string(consensus_dict):
    consensus = ""
    for pos in range(0, len(consensus_dict["C"])):
        consensus = consensus + max(consensus_dict, key=lambda k: consensus_dict[k][pos])
    return consensus


class FASTASequence:
    def __init__(self, fasta_sequence: str, fasta_id: str):
        self.fasta_sequence = fasta_sequence
        self.fasta_id = fasta_id

    def __str__(self):
        return self.fasta_sequence

    def __reversed__(self):
        return reverse_comp(self.fasta_sequence)

    def get_id(self):
        return self.fasta_id

    def to_protein(self):
        protein = ""
        protein_list = []
        start = self.fasta_sequence.find("AUG")
        for index in range(start, len(self.fasta_sequence), 3):
            protein_list.append(self.fasta_sequence[index:index + 3])
        for items in protein_list:
            if amino_acid_dict[items] == "STOP":
                break
            protein = protein + amino_acid_dict[items]
        return protein

    def to_rna(self):
        return str(self.fasta_sequence).replace("T", "U")

    def to_dna(self):
        return str(self.fasta_sequence).replace("U", "T")
