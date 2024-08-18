# Reverse translation of amino acid sequence to nucleotide sequence
def reverse_translate(aa_sequence):
    codon_table = {
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'C': ['TGT', 'TGC'],
        'D': ['GAC', 'GAT'],
        'E': ['GAA', 'GAG'],
        'F': ['TTC', 'TTT'],
        'G': ['GGA', 'GGC', 'GGG', 'GGT'],
        'H': ['CAC', 'CAT'],
        'I': ['ATA', 'ATC', 'ATT'],
        'K': ['AAA', 'AAG'],
        'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'M': ['ATG'],
        'N': ['AAC', 'AAT'],
        'P': ['CCA', 'CCC', 'CCG', 'CCT'],
        'Q': ['CAA', 'CAG'],
        'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'T': ['ACA', 'ACC', 'ACG', 'ACT'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT'],
    }

    nt_sequence = ''
    for aa in aa_sequence:
        codons = codon_table.get(aa, [''])  # Default to an empty string if amino acid not found
        nt_sequence += codons[0]  # Use the first codon (you can choose randomly if needed)

    return nt_sequence

# Example usage:
amino_acids = input("Enter the amino acid sequence: ").upper()
nucleotides = reverse_translate(amino_acids)
print(f"Nucleotide sequence: {nucleotides}")

from Bio.SeqUtils import MeltingTemp
from Bio.Seq import Seq

my_sequence = input("Enter the nucleotide sequence: ").upper()
my_seq_object = Seq(my_sequence)

# Calculate Tm using nearest neighbor method
tm_nn = MeltingTemp.Tm_NN(my_seq_object)
print(f"Melting temperature (Tm) using NN method: {tm_nn:.2f} °C")
