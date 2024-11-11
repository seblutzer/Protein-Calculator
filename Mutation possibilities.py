codons_to_amino_acids = {
    'UUU': ['Phe', 'F', 'Apolar'],
    'UUC': ['Phe', 'F', 'Apolar'],
    'UUA': ['Leu', 'L', 'Apolar'],
    'UUG': ['Leu', 'L', 'Apolar'],

    'UAU': ['Tyr', 'Y', 'Polar'],
    'UAC': ['Tyr', 'Y', 'Polar'],
    'UAA': ['Stop', '', ''],
    'UAG': ['Stop', '', ''],

    'UCU': ['Ser', 'S', 'Polar'],
    'UCC': ['Ser', 'S', 'Polar'],
    'UCA': ['Ser', 'S', 'Polar'],
    'UCG': ['Ser', 'S', 'Polar'],

    'UGU': ['Cys', 'C', 'Polar'],
    'UGC': ['Cys', 'C', 'Polar'],
    'UGA': ['Stop', '', ''],
    'UGG': ['Trp', 'W', 'Apolar'],

    'CUU': ['Leu', 'L', 'Apolar'],
    'CUC': ['Leu', 'L', 'Apolar'],
    'CUA': ['Leu', 'L', 'Apolar'],
    'CUG': ['Leu', 'L', 'Apolar'],

    'CCU': ['Pro', 'P', 'Apolar'],
    'CCC': ['Pro', 'P', 'Apolar'],
    'CCA': ['Pro', 'P', 'Apolar'],
    'CCG': ['Pro', 'P', 'Apolar'],

    'CAU': ['His', 'H', 'Polar'],
    'CAC': ['His', 'H', 'Polar'],
    'CAA': ['Gln', 'Q', 'Polar'],
    'CAG': ['Gln', 'Q', 'Polar'],

    'CGA': ['Arg', 'R', 'Polar'],
    'CGG': ['Arg', 'R', 'Polar'],
    'CGU': ['Arg', 'R', 'Polar'],
    'CGC': ['Arg', 'R', 'Polar'],

    'AUU': ['Ile', 'I', 'Apolar'],
    'AUC': ['Ile', 'I', 'Apolar'],
    'AUA': ['Ile', 'I', 'Apolar'],
    'AUG': ['Met', 'M', 'Apolar'],

    'ACU': ['Thr', 'T', 'Polar'],
    'ACC': ['Thr', 'T', 'Polar'],
    'ACA': ['Thr', 'T', 'Polar'],
    'ACG': ['Thr', 'T', 'Polar'],

    'AAU': ['Asn', 'N', 'Polar'],
    'AAC': ['Asn', 'N', 'Polar'],
    'AAA': ['Lys', 'K', 'Polar'],
    'AAG': ['Lys', 'K', 'Polar'],

    'AGU': ['Ser', 'S', 'Polar'],
    'AGC': ['Ser', 'S', 'Polar'],
    'AGA': ['Arg', 'R', 'Polar'],
    'AGG': ['Arg', 'R', 'Polar'],

    'GUU': ['Val', 'V', 'Apolar'],
    'GUC': ['Val', 'V', 'Apolar'],
    'GUA': ['Val', 'V', 'Apolar'],
    'GUG': ['Val', 'V', 'Apolar'],

    'GCU': ['Ala', 'A', 'Apolar'],
    'GCC': ['Ala', 'A', 'Apolar'],
    'GCA': ['Ala', 'A', 'Apolar'],
    'GCG': ['Ala', 'A', 'Apolar'],

    'GAU': ['Asp', 'D', 'Polar'],
    'GAC': ['Asp', 'D', 'Polar'],
    'GAA': ['Glu', 'E', 'Polar'],
    'GAG': ['Glu', 'E', 'Polar'],

    'GGU': ['Gly', 'G', 'Apolar'],
    'GGC': ['Gly', 'G', 'Apolar'],
    'GGA': ['Gly', 'G', 'Apolar'],
    'GGG': ['Gly', 'G', 'Apolar']
}


def generate_mutations(codon):
    nucleotides = ['A', 'U', 'C', 'G']
    mutations = {}

    for i in range(len(codon)):
        for nt in nucleotides:
            if nt != codon[i]:  # Evitar mutação para o mesmo nucleotídeo
                mutated_codon = codon[:i] + nt + codon[i + 1:]
                if mutated_codon in codons_to_amino_acids:
                    amino_acid = codons_to_amino_acids[mutated_codon]
                    mutations[mutated_codon] = f"{amino_acid[0]} - {amino_acid[2]}"

    return mutations


# Exemplo de uso
codon = 'GGG'
mutations = generate_mutations(codon)
print(f'O aminoácido sequenciado é: {codons_to_amino_acids[codon][0]}: {codons_to_amino_acids[codon][2]}')
for mutated_codon, description in mutations.items():
    if codons_to_amino_acids[codon][0] == description[:3]:
        print(f'{mutated_codon} não altera o aminoácido')
    else:
        print(f"{mutated_codon}: {description}")