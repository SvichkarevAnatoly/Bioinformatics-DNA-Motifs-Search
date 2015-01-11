import Bio.motifs as motifs
from Bio.Alphabet import IUPAC

dna_alf = sorted(list(IUPAC.unambiguous_dna.letters))

pwm_id_list = []
pwm_matrix_list = []
pwm_matrix = {n: [] for n in dna_alf}

with open("../samples/PWMs.txt", 'r') as pwms_plain_text:
    for line in pwms_plain_text:
        line_list = line.strip().split()
        if len(line_list) == 1:
            pwm_id_list.append(line_list[0][2:])
            pwm_matrix_list.append(pwm_matrix)
            pwm_matrix = {n: [] for n in dna_alf}
        else:
            if len(line_list) > 1 and line_list[1].isdigit():
                for i, nucleotide in enumerate(dna_alf):
                    pwm_matrix[nucleotide].append(int(line_list[i]))
pwm_matrix_list.pop(0)
pwm_matrix_list.append(pwm_matrix)

motif_list = []
for matrix in pwm_matrix_list:
    motif_list.append(motifs.Motif(alphabet=IUPAC.unambiguous_dna, counts=matrix))

with open("../samples/PWMs_TRANSFAC.txt", 'w') as motif_file:
    motif_file.write("VV  EXAMPLE January 15, 2013\nXX\n//\n")
    for i, motif_id in enumerate(pwm_id_list):
        motif_file.write("ID  " + motif_id + "\n")
        motif_file.write(motif_list[i].format("transfac"))
