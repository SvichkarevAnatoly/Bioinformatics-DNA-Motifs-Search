from time import gmtime, strftime
import Bio.motifs as motifs
from Bio.Alphabet import IUPAC
import lib

dna_alf = lib.dna_alphabet()

working_path_str = "../test_samples/"
with open(working_path_str + "PWMs.txt", 'r') as pwms_plain_text:
    pwm_id_list, pwm_matrix_list, pwm_matrix = lib.read_excel_motif_matrix_list_from_file(pwms_plain_text)

motif_list = []
for matrix in pwm_matrix_list:
    motif_list.append(motifs.Motif(alphabet=IUPAC.unambiguous_dna, counts=matrix))

now = gmtime()
time_str = strftime("%B %d, %Y %X", now)
with open(working_path_str + "PWMs_TRANSFAC.txt", 'w') as motif_file:
    motif_file.write("VV  " + time_str + "\nXX\n//\n")
    for motif_id, motif_matrix in zip(pwm_id_list, motif_list):
        motif_file.write("ID  " + motif_id + "\n")
        motif_file.write(motif_matrix.format("transfac"))
