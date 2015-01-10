import Bio.motifs as motifs

pwm_id_list = []
pwm_matrix_list = []
pfm3 = {'A': [],
        'C': [],
        'G': [],
        'T': []}

with open("../samples/PWMs.txt", 'r') as pwms_plain_text:
    for line in pwms_plain_text:
        line_list = line.strip().split()
        if len(line_list) == 1:
            pwm_id_list.append(line_list[0][2:])
            pwm_matrix_list.append(pfm3)
            pfm3 = {'A': [],
                    'C': [],
                    'G': [],
                    'T': []}
        else:
            if len(line_list) > 1 and line_list[1].isdigit():
                pfm3['A'].append(int(line_list[1]))
                pfm3['C'].append(int(line_list[2]))
                pfm3['G'].append(int(line_list[3]))
                pfm3['T'].append(int(line_list[4]))
pwm_matrix_list.pop(0)
pwm_matrix_list.append(pfm3)

motif_list = []
for matrix in pwm_matrix_list:
    motif_list.append(motifs.Motif(counts=matrix))

with open("../samples/PWMs_TRANSFAC.txt", 'w') as motif_file:
    motif_file.write("VV  EXAMPLE January 15, 2013\nXX\n//\n")
    for i, motif_id in enumerate(pwm_id_list):
        motif_file.write("ID  " + motif_id + "\n")
        motif_file.write(motif_list[i].format("transfac"))
