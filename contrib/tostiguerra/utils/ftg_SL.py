import numpy as np


def dH_dS(seq1, seq2):
    dH = 0
    dS = 0
    weight = 0.5

    # print(seq1, seq2)
    for k in range(len(seq1) - 1):
        if seq1[k : k + 2] + seq2[::-1][k : k + 2] in dH_stack.keys():
            if k == 0 or k == len(seq1) - 2:
                weight = 0.5
                dH += weight * dH_stack[seq1[k : k + 2] + seq2[::-1][k : k + 2]]
                dS += weight * dS_stack[seq1[k : k + 2] + seq2[::-1][k : k + 2]]

            else:
                weight = 1
                dH += dH_stack[seq1[k : k + 2] + seq2[::-1][k : k + 2]]
                dS += dS_stack[seq1[k : k + 2] + seq2[::-1][k : k + 2]]

        else:
            pass

    return dH, dS


def sliding_window(seqs, i, j):
    slide_time = len(seqs[i]) - len(seqs[j])
    HS = []
    seq_j = seqs[j + 1][0] + seqs[j][::-1] + seqs[j - 1][-1]
    for sl in range(slide_time + 1):
        if sl == 0:
            if sl == slide_time:
                HS.append(
                    dH_dS(seqs[i - 1][-1] + seqs[i] + seqs[i + 1][0], seq_j[::-1])
                )
            else:
                HS.append(
                    dH_dS(seqs[i - 1][-1] + seqs[i][: len(seqs[j]) + 1], seq_j[::-1])
                )

        elif sl == slide_time:
            HS.append(
                dH_dS(seqs[i][sl - 1 : sl + len(seqs[j])] + seqs[i + 1][0], seq_j[::-1])
            )

        else:
            HS.append(dH_dS(seqs[i][sl - 1 : sl + len(seqs[j]) + 1], seq_j[::-1]))

    HS.sort(key=lambda x: x[0])

    return HS[0][0], HS[0][1]


def ordering_seqs(seq1, seq2, idx1, idx2):
    if len(seq1) >= len(seq2):
        return idx1, idx2
    else:
        return idx2, idx1


def interaction_matrix(seqs, material, salt_conc):
    material = material.upper()
    if material == "DNA":
        seqs = [seqs[i].upper().replace("U", "T") for i in range(len(seqs))]
        from DNA_SL import dH_stack, dS_stack
    elif material == "RNA":
        seqs = [seqs[i].upper().replace("T", "U") for i in range(len(seqs))]
        from RNA_22 import dH_stack, dS_stack
    else:
        print("Error: material option must be specified! Choose between RNA and DNA")
        exit(0)

    global dH_stack
    global dS_stack

    dH = np.zeros((len(seqs), len(seqs)))
    dS = np.zeros((len(seqs), len(seqs)))

    for i in range(1, len(seqs) - 1):
        for j in range(i + 1, len(seqs) - 1):

            longer, shorter = ordering_seqs(seqs[i], seqs[j], i, j)

            H, S = sliding_window(seqs, longer, shorter)
            dH[i][j] = H
            dS[i][j] = S

            # salt_correction = 0.368 * (len(seqs[i]) - 1.0)  * math.log(salt_conc)
            # dS[i][j] += salt_correction

    return dH, dS
