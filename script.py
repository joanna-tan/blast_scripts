from Bio import AlignIO
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
from Bio.pairwise2 import format_alignment
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline

def blast_internet():
    input_file = "NATL2A_GCA_000012465.1_ASM1246v1_translated_cds.faa"
    seq1 = ""
    protein_sequences = SeqIO.parse(open(input_file), 'fasta')
    result_handle = NCBIWWW.qblast("blastn", "ptn", protein_sequences)

    output_file = open("output.txt", "w")
    output_file.writelines(result_handle.read)

def blast_local(query, database):
    for seq in query:
        cline = NcbiblastpCommandline(query=seq, db=database,
                                evalue=0.001, remote=True, ungapped=True, out = "optunia.xml")
        print(cline)
        # cline()
        # stdout, stderr = cline
        # print(stdout + stderr)
        # print(cline())


if __name__ == "__main__":
    print("Hello world!")
    # blast_internet()

    input_file = "NATL2A_GCA_000012465.1_ASM1246v1_translated_cds.faa"
    seq1 = ""
    query = []
    protein_sequences = SeqIO.parse(open(input_file), 'fasta')
    # print(protein_sequences)
    i = 0
    for seq in protein_sequences:
    # for i in range(10):
        # seq = protein_sequences[i]
        # print(seq.id) #--> retrives the header
        # print(seq.seq) #-->retrieves the sequence
        # print(len(seq))
        query.append(seq.seq)
        seq1 += (seq.seq)
        if i > 5:
            break
        i += 1

    seq1 = Seq(seq1)
    # print(seq1)

    input_file = "MED4_GCA_000011465.1_ASM1146v1_translated_cds.faa"
    seq2 = ""
    protein_sequences = SeqIO.parse(open(input_file), 'fasta')
    blast_local(query, "MED4_GCA_000011465.1_ASM1146v1_translated_cds.faa")
    # i = 0
    # for seq in protein_sequences:
    # # for i in range(10):
    #     # print(seq.id) #--> retrives the header
    #     # print(seq.seq) #-->retrieves the sequence
    #     # print(len(seq))
    #     # seq = protein_sequences[i]
    #     seq2 += (seq.seq)
    #     if i > 5:
    #         break
    #     i+= 1

    # seq2 = Seq(seq2)
    # print("seq1", len(seq1), seq1)
    # print("seq2", len(seq2), seq2)

    # output_file = open("out.txt", "w")
    # alignments = pairwise2.align.globalxx(seq1, seq2)
    # # print(format_alignment(*alignments[0]))
    # output_file.writelines(format_alignment(*alignments[0]))
    # for match in alignments:
    #     # output_file.writelines(format_alignment(*match))
    #     output_file.writelines(format_alignment(*alignments[0]))
    #     print(match)


    # aligner = Align.PairwiseAligner()
    # print(aligner)

    # alignment = aligner.align(seq1, seq2)
    # for align in alignment:
    #     print(align)