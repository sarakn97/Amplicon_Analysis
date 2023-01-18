from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """

    remainder = len(sequence) % 3

    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


seq_records = SeqIO.parse("/mnt/Active/Priya_Seq/mixed_v_infection/d12X4/S15_2.fa", "fasta")

amino_acids3 = []
# ids = []

for record in seq_records:
    record3 = record.seq[2:]
    record3trs = record3.transcribe()
    RecID = record.id
    # ids.append(RecID)
    record_final = SeqRecord(
        pad_seq(record3trs).translate(),
        id=RecID,
        name="S18",
        description="S15 envelope",
    )
    amino_acids3.append(record_final)

SeqIO.write(amino_acids3, "AA_S15.fasta", "fasta")
