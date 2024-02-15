import random

from Bio import SeqIO
from Bio import Entrez
from Bio.Blast.NCBIWWW import qblast
from Bio.Blast import NCBIXML

# 0. Wczytaj plik w formacie FastQ microbial_reads.fastq (jest już na naszym serwerze jupyter w katalogu WBO,
# nie trzeba go tam ładować), przy pomocy SeqIO.parse("..","fastq"). Są to odczyty z mikrobiomu jelitowego myszy.


microbial = SeqIO.parse("data/microbial_reads.fastq", "fastq")

# 1. Wybierz kilka losowych, dość długich sekwencji DNA i uruchom dla nich program BLAST online, obejrzyj wyniki (
# jeśli nic się nie “trafiło”, możesz wybrać inne, dłuższe sekwencje, np. powyżej 200 par zasad)

some_random_seqs = random.choices([seq for seq in microbial if len(seq) > 100], k=10)

print(some_random_seqs)

# 2. Wykonaj wyszukiwanie dla tych samych sekwencji przy pomocy interfejsu API biopython’a do blasta online (NCBIWWW)
# i parsera xml (NCBIXML). Pamiętaj o podaniu swojego adresu e-mail:

SeqIO.write(some_random_seqs, "data/random_seqs.fasta", "fasta")

# Get content of file with random sequences
with open("data/random_seqs.fasta", "r") as f:
    some_random_seqs_fasta = "".join(f.readlines())

    Entrez.email = "a.tudruj@student.uw.edu.pl"

    result_handle = qblast("blastn", "nt", some_random_seqs_fasta)

    blast_records = NCBIXML.read(result_handle)

    # Get the record with the lowest e-value
    record = min(blast_records.alignments, key=lambda x: x.hsps[0].expect)

    print(record)
