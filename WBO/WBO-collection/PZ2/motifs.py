import argparse
import random
import re

import Bio.motifs
import Bio.Seq
from Bio import SeqIO


# ======================================================================================================================
# Consensus method implementation
# ======================================================================================================================


class MotifData:
    def __init__(self, motif, value, idx):
        self.motif = motif
        self.value = value
        self.id = idx

    def __gt__(self, other):
        return self.value > other.value

    @staticmethod
    def of(motif, value):
        return MotifData(motif, value, None)

    def __str__(self):
        return f"Motif:\n{str(self.motif)}\nValue: {self.value}\n"


class BestKeeper:
    def __init__(self, size):
        self.size = size
        self.data = []

    def normalize(self):
        self.data.sort(reverse=True)
        self.data = self.data[: self.size]

    def add(self, motif_data):
        self.add_all([motif_data])

    def add_all(self, motif_datas: []):
        self.data.extend(motif_datas)
        self.normalize()

    def __iadd__(self, other: []):
        self.add_all(other)
        return self

    def into_list(self):
        return self.data


def sliding_window(elements, window_size):
    for i in range(len(elements) - window_size + 1):
        yield elements[i : i + window_size]


def without_ith_element(elements, i):
    return elements[:i] + elements[i + 1 :]


def consensus_with_unified_motif(sequences, motif_length, selected_motif):
    best = MotifData.of(selected_motif, 0)

    for seq_id, seq in enumerate(sequences):
        for window in sliding_window(seq, motif_length):
            motif = Bio.motifs.create([window] + selected_motif.instances)
            pwm = motif.counts.normalize(
                pseudocounts={"A": 0.6, "C": 0.4, "G": 0.4, "T": 0.6}
            )
            value = pwm.log_odds().std({"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3})
            if (current := MotifData(motif, value, seq_id)) > best:
                best = current

    return (
        consensus_with_unified_motif(
            without_ith_element(sequences, best.id), motif_length, best.motif
        )
        if len(sequences) > 1
        else best
    )


def consensus_single(sequences, motif_length, keep=5):
    motifs = BestKeeper(keep)

    first, *rest = sequences
    for window in sliding_window(first, motif_length):
        motif = Bio.motifs.create([window])
        motifs += [consensus_with_unified_motif(rest, motif_length, motif)]

    return motifs.into_list()


def consensus(sequences, motif_length=7, iterations=1, keep=5):
    best = BestKeeper(keep)

    for _ in range(iterations):
        permuted_sequences = list(sequences)
        random.shuffle(permuted_sequences)
        best += consensus_single(permuted_sequences, motif_length, keep=keep)

    return best.into_list()


# ======================================================================================================================
# End of consensus implementation
# ======================================================================================================================


def filter_by_group(seqs, group):
    return [seq for seq in seqs if seq.id.startswith(group)]


hidden_name = r"\(extended_from: (.*)\)"


def map_to_seq_name(seq):
    return re.search(hidden_name, seq.description).group(1)


def group_to_promoters(proteins, group, ids_to_promoters):
    group_items = filter_by_group(proteins, group)
    group_promoters_names = list(map(map_to_seq_name, group_items))
    return [ids_to_promoters[name] for name in group_promoters_names]


def find_motifs(seqs, motif_length, keep):
    return consensus(
        [str(seq.seq) for seq in seqs], motif_length=motif_length, keep=keep
    )


def do_run(extended_path, promoters_path, output, motif_length, keep):
    extended_proteins = list(SeqIO.parse(extended_path, "fasta"))
    promoters = list(SeqIO.parse(promoters_path, "fasta"))
    ids_to_promoters = {prom.id: prom for prom in promoters}

    group_A_promoters = group_to_promoters(
        extended_proteins, "groupA", ids_to_promoters
    )
    group_B_promoters = group_to_promoters(
        extended_proteins, "groupB", ids_to_promoters
    )

    print("Running consensus for group A")
    group_A_motifs = find_motifs(group_A_promoters, motif_length, keep)
    print("Running consensus for group B")
    group_B_motifs = find_motifs(group_B_promoters, motif_length, keep)

    print(f"Found {len(group_A_motifs)} motifs for group A")
    print(f"Found {len(group_B_motifs)} motifs for group B")

    with open(output, "w") as output_file:
        for motif in group_A_motifs:
            output_file.write(f">groupA\n{motif}\n")
        for motif in group_B_motifs:
            output_file.write(f">groupB\n{motif}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-e",
        "--extended",
        type=str,
        help="Extended proteins file name",
        default="translated_proteins.fa",
    )
    parser.add_argument(
        "-p",
        "--promoters",
        type=str,
        help="Promoters file name",
        default="proms_e_coli_fixed.fa",
    )
    parser.add_argument(
        "-o", "--output", type=str, help="Output file name", default="motifs.fa"
    )
    parser.add_argument(
        "-m", "--motif_length", type=int, help="Motif length", default=10
    )
    parser.add_argument("-k", "--keep", type=int, help="Keep best motifs", default=5)

    args = parser.parse_args()
    do_run(args.extended, args.promoters, args.output, args.motif_length, args.keep)


if __name__ == "__main__":
    main()
