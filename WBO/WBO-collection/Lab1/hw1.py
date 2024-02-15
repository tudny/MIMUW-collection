import Bio.SeqIO

FILE_NAME_TEST = "data/test_fasta.fa"
FILE_NAME_TEST_TARGET = "data/test_fasta_target.fa"

FILE_NAME = "data/yeast.fa"
FILE_NAME_TARGET = "data/yeast_target.fa"


def kmers(s: list, k: int) -> set:
    return set([s[i : i + k] for i in range(len(s) - k + 1)])


def check(seqs: list, k: int):
    seqs_as_kmers = [kmers(seq, k) for seq in seqs]
    all_unique_kmers = set()
    all_kmers = set()
    for seq_kmers in seqs_as_kmers:
        all_unique_kmers = (all_unique_kmers - seq_kmers) | (seq_kmers - all_kmers)
        all_kmers = all_kmers | seq_kmers

    solutions_intersections = [seq & all_unique_kmers for seq in seqs_as_kmers]
    seqs_matches = [
        next(iter(inter)) for inter in solutions_intersections if len(inter) > 0
    ]
    return seqs_matches if len(seqs_matches) == len(seqs) else None


def binary_serach_naive(seqs: list):
    begin = 0
    end = min([len(seq) for seq in seqs]) + 1
    saved_solution = None
    while end - begin > 1:
        k = (begin + end) // 2
        if maybe_solution := check(seqs, k):
            end = k
            saved_solution = maybe_solution
        else:
            begin = k
    return end, saved_solution


def read_seqs(filename, filetype) -> (list, list):
    return zip(
        *[
            (record.reverse_complement().seq, record.id)
            for record in Bio.SeqIO.parse(filename, filetype)
        ]
    )


def write_seqs(filename, probes, seq_names):
    with open(filename, "w") as f:
        for seq, name in zip(probes, seq_names):
            f.write(f">{name}\n{seq}\n")


def unique_complement_probes_of_min_length(
    filename_source, filename_target, filetype="fasta"
):
    seqs, names = read_seqs(filename_source, filetype)
    min_k, probes = binary_serach_naive(seqs)
    if probes is None:
        print(f"No solution found.")
    else:
        print(f"Minimal k is {min_k}")
        write_seqs(filename_target, probes, names)


if __name__ == "__main__":
    unique_complement_probes_of_min_length(FILE_NAME_TEST, FILE_NAME_TEST_TARGET)
    unique_complement_probes_of_min_length(FILE_NAME, FILE_NAME_TARGET)
