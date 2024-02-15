from Bio import Phylo, SeqIO

tree = Phylo.read("tree_nj", "newick")
Phylo.draw_ascii(tree)

tree = Phylo.read("tree_upgma", "newick")
Phylo.draw_ascii(tree)

tree = Phylo.read("bialka.dnd", "newick")
Phylo.draw_ascii(tree)

bialka = list(SeqIO.parse("bialka.fa", "fasta"))
origin = SeqIO.read("PAH.fa", "fasta")


def distance(s1: str, s2: str) -> float:
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2)) / len(s1)


for seq in bialka:
    print(seq.id)
    print(distance(seq, origin))
