import os

from Bio import SeqIO, Seq, AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def solution(path, source, clean=False):
    print("Running solution for", source, "in", path)

    def save_log(identifier, _stdout, _stderr):
        with open(f"{path}/{identifier}.out.log", "w") as log:
            log.write(_stdout)
        with open(f"{path}/{identifier}.err.log", "w") as log:
            log.write(_stderr)

    def draw_phylo(tree, name):
        Phylo.draw(tree)
        print(f"Tree {name}:")
        Phylo.draw_ascii(tree)
        Phylo.write(tree, f"{path}/{name}.nwk", "newick")

    # 1.Weźmy na początek naszą rodzinę sekwencji białek histonowych. Jest dostępna na naszym serwerze jupyter w pliku
    # histones.fa . Na początek potrzebujemy wczytać je z pliku (SeqIO.parse) przetłumaczyć je na sekwencje białkowe (
    # metoda translate()) i zapisać do pliku fasta (SeqIO.write) wraz z odpowiednimi identyfikatorami i usuniętymi
    # kodonami stopu (symbol "*" na końcu, który powstaje w wyniku translate)

    # ================ 1 =================

    def translate_seq(seq: Seq) -> Seq:
        seq.seq = seq.seq.translate(stop_symbol="")
        return seq

    sequences = list(SeqIO.parse(f"{path}/{source}.fa", "fasta"))
    translated_sequences = list(map(translate_seq, sequences))
    SeqIO.write(translated_sequences, f"{path}/{source}_translated.fa", "fasta")

    # 2. Wykorzystajmy program clustalw (Bio.Align.Applications.ClustalwCommandline) do wykonania mulituliniowienia tych
    # sekwencji na naszym serwerze plik wykonywalny jest zainstalowany w katalogu /usr/bin. Po wykonaniu zadania wczytaj
    # uliniowienie z pliku .aln (AlignIO.read ("..", format="clustal") oraz drzewo filogenetyczne z pliku .dnd (
    # Phylo.read("..", format="newick")

    # ================ 2 =================

    clustalw_cline = ClustalwCommandline(
        cmd="clustalw",
        infile=f"{path}/{source}_translated.fa",
        outfile=f"{path}/{source}_translated_clustalw.aln",
    )
    stdout, stderr = clustalw_cline()
    os.rename(  # rename files to match clustalw output
        f"{path}/{source}_translated.dnd", f"{path}/{source}_translated_clustalw.dnd"
    )
    save_log("clustalw", stdout, stderr)

    alignment = AlignIO.read(
        f"{path}/{source}_translated_clustalw.aln", format="clustal"
    )
    phylo_tree = Phylo.read(f"{path}/{source}_translated_clustalw.dnd", format="newick")

    draw_phylo(phylo_tree, "clustalw")

    # 3. Wykorzystaj teraz program muscle (Bio.Align.Applications.MuscleCommandline) do wykonania multiuliniowienia tych
    # samych sekwencji. Wczytaj je z pliku w formacie fasta

    # ================ 3 =================

    muscle_cline = MuscleCommandline(
        "muscle",
        input=f"{path}/{source}_translated.fa",
        out=f"{path}/{source}_translated_muscle.aln",
    )
    stdout, stderr = muscle_cline()
    save_log("muscle", stdout, stderr)

    muscle_alignment = AlignIO.read(
        f"{path}/{source}_translated_muscle.aln", format="fasta"
    )

    # 4. Stwórz drzewo filogenetyczne na podstawie uzyskanego uliniowienia metodą neighbor joining jak na poprzednich
    # zajęciach. Porównaj to drzewo do drzewa otrzymanego z programu clustalw

    # ================ 4 =================

    calculator = DistanceCalculator("blosum62")
    muscle_distance_matrix = calculator.get_distance(muscle_alignment)
    muscle_constructor = DistanceTreeConstructor()
    muscle_tree = muscle_constructor.nj(muscle_distance_matrix)

    draw_phylo(muscle_tree, "muscle")

    print("Muscle tree:")
    print(muscle_tree)

    print("Clustalw tree:")
    print(phylo_tree)

    # Cleanup
    if clean:
        files = [
            f"{path}/clustalw.nwk",
            f"{path}/muscle.nwk",
            f"{path}/{source}_translated.fa",
            f"{path}/{source}_translated_clustalw.aln",
            f"{path}/{source}_translated_clustalw.dnd",
            f"{path}/{source}_translated_muscle.aln",
            f"{path}/clustalw.out.log",
            f"{path}/clustalw.err.log",
            f"{path}/muscle.out.log",
            f"{path}/muscle.err.log",
        ]
        for file in files:
            os.remove(file)


clean = False
solution("data/histones", "histones", clean=clean)
solution("data/orto", "Human_PAH_orthologues", clean=clean)
