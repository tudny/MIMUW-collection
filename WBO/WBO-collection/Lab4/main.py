from Bio import Phylo
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


FILENAME1 = "data/Human_PAH_paralogues.nex"
FILENAME2 = "data/Human_H2BFS_paralogues.nex"


def read_nexus_file(filename: str):
    return AlignIO.read(filename, "nexus")


# 1. Na początek ustalmy listy sekwencji. Rozważmy sekwencje paralogiczne i ortologiczne ludzkiej hydroksylazy
# fenyloalaniny.  Załóżmy, że mamy już sekwencje aminokwasowe uliniowione globalnie.  Ze względu na czasochłonność
# procesu uliniowienia, użyjemy plików z uliniowieniami w formacie nexus:  Human_PAH_paralogues,
# Human_H2BFS_paralogues (ewentualnie Human_PAH_orthologues, ale to dość duży plik). Wczytaj te pliki przy pomocy
# metod z modułu Bio.AlignIO (pliki są już na serwerze jupyter w katalogu WBO)
def task1():
    content = read_nexus_file(FILENAME1)
    print(content)
    content = read_nexus_file(FILENAME2)
    print(content)


# 2.Wylicz macierze odległości dla tych grup sekwencji przy pomocy klasy
# Bio.Phylo.TreeConstruction.DistanceCalculator  (dla macierzy BLOSUM62, osobno dla paralogów i osobno dla ortologów
# genów PAH – to  zajmie chwilę).


def calculate_matrix(file_content):
    calculator = DistanceCalculator("blosum62")
    return calculator.get_distance(file_content)


def task2():
    para_content = read_nexus_file(FILENAME1)
    orto_content = read_nexus_file(FILENAME2)

    para_matrix = calculate_matrix(para_content)
    orto_matrix = calculate_matrix(orto_content)

    print(para_matrix)
    print(orto_matrix)


# 3. Stwórz drzewa filogenetyczne na podstawie macierzy przy pomocy klasy
# Bio.Phylo.TreeConstruction.DistanceTreeConstructor zarówno metodą UPGMA – hierarchiczną jak i nj (neighbor joining)
def generate_tree(matrix):
    constructor = DistanceTreeConstructor()
    return constructor.upgma(matrix), constructor.nj(matrix)


def task3():
    tree1_upgma, tree1_nj = generate_tree(calculate_matrix(read_nexus_file(FILENAME1)))
    tree2_upgma, tree2_nj = generate_tree(calculate_matrix(read_nexus_file(FILENAME2)))
    print(tree1_upgma, tree1_nj)
    print(tree2_upgma, tree2_nj)
    return tree1_upgma, tree1_nj, tree2_upgma, tree2_nj


# 4. Wyświetl uzyskane drzewa przy pomocy metody draw_ascii() i draw()
def print_tree(tree):
    Phylo.draw_ascii(tree)
    Phylo.draw(tree)


def task4():
    tree1_upgma, tree1_nj, tree2_upgma, tree2_nj = task3()
    print_tree(tree1_upgma)
    print_tree(tree1_nj)
    print_tree(tree2_upgma)
    print_tree(tree2_nj)


# 5. Zapisz uzyskane drzewa do formatów newick i phyloxml. obejrzyj wyniki.
def task5():
    tree1_upgma, tree1_nj, tree2_upgma, tree2_nj = task3()
    Phylo.write(tree1_upgma, "data/tree1_upgma.nwk", "newick")
    Phylo.write(tree1_nj, "data/tree1_nj.nwk", "newick")
    Phylo.write(tree2_upgma, "data/tree2_upgma.nwk", "newick")
    Phylo.write(tree2_nj, "data/tree2_nj.nwk", "newick")

    Phylo.write(tree1_upgma, "data/tree1_upgma.xml", "phyloxml")
    Phylo.write(tree1_nj, "data/tree1_nj.xml", "phyloxml")
    Phylo.write(tree2_upgma, "data/tree2_upgma.xml", "phyloxml")
    Phylo.write(tree2_nj, "data/tree2_nj.xml", "phyloxml")


def main():
    task1()
    task2()
    task3()
    task4()
    task5()


if __name__ == "__main__":
    main()
