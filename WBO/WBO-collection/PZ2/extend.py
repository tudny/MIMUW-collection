# Zadanie zaliczeniowe 2

# Tym razem nasze zadanie polegać będzie na napisaniu programu, który dla zadanej listy fragmentów białek uzupełnia
# je do pełnych sekwencji używając programu BLAST, a następnie wyszukuje w tych sekwencjach domen przy pomocy
# algorytmu HMMER, oraz wyszukuje motywów DNA w promotorach tych sekwencji (w grupach A i B).

# (6 pkt) Napisz program extend.py, który dla zadanej listy fragmentów białek w formacie fasta (
# protein_fragments.fa), znajdzie w bazie lokalnej stworzonej z podanego większego pliku fasta (genes_e_coli.fa)
# najbliższe białko i zwróci nowy plik fasta z wynikowymi sekwencjami genów już przetłumaczonych na białka. (jak
# uruchamiać blast lokalnie). Proszę zwrócić uwagę, że mamy fragmenty białek i sekwencje DNA genów.

# Rozwiązanie:
import os
import argparse

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline


# If len(arr) % n == 0 do nothing
# Else remove (len(arr) % n) elements from arr
def trim_to_kn_length(arr, n):
    return arr[: -(len(arr) % n)] if len(arr) % n != 0 else arr


def extend_proteins(protein_fragments_file, genes_file, output_file):
    all_genes = list(SeqIO.parse(genes_file, "fasta"))
    # Make a map of gene sequences by gene name
    gene_sequences = {}
    for gene in all_genes:
        gene_sequences[gene.id] = trim_to_kn_length(gene.seq, 3).translate()

    # If output file exists, remove it
    if os.path.exists(output_file):
        print(f"Plik {output_file} już istnieje. Nadpisuję.")
        os.remove(output_file)

    # Read protein fragments
    protein_fragments = list(SeqIO.parse(protein_fragments_file, "fasta"))

    # Query local BLAST database for each protein fragment
    for protein_fragment in protein_fragments:
        # Crate a temporary query file
        query_file = "query.fasta"
        SeqIO.write(protein_fragment, query_file, "fasta")

        # Call tblastn
        blast_result_file = "blast_results.xml"
        tblastn_cline = NcbitblastnCommandline(
            query=query_file, subject=genes_file, outfmt=5, out=blast_result_file
        )
        tblastn_cline()

        # Parse results
        result_handle = open(blast_result_file)
        blast_records = NCBIXML.parse(result_handle)

        # Find the best alignment
        bast_alignment = None
        best_gene = None
        best_name = None
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                if (
                    bast_alignment is None
                    or alignment.hsps[0].score > bast_alignment.hsps[0].score
                ):
                    bast_alignment = alignment
                    best_gene = gene_sequences[alignment.hit_id]
                    best_name = alignment.hit_id

        protein_fragment.description += f" (extended_from: {best_name})"
        protein_fragment.seq = best_gene
        # Save the best gene to the output file
        with open(output_file, "a") as output_handle:
            SeqIO.write(protein_fragment, output_handle, "fasta")

        # Remove temporary files
        os.remove(query_file)
        os.remove("blast_results.xml")

    print(
        "Przetłumaczono i zapisano rozszerzone sekwencje białek w pliku:", output_file
    )


def main():
    parser = argparse.ArgumentParser(
        description="Extend protein fragments to full sequences"
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output file name",
        default="translated_proteins.fa",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        help="Protein fragments file name",
        default="protein_fragments.fa",
    )
    parser.add_argument(
        "-g",
        "--genes",
        type=str,
        help="Genes database file name",
        default="genes_e_coli_new.fa",
    )

    args = parser.parse_args()
    files_to_check = [args.input, args.genes]
    for file in files_to_check:
        if not os.path.exists(file):
            print(f"Plik {file} nie istnieje. Przerywam.")
            return
    extend_proteins(args.input, args.genes, args.output)


if __name__ == "__main__":
    main()
