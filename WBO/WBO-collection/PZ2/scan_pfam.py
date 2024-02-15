# Zadanie zaliczeniowe 2

# Tym razem nasze zadanie polegać będzie na napisaniu programu, który dla zadanej listy fragmentów białek uzupełnia
# je do pełnych sekwencji używając programu BLAST, a następnie wyszukuje w tych sekwencjach domen przy pomocy
# algorytmu HMMER, oraz wyszukuje motywów DNA w promotorach tych sekwencji (w grupach A i B).

# (6 pkt) Napisz program extend.py, który dla zadanej listy fragmentów białek w formacie fasta (
# protein_fragments.fa), znajdzie w bazie lokalnej stworzonej z podanego większego pliku fasta (genes_e_coli.fa)
# najbliższe białko i zwróci nowy plik fasta z wynikowymi sekwencjami genów już przetłumaczonych na białka. (jak
# uruchamiać blast lokalnie). Proszę zwrócić uwagę, że mamy fragmenty białek i sekwencje DNA genów.

# 2.   (7 pkt) Napisz program scan_pfam.py, który na podstawie pliku z białkami wykonuje zapytanie do serwera hmmscan
# (serwer online tu, opis api) i pobiera pliki wynikowe w formacie tsv. Korzystając z plików pobranych z serwera
# hmmer poda nam w wyniku plik csv, w którym będziemy mieli w wierszach kolejne identyfikatory białek z pliku FASTA,
# zaś w kolumnach będzie miał kolejne identyfikatory domen białkowych PFAM. Na przecięciu wiersza i kolumny stawiamy
# 0, jeśli dana domena nie została znaleziona w danym białku, a 1 w przeciwnym wypadku. Niestety moduły do parsowania
# wyjścia z HMMera w biopythonie często sprawiają problemy, dlatego polecam ręcznie wczytywać pliki tsv.


import argparse
import csv
import os
import re
import xml.etree.ElementTree as ET

import requests
from Bio import SeqIO

# Extract uuid from xml like this
# <data name="results" algo="hmmscan" uuid="8D822B56-FD45-11ED-9657-88A2B569042E">
uuid_regex = r'uuid="([0-9A-F]{8}-([0-9A-F]{4}-){3}[0-9A-F]{12})"'


def extract_uuid_by_hand(xml_text):
    return re.findall(uuid_regex, xml_text)[0][0]


def query_hmmscan(filename):
    seq = open(filename, "r").read()
    url = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
    headers = {"Expect": "", "Accept": "text/xml"}
    parameters = {"hmmdb": "pfam", "seq": seq}

    response = requests.post(url, headers=headers, data=parameters)
    try:
        uuid = ET.fromstring(response.text).find("data").attrib["uuid"]
    except ET.ParseError:
        uuid = extract_uuid_by_hand(response.text)
    # tsv is not supported by the API :((
    # I could send data to the API as a file, but the wait time is about 9 days
    result_url = f"https://www.ebi.ac.uk/Tools/hmmer/results/{uuid}/score?output=json"
    result = requests.get(result_url).json()
    return [hit["acc"] for hit in result["results"]["hits"]]


def query_proteins_to_hmmscan(filename, output):
    proteins = list(SeqIO.parse(filename, "fasta"))
    temp_file = "query.fa"

    def query(protein):
        with open(temp_file, "w") as output_handle:
            SeqIO.write(protein, output_handle, "fasta")
        return query_hmmscan(temp_file)

    hits_dict = {protein.id: query(protein) for protein in proteins}
    csv_writer(output, hits_dict)
    os.remove(temp_file)


def csv_writer(filename, data):
    all_hits = set()
    for hits in data.values():
        all_hits.update(hits)
    all_hits = sorted(all_hits)
    header = ["protein"] + all_hits
    with open(filename, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for protein, hits in data.items():
            row = [1 if hit in hits else 0 for hit in all_hits]
            writer.writerow([protein] + row)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--input", help="Input file", default="translated_proteins.fa"
    )
    parser.add_argument("-o", "--output", help="Output file", default="pfam_hits.csv")
    args = parser.parse_args()
    query_proteins_to_hmmscan(args.input, args.output)


if __name__ == "__main__":
    main()
