import Bio.SeqIO
from itertools import permutations

seqs = Bio.SeqIO.parse("data/test_fasta.fa", "fasta")


def kmers(s, k):
    return list(set([s[i : i + k] for i in range(len(s) - k + 1)]))


def make_graph_edges(seq_kmers: list):
    graph = {}

    def new(a):
        if a not in graph:
            graph[a] = []

    def insert(a, b):
        new(a)
        new(b)
        graph[a].append(b)

    for node_a in seq_kmers:
        for node_b in seq_kmers:
            if node_a[:-1] == node_b[1:]:
                insert(node_b, node_a)
    return graph


def is_hamilton_path(perm, graph):
    for a, b in zip(perm, perm[1:]):
        if b not in graph[a]:
            return False
    return True


def find_hamilton_path(graph: dict):
    nodes = graph.keys()
    for perm in permutations(nodes):
        if is_hamilton_path(perm, graph):
            return perm


def make_graph_nodes(seq_kmers: list):
    graph = {kmer[:-1]: [] for kmer in seq_kmers}
    for kmer in seq_kmers:
        graph[kmer[:-1]].append(kmer[1:])
    return graph


def try_node(node, graph, n, path=None):
    if path is None:
        path = []
    for child in graph[node]:
        child_path = path.copy().append(child)
        child_graph = graph.copy()[node].pop(child)
        new_path = try_node(child, child_graph, n, child_path)
        if len(new_path) == n:
            return new_path
    return None


def find_euler_path(graph: dict):
    nodes = graph.keys()
    n = len(nodes)
    for node in nodes:
        path = try_node(node, graph, n, None)
        if path is not None:
            return path
    return None


for seq in seqs:
    print("Sequence: ", seq.id)
    seq_kmers = kmers(seq.seq, 3)
    # edge_graph = make_graph_edges(seq_kmers)
    # hamilton = find_hamilton_path(edge_graph)
    # print(edge_graph)
    # print(hamilton)
    node_graph = make_graph_nodes(seq_kmers)
    hamilton_path = find_hamilton_path(node_graph)
    print(node_graph)
    print(hamilton_path)
