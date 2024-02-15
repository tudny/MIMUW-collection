from itertools import zip_longest, islice

from Bio.SeqIO import parse


# Snippet from https://louisabraham.github.io/articles/suffix-arrays
def to_int_keys_best(l):
    """
    l: iterable of keys
    returns: a list with integer keys
    """
    seen = set()
    ls = []
    for e in l:
        if e not in seen:
            ls.append(e)
            seen.add(e)
    ls.sort()
    index = {v: i for i, v in enumerate(ls)}
    return [index[v] for v in l]


def suffix_array_best(s):
    """
    suffix array of s
    O(n * log(n)^2)
    """
    n = len(s)
    k = 1
    line = to_int_keys_best(s)
    while max(line) < n - 1:
        line = to_int_keys_best(
            [
                a * (n + 1) + b + 1
                for (a, b) in zip_longest(line, islice(line, k, None), fillvalue=-1)
            ]
        )
        k <<= 1
    return line


# End of snippet https://louisabraham.github.io/articles/suffix-arrays


def reverse_sa(sa: [int]) -> [int]:
    return [x[0] for x in sorted(enumerate(sa), key=lambda x: x[1])]


def cyclic(txt: str):
    for d in range(len(txt)):
        yield txt[d:] + txt[0:d]


def compute_BWT_naive(txt: str) -> str:
    """
    Compute Burrows-Wheeler transform of input string using naive algorithm.
    :param txt: input string
    :return: BWT(txt)
    """
    txt = txt + "~"
    all_rotations = list(cyclic(txt))
    all_rotations.sort()
    bwt = [x[-1] for x in all_rotations]
    return "".join(bwt)


def construct_suffix_array(txt: str) -> [int]:
    """
    Construct suffix array of input string.
    :param txt: input string
    :return: suffix array of txt
    """
    suffix_rev = suffix_array_best(txt)
    suffix = [x[0] for x in sorted(enumerate(suffix_rev), key=lambda x: x[1])]
    return suffix


def compute_BWT(txt: str) -> str:
    """
    Compute Burrows-Wheeler transform of input string using algorithm based on suffix array.
    :param txt: input string
    :return: BWT(txt)
    """
    txt = txt + "~"
    suffix_array = construct_suffix_array(txt)
    bwt = [txt[i - 1] for i in suffix_array]
    return "".join(bwt)


def construct_c_array(bwt: str) -> dict[str, int]:
    """
    Construct C array for input string.
    C array is a dictionary where keys are characters from input string and values are number of occurrences of
    characters that are lexically smaller than key.
    :param bwt: input string
    :return: C array for bwt
    """
    count_map = {}
    for letter in bwt:
        count_map[letter] = count_map.get(letter, 0) + 1

    alphabet = sorted(set(bwt))
    c_array = {alphabet[0]: 0}
    for i in range(1, len(alphabet)):
        c_array[alphabet[i]] = c_array[alphabet[i - 1]] + count_map[alphabet[i - 1]]

    return c_array


def construct_occ_array(bwt: str) -> dict[tuple[str, int], int]:
    alphabet = sorted(set(bwt))
    occ_dict = {(letter, -1): 0 for letter in alphabet}
    for pos in range(len(bwt)):
        for letter in alphabet:
            is_this_letter = 1 if bwt[pos] == letter else 0
            occ_dict[(letter, pos)] = occ_dict[(letter, pos - 1)] + is_this_letter
    return {k: v for k, v in occ_dict.items() if k[1] >= 0}


def construct_lf_array(
    bwt: str, c_array: dict[str, int], occ_array: dict[tuple[str, int], int]
) -> dict[int, int]:
    d = {}
    for i in range(len(bwt)):
        d[i] = c_array[bwt[i]] + occ_array[(bwt[i], i)] - 1
    return d


def find_occurrences(subtxt: str, bwt: str, lf_array: dict[int, int]) -> [int]:
    """
    Find all occurrences of subtxt in txt.
    :param txt: input string
    :param subtxt: string to find
    :param bwt: Burrows-Wheeler transformata (txt)
    :param lf_array: Last-to-First array
    :return: list of positions of subtxt in txt
    """

    def traverse(idx: int, local_substr: str):
        if bwt[idx] == local_substr[-1]:
            if len(local_substr) == 1:
                return True
            else:
                return traverse(lf_array[idx], local_substr[:-1])
        else:
            return False

    matches = []
    for idx, letter in enumerate(bwt):
        if traverse(idx, subtxt):
            matches.append(idx)
    return matches


def test(txt: str) -> bool:
    good_res = compute_BWT_naive(txt)
    maybe_res = compute_BWT(txt)
    if good_res == maybe_res:
        print("Passed!")
        return True
    else:
        print("Wrong!")
        print(f"Good result: {good_res}")
        print(f"Maybe result: {maybe_res}")
        exit(1)


def find_substring(text: str, subtext: str, do_print: bool = True):
    """
    Find all occurrences of subtext in text.
    :param do_print: Print debug info
    :param text: Text to search in
    :param subtext: Subtext to search for
    :return: List of positions [begin, end) of subtext in text
    """

    def local_print(*args, **kwargs):
        if do_print:
            print(*args, **kwargs)

    local_print("=" * 80)
    local_print("text", text)
    bwt = compute_BWT(text)
    local_print("bwt:", bwt)
    c = construct_c_array(bwt)
    local_print("c:", c)
    occ = construct_occ_array(bwt)
    local_print("occ", occ)
    lf = construct_lf_array(bwt, c, occ)
    local_print("lf", lf)

    local_print("subtext:", subtext)
    local_print("=" * 80)

    matches = find_occurrences(subtext, bwt, lf)
    local_print("matches:", matches)
    rev_sa = construct_suffix_array(text + "~")
    mapped_matches = [rev_sa[x] for x in matches]
    local_print("mapped_matches:", mapped_matches)

    mark_array = ["_"] * len(text)

    for idx in mapped_matches:
        for i in range(len(subtext)):
            mark_array[idx - i - 1] = "^"

    local_print(text)
    local_print("".join(mark_array))

    result = [(x - len(subtext), x) for x in mapped_matches]
    for b, e in result:
        assert text[b:e] == subtext

    return result


if __name__ == "__main__":
    test("ababba")
    test("ACTG")
    test("ATGCTAGCATCTAGCTAC")
    test("AGTCGTACGTA")
    test("GCATGTCTAGCCAG")
    test("GCATGTCAGCCAG")

    text = "CGAGCCGCTTTCCATATCTATTAACGCATAAAAAACTCTGCTGGCATTCACAAATGCGCAGGGGTAAAACGTTTCCTGTAGCACCGTGAGTTATACTTTGT"
    subtext = "AAAC"
    substrs = find_substring(text, subtext)
    print("substrs:", substrs)
    #
    # big_file = 'data/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna_.chromosome.Chromosome.fa'
    # big_subfile = 'data/ecoli_proms.fa'
    #
    # big_text = list(parse(big_file, 'fasta'))[0].seq
    #
    # for big_subtext in parse(big_subfile, 'fasta'):
    #     print('big_subtext:', big_subtext)
    #     big_subtext = str(big_subtext.seq)
    #     substrs = find_substring(big_text, big_subtext, do_print=False)
    #     print('substrs:', substrs)
    #
