#include <bits/stdc++.h>

using namespace std;

class SuperString {

    using hash_t = long long;
    using sub_str_t = pair<size_t, size_t>;

    const hash_t mod = 1e9 + 7;
    const hash_t p = 313;

    string s;
    vector<hash_t> hash;
    vector<hash_t> pows;

    hash_t sub_hash(sub_str_t a) {
        size_t i = a.first,
               j = a.second;

        return (((hash[j] - hash[i - 1] * pows[j - i + 1]) % mod) + mod) % mod;
    }

    size_t length(sub_str_t a) {
        return a.second - a.first + 1;
    }

public:
    SuperString(const string &_s) : s('$' + _s) {
        size_t len = s.length();

        hash.resize(len);
        pows.resize(len);

        hash[0] = 0;
        pows[0] = 1;

        for (size_t i = 1; i < len; ++i) {
            hash[i] = (hash[i - 1] * p + s[i]) % mod;
            pows[i] = (pows[i - 1] * p) % mod;
        }
    }

    bool are_subs_equal(sub_str_t a, sub_str_t b) {
        return sub_hash(a) == sub_hash(b);
    }

    char compare(sub_str_t a, sub_str_t b) {
        size_t lenA = length(a);
        size_t lenB = length(b);
        bool swapped = false;

        if (lenA < lenB) {
            swap(lenA, lenB);
            swap(a, b);
            swapped = true;
        }

        // a  b  a  a  b  a  b  a  a  b  a  a  b
        // 1  2  3  4  5  6  7  8  9  10 11 12 13
        //                   7  8  -  -  -  -  13

        size_t n = min(lenA, lenB);
        size_t x1 = a.first, y1 = a.second;
        size_t x2 = b.first, y2 = b.second;

        if (lenA == lenB && are_subs_equal({ x1, y1 }, { x2, y2 })) {
            return '=';
        }

        if (are_subs_equal({ x1, x1 + n - 1 }, { x2, x2 + n - 1 })) {
            return (swapped) ? '<' : '>';
        }

        size_t pocz = -1, kon = n;

        while (kon - pocz > 1) {
            size_t midd = (pocz + kon) / 2;

            if (are_subs_equal({ x1, x1 + midd }, { x2, x2 + midd })) {
                pocz = midd;
            } else {
                kon = midd;
            }
        }

        return (s[x1 + pocz + 1] < s[x2 + pocz + 1]) ? (!swapped ? '<' : '>') : (!swapped ? '>' : '<');
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    size_t a, b, c, d, N, Q;
    cin >> N >> Q;

    string s;
    cin >> s;

    SuperString ss = s;

    while (Q --> 0) {
        cin >> a >> b >> c >> d;
        cout << ss.compare({ a, b }, { c, d }) << "\n";
    }

    return 0;
}