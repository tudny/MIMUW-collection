#include <bits/stdc++.h>

using namespace std;

using type_t = uint64_t;
using score_t = uint64_t;

const score_t MOD = 1e9;
const size_t _MAX_N = 1010;

score_t dpp[_MAX_N][_MAX_N][2];

score_t add(score_t s, score_t val) {
    s += val;
    s %= MOD;
    return s;
}

// p, k, b - 0 jeżeli poprzedni dłuższy w lewo, 1 jeżeli w prawo
score_t dp(vector<type_t> &tab, size_t p, size_t k, bool koniec) {
    // cout << "p=" << (p + 1) << ", " << "k=" << (k + 1) << "\n";

    if (p > k) return 0;

    if (p == k) {
        score_t s = 0;

        if (koniec && k + 1 < tab.size()) {
            type_t b = tab[k + 1];
            
            if (tab[k] < b) {
                s = add(s, 1);
            }
        } else if (!koniec && p - 1 >= 0) {
            type_t a = tab[p - 1];

            if (tab[k] > a) {
                s = add(s, 1);
            }
        }

        dpp[p][k][koniec] = s;
    }

    if (dpp[p][k][koniec] != -1) {
        return dpp[p][k][koniec];
    }

    score_t s = 0;

    if (koniec && k + 1 < tab.size()) {
        type_t b = tab[k + 1];

        if (tab[k] < b) {
            s = add(s, dp(tab, p, k - 1, true));
        }
        if (tab[p] < b) {
            s = add(s, dp(tab, p + 1, k, false));
        }
    } else if (!koniec && p - 1 >= 0) {
        type_t a = tab[p - 1];

        if (tab[k] > a) {
            s = add(s, dp(tab, p, k - 1, true));
        }
        if (tab[p] > a) {
            s = add(s, dp(tab, p + 1, k, false));
        }
    }

    return dpp[p][k][koniec] = s;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    size_t N;
    cin >> N;

    vector<type_t> tab(N);

    for (auto &ele : tab)
        cin >> ele;

    for (size_t i = 0; i < _MAX_N; ++i) {
        for (size_t j = 0; j < _MAX_N; ++j) {
            for (size_t k = 0; k < 2; ++k) {
                dpp[i][j][k] = -1;
            }
        }
    }

    score_t s = 
        N == 1 ? 
            1 : add(dp(tab, 0, N - 2, true), dp(tab, 1, N - 1, false));

    cout << s << '\n';

    return 0;
}