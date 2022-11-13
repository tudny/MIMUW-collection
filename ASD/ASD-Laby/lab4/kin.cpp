#include <bits/stdc++.h>

using namespace std;

using ll = long long;
const ll MOD = 1e9;

class PartialTree {
    static const size_t log = 15;
    static const size_t p2 = (1 << log);
    vector<ll> drz;
    
public:
    PartialTree() {
        drz.resize(p2 << 1);
    }

    void add(size_t nr, ll val) {
        nr = nr + p2 - 1;
        drz[nr] += val;
        drz[nr] %= MOD;
        nr /= 2;

        while (nr > 0) {
            drz[nr] = (drz[2 * nr] + drz[2 * nr + 1]) % MOD;
            nr /= 2;
        }
    }

    ll read(size_t p, size_t k, size_t id = 1, size_t aktP = 1, size_t aktK = p2) {
        if (p <= aktP && aktK <= k) {
            return drz[id];
        }
        if (aktK < p || k < aktP) {
            return 0;
        }

        size_t midd = (aktP + aktK) / 2;
        return (read(p, k, 2 * id, aktP, midd) 
                + read(p, k, 2 * id + 1, midd + 1, aktK)) % MOD;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    size_t N, K;
    cin >> N >> K;

    vector<ll> tab(N + 1);
    for (size_t i = 1; i <= N; i++) {
        cin >> tab[i];
    }

    vector<vector<ll>> dp(K + 1, vector<ll>(N + 1, 0));

    for (size_t n = 1; n <= N; ++n) {
        dp[1][n] = 1;
    }

    for (size_t k = 2; k <= K; ++k) {
        PartialTree pt;
        for (size_t n = 1; n <= N; ++n) {
            dp[k][n] = pt.read(tab[n], N);
            pt.add(tab[n], dp[k - 1][n]);
        }
    }

    ll sum = 0L;
    for (size_t n = 1; n <= N; ++n) {
        sum = (sum + dp[K][n]) % MOD;
    }

    cout << sum << "\n";

    return 0;
}