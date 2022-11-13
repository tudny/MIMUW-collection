#include <bits/stdc++.h>

using namespace std;

using ll = unsigned long long;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    ll N, M;
    cin >> N >> M;

    vector<ll> masks(N, 0L);

    for (ll m = 0L; m < M; ++m) {
        for (ll a, i = 0L; i < N; ++i) {
            cin >> a;
            masks[a - 1] |= (((ll)(i < (N >> 1L))) << m);
        }
    }

    unordered_set<ll> s(masks.begin(), masks.end());

    cout << (s.size() == N ? "TAK":"NIE") << "\n"; 

    return 0;
}