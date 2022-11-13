#include <bits/stdc++.h>

using namespace std;

const int _MAX_N = 1 << 17;

vector<int> G[_MAX_N];
vector<int> T[_MAX_N];
vector<int> ile_do;
vector<int> P(_MAX_N);
vector<int> topo;
vector<bool> biore(_MAX_N);

int check(int k) {
    fill(biore.begin(), biore.end(), false);
    
    int sum = 0;
    for (auto v : topo) {
        biore[v] = P[v] <= k;
        for (auto w : G[v]) {
            biore[v] = biore[v] && biore[w];
        }
        if (biore[v]) {
            ++sum;
        }
    }

    return sum;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    int N, M, K;
    cin >> N >> M >> K;

    for (int i = 1; i <= N; ++i) {
        cin >> P[i];
    }

    ile_do.resize(N + 1);
    for (int a, b, i = 0; i < M; ++i) {
        cin >> a >> b;
        T[b].push_back(a);
        G[a].push_back(b);
        ile_do[a]++;
    }

    queue<int> q;
    for (int i = 1; i <= N; ++i) {
        if (ile_do[i] == 0) {
            q.push(i);
        }
    }

    while (!q.empty()) {
        int v = q.front(); q.pop();
        topo.push_back(v);

        for (auto w : T[v]) {
            if (--ile_do[w] == 0) {
                q.push(w);
            }
        }
    }

    int p = 0, k = 100000001;
    while (k - p > 1) {
        int midd = (p + k) / 2;

        if (check(midd) >= K) {
            k = midd;
        } else {
            p = midd;
        }
    }

    cout << k << "\n";

    return 0;
}