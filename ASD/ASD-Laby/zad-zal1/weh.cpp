#include <bits/stdc++.h>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    // Variables
    vector<bool> zw;
    vector<int> e;
    vector<vector<int>> G;
    int p, k;
    size_t z, n, m;

    // Read
    cin >> p >> k >> z;
    zw.resize(p + 1, false);
    for (size_t val, i = 0; i < z; i++) {
        cin >> val;
        zw[val] = true;
    }

    cin >> n >> m;
    G.resize(n + 1);
    for (size_t a, b, i = 0; i < m; ++i) {
        cin >> a >> b;
        G[a].push_back(b);
        G[b].push_back(a);
    }

    e.resize(n);
    for (auto &ele : e) {
        cin >> ele;
    }

    // Solution for distance O(n + m)
    vector<int> d(n + 1, -1);
    queue<size_t> q;
    q.push(1);
    d[1] = 0;

    while (!q.empty()) {
        int v = q.front();
        q.pop();

        for (auto w : G[v]) {
            if (d[w] < 0) {
                d[w] = d[v] + 1;
                q.push(w);
            }
        }
    }

    // Solution for road with given length O(n + m)
    vector<size_t> road;
    size_t act = n;

    while ((int) road.size() != d[n] + 1) {
        road.push_back(act);

        for (auto w : G[act]) {
            if (d[w] < d[act]) {
                act = w;
                break;
            }
        }
    }

    reverse(road.begin(), road.end());

    // Solution for fill-ups O(np)
    // 0 - cannot, 1 - can w/o filling up, 2 - can only with a fill up
    int can[n][p + 1];
    for (size_t v = 0; v < n; ++v) 
        for (int j = 0; j <= p; ++j)
            can[v][j] = 0;

    can[0][p] = 1;
    for (size_t v = 1; v < n; ++v) {
        for (int j = k; j <= p; ++j) {
            can[v][j - k] = (can[v - 1][j] > 0);
        }

        for (int j = 0; j <= p; ++j) {
            if (can[v][j] == 1) {
                if (j + e[v] <= p && !zw[j + e[v]] && can[v][j + e[v]] != 1) {
                    can[v][j + e[v]] = 2;
                }
            }
        }
    }

    int best = -1;

    for (int j = 0; j <= p; ++j) {
        if (can[d[n]][j]) {
            best = j;
        }
    }

    // No solution
    if (d[n] < 0 || best < 0) {
        cout << -1 << "\n";
        return 0;
    }

    // Finding charge stations O(n)
    vector<int> charge;

    int actN = d[n];
    int actP = best;

    while (actN != 0) {
        if (can[actN][actP] == 1) {
            actP = actP + k;
        } else if (can[actN][actP] == 2) {
            charge.push_back(actN);
            actP = actP - e[actN] + k;
        } else if (can[actN][actP] == 0) {
            assert(false);
        }
        actN--;
    }

    reverse(charge.begin(), charge.end());

    // Print the solution
    cout << d[n] + 1 << " " << best << " " << charge.size() << "\n";

    for (auto r : road) {
        cout << r << " ";
    } cout << "\n";

    for (auto el : charge) {
        cout << road[el] << " ";
    } cout << "\n";

    return 0;
}
