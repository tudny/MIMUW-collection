#include <bits/stdc++.h>

using namespace std;

const int LOG = 19;
const int _MAX_N = (1 << LOG);
vector<int> G[_MAX_N];

void dfs(int v, vector<int> &d) {
    for (auto w : G[v]) {
        if (d[w] == -1) {
            d[w] = d[v] + 1;
            dfs(w, d);
        }
    }
}

pair<int, int> diameter() {
    vector<int> d(_MAX_N, -1);

    d[1] = 0;
    dfs(1, d);

    int max_v = 1;
    for (int i = 1; i < _MAX_N; ++i) {
        if (d[i] > d[max_v]) {
            max_v = i;
        }
    }

    fill(d.begin(), d.end(), -1);
    d[max_v] = 0;
    dfs(max_v, d);
    
    int max_w = 1;
    for (int i = 1; i < _MAX_N; ++i) {
        if (d[i] > d[max_w]) {
            max_w = i;
        }
    }

    return { max_v, max_w };
}

void calcluate_jump(int v, int w, vector<int> &d, vector<vector<int>> &jump) {
    jump[v][0] = w;
    for (int l = 1; l < LOG; ++l) {
        jump[v][l] = jump[jump[v][l - 1]][l - 1];
    }

    for (auto z : G[v]) {
        if (d[z] == -1) {
            d[z] = d[v] + 1;
            calcluate_jump(z, v, d, jump);
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    int N;
    cin >> N;

    for (int i = 1; i <= N; i++) {
        for (int a, k = 0; k < 2; ++k) {
            cin >> a;
            if (a > 0) {
                G[i].push_back(a);
                G[a].push_back(i);
            }
        }
    }

    auto diam = diameter();

    vector<vector<int>> jump1(_MAX_N, vector<int>(LOG, 0));
    vector<int> vis1(_MAX_N, -1);

    vis1[diam.first] = 0;
    calcluate_jump(diam.first, 0, vis1, jump1);

    vector<vector<int>> jump2(_MAX_N, vector<int>(LOG, 0));
    vector<int> vis2(_MAX_N, -1);

    vis2[diam.second] = 0;
    calcluate_jump(diam.second, 0, vis2, jump2);

    int M;
    cin >> M;

    while (M --> 0) {
        int v, d;
        cin >> v >> d;

        int w = v;
        for (int l = 0; l < LOG; ++l) {
            if (d & (1 << l)) {
                w = jump1[w][l];
            }
        }

        if (w == 0) {
            w = v;
            for (int l = 0; l < LOG; ++l) {
                if (d & (1 << l)) {
                    w = jump2[w][l];
                }
            }
        }

        if (w != 0) {
            cout << w << "\n";
        } else {
            cout << -1 << "\n";
        }
    }

    return 0;
}