#include <bits/stdc++.h>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    int N;
    cin >> N;

    stringstream ss;

    vector<int> tab(N);

    for (auto &el : tab)
        cin >> el;

    int j = 0, maks = -1000000, cnt = 0;

    for (int i = 1; i <= N; ++i) {
        maks = max(maks, tab[i - 1]);

        if (i == maks) {
            ++cnt;
            ss << (i - j) << " ";
            for (int l = j + 1; l <= i; ++l) {
                ss << l << " ";
            }
            ss << "\n";

            j = i;
        }
    }

    cout << cnt << "\n";
    cout << ss.str();

    return 0;
}