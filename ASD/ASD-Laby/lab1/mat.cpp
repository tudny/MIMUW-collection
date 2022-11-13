#include <bits/stdc++.h>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    string str;
    cin >> str;

    string empty = "$";
    auto prev = empty.begin();

    long best = str.length() + 1;
    for (auto it = str.begin(); it != str.end(); it++) {
        if (*it != '*') {
            if (*prev != '$' && *prev != *it) {
                best = min(best, ((long) distance(prev, it)) + 1);
            }

            prev = it;
        }
    }

    cout << (str.length() - best + 2) << "\n";

    return 0;
}