#include <iostream>
#include <vector>

using namespace std;

using type_t = long long;

int main() {

    ios_base::sync_with_stdio(false);
    cout.tie(0);
    cin.tie(0);

    size_t N;
    cin >> N;

    vector<type_t> tab(N + 1);
    vector<type_t> prefSum(N + 1);

    prefSum[0] = 0;

    for (size_t i = 1; i <= N; ++i) {
        cin >> tab[i];
        prefSum[i] = prefSum[i - 1] + tab[i];
    }

    const type_t INF = 1000000000000009;

    vector<type_t> maxLeftOdd(N + 1), maxLeftEven(N + 1),
                    minRightOdd(N + 2), minRightEven(N + 2);
    
    maxLeftEven[0] = -INF;
    maxLeftOdd[0] = -INF;

    minRightEven[N + 1] = INF;
    minRightOdd[N + 1] = INF;

    for (size_t i = 1; i <= N; i++) {
        if (tab[i] & 1) {
            maxLeftOdd[i] = max(maxLeftOdd[i - 1], tab[i]);
            maxLeftEven[i] = maxLeftEven[i - 1];
        } else {
            maxLeftOdd[i] = maxLeftOdd[i - 1];
            maxLeftEven[i] = max(maxLeftEven[i - 1], tab[i]);
        }
    }

    for (size_t i = N; i > 0; i--) {
        if (tab[i] & 1) {
            minRightOdd[i] = min(minRightOdd[i + 1], tab[i]);
            minRightEven[i] = minRightEven[i + 1];
        } else {
            minRightOdd[i] = minRightOdd[i + 1];
            minRightEven[i] = min(minRightEven[i + 1], tab[i]);
        }
    }

    int K;
    cin >> K;

    while (K --> 0) {
        int k;
        cin >> k;

        int m = N - k;

        type_t sum = prefSum[N] - prefSum[m];

        if (sum & 1) {
            cout << sum << "\n";
            continue;
        }

        type_t l1 = maxLeftEven[m],
                r1 = minRightOdd[m + 1],
                l2 = maxLeftOdd[m],
                r2 = minRightEven[m + 1];

        type_t s1 = sum - r1 + l1;
        type_t s2 = sum - r2 + l2;

        bool s1good = abs(l1) != INF && abs(r1) != INF;
        bool s2good = abs(l2) != INF && abs(r2) != INF;

        type_t maks = 0;
        if (s1good && s2good) {
            maks = max(s1, s2);
        } else if (s1good) {
            maks = s1;
        } else if (s2good) {
            maks = s2;
        } else {
            maks = -1;
        }

        cout << maks << "\n";
    }

    return 0;
}
