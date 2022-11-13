#include <bits/stdc++.h>

#define L(x) (2 * x)
#define P(x) (2 * x + 1)

using namespace std;

using ll = uint32_t;

enum Update {
    SET_0 = 'Z',
    SET_1 = 'J',
    NOTHING = 'N'
};

class PartialTree {
    static const size_t log = 20;
    static const size_t p2 = (1 << log);
    vector<ll> drz;
    vector<Update> push;
    
public:
    PartialTree() {
        drz.resize(p2 << 1);
        fill(drz.begin(), drz.end(), 0);

        push.resize(p2 << 1);
        fill(push.begin(), push.end(), NOTHING);
    }

private:
    void update(size_t id, size_t aktP, size_t aktK) {
        if (push[id] == NOTHING)
            return;

        Update update = push[id];
        push[id] = NOTHING;

        push[L(id)] = update;
        push[P(id)] = update;
        size_t cnt = (aktK - aktP + 1) / 2;

        switch (update) {
            case SET_0:
                drz[L(id)] = 0;
                drz[P(id)] = 0;
                break;
            case SET_1:
                drz[L(id)] = cnt;
                drz[P(id)] = cnt;
                break;
            case NOTHING:
                break;
        }
    }

public:
    void set(size_t p, size_t k, Update val, size_t id = 1, size_t aktP = 1, size_t aktK = p2){

        if (p <= aktP && aktK <= k) {
            switch (val) {
                case SET_0:
                    drz[id] = 0;
                    push[id] = SET_0;
                    break;
                case SET_1:
                    drz[id] = (aktK - aktP + 1);
                    push[id] = SET_1;
                    break;
                case NOTHING:
                    assert(false);
            }

            return;
        }

        if (aktK < p || k < aktP) {
            return;
        }

        update(id, aktP, aktK);

        size_t midd = (aktP + aktK) / 2;
        set(p, k, val, L(id), aktP, midd);
        set(p, k, val, P(id), midd + 1, aktK);

        drz[id] = drz[L(id)] + drz[P(id)];
    }

    ll read(size_t p, size_t k, size_t id = 1, size_t aktP = 1, size_t aktK = p2) {
        if (p <= aktP && aktK <= k) {
            return drz[id];
        }
        if (aktK < p || k < aktP) {
            return 0;
        }

        update(id, aktP, aktK);

        size_t midd = (aktP + aktK) / 2;
        return read(p, k, L(id), aktP, midd) 
                + read(p, k, P(id), midd + 1, aktK);
    }

    void print() {
        for (size_t id = 1; id < p2 * 2; ++id) {
            if (!(id & (id - 1))) cout << "\n";
            cout << "[" << drz[id] << ", " << (char) push[id] << "] ";
        } cout << "\n\n";
    }
};


int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    size_t n, m;
    cin >> n >> m;

    PartialTree pt;

    // pt.print();

    while (m --> 0) {
        size_t a, b;
        char c;
        cin >> a >> b >> c;

        switch (c) {
            case 'C':
                pt.set(a, b, SET_0);
                break;
            case 'B':
                pt.set(a, b, SET_1);
                break;
        }

        ll read = pt.read(1, n);
        cout << read << "\n";

        // pt.print();
    }

    return 0;
}