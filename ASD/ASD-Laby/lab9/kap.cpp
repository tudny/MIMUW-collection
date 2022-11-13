#include <bits/stdc++.h>

using namespace std;

struct Point {
    int x, y;
    int id;

    Point() = default;
    Point(int x, int y) : x(x), y(y) {}

    int dis(const Point &other) const {
        int dx = abs(this->x - other.x);
        int dy = abs(this->y - other.y);

        return min(dx, dy);
    }
};

istream &operator>>(istream &is, Point &point) {
    return is >> point.x >> point.y;
}

ostream &operator<<(ostream &os, const Point &point) {
    return os << "[" << point.id << "]" << "(" << point.x << ", " << point.y << ")";
}

bool cmpX(const Point &a, const Point &b) {
    return a.x < b.x;
}

bool cmpY(const Point &a, const Point &b) {
    return a.y < b.y;
}

const long long INF = 1000000000000000000;

void DJ(int v, vector<long long> &dis, const vector<vector<pair<int, int>>> &G) {
    fill(dis.begin(), dis.end(), INF);
    dis[v] = 0;

    priority_queue<pair<long long, int>> q;
    q.push({ -dis[v], v });

    while (!q.empty()) {
        auto val = q.top(); q.pop();
        int v = val.second;
        long long d = -val.first;

        if (dis[v] != d)
            continue;

        for (const auto &vall : G[v]) {
            int w = vall.first;
            long long dd = vall.second;
            if (dis[w] > dis[v] + dd) {
                dis[w] = dis[v] + dd;
                q.push({ -dis[w], w });
            }
        }
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

    int N;
    cin >> N;

    vector<Point> points(N);
    vector<vector<pair<int, int>>> G(N, vector<pair<int, int>>(0));

    for (int i = 0; i < N; ++i) {
        cin >> points[i];
        points[i].id = i;
    }

    auto cmps = { cmpX, cmpY };

    for (const auto &cmp : cmps) {
        sort(points.begin(), points.end(), cmp);

        for (int i = 1; i < N; ++i) {
            G[points[i].id].push_back({ points[i - 1].id, points[i - 1].dis(points[i]) });
            G[points[i - 1].id].push_back({ points[i].id, points[i].dis(points[i - 1]) });
        }
    }

    vector<long long> dis(N);

    DJ(0, dis, G);

    cout << dis[N - 1] << "\n";

    return 0;
}