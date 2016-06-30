#include <bits/stdc++.h>
using namespace std;
#define int long long   // <-----!!!!!!!!!!!!!!!!!!!

#define rep(i,n) for (int i=0;i<(n);i++)
#define rep2(i,a,b) for (int i=(a);i<(b);i++)
#define rrep(i,n) for (int i=(n)-1;i>=0;i--)
#define rrep2(i,a,b) for (int i=(b)-1;i>=(a);i--)
#define all(a) (a).begin(),(a).end()

typedef long long ll;
typedef pair<int, int> Pii;
typedef tuple<int, int, int> TUPLE;
typedef vector<int> V;
typedef vector<V> VV;
typedef vector<VV> VVV;
const int inf = 1e9;
const int mod = 1e9 + 7;

struct edge {
    int to, cost;
    edge(){}
    edge(int _to, int _cost) : to(_to), cost(_cost) {}
};
typedef vector<vector<edge>> Graph;

int dijkstra(const Graph& G, int s, int g) {
    int N = G.size();
    priority_queue<Pii, vector<Pii>, greater<Pii>> pq;   // cost, vertex
    V d(N, inf);
    d[s] = 0;
    pq.push(make_pair(0, s));

    while (!pq.empty()) {
        auto p = pq.top(); pq.pop();
        int v = p.second;
        if (v == g) break;
        if (d[v] < p.first) continue;
        for (const auto& e : G[v]) {
            if (d[e.to] > d[v] + e.cost) {
                d[e.to] = d[v] + e.cost;
                pq.push(make_pair(d[e.to], e.to));
            }
        }
    }
    return d[g];
}
