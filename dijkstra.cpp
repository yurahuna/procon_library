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
