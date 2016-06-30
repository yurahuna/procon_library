template <typename T>
struct edge {
    int from, to;
    T cost;
    edge(){}
    edge(int _from, int _to, T _cost) : from(_from), to(_to), cost(_cost) {}
    bool operator< (const edge& e) const {
        return cost == e.cost ? (from == e.from ? to < e.to : from < e.from) : cost < e.cost;
    }
    bool operator> (const edge& e) const {
        return cost == e.cost ? (from == e.from ? to > e.to : from > e.from) : cost > e.cost;
    }
};

template <typename T>
class Prim {
private:
    int n;
    vector<vector<edge<T>>> G;
public:
    Prim(int _n) : n(_n) {
        G.resize(n);
    }
    // undirected
    void addEdge(int u, int v, T c) {
        G[u].emplace_back(u, v, c);
        G[v].emplace_back(v, u, c);
    }
    T calc(int s = 0) {
        T total_weight = 0;
        vector<bool> visited(n);
        int num_visited = 0;
        priority_queue<edge<T>, vector<edge<T>>, greater<edge<T>>> pq;
        pq.push(edge<T>(-1, s, 0));
        while (!pq.empty() && num_visited < n) {
            auto e = pq.top(); pq.pop();
            if (visited[e.to]) continue;
            total_weight += e.cost;
            visited[e.to] = true;
            num_visited++;
            for (const auto& ne: G[e.to]) {
                if (!visited[ne.to]) pq.push(ne);
            }
        }
        return total_weight;
    }
};
/* !!!!! T: 辺のコストの型(int, double) !!!!! */
