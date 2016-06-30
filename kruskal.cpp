class UnionFind {
private:
    const int n;
    vector<int> uni;
public:
    UnionFind(int _n) : n(_n), uni(_n, -1) {}
    int root(int x) {
        if (uni[x] < 0) return x;
        return uni[x] = root(uni[x]);
    }
    bool same(int x, int y) {
        return root(x) == root(y);
    }
    bool unite(int x, int y) {
        x = root(x);
        y = root(y);
        if (x == y) return false;
        if (uni[x] > uni[y]) swap(x, y);
        uni[x] += uni[y];
        uni[y] = x;
        return true;
    }
    void print() {
        for (auto x : uni) cout << x << " ";
        cout << endl;
    }
};

template <typename T>
struct edge {
    int from, to;
    T cost;
    edge(){}
    edge(int _from, int _to, T _cost) : from(_from), to(_to), cost(_cost) {}
    bool operator< (const edge& e) const {
        return cost == e.cost ? (from == e.from ? to < e.to : from < e.from) : cost < e.cost;
    }
};

template <typename T>
class Kruskal {
private:
    int n;
    vector<edge<T>> edges;
    UnionFind uf;
public:
    Kruskal(int _n) : n(_n), uf(_n) {}
    void addEdge(int _from, int _to, T _cost) {
        edges.emplace_back(_from, _to, _cost);
    }
    T calc() {
        sort(all(edges));
        T res = 0;
        for (auto &e : edges) {
            if (uf.unite(e.from, e.to)) {
                res += e.cost;
            }
        }
        return res;
    }
};
/* !!!!! T: 辺のコストの型(int, double) !!!!! */
