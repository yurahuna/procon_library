typedef vector<vector<int>> Graph;
typedef pair<int, int>      Edge;   // (a < b: undirected)

class BICC {
private:
    const int n;

    Graph G;
    vi depth;
    vi par;
    map<Edge, int> imosEdge;
    map<Edge, int> EdgeType;
    enum {UNUSED, USED_DFS, BRIDGE};
    vector<Edge> bridges;

	vi cmp;
	int num_cc;
	vi size_of_vertex;
	Graph G_cc;
public:
    BICC(int _n) : n(_n), G(_n), depth(_n, -1), par(_n, -1), cmp(_n, -1), num_cc(0) {}
    Edge getEdge(int a, int b) {
        if (a > b) swap(a, b);
        return Edge(a, b);
    }
    void updateEdgeType(int a, int b, int type) {
        if (a < 0 || b < 0) return;
        EdgeType[getEdge(a, b)] = type;
    }
    void addEdge(int a, int b) {
        G[a].emplace_back(b);
        G[b].emplace_back(a);
        updateEdgeType(a, b, UNUSED);
    }
    void dfsTreeConstruct(int v, int pre) {
        if (depth[v] != -1) return;
        depth[v] = (pre == -1 ? 0 : depth[pre] + 1);
        par[v] = pre;
        updateEdgeType(pre, v, USED_DFS);
        for (auto&& nxt : G[v]) {
            if (nxt != pre) dfsTreeConstruct(nxt, v);
        }
    }
    void updateImos(int a, int b) {
        if (depth[a] < depth[b]) swap(a, b);

        if (par[a] != -1) {
            imosEdge[getEdge(a, par[a])]++;
        }
        if (par[b] != -1) {
            imosEdge[getEdge(b, par[b])]--;
        }
    }
    int imosFinal(int v, int pre) {
        int t = 0;
        for (auto&& nxt : G[v]) {
            if (nxt != pre && EdgeType[getEdge(nxt, v)] == USED_DFS) {
                t += imosFinal(nxt, v);
            }
        }
        if (pre != -1) imosEdge[getEdge(v, pre)] += t;
        return pre == -1 ? 0 : imosEdge[getEdge(v, pre)];
    }
    int extractCC(int v, int color) {
    	if (cmp[v] != -1) return 0;
    	cmp[v] = color;
    	int t = 1;
    	for (auto&& nxt : G[v]) {
    		if (EdgeType[getEdge(v, nxt)] != BRIDGE) {
    			t += extractCC(nxt, color);
    		}
    	}
    	return t;
    }
    void bicc() {
        dfsTreeConstruct(0, -1);
        for (auto&& p : EdgeType) {
            Edge e;
            int type;
            tie(e, type) = p;
            if (type == UNUSED) {
                updateImos(e.first, e.second);
            }
        }
        imosFinal(0, -1);
        for (auto&& p : EdgeType) {
            Edge e;
            int type;
            tie(e, type) = p;
            if (type == USED_DFS) {
                if (imosEdge[e] == 0) {
                    EdgeType[e] = BRIDGE;
                    bridges.emplace_back(e);
                }
            }
        }

		rep(i, n) {
			int size_cc = extractCC(i, num_cc);
			if (size_cc > 0) {
				size_of_vertex.emplace_back(size_cc);
				num_cc++;
			}
		}

	 	vector<set<int>> G_cc_st(num_cc);
		for (auto&& p : EdgeType) {
            Edge e;
            int type;
            tie(e, type) = p;
            if (type == BRIDGE) {
				G_cc_st[cmp[e.first]].insert(cmp[e.second]);
				G_cc_st[cmp[e.second]].insert(cmp[e.first]);
            }
        }

		rep(i, num_cc) {
			G_cc.emplace_back(vector<int>(all(G_cc_st[i])));
		}
    }
};
