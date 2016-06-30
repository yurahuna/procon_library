class WarshallFloyd {
private:
    int n;
    VV d;
public:
    WarshallFloyd(int _n) : n(_n) {
        d.resize(n, V(n));
        rep(i, n) rep(j, n) {
            d[i][j] = (i == j ? 0 : inf);
        }
    }
    // directed
    void setDist(int i, int j, int c) {
        d[i][j] = c;
    }
    int getDist(int i, int j) {
        return d[i][j];
    }
    void calc() {
        rep(k, n) rep(i, n) rep(j, n) {
            d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
        }
    }
};
