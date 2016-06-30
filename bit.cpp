class Bit {
private:
    vector<int> bit;
    int n;

public:
    Bit(int _n) {
        n = _n;
        bit.clear();
        bit.resize(n + 1, 0);
    }

    // get sum in [1, i]
    // sum{[i, j]} = sum{[1, j]} - sum{[1, i-1]}
    int sum(int i) {
        int s = 0;
        while (i > 0) {
            s += bit[i];
            i -= i & -i;
        }
        return s;
    }

    int sum(int l, int r) {
        return sum(r) - sum(l - 1);
    }

    // add x to bit[i]
    void add(int i, int x) {
        while (i <= n) {
            bit[i] += x;
            i += i & -i;
        }
    }
};
