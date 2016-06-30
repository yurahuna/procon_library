class SegTreeMax {
            private:
                const int INF = (1LL << 31) - 1;
                int n;  // width of bottom
                vector<int> dat;

            public:
                SegTreeMax(int _n) {
                    n = 1;
                    while (n < _n) n *= 2;
                    dat.clear();
                    dat.resize(2 * n - 1, -INF);
                }

                void update(int k, int a) {
                    k += n - 1;
                    dat[k] = a;
                    while (k > 0) {
                        k = (k - 1) / 2;
                        dat[k] = max(dat[2 * k + 1], dat[2 * k + 2]);
                    }
                }

                // get max value in [a, b)
                // dat[k] corresponds to [l, r)
                int query(int a, int b, int k, int l, int r) {
                    if (r <= a || b <= l) return -INF;
                    if (a <= l && r <= b) return dat[k];
                    return max(query(a, b, k * 2 + 1, l, (l + r) / 2), query(a, b, k * 2 + 2, (l + r) / 2, r));
                }

                int query(int a, int b) { return query(a, b, 0, 0, n); }
            };
