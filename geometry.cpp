#include <bits/stdc++.h>
using namespace std;
#define int long long   // <-----!!!!!!!!!!!!!!!!!!!

#define rep(i,n) for (int i=0;i<(n);i++)
#define rep2(i,a,b) for (int i=(a);i<(b);i++)
#define rrep(i,n) for (int i=(n)-1;i>=0;i--)
#define rrep2(i,a,b) for (int i=(a)-1;i>=b;i--)
#define all(a) (a).begin(),(a).end()
#define rall(a) (a).rbegin(),(a).rend()
#define printV(_v) for(auto _x:_v){cout<<_x<<" ";}cout<<endl
#define printVS(vs) for(auto x : vs){cout << x << endl;}
#define printVV(_vv) for(auto _v:_vv){for(auto _x:_v){cout<<_x<<" ";}cout<<endl;}
#define printP(p) cout << p.first << " " << p.second << endl
#define printVP(vp) for(auto p : vp) printP(p);

typedef long long ll;
typedef pair<int, int> Pii;
typedef tuple<int, int, int> TUPLE;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<vvi> vvvi;
typedef vector<Pii> vp;
const int inf = 1e9;
const int mod = 1e9 + 7;

typedef complex<double> P;
typedef vector<P> G;
#define here(g, i) g[i]
#define next(g, i) g[(i + 1) % g.size()]
#define prev(g, i) g[(i - 1 + g.size()) % g.size()]
const double EPS = 1e-10;
const double INF = 1e12;
const double PI = acos(-1);

struct L {
    P a, b, v;
    L(){}
    L(P _a, P _b) : a(_a), b(_b), v(b - a) {}
    L(double _ax, double _ay, double _bx, double _by) : L(P(_ax, _ay), P(_bx, _by)) {}
};

struct C {
    P p;
    double r;
    C(){}
    C(P _p, double _r) : p(_p), r(_r) {}
    C(double _x, double _y, double _r) : p(_x, _y), r(_r) {}
    void print() {
        cerr << p.real() << " " << p.imag() << " " << r << endl;
    }
};

bool cmp_x(const P& p, const P& q) {
    if (p.real() != q.real()) return p.real() < q.real();
    return p.imag() < q.imag();
}

double cross(P a, P b) {
    return imag(conj(a) * b);
}

double dot(P a, P b) {
    return real(conj(a) * b);
}

int ccw(P p0, P p1, P p2) {
    if (cross(p1 - p0, p2 - p0) > 0) return +1; // counter-clockwise
    if (cross(p1 - p0, p2 - p0) < 0) return -1; // clockwise
    if (dot(p1 - p0, p2 - p0) < 0) return +2;   // online_back
    if (dot(p0 - p1, p2 - p1) < 0) return -2;   // online_front
    return 0;                                   // on_segment
}

// !! 端点を含まないときは "<" !!
bool intersectSS(L l1, L l2) {
    return (ccw(l1.a, l1.b, l2.a) * ccw(l1.a, l1.b, l2.b) <= 0 &&
            ccw(l2.a, l2.b, l1.a) * ccw(l2.a, l2.b, l1.b) <= 0);
}


P crosspointSS(L l1, L l2) {
    assert(intersectSS(l1, l2));
    double d1 = abs(cross(l2.v, l1.a - l2.a));
    double d2 = abs(cross(l2.v, l1.b - l2.a));
    double t = d1 / (d1 + d2);
    return l1.a + t * l1.v;
}

// 2: parallel
// 1: orthogonal
// 0: otherwise
int relationLL(L l1, L l2) {
    if (cross(l1.v, l2.v) == 0) return 2;
    if (dot(l1.v, l2.v) == 0) return 1;
    return 0;
}


bool intersectLL(L l1, L l2) {
    return cross(l1.v, l2.v) != 0;
}

P crosspointLL(L l1, L l2) {
    return l1.a + l1.v * cross(l2.v, l2.a - l1.a) / cross(l2.v, l1.v);
}

bool intersectSG(L l, G g) {
    int n = g.size();
    rep(i, n) {
        if (intersectSS(l, L(here(g, i), next(g, i)))) {
            return true;
        }
    }
    return false;
}


double distanceLP(L l, P p) {
    return abs(cross(l.v, p - l.a)) / abs(l.v);
}

double distanceSP(L l, P p) {
    if (dot(l.v, p - l.a) < 0) return abs(p - l.a);
    if (dot(-l.v, p - l.b) < 0) return abs(p - l.b);
    return distanceLP(l, p);
}

double distanceSS(L l1, L l2) {
    if (intersectSS(l1, l2)) return 0;
    double d = INF;
    d = min(d, distanceSP(l1, l2.a));
    d = min(d, distanceSP(l1, l2.b));
    d = min(d, distanceSP(l2, l1.a));
    d = min(d, distanceSP(l2, l1.b));
    return d;
}

double distanceCC(C c1, C c2) {
    double d = abs(c1.p - c2.p);
    return max(0., d - (c1.r + c2.r));
}

P projection(P x, P d) {
    return dot(x, d) * d / abs(d) / abs(d);
}

P reflection(P x, P d) {
    P t = projection(x, d);
    return 2. * t - x;
}

double deg2rad(double deg) {
    return deg * 2. * PI / 360;
}

// rot p around q by theta (counter-clockwise)
P rotP(P p, P q, double theta) {
    p -= q;
    double x = p.real(), y = p.imag();
    p = P(x * cos(theta) - y * sin(theta), x * sin(theta) + y * cos(theta));
    p += q;
    return p;
}

G convex_hull(vector<P> ps) {
    int n = ps.size();
    sort(all(ps), cmp_x);
    G qs(n * 2);
    int k = 0;
    rep(i, n) {
        // 一直線上の3つ以上の頂点を残したい場合は「<=」を「<」に
        while (k > 1 && cross(qs[k - 1] - qs[k - 2], ps[i] - qs[k - 1]) <= 0) k--;
        qs[k++] = ps[i];
    }
    for (int i = n - 2, t = k; i >= 0; i--) {
        // 一直線上の3つ以上の頂点を残したい場合は「<=」を「<」に
        while (k > t && cross(qs[k - 1] - qs[k - 2], ps[i] - qs[k - 1]) <= 0) k--;
        qs[k++] = ps[i];
    }
    qs.resize(k - 1);
    return qs;
}

// ベクトルp, qのなす角(q->pの方向が正)
double thetaPP(P p, P q) {
    int sgn = cross(q, p) > 0 ? +1 : -1;
    return sgn * acos(dot(p, q) / abs(p) / abs(q));
}

// lがx軸となす角
double thetaL(L l) {
    return thetaPP(l.v, P(1, 0));
}

G rotG(G g, P p, double theta) {
    rep(i, g.size()) {
        g[i] = rotP(g[i], p, theta);
    }
    return g;
}

P centroidG(G g) {
    int n = g.size();
    double x = 0, y = 0;
    rep(i, n) {
        x += g[i].real();
        y += g[i].imag();
    }
    return P(x / n, y / n);
}

double areaG(G g) {
    int n = g.size();
    double ret;
    rep(i, n) {
        ret += cross(here(g, i), next(g, i));
    }
    return ret / 2;
}

bool isConvex(G g) {
    int n = g.size();
    rep(i, n) {
        // when 180 deg is permmitted
        if (ccw(prev(g, i), here(g, i), next(g, i)) == -1) {
            return false;
        }
    }
    return true;
}

enum { OUT, ON, IN };
int containPG(const P& p, const G& g) {
    int n = g.size();
    bool in = false;
    rep(i, n) {
        P a = here(g, i) - p, b = next(g, i) - p;
        if (a.imag() > b.imag()) swap(a, b);
        if (a.imag() <= 0 && 0 < b.imag() && cross(a, b) > 0) in = !in;
        if (cross(a, b) == 0 && dot(a, b) <= 0) return ON;
    }
    return in ? IN : OUT;
}

P readP() {
    double x, y;
    cin >> x >> y;
    return P(x, y);
}

L readL() {
    double xa, ya, xb, yb;
    cin >> xa >> ya >> xb >> yb;
    return L(xa, ya, xb, yb);
}

G readG() {
    int n;
    cin >> n;
    G g(n);
    rep(i, n) g[i] = readP();
    return g;
}

C readC() {
    double x, y, r;
    cin >> x >> y >> r;
    return C(x, y, r);
}

bool EQ(double a, double b) {
    return abs(a - b) < EPS;
}

bool EqP(P p, P q) {
    return EQ(p.real(), q.real()) && EQ(p.imag(), q.imag());
}

bool EqG(G g, G h) {
    if (g.size() != h.size()) return false;
    rep(k, g.size()) {
        bool flag = true;
        rep(i, g.size()) {
            if (!EqP(g[(i + k) % g.size()], h[i])) {
                flag = false;
                break;
            }
        }
        if (flag) {
            return true;
        }
    }
    return false;
}



void printG(G g) {
    rep(i, g.size()) {
        cout << g[i].real() << " " << g[i].imag() << endl;
    }
    cout << endl;
}


signed main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);

    double x0, y0, x1, y1;
    cin >> x0 >> y0 >> x1 >> y1;
    P p0(x0, y0), p1(x1, y1);
    int q;
    cin >> q;
    while (q--) {
        double x2, y2;
        cin >> x2 >> y2;
        P p2(x2, y2);
        P a = project(p2 - p0, p1 - p0) + p0;
        printf("%.10f %.10f\n", a.real(), a.imag());
    }

}
