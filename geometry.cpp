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

using P = complex<double>;
using G = vector<P>;
#define here(g, i) g[i]
#define next(g, i) g[(i + 1) % g.size()]
#define prev(g, i) g[(i - 1 + g.size()) % g.size()]
const double EPS = 1e-8;
const double INF = 1e12;
const double PI = acos(-1);
inline int sgn(double a, double b = 0) { return a < b - EPS ? -1 : a > b + EPS ? 1 : 0; }

struct L {
    P a, b, v, h;
    L(){}
    L(P _a, P _b) : a(_a), b(_b), v(b - a), h(v / abs(v) * P(0, 1)) {}
    L(double _ax, double _ay, double _bx, double _by) : L(P(_ax, _ay), P(_bx, _by)) {}
    void trans(double d) {
        a += d * h;
        b += d * h;
    }
    void print() {
        cerr << "{(" << a.real() << ", " << a.imag() << "), (" << b.real() << ", " << b.imag() << ")}" << endl;
    }
};
using S = L;

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

double cross(P a, P b) {
    return imag(conj(a) * b);
}

double dot(P a, P b) {
    return real(conj(a) * b);
}

int ccw(P p0, P p1, P p2) {
    p1 -= p0; p2 -= p0;
    if (sgn(cross(p1, p2)) == 1) return +1;         // counter-clockwise
    if (sgn(cross(p1, p2)) == -1) return -1;        // clockwise
    if (sgn(dot(p1, p2)) == -1) return +2;          // p2 -- p0 -- p1
    if (sgn(norm(p1), norm(p2)) == -1) return -2;   // p0 -- p1 -- p2
    return 0;                                       // p0 -- p2 -- p1
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

// Assume that there are two crosspoints
vector<P> crosspointLC(L l, C c) {
    double d = distanceLP(l, c.p);
    double len = sqrt(c.r * c.r - d * d);
    P h = l.h * d;
    if (!intersectPL(c.p + h, l)) h *= -1;
    assert(intersectPL(c.p + h, l));

    vector<P> ret;
    if (intersectPL(c.p + h + l.v / abs(l.v) * len, l)) ret.emplace_back(c.p + h + l.v / abs(l.v) * len);
    if (intersectPL(c.p + h - l.v / abs(l.v) * len, l)) ret.emplace_back(c.p + h - l.v / abs(l.v) * len);
    return ret;
}

P crosspointLL(L l1, L l2) {
    return l1.a + l1.v * cross(l2.v, l2.a - l1.a) / cross(l2.v, l1.v);
}

L intersectSegmentCC(C c1, C c2) {
    P v = c2.p - c1.p;
    double ac = (norm(v) + c1.r * c1.r - c2.r * c2.r) / (2 * abs(v));
    double as = sqrt(c1.r * c1.r - ac * ac);
    P u = v / abs(v);
    P h = u * P(0, as);
    u *= ac;
    return L(c1.p + u + h, c1.p + u - h);
}

// 2: parallel
// 1: orthogonal
// 0: otherwise
int relationLL(L l1, L l2) {
    if (cross(l1.v, l2.v) == 0) return 2;
    if (dot(l1.v, l2.v) == 0) return 1;
    return 0;
}

int relationCC(C c1, C c2) {
    double d = abs(c1.p - c2.p);
    if (sgn(d, c1.r + c2.r) == 1)       return 4; // distant
    if (sgn(d, c1.r + c2.r) == 0)       return 3; // touch outside
    if (sgn(d, abs(c1.r - c2.r)) == -1) return 0; // c1 in c2
    if (sgn(d, abs(c1.r - c2.r)) == 0)  return 1; // c1 touches in c2
    return 2; // two crosspoints
}

bool intersectLL(L l1, L l2) {
    return cross(l1.v, l2.v) != 0;
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

// exclude endpoints
bool intersectSC(L s, C c) {
    return sgn(distanceSP(s, c.p), c.r) == -1;
}

bool intersectPL(P p, L l) {
    return abs(ccw(l.a, l.b, p)) != 1;
}

bool intersectGC(G g, C c) {
    int n = g.size();
    rep(i, n) {
        L s(here(g, i), next(g, i));
        if (intersectSC(s, c)) return true;
    }
    return false;
}

double distanceLP(L l, P p) {
    return abs(cross(l.v, p - l.a)) / abs(l.v);
}

double distanceSP(L l, P p) {
    if (sgn(dot(l.v, p - l.a)) == -1) return abs(p - l.a);
    if (sgn(dot(-l.v, p - l.b)) == -1) return abs(p - l.b);
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

P projection(P p, P a, P b) {
    p -= a; b -= a;
    return dot(p, b) * b / abs(b) / abs(b) + a;
}

P reflection(P p, P a, P b) {
    P t = projection(p, a, b);
    return 2. * t - p;
}

L PerpendicularBisector(P p, P q) {
    P c = (p + q) / 2.;
    P d = (q - p) / 2.;
    P h = d * P(0, 1);
    return L(c + h, c - h);
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

G convex_cut(G g, L l) {
    G h;
    rep(i, (int)g.size()) {
        P p = here(g, i), q = next(g, i);
        if (ccw(p, q, l.a) == 0 && ccw(p, q, l.b) == 0) {
            if (ccw(p, l.b, l.a) == 0) return g;    // p -- l.a -- l.b -- q
            else return G{};                        // p -- l.b -- l.a -- q
        }
        if (ccw(l.a, l.b, p) != -1) h.emplace_back(p);
        if (ccw(l.a, l.b, p) * ccw(l.a, l.b, q) < 0)
            h.emplace_back(crosspointLL(L(p, q), l));
    }
    return h;
}

vector<G> Voronoi(G g, vector<P> p) {
    vector<G> ret;
    rep(i, p.size()) {
        G h = g;
        rep(j, p.size()) {
            if (i == j) continue;
            L l = PerpendicularBisector(p[i], p[j]);
            if (ccw(l.a, l.b, p[i]) == -1) l.reverse();
            h = convex_cut(h, l);
        }
        ret.emplace_back(h);
    }
    return ret;
}

bool cmp_x(const P& p, const P& q) {
    if (sgn(p.real(), q.real()) != 0) return sgn(p.real(), q.real()) == -1;
    return sgn(p.imag(), q.imag()) == -1;
}

G convex_hull(vector<P> ps) {
    int n = ps.size();
    sort(all(ps), cmp_x);
    G qs(n * 2);
    int k = 0;
    rep(i, n) {
        // 一直線上の3つ以上の頂点を残したい場合は「<=」を「<」に
        while (k > 1 && sgn(cross(qs[k - 1] - qs[k - 2], ps[i] - qs[k - 1])) < 0) k--;
        qs[k++] = ps[i];
    }
    for (int i = n - 2, t = k; i >= 0; i--) {
        // 一直線上の3つ以上の頂点を残したい場合は「<=」を「<」に
        while (k > t && sgn(cross(qs[k - 1] - qs[k - 2], ps[i] - qs[k - 1])) < 0) k--;
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

// signed!!!
double areaG(G g) {
    int n = g.size();
    double ret = 0;
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

void printG(G g) {
    rep(i, g.size()) {
        cout << g[i].real() << " " << g[i].imag() << endl;
    }
    cout << endl;
}

void printPoint(P p) {
    cout << fixed << setprecision(10) << p.real() << " " << p.imag() << endl;
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
