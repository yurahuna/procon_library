const double EPS = 1e-10;

struct P {
    double x, y, z;
    P(){}
    P(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    const P operator+ (const P& p) {
        return P(this->x + p.x, this->y + p.y , this->z + p.z);
    }
    const P operator- (const P &p) {
        return P(this->x - p.x, this->y - p.y , this->z - p.z);
    }
    void print() {
        printf("%.1f %.1f %.1f\n", x, y, z);
    }
};

struct C {
    P c;
    double r;
    C(){}
    C(P _c, double _r) : c(_c), r(_r) {}
    C(double _x, double _y, double _z, double _r) : C(P(_x, _y, _z), _r) {}
};

struct L {
    P a, b;
    L(){}
    L(P _a, P _b) : a(_a), b(_b) {}
};

double dot(const P& p, const P& q) {
    return p.x * q.x + p.y * q.y + p.z * q.z;
}

P cross(P p, P q) {
    return P(p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x);
}

double distance(P p) {
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

double distancePP(P p, P q) {
    return sqrt((p.x - q.x) * (p.x - q.x)
            + (p.y - q.y) * (p.y - q.y)
            + (p.z - q.z) * (p.z - q.z));
}

double distanceLP(L l, P p) {
    return distance(cross(l.v, p - l.a)) / distance(l.v);
}

double distanceSP(L l, P p) {
    if (dot(l.b - l.a, p - l.a) < 0) return distancePP(p, l.a);
    if (dot(l.a - l.b, p - l.b) < 0) return distancePP(p, l.b);
    return distanceLP(l, p);
}

double distanceCC(C c1, C c2) {
    return max(0., distancePP(c1.c, c2.c) - c1.r - c2.r);
}

bool intersectSC(L s, C c) {
    // return distanceSP(s, c.c) <= c.r;
    return distanceSP(s, c.c) < c.r + EPS;
}

C readC() {
    double x, y, z, r;
    cin >> x >> y >> z >> r;
    return C(P(x, y, z), r);
}

L readL() {
    double x1, y1, z1, x2, y2, z2;
    cin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
    return L(P(x1, y1, z1), P(x2, y2, z2));
}
