#include <bits/stdc++.h>
using namespace std;
#define int long long   // <-----!!!!!!!!!!!!!!!!!!!

#define rep(i,n) for (int i=0;i<(n);i++)
#define rep2(i,a,b) for (int i=(a);i<(b);i++)
#define rrep(i,n) for (int i=(n)-1;i>=0;i--)
#define rrep2(i,a,b) for (int i=(b)-1;i>=(a);i--)
#define all(a) (a).begin(),(a).end()

typedef long long ll;
typedef pair<int, int> Pii;
typedef tuple<int, int, int> TUPLE;
typedef vector<int> V;
typedef vector<V> VV;
typedef vector<VV> VVV;
const int inf = 1e9;
const int mod = 1e9 + 7;

struct Point3d {
    int x, y, z;
    Point3d(){}
    Point3d(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {}
    Point3d cross(Point3d p, Point3d q) {
        return Point3d(p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x);
    }
    // axis: 0:x, 1:y, 2:z
    // s: counter-clockwise-> +1, else-> -1
    Point3d rot(int axis, int s) {
        int c = 0;
        VV a;
        switch(axis) {
            case 0: // x
                a = {{1, 0, 0},
                     {0, c, -s},
                     {0, s, c}};
                break;
            case 1:  // y
                a = {{c, 0, s},
                     {0, 1, 0},
                     {-s, 0, c}};
                break;
            case 2:  // z
                a = {{c, -s, 0},
                     {s, c, 0},
                     {0, 0, 1}};
                break;
        }
        Point3d b;
        b.x = a[0][0] * x + a[0][1] * y + a[0][2] * z;
        b.y = a[1][0] * x + a[1][1] * y + a[1][2] * z;
        b.z = a[2][0] * x + a[2][1] * y + a[2][2] * z;
        return b;
    }
    void print() {
        cout << x << " " << y << " " << z << endl;
    }
    bool operator<(const Point3d &p) const {
       return x == p.x ? (y == p.y ? z < p.z : y < p.y) : x < p.x;
    }
    // bool operator>(const Point3d &p) const {
    //   return x == p.x ? (y == p.y ? z > p.z : y > p.y) : x > p.x;
    // }
};

struct Dice {
    map<Point3d, int> dice;
    Dice() {}
    Dice(V v) {
        dice[Point3d(1, 0, 0)] = v[0];
        dice[Point3d(0, 1, 0)] = v[1];
        dice[Point3d(0, 0, 1)] = v[2];
        dice[Point3d(-1, 0, 0)] = v[5];
        dice[Point3d(0, -1, 0)] = v[4];
        dice[Point3d(0, 0, -1)] = v[3];
    }
    void rot(int axis, int s) {
        map<Point3d, int> dice_;
        dice_[Point3d(1, 0, 0).rot(axis, s)] = dice[Point3d(1, 0, 0)];
        for (/*const auto&*/ auto p : dice) {
            p.first.rot(axis, s);
            // dice_[b] = p.second;
            // dice_[p.first.rot(axis, s)] = p.second;
        }
        dice = dice_;
    }
    int get(int x, int y, int z) {
        return dice[Point3d(x, y, z)];
    }
};

signed main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(0);

    // Point3d p(1, 0, 0);
    // p.rot(2, 1);
    // p.print();
    // p.rot(0, 1);
    // p.print();
    // p.rot(1, 1);
    // p.print();

    // Dice d(V{1, 2, 3, 4, 5, 6});
    // cout << d.get(1, 0, 0) << endl;

}
