// Written by Qiwei Feng, 2015.05.05
#include <algorithm>
#include <iostream>
#include <cassert>
#include <random>
#include <vector>
#include <set>
#include <map>

//-----------------------------------------------
//-- Global Variables
//-----------------------------------------------
// {{{

template<class T> T sqr(T x) { return x * x; }

struct Point {
    double x, y;
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}
    Point operator +(const Point p) const {
        return Point(x + p.x, y + p.y);
    }
    Point operator -(const Point p) const {
        return Point(x - p.x, y - p.y);
    }
    friend double sqr_dist(Point p1, Point p2) {
        return sqr(p1.x - p2.x) + sqr(p1.y - p2.y);
    }
    friend double dist(Point p1, Point p2) {
        return sqrt(sqr_dist(p1, p2));
    }
    friend double mult(Point p1, Point p2) {
        return p1.x * p2.y - p1.y * p2.x;
    }
    friend double mult(Point p0, Point p1, Point p2) {
        return mult(p1 - p0, p2 - p0);
    }
    bool operator <(const Point &p) const {
        return y < p.y;
    }
    void normalize() {
        double len = dist(*this, Point());
        x /= len, y /= len;
    }
};

const int MAXN = 100000;

Point sites[MAXN];
int n;
double sweepline;

void read_sites() {
    n = 0;
    while (true) {
        int i = n + 1;
        if (scanf("%lf%lf", &sites[i].x, &sites[i].y) < 2)
            break;
        n = i;
    }
}

template<class T> struct Maybe {
    bool ok;
    T val;
    Maybe() : ok(false) {}
    Maybe(T v) : ok(true), val(v) {}
    static Maybe nothing() {
        return Maybe();
    }
    T get() const {
        if (!ok)
            throw "error";
        return val;
    }
};

// }}}

//-----------------------------------------------
//-- Basic Calculations
//-----------------------------------------------
// {{{

Maybe<Point> raw_intersect(int i, int j, int k) {
    // make sure i, j, k be valid and different
    assert(i != j && j!= k && j != 0);
    if (i == 0 || k == 0 || i == k)
        return Maybe<Point>::nothing();

    if (mult(sites[i], sites[j], sites[k]) < 0)
        return Maybe<Point>::nothing();

    Point p0 = sites[i], p1 = sites[j], p2 = sites[k];
    double a1 = p1.x - p0.x, b1 = p1.y - p0.y, c1 = (sqr(a1) + sqr(b1)) / 2.;
    double a2 = p2.x - p0.x, b2 = p2.y - p0.y, c2 = (sqr(a2) + sqr(b2)) / 2.;
    double d = a1 * b2 - a2 * b1;
    double x = p0.x + (c1 * b2 - c2 * b1) / d;
    double y = p0.y + (a1 * c2 - a2 * c1) / d;
    return Maybe<Point>(Point(x, y));
}

Maybe<Point> intersect(int i, int j, int k) {
    Maybe<Point> p = raw_intersect(i, j, k);
    if (p.ok) {
        double d = dist(sites[i], p.val);
        p.val.y += d;
    }
    return p;
}

Maybe<std::pair<double, double>> solve_equation(double a, double b, double c) {
    double delta = b * b - 4. * a * c;
    if (delta < -1e-10)
        return Maybe<std::pair<double, double>>::nothing();
    double sqrt_delta = delta < 0 ? 0. : sqrt(delta);
    double x1 = (-b + sqrt_delta) / (2. * a);
    double x2 = (-b - sqrt_delta) / (2. * a);
    if (x1 > x2)
        std::swap(x1, x2);
    return Maybe<std::pair<double, double>>(std::make_pair(x1, x2));
}

// }}}

//-----------------------------------------------
//-- Tree for List
//-----------------------------------------------
// {{{

class Bisector {
    double xpos(double y) const {
        assert(i != 0 && j != 0);
        double mx = (sites[i].x + sites[j].x) / 2.;
        double my = (sites[i].y + sites[j].y) / 2.;
        if (sites[i].y == sites[j].y)
            return mx;
        double k = (sites[j].x - sites[i].x) / (sites[i].y - sites[j].y);
        double b = my - k * mx;
        double x0 = sites[i].x, y0 = sites[i].y;
        double t = y - b, t0 = y0 - b;
        auto sol = solve_equation(1,
            2 * k * t - 2 * k * t0 - 2 * x0,
            x0 * x0 + t0 * t0 - t * t);
        assert(sol.ok);
        double ans = 0;
        if (sites[i].y < sites[j].y)
            ans = sol.val.first;
        else
            ans = sol.val.second;
        //std::cerr.precision(10);
        //std::cerr << "xpos [" << i << " " << j << "] at " << y << " -> " << ans << std::endl;
        return ans;
    }
    double yval(double x) const {
        double mx = (sites[i].x + sites[j].x) / 2.;
        double my = (sites[i].y + sites[j].y) / 2.;
        double k = (sites[j].x - sites[i].x) / (sites[i].y - sites[j].y);
        double b = my - k * mx;
        double y = k * x + b;
        y += dist(Point(x, y), sites[i]);
        return y;
    }
public:
    int i, j;
    Bisector(int i, int j) : i(i), j(j) {}
    Point bottom() const {
        if (sites[i].y == sites[j].y)
            return Point((sites[i].x + sites[j].x) / 2., (sites[i].y + sites[j].y) / 2.);
        else {
            if (sites[i].y > sites[j].y)
                return sites[i];
            else
                return sites[j];
        }
    }
    int bottom_i() const {
        if (sites[i].y == sites[j].y)
            return std::min(i, j);
        else
            return (sites[i].y > sites[j].y) ? i : j;
    }
    bool inside(Point p) const {
        //std::cerr.precision(3);
        //std::cerr << "inside [" << i << " " << j << "] (" << p.x << ", " << p.y << " ";
        //std::cerr << " yval=" << yval(p.x) << std::endl;
        return p.y >= yval(p.x);
    }
    bool is_left() const {
        if (sites[i].y < sites[j].y)
            return true;
        return false;
    }
    bool operator <(const Bisector &b) const {
        if (i == b.i && j == b.j)
            return false;
        if (i == 0 || b.i == 0) // one is the leftmost one
            return i == 0;
        if (j == 0 || b.j == 0) // one is the rightmost one
            return b.j == 0;
        if (i == b.j && j == b.i) // of same raw_bisector
            return sites[i].y < sites[j].y;
        //std::cerr << "comparing [" << i << " " << j << "] [" << b.i << " " << b.j << "]" << std::endl;

        if (bottom_i() == b.bottom_i()) {
            if (is_left() && !b.is_left())
                return true;
            if (!is_left() && b.is_left())
                return false;
            int k1 = i + j - bottom_i();
            int k2 = b.i + b.j - b.bottom_i();
            //std::cerr << " -- same bottom" << std::endl;
            return sites[k1].x < sites[k2].x;
        }
        
        if (inside(b.bottom())) {
            if (sites[i].y < sites[j].y)
                return true;
            return false;
        }
        if (b.inside(bottom())) {
            return !(b < *this);
        }
        return bottom().x < b.bottom().x;
    }
    bool operator ==(const Bisector &b) const {
        if (i == b.i && j == b.j)
            return true;
        return false;
    }
    bool after(double x) const {
        if (i == 0) // leftmost
            return false;
        if (j == 0) // rightmost
            return true;
        return x < xpos(sweepline);
    }
};

// Balanced Binary Search Tree
// Implementation: Treap

struct Node {
    Bisector val;
    unsigned randkey;
    Node *prev, *next;
    Node *l, *r;
    
    Node(Bisector v) : val(v) {
        static std::default_random_engine gen;
        static std::uniform_int_distribution<unsigned> distribution(0u, ~0u);
        randkey = distribution(gen);
        prev = next = nullptr;
        l = r = nullptr;
    }
};

struct Tree {
    Node *root;
    Node *_le, *_ri;

    Node *_lrotate(Node *x) {
        assert(x->l != nullptr);
        Node *y = x->l;
        x->l = y->r;
        y->r = x;
        return y;
    }
    Node *_rrotate(Node *x) {
        assert(x->r != nullptr);
        Node *y = x->r;
        x->r = y->l;
        y->l = x;
        return y;
    }
    
    Node *_insert(Node *cur, Node *tmp) {
        if (cur == nullptr)
            return tmp;
        if (tmp->val < cur->val) {
            _ri = cur;
            cur->l = _insert(cur->l, tmp);
            if (cur->l->randkey > cur->randkey)
                cur = _lrotate(cur);
        }
        else {
            _le = cur;
            cur->r = _insert(cur->r, tmp);
            if (cur->r->randkey > cur->randkey)
                cur = _rrotate(cur);
        }
        return cur;
    }
    Node *_remove(Node *cur, Bisector val) {
        assert(cur != nullptr);
        if (cur->val == val) {
            if (cur->l == nullptr && cur->r == nullptr) {
                _le = cur->prev, _ri = cur->next;
                //std::cerr << "removed [" << cur->val.i << " " << cur->val.j << "]" << std::endl;
                delete cur;
                cur = nullptr;
            }
            else {
                if (cur->l != nullptr && (cur->r == nullptr || cur->l->randkey > cur->r->randkey)) {
                    cur = _lrotate(cur);
                    cur->r = _remove(cur->r, val);
                }
                else {
                    cur = _rrotate(cur);
                    cur->l = _remove(cur->l, val);
                }
            }
        }
        else {
            if (val < cur->val)
                cur->l = _remove(cur->l, val);
            else
                cur->r = _remove(cur->r, val);
        }
        return cur;
    }
    void _print(Node *cur) const {
        if (cur == nullptr)
            return;
        _print(cur->l);
        //std::cerr << "-- bisector: [" << cur->val.i << " " << cur->val.j << "] ";
        //std::cerr << std::endl;
        _print(cur->r);
    }
public:
    Tree() {
        clear();
    }
    void clear() {
        root = nullptr; // FIXME: free previous nodes
    }
    void insert(Bisector val) {
        Node *tmp = new Node(val);
        _le = _ri = nullptr;
        root = _insert(root, tmp);
        tmp->prev = _le;
        tmp->next = _ri;
        if (_le != nullptr) _le->next = tmp;
        if (_ri != nullptr) _ri->prev = tmp;
    }
    std::pair<Node *, Node *> remove(Bisector val) {
        // [return] its neighbors before removing
        _le = _ri = nullptr;
        root = _remove(root, val);
        if (_le != nullptr) _le->next = _ri;
        if (_ri != nullptr) _ri->prev = _le;
        return std::make_pair(_le, _ri);
    }
    Node *find_left(double xpos) {
        Node *result = nullptr;
        for (Node *cur = root; cur; ) {
            if (cur->val.after(xpos)) {
                cur = cur->l;
            }
            else {
                result = cur;
                cur = cur->r;
            }
        }
        return result;
    }
    void print() const {
        _print(root);
    }
};

Tree tree;

// }}}

//-----------------------------------------------
//-- Queue for Events
//-----------------------------------------------
// {{{

struct Event {
    // if is_site:
    //   p = site position
    //   i = site id
    //   j, k = nothing
    // if !is_site, i.e., is intersection
    //   p = intersection
    //   i, j, k = intersecting sites, in that order as in list
    bool is_site;
    Point p;
    int i, j, k;
    static Event from_site(int i) {
        Event e;
        e.is_site = true;
        e.p = sites[i];
        e.i = i;
        return e;
    }
    static Event from_intersection(int i, int j, int k, Point p) {
        Event e;
        e.is_site = false;
        e.p = p;
        e.i = i, e.j = j, e.k = k;
        return e;
    }
    bool operator <(const Event &e) const {
        if (p.y != e.p.y)
            return p.y < e.p.y;
        if (is_site != e.is_site) // intersection before site
            return !is_site;
        if (is_site) {
            if (i != e.i)
                return i < e.i;
        }
        else {
            if (i != e.i)
                return i < e.i;
            if (j != e.j)
                return j < e.j;
            if (k != e.k)
                return k < e.k;
        }
        return false;
    }
};

std::multiset<Event> events;

Event extract_event() {
    Event e = *events.begin();
    events.erase(events.begin());
    return e;
}

void add_intersection(int i, int j, int k) {
    Maybe<Point> p = intersect(i, j, k);
    if (p.ok) {
        events.insert(Event::from_intersection(i, j, k, p.get()));
    }
}


void remove_intersection(int i, int j, int k) {
    Maybe<Point> p = intersect(i, j, k);
    if (p.ok) {
        events.erase(events.find(Event::from_intersection(i, j, k, p.get())));
    }
}

// }}}

//-----------------------------------------------
//-- Main
//-----------------------------------------------

struct VoronoiVertex {
    Point p;
    int i, j, k;
    VoronoiVertex(int i, int j, int k) : i(i), j(j), k(k) {
        p = raw_intersect(i, j, k).get();
    }
};
std::vector<VoronoiVertex> voronoi_vertices;

// {{{
void solve() {
    voronoi_vertices.clear();

    std::sort(sites + 1, sites + n + 1);
    events.clear();
    for (int i = 1; i <= n; ++i)
        events.insert(Event::from_site(i));

    int first_site = extract_event().i;
    sweepline = sites[first_site].y;
    tree.clear();
    tree.insert(Bisector(0, first_site));
    tree.insert(Bisector(first_site, 0));

    while (!events.empty()) {
        Event e = extract_event();
        sweepline = e.p.y;
        //std::cerr.precision(15);
        //std::cerr << "sweepline: " << sweepline << std::endl;
        //std::cerr << "event point: " << e.p.x << " " << e.p.y << std::endl;

        if (e.is_site) {
            //std::cerr << "site event\n";
            int site_id = e.i;
            //std::cerr << "site: " << site_id << std::endl;
            Node *c = tree.find_left(sites[site_id].x);
            assert(c->next != nullptr);
            assert(c->val.j == c->next->val.i);
            remove_intersection(c->val.i, c->val.j, c->next->val.j);
            add_intersection(c->val.i, c->val.j, site_id);
            add_intersection(site_id, c->next->val.i, c->next->val.j);
            tree.insert(Bisector(c->val.j, site_id));
            tree.insert(Bisector(site_id, c->next->val.i));
        }
        else {
            //std::cerr << "intersection event\n";
            voronoi_vertices.push_back(VoronoiVertex(e.i, e.j, e.k));

            Node *x = tree.remove(Bisector(e.i, e.j)).first;
        //std::cerr << "removed [" << e.i << " " << e.j << "]; tree:\n"; tree.print();
            Node *y = tree.remove(Bisector(e.j, e.k)).second;
        //std::cerr << "removed [" << e.j << " " << e.k << "]; tree:\n"; tree.print();
            assert(x != nullptr && x->val.j == e.i);
            assert(y != nullptr && y->val.i == e.k);
            remove_intersection(x->val.i, e.i, e.j);
            remove_intersection(e.j, e.k, y->val.j);
            add_intersection(x->val.i, e.i, e.k);
            add_intersection(e.i, e.k, y->val.j);
            tree.insert(Bisector(e.i, e.k));
        }
        //std::cerr << ""; tree.print();
        //std::cerr << "===========================" << std::endl;
    }

    //for (auto &vp : voronoi_vertices)
    //    printf("%.9lf %.9lf\n", vp.p.x, vp.p.y);
}
// }}}

// {{{

const double boundary_l = -9.2;
const double boundary_r = 1.2;
const double boundary_x[4] = { boundary_l, boundary_l, boundary_r, boundary_r };
const double boundary_y[4] = { boundary_l, boundary_r, boundary_r, boundary_l };

bool inside_boundary(Point p) {
    const double eps = 1e-7;
    if (p.x + eps < boundary_l || p.y + eps < boundary_l)
        return false;
    if (p.x - eps > boundary_r || p.y - eps > boundary_r)
        return false;
    return true;
}

Point translate_to_boundary(Point p0, Point dir) {
    Point target = p0;
    double smallest_c = 1e50;
    for (int t = 0; t < 4; ++t) {
        double boundary = (t & 1) ? boundary_l : boundary_r;
        double diff, velocity;
        if (t & 2) {
            diff = boundary - p0.x;
            velocity = dir.x;
        }
        else {
            diff = boundary - p0.y;
            velocity = dir.y;
        }
        if (fabs(velocity) > 1e-10) {
            double c = diff / velocity;
            if (c > 1e-10 && c < smallest_c) {
                Point tmp = p0 + Point(dir.x * c, dir.y * c);
                bool flag = inside_boundary(tmp);
                if (flag) {
                    target = tmp;
                    smallest_c = c;
                }
            }
        }
    }
    return target;
}

void translate_results() {
    for (int t = 0; t < 4; ++t) {
        printf("%.2lf %.2lf\n", boundary_x[t], boundary_y[t]);
        int z = (t + 1) % 4;
        printf("%.2lf %.2lf\n", boundary_x[z], boundary_y[z]);
    }
    std::map<std::pair<int, int>, Point> s;
    for (auto &vp : voronoi_vertices) {
        int ver[3] = { vp.i, vp.j, vp.k };
        for (int cur = 0; cur < 3; ++cur) {
            int i = ver[cur], j = ver[(cur + 1) % 3];
            auto iter = s.find(std::make_pair(j, i));
            if (iter == s.end()) {
                s[std::make_pair(i, j)] = vp.p;
            }
            else {
                Point p = vp.p, q = iter->second;
                s.erase(iter);
                while (true) {
                    bool f = false;
                    if (!inside_boundary(p)) {
                        p = translate_to_boundary(p, q - p);
                        f = true;
                    }
                    if (!inside_boundary(q)) {
                        q = translate_to_boundary(q, p - q);
                        f = true;
                    }
                    if (!f) break;
                }
                printf("%.9lf %.9lf\n%.9lf %.9lf\n", p.x, p.y, q.x, q.y);
            }
        }
    }
    if (n == 2) {
        Point pi = sites[1];
        Point pj = sites[2];
        Point p = pi + pj; p.x /= 2, p.y /= 2;
        Point dir(pj.y - pi.y, pi.x - pj.x);
        Point q = translate_to_boundary(p, dir);
        p = translate_to_boundary(p, p - q);
        printf("%.9lf %.9lf\n%.9lf %.9lf\n", p.x, p.y, q.x, q.y);
    }
    for (auto &sp : s) {
        Point p = sp.second;
        if (!inside_boundary(p))
            continue;
        Point pi = sites[sp.first.first];
        Point pj = sites[sp.first.second];
        Point dir(pj.y - pi.y, pi.x - pj.x);
        Point q = translate_to_boundary(p, dir);
        printf("%.9lf %.9lf\n%.9lf %.9lf\n", p.x, p.y, q.x, q.y);
    }
}
// }}}

int main(int argc, char *argv[]) {
    if (argc >= 2) {
        freopen(argv[1], "r", stdin);
    }
    if (argc >= 3) {
        freopen(argv[2], "w", stdout);
    }
    read_sites();
    solve();
    translate_results();
    return 0;
}

