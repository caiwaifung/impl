#include <algorithm>
#include <iostream>
#include <string>
#include "comm.h"
using namespace std;

struct Node {
    Point vantage;
    Node *l, *r;
    double l1, l2, r1, r2;
};

Node *make_node(vector<Point> arr) {
    if (arr.empty())
        return nullptr;
    Node *cur = new Node();
    cur->vantage = arr.back();
    arr.pop_back();

    vector<double> ds;
    for (auto &p : arr) ds.push_back(p.dis(cur->vantage));
    double median = 0;
    if (!ds.empty()) {
        auto kth = ds.begin() + ds.size() / 2;
        nth_element(ds.begin(), kth, ds.end());
        median = *kth;
    }

    vector<Point> la, ra;
    cur->l1 = cur->r1 = +1e50;
    cur->l2 = cur->r2 = -1e50;
    for (auto &p : arr) {
        double d = p.dis(cur->vantage);
        if (d < median) {
            la.push_back(p);
            cur->l1 = min(cur->l1, d);
            cur->l2 = max(cur->l2, d);
        }
        else {
            ra.push_back(p);
            cur->r1 = min(cur->r1, d);
            cur->r2 = max(cur->r2, d);
        }
    }
    cur->l = make_node(la);
    cur->r = make_node(ra);
    return cur;
}

double ans;
Point query_point;

void query(Node *cur) {
    if (cur == nullptr)
        return;
    double tmp = query_point.dis(cur->vantage);
    ans = min(ans, tmp);

    if (cur->l != nullptr && cur->l1 - ans < tmp && cur->l2 + ans > tmp)
        query(cur->l);
    if (cur->r != nullptr && cur->r1 - ans < tmp && cur->r2 + ans > tmp)
        query(cur->r);
}

int main() {
    vector<Place> cities = read_cities();
    vector<Place> queries = read_queries();

    vector<Point> points(cities.size());
    transform(cities.begin(), cities.end(), points.begin(),
        [](const Place &p) { return p.point(); });
    random_shuffle(points.begin(), points.end());
    Node *root = make_node(points);

    cout.precision(20);
    for (auto &q : queries) {
        ans = 1e50;
        query_point = q.point();

        query(root);

        cout << ans << endl;
    }
    return 0;
}
