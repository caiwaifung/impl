// Example problem: URAL 1076
// Time complexity: O(n^3)

#include <iostream>
#include <algorithm>
#include <cassert>
#include <queue>

int w[150][150];
int x[150], y[150], slack[150];
int n;

std::queue<int> queue;
bool vx[150], vy[150];
int lx[150], ly[150], prev[150];

void visit(int i) {
    vx[i] = true;
    for (int j = 0; j < n; ++j) {
        if (vy[j]) continue;
        if (w[i][j] == x[i] + y[j]) {
            vy[j] = true;
            prev[j] = i;
            queue.push(j);
        }
        else {
            int d = x[i] + y[j] - w[i][j];
            if (d < slack[j]) {
                slack[j] = d;
                prev[j] = i;
            }
        }
    }
}

void augment(int j) {
    while (j >= 0) {
        int i = prev[j];
        int j2 = lx[i];
        lx[i] = j, ly[j] = i;
        j = j2;
    }
}

bool bfs() {
    while (!queue.empty()) {
        int j = queue.front();
        queue.pop();
        int i = ly[j];
        if (i < 0) {
            augment(j);
            return true;
        }
        visit(i);
    }
    return false;
}

int solve() {
    std::fill(lx, lx + n, -1);
    std::fill(ly, ly + n, -1);
    std::fill(x, x + n, 0);
    std::fill(y, y + n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            x[i] = std::max(x[i], w[i][j]);
        }
    }
    for (int k = 0; k < n; ++k) {
        std::fill(vx, vx + n, false);
        std::fill(vy, vy + n, false);
        std::fill(slack, slack + n, 0x7fffffff);
        queue = std::queue<int>();
        visit(k);
        while (!bfs()) {
            int d = 0x7fffffff;
            for (int j = 0; j < n; j++) {
                if(!vy[j])
                    d = std::min(d, slack[j]);
            }
            for (int i = 0; i < n; i++) {
                if (vx[i]) x[i] -= d;
            }
            for (int j = 0; j < n; j++) {
                if (vy[j]) y[j] += d;
            }
            for (int j = 0; j < n; j++) {
                slack[j] -= d;
                if (!vy[j] && slack[j] == 0) {
                    vy[j] = true;
                    queue.push(j);
                }
            }
        }
        assert(lx[k] >= 0);
    }
    int ans = 0;
    for (int i = 0; i < n; i++) {
        ans += w[i][lx[i]];
    }
    return ans;
}

int main() {
    int ans = 0;
    std::cin >> n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> w[i][j];
            ans += w[i][j];
        }
    }
    ans -= solve();
    std::cout << ans << std::endl;
    return 0;
}
