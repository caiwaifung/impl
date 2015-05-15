// Example problem: URAL 1076
// Time complexity: O(n^4)

#include <iostream>
#include <algorithm>

int w[150][150];
int x[150], y[150];
int n;

bool vx[150], vy[150];
int lx[150], ly[150];

bool find(int i) {
    vx[i] = true;
    for (int j = 0; j < n; j++) {
        if (w[i][j] == x[i] + y[j] && !vy[j]) {
            vy[j] = true;
            if (ly[j] < 0 || find(ly[j])) {
                lx[i] = j; ly[j] = i;
                return true;
            }
        }
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
    for (int k = 0; k < n; ) {
        std::fill(vx, vx + n, false);
        std::fill(vy, vy + n, false);
        if (find(k)) {
            ++k;
            continue;
        }
        int d = 0x7fffffff;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (vx[i] && !vy[j])
                    d = std::min(d, x[i] + y[j] - w[i][j]);
            }
        }
        for (int i = 0; i < n; i++) {
            if (vx[i])
                x[i] -= d;
        }
        for (int i = 0; i < n; i++) {
            if (vy[i])
                y[i] += d;
        }
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
