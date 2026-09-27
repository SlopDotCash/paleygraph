#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Exact integer certificate data. Analytic interpretation is in the pass-36 proof.
using i64 = std::int64_t;
using i128 = __int128_t;
static std::string decimal(i128 x) {
    if (x == 0) return "0";
    bool negative = x < 0;
    if (negative) x = -x;
    std::string s;
    while (x) { s.push_back(char('0' + x % 10)); x /= 10; }
    if (negative) s.push_back('-');
    std::reverse(s.begin(), s.end()); return s;
}
static i64 power(i64 a, i64 b, i64 p) {
    i64 r = 1;
    for (; b; b >>= 1, a = a * a % p) if (b & 1) r = r * a % p;
    return r;
}
static bool prime(i64 p) {
    if (p < 2) return false;
    for (i64 d = 2; d * d <= p; ++d) if (p % d == 0) return false;
    return true;
}
static std::vector<i64> factors(i64 v) {
    std::vector<i64> out;
    for (i64 d = 2; d * d <= v; ++d) if (v % d == 0) {
        out.push_back(d); while (v % d == 0) v /= d;
    }
    if (v > 1) out.push_back(v);
    return out;
}
static int mobius(int v) {
    int out = 1;
    for (int d = 2; d * d <= v; ++d) if (v % d == 0) {
        v /= d; out = -out; if (v % d == 0) return 0;
    }
    return v > 1 ? -out : out;
}
int main(int argc, char **argv) {
    try {
        const auto started = std::chrono::steady_clock::now();
        if (argc != 4) throw std::runtime_error("usage: program prime even_order truncation");
        const i64 p = std::stoll(argv[1]), n = std::stoll(argv[2]);
        const int L = std::stoi(argv[3]);
        if (p > 10000000 || !prime(p) || p < 3 || n < 2 || n % 2 ||
            (p - 1) % n || L < 1 || L > 4096 || (p - 1) / n > 500000)
            throw std::runtime_error("invalid input or explicit resource cap exceeded");
        if (i128(n) * p * (p + 1) > std::numeric_limits<i64>::max() / 4)
            throw std::runtime_error("exact accumulator range exceeded");
        const i64 q = (p - 1) / n, scale = i64(1) << 40;
        const auto fs = factors(p - 1);
        i64 g = 2;
        for (;; ++g) {
            bool good = true;
            for (i64 f : fs) if (power(g, (p - 1) / f, p) == 1) good = false;
            if (good) break;
        }
        std::vector<i64> C(q, 0), d(q), small_index(L + 1, -1);
        i64 x = 1, j = 0;
        for (i64 t = 0; t < p - 1; ++t) {
            i64 r = std::min(x, p - x);
            C[j] += r * r;
            if (x <= L) small_index[x] = j;
            x = x * g % p;
            if (++j == q) j = 0;
        }
        if (x != 1 || j != 0) throw std::runtime_error("orbit did not close");
        i128 mass = 0, centered_mass = 0;
        i64 D = 0, min_j = 0, max_j = 0;
        for (j = 0; j < q; ++j) {
            mass += C[j];
            d[j] = 12 * C[j] - n * p * (p + 1);
            centered_mass += d[j];
            D = std::max(D, std::abs(d[j]));
            if (C[j] < C[min_j]) min_j = j;
            if (C[j] > C[max_j]) max_j = j;
        }
        if (12 * mass != i128(p) * (p * p - 1) || centered_mass != 0)
            throw std::runtime_error("exact mass identity failed");
        std::vector<i128> sums(q, 0);
        i64 weight_mass = 0, terms = 0;
        for (int k = 1; k <= L; ++k) {
            if (k % p == 0) continue;
            int m = k, t = 0;
            while (m % 2 == 0) { ++t; m /= 2; }
            const int mu = mobius(m);
            const i64 denominator = i64(m) * m * (t ? (i64(1) << (t + 1)) : 1);
            // Round toward zero: the omitted coefficient has the same sign.
            const i64 w = -mu * (scale / denominator);
            if (!w) continue;
            const i64 shift = small_index[k % p];
            if (shift < 0) throw std::runtime_error("missing small-residue coset");
            weight_mass += std::abs(w); ++terms;
            const i64 split = q - shift;
            for (j = 0; j < split; ++j) sums[j] += i128(w) * d[j + shift];
            for (j = split; j < q; ++j) sums[j] += i128(w) * d[j - split];
        }
        i128 largest = 0;
        i64 largest_j = 0;
        for (j = 0; j < q; ++j) {
            i128 v = sums[j] < 0 ? -sums[j] : sums[j];
            if (v > largest) { largest = v; largest_j = j; }
        }
        i128 runner_up = 0;
        i64 runner_up_j = -1;
        for (j = 0; j < q; ++j) if (j != largest_j) {
            i128 v = sums[j] < 0 ? -sums[j] : sums[j];
            if (v > runner_up) { runner_up = v; runner_up_j = j; }
        }
        const i128 benefit = i128(D) * weight_mass - largest;
        if (benefit < 0) throw std::runtime_error("weighted triangle inequality failed");
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        std::cout << "{\"p\":" << p << ",\"n\":" << n << ",\"cosets\":" << q
            << ",\"primitive_generator\":" << g
            << ",\"subgroup_generator\":" << power(g, q, p)
            << ",\"contains_two\":" << (power(2, n, p) == 1 ? "true" : "false")
            << ",\"truncation\":" << L << ",\"weight_scale\":" << scale
            << ",\"nonzero_weights\":" << terms << ",\"weight_mass\":" << weight_mass
            << ",\"distance_denominator\":" << 24 * p * p
            << ",\"max_deviation_numerator\":" << D
            << ",\"minimum_all_H_square_sum\":" << C[min_j]
            << ",\"maximum_all_H_square_sum\":" << C[max_j]
            << ",\"minimum_coset_index\":" << min_j
            << ",\"maximum_coset_index\":" << max_j
            << ",\"weighted_max_coset_index\":" << largest_j
            << ",\"weighted_max_abs_numerator\":\"" << decimal(largest) << "\""
            << ",\"weighted_runner_up_abs_numerator\":\"" << decimal(runner_up) << "\""
            << ",\"weighted_runner_up_coset_index\":" << runner_up_j
            << ",\"benefit_numerator\":\"" << decimal(benefit) << "\""
            << ",\"elapsed_seconds\":" << elapsed << ",\"mass_checked\":true}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
