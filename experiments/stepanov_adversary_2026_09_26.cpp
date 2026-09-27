// Stepanov wave, worker `adversary` (2026-09-26): search helper for robust Hanson--Petridis.
//
// Build:  clang++ -O2 -std=c++17 -o stepanov_adversary experiments/stepanov_adversary_2026_09_26.cpp
//
// Notation (research/stepanov-brief-2026-09-26.md): chi = quadratic character, chi(0) = 0,
// for A subset of the field, e_b = #{a in A : chi(a+b) = -1}, B_e(A) = {b : e_b <= e},
// r_e = |B_e(A) cap (-A)|, obj_e(A) = |A|*|B_e(A)| - r_e, R_e(A) = obj_e(A)/((q-1)/2).
//
// Modes (all output JSON, one object per line, on stdout):
//   exh  <p> <mmax> <threads>       exact max of obj_e over all A with |A| = m <= mmax, every e,
//                                   over F_p (p <= 61); A normalised to contain {0,1} or {0,nu}.
//   exhq <p> <mmax> <threads>       the same over F_{p^2} = F_p[x]/(x^2 - nu) (p^2 <= 64).
//   rhp  <fp|fq> <p> <m> <e> <family> <seed> <seconds> [A=a1,a2,...]
//                                   local search (annealing + exact steepest polish) maximising
//                                   obj_e(A), |A| = m, started from a structured seed family.
//   rect <fp|fq> <p> <m> <n> <family> <seed> <seconds>
//                                   alternating maximisation + iterated local search of
//                                   S(A,B) = sum chi(a+b) over |A| = m, |B| = n.
// Field elements of F_{p^2} are indexed u + p*v for u + v*x.
// Every witness printed here is re-checked by exact arithmetic in
// experiments/stepanov_adversary_2026_09_26.py; nothing here is trusted by itself.

#include <cstdio>
#include <climits>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <sstream>

using namespace std;
typedef long long ll;
typedef unsigned long long ull;

// Budgets are in process CPU seconds (the machine is heavily shared; wall time is not a work measure).
static double now_sec() {
    return (double)clock() / CLOCKS_PER_SEC;
}
static double wall_sec() {
    return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();
}

struct Field {
    int p = 0, q = 0;
    bool ext = false;
    int nu = 0;      // non-residue mod p (defines F_{p^2})
    int nonsq = 0;   // index of a non-square of the field
    vector<int8_t> chi, chi2;
    vector<int> negt;
    vector<int> addt;  // only for ext
    inline int add(int a, int b) const {
        if (!ext) { int s = a + b; return s >= p ? s - p : s; }
        return addt[(size_t)a * q + b];
    }
    inline int chisum(int a, int b) const { return ext ? chi[addt[(size_t)a * q + b]] : chi2[a + b]; }
    inline int neg(int a) const { return negt[a]; }
    inline int mul(int a, int b) const {  // only used for F_p
        return (int)((ll)a * b % p);
    }
};

static Field make_fp(int p) {
    Field F; F.p = p; F.q = p; F.ext = false;
    F.chi.assign(p, -1); F.chi[0] = 0;
    for (ll x = 1; x < p; x++) F.chi[(x * x) % p] = 1;
    F.chi2.resize(2 * p);
    for (int k = 0; k < 2 * p; k++) F.chi2[k] = F.chi[k % p];
    F.negt.resize(p);
    for (int a = 0; a < p; a++) F.negt[a] = a ? p - a : 0;
    for (int x = 2; x < p; x++) if (F.chi[x] == -1) { F.nonsq = x; break; }
    F.nu = F.nonsq;
    return F;
}

// Null model: F_p with the Legendre symbol replaced by a balanced random +-1 function
// (value 0 at 0) read from a file of p integers; the driver generates the table from a
// Python seed so that the verifier can regenerate it exactly.
static Field make_fr(int p, const string& path) {
    Field F = make_fp(p);
    FILE* fh = fopen(path.c_str(), "r");
    if (!fh) { fprintf(stderr, "cannot open %s\n", path.c_str()); exit(1); }
    for (int x = 0; x < p; x++) { int v; if (fscanf(fh, "%d", &v) != 1) { fprintf(stderr, "bad table\n"); exit(1); } F.chi[x] = (int8_t)v; }
    fclose(fh);
    for (int k = 0; k < 2 * p; k++) F.chi2[k] = F.chi[k % p];
    for (int x = 1; x < p; x++) if (F.chi[x] == -1) { F.nonsq = x; break; }
    F.nu = F.nonsq;
    return F;
}

static Field make_fq(int p);

static Field make_field(const string& fk, int p) {
    if (fk == "fp") return make_fp(p);
    if (fk == "fq") return make_fq(p);
    if (fk.rfind("fr:", 0) == 0) return make_fr(p, fk.substr(3));
    fprintf(stderr, "unknown field kind\n"); exit(1);
}

static Field make_fq(int p) {
    Field F; F.p = p; F.q = p * p; F.ext = true;
    vector<int8_t> chip(p, -1); chip[0] = 0;
    for (ll x = 1; x < p; x++) chip[(x * x) % p] = 1;
    int nu = 0;
    for (int x = 2; x < p; x++) if (chip[x] == -1) { nu = x; break; }
    F.nu = nu;
    int q = p * p;
    F.chi.resize(q);
    for (int v = 0; v < p; v++)
        for (int u = 0; u < p; u++) {
            int idx = u + p * v;
            if (idx == 0) { F.chi[idx] = 0; continue; }
            ll N = ((ll)u * u - (ll)nu * v % p * v) % p;
            if (N < 0) N += p;
            F.chi[idx] = chip[N];  // chi_q(z) = chi_p(Norm z)
        }
    F.addt.resize((size_t)q * q);
    for (int a = 0; a < q; a++) {
        int ua = a % p, va = a / p;
        for (int b = 0; b < q; b++) {
            int ub = b % p, vb = b / p;
            F.addt[(size_t)a * q + b] = (ua + ub) % p + p * ((va + vb) % p);
        }
    }
    F.negt.resize(q);
    for (int a = 0; a < q; a++) { int u = a % p, v = a / p; F.negt[a] = ((p - u) % p) + p * ((p - v) % p); }
    for (int x = 1; x < q; x++) if (F.chi[x] == -1) { F.nonsq = x; break; }
    return F;
}

static string vec_json(const vector<int>& v) {
    string s = "[";
    for (size_t i = 0; i < v.size(); i++) { if (i) s += ","; s += to_string(v[i]); }
    return s + "]";
}

// ---------------------------------------------------------------- exhaustive (q <= 64)

struct ExhBest { ll obj = -1; int Bsz = 0, r = 0; vector<int> A; ll count = 0; };
static const int NPL = 5;

struct ExhWorker {
    const Field* F; int q, mmax; ull full;
    const vector<ull>* bad; const vector<ull>* negbit;
    vector<vector<ExhBest>> best;
    vector<int> cand;
    int A[64]; int msz = 0;
    ll nodes = 0;
    void eval(const ull* P, ull negA) {
        nodes++;
        int m = msz; ull le = 0;
        for (int e = 0; e <= m; e++) {
            ull eq = full;
            for (int k = 0; k < NPL; k++) eq &= ((e >> k) & 1) ? P[k] : ~P[k];
            le |= eq;
            int Bsz = __builtin_popcountll(le), r = __builtin_popcountll(le & negA);
            ll obj = (ll)m * Bsz - r;
            ExhBest& bb = best[m][e];
            if (obj > bb.obj) { bb.obj = obj; bb.Bsz = Bsz; bb.r = r; bb.A.assign(A, A + m); bb.count = 1; }
            else if (obj == bb.obj) bb.count++;
        }
    }
    void dfs(int start, const ull* P, ull negA) {
        eval(P, negA);
        if (msz == mmax) return;
        for (int i = start; i < (int)cand.size(); i++) {
            int c = cand[i];
            ull Q[NPL]; memcpy(Q, P, sizeof(Q));
            ull carry = (*bad)[c];
            for (int k = 0; k < NPL; k++) { ull t = Q[k] & carry; Q[k] ^= carry; carry = t; }
            A[msz++] = c;
            dfs(i + 1, Q, negA | (*negbit)[c]);
            msz--;
        }
    }
};

static void run_exh(const Field& F, int mmax, int nthreads, const char* tag) {
    int q = F.q;
    if (q > 64 || mmax > 31) { fprintf(stderr, "exh: q<=64, mmax<=31 required\n"); exit(1); }
    double t0 = wall_sec(), c0 = now_sec();
    ull full = (q == 64) ? ~0ULL : ((1ULL << q) - 1);
    vector<ull> bad(q, 0), negbit(q, 0);
    for (int c = 0; c < q; c++) {
        for (int b = 0; b < q; b++) if (F.chisum(c, b) == -1) bad[c] |= 1ULL << b;
        negbit[c] = 1ULL << F.neg(c);
    }
    // global best incl. m = 1 (A = {0}; all singletons equivalent under translation)
    vector<vector<ExhBest>> gbest(mmax + 1, vector<ExhBest>(mmax + 1));
    {
        ExhWorker w; w.F = &F; w.q = q; w.mmax = mmax; w.full = full; w.bad = &bad; w.negbit = &negbit;
        w.best.assign(mmax + 1, vector<ExhBest>(mmax + 1));
        ull P[NPL] = {0, 0, 0, 0, 0};
        ull carry = bad[0];
        for (int k = 0; k < NPL; k++) { ull t = P[k] & carry; P[k] ^= carry; carry = t; }
        w.A[0] = 0; w.msz = 1; w.eval(P, negbit[0]);
        for (int e = 0; e <= 1 && e <= mmax; e++) gbest[1][e] = w.best[1][e];
    }
    int svals[2] = {1, F.nonsq};
    // tasks: (s, i) with first extra element cand[i]; plus (s, -1) = the pair itself
    vector<pair<int, int>> tasks;
    vector<vector<int>> cands(2);
    for (int si = 0; si < 2; si++) {
        for (int c = 0; c < q; c++) if (c != 0 && c != svals[si]) cands[si].push_back(c);
        tasks.push_back({si, -1});
        if (mmax >= 3) for (int i = 0; i < (int)cands[si].size(); i++) tasks.push_back({si, i});
    }
    atomic<int> next(0);
    mutex mu;
    ll total_nodes = 0;
    auto worker = [&]() {
        ExhWorker w; w.F = &F; w.q = q; w.mmax = mmax; w.full = full; w.bad = &bad; w.negbit = &negbit;
        w.best.assign(mmax + 1, vector<ExhBest>(mmax + 1));
        while (true) {
            int t = next.fetch_add(1);
            if (t >= (int)tasks.size()) break;
            int si = tasks[t].first, i = tasks[t].second, s = svals[si];
            w.cand = cands[si];
            ull P[NPL] = {0, 0, 0, 0, 0};
            for (int c : {0, s}) {
                ull carry = bad[c];
                for (int k = 0; k < NPL; k++) { ull tt = P[k] & carry; P[k] ^= carry; carry = tt; }
            }
            ull negA = negbit[0] | negbit[s];
            w.A[0] = 0; w.A[1] = s; w.msz = 2;
            if (i < 0) { w.eval(P, negA); continue; }
            int c = w.cand[i];
            ull carry = bad[c];
            for (int k = 0; k < NPL; k++) { ull tt = P[k] & carry; P[k] ^= carry; carry = tt; }
            w.A[w.msz++] = c;
            w.dfs(i + 1, P, negA | negbit[c]);
            w.msz--;
        }
        lock_guard<mutex> lk(mu);
        total_nodes += w.nodes;
        for (int m = 2; m <= mmax; m++)
            for (int e = 0; e <= m; e++) {
                ExhBest& g = gbest[m][e]; ExhBest& b = w.best[m][e];
                if (b.obj < 0) continue;
                if (b.obj > g.obj) { ll c = b.count; g = b; g.count = c; }
                else if (b.obj == g.obj) { g.count += b.count; if (b.A < g.A) { g.A = b.A; g.Bsz = b.Bsz; g.r = b.r; } }
            }
    };
    vector<thread> th;
    for (int i = 0; i < nthreads; i++) th.emplace_back(worker);
    for (auto& t : th) t.join();
    double secs = wall_sec() - t0, csecs = now_sec() - c0;
    for (int m = 1; m <= mmax; m++)
        for (int e = 0; e <= m; e++) {
            ExhBest& g = gbest[m][e];
            if (g.obj < 0) continue;
            printf("{\"mode\":\"%s\",\"p\":%d,\"q\":%d,\"m\":%d,\"e\":%d,\"obj\":%lld,\"Bsize\":%d,\"r\":%d,"
                   "\"A\":%s,\"n_normalised_attaining\":%lld}\n",
                   tag, F.p, q, m, e, g.obj, g.Bsz, g.r, vec_json(g.A).c_str(), g.count);
        }
    printf("{\"mode\":\"%s_done\",\"p\":%d,\"q\":%d,\"mmax\":%d,\"nodes\":%lld,\"seconds\":%.2f,\"cpu_seconds\":%.2f,\"nonsq\":%d}\n",
           tag, F.p, q, mmax, total_nodes, secs, csecs, F.nonsq);
    fflush(stdout);
}

// ---------------------------------------------------------------- helpers for seeds

static vector<int> random_subset(int q, int m, mt19937_64& rng, const vector<int>* pool = nullptr) {
    vector<int> out;
    if (pool) {
        vector<int> v = *pool;
        shuffle(v.begin(), v.end(), rng);
        for (int i = 0; i < m && i < (int)v.size(); i++) out.push_back(v[i]);
        return out;
    }
    vector<char> used(q, 0);
    while ((int)out.size() < m) { int x = rng() % q; if (!used[x]) { used[x] = 1; out.push_back(x); } }
    return out;
}

static void fill_to(vector<int>& A, int q, int m, mt19937_64& rng) {
    vector<char> used(q, 0);
    vector<int> B;
    for (int x : A) if (!used[x]) { used[x] = 1; B.push_back(x); }
    A = B;
    if ((int)A.size() > m) { shuffle(A.begin(), A.end(), rng); A.resize(m); return; }
    while ((int)A.size() < m) { int x = rng() % q; if (!used[x]) { used[x] = 1; A.push_back(x); } }
}

static int mulorder(ll g, int p) {
    ll x = g % p; int k = 1;
    while (x != 1) { x = x * g % p; k++; if (k > p) return -1; }
    return k;
}

static vector<int> greedy_clique(const Field& F, int kmax, mt19937_64& rng) {
    // clique of the difference graph: chi(x - y) = 1; greedy over random order
    int q = F.q;
    vector<int> order(q); iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);
    vector<int> K;
    for (int x : order) {
        bool ok = true;
        for (int y : K) if (F.chisum(x, F.neg(y)) != 1) { ok = false; break; }
        if (ok) { K.push_back(x); if ((int)K.size() >= kmax) break; }
    }
    return K;
}

static vector<int> seed_m4(const Field& F, int m, mt19937_64& rng, int iters);
static vector<int> seed_rect(const Field& F, int m, int n, mt19937_64& rng, double seconds);

static vector<int> make_seed(const Field& F, const string& fam, int m, mt19937_64& rng, const vector<int>& given) {
    int q = F.q, p = F.p;
    vector<int> A;
    if (fam == "random" || (F.ext && fam != "given" && fam != "m4" && fam != "rect" && fam != "clique" && fam != "nbhd" && fam != "b01")) {
        return random_subset(q, m, rng);
    }
    if (fam == "given") { A = given; fill_to(A, q, m, rng); return A; }
    if (fam == "interval") {
        int c = (rng() & 1) ? 1 : F.nonsq; int t = rng() % p;
        for (int i = 0; i < m; i++) A.push_back((int)(((ll)c * i + t) % p));
        return A;
    }
    if (fam == "gp") {
        for (int tries = 0; tries < 1000; tries++) {
            int g = 2 + rng() % (p - 3);
            int o = mulorder(g, p);
            if (o < m) continue;
            int c = (rng() & 1) ? 1 : F.nonsq;
            ll x = c;
            A.clear();
            for (int i = 0; i < m; i++) { A.push_back((int)x); x = x * g % p; }
            return A;
        }
        return random_subset(q, m, rng);
    }
    if (fam == "coset") {
        // union of cosets of a subgroup H of order h | p-1, h <= m (prefer h = m, then large h)
        vector<int> divs;
        for (int h = 2; h <= m && h < p - 1; h++) if ((p - 1) % h == 0) divs.push_back(h);
        if (divs.empty()) return random_subset(q, m, rng);
        int h;
        if (find(divs.begin(), divs.end(), m) != divs.end() && (rng() % 2 == 0)) h = m;
        else h = divs[rng() % divs.size()];
        // generator of the order-h subgroup: g^((p-1)/h) for a primitive root g
        int g = 2;
        while (mulorder(g, p) != p - 1) g++;
        ll gen = 1, base = g; ll ex = (p - 1) / h;
        { ll b = base, e = ex; gen = 1; while (e) { if (e & 1) gen = gen * b % p; b = b * b % p; e >>= 1; } }
        vector<char> used(p, 0);
        while ((int)A.size() < m) {
            int c = 1 + rng() % (p - 1);
            if (used[c]) continue;
            ll x = c;
            for (int i = 0; i < h; i++) { if (!used[x]) { used[x] = 1; A.push_back((int)x); } x = x * gen % p; }
        }
        if ((int)A.size() > m) { shuffle(A.begin(), A.end(), rng); A.resize(m); }
        return A;
    }
    if (fam == "clique") {  // near-clique: greedy clique, filled from its common neighbourhood
        vector<int> K = greedy_clique(F, m, rng);
        A = K;
        if ((int)A.size() < m) {
            vector<int> pool;
            vector<char> inK(q, 0); for (int x : K) inK[x] = 1;
            for (int x = 0; x < q; x++) {
                if (inK[x]) continue;
                int bad = 0;
                for (int y : K) if (F.chisum(x, F.neg(y)) != 1) bad++;
                if (bad <= 1) pool.push_back(x);
            }
            shuffle(pool.begin(), pool.end(), rng);
            for (int x : pool) { if ((int)A.size() >= m) break; A.push_back(x); }
        }
        fill_to(A, q, m, rng);
        return A;
    }
    if (fam == "nbhd") {  // A inside N[K] for a clique K, so that -K lies in B_0(A)
        int kmax = max(1, (int)floor(log2((double)q / m)));  // keeps |N[K]| ~ q/2^k >= m
        int k = 1 + rng() % kmax;
        vector<int> K = greedy_clique(F, k, rng);
        vector<int> pool;
        for (int x = 0; x < q; x++) {
            bool ok = true;
            for (int y : K) if (F.chisum(x, F.neg(y)) == -1) { ok = false; break; }
            if (ok) pool.push_back(x);
        }
        if ((int)pool.size() >= m) return random_subset(q, m, rng, &pool);
        A = pool; fill_to(A, q, m, rng); return A;
    }
    if (fam == "b01") {  // A inside B_0({0,1}) = {x : chi(x), chi(x+1) != -1}
        vector<int> pool;
        for (int x = 0; x < q; x++) if (F.chisum(x, 0) != -1 && F.chisum(x, 1) != -1) pool.push_back(x);
        return random_subset(q, m, rng, &pool);
    }
    if (fam == "b01core") {  // {0,1} plus elements of B_0({0,1})
        vector<int> pool;
        for (int x = 2; x < q; x++) if (F.chisum(x, 0) != -1 && F.chisum(x, 1) != -1) pool.push_back(x);
        A = {0, 1};
        shuffle(pool.begin(), pool.end(), rng);
        for (int x : pool) { if ((int)A.size() >= m) break; A.push_back(x); }
        fill_to(A, q, m, rng); return A;
    }
    if (fam == "squares") {
        vector<char> used(q, 0);
        for (ll k = 1; (int)A.size() < m && k < p; k++) { int x = (int)(k * k % p); if (!used[x]) { used[x] = 1; A.push_back(x); } }
        return A;
    }
    if (fam == "m4") return seed_m4(F, m, rng, 300);
    if (fam == "rect") return seed_rect(F, m, max(m, (q - 1) / (2 * m)), rng, 0.25);
    return random_subset(q, m, rng);
}

// ---------------------------------------------------------------- fourth-moment seed

static vector<int> seed_m4(const Field& F, int m, mt19937_64& rng, int iters) {
    int q = F.q;
    vector<int> A = random_subset(q, m, rng);
    vector<char> inA(q, 0); for (int a : A) inA[a] = 1;
    vector<int> FA(q, 0);
    for (int a : A) for (int b = 0; b < q; b++) FA[b] += F.chisum(a, b);
    for (int it = 0; it < iters; it++) {
        int io = rng() % m; int ao = A[io];
        int bestin = -1; double bestd = 0;
        for (int c = 0; c < 8; c++) {
            int ai = rng() % q; if (inA[ai]) continue;
            double d = 0;
            for (int b = 0; b < q; b++) {
                ll f = FA[b], g = f - F.chisum(ao, b) + F.chisum(ai, b);
                d += (double)(g * g * g * g - f * f * f * f);
            }
            if (d > bestd) { bestd = d; bestin = ai; }
        }
        if (bestin >= 0) {
            for (int b = 0; b < q; b++) FA[b] += F.chisum(bestin, b) - F.chisum(ao, b);
            inA[ao] = 0; inA[bestin] = 1; A[io] = bestin;
        }
    }
    return A;
}

// ---------------------------------------------------------------- rectangle bias search

struct RectRes { ll S = LLONG_MIN; vector<int> A, B; ll rounds = 0; };

static void topk(const vector<int>& val, int k, vector<int>& out, mt19937_64& rng) {
    int q = val.size();
    vector<pair<ll, int>> key(q);
    for (int x = 0; x < q; x++) key[x] = {-(((ll)val[x]) << 20) - (ll)(rng() & 0xFFFFF), x};
    nth_element(key.begin(), key.begin() + k, key.end());
    out.resize(k);
    for (int i = 0; i < k; i++) out[i] = key[i].second;
}

static RectRes rect_search(const Field& F, int m, int n, mt19937_64& rng, double seconds,
                           const vector<int>* initA = nullptr) {
    int q = F.q;
    double t0 = now_sec();
    RectRes best;
    vector<int> A, B, FA(q), GB(q);
    auto compFA = [&]() { fill(FA.begin(), FA.end(), 0); for (int a : A) for (int b = 0; b < q; b++) FA[b] += F.chisum(a, b); };
    auto compGB = [&]() { fill(GB.begin(), GB.end(), 0); for (int b : B) for (int a = 0; a < q; a++) GB[a] += F.chisum(a, b); };
    int restart = 0;
    while (true) {
        if (restart == 0 && initA) A = *initA, fill_to(A, q, m, rng);
        else if (restart == 0 || best.A.empty() || rng() % 3 == 0) A = random_subset(q, m, rng);
        else {  // perturb the best
            A = best.A;
            int t = max(1, m / 5);
            vector<char> used(q, 0); for (int a : A) used[a] = 1;
            for (int i = 0; i < t; i++) {
                int j = rng() % m; int x;
                do { x = rng() % q; } while (used[x]);
                used[A[j]] = 0; used[x] = 1; A[j] = x;
            }
        }
        restart++;
        compFA();
        ll S = LLONG_MIN;
        for (int it = 0; it < 60; it++) {
            topk(FA, n, B, rng);
            compGB();
            vector<int> A2; topk(GB, m, A2, rng);
            A = A2; compFA();
            ll S2 = 0; for (int b : B) S2 += FA[b];
            // B is optimal for the new A only after re-selection; evaluate with re-selected B
            vector<int> B2; topk(FA, n, B2, rng);
            ll S3 = 0; for (int b : B2) S3 += FA[b];
            B = B2;
            best.rounds++;
            if (S3 > best.S) { best.S = S3; best.A = A; best.B = B; }
            if (S3 <= S) break;
            S = S3;
        }
        if (now_sec() - t0 > seconds) break;
    }
    return best;
}

static vector<int> seed_rect(const Field& F, int m, int n, mt19937_64& rng, double seconds) {
    RectRes r = rect_search(F, m, min(n, F.q - 1), rng, seconds);
    return r.A;
}

// ---------------------------------------------------------------- RHP local search

struct RHPSearch {
    const Field& F; int q, m, e; mt19937_64& rng;
    vector<int> A; vector<int> pos; vector<int> cnt; int good = 0;
    int K = 2; vector<double> w; vector<int> L, bd;
    vector<int> R0;  // elements s with chi(s) != -1
    RHPSearch(const Field& F_, int m_, int e_, mt19937_64& rng_) : F(F_), q(F_.q), m(m_), e(e_), rng(rng_) {
        for (int s = 0; s < q; s++) if (F.chi[s] != -1) R0.push_back(s);
    }
    inline int bad(int a, int b) const { return F.chisum(a, b) == -1; }
    void load(const vector<int>& A0) {
        A = A0; pos.assign(q, -1);
        for (int i = 0; i < m; i++) pos[A[i]] = i;
        cnt.assign(q, 0);
        for (int a : A) for (int b = 0; b < q; b++) cnt[b] += bad(a, b);
        good = 0; for (int b = 0; b < q; b++) good += (cnt[b] <= e);
    }
    int rcount() const { int r = 0; for (int a : A) r += (cnt[F.neg(a)] <= e); return r; }
    ll obj() const { return (ll)m * good - rcount(); }
    void choose_K(int Ltarget) {
        vector<int> hist(m + 2, 0);
        for (int b = 0; b < q; b++) hist[min(cnt[b], m + 1)]++;
        int acc = 0; for (int j = 0; j <= e; j++) acc += hist[j];
        K = 1;
        int Kcap = max(10, m / 2 - e);  // large m, small e: widen the window so the surrogate is never flat
        while (e + K < m && K < Kcap) { acc += hist[e + K]; if (acc >= Ltarget && K >= 2) break; K++; }
        if (K < 2) K = 2;
        w.assign(m + 3, 0.0);
        double lam = 4.0 / K;
        for (int j = 0; j <= m + 2; j++) {
            if (j <= e) w[j] = 1.0;
            else if (j <= e + K - 1) w[j] = (K <= 10) ? exp(-lam * (j - e)) : 0.5 * exp(-lam * (j - e - 1));
            else w[j] = 0.0;
        }
    }
    void rebuild() {
        L.clear(); bd.clear();
        for (int b = 0; b < q; b++) {
            if (cnt[b] <= e + K) L.push_back(b);
            if (cnt[b] == e + 1) bd.push_back(b);
        }
    }
    double sur() const { double s = 0; for (int b = 0; b < q; b++) s += w[min(cnt[b], m + 2)]; return s; }
    inline double delta(int ao, int ai) const {
        double d = 0;
        for (int b : L) {
            int c = cnt[b], c2 = c - bad(ao, b) + bad(ai, b);
            if (c2 != c) d += w[c2] - w[c];
        }
        return d;
    }
    void apply(int io, int ai) {
        int ao = A[io];
        for (int b = 0; b < q; b++) {
            int dd = bad(ai, b) - bad(ao, b);
            if (dd) {
                int c = cnt[b];
                if (c <= e) good--;
                cnt[b] = c + dd;
                if (c + dd <= e) good++;
            }
        }
        pos[ao] = -1; pos[ai] = io; A[io] = ai;
    }
    // exact steepest ascent on |B_e| (ties broken by the exact objective)
    bool polish_once(double deadline) {
        vector<int> order(m); iota(order.begin(), order.end(), 0);
        shuffle(order.begin(), order.end(), rng);
        vector<int> gcount(q);
        for (int io : order) {
            if (now_sec() > deadline) return false;
            int ao = A[io];
            vector<int> T; int base = 0;
            for (int b = 0; b < q; b++) {
                int c = cnt[b] - bad(ao, b);
                if (c <= e) base++;
                if (c == e) T.push_back(b);
            }
            if ((ll)T.size() * (ll)R0.size() > 20000000LL) continue;
            fill(gcount.begin(), gcount.end(), 0);
            for (int b : T) { int nb = F.neg(b); for (int s : R0) gcount[F.add(s, nb)]++; }
            int bestx = -1; ll bestgood = good;
            for (int x = 0; x < q; x++) {
                if (pos[x] >= 0) continue;
                ll ng = base - (ll)T.size() + gcount[x];
                if (ng > bestgood) { bestgood = ng; bestx = x; }
            }
            if (bestx >= 0) { apply(io, bestx); return true; }
        }
        return false;
    }
};

struct RHPResult { ll obj = -1; int good = 0, r = 0; vector<int> A; ll steps = 0; int restarts = 0; };

static RHPResult rhp_search(const Field& F, int m, int e, const string& fam, mt19937_64& rng, double seconds,
                            const vector<int>& given) {
    double t0 = now_sec(), deadline = t0 + seconds;
    RHPSearch S(F, m, e, rng);
    RHPResult best;
    int restart = 0;
    auto record = [&]() {
        ll o = S.obj();
        if (o > best.obj) { best.obj = o; best.good = S.good; best.r = S.rcount(); best.A = S.A; }
    };
    while (now_sec() < deadline) {
        vector<int> A0;
        if (restart == 0 || best.A.empty() || rng() % 2 == 0) A0 = make_seed(F, fam, m, rng, given);
        else {
            A0 = best.A;
            int t = max(1, m / 4);
            vector<char> used(F.q, 0); for (int a : A0) used[a] = 1;
            for (int i = 0; i < t; i++) {
                int j = rng() % m; int x;
                do { x = rng() % F.q; } while (used[x]);
                used[A0[j]] = 0; used[x] = 1; A0[j] = x;
            }
        }
        if ((int)A0.size() != m) fill_to(A0, F.q, m, rng);
        restart++;
        S.load(A0);
        record();
        S.choose_K(400);
        S.rebuild();
        // annealing on the surrogate
        double rem = deadline - now_sec();
        double budget = min(rem, max(0.5, seconds / 6.0));
        double tstart = now_sec(), T0 = 1.0, T1 = 0.02;
        ll step = 0;
        while (true) {
            double el = now_sec() - tstart;
            if (el > budget * 0.8) break;
            double frac = el / (budget * 0.8);
            double T = T0 * pow(T1 / T0, frac);
            for (int inner = 0; inner < 16; inner++) {
                step++;
                int io = rng() % m, ao = S.A[io];
                int bestai = -1; double bestd = -1e18;
                for (int c = 0; c < 24; c++) {
                    int ai;
                    if (c % 2 == 1 && !S.bd.empty()) {
                        int b = S.bd[rng() % S.bd.size()];
                        int s = S.R0[rng() % S.R0.size()];
                        ai = F.add(s, F.neg(b));
                    } else ai = rng() % F.q;
                    if (S.pos[ai] >= 0) continue;
                    double d = S.delta(ao, ai);
                    if (d > bestd) { bestd = d; bestai = ai; }
                }
                if (bestai < 0) continue;
                if (bestd >= 0 || (double)(rng() % 1000000) / 1e6 < exp(bestd / T)) {
                    S.apply(io, bestai);
                    S.rebuild();
                    record();
                }
            }
        }
        best.steps += step;
        // exact polish
        while (now_sec() < deadline && S.polish_once(deadline)) { record(); }
        record();
    }
    best.restarts = restart;
    return best;
}

// ---------------------------------------------------------------- second-moment share search
// maximise P_thr(A) = sum_{b : F_A(b) >= thr} F_A(b)^2 over |A| = m; share = P/(m(q-m)).
// thr = 1: positive part of the second moment; thr = m - 2e - 1: rows with at most e bad sums.

struct ShareRes { ll obj = -1; vector<int> A; ll steps = 0; };

static ShareRes share_search(const Field& F, int m, int thr, const string& fam, mt19937_64& rng, double seconds) {
    int q = F.q;
    double t0 = now_sec(), deadline = t0 + seconds;
    ShareRes best;
    vector<int> FA(q);
    auto val = [&](int f) -> ll { return f >= thr ? (ll)f * f : 0; };
    int restart = 0;
    while (now_sec() < deadline) {
        vector<int> A;
        if (restart == 0 || best.A.empty() || rng() % 2 == 0) A = make_seed(F, fam, m, rng, {});
        else {
            A = best.A;
            int t = max(1, m / 4);
            vector<char> used(q, 0); for (int a : A) used[a] = 1;
            for (int i = 0; i < t; i++) { int j = rng() % m; int x; do { x = rng() % q; } while (used[x]); used[A[j]] = 0; used[x] = 1; A[j] = x; }
        }
        if ((int)A.size() != m) fill_to(A, q, m, rng);
        restart++;
        vector<char> inA(q, 0); for (int a : A) inA[a] = 1;
        fill(FA.begin(), FA.end(), 0);
        for (int a : A) for (int b = 0; b < q; b++) FA[b] += F.chisum(a, b);
        ll cur = 0; for (int b = 0; b < q; b++) cur += val(FA[b]);
        if (cur > best.obj) { best.obj = cur; best.A = A; }
        double budget = min(deadline - now_sec(), max(0.3, seconds / 4.0));
        double ts = now_sec(), T0 = 0.02 * (double)m * q, T1 = 0.5;
        while (true) {
            double el = now_sec() - ts;
            if (el > budget) break;
            double T = T0 * pow(T1 / T0, el / budget);
            for (int inner = 0; inner < 8; inner++) {
                best.steps++;
                int io = rng() % m, ao = A[io];
                int bi = -1; ll bd = LLONG_MIN;
                for (int c = 0; c < 6; c++) {
                    int ai = rng() % q; if (inA[ai]) continue;
                    ll d = 0;
                    for (int b = 0; b < q; b++) {
                        int f0 = FA[b], f1 = f0 - F.chisum(ao, b) + F.chisum(ai, b);
                        if (f1 != f0) d += val(f1) - val(f0);
                    }
                    if (d > bd) { bd = d; bi = ai; }
                }
                if (bi < 0) continue;
                if (bd >= 0 || (double)(rng() % 1000000) / 1e6 < exp((double)bd / T)) {
                    for (int b = 0; b < q; b++) FA[b] += F.chisum(bi, b) - F.chisum(ao, b);
                    inA[ao] = 0; inA[bi] = 1; A[io] = bi; cur += bd;
                    if (cur > best.obj) { best.obj = cur; best.A = A; }
                }
            }
        }
    }
    return best;
}

// ---------------------------------------------------------------- main

static vector<int> parse_list(const string& s) {
    vector<int> v; string cur;
    for (char c : s) { if (c == ',') { if (!cur.empty()) v.push_back(atoi(cur.c_str())); cur.clear(); } else cur += c; }
    if (!cur.empty()) v.push_back(atoi(cur.c_str()));
    return v;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "see header for usage\n"); return 1; }
    string mode = argv[1];
    if (mode == "exh" || mode == "exhq") {
        int p = atoi(argv[2]), mmax = atoi(argv[3]), th = atoi(argv[4]);
        Field F = (mode == "exh") ? make_fp(p) : make_fq(p);
        run_exh(F, mmax, th, mode.c_str());
        return 0;
    }
    if (mode == "rhp") {
        string fk = argv[2]; int p = atoi(argv[3]), m = atoi(argv[4]), e = atoi(argv[5]);
        string fam = argv[6]; ull seed = strtoull(argv[7], 0, 10); double secs = atof(argv[8]);
        vector<int> given;
        if (argc > 9 && strncmp(argv[9], "A=", 2) == 0) given = parse_list(argv[9] + 2);
        Field F = make_field(fk, p);
        mt19937_64 rng(seed * 1000003ULL + 17);
        double t0 = now_sec();
        RHPResult r = rhp_search(F, m, e, fam, rng, secs, given);
        vector<int> A = r.A; sort(A.begin(), A.end());
        printf("{\"mode\":\"rhp\",\"field\":\"%s\",\"p\":%d,\"q\":%d,\"m\":%d,\"e\":%d,\"family\":\"%s\",\"seed\":%llu,"
               "\"obj\":%lld,\"Bsize\":%d,\"r\":%d,\"A\":%s,\"steps\":%lld,\"restarts\":%d,\"seconds\":%.2f}\n",
               fk.c_str(), p, F.q, m, e, fam.c_str(), seed, r.obj, r.good, r.r, vec_json(A).c_str(), r.steps,
               r.restarts, now_sec() - t0);
        return 0;
    }
    if (mode == "rect") {
        string fk = argv[2]; int p = atoi(argv[3]), m = atoi(argv[4]), n = atoi(argv[5]);
        string fam = argv[6]; ull seed = strtoull(argv[7], 0, 10); double secs = atof(argv[8]);
        Field F = make_field(fk, p);
        mt19937_64 rng(seed * 1000003ULL + 29);
        double t0 = now_sec();
        vector<int> init;
        vector<int>* ip = nullptr;
        if (fam != "random") { init = make_seed(F, fam, m, rng, {}); ip = &init; }
        RectRes r = rect_search(F, m, n, rng, secs, ip);
        vector<int> A = r.A, B = r.B; sort(A.begin(), A.end()); sort(B.begin(), B.end());
        printf("{\"mode\":\"rect\",\"field\":\"%s\",\"p\":%d,\"q\":%d,\"m\":%d,\"n\":%d,\"family\":\"%s\",\"seed\":%llu,"
               "\"S\":%lld,\"A\":%s,\"B\":%s,\"rounds\":%lld,\"seconds\":%.2f}\n",
               fk.c_str(), p, F.q, m, n, fam.c_str(), seed, r.S, vec_json(A).c_str(), vec_json(B).c_str(), r.rounds,
               now_sec() - t0);
        return 0;
    }
    if (mode == "share") {
        string fk = argv[2]; int p = atoi(argv[3]), m = atoi(argv[4]), thr = atoi(argv[5]);
        string fam = argv[6]; ull seed = strtoull(argv[7], 0, 10); double secs = atof(argv[8]);
        Field F = make_field(fk, p);
        mt19937_64 rng(seed * 1000003ULL + 41);
        double t0 = now_sec();
        ShareRes r = share_search(F, m, thr, fam, rng, secs);
        vector<int> A = r.A; sort(A.begin(), A.end());
        printf("{\"mode\":\"share\",\"field\":\"%s\",\"p\":%d,\"q\":%d,\"m\":%d,\"thr\":%d,\"family\":\"%s\",\"seed\":%llu,"
               "\"obj\":%lld,\"share\":%.6f,\"A\":%s,\"steps\":%lld,\"seconds\":%.2f}\n",
               fk.c_str(), p, F.q, m, thr, fam.c_str(), seed, r.obj, (double)r.obj / ((double)m * (F.q - m)),
               vec_json(A).c_str(), r.steps, now_sec() - t0);
        return 0;
    }
    fprintf(stderr, "unknown mode\n");
    return 1;
}
