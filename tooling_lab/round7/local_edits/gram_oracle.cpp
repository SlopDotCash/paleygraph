// Independent exact contraction backend: synthetic division and weighted Gram.
// This is not neighbour enumeration. Python's direct tiny oracle is separate.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I = __int128_t;

std::string decimal(I x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    std::string out;
    while (x) { out.push_back(char('0' + x % 10)); x /= 10; }
    if (neg) out.push_back('-');
    std::reverse(out.begin(), out.end());
    return out;
}

int main() {
    try {
        int q,n,d;
        if (!(std::cin >> q >> n >> d) || q < 5 || q > 10000000 || q%4 != 1 || n < 1 || n >= q || n > 64 || d < 0 || d > 6 || d > n)
            throw std::runtime_error("invalid parameters");
        for (int a=2; int64_t(a)*a<=q; ++a) if (q%a==0) throw std::runtime_error("q not prime");
        std::vector<int> c(n), mark(q,-1);
        for (int i=0;i<n;++i) {
            if (!(std::cin >> c[i]) || c[i]<0 || c[i]>=q || mark[c[i]]!=-1) throw std::runtime_error("invalid selected set");
            mark[c[i]]=i;
        }
        std::sort(c.begin(),c.end());
        std::fill(mark.begin(),mark.end(),-1);
        for (int i=0;i<n;++i) mark[c[i]]=i;
        std::vector<int8_t> chi(q,-1);
        chi[0]=0;
        for (int64_t x=1;x<q;++x) chi[(x*x)%q]=1;
        std::vector<int64_t> sa(n), correction(n), sum_r(n);
        std::vector<I> norm_r(n);
        std::vector<std::vector<int64_t>> gram(n,std::vector<int64_t>(n));
        int64_t target=0;
        std::vector<int> signs(n);
        for (int x=0;x<q;++x) {
            int64_t e[7]={1,0,0,0,0,0,0};
            for (int a=0;a<n;++a) {
                int diff=x-c[a]; if (diff<0) diff+=q;
                signs[a]=chi[diff];
                for (int j=d;j>0;--j) e[j]+=signs[a]*e[j-1];
            }
            target+=e[d];
            int64_t rp=d?1:0, rm=d?1:0, rz=d?e[d-1]:0;
            for (int j=1;j<d;++j) { rp=e[j]-rp; rm=e[j]+rm; }
            if ((rp+rm)%2 || (rp-rm)%2) throw std::runtime_error("parity error");
            int64_t A=(rp+rm)/2, B=(rp-rm)/2;
            if (mark[x]>=0) correction[mark[x]]=rz-A;
            for (int a=0;a<n;++a) {
                int s=signs[a];
                int64_t r=s>0?rp:s<0?rm:rz;
                sum_r[a]+=r; norm_r[a]+=I(r)*r;
                sa[a]+=A*s;
                int64_t signed_B=B*s;
                for (int b=a;b<n;++b) gram[a][b]+=signed_B*signs[b];
            }
        }
        std::cout << "{\"q\":" << q << ",\"n\":" << n << ",\"degree\":" << d << ",\"target\":" << target << ",\"internal_contraction\":[";
        for (int a=0;a<n;++a) {
            if (a) std::cout << ',';
            std::cout << '[';
            for (int b=0;b<n;++b) {
                if (b) std::cout << ',';
                int diff=c[a]-c[b];if(diff<0) diff+=q;
                std::cout << sa[b]+gram[std::min(a,b)][std::max(a,b)]+correction[a]*chi[diff];
            }
            std::cout << ']';
        }
        std::cout << "],\"insertion_sums\":[";
        for (int a=0;a<n;++a) { if(a) std::cout << ','; std::cout << sum_r[a]; }
        std::cout << "],\"insertion_norms_squared\":[";
        for (int a=0;a<n;++a) { if(a) std::cout << ','; std::cout << decimal(norm_r[a]); }
        std::cout << "]}\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
