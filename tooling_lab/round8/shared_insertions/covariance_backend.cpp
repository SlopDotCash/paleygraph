// Shared-insertion covariance backend, extending the frozen round7 Gram method.
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
        std::vector<std::vector<int64_t>> gram_b2(n,std::vector<int64_t>(n));
        std::vector<I> sum_ab_sign(n);
        std::vector<int64_t> at_A(n), at_B(n);
        I sum_a2=0;
        // |B| <= sum_{j=d-2,d-4,...} binom(n,j), at most637393 here.
        // Hence q*B^2 <=4.063e18 in the accepted domain, within int64.
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
            sum_a2+=I(A)*A;
            if (mark[x]>=0) { correction[mark[x]]=rz-A; at_A[mark[x]]=A; at_B[mark[x]]=B; }
            for (int a=0;a<n;++a) {
                int s=signs[a];
                int64_t r=s>0?rp:s<0?rm:rz;
                sum_r[a]+=r; norm_r[a]+=I(r)*r;
                sa[a]+=A*s; sum_ab_sign[a]+=I(A)*B*s;
                int64_t signed_B=B*s;
                int64_t signed_b2=B*B*s;
                for (int b=a;b<n;++b) {
                    gram[a][b]+=signed_B*signs[b];
                    gram_b2[a][b]+=signed_b2*signs[b];
                }
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
        std::cout << "],\"derivative_gram\":[";
        for (int a=0;a<n;++a) {
            if (a) std::cout << ',';
            std::cout << '[';
            for (int b=0;b<n;++b) {
                if (b) std::cout << ',';
                int diff=c[a]-c[b];if(diff<0) diff+=q;
                int sign=chi[diff];
                I value=sum_a2+sum_ab_sign[a]+sum_ab_sign[b]+gram_b2[std::min(a,b)][std::max(a,b)]
                    +I(correction[a])*(at_A[a]+at_B[a]*sign)
                    +I(correction[b])*(at_A[b]+at_B[b]*sign);
                if(a==b) { value+=I(correction[a])*correction[a]; if(value!=norm_r[a]) throw std::runtime_error("Gram diagonal mismatch"); }
                std::cout << decimal(value);
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
