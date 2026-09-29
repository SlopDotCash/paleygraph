// Exact selected contractions and one character convolution for two insertions.
// Iterative radix-two NTT follows the earlier lab backend; two new primes extend
// its length and signed reconstruction range. No floating-point convolution.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I=__int128_t;
using U=uint32_t;
using W=uint64_t;

void require(bool b,const std::string& s){if(!b)throw std::runtime_error(s);}
bool prime(W n){if(n<2)return false;for(W i=2;i*i<=n;++i)if(n%i==0)return false;return true;}
std::string dec(I x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;std::string s;while(x){s.push_back('0'+x%10);x/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
U power(U a,W n,U mod){W r=1;while(n){if(n&1)r=r*a%mod;a=W(a)*a%mod;n>>=1;}return U(r);}
template<U MOD,U ROOT> void ntt(std::vector<U>& a,bool invert){
    size_t n=a.size();require(n && !(n&(n-1)) && (MOD-1)%n==0,"invalid NTT length");
    for(size_t i=1,j=0;i<n;++i){size_t bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j)std::swap(a[i],a[j]);}
    for(size_t len=2;len<=n;len<<=1){U step=power(ROOT,(MOD-1)/len,MOD);if(invert)step=power(step,MOD-2,MOD);
        for(size_t i=0;i<n;i+=len){W w=1;for(size_t j=0;j<len/2;++j){
            U u=a[i+j],v=W(a[i+j+len/2])*w%MOD;
            U sum=u+v; if(sum>=MOD)sum-=MOD;
            a[i+j]=sum;a[i+j+len/2]=u>=v?u-v:u+MOD-v;w=w*step%MOD;
        }}
    }
    if(invert){W ni=power(n,MOD-2,MOD);for(auto& x:a)x=x*ni%MOD;}
}
template<U MOD,U ROOT> std::vector<U> convolution(const std::vector<int8_t>& chi,const std::vector<int64_t>& h){
    size_t q=h.size(),length=1;while(length<2*q-1)length<<=1;
    std::vector<U> a(length),b(length);
    for(size_t x=0;x<q;++x){a[x]=chi[x]<0?MOD-1:chi[x];int64_t r=h[x]%MOD;if(r<0)r+=MOD;b[x]=U(r);}
    ntt<MOD,ROOT>(a,false);ntt<MOD,ROOT>(b,false);
    for(size_t x=0;x<length;++x)a[x]=W(a[x])*b[x]%MOD;
    b.clear();b.shrink_to_fit();ntt<MOD,ROOT>(a,true);
    std::vector<U> out(q);
    for(size_t x=0;x<q;++x){W v=a[x];if(x+q<2*q-1)v+=a[x+q];out[x]=v%MOD;}
    return out;
}
template<class T> void array(const std::vector<T>& a){std::cout<<'[';for(size_t i=0;i<a.size();++i){if(i)std::cout<<',';std::cout<<dec(a[i]);}std::cout<<']';}

int main(){try{
    int q,n,d,a0,a1;
    require(bool(std::cin>>q>>n>>d>>a0>>a1),"missing parameters");
    require(q>=5 && q<=10000000 && q%4==1 && prime(q) && n>=2 && n<=64 && n<=q-2 && d>=0 && d<=6 && d<=n,"invalid parameters");
    std::vector<int> c(n),mark(q,-1);
    for(int& x:c){require(bool(std::cin>>x),"missing set entry");require(x>=0 && x<q && mark[x]<0,"invalid set entry");mark[x]=0;}
    require(a0!=a1 && a0>=0 && a0<q && a1>=0 && a1<q && mark[a0]>=0 && mark[a1]>=0,"invalid deletion pair");
    std::sort(c.begin(),c.end());if(a0>a1)std::swap(a0,a1);
    for(int i=0;i<n;++i)mark[c[i]]=i;
    std::vector<int8_t> chi(q,-1);chi[0]=0;for(int64_t x=1;x<q;++x)chi[x*x%q]=1;
    std::vector<int64_t> g(q),h(q),gs(n),hs(n),f(n),v(n);
    std::vector<I> w(n);
    std::vector<std::vector<int64_t>> K(n,std::vector<int64_t>(n));
    I G0=0,G2=0,H0=0,H2=0,t0=0;int64_t maxh=0;bool zero_g=true,constant_h=true;
    std::vector<int> signs(n);
    auto begin=std::chrono::steady_clock::now();
    for(int x=0;x<q;++x){
        int64_t e[7]={1,0,0,0,0,0,0};
        for(int i=0;i<n;++i){int diff=x-c[i];if(diff<0)diff+=q;signs[i]=chi[diff];
            if(c[i]!=a0 && c[i]!=a1)for(int j=d;j>0;--j)e[j]+=signs[i]*e[j-1];
        }
        g[x]=d>=1?e[d-1]:0;h[x]=d>=2?e[d-2]:0;t0+=e[d];
        G0+=g[x];G2+=I(g[x])*g[x];H0+=h[x];H2+=I(h[x])*h[x];
        maxh=std::max(maxh,std::abs(h[x]));zero_g &= g[x]==0;if(x)constant_h &= h[x]==h[0];
        if(mark[x]>=0){gs[mark[x]]=g[x];hs[mark[x]]=h[x];}
        for(int i=0;i<n;++i){f[i]+=g[x]*signs[i];v[i]+=h[x]*signs[i];w[i]+=I(g[x])*h[x]*signs[i];
            int64_t sh=h[x]*signs[i];for(int j=i;j<n;++j)K[i][j]+=sh*signs[j];
        }
    }
    auto streamed=std::chrono::steady_clock::now();
    constexpr U P1=2013265921,R1=31,P2=1811939329,R2=13;
    require(prime(P1)&&prime(P2),"invalid NTT prime");
    for(U factor:{2,3,5})require(power(R1,(P1-1)/factor,P1)!=1,"invalid first primitive root");
    for(U factor:{2,3})require(power(R2,(P2-1)/factor,P2)!=1,"invalid second primitive root");
    require(P1-1==W(15)*(1u<<27) && P2-1==W(27)*(1u<<26),"invalid prime factorization");
    int64_t modulus=int64_t(P1)*P2,bound=int64_t(q-1)*maxh;
    require(modulus>2*bound,"CRT range insufficient");
    I Q=0,transform_sum=0,transform_norm=0;int checked=0;
    bool skip=zero_g||constant_h;
    if(!skip){
        auto r1=convolution<P1,R1>(chi,h);auto r2=convolution<P2,R2>(chi,h);
        U inv=power(P1%P2,P2-2,P2);
        for(int x=0;x<q;++x){
            int64_t diff=int64_t(r2[x])-int64_t(r1[x]%P2);if(diff<0)diff+=P2;
            W lift=W(diff)*inv%P2;int64_t y=int64_t(r1[x])+int64_t(P1)*int64_t(lift);if(y>modulus/2)y-=modulus;
            require(std::abs(y)<=bound,"CRT coordinate exceeds proven bound");
            if(mark[x]>=0){require(y==v[mark[x]],"selected direct convolution mismatch");++checked;}
            Q+=I(g[x])*y;transform_sum+=y;transform_norm+=I(y)*y;
        }
        require(transform_sum==0 && transform_norm==I(q)*H2-H0*H0,"full transform sum or norm mismatch");
    }
    auto finish=std::chrono::steady_clock::now();
    std::cout<<"{\"q\":"<<q<<",\"selected\":";array(c);std::cout<<",\"deleted\":["<<a0<<','<<a1<<"],\"degree\":"<<d;
    std::cout<<",\"t0\":"<<dec(t0)<<",\"G0\":"<<dec(G0)<<",\"G2\":"<<dec(G2)<<",\"H0\":"<<dec(H0)<<",\"H2\":"<<dec(H2)<<",\"Q\":"<<dec(Q);
    std::cout<<",\"g_selected\":";array(gs);std::cout<<",\"h_selected\":";array(hs);
    std::cout<<",\"f_selected\":";array(f);std::cout<<",\"v_selected\":";array(v);std::cout<<",\"w_selected\":";array(w);
    std::cout<<",\"K_selected\":[";for(int i=0;i<n;++i){if(i)std::cout<<',';std::cout<<'[';for(int j=0;j<n;++j){if(j)std::cout<<',';std::cout<<K[std::min(i,j)][std::max(i,j)];}std::cout<<']';}std::cout<<']';
    std::cout<<",\"convolution\":{\"skipped_exact_zero_Q\":"<<(skip?"true":"false")<<",\"primes\":["<<P1<<','<<P2<<"],\"roots\":["<<R1<<','<<R2<<"],\"coordinate_bound\":"<<bound<<",\"crt_modulus\":"<<modulus
             <<",\"selected_coordinates_compared\":"<<checked<<",\"sum\":"<<dec(transform_sum)<<",\"norm_squared\":"<<dec(transform_norm)<<"}";
    std::cout<<",\"timing\":{\"stream_seconds\":"<<std::chrono::duration<double>(streamed-begin).count()<<",\"convolution_seconds\":"<<std::chrono::duration<double>(finish-streamed).count()<<"}}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
