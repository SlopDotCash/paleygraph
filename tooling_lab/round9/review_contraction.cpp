// Separate Q verifier: coefficient counts, DIF transforms, Sg instead of Sh,
// and three-prime reconstruction of the final bilinear scalar in Python.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
using U=uint32_t;using W=uint64_t;
template<U P> U power(U x,W n){W r=1;while(n){if(n&1)r=r*x%P;x=W(x)*x%P;n/=2;}return r;}
template<U P,U R> void transform(std::vector<U>& a,bool inverse){
    size_t n=a.size();if((P-1)%n)throw std::runtime_error("bad transform length");
    for(size_t len=n;len>1;len/=2){U step=power<P>(R,(P-1)/len);if(inverse)step=power<P>(step,P-2);
        for(size_t i=0;i<n;i+=len){W w=1;for(size_t j=0;j<len/2;++j){U x=a[i+j],y=a[i+j+len/2];
            U s=x+y;if(s>=P)s-=P;a[i+j]=s;
            U delta=x>=y?x-y:x+P-y;a[i+j+len/2]=delta*w%P;w=w*step%P;
        }}
    }
    for(size_t i=0,j=0;i<n;++i){if(i<j)std::swap(a[i],a[j]);size_t b=n/2;while(b && (j&b)){j^=b;b/=2;}j^=b;}
    if(inverse){W scale=power<P>(n,P-2);for(auto& x:a)x=x*scale%P;}
}
template<U P,U R> U bilinear(const std::vector<int8_t>& chi,const std::vector<int64_t>& g,const std::vector<int64_t>& h){
    size_t q=g.size(),N=1;while(N<2*q-1)N*=2;std::vector<U> a(N),b(N);
    for(size_t i=0;i<q;++i){a[i]=chi[i]<0?P-1:chi[i];int64_t v=g[i]%P;if(v<0)v+=P;b[i]=v;}
    transform<P,R>(a,false);transform<P,R>(b,false);for(size_t i=0;i<N;++i)a[i]=W(a[i])*b[i]%P;
    b.clear();b.shrink_to_fit();transform<P,R>(a,true);W result=0;
    for(size_t i=0;i<q;++i){W v=a[i];if(i+q<2*q-1)v+=a[i+q];v%=P;int64_t w=h[i]%P;if(w<0)w+=P;result=(result+v*W(w))%P;}
    return U(result);
}
int64_t choose(int n,int k){if(k<0 || k>n)return 0;int64_t r=1;for(int j=1;j<=k;++j)r=r*(n-j+1)/j;return r;}
int64_t coeff(int positives,int negatives,int d){if(d<0)return 0;int64_t x=0;for(int j=0;j<=d;++j)x+=(j%2?-1:1)*choose(negatives,j)*choose(positives,d-j);return x;}
int main(){try{
    int q,n,d,a,b;if(!(std::cin>>q>>n>>d>>a>>b))return 1;
    if(q<5 || q>10000000 || n<2 || n>64 || d<0 || d>6)throw std::runtime_error("bad input");
    std::vector<int> base;for(int i=0,x;i<n;++i){std::cin>>x;if(x!=a && x!=b)base.push_back(x);}
    std::vector<int8_t> chi(q,-1);chi[0]=0;for(int64_t x=1;x<=q/2;++x)chi[x*x%q]=1;
    int64_t table_g[2][65],table_h[2][65];int s=base.size();
    for(int z=0;z<=1;++z)for(int p=0;p<=s;++p){table_g[z][p]=coeff(p,s-z-p,d-1);table_h[z][p]=coeff(p,s-z-p,d-2);}
    std::vector<int64_t> g(q),h(q);int64_t maxg=0,maxh=0;
    for(int x=0;x<q;++x){int positive=0,zero=0;for(int j:base){int v=x-j;if(v<0)v+=q;positive+=chi[v]==1;zero+=chi[v]==0;}
        g[x]=table_g[zero][positive];h[x]=table_h[zero][positive];maxg=std::max(maxg,std::abs(g[x]));maxh=std::max(maxh,std::abs(h[x]));
    }
    U r1=bilinear<2013265921,31>(chi,g,h),r2=bilinear<1811939329,13>(chi,g,h),r3=bilinear<469762049,3>(chi,g,h);
    std::cout<<"{\"residues\":["<<r1<<','<<r2<<','<<r3<<"],\"primes\":[2013265921,1811939329,469762049],\"max_abs_g\":"<<maxg<<",\"max_abs_h\":"<<maxh<<"}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
