// Independent mean verifier inputs: directly count inward elementary targets.
// No forward-insertion weights, Boolean interpolation, or weighted Gram sums.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>
int64_t choose(int n,int k){if(n<0 || k<0 || k>n)return 0;int64_t x=1;for(int j=1;j<=k;++j)x=x*(n-j+1)/j;return x;}
int64_t elementary(int n,int p,int z,int d){if(d<0)return 0;int64_t v=0;for(int j=0;j<=d;++j)v+=(j%2?-1:1)*choose(n-p-z,j)*choose(p,d-j);return v;}
int main(){try{
    int q,n,d;if(!(std::cin>>q>>n>>d))throw std::runtime_error("missing parameters");
    if(q<5 || q>10000000 || q%4!=1 || n<2 || n>64 || n>q-2 || d<0 || d>6 || d>n)throw std::runtime_error("invalid parameters");
    for(int j=2;int64_t(j)*j<=q;++j)if(q%j==0)throw std::runtime_error("nonprime field");
    std::vector<int> c(n);for(int& x:c){std::cin>>x;if(x<0 || x>=q)throw std::runtime_error("invalid point");}std::sort(c.begin(),c.end());
    if(std::adjacent_find(c.begin(),c.end())!=c.end())throw std::runtime_error("duplicate point");
    int64_t table[3][3][2][65]={};
    for(int omit=0;omit<3;++omit)for(int e=0;e<3;++e)for(int z=0;z<2;++z)for(int p=0;p<=n-omit-z;++p)
        table[omit][e][z][p]=elementary(n-omit,p,z,d-2*e);
    std::vector<int8_t> chi(q,-1);chi[0]=0;for(int64_t x=1;x<=(q-1)/2;++x)chi[x*x%q]=1;
    int64_t full[3]={};std::vector<std::vector<int64_t>> singles(n,std::vector<int64_t>(3));
    std::vector<std::vector<std::vector<int64_t>>> pairs(n,std::vector<std::vector<int64_t>>(n,std::vector<int64_t>(3)));
    std::vector<int> positive(n),zero(n);
    for(int x=0;x<q;++x){int P=0,Z=0;for(int i=0;i<n;++i){int diff=x-c[i];if(diff<0)diff+=q;positive[i]=chi[diff]>0;zero[i]=chi[diff]==0;P+=positive[i];Z+=zero[i];}
        for(int e=0;e<3;++e)full[e]+=table[0][e][Z][P];
        for(int i=0;i<n;++i){int p=P-positive[i],z=Z-zero[i];
            for(int e=0;e<3;++e)singles[i][e]+=table[1][e][z][p];
            for(int j=i+1;j<n;++j){int pp=p-positive[j],zz=z-zero[j];for(int e=0;e<3;++e)pairs[i][j][e]+=table[2][e][zz][pp];}
        }
    }
    auto array=[](auto& a){std::cout<<'[';bool comma=false;for(auto x:a){if(comma)std::cout<<',';comma=true;std::cout<<x;}std::cout<<']';};
    std::cout<<"{\"q\":"<<q<<",\"degree\":"<<d<<",\"selected\":";array(c);std::cout<<",\"full\":";array(full);
    std::cout<<",\"singles\":[";for(int i=0;i<n;++i){if(i)std::cout<<',';array(singles[i]);}std::cout<<"],\"pairs\":[";bool comma=false;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){if(comma)std::cout<<',';comma=true;std::cout<<"{\"indices\":["<<i<<','<<j<<"],\"values\":";array(pairs[i][j]);std::cout<<'}';}
    std::cout<<"],\"induced\":[";for(int i=0;i<n;++i){if(i)std::cout<<',';std::cout<<'[';for(int j=0;j<n;++j){if(j)std::cout<<',';int diff=c[i]-c[j];if(diff<0)diff+=q;std::cout<<int(chi[diff]);}std::cout<<']';}std::cout<<"]}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
