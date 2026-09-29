// All fixed-deletion-pair two-insertion means. One row pass; no convolution.
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using I=__int128_t;
std::string dec(I x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;std::string s;while(x){s.push_back('0'+x%10);x/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
int main(){try{
    int q,n,d;require(bool(std::cin>>q>>n>>d),"missing parameters");
    require(q>=5 && q<=10000000 && q%4==1 && n>=2 && n<=64 && n<=q-2 && d>=0 && d<=6 && d<=n,"invalid parameters");
    for(int j=2;int64_t(j)*j<=q;++j)require(q%j!=0,"q is not prime");
    std::vector<int> c(n),mark(q,-1);
    for(int& a:c){require(bool(std::cin>>a),"missing set entry");require(a>=0 && a<q && mark[a]<0,"invalid selected set");mark[a]=0;}
    std::sort(c.begin(),c.end());for(int i=0;i<n;++i)mark[c[i]]=i;
    std::vector<int8_t> chi(q,-1);chi[0]=0;for(int64_t x=1;x<q;++x)chi[x*x%q]=1;
    int64_t m=q-n;I N=I(m)*(m-1),sum0=0;
    std::vector<I> linear(n);std::vector<std::vector<I>> gram(n,std::vector<I>(n)),correction(n,std::vector<I>(n));
    std::vector<int> signs(n);
    for(int x=0;x<q;++x){
        int64_t e[7]={1,0,0,0,0,0,0};int R=0,inside=mark[x]>=0;
        for(int i=0;i<n;++i){int diff=x-c[i];if(diff<0)diff+=q;int sign=chi[diff];signs[i]=sign;R+=sign;
            for(int j=d;j>0;--j)e[j]+=sign*e[j-1];
        }
        auto phi=[&](int u,int v){int64_t e1[7]={1},e2[7]={1};
            for(int j=1;j<=d;++j){e1[j]=e[j]-u*e1[j-1];e2[j]=e1[j]-v*e2[j-1];}
            auto get=[&](int j)->I{return j>=0?e2[j]:0;};
            return N*get(d)-2*I(m-1)*R*get(d-1)+(I(R)*R-(m-1)-inside)*get(d-2);
        };
        I pp=phi(1,1),mm=phi(-1,-1),pm=phi(1,-1);
        I w0=pp+mm+2*pm,w1=pp-mm,w2=pp+mm-2*pm;sum0+=w0;
        for(int i=0;i<n;++i){linear[i]+=w1*signs[i];I weight=w2*signs[i];
            for(int j=i+1;j<n;++j)gram[i][j]+=weight*signs[j];
        }
        if(inside){int i=mark[x];for(int j=0;j<n;++j)if(i!=j){int s=signs[j];correction[std::min(i,j)][std::max(i,j)]+=4*phi(0,s)-w0-w1*s;}}
    }
    std::cout<<"{\"q\":"<<q<<",\"degree\":"<<d<<",\"selected\":[";
    for(int i=0;i<n;++i){if(i)std::cout<<',';std::cout<<c[i];}
    std::cout<<"],\"ordered_insertion_pairs_per_deletion\":"<<dec(N)<<",\"pairs\":[";bool comma=false;
    for(int i=0;i<n;++i)for(int j=i+1;j<n;++j){I four=sum0+linear[i]+linear[j]+gram[i][j]+correction[i][j];require(four%4==0,"nonintegral numerator");
        if(comma)std::cout<<',';comma=true;std::cout<<"{\"deleted\":["<<c[i]<<','<<c[j]<<"],\"target_sum\":"<<dec(four/4)<<'}';
    }
    std::cout<<"]}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
