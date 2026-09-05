// Standard exact NTT, adapted from the frozen round4 marked_counts backend.
#pragma once
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
namespace exact {
constexpr int64_t modulus=998244353,root=3;
constexpr size_t max_length=1u<<23;
inline void require(bool ok,const std::string&why){if(!ok)throw std::runtime_error(why);}
inline bool prime(uint64_t n){if(n<2)return false;for(uint64_t d=2;d*d<=n;d++)if(n%d==0)return false;return true;}
inline int64_t power(int64_t a,int64_t n){int64_t r=1;while(n){if(n&1)r=r*a%modulus;a=a*a%modulus;n>>=1;}return r;}
inline void validate(){
 require(prime(modulus) && modulus-1==(int64_t(1)<<23)*7*17,"invalid NTT prime or factorization");
 for(int d:{2,7,17})require(power(root,(modulus-1)/d)!=1,"invalid primitive root");
}
inline void ntt(std::vector<int64_t>&a,bool invert){
 size_t n=a.size();require(n && !(n&(n-1)) && n<=max_length && (modulus-1)%n==0,"unsupported transform length");
 for(size_t i=1,j=0;i<n;i++){size_t bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j)std::swap(a[i],a[j]);}
 for(size_t len=2;len<=n;len<<=1){int64_t wlen=power(root,(modulus-1)/len);if(invert)wlen=power(wlen,modulus-2);
  for(size_t i=0;i<n;i+=len){int64_t w=1;for(size_t j=0;j<len/2;j++){
   int64_t u=a[i+j],v=a[i+j+len/2]*w%modulus;
   a[i+j]=u+v;if(a[i+j]>=modulus)a[i+j]-=modulus;
   a[i+j+len/2]=u-v;if(a[i+j+len/2]<0)a[i+j+len/2]+=modulus;w=w*wlen%modulus;
  }}
 }
 if(invert){int64_t ni=power(n,modulus-2);for(auto&x:a)x=x*ni%modulus;}
}
}
