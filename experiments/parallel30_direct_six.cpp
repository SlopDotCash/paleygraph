// Independent complete enumeration with the first coordinate normalized to one.
// Pass 30 adds direct detection of proper zero-sum triple partitions.
// Uses no triangle-orbit lower bound to compute the count.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <unordered_map>
#include <vector>

using U = std::uint64_t;
U mul(U a, U b, U p) { return static_cast<U>((static_cast<__uint128_t>(a)*b)%p); }
U power(U a, U b, U p) {
  U out=1;
  while(b) { if(b&1) out=mul(out,a,p); a=mul(a,a,p); b>>=1; }
  return out;
}
struct Pair { U a,b,weight; };

int main() {
  U p,n;
  while(std::cin >> p >> n) {
    assert(n>=4 && n<=1024 && (n&(n-1))==0 && (p-1)%n==0 && p<(U(1)<<50));
    U g=0;
    for(U a=3; a<p; ++a) {
      U candidate=power(a,(p-1)/n,p);
      if(power(candidate,n/2,p)!=1) { g=candidate; break; }
    }
    assert(g && power(g,n,p)==1);
    std::vector<U> H(n);
    H[0]=1;
    for(U i=1;i<n;++i) H[i]=mul(H[i-1],g,p);
    std::sort(H.begin(),H.end());
    assert(std::adjacent_find(H.begin(),H.end())==H.end());
    std::unordered_map<U,std::vector<Pair>> pairs;
    pairs.reserve(n*n);
    for(U i=0;i<n;++i) for(U j=i;j<n;++j)
      pairs[(1+H[i]+H[j])%p].push_back({H[i],H[j],i==j?U(1):U(2)});
    U E3=0,T6=0,J6=0,repeated=0,balanced=0,D6=0,triangleD6=0,candidates=0;
    for(U i=0;i<n;++i) for(U j=i;j<n;++j) for(U k=j;k<n;++k) {
      U sum=(H[i]+H[j]+H[k])%p;
      auto it=pairs.find(sum? p-sum:0);
      if(it==pairs.end()) continue;
      U wt=i==k?1:(i==j || j==k?3:6);
      for(const Pair& pair:it->second) {
        ++candidates;
        U mass=n*wt*pair.weight;
        std::array<U,6> word={1,pair.a,pair.b,H[i],H[j],H[k]};
        bool intrinsic=true,has_opposite=false,has_repeat=false;
        for(int a=0;a<6;++a) {
          int equal_count=0,negative_count=0;
          for(int b=0;b<6;++b) {
            equal_count+=(word[a]==word[b]);
            negative_count+=(word[a]+word[b]==p);
            if(a<b && word[a]==word[b]) has_repeat=true;
            if(word[a]+word[b]==p) has_opposite=true;
          }
          intrinsic &= equal_count==negative_count;
        }
        E3+=mass;
        if(intrinsic) { T6+=mass; continue; }
        if(has_opposite) { J6+=mass; continue; }
        if(has_repeat) { repeated+=mass; continue; }
        bool any_balanced=false;
        for(int a=1;a<6;++a) for(int b=a+1;b<6;++b) {
          U left=1,right=1;
          for(int j2=0;j2<6;++j2) {
            if(j2==0 || j2==a || j2==b) left=mul(left,word[j2],p);
            else right=mul(right,word[j2],p);
          }
          any_balanced |= (left+right==p);
        }
        if(any_balanced) balanced+=mass;
        else {
          D6+=mass;
          int zeroPartitions=0;
          for(int a=1;a<6;++a) for(int b=a+1;b<6;++b)
            zeroPartitions += (word[0]+word[a]+word[b])%p==0;
          assert(zeroPartitions<=1);
          if(zeroPartitions) triangleD6+=mass;
        }
      }
    }
    std::cout << "{\"p\":" << p << ",\"n\":" << n
      << ",\"E3\":" << E3 << ",\"T6\":" << T6 << ",\"J6\":" << J6
      << ",\"repeated_R6\":" << repeated << ",\"distinct_balanced_R6\":" << balanced
      << ",\"distinct_unbalanced_R6\":" << D6
      << ",\"triangle_D6\":" << triangleD6
      << ",\"primitive_D6\":" << D6-triangleD6
      << ",\"normalized_multiset_candidates\":" << candidates << "}" << std::endl;
  }
}
