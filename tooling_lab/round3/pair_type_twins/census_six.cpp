#include <array>
#include <cstdint>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

// Exact bit-parity census. Inputs are negative-entry masks of two symmetric
// sign matrices with zero diagonal. No floating point or graph library.
int main() {
  int q; std::cin >> q;
  if(q<6 || q>63) return 2;
  std::vector<uint64_t> P(q), Q(q);
  for(auto &x:P) std::cin>>x;
  for(auto &x:Q) std::cin>>x;
  std::map<int,uint64_t> hp,hq;
  std::map<std::pair<int,int>,uint64_t> joint;
  std::map<int,std::array<int,6>> wp,wq;
  uint64_t total=0;
  int biggest=-1; std::array<int,6> difference{}; int dp=0,dq=0;
  const uint64_t all=(uint64_t(1)<<q)-1;
  for(int a=0;a<q-5;a++) for(int b=a+1;b<q-4;b++)
  for(int c=b+1;c<q-3;c++) for(int d=c+1;d<q-2;d++)
  for(int e=d+1;e<q-1;e++) {
    uint64_t pp=P[a]^P[b]^P[c]^P[d]^P[e];
    uint64_t qq=Q[a]^Q[b]^Q[c]^Q[d]^Q[e];
    uint64_t marked=(uint64_t(1)<<a)|(uint64_t(1)<<b)|(uint64_t(1)<<c)|
                    (uint64_t(1)<<d)|(uint64_t(1)<<e);
    for(int f=e+1;f<q;f++) {
      uint64_t live=all&~(marked|(uint64_t(1)<<f));
      int vp=q-6-2*__builtin_popcountll((pp^P[f])&live);
      int vq=q-6-2*__builtin_popcountll((qq^Q[f])&live);
      std::array<int,6> set={a,b,c,d,e,f};
      ++hp[vp]; ++hq[vq]; ++joint[{vp,vq}]; ++total;
      if(!wp.count(vp))wp[vp]=set;
      if(!wq.count(vq))wq[vq]=set;
      int delta=vp>vq?vp-vq:vq-vp;
      if(delta>biggest){biggest=delta;difference=set;dp=vp;dq=vq;}
    }
  }
  std::cout<<"COUNT "<<total<<"\n";
  for(auto [value,count]:hp){std::cout<<"P "<<value<<" "<<count;for(int i:wp[value])std::cout<<" "<<i;std::cout<<"\n";}
  for(auto [value,count]:hq){std::cout<<"Q "<<value<<" "<<count;for(int i:wq[value])std::cout<<" "<<i;std::cout<<"\n";}
  for(auto [values,count]:joint)std::cout<<"J "<<values.first<<" "<<values.second<<" "<<count<<"\n";
  std::cout<<"D "<<dp<<" "<<dq;for(int i:difference)std::cout<<" "<<i;std::cout<<"\n";
}
