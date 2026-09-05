// Independent histogram oracle: ordinary integer products, six-sets containing 0.
// Translation invariance permits exact conversion to the full subset census.
#include <array>
#include <cstdint>
#include <iostream>

int main() {
  constexpr int q=49;
  int P[q][q], Q[q][q];
  for(auto &row:P) for(int &v:row) if(!(std::cin>>v)) return 2;
  for(auto &row:Q) for(int &v:row) if(!(std::cin>>v)) return 2;
  std::array<uint64_t,99> hp{},hq{};
  std::array<std::array<uint64_t,99>,99> joint{};
  uint64_t count=0;
  for(int a=1;a<q-4;++a) for(int b=a+1;b<q-3;++b)
  for(int c=b+1;c<q-2;++c) for(int d=c+1;d<q-1;++d) {
    int partialP[q],partialQ[q];
    for(int y=0;y<q;++y) {
      partialP[y]=P[y][0]*P[y][a]*P[y][b]*P[y][c]*P[y][d];
      partialQ[y]=Q[y][0]*Q[y][a]*Q[y][b]*Q[y][c]*Q[y][d];
    }
    for(int e=d+1;e<q;++e) {
      int vp=0,vq=0;
      for(int y=0;y<q;++y) {vp+=partialP[y]*P[y][e];vq+=partialQ[y]*Q[y][e];}
      ++hp[vp+49];++hq[vq+49];++joint[vp+49][vq+49];++count;
    }
  }
  std::cout<<"COUNT "<<count<<"\n";
  for(int i=0;i<99;++i) {
    if(hp[i])std::cout<<"P "<<i-49<<" "<<hp[i]<<"\n";
    if(hq[i])std::cout<<"Q "<<i-49<<" "<<hq[i]<<"\n";
    for(int j=0;j<99;++j) if(joint[i][j])std::cout<<"J "<<i-49<<" "<<j-49<<" "<<joint[i][j]<<"\n";
  }
}
