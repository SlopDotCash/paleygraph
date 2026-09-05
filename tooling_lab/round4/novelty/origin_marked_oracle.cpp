// Independent T6 sums for selected marks in Paley49/Peisert49.
// Ordinary row products; origin-containing subsets only. No cell-pair compiler.
#include <cstdint>
#include <iostream>
#include <vector>
int main(){
 const int q=49; int P[q][q],Q[q][q];
 for(auto &r:P)for(int &x:r)if(!(std::cin>>x))return 2;
 for(auto &r:Q)for(int &x:r)if(!(std::cin>>x))return 2;
 int h;std::cin>>h;std::vector<uint64_t> mark(h);for(auto &v:mark)std::cin>>v;
 std::vector<int64_t> count(h),sp(h),sq(h),pp(h),qq(h);
 for(int a=1;a<q-4;++a)for(int b=a+1;b<q-3;++b)
 for(int c=b+1;c<q-2;++c)for(int d=c+1;d<q-1;++d){
  int p[q],r[q];for(int y=0;y<q;++y){p[y]=P[y][0]*P[y][a]*P[y][b]*P[y][c]*P[y][d];r[y]=Q[y][0]*Q[y][a]*Q[y][b]*Q[y][c]*Q[y][d];}
  uint64_t partial=1|(uint64_t(1)<<a)|(uint64_t(1)<<b)|(uint64_t(1)<<c)|(uint64_t(1)<<d);
  for(int e=d+1;e<q;++e){
   int vp=0,vq=0;for(int y=0;y<q;++y){vp+=p[y]*P[y][e];vq+=r[y]*Q[y][e];}
   uint64_t chosen=partial|(uint64_t(1)<<e);
   for(int i=0;i<h;++i)if((chosen&mark[i])==mark[i]){++count[i];sp[i]+=vp;sq[i]+=vq;pp[i]+=vp*vp;qq[i]+=vq*vq;}
  }
 }
 for(int i=0;i<h;++i)std::cout<<i<<" "<<count[i]<<" "<<sp[i]<<" "<<pp[i]<<" "<<sq[i]<<" "<<qq[i]<<"\n";
}
