// Direct elementary-sign coefficient sums on all origin-containing subsets.
// No character-trace identities, row-triple compiler, or moment interpolation.
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <stdexcept>

constexpr int q=49;
uint64_t choose(int n,int k){uint64_t a=1;for(int i=1;i<=k;++i)a=a*(n-k+i)/i;return a;}
int main(int argc,char**argv){
 const int n=argc>1?std::stoi(argv[1]):7;
 if(n<6||n>8)throw std::runtime_error("only n=6,7,8 supported by this audited integer budget");
 const auto started=std::chrono::steady_clock::now();
 std::array<std::array<uint64_t,q>,2> positive{};
 for(int g=0;g<2;++g)for(int x=0;x<q;++x)for(int y=0;y<q;++y){
  int sign;if(!(std::cin>>sign))return 2;
  if(sign==1)positive[g][x]|=uint64_t(1)<<y;
  if(sign < -1||sign>1||(x==y)!=(sign==0))return 3;
 }
 int coefficients[2][9]{};
 for(int inside=0;inside<=1;++inside)for(int a=0;a<=n-inside;++a){
  int b=n-inside-a,v=0;
  for(int j=0;j<=6;++j)if(j<=b&&6-j<=a)
   v+=(j%2?-1:1)*int(choose(b,j)*choose(a,6-j));
  coefficients[inside][a]=v;
 }
 // Even the global worst-case sum of absolute cubes is below 2^63 for n<=8:
 // C(49,8)*(49*C(8,6))^3 = 1164709865022979968.
 std::array<std::array<int64_t,3>,2> moments{};
 std::array<int,2> minimum{1000000,1000000},maximum{-1000000,-1000000};
 uint64_t count=0,selection=(uint64_t(1)<<(n-1))-1;
 const uint64_t limit=uint64_t(1)<<48;
 while(selection<limit){
  const uint64_t mask=(selection<<1)|1;
  for(int g=0;g<2;++g){
   int value=0;
   for(int x=0;x<q;++x){
    int a=__builtin_popcountll(positive[g][x]&mask);
    value+=coefficients[(mask>>x)&1][a];
   }
   const int64_t v=value;
   moments[g][0]+=v;moments[g][1]+=v*v;moments[g][2]+=v*v*v;
   if(value<minimum[g])minimum[g]=value;if(value>maximum[g])maximum[g]=value;
  }
  ++count;
  const uint64_t low=selection&-selection,next=selection+low;
  selection=next|(((selection^next)>>2)/low);
 }
 if(count!=choose(48,n-1))return 4;
 std::cout<<"COUNT "<<count<<"\n";
 for(int g=0;g<2;++g)std::cout<<g<<" "<<moments[g][0]<<" "<<moments[g][1]<<" "<<moments[g][2]<<" "<<minimum[g]<<" "<<maximum[g]<<"\n";
 std::cout<<"SECONDS "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<"\n";
}
