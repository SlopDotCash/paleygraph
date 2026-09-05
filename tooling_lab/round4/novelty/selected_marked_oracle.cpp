// Direct elementary-coefficient oracle on all completions of selected marks.
#include <cstdint>
#include <functional>
#include <iostream>
#include <vector>
int binom(int n,int k){if(k<0||k>n)return 0;int out=1;for(int j=1;j<=k;++j)out=out*(n-j+1)/j;return out;}
int main(){
 const int q=49;int P[q][q],Q[q][q];
 for(auto&r:P)for(int&v:r)if(!(std::cin>>v))return 2;
 for(auto&r:Q)for(int&v:r)if(!(std::cin>>v))return 2;
 int requests;std::cin>>requests;
 for(int request=0;request<requests;++request){
  int n,m;std::cin>>n>>m;if(n<6||n>8||m>n)return 3;
  std::vector<int>C(m),outside;bool used[q]{};
  for(int&x:C){std::cin>>x;if(x<0||x>=q||used[x])return 4;used[x]=true;}
  for(int x=0;x<q;++x)if(!used[x])outside.push_back(x);
  int table[9][9]{};for(int a=0;a<=n;++a)for(int b=0;a+b<=n;++b)
    for(int j=0;j<=6;++j)table[a][b]+=(j%2?-1:1)*binom(b,j)*binom(a,6-j);
  int64_t count=0,sp=0,pp=0,sq=0,qq=0;
  std::function<void(int,int)>visit=[&](int pos,int left){
   if(!left){
    int vp=0,vq=0;
    for(int y=0;y<q;++y){int ap=0,bp=0,aq=0,bq=0;for(int x:C){ap+=P[y][x]==1;bp+=P[y][x]==-1;aq+=Q[y][x]==1;bq+=Q[y][x]==-1;}
      vp+=table[ap][bp];vq+=table[aq][bq];}
    ++count;sp+=vp;pp+=vp*vp;sq+=vq;qq+=vq*vq;return;
   }
   for(int j=pos;j<=int(outside.size())-left;++j){C.push_back(outside[j]);visit(j+1,left-1);C.pop_back();}
  };
  visit(0,n-m);
  std::cout<<request<<" "<<count<<" "<<sp<<" "<<pp<<" "<<sq<<" "<<qq<<"\n";
 }
}
