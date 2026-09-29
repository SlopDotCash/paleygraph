// Independent exact pointwise oracle. Enumerates all completions; no cell counts.
#include <cstdint>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace std;
void need(bool x){if(!x)throw runtime_error("oracle input/check failed");}
int main(){try{
 int graphs;cin>>graphs;vector<vector<vector<int>>>S(graphs);vector<vector<uint64_t>>negative(graphs);
 for(int g=0;g<graphs;g++){
  int q;cin>>q;need(q<=63 && q>=5);S[g]=vector<vector<int>>(q,vector<int>(q));negative[g].resize(q);
  for(int x=0;x<q;x++)for(int y=0;y<q;y++){int s;cin>>s;need(s>=-1 && s<=1);S[g][x][y]=s;if(s==-1)negative[g][x]|=uint64_t(1)<<y;}
  for(int x=0;x<q;x++){
   int rowsum=0;for(int y=0;y<q;y++){need(S[g][x][y]==S[g][y][x]);need((S[g][x][y]==0)==(x==y));rowsum+=S[g][x][y];}need(rowsum==0);
   for(int y=0;y<q;y++){int dot=0;for(int z=0;z<q;z++)dot+=S[g][x][z]*S[g][y][z];need(dot==(x==y?q-1:-1));}
  }
 }
 int queries;cin>>queries;
 for(int test=0;test<queries;test++){
  int g,n,m;cin>>g>>n>>m;int q=S[g].size();need(6<=n && n<=8 && 0<=m && m<=n);uint64_t fixed=0;
  for(int k=0;k<m;k++){int v;cin>>v;need(0<=v && v<q && !(fixed>>v&1));fixed|=uint64_t(1)<<v;}
  vector<int>outside;for(int x=0;x<q;x++)if(!(fixed>>x&1))outside.push_back(x);
  // e6 of a positives and b negatives, computed by direct polynomial update.
  int table[9][9]={};for(int a=0;a<=n;a++)for(int b=0;a+b<=n;b++){
   int coef[7]={1,0,0,0,0,0,0};for(int k=0;k<a+b;k++){int sign=k<a?1:-1;for(int j=6;j>=1;j--)coef[j]+=sign*coef[j-1];}table[a][b]=coef[6];
  }
  int64_t count=0,total=0,squares=0;map<int,int64_t>hist;
  auto recurse=[&](auto&&self,int at,int left,uint64_t mask)->void{
   if(left==0){
    int value=0;for(int x=0;x<q;x++){int b=__builtin_popcountll(mask&negative[g][x]);int a=n-b-int(mask>>x&1);value+=table[a][b];}
    count++;total+=value;squares+=int64_t(value)*value;hist[value]++;return;
   }
   for(int i=at;i+left<=int(outside.size());i++)self(self,i+1,left-1,mask|(uint64_t(1)<<outside[i]));
  };
  recurse(recurse,0,n-m,fixed);
  cout<<"{\"query\":"<<test<<",\"graph_index\":"<<g<<",\"n\":"<<n<<",\"marks_mask\":"<<fixed<<",\"completions\":"<<count<<",\"sum_T6\":"<<total<<",\"sum_T6_squared\":"<<squares<<",\"histogram\":[";
  bool first=true;for(auto item:hist){if(!first)cout<<',';first=false;cout<<'['<<item.first<<','<<item.second<<']';}cout<<"]}\n";
 }
 return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
