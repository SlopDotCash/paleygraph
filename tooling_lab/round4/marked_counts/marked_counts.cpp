// Exact marked-cell relation counts for prime Paley sign matrices.
// Standard radix-2 number-theoretic transforms, adapted from the round2 backend.
// No novel transform algorithm or spectral bound is claimed.
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
constexpr int64_t MOD=998244353,ROOT=3;
constexpr size_t MAX_LENGTH=1u<<23;
void require(bool ok,const string &why){if(!ok)throw runtime_error(why);}
bool is_prime(uint64_t n){if(n<2)return false;for(uint64_t d=2;d*d<=n;d++)if(n%d==0)return false;return true;}
int64_t power(int64_t a,int64_t n){int64_t r=1;while(n){if(n&1)r=r*a%MOD;a=a*a%MOD;n>>=1;}return r;}
void ntt(vector<int64_t>&a,bool invert){
 size_t n=a.size();require(n && (n&(n-1))==0 && n<=MAX_LENGTH && (MOD-1)%n==0,"unsupported transform length");
 for(size_t i=1,j=0;i<n;i++){size_t bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j)swap(a[i],a[j]);}
 for(size_t len=2;len<=n;len<<=1){int64_t wlen=power(ROOT,(MOD-1)/len);if(invert)wlen=power(wlen,MOD-2);
  for(size_t i=0;i<n;i+=len){int64_t w=1;for(size_t j=0;j<len/2;j++){
   int64_t u=a[i+j],v=a[i+j+len/2]*w%MOD;
   a[i+j]=u+v;if(a[i+j]>=MOD)a[i+j]-=MOD;
   a[i+j+len/2]=u-v;if(a[i+j+len/2]<0)a[i+j+len/2]+=MOD;
   w=w*wlen%MOD;
  }}
 }
 if(invert){int64_t ni=power(n,MOD-2);for(auto &x:a)x=x*ni%MOD;}
}
void pattern_json(const vector<int>&pattern){cout<<'[';for(size_t i=0;i<pattern.size();i++){if(i)cout<<',';cout<<pattern[i];}cout<<']';}
int main(int argc,char**argv){try{
 auto start=chrono::steady_clock::now();require(argc>=2,"usage: marked_counts p [distinct marks]");
 size_t used=0;uint64_t largep=stoull(argv[1],&used);require(used==string(argv[1]).size(),"invalid p");
 require(largep>=5 && largep%4==1 && largep<=(MAX_LENGTH+1)/2,"p must be 1 mod 4 and fit exact transform domain");
 require(is_prime(largep),"field parameter is not prime");uint32_t p=largep;
 require(argc-2<=6,"at most six distinct marks supported");vector<int>marks;
 for(int i=2;i<argc;i++){int64_t m=stoll(argv[i],&used);require(used==string(argv[i]).size() && m>=0 && m<p,"marks must be canonical field elements");marks.push_back(m);}
 require(set<int>(marks.begin(),marks.end()).size()==marks.size(),"marks must be distinct");
 require(is_prime(MOD) && (MOD-1)==(int64_t(1)<<23)*7*17,"invalid NTT modulus or factorization");
 for(int d:{2,7,17})require(power(ROOT,(MOD-1)/d)!=1,"NTT root is not primitive");
 size_t L=1;while(L<2*size_t(p)-1)L<<=1;require(L<=MAX_LENGTH,"linear convolution exceeds transform domain");
 vector<int8_t>chi(p,-1);chi[0]=0;for(uint64_t x=1;x<=(p-1)/2;x++)chi[x*x%p]=1;
 require(count(chi.begin(),chi.end(),int8_t(1))==(p-1)/2,"incorrect quadratic character");
 map<vector<int>,vector<uint32_t>>groups;
 for(uint32_t x=0;x<p;x++){vector<int>pattern;for(int m:marks)pattern.push_back(chi[(x+p-m)%p]);groups[pattern].push_back(x);}
 vector<vector<int>>patterns;vector<vector<uint32_t>>cells;
 for(auto &kv:groups){patterns.push_back(kv.first);cells.push_back(std::move(kv.second));}
 size_t c=cells.size();vector<size_t>cell_of(p);for(size_t i=0;i<c;i++)for(auto x:cells[i])cell_of[x]=i;
 vector<uint32_t>samples{0,1,2,p-1};for(int m:marks)samples.push_back(m);
 uint64_t state=20260905;
 while(samples.size()<36){state=(1664525*state+1013904223)%4294967296ULL;samples.push_back(state%p);}
 sort(samples.begin(),samples.end());samples.erase(unique(samples.begin(),samples.end()),samples.end());
 vector<vector<int64_t>>sample_values(samples.size(),vector<int64_t>(c));
 vector<vector<int64_t>>signed_sums(c,vector<int64_t>(c));
 vector<int64_t>kernel(L),work(L);for(uint32_t x=0;x<p;x++)kernel[x]=(chi[x]+MOD)%MOD;
 ntt(kernel,false);size_t transforms=1,singletons=0,fullcells=0;int64_t minmargin=MOD;
 for(size_t v=0;v<c;v++){
  const auto&V=cells[v];require(2*int64_t(V.size())<MOD,"ambiguous centered integer recovery");
  minmargin=min(minmargin,MOD-2*int64_t(V.size()));bool singleton=V.size()==1,full=V.size()==p;
  if(singleton)singletons++;else if(full)fullcells++;else{
   fill(work.begin(),work.end(),0);for(auto y:V)work[y]=1;
   ntt(work,false);for(size_t i=0;i<L;i++)work[i]=work[i]*kernel[i]%MOD;ntt(work,true);transforms+=2;
  }
  size_t sample_idx=0;int64_t global_sum=0;
  for(uint32_t x=0;x<p;x++){
   int64_t value;
   if(singleton)value=chi[(x+p-V[0])%p];else if(full)value=0;else{
    value=work[x];if(x+p<2*p-1)value=(value+work[x+p])%MOD;if(value>MOD/2)value-=MOD;
   }
   require(abs(value)<=int64_t(V.size()),"recovered coefficient exceeds exact a priori bound");
   signed_sums[cell_of[x]][v]+=value;global_sum+=value;
   if(sample_idx<samples.size() && x==samples[sample_idx])sample_values[sample_idx++][v]=value;
  }
  require(global_sum==0,"character convolution does not have zero global sum");
 }
 vector<vector<array<int64_t,3>>>counts(c,vector<array<int64_t,3>>(c));
 int64_t global[3]={0,0,0};
 for(size_t u=0;u<c;u++)for(size_t v=0;v<c;v++){
  int64_t diagonal=(u==v)?cells[u].size():0,total=int64_t(cells[u].size())*cells[v].size()-diagonal,sg=signed_sums[u][v];
  require(sg==signed_sums[v][u],"ordered relation symmetry failed");
  require(abs(sg)<=total && (total+sg)%2==0,"sign totals fail parity or positivity");
  counts[u][v]={(total-sg)/2,diagonal,(total+sg)/2};
  for(int k=0;k<3;k++)global[k]+=counts[u][v][k];
 }
 require(global[1]==p && global[0]==int64_t(p)*(p-1)/2 && global[2]==global[0],"global Paley pair counts fail");
 cout<<"{\"p\":"<<p<<",\"marks\":";pattern_json(marks);cout<<",\"cells\":[";
 for(size_t i=0;i<c;i++){if(i)cout<<',';cout<<"{\"pattern\":";pattern_json(patterns[i]);cout<<",\"size\":"<<cells[i].size()<<'}';}
 cout<<"],\"relations\":[";bool first=true;
 for(size_t u=0;u<c;u++)for(size_t v=0;v<c;v++)for(int sign=-1;sign<=1;sign++)if(counts[u][v][sign+1]){
  if(!first)cout<<',';first=false;cout<<"{\"left_pattern\":";pattern_json(patterns[u]);cout<<",\"right_pattern\":";pattern_json(patterns[v]);cout<<",\"sign\":"<<sign<<",\"count\":"<<counts[u][v][sign+1]<<'}';
 }
 cout<<"],\"sampled_rows\":[";
 for(size_t i=0;i<samples.size();i++){if(i)cout<<',';cout<<"{\"x\":"<<samples[i]<<",\"cell_character_sums\":[";for(size_t j=0;j<c;j++){if(j)cout<<',';cout<<sample_values[i][j];}cout<<"]}";}
 double elapsed=chrono::duration<double>(chrono::steady_clock::now()-start).count();
 cout<<"],\"metadata\":{\"ordered_row_pairs\":true,\"field_primality_method\":\"trial division through sqrt(p)\",\"ntt_modulus\":"<<MOD<<",\"ntt_primitive_root\":"<<ROOT<<",\"ntt_length\":"<<L<<",\"ntt_transforms\":"<<transforms<<",\"singleton_cells_direct\":"<<singletons<<",\"full_cells_direct\":"<<fullcells<<",\"minimum_integer_recovery_margin\":"<<minmargin<<",\"pair_count_bit_bound\":";
 int bits=0;for(uint64_t val=uint64_t(p)*p;val;val>>=1)bits++;cout<<bits<<",\"backend_elapsed_seconds\":"<<elapsed<<"}}\n";
 return 0;
}catch(const exception &e){cerr<<e.what()<<'\n';return 1;}}
