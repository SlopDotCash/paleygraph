// One exact character convolution for Q=h2^T S h3, masked off three marks.
#include "ntt_exact.hpp"
#include <array>
#include <chrono>
#include <iostream>
#include <map>
#include <set>
using namespace std;
int main(int argc,char**argv){try{
 auto start=chrono::steady_clock::now();exact::require(argc==5,"usage: single_contraction p mark1 mark2 mark3");
 size_t used;uint64_t q=stoull(argv[1],&used);exact::require(used==string(argv[1]).size() && q>=5 && q%4==1 && q<=(exact::max_length+1)/2,"field exceeds supported exact prime domain");
 exact::require(exact::prime(q),"field parameter is not prime");uint32_t p=q;array<int,3>marks;
 for(int i=0;i<3;i++){int64_t m=stoll(argv[i+2],&used);exact::require(used==string(argv[i+2]).size() && 0<=m && m<p,"noncanonical mark");marks[i]=m;}
 exact::require(set<int>(marks.begin(),marks.end()).size()==3,"marks must be distinct");exact::validate();
 size_t L=1;while(L<2*size_t(p)-1)L<<=1;exact::require(L<=exact::max_length,"transform too long");
 vector<int8_t>chi(p,-1),h2(p),h3(p);chi[0]=0;for(uint64_t x=1;x<=(p-1)/2;x++)chi[x*x%p]=1;
 map<array<int,3>,uint32_t>cells;int64_t l1=0,h2l1=0;
 for(uint32_t x=0;x<p;x++){
  array<int,3>a;for(int i=0;i<3;i++)a[i]=chi[(x+p-marks[i])%p];cells[a]++;
  if(a[0] && a[1] && a[2]){h3[x]=a[0]*a[1]*a[2];h2[x]=a[0]*a[1]+a[0]*a[2]+a[1]*a[2];}
  l1+=abs(int(h3[x]));h2l1+=abs(int(h2[x]));
 }
 exact::require(l1==p-3 && 2*l1<exact::modulus,"ambiguous signed convolution recovery");
 vector<int64_t>a(L),b(L);for(uint32_t x=0;x<p;x++){a[x]=(chi[x]+exact::modulus)%exact::modulus;b[x]=(h3[x]+exact::modulus)%exact::modulus;}
 exact::ntt(a,false);exact::ntt(b,false);for(size_t x=0;x<L;x++)a[x]=a[x]*b[x]%exact::modulus;exact::ntt(a,true);
 int64_t Q=0,global=0;vector<pair<uint32_t,int64_t>>samples;set<uint32_t>samplex{0,1,2,p-1};for(int m:marks)samplex.insert(m);
 uint64_t state=20260905;for(int i=0;i<29;i++){state=(1664525*state+1013904223)%4294967296ULL;samplex.insert(state%p);}
 for(uint32_t x=0;x<p;x++){
  int64_t v=a[x];if(x+p<2*p-1)v=(v+a[x+p])%exact::modulus;if(v>exact::modulus/2)v-=exact::modulus;
  exact::require(abs(v)<=l1,"coefficient exceeds exact signed bound");Q+=int64_t(h2[x])*v;global+=v;
  if(samplex.count(x))samples.emplace_back(x,v);
 }
 exact::require(global==0 && abs(Q)<=h2l1*l1,"global convolution or contraction bound failed");
 cout<<"{\"p\":"<<p<<",\"marks\":["<<marks[0]<<','<<marks[1]<<','<<marks[2]<<"],\"Q\":"<<Q<<",\"cells\":[";bool first=true;
 for(auto&c:cells){if(!first)cout<<',';first=false;cout<<"{\"pattern\":["<<c.first[0]<<','<<c.first[1]<<','<<c.first[2]<<"],\"size\":"<<c.second<<'}';}
 cout<<"],\"sampled_convolution_rows\":[";first=true;for(auto&kv:samples){if(!first)cout<<',';first=false;cout<<"{\"x\":"<<kv.first<<",\"value\":"<<kv.second<<'}';}
 cout<<"],\"metadata\":{\"ntt_modulus\":"<<exact::modulus<<",\"ntt_length\":"<<L<<",\"transforms\":3,\"h3_l1\":"<<l1<<",\"h2_l1\":"<<h2l1<<",\"integer_recovery_margin\":"<<exact::modulus-2*l1<<",\"Q_absolute_bound\":"<<h2l1*l1<<",\"backend_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}}\n";
 return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
