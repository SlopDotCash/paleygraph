// Exact normalized prime-Paley triple-row trace inventory.
// Standard NTT, copied with provenance from the frozen round4 backend.
#include "ntt_exact.hpp"
#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
using namespace std;
int main(int argc,char**argv){try{
 auto start=chrono::steady_clock::now();exact::require(argc==3,"usage: trace_inventory p output_tau.bin");
 size_t used=0;uint64_t q=stoull(argv[1],&used);exact::require(used==string(argv[1]).size() && q>=13 && q%4==1 && q<=(exact::max_length+1)/2,"p outside prime Paley degree6 transform domain");
 exact::require(exact::prime(q),"field parameter is not prime");exact::validate();uint32_t p=q;
 size_t L=1;while(L<2*size_t(p)-1)L<<=1;exact::require(L<=exact::max_length,"unsupported transform length");
 vector<int8_t>chi(p,-1),input(p);chi[0]=0;for(uint64_t x=1;x<=(p-1)/2;x++)chi[x*x%p]=1;
 vector<int64_t>a(L),b(L);int64_t l1=0;
 for(uint32_t x=0;x<p;x++){input[x]=chi[x]*chi[(x+p-1)%p];l1+=abs(int(input[x]));a[x]=(chi[x]+exact::modulus)%exact::modulus;b[x]=(input[x]+exact::modulus)%exact::modulus;}
 exact::require(l1==p-2 && 2*l1<exact::modulus,"ambiguous centered integer recovery");
 exact::ntt(a,false);exact::ntt(b,false);for(size_t i=0;i<L;i++)a[i]=a[i]*b[i]%exact::modulus;exact::ntt(a,true);
 vector<int32_t>tau(p);map<array<int,4>,uint32_t>inventory;int64_t sum=0;int maxabs=0;
 for(uint32_t t=0;t<p;t++){
  int64_t value=a[t];if(t+p<2*p-1)value=(value+a[t+p])%exact::modulus;if(value>exact::modulus/2)value-=exact::modulus;
  exact::require(abs(value)<=l1,"coefficient exceeds exact input l1 bound");tau[t]=value;sum+=value;
  if(t<2){exact::require(value==-1,"singular parameter identity failed");continue;}
  // Hasse is a consistency check AFTER exact recovery; not used to recover it.
  exact::require(value*value<=4*int64_t(p),"Legendre Hasse consistency failed");
  maxabs=max(maxabs,int(abs(value)));int e01=1,e02=chi[t],e12=chi[t-1];
  // Direct row ordering: (S[0,1],S[0,t],S[1,t]).
  exact::require(e01==chi[p-1] && e02==chi[p-t] && e12==chi[(p+1-t)%p],"normalized edge orientation mismatch");
  int64_t bulk=0;
  for(int u:{-1,1})for(int v:{-1,1})for(int w:{-1,1}){
   int64_t num=int64_t(p)-3+u*(-e01-e02)+v*(-e01-e12)+w*(-e02-e12)
    +u*v*(-1-e02*e12)+u*w*(-1-e01*e12)+v*w*(-1-e01*e02)+u*v*w*value;
   exact::require(num>=0 && num%8==0,"nonintegral or negative triple-type inversion");bulk+=num/8;
  }
  exact::require(bulk==p-3,"bulk type counts fail total");inventory[{e01,e02,e12,int(value)}]++;
 }
 exact::require(sum==0,"total convolution sum failed");
 ofstream file(argv[2],ios::binary);exact::require(bool(file),"cannot open tau output");
 // Explicit little-endian encoding, including signed int32 values via uint32.
 auto put=[&](uint32_t n){for(int k=0;k<4;k++)file.put(char((n>>(8*k))&255));};put(p);for(int32_t v:tau)put(uint32_t(v));file.close();exact::require(bool(file),"tau output write failed");
 set<uint32_t>samples{0,1,2,p-1};uint64_t state=20260905;
 for(int i=0;i<29;i++){state=(1664525*state+1013904223)%4294967296ULL;samples.insert(state%p);}
 cout<<"{\"p\":"<<p<<",\"q\":"<<p<<",\"records\":[";bool first=true;
 for(auto&kv:inventory){if(!first)cout<<',';first=false;cout<<"{\"edges\":["<<kv.first[0]<<','<<kv.first[1]<<','<<kv.first[2]<<"],\"tau\":"<<kv.first[3]<<",\"count\":"<<kv.second<<'}';}
 cout<<"],\"normalized_parameter_count\":"<<p-2<<",\"ordered_distinct_row_weight\":"<<uint64_t(p)*(p-1)<<",\"sampled_tau_rows\":[";first=true;
 for(auto t:samples){if(!first)cout<<',';first=false;cout<<"{\"t\":"<<t<<",\"tau\":"<<tau[t]<<'}';}
 cout<<"],\"metadata\":{\"ntt_modulus\":"<<exact::modulus<<",\"ntt_length\":"<<L<<",\"transforms\":3,\"input_l1\":"<<l1<<",\"integer_recovery_margin\":"<<exact::modulus-2*l1<<",\"singular_tau_0\":"<<tau[0]<<",\"singular_tau_1\":"<<tau[1]<<",\"max_abs_nonsingular_tau\":"<<maxabs<<",\"Hasse_checked_after_recovery\":true,\"all_normalized_triple_type_inversions_checked\":"<<p-2<<",\"backend_seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}}\n";
 return 0;
}catch(const exception&e){cerr<<e.what()<<'\n';return 1;}}
