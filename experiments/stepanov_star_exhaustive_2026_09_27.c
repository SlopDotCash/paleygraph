/* Independent exhaustive check of the robust Hanson-Petridis inequality (star):
   for every A subset F_p with 0 in A, 1<=m<=(p+1)/2, and every 0<=e<=(m-1)/2,
   sum_{b: e_b<=e} (e+1-e_b)(2m-3e-e_b-2 delta_b) <= 2(e+1)(d-e).
   Also records the tightest ratio for e>=1. No shared code with the worker's verifier. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc,char**argv){
  int p=atoi(argv[1]); int d=(p-1)/2; int chi[64]={0};
  for(int x=1;x<p;x++) chi[(x*x)%p]=1;
  for(int x=1;x<p;x++) if(!chi[x]) chi[x]=-1;
  uint64_t Nb[64]; for(int b=0;b<p;b++){Nb[b]=0; for(int a=0;a<p;a++) if(chi[(a+b)%p]==-1) Nb[b]|=1ULL<<a;}
  long long tested=0, fails=0; double worst=0; uint64_t worstA=0; int worste=-1;
  uint64_t full=1ULL<<(p-1);
  for(uint64_t mask=0; mask<full; mask++){
    uint64_t A=(mask<<1)|1ULL; int m=__builtin_popcountll(A);
    if(2*m-1>p) continue;
    int eb[64], del[64];
    for(int b=0;b<p;b++){ eb[b]=__builtin_popcountll(A&Nb[b]); del[b]=(A>>((p-b)%p))&1ULL; }
    for(int e=0; 2*e<=m-1; e++){
      long long lhs=0;
      for(int b=0;b<p;b++) if(eb[b]<=e) lhs+=(long long)(e+1-eb[b])*(2*m-3*e-eb[b]-2*del[b]);
      long long rhs=2LL*(e+1)*(d-e); tested++;
      if(lhs>rhs){ fails++; if(fails<5) printf("FAIL p=%d A=%llx e=%d lhs=%lld rhs=%lld\n",p,(unsigned long long)A,e,lhs,rhs);}
      if(e>=1 && rhs>0){ double r=(double)lhs/rhs; if(r>worst){worst=r;worstA=A;worste=e;} }
    }
  }
  printf("p=%d tested=%lld fails=%lld worst_ratio_e>=1=%.4f at A=%llx e=%d\n",p,tested,fails,worst,(unsigned long long)worstA,worste);
  return fails?1:0;
}
