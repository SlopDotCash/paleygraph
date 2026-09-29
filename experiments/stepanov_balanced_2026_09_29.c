// stepanov_balanced_2026_09_29.c  (worker `balanced`, Stepanov wave, 2026-09-29)
//
// Exact biclique data for complete bicliques A + B subset of Q u {0} in F_p (p = 1 mod 4 prime).
// Relation: b in N'(a)  iff  a + b in Q u {0}.  B(A) = intersection of N'(a), a in A.
//
// Modes:
//   profile p k        : exact M_k(p) = max_{|A| = k} |B(A)| with a witness A.
//   certify p m1 T     : exhaustive branch and bound; decides whether some complete biclique
//                        with m1 <= |A| <= |B| and |A||B| > T exists.  Prints a violator if found.
//
// Normalisation (valid for |A| >= 2): the maps (A,B) -> (uA + t, uB - t), u in Q, preserve the
// relation, so A may be assumed to contain {0,1} (if two elements of A differ by a residue) or
// {0,n0} (n0 the least non-residue) otherwise.  The search enumerates supersets of {0,s},
// s in {1,n0}, adding further elements in increasing order, so every k-set is covered up to the
// symmetries.  Pruning (both modes) only uses: B(A u X) is contained in B(A u {x}) for x in X.
//
// Compile:  clang -O3 -o sb experiments/stepanov_balanced_2026_09_29.c
// Output: one JSON line on stdout.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

typedef uint64_t u64;
static int p, W, n0;
static u64 *N;          // N[x*W .. ] bitset of N'(x)
static signed char *chi;
static long long nodes;

static int popc(const u64 *b){ int c=0; for(int w=0;w<W;w++) c+=__builtin_popcountll(b[w]); return c; }
static int popand(const u64 *a, const u64 *b){ int c=0; for(int w=0;w<W;w++) c+=__builtin_popcountll(a[w]&b[w]); return c; }

static void setup(int pp){
  p=pp; W=(p+63)/64;
  chi=calloc(p,1);
  for(int x=1;x<p;x++) chi[(long long)x*x%p]=1;
  for(int x=1;x<p;x++) if(!chi[x]) chi[x]=-1;
  n0=0; for(int x=2;x<p;x++) if(chi[x]==-1){n0=x;break;}
  N=calloc((size_t)p*W,sizeof(u64));
  for(int a=0;a<p;a++) for(int b=0;b<p;b++){ int s=(a+b)%p; if(s==0||chi[s]==1) N[(size_t)a*W+b/64]|=1ULL<<(b%64); }
}

// ---------------- profile ----------------
static int K, bestM, bestA[64], curA[64];
static void dfs_profile(int j, int last, const u64 *B, const int *cand, int nc){
  nodes++;
  if(j==K){ int c=popc(B); if(c>bestM){ bestM=c; memcpy(bestA,curA,sizeof(int)*K);} return; }
  int need=K-j;
  // counts for candidates > last
  int *cnt=malloc(sizeof(int)*(nc+1)); int *cl=malloc(sizeof(int)*(nc+1)); int m=0;
  for(int i=0;i<nc;i++){ int x=cand[i]; if(x<=last) continue; int c=popand(B,N+(size_t)x*W); if(c>bestM){ cl[m]=x; cnt[m]=c; m++; } }
  if(m<need){ free(cnt); free(cl); return; }
  // bound: the need-th largest count must exceed bestM
  { int *tmp=malloc(sizeof(int)*m); memcpy(tmp,cnt,sizeof(int)*m);
    // partial selection of need-th largest
    for(int a=0;a<need;a++){ int mi=a; for(int b=a+1;b<m;b++) if(tmp[b]>tmp[mi]) mi=b; int t=tmp[a]; tmp[a]=tmp[mi]; tmp[mi]=t; }
    int vr=tmp[need-1]; free(tmp); if(vr<=bestM){ free(cnt); free(cl); return; } }
  u64 *NB=malloc(sizeof(u64)*W);
  for(int i=0;i<m;i++){
    if(cnt[i]<=bestM) continue;
    int x=cl[i];
    for(int w=0;w<W;w++) NB[w]=B[w]&N[(size_t)x*W+w];
    curA[j]=x;
    // children candidates: those after i with count > bestM (monotone filter)
    dfs_profile(j+1,x,NB,cl+i+1,m-i-1);
  }
  free(NB); free(cnt); free(cl);
}

static void run_profile(int k){
  K=k; bestM=0; nodes=0;
  if(k==1){ bestM=(p+1)/2; bestA[0]=0; return; }
  int *cand=malloc(sizeof(int)*p);
  for(int si=0;si<2;si++){
    int s= si==0?1:n0;
    u64 *B=malloc(sizeof(u64)*W);
    for(int w=0;w<W;w++) B[w]=N[w]&N[(size_t)s*W+w];
    int nc=0; for(int x=1;x<p;x++) if(x!=s) cand[nc++]=x;
    curA[0]=0; curA[1]=s;
    if(k==2){ int c=popc(B); if(c>bestM){bestM=c; bestA[0]=0; bestA[1]=s;} }
    else dfs_profile(2,0,B,cand,nc);
    free(B);
  }
  free(cand);
}

// ---------------- certify ----------------
static int M1; static long long T; static int found, foundA[256], foundm;
static int Mcap[512]; // Mcap[m] >= M_m(p) (exact profile values for m <= mk, then M_mk)
static int thr, BMAX; static double tlimit, tstart; static int timed_out;
static double now(void);
static void dfs_cert(int j, int last, const u64 *B, const int *cand, int nc){
  if(timed_out) return;
  nodes++;
  if((nodes & 1023)==0 && now()-tstart>tlimit){ timed_out=1; return; }
  int nb=popc(B);
  if(j>=M1 && j<=BMAX && nb>=j && (long long)j*nb>T){ found=1; foundm=j; memcpy(foundA,curA,sizeof(int)*j); T=(long long)j*nb; }
  if(j>=BMAX) return;
  int *cnt=malloc(sizeof(int)*(nc+1)); int *cl=malloc(sizeof(int)*(nc+1)); int m=0;
  for(int i=0;i<nc;i++){ int x=cand[i]; if(x<=last) continue; int c=popand(B,N+(size_t)x*W); if(c>=thr){ cl[m]=x; cnt[m]=c; m++; } }
  int *s=malloc(sizeof(int)*(m+1)); memcpy(s,cnt,sizeof(int)*m);
  for(int a=1;a<m;a++){ int v=s[a], b=a-1; while(b>=0 && s[b]<v){ s[b+1]=s[b]; b--; } s[b+1]=v; }
  int feasible=0;
  for(int r=1;r<=m;r++){ int mm=j+r; if(mm>BMAX) break; if(mm<M1) continue; int v=s[r-1]; if(mm<512 && Mcap[mm]<v) v=Mcap[mm]; if(v<mm) break; if((long long)mm*v>T){ feasible=1; break; } }
  free(s);
  if(!feasible){ free(cnt); free(cl); return; }
  u64 *NB=malloc(sizeof(u64)*W);
  for(int i=0;i<m && !timed_out;i++){
    int x=cl[i];
    for(int w=0;w<W;w++) NB[w]=B[w]&N[(size_t)x*W+w];
    curA[j]=x;
    dfs_cert(j+1,x,NB,cl+i+1,m-i-1);
  }
  free(NB); free(cnt); free(cl);
}

static void run_certify(int m1, long long t){
  M1=m1; T=t; found=0; nodes=0; timed_out=0; tstart=now();
  // final n satisfies n >= m >= m1 and m n > T, so n >= max(m1, floor(sqrt T)+1)
  long long q=0; while((q+1)*(q+1)<=T) q++;
  thr = m1 > q+1 ? m1 : (int)(q+1);
  int *cand=malloc(sizeof(int)*p);
  for(int si=0;si<2 && !found;si++){
    int s= si==0?1:n0;
    u64 *B=malloc(sizeof(u64)*W);
    for(int w=0;w<W;w++) B[w]=N[w]&N[(size_t)s*W+w];
    int nc=0; for(int x=1;x<p;x++) if(x!=s) cand[nc++]=x;
    curA[0]=0; curA[1]=s;
    dfs_cert(2,0,B,cand,nc);
    free(B);
  }
  free(cand);
}

static double now(void){ struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec+1e-9*ts.tv_nsec; }

int main(int argc, char **argv){
  if(argc<3){ fprintf(stderr,"usage: sb p profile k | sb p pk K mk bmax tlimit\n"); return 1; }
  int pp=atoi(argv[1]); setup(pp);
  double t0=now();
  if(!strcmp(argv[2],"profile")){ (void)0;
    int k=atoi(argv[3]); run_profile(k);
    printf("{\"p\":%d,\"mode\":\"profile\",\"k\":%d,\"M\":%d,\"A\":[",p,k,bestM);
    for(int i=0;i<k;i++) printf("%s%d",i?",":"",bestA[i]);
    printf("],\"nodes\":%lld,\"time\":%.3f}\n",nodes,now()-t0);
  } else if(!strcmp(argv[2],"cert")){
    // cert p m1 T bmax tlimit c_1 ... c_L : maximise |A||B| over complete bicliques with
    // m1 <= |A| <= |B|, |A| <= bmax, starting from the value T (only strictly larger products are
    // searched); c_m are valid upper bounds for M_m(p) (c_m for m > L is c_L).
    int m1=atoi(argv[3]); long long t=atoll(argv[4]); BMAX=atoi(argv[5]); tlimit=atof(argv[6]);
    int L=argc-7; for(int m=0;m<512;m++) Mcap[m]=p;
    for(int i=0;i<L;i++) Mcap[i+1]=atoi(argv[7+i]);
    for(int m=L+1;m<512;m++) Mcap[m]=Mcap[L];
    run_certify(m1,t);
    printf("{\"p\":%d,\"mode\":\"cert\",\"m1\":%d,\"T_in\":%lld,\"bmax\":%d,\"P\":%lld,\"exact\":%s,\"improved\":%s,\"A\":[",p,m1,t,BMAX,T,timed_out?"false":"true",found?"true":"false");
    if(found) for(int i=0;i<foundm;i++) printf("%s%d",i?",":"",foundA[i]);
    printf("],\"nodes\":%lld,\"time\":%.3f}\n",nodes,now()-t0);
  } else if(!strcmp(argv[2],"pk")){
    // P_K(p) = max{|A||B| : complete biclique, min(|A|,|B|) >= K}.
    // Step 1: exact M_m for K <= m <= mk.  T = max m*M_m over those m with M_m >= K.
    // Step 2: certify that no biclique with mk+1 <= |A| <= |B| has |A||B| > T
    //         (bound uses M_m <= Mcap[m]; M_m is nonincreasing in m).
    int KK=atoi(argv[3]), mk=atoi(argv[4]); BMAX=atoi(argv[5]); tlimit=atof(argv[6]);
    int Ms[64]; int As[64][64]; long long Tbest=0; int argm=0;
    printf("{\"p\":%d,\"mode\":\"pk\",\"K\":%d,\"mk\":%d,\"profile\":{",p,KK,mk);
    for(int m=KK;m<=mk;m++){ run_profile(m); Ms[m]=bestM; memcpy(As[m],bestA,sizeof(int)*m);
      printf("%s\"%d\":{\"M\":%d,\"A\":[",m>KK?",":"",m,bestM); for(int i=0;i<m;i++) printf("%s%d",i?",":"",bestA[i]);
      printf("],\"nodes\":%lld}",nodes);
      if(bestM>=KK && (long long)m*bestM>Tbest){ Tbest=(long long)m*bestM; argm=m; } }
    for(int m=0;m<512;m++) Mcap[m]= m<=mk ? (m>=KK?Ms[m]:p) : Ms[mk];
    double t1=now();
    run_certify(mk+1,Tbest);
    printf("},\"T_profile\":%lld,\"argm\":%d,\"bmax\":%d,\"P\":%lld,\"exact\":%s,\"improved_by_large_m\":%s,\"A_large\":[",Tbest,argm,BMAX,T,timed_out?"false":"true",found?"true":"false");
    if(found) for(int i=0;i<foundm;i++) printf("%s%d",i?",":"",foundA[i]);
    printf("],\"cert_nodes\":%lld,\"cert_time\":%.3f,\"time\":%.3f}\n",nodes,now()-t1,now()-t0);
  }
  return 0;
}
