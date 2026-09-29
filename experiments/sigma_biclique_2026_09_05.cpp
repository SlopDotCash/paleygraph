// sigma_biclique_2026_09_05.cpp -- exploration helper for the sigma 'biclique' direction.
// Exact clique numbers of Paley graphs, exact profiles M_k(p) = max_{|A|=k} |B(A)| where
// B(A) = {b : a+b in Q ∪ {0} for all a in A}, exact balanced biclique numbers for small p,
// threshold searches, geometric-progression comparisons, and randomized lower bounds.
// Output: one JSON object per line.  Compile: clang++ -O3 -march=native -std=c++17 -o sb sigma_biclique_2026_09_05.cpp
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <chrono>
#include <random>
#include <string>
using namespace std;
typedef uint64_t u64;
static int p, W;
static vector<int> chi;
static vector<vector<u64>> adj;   // Paley adjacency: adj[x] = {y : chi(x-y) = 1}
static vector<vector<u64>> Np;    // Np[a] = {b : chi(a+b) = 1} ∪ {-a}
static inline int popc(const vector<u64>& v){ int c=0; for(int w=0;w<W;w++) c+=__builtin_popcountll(v[w]); return c; }
static inline int popc_and(const vector<u64>& a, const vector<u64>& b){ int c=0; for(int w=0;w<W;w++) c+=__builtin_popcountll(a[w]&b[w]); return c; }
static inline bool empty_bs(const vector<u64>& v){ for(int w=0;w<W;w++) if(v[w]) return false; return true; }
static inline int first_bit(const vector<u64>& v){ for(int w=0;w<W;w++) if(v[w]) return w*64+__builtin_ctzll(v[w]); return -1; }
static inline void clr(vector<u64>& v,int i){ v[i>>6] &= ~(1ULL<<(i&63)); }
static inline void setb(vector<u64>& v,int i){ v[i>>6] |= (1ULL<<(i&63)); }
static inline bool get(const vector<u64>& v,int i){ return (v[i>>6]>>(i&63))&1ULL; }
static vector<int> bits(const vector<u64>& v){ vector<int> r; for(int i=0;i<p;i++) if(get(v,i)) r.push_back(i); return r; }
static string jlist(const vector<int>& v){ string s="["; for(size_t i=0;i<v.size();i++){ if(i) s+=","; s+=to_string(v[i]); } return s+"]"; }
static double now(){ return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); }

static void setup(int pp){
  p=pp; W=(p+63)/64; chi.assign(p,-1); chi[0]=0;
  for(long long x=1;x<p;x++) chi[(x*x)%p]=1;
  adj.assign(p, vector<u64>(W,0)); Np.assign(p, vector<u64>(W,0));
  for(int x=0;x<p;x++) for(int y=0;y<p;y++){ if(chi[((x-y)%p+p)%p]==1) setb(adj[x],y); int s=(x+y)%p; if(chi[s]==1||s==0) setb(Np[x],y); }
}

// ---------------- maximum clique (BBMC: bitset branch and bound with greedy colouring) ----------------
static int best_omega; static vector<int> best_clique, cur_clique; static long long clique_nodes;
static void expand(vector<u64>& P, int depth){
  clique_nodes++;
  vector<int> order, color; vector<u64> U=P; int col=0;
  while(!empty_bs(U)){ col++; vector<u64> Qc=U; while(!empty_bs(Qc)){ int v=first_bit(Qc); clr(Qc,v); clr(U,v); for(int w=0;w<W;w++) Qc[w]&=~adj[v][w]; order.push_back(v); color.push_back(col);} }
  for(int i=(int)order.size()-1;i>=0;i--){
    if(depth+color[i]<=best_omega) return;
    int v=order[i]; vector<u64> NP(W); for(int w=0;w<W;w++) NP[w]=P[w]&adj[v][w];
    cur_clique.push_back(v);
    if(empty_bs(NP)){ if(depth+1>best_omega){ best_omega=depth+1; best_clique=cur_clique; } }
    else expand(NP,depth+1);
    cur_clique.pop_back(); clr(P,v);
  }
}
static int omega_exact(){
  // every clique of size >= 2 is equivalent under x -> (x-a)/(b-a) (b-a a square, since a~b) to one containing 0 and 1
  vector<u64> P(W,0); for(int x=2;x<p;x++) if(chi[x]==1 && chi[x-1]==1) setb(P,x);
  best_omega=2; best_clique={0,1}; cur_clique={0,1}; clique_nodes=0;
  // heuristic warm start
  mt19937_64 rng(12345);
  for(int it=0;it<200;it++){ vector<u64> C=P; vector<int> cl={0,1}; while(!empty_bs(C)){ vector<int> b=bits(C); int v=b[rng()%b.size()]; cl.push_back(v); for(int w=0;w<W;w++) C[w]&=adj[v][w]; } if((int)cl.size()>best_omega){best_omega=cl.size(); best_clique=cl;} }
  expand(P,2);
  return best_omega;
}

// ---------------- profile M_k(p): max |B(A)| over |A| = k  ----------------
static int Kt, bestM; static vector<int> bestA, curA; static long long prof_nodes; static double deadline; static bool timed_out;
static void dfs_profile(int j, int last, const vector<u64>& B, int e){
  prof_nodes++;
  if((prof_nodes&1023)==0 && now()>deadline){ timed_out=true; }
  if(timed_out) return;
  int cnt=popc(B);
  if(j==Kt){ if(cnt>bestM){ bestM=cnt; bestA=curA; } return; }
  if(cnt<=bestM) return;
  int r=Kt-j;
  vector<pair<int,int>> cand; cand.reserve(p);
  for(int a=last+1;a<p;a++){ if(a==e) continue; int v=popc_and(B,Np[a]); if(v>bestM) cand.push_back({a,v}); }
  if((int)cand.size()<r) return;
  { vector<int> vs; vs.reserve(cand.size()); for(auto&c:cand) vs.push_back(c.second); nth_element(vs.begin(), vs.begin()+(r-1), vs.end(), greater<int>()); if(vs[r-1]<=bestM) return; }
  for(auto&c:cand){ if(c.second<=bestM) continue; vector<u64> NB(W); for(int w=0;w<W;w++) NB[w]=B[w]&Np[c.first][w]; curA.push_back(c.first); dfs_profile(j+1,c.first,NB,e); curA.pop_back(); if(timed_out) return; }
}
// returns exact M_k (bestM) unless timed out; start value lb (search finds > lb only; pass lb = threshold-1 for threshold search)
static int profile(int k, int lb, double tlimit, bool& exact){
  int n0=2; while(chi[n0]!=-1) n0++;
  Kt=k; bestM=lb; bestA.clear(); prof_nodes=0; deadline=now()+tlimit; timed_out=false;
  for(int e: {1,n0}){
    vector<u64> B(W); for(int w=0;w<W;w++) B[w]=Np[0][w]&Np[e][w];
    curA={0,e};
    if(k==2){ int c=popc(B); if(c>bestM){bestM=c; bestA=curA;} continue; }
    dfs_profile(2,0,B,e);
  }
  exact=!timed_out; return bestM;
}

// ---------------- balanced biclique number b(p) = max min(|A|,|B(A)|) ----------------
static int bestBal; static vector<int> bestBalA; static long long bal_nodes;
static void dfs_bal(int j, int last, const vector<u64>& B, int e){
  bal_nodes++;
  if((bal_nodes&1023)==0 && now()>deadline) timed_out=true;
  if(timed_out) return;
  int cnt=popc(B);
  if(min(j,cnt)>bestBal){ bestBal=min(j,cnt); bestBalA=curA; }
  if(cnt<=bestBal) return;
  vector<pair<int,int>> cand; cand.reserve(p);
  for(int a=last+1;a<p;a++){ if(a==e) continue; int v=popc_and(B,Np[a]); if(v>bestBal) cand.push_back({a,v}); }
  if(cand.empty()) return;
  { vector<int> vs; for(auto&c:cand) vs.push_back(c.second); sort(vs.begin(),vs.end(),greater<int>()); int bound=0; for(size_t r=1;r<=vs.size();r++) bound=max(bound,min(j+(int)r,vs[r-1])); if(bound<=bestBal) return; }
  for(auto&c:cand){ if(c.second<=bestBal) continue; vector<u64> NB(W); for(int w=0;w<W;w++) NB[w]=B[w]&Np[c.first][w]; curA.push_back(c.first); dfs_bal(j+1,c.first,NB,e); curA.pop_back(); if(timed_out) return; }
}
static int balanced(int lb, double tlimit, bool& exact){
  int n0=2; while(chi[n0]!=-1) n0++;
  bestBal=lb; bestBalA.clear(); bal_nodes=0; deadline=now()+tlimit; timed_out=false;
  for(int e: {1,n0}){ vector<u64> B(W); for(int w=0;w<W;w++) B[w]=Np[0][w]&Np[e][w]; curA={0,e}; dfs_bal(2,0,B,e); }
  exact=!timed_out; return bestBal;
}

// ---------------- global product max over min(|A|,|B|) >= 2 ----------------
static long long bestProd; static vector<int> bestProdA; static long long prod_nodes;
static void dfs_prod(int j, int last, const vector<u64>& B, int e){
  prod_nodes++;
  if((prod_nodes&1023)==0 && now()>deadline) timed_out=true;
  if(timed_out) return;
  int cnt=popc(B);
  if(cnt>=2 && (long long)j*cnt>bestProd){ bestProd=(long long)j*cnt; bestProdA=curA; }
  if(cnt<2) return;
  vector<pair<int,int>> cand; cand.reserve(p);
  for(int a=last+1;a<p;a++){ if(a==e) continue; int v=popc_and(B,Np[a]); if(v>=2) cand.push_back({a,v}); }
  if(cand.empty()) return;
  { vector<int> vs; for(auto&c:cand) vs.push_back(c.second); sort(vs.begin(),vs.end(),greater<int>()); long long bound=0; for(size_t r=1;r<=vs.size();r++) bound=max(bound,(long long)(j+(int)r)*vs[r-1]); if(bound<=bestProd) return; }
  { vector<int> vs; for(auto&c:cand) vs.push_back(c.second); sort(vs.begin(),vs.end(),greater<int>());
    for(auto&c:cand){ long long cb=0; for(size_t r=1;r<=vs.size();r++) cb=max(cb,(long long)(j+(int)r)*min(c.second,vs[r-1])); if(cb<=bestProd) continue; vector<u64> NB(W); for(int w=0;w<W;w++) NB[w]=B[w]&Np[c.first][w]; curA.push_back(c.first); dfs_prod(j+1,c.first,NB,e); curA.pop_back(); if(timed_out) return; } }
}
static long long product_max(double tlimit, bool& exact){
  int n0=2; while(chi[n0]!=-1) n0++;
  bestProd=0; bestProdA.clear(); prod_nodes=0; deadline=now()+tlimit; timed_out=false;
  for(int e: {1,n0}){ vector<u64> B(W); for(int w=0;w<W;w++) B[w]=Np[0][w]&Np[e][w]; curA={0,e}; dfs_prod(2,0,B,e); }
  exact=!timed_out; return bestProd;
}



// ---------------- sum-clique number s(p) = max |A| with A + A ⊆ Q ∪ {0} ----------------
// (A, A) is then a complete biclique, so s(p) <= b(p).  Vertices x with 2x ∈ Q ∪ {0}; x ~ y iff chi(x+y) = 1.
// 0 is adjacent to every x with chi(x)=1.  Scaling by squares preserves the relation, so a sum-clique with a
// nonzero element can be normalised to contain v0 = 1 (if 2 ∈ Q, all vertices lie in Q ∪ {0}) or v0 = n0 (if 2 ∉ Q).
static vector<vector<u64>> sadj;
static void sexpand(vector<u64>& P, int depth){
  clique_nodes++;
  vector<int> order, color; vector<u64> U=P; int col=0;
  while(!empty_bs(U)){ col++; vector<u64> Qc=U; while(!empty_bs(Qc)){ int v=first_bit(Qc); clr(Qc,v); clr(U,v); for(int w=0;w<W;w++) Qc[w]&=~sadj[v][w]; order.push_back(v); color.push_back(col);} }
  for(int i=(int)order.size()-1;i>=0;i--){
    if(depth+color[i]<=best_omega) return;
    int v=order[i]; vector<u64> NP(W); for(int w=0;w<W;w++) NP[w]=P[w]&sadj[v][w];
    cur_clique.push_back(v);
    if(empty_bs(NP)){ if(depth+1>best_omega){ best_omega=depth+1; best_clique=cur_clique; } }
    else sexpand(NP,depth+1);
    cur_clique.pop_back(); clr(P,v);
  }
}
static int sumclique_exact(){
  int n0=2; while(chi[n0]!=-1) n0++;
  sadj.assign(p, vector<u64>(W,0));
  for(int x=0;x<p;x++){ if(!(x==0 || chi[(2*x)%p]==1)) continue; for(int y=0;y<p;y++){ if(y==x) continue; if(!(y==0 || chi[(2*y)%p]==1)) continue; if((x+y)%p==0 || chi[(x+y)%p]==1) setb(sadj[x],y); } }
  int v0 = (chi[2]==1) ? 1 : n0;   // normalised nonzero element
  best_omega=1; best_clique={0}; cur_clique={v0}; clique_nodes=0;
  vector<u64> P=sadj[v0];
  if(best_omega<1){best_omega=1;}
  best_omega=1; best_clique={v0};
  sexpand(P,1);
  return best_omega;
}

// ---------------- HP-reduced product check ----------------
// By Hanson-Petridis Thm 1.2 (d=(p-1)/2): |A||B| <= (p-1)/2 + |B ∩ (-A)|.  A biclique with |A||B| >= (p+5)/2
// has |B ∩ (-A)| >= 3, and C = A ∩ (-B) is a clique; normalising three of its elements to {0,1,c} ⊆ A gives
// {0,-1,-c} ⊆ B, so A ⊆ CA := Np[0]∩Np[-1]∩Np[-c] and B ⊆ CB := Np[0]∩Np[1]∩Np[c].
static long long p3_best; static vector<int> p3_A; static long long p3_nodes; static vector<int> CAlist;
static int p3_minB;   // a violator has |A| <= |B| and |A||B| >= (p+5)/2, so |B| >= ceil(sqrt((p+5)/2))
static void dfs_p3(int idx, const vector<u64>& B, int j){
  p3_nodes++;
  int cnt=popc(B);
  if(cnt<p3_minB || j>cnt) return;          // B is the larger side of the violator
  if(cnt>=3 && (long long)j*cnt>p3_best){ p3_best=(long long)j*cnt; p3_A=curA; }
  vector<pair<int,int>> cand;
  for(int i=idx;i<(int)CAlist.size();i++){ int v=popc_and(B,Np[CAlist[i]]); if(v>=p3_minB) cand.push_back({i,v}); }
  if(cand.empty()) return;
  vector<int> vs; for(auto&c:cand) vs.push_back(c.second); sort(vs.begin(),vs.end(),greater<int>());
  long long bound=0; for(size_t r=1;r<=vs.size();r++) bound=max(bound,(long long)(j+(int)r)*vs[r-1]);
  if(bound<=p3_best) return;
  for(auto&c:cand){ long long cb=0; for(size_t r=1;r<=vs.size();r++) cb=max(cb,(long long)(j+(int)r)*min(c.second,vs[r-1])); if(cb<=p3_best) continue;
    vector<u64> NB(W); for(int w=0;w<W;w++) NB[w]=B[w]&Np[CAlist[c.first]][w]; curA.push_back(CAlist[c.first]); dfs_p3(c.first+1,NB,j+1); curA.pop_back(); }
}
static long long product3(long long threshold, int& ntri){
  p3_best=threshold; p3_A.clear(); p3_nodes=0; ntri=0;
  p3_minB=1; while((long long)p3_minB*p3_minB<(p+5)/2) p3_minB++;
  for(int c=2;c<p;c++){ if(!(chi[c]==1 && chi[c-1]==1)) continue; ntri++;
    vector<u64> CA(W), CB(W); for(int w=0;w<W;w++){ CA[w]=Np[0][w]&Np[p-1][w]&Np[p-c][w]; CB[w]=Np[0][w]&Np[1][w]&Np[c][w]; }
    clr(CA,0); clr(CA,1); clr(CA,c); CAlist=bits(CA); curA={0,1,c}; dfs_p3(0,CB,3); }
  return p3_best;
}

// ---------------- geometric progressions and random greedy ----------------
static int gp_best(int k, int& bestr, int& bestc){
  int n0=2; while(chi[n0]!=-1) n0++;
  int best=-1; bestr=-1; bestc=-1;
  for(int r=2;r<p;r++){
    // powers r^0..r^{k-1} must be distinct
    vector<int> U; long long x=1; bool ok=true; for(int i=0;i<k;i++){ for(int u:U) if(u==x){ok=false;break;} if(!ok) break; U.push_back((int)x); x=(x*r)%p; }
    if(!ok) continue;
    for(int c: {1,n0}){ vector<u64> B(W,~0ULL); for(int u:U){ int cu=(int)(((long long)c*u)%p); for(int w=0;w<W;w++) B[w]&=Np[cu][w]; } int v=popc(B); if(v>best){best=v; bestr=r; bestc=c;} }
  }
  return best;
}
static int random_greedy_balanced(int iters, vector<int>& outA, unsigned seed){
  mt19937_64 rng(seed); int best=0;
  for(int it=0;it<iters;it++){
    vector<u64> B(W,0); for(int i=0;i<p;i++) setb(B,i);
    vector<int> A; vector<u64> used(W,0);
    while(true){
      int cnt=popc(B); if(min((int)A.size(),cnt)>best){best=min((int)A.size(),cnt); outA=A;}
      if(cnt<= (int)A.size()) break;
      // pick among top-3 candidates by |B ∩ Np[a]| with randomisation
      vector<pair<int,int>> cand; for(int a=0;a<p;a++){ if(get(used,a)) continue; cand.push_back({-popc_and(B,Np[a]),a}); }
      if(cand.empty()) break; int m=min((size_t)4,cand.size()); partial_sort(cand.begin(),cand.begin()+m,cand.end());
      int pick=cand[rng()%m].second; A.push_back(pick); setb(used,pick); for(int w=0;w<W;w++) B[w]&=Np[pick][w];
    }
  }
  return best;
}

int main(int argc, char** argv){
  if(argc<3){ fprintf(stderr,"usage: sb p mode [args]\n modes: omega | profile kmin kmax tlimit | threshold k T tlimit | balanced tlimit | product tlimit | gp kmax | greedy iters\n"); return 1; }
  int pp=atoi(argv[1]); string mode=argv[2]; setup(pp);
  if(mode=="omega"){ double t0=now(); int om=omega_exact(); vector<u64> B(W,~0ULL); for(int c:best_clique) for(int w=0;w<W;w++) B[w]&=Np[c][w];
    printf("{\"p\":%d,\"omega\":%d,\"clique\":%s,\"nodes\":%lld,\"B_of_clique\":%d,\"time\":%.3f}\n",p,om,jlist(best_clique).c_str(),clique_nodes,popc(B),now()-t0); }
  else if(mode=="profile"){ int kmin=atoi(argv[3]), kmax=atoi(argv[4]); double tl=atof(argv[5]);
    for(int k=kmin;k<=kmax;k++){ double t0=now(); bool ex; int m=profile(k,k-1>0?0:0,tl,ex); // lb 0: exact
      printf("{\"p\":%d,\"k\":%d,\"M\":%d,\"exact\":%s,\"A\":%s,\"nodes\":%lld,\"time\":%.3f}\n",p,k,m,ex?"true":"false",jlist(bestA).c_str(),prof_nodes,now()-t0); fflush(stdout); if(!ex) break; } }
  else if(mode=="threshold"){ int k=atoi(argv[3]), T=atoi(argv[4]); double tl=atof(argv[5]); double t0=now(); bool ex; int m=profile(k,T-1,tl,ex);
    printf("{\"p\":%d,\"k\":%d,\"T\":%d,\"found\":%s,\"M_ge_T\":%d,\"exact\":%s,\"A\":%s,\"nodes\":%lld,\"time\":%.3f}\n",p,k,T,(m>=T)?"true":"false",m,ex?"true":"false",jlist(bestA).c_str(),prof_nodes,now()-t0); }
  else if(mode=="balanced"){ double tl=atof(argv[3]); double t0=now(); bool ex; vector<int> A0; int lb=random_greedy_balanced(200,A0,7); int b=balanced(lb-1,tl,ex);
    printf("{\"p\":%d,\"balanced\":%d,\"exact\":%s,\"A\":%s,\"greedy_lb\":%d,\"greedy_A\":%s,\"nodes\":%lld,\"time\":%.3f}\n",p,b,ex?"true":"false",jlist(bestBalA.empty()?A0:bestBalA).c_str(),lb,jlist(A0).c_str(),bal_nodes,now()-t0); }
  else if(mode=="product"){ double tl=atof(argv[3]); double t0=now(); bool ex; long long pr=product_max(tl,ex);
    printf("{\"p\":%d,\"product\":%lld,\"exact\":%s,\"A\":%s,\"nodes\":%lld,\"time\":%.3f}\n",p,pr,ex?"true":"false",jlist(bestProdA).c_str(),prod_nodes,now()-t0); }
  else if(mode=="sumclique"){ double t0=now(); int s=sumclique_exact();
    printf("{\"p\":%d,\"sumclique\":%d,\"A\":%s,\"nodes\":%lld,\"time\":%.3f}\n",p,s,jlist(best_clique).c_str(),clique_nodes,now()-t0); }
  else if(mode=="product3"){ double t0=now(); int ntri; long long thr=(p+3)/2; long long r=product3(thr,ntri);
    printf("{\"p\":%d,\"threshold\":%lld,\"violator_found\":%s,\"best\":%lld,\"A\":%s,\"triangles\":%d,\"nodes\":%lld,\"time\":%.3f}\n",p,thr,(r>thr)?"true":"false",r,jlist(p3_A).c_str(),ntri,p3_nodes,now()-t0); }
  else if(mode=="gp"){ int kmax=atoi(argv[3]); printf("{\"p\":%d,\"gp\":[",p); for(int k=2;k<=kmax;k++){ int r,c; int v=gp_best(k,r,c); printf("%s[%d,%d,%d,%d]",k>2?",":"",k,v,r,c);} printf("]}\n"); }
  else if(mode=="greedy"){ int iters=atoi(argv[3]); vector<int> A; int b=random_greedy_balanced(iters,A,99); vector<u64> B(W,~0ULL); for(int a:A) for(int w=0;w<W;w++) B[w]&=Np[a][w];
    printf("{\"p\":%d,\"greedy_balanced\":%d,\"A\":%s,\"B\":%s}\n",p,b,jlist(A).c_str(),jlist(bits(B)).c_str()); }
  return 0;
}
