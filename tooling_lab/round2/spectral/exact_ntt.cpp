// Exact cyclic character convolution for integer witnesses.
// Established radix-2 NTT; this implementation is a certificate backend.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using namespace std;
constexpr int64_t MOD=998244353,ROOT=3;
int64_t power(int64_t a,int64_t n){int64_t r=1;while(n){if(n&1)r=r*a%MOD;a=a*a%MOD;n>>=1;}return r;}
void ntt(vector<int64_t>&a,bool invert){
 size_t n=a.size();assert(n<=(1u<<23) && (MOD-1)%n==0);
 for(size_t i=1,j=0;i<n;i++){size_t bit=n>>1;for(;j&bit;bit>>=1)j^=bit;j^=bit;if(i<j)swap(a[i],a[j]);}
 for(size_t len=2;len<=n;len<<=1){int64_t wlen=power(ROOT,(MOD-1)/len);if(invert)wlen=power(wlen,MOD-2);
  for(size_t i=0;i<n;i+=len){int64_t w=1;for(size_t j=0;j<len/2;j++){
   int64_t u=a[i+j],v=a[i+j+len/2]*w%MOD;a[i+j]=u+v;if(a[i+j]>=MOD)a[i+j]-=MOD;
   a[i+j+len/2]=u-v;if(a[i+j+len/2]<0)a[i+j+len/2]+=MOD;w=w*wlen%MOD;
  }}
 }
 if(invert){int64_t ni=power(n,MOD-2);for(auto &x:a)x=x*ni%MOD;}
}
int main(int argc,char**argv){
 assert(power(ROOT,(MOD-1)/2)!=1 && power(ROOT,(MOD-1)/7)!=1 && power(ROOT,(MOD-1)/17)!=1);
 for(int64_t d=2;d*d<=MOD;d++)assert(MOD%d);
 assert(argc==2);ifstream in(argv[1],ios::binary);assert(in);uint32_t p,m;in.read((char*)&p,4);in.read((char*)&m,4);
 assert(in && p%4==1 && p>5);for(uint32_t d=2;uint64_t(d)*d<=p;d++)assert(p%d);
 vector<int>ch(p,-1);ch[0]=0;for(uint64_t x=1;x<=(p-1)/2;x++)ch[x*x%p]=1;
 vector<int32_t>C(m),z(m);int64_t n=0,sm=0,absum=0;uint32_t expected=0;
 for(uint32_t x=0;x<p;x++)if(ch[x]==1 && ch[(x+p-1)%p]==1)expected++;
 assert(expected==m);
 for(uint32_t i=0;i<m;i++){in.read((char*)&C[i],4);in.read((char*)&z[i],4);assert(in);assert(C[i]>=0 && C[i]<(int32_t)p && z[i]>=-1024 && z[i]<=1024);
  assert(!i || C[i]>C[i-1]);assert(ch[C[i]]==1 && ch[(C[i]+p-1)%p]==1);
  n+=int64_t(z[i])*z[i];sm+=z[i];absum+=abs(z[i]);
 }
 assert(in.peek()==EOF && n>0 && 2*absum<MOD);
 size_t L=1;while(L<2*p-1)L<<=1;assert(L<=(1u<<23));
 vector<int64_t>a(L),b(L);for(uint32_t x=0;x<p;x++)a[x]=(ch[x]+MOD)%MOD;for(uint32_t i=0;i<m;i++)b[C[i]]=(z[i]+MOD)%MOD;
 ntt(a,false);ntt(b,false);for(size_t i=0;i<L;i++)a[i]=a[i]*b[i]%MOD;ntt(a,true);
 int64_t A=0;uint32_t direct=0;
 for(uint32_t i=0;i<m;i++){
  uint32_t x=C[i];int64_t conv=a[x];if(x+p<2*p-1)conv=(conv+a[x+p])%MOD;if(conv>MOD/2)conv-=MOD;
  assert(abs(conv)<=absum);A+=int64_t(z[i])*conv;
  if(p<=4001){int64_t naive=0;for(uint32_t j=0;j<m;j++)naive+=int64_t(ch[(x+p-C[j])%p])*z[j];assert(naive==conv);direct++;}
 }
 int64_t maxcorr=-1;uint32_t maxt=0,largecount=0;
 if(2*n<MOD){
  fill(a.begin(),a.end(),0);fill(b.begin(),b.end(),0);
  for(uint32_t i=0;i<m;i++){a[C[i]]=(z[i]+MOD)%MOD;b[(p-C[i])%p]=(z[i]+MOD)%MOD;}
  ntt(a,false);ntt(b,false);for(size_t i=0;i<L;i++)a[i]=a[i]*b[i]%MOD;ntt(a,true);
  for(uint32_t t=0;t<p;t++){int64_t corr=a[t];if(t+p<2*p-1)corr=(corr+a[t+p])%MOD;if(corr>MOD/2)corr-=MOD;
   assert(abs(corr)<=n);if(!t)assert(corr==n);
   if(t && abs(corr)>maxcorr){maxcorr=abs(corr);maxt=t;}if(5*abs(corr)>=3*n)largecount++;
  }
 }
 long double h=(long double)A/(sqrt((long double)p)*n)-(long double)sm*sm/((long double)p*n);
 cout<<"{\"p\":"<<p<<",\"m\":"<<m<<",\"norm_squared\":"<<n<<",\"sum_entries\":"<<sm<<",\"sum_absolute_entries\":"<<absum<<",\"signed_quadratic_form\":"<<A<<",\"ntt_modulus\":"<<MOD<<",\"ntt_length\":"<<L<<",\"integer_recovery_margin\":"<<MOD-2*absum<<",\"direct_small_rows_checked\":"<<direct<<",\"autocorrelation_computed\":"<<(2*n<MOD?"true":"false")<<",\"max_abs_nonzero_translation_numerator\":"<<maxcorr<<",\"max_translation_shift\":"<<maxt<<",\"count_translation_ge_3_5\":"<<largecount<<",\"rayleigh_H_display\":"<<setprecision(20)<<h<<"}\n";
}
