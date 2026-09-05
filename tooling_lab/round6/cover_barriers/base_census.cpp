// Exact scalar and vector-RS maximum agreement census; n<=32.
// Generate by barycentric evaluation; independently replay by Newton differences.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using U=std::uint64_t;
U p;
U add(U a,U b){return (a+b)%p;}
U sub(U a,U b){return a>=b?a-b:a+p-b;}
U mul(U a,U b){return a*b%p;}
U power(U a,U b){U y=1;while(b){if(b&1)y=mul(y,a);a=mul(a,a);b>>=1;}return y;}
bool prime(U v){if(v<2)return false;if(v%2==0)return v==2;for(U d=3;d<=v/d;d+=2)if(v%d==0)return false;return true;}
void write32(std::ostream& out,std::uint32_t v){char a[4];for(int j=0;j<4;j++)a[j]=char((v>>(8*j))&255);out.write(a,4);}
std::uint32_t read32(std::istream& in){unsigned char a[4];in.read(reinterpret_cast<char*>(a),4);if(!in)throw std::runtime_error("short ledger");std::uint32_t v=0;for(int j=0;j<4;j++)v|=std::uint32_t(a[j])<<(8*j);return v;}
struct Sample{std::array<std::vector<U>,2> word;std::array<int,3> maximum{0,0,0};std::array<std::vector<int>,3> best_base;std::array<std::uint32_t,3> best_support{};};
int main(int argc,char**argv){try{
 if(argc!=5)throw std::runtime_error("usage: base_census generate|verify input.txt ledger.bin summary.json");
 bool verify=std::string(argv[1])=="verify";if(!verify&&std::string(argv[1])!="generate")throw std::runtime_error("mode");
 std::ifstream input(argv[2]);int n,k,m;input>>p>>n>>k>>m;
 if(!input||n<1||n>32||k<1||k>n||m<1||p<3||p>2147483713ULL)throw std::runtime_error("unsupported input");
 if(!prime(p))throw std::runtime_error("prime characteristic required");
 std::vector<U> dom(n);for(auto&x:dom){input>>x;if(x>=p)throw std::runtime_error("domain");}
 auto unique=dom;std::sort(unique.begin(),unique.end());if(std::unique(unique.begin(),unique.end())!=unique.end())throw std::runtime_error("repeated domain");
 std::vector<Sample> samples(m);for(auto&s:samples)for(auto&w:s.word){w.resize(n);for(auto&x:w){input>>x;if(x>=p)throw std::runtime_error("word");}}
 if(!input)throw std::runtime_error("short input");
 std::vector<std::vector<U>> dif(n,std::vector<U>(n)),inv(n,std::vector<U>(n));
 for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(i!=j){dif[i][j]=sub(dom[i],dom[j]);inv[i][j]=power(dif[i][j],p-2);if(mul(dif[i][j],inv[i][j])!=1)throw std::runtime_error("not a field inverse");}
 std::fstream ledger(argv[3],std::ios::binary|(verify?std::ios::in:std::ios::out|std::ios::trunc));if(!ledger)throw std::runtime_error("ledger open");
 std::vector<int> base(k);for(int i=0;i<k;i++)base[i]=i;U count=0;
 do{
  std::uint32_t base_mask=0;for(int b:base)base_mask|=std::uint32_t(1)<<b;
  std::vector<std::array<std::uint32_t,2>> masks(m,{base_mask,base_mask});
  if(!verify){
   std::vector<U> weights(k,1);for(int j=0;j<k;j++)for(int ell=0;ell<k;ell++)if(j!=ell)weights[j]=mul(weights[j],inv[base[j]][base[ell]]);
   for(int i=0;i<n;i++)if(!(base_mask&(std::uint32_t(1)<<i))){
    U prod=1;for(int b:base)prod=mul(prod,dif[i][b]);
    std::vector<U> lambda(k);for(int j=0;j<k;j++)lambda[j]=mul(mul(prod,weights[j]),inv[i][base[j]]);
    for(int a=0;a<m;a++)for(int c=0;c<2;c++){
     U y=0;for(int j=0;j<k;j++)y=add(y,mul(lambda[j],samples[a].word[c][base[j]]));
     if(y==samples[a].word[c][i])masks[a][c]|=std::uint32_t(1)<<i;
    }
   }
  }else{
   for(int a=0;a<m;a++)for(int c=0;c<2;c++){
    std::vector<U> dd(k);for(int j=0;j<k;j++)dd[j]=samples[a].word[c][base[j]];
    for(int step=1;step<k;step++)for(int j=k-1;j>=step;j--)dd[j]=mul(sub(dd[j],dd[j-1]),inv[base[j]][base[j-step]]);
    for(int i=0;i<n;i++)if(!(base_mask&(std::uint32_t(1)<<i))){
     U y=dd.back();for(int j=k-2;j>=0;j--)y=add(mul(y,sub(dom[i],dom[base[j]])),dd[j]);
     if(y==samples[a].word[c][i])masks[a][c]|=std::uint32_t(1)<<i;
    }
   }
  }
  for(int a=0;a<m;a++){
   std::array<std::uint32_t,3> all{masks[a][0],masks[a][1],masks[a][0]&masks[a][1]};
   for(int c=0;c<3;c++){
    if(verify){if(read32(ledger)!=all[c])throw std::runtime_error("Newton support mismatch at base "+std::to_string(count));}
    else write32(ledger,all[c]);
    int size=__builtin_popcount(all[c]);if(size>samples[a].maximum[c]){samples[a].maximum[c]=size;samples[a].best_base[c]=base;samples[a].best_support[c]=all[c];}
   }
  }
  count++;if(count>1000000)throw std::runtime_error("million-base budget exceeded");
  int j=k-1;while(j>=0&&base[j]==n-k+j)j--;if(j<0)break;base[j]++;for(int ell=j+1;ell<k;ell++)base[ell]=base[ell-1]+1;
 }while(true);
 if(verify&&ledger.peek()!=std::char_traits<char>::eof())throw std::runtime_error("trailing ledger bytes");
 std::ofstream out(argv[4]);out<<"{\"status\":\"passed\",\"method\":\""<<(verify?"Newton replay":"barycentric enumeration")<<"\",\"basis_count\":"<<count<<",\"samples\":[";
 for(int a=0;a<m;a++){
  if(a)out<<',';out<<"{\"maxima_u0_u1_joint\":[";for(int c=0;c<3;c++){if(c)out<<',';out<<samples[a].maximum[c];}out<<"],\"witness_bases\":[";
  for(int c=0;c<3;c++){if(c)out<<',';out<<'[';for(int j=0;j<k;j++){if(j)out<<',';out<<samples[a].best_base[c][j];}out<<']';}out<<"],\"witness_support_masks\":[";
  for(int c=0;c<3;c++){if(c)out<<',';out<<samples[a].best_support[c];}out<<"]}";
 }out<<"]}\n";std::cout<<"passed "<<count<<" bases, "<<m<<" samples, "<<(verify?"Newton verification":"barycentric enumeration")<<"\n";
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
