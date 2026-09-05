#include <boost/multiprecision/cpp_int.hpp>
#include <array>
#include <chrono>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;
using Hist = std::vector<std::array<long long,4>>;

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

long long choose(int n, int k) {
    if (k<0 || k>n) return 0;
    long long out=1;
    for (int j=1;j<=k;++j) out=out*(n-j+1)/j;
    return out;
}

int sign_power(long long value, int exponent) {
    if (!exponent) return 1;
    if (!value) return 0;
    return exponent%2 ? int(value) : 1;
}

struct Algebra {
    int degree, limit, side, size;
    std::vector<std::array<int,3>> exponent;
    std::vector<std::vector<long long>> cover;
    std::vector<std::vector<std::pair<int,int>>> transitions;

    int index(int a,int b,int c) const { return (a*side+b)*side+c; }

    Algebra(int d,int k):degree(d),limit(k),side(d+1),size(side*side*side),
        exponent(size),cover(k+1,std::vector<long long>(size)),transitions(size) {
        for(int a=0;a<=d;++a) for(int b=0;b<=d;++b) for(int c=0;c<=d;++c) {
            int alpha=index(a,b,c);
            exponent[alpha]={a,b,c};
            for(int u=1;u<=k;++u) {
                if (u<a || u<b || u<c || u>a+b+c) continue;
                // Inclusion-exclusion: three subsets of a u-element set,
                // with respective sizes a,b,c, whose union is the whole set.
                cpp_int value=0;
                for(int h=0;h<=u;++h) {
                    cpp_int term=cpp_int(choose(u,h))*choose(h,a)*choose(h,b)*choose(h,c);
                    value += (u-h)%2 ? -term : term;
                }
                require(value>=0 && value<=std::numeric_limits<long long>::max(),"cover coefficient overflow");
                cover[u][alpha]=value.convert_to<long long>();
            }
            for(int x=0;x+a<=d;++x) for(int y=0;y+b<=d;++y) for(int z=0;z+c<=d;++z)
                transitions[alpha].push_back({index(x,y,z),index(a+x,b+y,c+z)});
        }
    }
};

void run(int case_id,const Algebra& alg,const Hist& hist) {
    auto started=std::chrono::steady_clock::now();
    int d=alg.degree,K=alg.limit,N=alg.size;
    std::vector<long long> power_sums(N);
    for(int alpha=0;alpha<N;++alpha) {
        auto exp=alg.exponent[alpha];
        for(auto row:hist)
            power_sums[alpha] += row[3]*sign_power(row[0],exp[0])*sign_power(row[1],exp[1])*sign_power(row[2],exp[2]);
    }
    // E_k is the k-th elementary polynomial in H_y,
    // H_y=(1+a_y*t1)(1+b_y*t2)(1+c_y*t3)-1.
    // Newton: k E_k = sum_{u=1}^k (-1)^(u-1) E_{k-u} sum_y H_y^u.
    std::vector<std::vector<cpp_int>> E(K+1,std::vector<cpp_int>(N));
    E[0][0]=1;
    unsigned long long products=0;
    std::size_t max_bits=1;
    for(int k=1;k<=K;++k) {
        for(int u=1;u<=k;++u) {
            const auto& previous=E[k-u];
            for(int alpha=1;alpha<N;++alpha) {
                if (!alg.cover[u][alpha] || !power_sums[alpha]) continue;
                cpp_int factor=cpp_int(alg.cover[u][alpha])*power_sums[alpha];
                if (u%2==0) factor=-factor;
                for(auto [beta,target]:alg.transitions[alpha]) {
                    if (previous[beta]==0) continue;
                    E[k][target] += factor*previous[beta];
                    ++products;
                }
            }
        }
        for(auto& value:E[k]) {
            require(value%k==0,"Newton division was not exact");
            value/=k;
            if (value!=0) {
                cpp_int positive=value<0 ? -value : value;
                max_bits=std::max(max_bits,boost::multiprecision::msb(positive)+1);
            }
        }
    }
    double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::cout << "{\"case_id\":"<<case_id<<",\"degree\":"<<d<<",\"max_union\":"<<K
              <<",\"integer_products\":"<<products<<",\"maximum_intermediate_bits\":"<<max_bits
              <<",\"backend_seconds\":"<<elapsed<<",\"coefficients\":[";
    for(int k=0;k<=K;++k) {
        if(k) std::cout<<',';
        std::cout<<E[k][alg.index(d,d,d)];
    }
    std::cout<<"]}\n";
}

int main() {
    try {
        int cases;
        require(bool(std::cin>>cases) && cases>=1 && cases<=10000,"invalid case count");
        std::map<std::pair<int,int>,std::shared_ptr<Algebra>> cached;
        for(int id=0;id<cases;++id) {
            int d,K,groups;
            require(bool(std::cin>>d>>K>>groups),"missing case header");
            require(d>=0 && d<=6 && K>=0 && K<=3*d && groups>=1 && groups<=27,"unsupported truncation or groups");
            Hist hist(groups);
            long long population=0;
            for(auto& row:hist) {
                require(bool(std::cin>>row[0]>>row[1]>>row[2]>>row[3]),"missing histogram row");
                require(row[0]>=-1 && row[0]<=1 && row[1]>=-1 && row[1]<=1 && row[2]>=-1 && row[2]<=1,"invalid sign");
                require(row[3]>=0 && row[3]<=1000000000000LL,"invalid multiplicity");
                population+=row[3];
            }
            require(population>0 && population<=1000000000000LL && K<=population,"invalid population");
            auto key=std::make_pair(d,K);
            if (!cached.count(key)) cached[key]=std::make_shared<Algebra>(d,K);
            run(id,*cached[key],hist);
        }
        std::string extra;
        require(!(std::cin>>extra),"trailing input");
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
