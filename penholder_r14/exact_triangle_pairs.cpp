#include "exact_predicate.hpp"
struct Raw {uint64_t mant;int exp;bool neg;};
int main(int argc,char**argv){
    if(argc!=3){std::cerr<<"usage: exact_triangle_pairs input.bin output.bin\n";return 2;}
    std::ifstream in(argv[1],std::ios::binary);uint64_t n=0,np=0;in.read((char*)&n,8);in.read((char*)&np,8);
    if(!in||n>10000000||np>1000000000){std::cerr<<"invalid input header\n";return 2;}
    std::vector<Raw> raw(n*9);int minexp=std::numeric_limits<int>::max();
    for(auto&r:raw){
        uint64_t bits;in.read((char*)&bits,8);uint64_t e=(bits>>52)&2047,f=bits&((uint64_t(1)<<52)-1);
        if(e==2047){std::cerr<<"nonfinite coordinate\n";return 2;}
        r.neg=bits>>63;r.mant=e?(f|(uint64_t(1)<<52)):f;r.exp=e?int(e)-1075:-1074;
        if(r.mant){while((r.mant&1)==0){r.mant>>=1;r.exp++;}minexp=std::min(minexp,r.exp);}
    }
    if(minexp==std::numeric_limits<int>::max())minexp=0;
    std::vector<Tri> ts(n);
    for(uint64_t i=0;i<n*9;i++){
        const Raw&r=raw[i];Z z=r.mant;if(r.mant)z<<=(r.exp-minexp);if(r.neg)z=-z;ts[i/9][(i%9)/3].x[i%3]=std::move(z);
    }
    std::ofstream out(argv[2],std::ios::binary);std::array<uint64_t,7> counts{};
    for(uint64_t i=0;i<np;i++){
        uint32_t p[2];in.read((char*)p,8);if(!in||p[0]>=n||p[1]>=n){std::cerr<<"bad pair\n";return 2;}
        uint8_t c=intersect(ts[p[0]],ts[p[1]]);counts[c]++;out.write((char*)&c,1);
    }
    std::cout<<"triangles="<<n<<" pairs="<<np<<" common_power_of_two="<<minexp<<" counts=";
    for(auto c:counts)std::cout<<c<<",";std::cout<<"\n";
    return 0;
}
