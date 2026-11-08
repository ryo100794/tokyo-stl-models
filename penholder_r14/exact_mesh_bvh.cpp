// Exhaustive self-intersection of an indexed mesh with exact binary-coordinate predicates.
#include "exact_predicate.hpp"
#include <cstring>
#include <chrono>
#include <atomic>
#include <mutex>
#include <omp.h>
struct D {double x[3];};
struct Box {double lo[3],hi[3];};
struct Face {uint32_t v[3];};
struct Node {Box box;uint32_t begin,count,left=0,right=0;bool leaf()const{return left==0;}};
struct RawNum {uint64_t mant;int exp;bool neg;};
struct Hit {uint32_t a,b;uint8_t code;};
std::vector<V> verts;
std::vector<D> coords;
std::vector<Face> faces;
std::vector<Box> boxes;
std::vector<uint32_t> perm;
std::vector<Node> nodes;
std::vector<uint32_t> degenerate_faces;
std::atomic<uint64_t> processed{0};std::mutex printmutex;
Box emptybox(){Box b;for(int k=0;k<3;k++){b.lo[k]=std::numeric_limits<double>::infinity();b.hi[k]=-b.lo[k];}return b;}
void addbox(Box&a,const Box&b){for(int k=0;k<3;k++){a.lo[k]=std::min(a.lo[k],b.lo[k]);a.hi[k]=std::max(a.hi[k],b.hi[k]);}}
bool overlap(const Box&a,const Box&b){for(int k=0;k<3;k++)if(a.lo[k]>b.hi[k]||b.lo[k]>a.hi[k])return false;return true;}
uint32_t build(uint32_t begin,uint32_t count){
    Node n;n.begin=begin;n.count=count;n.box=emptybox();Box centers=emptybox();
    for(uint32_t j=begin;j<begin+count;j++){const Box&b=boxes[perm[j]];addbox(n.box,b);for(int k=0;k<3;k++){double c=(b.lo[k]+b.hi[k])*.5;centers.lo[k]=std::min(centers.lo[k],c);centers.hi[k]=std::max(centers.hi[k],c);}}
    uint32_t idx=nodes.size();nodes.push_back(n);
    if(count<=8)return idx;
    int axis=0;for(int k=1;k<3;k++)if(centers.hi[k]-centers.lo[k]>centers.hi[axis]-centers.lo[axis])axis=k;
    uint32_t half=count/2;
    std::nth_element(perm.begin()+begin,perm.begin()+begin+half,perm.begin()+begin+count,[axis](uint32_t a,uint32_t b){return boxes[a].lo[axis]+boxes[a].hi[axis]<boxes[b].lo[axis]+boxes[b].hi[axis];});
    uint32_t left=build(begin,half),right=build(begin+half,count-half);nodes[idx].left=left;nodes[idx].right=right;return idx;
}
uint8_t check(uint32_t ia,uint32_t ib){
    const Face&a=faces[ia];const Face&b=faces[ib];
    int common=0,other=-1;
    for(int j=0;j<3;j++){bool same=false;for(int i=0;i<3;i++)if(a.v[i]==b.v[j])same=true;if(same)common++;else other=j;}
    // For two distinct triangles sharing an exact edge, different planes can only meet on that edge.
    if(common==2){V n=cross(sub(verts[a.v[1]],verts[a.v[0]]),sub(verts[a.v[2]],verts[a.v[0]]));if(dot(n,sub(verts[b.v[other]],verts[a.v[0]]))!=0)return 1;}
    Tri ta{{verts[a.v[0]],verts[a.v[1]],verts[a.v[2]]}},tb{{verts[b.v[0]],verts[b.v[1]],verts[b.v[2]]}};
    return intersect(ta,tb);
}
struct Acc {std::array<uint64_t,7> count{};std::vector<Hit> hits;uint64_t pending=0;};
void paircheck(uint32_t a,uint32_t b,Acc&acc){
    if(!overlap(boxes[a],boxes[b]))return;
    uint8_t c=check(a,b);acc.count[c]++;if(c>=2)acc.hits.push_back(Hit{std::min(a,b),std::max(a,b),c});
    if(++acc.pending==100000){uint64_t now=processed.fetch_add(acc.pending)+acc.pending;acc.pending=0;if(now%1000000==0){std::lock_guard<std::mutex>guard(printmutex);std::cerr<<"exact_pairs_processed="<<now<<"\n";}}
}
void visit(uint32_t ai,uint32_t bi,Acc&acc){
    const Node&a=nodes[ai];const Node&b=nodes[bi];if(!overlap(a.box,b.box))return;
    if(ai==bi){
        if(a.leaf()){for(uint32_t i=0;i<a.count;i++)for(uint32_t j=i+1;j<a.count;j++)paircheck(perm[a.begin+i],perm[a.begin+j],acc);}
        else{visit(a.left,a.left,acc);visit(a.left,a.right,acc);visit(a.right,a.right,acc);}return;
    }
    if(a.leaf()&&b.leaf()){for(uint32_t i=0;i<a.count;i++)for(uint32_t j=0;j<b.count;j++)paircheck(perm[a.begin+i],perm[b.begin+j],acc);return;}
    if(!a.leaf()&&(b.leaf()||a.count>=b.count)){visit(a.left,bi,acc);visit(a.right,bi,acc);}else{visit(ai,b.left,acc);visit(ai,b.right,acc);}
}
void frontier(uint32_t ai,uint32_t bi,int depth,std::vector<std::pair<uint32_t,uint32_t>>&tasks){
    const Node&a=nodes[ai];const Node&b=nodes[bi];if(!overlap(a.box,b.box))return;
    if(depth==0||(a.leaf()&&b.leaf())){tasks.emplace_back(ai,bi);return;}
    if(ai==bi){frontier(a.left,a.left,depth-1,tasks);frontier(a.left,a.right,depth-1,tasks);frontier(a.right,a.right,depth-1,tasks);}
    else if(!a.leaf()&&(b.leaf()||a.count>=b.count)){frontier(a.left,bi,depth-1,tasks);frontier(a.right,bi,depth-1,tasks);}
    else{frontier(ai,b.left,depth-1,tasks);frontier(ai,b.right,depth-1,tasks);}
}
int main(int argc,char**argv){
    if(argc<3){std::cerr<<"usage: exact_mesh_bvh indexed_input.bin bad_pairs.bin [threads]\n";return 2;}
    auto t0=std::chrono::steady_clock::now();int threads=argc>3?std::max(1,atoi(argv[3])):3;
    std::ifstream in(argv[1],std::ios::binary);uint64_t nv,nf;in.read((char*)&nv,8);in.read((char*)&nf,8);if(!in||nv>100000000||nf>100000000)return 2;
    coords.resize(nv);std::vector<RawNum> raw(nv*3);int minexp=std::numeric_limits<int>::max();
    for(uint64_t i=0;i<nv*3;i++){uint64_t bits;in.read((char*)&bits,8);std::memcpy(&coords[i/3].x[i%3],&bits,8);uint64_t e=(bits>>52)&2047,f=bits&((uint64_t(1)<<52)-1);if(e==2047)return 3;
        auto&r=raw[i];r.neg=bits>>63;r.mant=e?(f|(uint64_t(1)<<52)):f;r.exp=e?int(e)-1075:-1074;if(r.mant){while((r.mant&1)==0){r.mant>>=1;r.exp++;}minexp=std::min(minexp,r.exp);}}
    if(minexp==std::numeric_limits<int>::max())minexp=0;verts.resize(nv);
    for(uint64_t i=0;i<nv*3;i++){auto&r=raw[i];Z z=r.mant;if(r.mant)z<<=(r.exp-minexp);if(r.neg)z=-z;verts[i/3].x[i%3]=std::move(z);}raw.clear();raw.shrink_to_fit();
    faces.resize(nf);in.read((char*)faces.data(),nf*sizeof(Face));if(!in)return 4;boxes.resize(nf);perm.resize(nf);
    for(uint64_t i=0;i<nf;i++){Box b=emptybox();for(int j=0;j<3;j++){if(faces[i].v[j]>=nv)return 5;for(int k=0;k<3;k++){double x=coords[faces[i].v[j]].x[k];b.lo[k]=std::min(b.lo[k],x);b.hi[k]=std::max(b.hi[k],x);}}boxes[i]=b;perm[i]=i;
        if(zero(cross(sub(verts[faces[i].v[1]],verts[faces[i].v[0]]),sub(verts[faces[i].v[2]],verts[faces[i].v[0]]))))degenerate_faces.push_back(i);
    }
    nodes.reserve(nf/3+1);build(0,nf);std::vector<std::pair<uint32_t,uint32_t>>tasks;frontier(0,0,7,tasks);
    std::cerr<<"bvh_nodes="<<nodes.size()<<" tasks="<<tasks.size()<<" threads="<<threads<<"\n";
    std::vector<Acc> accum(threads);
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for(size_t i=0;i<tasks.size();i++){int tid=omp_get_thread_num();visit(tasks[i].first,tasks[i].second,accum[tid]);}
    std::array<uint64_t,7> counts{};std::vector<Hit> hits;
    for(auto&a:accum){for(int k=0;k<7;k++)counts[k]+=a.count[k];hits.insert(hits.end(),a.hits.begin(),a.hits.end());}
    std::sort(hits.begin(),hits.end(),[](const Hit&a,const Hit&b){return a.a<b.a||(a.a==b.a&&a.b<b.b);});
    std::ofstream out(argv[2],std::ios::binary);for(const auto&h:hits){out.write((char*)&h.a,4);out.write((char*)&h.b,4);out.write((char*)&h.code,1);}
    std::ofstream degout(std::string(argv[2])+".degenerate.bin",std::ios::binary);degout.write((char*)degenerate_faces.data(),4*degenerate_faces.size());
    uint64_t total=0;for(auto c:counts)total+=c;double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    std::cout<<"{\"vertices\":"<<nv<<",\"faces\":"<<nf<<",\"aabb_candidate_pairs\":"<<total<<",\"common_power_of_two\":"<<minexp<<",\"counts\":[";
    for(int i=0;i<7;i++){if(i)std::cout<<",";std::cout<<counts[i];}
    std::cout<<"],\"abnormal_pairs\":"<<hits.size()<<",\"degenerate_faces\":"<<degenerate_faces.size()<<",\"threads\":"<<threads<<",\"elapsed_seconds\":"<<sec<<"}\n";
}
