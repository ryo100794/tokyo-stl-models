#pragma once
// Exact triangle-pair intersection for binary floating-point input coordinates.
// Every double is decoded as a signed integer times a common power of two.
// GMP arbitrary-precision integer signs and rational endpoint comparisons determine the result.
// Input: uint64 ntri,npairs; ntri*9 little-endian doubles; npairs*2 uint32 IDs.
// Output: one uint8 classification per pair.
#include <gmpxx.h>
#include <array>
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include <stdexcept>
using Z=mpz_class;
struct V {std::array<Z,3> x;};
V sub(const V&a,const V&b){V c;for(int k=0;k<3;k++)c.x[k]=a.x[k]-b.x[k];return c;}
V cross(const V&a,const V&b){return V{{a.x[1]*b.x[2]-a.x[2]*b.x[1],a.x[2]*b.x[0]-a.x[0]*b.x[2],a.x[0]*b.x[1]-a.x[1]*b.x[0]}};}
Z dot(const V&a,const V&b){return a.x[0]*b.x[0]+a.x[1]*b.x[1]+a.x[2]*b.x[2];}
int sign(const Z&z){return z==0?0:(z>0?1:-1);}
bool eq(const V&a,const V&b){return a.x==b.x;}
bool zero(const V&a){return a.x[0]==0&&a.x[1]==0&&a.x[2]==0;}
using Tri=std::array<V,3>;
struct P {V n; Z d=1;};
P point(const V&v){return P{v,Z(1)};}
P normalize(P p){if(p.d<0){p.d=-p.d;for(auto&x:p.n.x)x=-x;}return p;}
bool eqp(const P&a,const P&b){for(int k=0;k<3;k++)if(a.n.x[k]*b.d!=b.n.x[k]*a.d)return false;return true;}
int cmp(const P&a,const P&b,int k){return sign(a.n.x[k]*b.d-b.n.x[k]*a.d);}
P edge_plane(const V&a,const V&b,const Z&da,const Z&db){
    P p;p.d=da-db;for(int k=0;k<3;k++)p.n.x[k]=da*b.x[k]-db*a.x[k];return normalize(p);
}
std::vector<P> cut(const Tri&t,const std::array<Z,3>&d){
    std::vector<P> p;
    for(int i=0;i<3;i++){
        if(d[i]==0)p.push_back(point(t[i]));
        int j=(i+1)%3;
        if(sign(d[i])*sign(d[j])<0)p.push_back(edge_plane(t[i],t[j],d[i],d[j]));
    }
    return p;
}
bool strict_side(const std::array<Z,3>&d){return (d[0]>0&&d[1]>0&&d[2]>0)||(d[0]<0&&d[1]<0&&d[2]<0);}
bool spans(const std::array<Z,3>&d){bool pos=false,neg=false;for(const auto&x:d){pos|=x>0;neg|=x<0;}return pos&&neg;}
Z orient2(const V&a,const V&b,const V&c,int u,int v){return (b.x[u]-a.x[u])*(c.x[v]-a.x[v])-(b.x[v]-a.x[v])*(c.x[u]-a.x[u]);}
bool inside(const V&p,const Tri&t,int u,int v){
    int a=sign(orient2(t[0],t[1],p,u,v)),b=sign(orient2(t[1],t[2],p,u,v)),c=sign(orient2(t[2],t[0],p,u,v));
    return (a>=0&&b>=0&&c>=0)||(a<=0&&b<=0&&c<=0);
}
P cross2edges(const V&a,const V&b,const V&c,const V&d,int u,int v){
    V ab=sub(b,a),cd=sub(d,c),ca=sub(c,a);
    Z den=ab.x[u]*cd.x[v]-ab.x[v]*cd.x[u];
    Z num=ca.x[u]*cd.x[v]-ca.x[v]*cd.x[u];
    P p;p.d=den;for(int k=0;k<3;k++)p.n.x[k]=a.x[k]*den+ab.x[k]*num;return normalize(p);
}
bool on_shared_segment(const P&p,const V&a,const V&b,int u,int v){
    Z px=p.n.x[u]-a.x[u]*p.d,py=p.n.x[v]-a.x[v]*p.d;
    if((b.x[u]-a.x[u])*py-(b.x[v]-a.x[v])*px!=0)return false;
    int k=(a.x[u]!=b.x[u]?u:v);
    Z lo=std::min(a.x[k],b.x[k])*p.d,hi=std::max(a.x[k],b.x[k])*p.d;
    return p.n.x[k]>=lo&&p.n.x[k]<=hi;
}
int dimension(const std::vector<P>&p,int u,int v){
    if(p.empty())return -1;
    size_t j=1;while(j<p.size()&&eqp(p[0],p[j]))j++;
    if(j==p.size())return 0;
    Z ax=p[j].n.x[u]*p[0].d-p[0].n.x[u]*p[j].d;
    Z ay=p[j].n.x[v]*p[0].d-p[0].n.x[v]*p[j].d;
    for(size_t k=j+1;k<p.size();k++){
        Z bx=p[k].n.x[u]*p[0].d-p[0].n.x[u]*p[k].d;
        Z by=p[k].n.x[v]*p[0].d-p[0].n.x[v]*p[k].d;
        if(ax*by-ay*bx!=0)return 2;
    }
    return 1;
}
// 0 disjoint; 1 expected shared-simplex contact; 2 other point contact;
// 3 noncoplanar proper crossing segment; 4 coplanar area overlap;
// 5 other segment contact; 6 degenerate input triangle.
uint8_t intersect(const Tri&a,const Tri&b){
    V na=cross(sub(a[1],a[0]),sub(a[2],a[0]));
    V nb=cross(sub(b[1],b[0]),sub(b[2],b[0]));
    if(zero(na)||zero(nb))return 6;
    std::array<Z,3> da,db;
    for(int i=0;i<3;i++){da[i]=dot(nb,sub(a[i],b[0]));db[i]=dot(na,sub(b[i],a[0]));}
    if(strict_side(da)||strict_side(db))return 0;
    std::vector<V> shared;
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(eq(a[i],b[j]))shared.push_back(a[i]);
    bool cop=(db[0]==0&&db[1]==0&&db[2]==0);
    if(cop){
        int drop=0;while(drop<3&&na.x[drop]==0)drop++;
        int u=(drop+1)%3,v=(drop+2)%3;
        std::vector<P> p;
        for(int i=0;i<3;i++){if(inside(a[i],b,u,v))p.push_back(point(a[i]));if(inside(b[i],a,u,v))p.push_back(point(b[i]));}
        for(int i=0;i<3;i++)for(int j=0;j<3;j++){
            const V&a0=a[i],&a1=a[(i+1)%3],&b0=b[j],&b1=b[(j+1)%3];
            int o1=sign(orient2(a0,a1,b0,u,v)),o2=sign(orient2(a0,a1,b1,u,v)),o3=sign(orient2(b0,b1,a0,u,v)),o4=sign(orient2(b0,b1,a1,u,v));
            if(o1*o2<0&&o3*o4<0)p.push_back(cross2edges(a0,a1,b0,b1,u,v));
        }
        if(p.empty())return 0;
        bool allowed=false;
        if(shared.size()==1){allowed=true;P s=point(shared[0]);for(const P&q:p)if(!eqp(q,s)){allowed=false;break;}}
        else if(shared.size()==2){allowed=true;for(const P&q:p)if(!on_shared_segment(q,shared[0],shared[1],u,v)){allowed=false;break;}}
        if(allowed)return 1;
        int dim=dimension(p,u,v);return dim==2?4:(dim==1?5:2);
    }
    V direction=cross(na,nb);int k=0;while(k<3&&direction.x[k]==0)k++;
    if(k==3)return 0;
    auto pa=cut(a,da),pb=cut(b,db);if(pa.empty()||pb.empty())return 0;
    auto bounds=[k](const std::vector<P>&p){P lo=p[0],hi=p[0];for(const P&q:p){if(cmp(q,lo,k)<0)lo=q;if(cmp(q,hi,k)>0)hi=q;}return std::pair<P,P>(lo,hi);};
    auto ia=bounds(pa),ib=bounds(pb);
    P lo=cmp(ia.first,ib.first,k)>=0?ia.first:ib.first;
    P hi=cmp(ia.second,ib.second,k)<=0?ia.second:ib.second;
    int c=cmp(lo,hi,k);if(c>0)return 0;
    bool allowed=false;
    if(shared.size()==1){P s=point(shared[0]);allowed=eqp(lo,s)&&eqp(hi,s);}
    else if(shared.size()==2){
        P l=point(shared[0]),h=point(shared[1]);if(cmp(l,h,k)>0)std::swap(l,h);
        allowed=cmp(lo,l,k)>=0&&cmp(hi,h,k)<=0;
    }
    if(allowed)return 1;
    if(c==0)return 2;
    return spans(da)&&spans(db)?3:5;
}
