// Orthographic, opaque, depth-buffered rendering of the exact saved triangles.
// No graphics driver. Input contains screen-space vertices and camera normals.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <vector>

template<class T> void read(std::ifstream& in,T& x){in.read(reinterpret_cast<char*>(&x),sizeof(T));}
struct Vertex {float x,y,z,nx,ny,nz;};
struct Face {uint32_t a,b,c;};
float edge(float ax,float ay,float bx,float by,float px,float py){return (px-ax)*(by-ay)-(py-ay)*(bx-ax);}
int main(int argc,char**argv){
  if(argc!=3)return 2;
  std::ifstream in(argv[1],std::ios::binary);
  uint32_t width,height,nparts;read(in,width);read(in,height);read(in,nparts);
  std::array<uint8_t,3> bg;in.read(reinterpret_cast<char*>(bg.data()),3);
  size_t count=size_t(width)*height;
  std::vector<float> zbuf(count,-std::numeric_limits<float>::infinity());
  std::vector<uint8_t> rgb(count*3);
  for(size_t i=0;i<count;i++)for(int c=0;c<3;c++)rgb[3*i+c]=bg[c];
  float lx=.32f,ly=.46f,lz=.828f;
  float ll=std::sqrt(lx*lx+ly*ly+lz*lz);lx/=ll;ly/=ll;lz/=ll;
  float hx=lx,hy=ly,hz=lz+1,hl=std::sqrt(hx*hx+hy*hy+hz*hz);hx/=hl;hy/=hl;hz/=hl;
  for(uint32_t part=0;part<nparts;part++){
    uint32_t nv,nf;read(in,nv);read(in,nf);
    float color[3];in.read(reinterpret_cast<char*>(color),sizeof(color));
    std::vector<Vertex> vertices(nv);std::vector<Face> faces(nf);
    in.read(reinterpret_cast<char*>(vertices.data()),sizeof(Vertex)*nv);
    in.read(reinterpret_cast<char*>(faces.data()),sizeof(Face)*nf);
    if(!in)throw std::runtime_error("Incomplete mesh input");
    for(auto face:faces){
      const auto& a=vertices[face.a];const auto& b=vertices[face.b];const auto& c=vertices[face.c];
      float ar=edge(a.x,a.y,b.x,b.y,c.x,c.y);
      if(std::abs(ar)<1e-8f)continue;
      int x0=std::max(0,int(std::floor(std::min({a.x,b.x,c.x}))));
      int x1=std::min(int(width)-1,int(std::ceil(std::max({a.x,b.x,c.x}))));
      int y0=std::max(0,int(std::floor(std::min({a.y,b.y,c.y}))));
      int y1=std::min(int(height)-1,int(std::ceil(std::max({a.y,b.y,c.y}))));
      if(x0>x1||y0>y1)continue;
      for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++){
        float wa=edge(b.x,b.y,c.x,c.y,x+.5f,y+.5f)/ar;
        float wb=edge(c.x,c.y,a.x,a.y,x+.5f,y+.5f)/ar;
        float wc=1.f-wa-wb;
        if(wa<0||wb<0||wc<0)continue;
        float z=wa*a.z+wb*b.z+wc*c.z;
        size_t pos=size_t(y)*width+x;
        if(z<=zbuf[pos])continue;
        float nx=wa*a.nx+wb*b.nx+wc*c.nx;
        float ny=wa*a.ny+wb*b.ny+wc*c.ny;
        float nz=wa*a.nz+wb*b.nz+wc*c.nz;
        float nl=std::sqrt(nx*nx+ny*ny+nz*nz);
        if(nl>1e-20f){nx/=nl;ny/=nl;nz/=nl;}
        float diffuse=std::max(0.f,nx*lx+ny*ly+nz*lz);
        float specular=std::pow(std::max(0.f,nx*hx+ny*hy+nz*hz),30.f);
        float shade=.28f+.70f*diffuse;
        for(int k=0;k<3;k++)rgb[3*pos+k]=uint8_t(std::clamp(255.f*(color[k]*shade+.10f*specular),0.f,255.f));
        zbuf[pos]=z;
      }
    }
  }
  std::ofstream out(argv[2],std::ios::binary);out.write(reinterpret_cast<char*>(rgb.data()),rgb.size());
  return 0;
}
