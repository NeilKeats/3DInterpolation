#include "volume.h"
#include "spline.h"

#define SCALE_RATE 2.0

template <class real>
int scale(Volume<real> *vc, Volume<real> *vo){
    int width = int(vc->VolWidth*SCALE_RATE);
    int height = int(vc->VolHeight*SCALE_RATE);
    int depth = int(vc->VolDepth*SCALE_RATE);
    vo->Init(width, height, depth);

    real x,y,z;

#pragma omp parallel for schedule(dynamic) num_threads(4)
    for(int k=0; k<depth; ++k){
        z = real(k)/real(depth-1) * real(vc->VolDepth-1);
        for(int j=0; j<height; ++j){
            y = real(j)/real(height-1) * real(vc->VolHeight-1);
            for(int i=0; i<width; ++i){
                x = real(i)/real(width-1) * real(vc->VolWidth-1);
                vo->VolData[k*(width*height) + j*width + i] = interpolation(vc->VolData, vc->VolWidth, vc->VolHeight, vc->VolDepth, x, y, z);
            }
        }
    }

	return 0;
}

int main(){

    Volume<float> *vr , *vc , *vo= nullptr;    
    vr = new Volume<float>;
    vc = new Volume<float>;
    vo = new Volume<float>;

    vr->ReadFromDisk("D:\\Projects\\SIFT&DVC\\SIFT_DATA\\Voxel_512_ori_GN.bin");
    vc->Init(vr->VolWidth,vr->VolHeight,vr->VolDepth);
    Prefilter(vr,vc);
    scale(vc,vo);
    vo->WriteToDisk("D:\\Projects\\SIFT&DVC\\SIFT_DATA\\Voxel_512_ori_GN_INTE.bin");

    delete vr;
    delete vc;
    delete vo;

    return 0;
};