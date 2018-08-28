#include "spline.h"
#include <omp.h>

#define MIN(x,y) ( ((x)<(y)) ? (x):(y) )
#define MAX(x,y) ( ((x)>(y)) ? (x):(y) )

#define ELT(height,width,x,y,z) ((((long)z)*((long)height)+((long)y))*((long)width)+(long)x)

template <class real> 
int FIR_1D(real *Vinput, real *Voutput, const int length, const int num){
#pragma omp parallel for schedule(dynamic) num_threads(19)
    for(int i=0; i<num; ++i){
        real *base = Vinput + i*length;
        for(int j=7; j<length-7; ++j){
            Voutput[i*length + j] = 
                BSplinePreFilter[7] * (base[j-7] + base[j+7]) +
                BSplinePreFilter[6] * (base[j-6] + base[j+6]) +
                BSplinePreFilter[5] * (base[j-5] + base[j+5]) +
                BSplinePreFilter[4] * (base[j-4] + base[j+4]) +
                BSplinePreFilter[3] * (base[j-3] + base[j+3]) +
                BSplinePreFilter[2] * (base[j-2] + base[j+2]) +
                BSplinePreFilter[1] * (base[j-1] + base[j+1]) +
                BSplinePreFilter[0] * (base[j]);
        }

        for(int j=0; j<7; ++j){
            Voutput[i*length + j] = 
                BSplinePreFilter[7] * (base[MAX(j-7,0)] + base[j+7]) +
                BSplinePreFilter[6] * (base[MAX(j-6,0)] + base[j+6]) +
                BSplinePreFilter[5] * (base[MAX(j-5,0)] + base[j+5]) +
                BSplinePreFilter[4] * (base[MAX(j-4,0)] + base[j+4]) +
                BSplinePreFilter[3] * (base[MAX(j-3,0)] + base[j+3]) +
                BSplinePreFilter[2] * (base[MAX(j-2,0)] + base[j+2]) +
                BSplinePreFilter[1] * (base[MAX(j-1,0)] + base[j+1]) +
                BSplinePreFilter[0] * (base[j]);
        }
        
        for(int j=length-7; j<length; ++j){
            Voutput[i*length + j] = 
                BSplinePreFilter[7] * (base[j-7] + base[MIN(j+7,length-1)]) +
                BSplinePreFilter[6] * (base[j-6] + base[MIN(j+6,length-1)]) +
                BSplinePreFilter[5] * (base[j-5] + base[MIN(j+5,length-1)]) +
                BSplinePreFilter[4] * (base[j-4] + base[MIN(j+4,length-1)]) +
                BSplinePreFilter[3] * (base[j-3] + base[MIN(j+3,length-1)]) +
                BSplinePreFilter[2] * (base[j-2] + base[MIN(j+2,length-1)]) +
                BSplinePreFilter[1] * (base[j-1] + base[MIN(j+1,length-1)]) +
                BSplinePreFilter[0] * (base[j]);

        }
    }
    return 0;
}


template <class real>
int Prefilter(Volume<real> *VData, Volume<real> *VCoeffi){

    Volume<real> Vtmp;
    Vtmp.Init(VData->VolWidth, VData->VolHeight, VData->VolDepth);

    //prefilter along x-axis direction
    FIR_1D(VData->VolData, VCoeffi->VolData, VData->VolWidth, VData->VolHeight*VData->VolDepth);
    
    
    //transpose xy
    transpose_vol_xy(VCoeffi, &Vtmp);
    //prefilter along y-axis direction(transposed)
    FIR_1D(Vtmp.VolData, VCoeffi->VolData, Vtmp.VolWidth, Vtmp.VolHeight*Vtmp.VolDepth);
    //transpose yx
    transpose_vol_xy(VCoeffi, &Vtmp);

    //transpose xz
    transpose_vol_xz(&Vtmp, VCoeffi);
    //prefilter along z-axis direction(transposed)
    FIR_1D(VCoeffi->VolData, Vtmp.VolData, VCoeffi->VolWidth, VCoeffi->VolHeight*VCoeffi->VolDepth);
    //transpose zx
    transpose_vol_xz(&Vtmp, VCoeffi);

    return 0;
}

template <class real>
inline real interpolation(real *VCoeffi, const int Width ,const int Height, const int Depth, const real x, const real y, const real z)
{
    int ix=floor();
    int iy=floor();
    int iz=floor();

    real fx = x - real(ix);
    real fy = y - real(iy);
    real fz = z - real(iz);
    
    real w_x[4];
	real w_y[4];
	real w_z[4];
	real sum_x[4];
	real sum_y[4];
    real r = 0;

	w_x[0] = w0_c(fx);
	w_x[1] = w1_c(fx);
	w_x[2] = w2_c(fx);
	w_x[3] = w3_c(fx);

	w_y[0] = w0_c(fy);
	w_y[1] = w1_c(fy);
	w_y[2] = w2_c(fy);
	w_y[3] = w3_c(fy);

    w_z[0] = w0_c(fz);
    w_z[1] = w1_c(fz);
    w_z[2] = w2_c(fz);
    w_z[3] = w3_c(fz);

    //z
    for(int j=0; j<4; ++j){
        //y
        for(int i=0; i<4; ++i){
            sum_x[i] = 
               w_x[0] * VCoeffi[ELT(Height,Width, MIN( MAX(ix-1,0),Width-1), MIN( MAX(iy-1+i,0),Height-1), MIN( MAX(iz-1+j,0),Depth-1))]
            +  w_x[1] * VCoeffi[ELT(Height,Width, MIN( MAX(ix,0),Width-1), MIN( MAX(iy-1+i,0),Height-1), MIN( MAX(iz-1+j,0),Depth-1))]
            +  w_x[2] * VCoeffi[ELT(Height,Width, MIN( MAX(ix+1,0),Width-1), MIN( MAX(iy-1+i,0),Height-1), MIN( MAX(iz-1+j,0),Depth-1))]
            +  w_x[3] * VCoeffi[ELT(Height,Width, MIN( MAX(ix+2,0),Width-1), MIN( MAX(iy-1+i,0),Height-1), MIN( MAX(iz-1+j,0),Depth-1))];

        }
        sum_y[j] = 		
          w_y[0] * sum_x[0] 
		+ w_y[1] * sum_x[1]
		+ w_y[2] * sum_x[2]
		+ w_y[3] * sum_x[3];
    }

    r = 
      w_z[0] * sum_y[0] 
    + w_z[1] * sum_y[1]
    + w_z[2] * sum_y[2]
    + w_z[3] * sum_y[3];

    return r;
}