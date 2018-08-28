#pragma once

#include "volume.h"
#include <math.h>

float BSplinePreFilter[8] = {
	1.732176555412859f,  //b0
	-0.464135309171000f, //b1
	0.124364681271139f,
	-0.033323415913556f,
	0.008928982383084f,
	-0.002392513618779f,
	0.000641072092032f,
	-0.000171774749350f,  //b7
};

// w0, w1, w2, and w3 are the four cubic B-spline basis functions
template <class real>
inline real w0_c(real a)
{
	//    return (1.0f/6.0f)*(-a*a*a + 3.0f*a*a - 3.0f*a + 1.0f);
	return (1.0f / 6.0f)*(a*(a*(-a + 3.0f) - 3.0f) + 1.0f);   // optimized
}

template <class real>
inline real w1_c(real a)
{
	//    return (1.0f/6.0f)*(3.0f*a*a*a - 6.0f*a*a + 4.0f);
	return (1.0f / 6.0f)*(a*a*(3.0f*a - 6.0f) + 4.0f);
}

template <class real>
inline real w2_c(real a)
{
	//    return (1.0f/6.0f)*(-3.0f*a*a*a + 3.0f*a*a + 3.0f*a + 1.0f);
	return (1.0f / 6.0f)*(a*(a*(-3.0f*a + 3.0f) + 3.0f) + 1.0f);
}

template <class real>
inline real w3_c(real a)
{
	return (1.0f / 6.0f)*(a*a*a);
}

template <class real> 
int FIR_1D(real *Vinput, real *Voutput, const int length, const int num);

template <class real>
int Prefilter(Volume<real> *VData, Volume<real> *VCoeffi);

template <class real>
real interpolation(real *VCoeffi, const int Width ,const int Height, const int Depth, const real x, const real y, const real z);