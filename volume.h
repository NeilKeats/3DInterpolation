#pragma once

#include <inttypes.h>
#include <stdio.h>

template <class real>
class Volume{

public:
    Volume(){
        VolWidth = VolHeight = VolDepth = 0 ;
        VolData = nullptr;
        VolCoeffi = nullptr;
    }

    ~Volume(){
        if (VolData != nullptr )
            free(VolData);
        if (VolCoeffi != nullptr )
            free(VolCoeffi);
        
    }

    int VolWidth;
    int VolHeight;
    int VolDepth;

    int Init(const int Width, const int Height, const int Depth);
    int ReadFromDisk(const char* filename);
    int WriteToDisk(const char* filename);
    

    real* VolData = nullptr;


};

int ReadMatrixSizeFromStream(FILE * file, int * m, int * n, int *p);

template <class real>
int ReadMatrixFromStream(FILE * file, const int M, const int N, const int P, real * Voldata);

int WriteMatrixHeaderToStream(FILE * file,  const int M, const int N, const int P);

template <class real>
int WriteMatrixToStream(FILE * file, const int M, const int N, const int P, const real * Voldata);

template <class real> 
int transpose_xy(real *Vdata, real *VTdata, const int Width, const int Height , const int Depth);

template <class real> 
int transpose_xz(real *Vdata, real *VTdata, const int Width, const int Height , const int Depth);

template <class real>
int transpose_vol_xy(Volume<real> *VData, Volume<real> *VTData);

template <class real>
int transpose_vol_xz(Volume<real> *VData, Volume<real> *VTData);