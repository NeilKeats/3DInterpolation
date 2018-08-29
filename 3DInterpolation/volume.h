#pragma once

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>


int ReadMatrixSizeFromStream(FILE * file, int * m, int * n, int *p);

int WriteMatrixHeaderToStream(FILE * file, const int M, const int N, const int P);


template <class real>
int ReadMatrixFromStream(FILE * file, const  int M, const int N, const int P, real * Voldata)
{
	unsigned int readBytes;
	if ((readBytes = fread(Voldata, sizeof(real), M*N*P, file)) < (unsigned int)M*N*P)
	{
		perror("Reading unmatched size of data.\n");
		return -1;
	}

	return 0;
}


template <class real>
int WriteMatrixToStream(FILE * file, const int M, const int N, const int P, const real * Voldata) {

	unsigned int readBytes;
	if ((readBytes = fwrite(Voldata, sizeof(real), M*N*P, file)) < (unsigned int)M*N*P)
	{
		perror("Writing unmatched size of data.\n");
		return -1;
	}

	return 0;
};



template <class real>
class Volume{

public:
    Volume(){
        VolWidth = VolHeight = VolDepth = 0 ;
        VolData = nullptr;
    }

    ~Volume(){
        if (VolData != nullptr )
            free(VolData);
        
    }

    int VolWidth;
    int VolHeight;
    int VolDepth;

	//
    int Init(const int Width, const int Height, const int Depth) {

		if (Width<1 || Height<1 || Depth<1) {
			perror("Error for illegal size of volume");
			return -1;
		}

		VolWidth = Width;
		VolHeight = Height;
		VolDepth = Depth;

		if (VolData != nullptr)
			free(VolData);

		VolData = (real *)malloc(sizeof(real) * VolWidth * VolHeight * VolDepth);

		return 0;
	};

    //
	int ReadFromDisk(const char* filename) {

		FILE *file = fopen(filename, "rb");
		if (!file) {
			perror("File opening failed \n");
			return -1;
		}

		if (ReadMatrixSizeFromStream(file, &VolWidth, &VolHeight, &VolDepth) != 0)
		{
			perror("Error reading matrix header from disk file.\n");
			return -1;
		}

		if (Init(VolWidth, VolHeight, VolDepth) != 0) {
			perror("");
			return -1;
		}

		if (ReadMatrixFromStream(file, VolWidth, VolHeight, VolDepth, VolData) != 0) {
			perror("Error reading matrix data from disk file.\n");
			return -1;
		}

		fclose(file);

		return 0;

	};

	//
    int WriteToDisk(const char* filename) {

		FILE *file = fopen(filename, "wb");
		if (!file) {
			perror("File opening failed \n");
			return -1;
		}

		if (WriteMatrixHeaderToStream(file, VolWidth, VolHeight, VolDepth) != 0)
		{
			perror("Error writing matrix header to disk file.\n");
			return -1;
		}

		if (WriteMatrixToStream(file, VolWidth, VolHeight, VolDepth, VolData) != 0) {
			perror("Error writing matrix data to disk file.\n");
			return -1;
		}

		fclose(file);

		return 0;

	};
    

    real* VolData = nullptr;


};

/*

template <class real>
int ReadMatrixFromStream(FILE * file, const int M, const int N, const int P, real * Voldata);

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
*/


//Definition
template <class real>
int transpose_xy(real *Vdata, real *VTdata, const int Width, const int Height, const int Depth) {

	for (int z = 0; z<Depth; ++z) {
		for (int y = 0; y<Height; ++y) {
			int base = z*(Width*Height) + y*Width;
			for (int x = 0; x<Width; ++x) {
				int trans_id = z*(Width*Height) + x*Height + y;
				VTdata[trans_id] = Vdata[base + x];
			}
		}
	}

	return 0;
}

template <class real>
int transpose_xz(real *Vdata, real *VTdata, const int Width, const int Height, const int Depth) {

	for (int z = 0; z<Depth; ++z) {
		for (int y = 0; y<Height; ++y) {
			for (int x = 0; x<Width; ++x) {
				int ori_id = z*(Width*Height) + y*Width + x;
				int trans_id = x*(Height*Depth) + y*Depth + z;
				VTdata[trans_id] = Vdata[ori_id];
			}
		}
	}

	return 0;
}

template <class real>
int transpose_vol_xy(Volume<real> *VData, Volume<real> *VTData) {

	if (VData->VolWidth * VData->VolHeight * VData->VolDepth != VTData->VolWidth * VTData->VolHeight * VTData->VolDepth) {
		perror("Failed in transposing volume because unequal size.\n");
		return -1;
	}

	transpose_xy(VData->VolData, VTData->VolData, VData->VolWidth, VData->VolHeight, VData->VolDepth);

	VTData->VolDepth = VData->VolDepth;
	VTData->VolHeight = VData->VolWidth;
	VTData->VolWidth = VData->VolHeight;

	return 0;
}

template <class real>
int transpose_vol_xz(Volume<real> *VData, Volume<real> *VTData) {

	if (VData->VolWidth * VData->VolHeight * VData->VolDepth != VTData->VolWidth * VTData->VolHeight * VTData->VolDepth) {
		perror("Failed in transposing volume because unequal size.\n");
		return -1;
	}

	transpose_xz(VData->VolData, VTData->VolData, VData->VolWidth, VData->VolHeight, VData->VolDepth);

	VTData->VolDepth = VData->VolWidth;
	VTData->VolHeight = VData->VolHeight;
	VTData->VolWidth = VData->VolDepth;

	return 0;

}
