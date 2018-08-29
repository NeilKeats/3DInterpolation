#include "volume.h"


int ReadMatrixSizeFromStream(FILE * file, int * m, int * n, int *p) {
	if (fread(m, sizeof(int), 1, file) < 1)
		return -1;
	if (fread(n, sizeof(int), 1, file) < 1)
		return -1;
	if (fread(p, sizeof(int), 1, file) < 1)
		return -1;

	return 0;
};


int WriteMatrixHeaderToStream(FILE * file, const int M, const int N, const int P) {

	if (fwrite(&M, sizeof(int), 1, file) < 1)
		return -1;
	if (fwrite(&N, sizeof(int), 1, file) < 1)
		return -1;
	if (fwrite(&P, sizeof(int), 1, file) < 1)
		return -1;

	return 0;
}