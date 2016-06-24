// Beispielfunktion


#include <R.h>
// #include <Rmath.h>
// #include <R_ext/Lapack.h>
// #include <R_ext/BLAS.h>


extern "C" {
void pseudoflood(double *bild, int *nrow, int *ncol, double *oldColor, double *newColor)
{

	int nrows = *nrow;
	int ncols = *ncol;
	
	for(int i = 0; i < ncols; i++)
		for(int j = 0; j < nrows; j++)
			if(bild[i * nrows + j] == *oldColor){
				bild[i * nrows + j] = *newColor; 
			}
			    
}

}

/*

extern "C" {
	
void convolve(double *a, int *na, double *b, int *nb, double *ab)
{
    int nab = *na + *nb - 1;

    for(int i = 0; i < nab; i++)
        ab[i] = 0.0;
    for(int i = 0; i < *na; i++)
        for(int j = 0; j < *nb; j++)
            ab[i + j] += a[i] * b[j];
}

}

*/