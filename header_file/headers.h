/*Header file for the source code*/

typedef struct{
    fftw_complex *c;
}composition;

extern void modelB(int NC, int N, double dr, int time, double dt, double kappa, double **chi, double lambda, composition *comp);
