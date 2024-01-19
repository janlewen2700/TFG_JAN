
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "fftw3.h"
#include<gsl/gsl_rng.h>
#include "../header_file/headers.h"

int main(void){
    
    FILE *fr, *fw;
    int N, NC, time;
    double dr, dt, lambda, kappa, c_noise;
    double **chi;
    gsl_rng * ran_num;
    const gsl_rng_type * Taus;

    /*Clearing the output directory in order to print fresh result from every new simulation run.*/

    (void) system("rm -rf output/*");

    
    fw = fopen("output/simulation_data","w");
    fr = fopen("input/system_time","r");
    if (fw == NULL || fr == NULL){
        printf("Unable to open file.");
        exit(0);
    }
    /*let's take the time and space info*/
    fscanf(fr,"%d%d%le%d%le",&NC,&N,&dr,&time,&dt);
    fclose(fr);
    fflush(fr);
    fprintf(fw,"Number of components = %d\n Size grid = %d\n Grid spacing = %le\n Time steps = %d\n time jumps = %le\n",NC,N,dr,time,dt);
    
    double c_zero[NC];
    
    /*model parameters data*/
    fr = fopen("input/model_para","r");
    if (fr == NULL){
        printf("Unable to open file.");
        exit(0);
    }
    
    fscanf(fr,"%le%le",&kappa,&lambda);
    fprintf(fw,"kappa = %le\n lambda = %le\n",kappa,lambda);
    
    /*saving memory for chi, it would'nt allow to pass it as argument directly*/
    chi = (double **)malloc(NC*sizeof(double*));
    if(chi==NULL){
        printf("Error saving memory for chi\n");
        exit(19);
    }
    for(int i=0; i<NC; i++){
        chi[i] = (double *)malloc(NC*sizeof(double));
        if(chi[i]==NULL){
            printf("Error saving memory for chi\n");
            exit(19);
        }
    }
    
    for(int i=0; i<NC; i++){
        for(int j=i+1; j<NC; j++){
            fscanf(fr,"%le",&chi[i][j]);
            chi[j][i]=chi[i][j];
        }
        chi[i][i]=0;
    }
    fclose(fr);
    fflush(fr);
    
    fprintf(fw,"chi:\n");
    for(int i=0; i<NC; i++){
        for(int j=0; j<NC; j++){
            fprintf(fw,"%10.3e",chi[i][j]);
        }
        fprintf(fw,"\n");
    }
    

    /*Lastly initial composition profile and noise.*/
    fr = fopen("input/initial_comp","r");
    if (fr == NULL){
        printf("Unable to open file.");
        exit(0);
    }
    
    for(int i=0; i<NC; i++){
        fscanf(fr,"%le",&c_zero[i]);
    }
    
    fscanf(fr,"%le", &c_noise);
    fclose(fr);
    fflush(fr);
    fprintf(fw,"\n");
    fprintf(fw,"c_noise = %le\n", c_noise);

    /*New structure defined using the FFTW3 base variable.*/
    
    composition comp[NC];

    for(int i=0; i<NC; i++){
        comp[i].c = (fftw_complex *)fftw_malloc(N*N* sizeof(fftw_complex));
        if(comp[i].c==NULL){
            printf("Error saving memory\n");
            exit(1);
        }
    }

    /*Random generator seed, anything could have been used..*/

    (void) gsl_rng_env_setup();
    Taus = gsl_rng_taus;
    ran_num = gsl_rng_alloc (Taus);

    /*Setting the initial composition profile.*/

    double average[NC];
    for(int i=0; i<NC; i++){
        average[i]=0;
    }
    for(int i1=0; i1 < N; i1++){
        for(int i2=0; i2 < N; i2++){
            comp[NC-1].c[i2+N*i1] = 1;
            for(int i=0; i<NC-1; i++){
                __real__(comp[i].c[i2+N*i1]) = c_zero[i] + c_noise*(0.5 - gsl_rng_uniform_pos(ran_num));
                __imag__(comp[i].c[i2+N*i1]) = 0.0;
                average[i] += comp[i].c[i2+N*i1];
                /*For the last molecule, just 1- all the others*/
                comp[NC-1].c[i2+N*i1] -= comp[i].c[i2+N*i1];
            }
            if(creal(comp[NC-1].c[i2+N*i1])<0){
                printf("Error negative composition profile\n");
            }
            average[NC-1] += comp[NC-1].c[i2+N*i1];
        }
    }

    /*Writing the average composition, which */
    fprintf(fw,"average_comp:");
    for(int i=0; i<NC; i++){
        fprintf(fw,"%10.3le", average[i]/(N*N));
    }
    
    /*making sure all files are closed*/
    fclose(fw);
    fflush(fw);
   
   
    /*Performing the evolution using Variable Mobility Cahn-Hilliard equation.*/

    modelB(NC, N, dr, time, dt, kappa, chi, lambda, comp);

    /*free the composition and the random number variable.*/
    for(int i=0; i<NC; i++){
        fftw_free(comp[i].c);
    }
    gsl_rng_free(ran_num);

    return 0;
}
