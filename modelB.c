/*The program has some commented parts, feel free to modify them. Some are parts that helped in th eproblem solving and others are just there in case some changes wanted to be done. */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include "fftw3.h"
#include<gsl/gsl_math.h>
#include "../header_file/headers.h"




typedef struct{
    fftw_plan plan;
}fplan;

int delta(int , int );



void modelB(int NC, int N, double dr, int time, double dt, double kappa, double **chi, double lambda, composition *comp){
    
    /*for the time beeing, let us just consider the NC=2 and dim=2 cases.*/
    
    FILE *fw;
    char file_name[50];
    /*char ps_file_name[50];*/
    double kv[2];
    double k2, k4, delk, denominator;
    int i, j, k, halfN, temp=0, dimtot=N*N, over=0;
    int *c;
    
    composition g[NC], chem[NC][2], y[NC-1][2], copy_comp[NC];
    fplan plan1[NC], plan2[NC], plan3[NC][2], plan4[NC-1][2], plan5[NC-1];
    
    /*We save memory for our structures. A fucntion would allow us to make it simpler and nicer, but we might have problems to free the space later.*/
    
    
    for(i=0; i<NC; i++){
        g[i].c = (fftw_complex *)fftw_malloc(dimtot*sizeof(fftw_complex));
        if(g[i].c == NULL){
            printf("Unexpected error while saving memory for the composition inside the structure\n");
            exit(1);
        }
        copy_comp[i].c = (fftw_complex *)fftw_malloc(dimtot*sizeof(fftw_complex));
        if(copy_comp[i].c == NULL){
            printf("Unexpected error while saving memory for the composition inside the structure\n");
            exit(1);
        }
        
    }
    
    //dimensional variables
    
    for(i=0; i<NC; i++){
        for(j=0; j<2; j++){
            chem[i][j].c = (fftw_complex*)fftw_malloc(dimtot*sizeof(fftw_complex));
            if(chem[i][j].c==NULL){
                printf("Unexpected error while saving memory for the composition inside the structure\n");
                exit(1);
            }
            if(i<NC-1){
                y[i][j].c = (fftw_complex*)fftw_malloc(dimtot*sizeof(fftw_complex));
                if(y[i][j].c==NULL){
                    printf("Unexpected error while saving memory for the composition inside the structure\n");
                    exit(1);
                }
            }
        }
    }
    
    /*Will tell which is the max concentration per place*/
    c = (int *)malloc(dimtot*sizeof(int));
    if(c==NULL){
        printf("Unexpected error while saving memory for the composition\n");
        exit(1);
    }
    
    
    /*We define A as in the reference for soft matter, being A=0.5*max(chi):*/
    
    double A=0;
    for(i=0; i<NC; i++){
        for(j=i+1; j<NC; j++){
            if (A<0.5*chi[i][j]){
                A=0.5*chi[i][j];
            }
        }
    }
    
    /*Same for the values of delkx/ky, WHY ARE WE DIVING BY DR IN DELK?*/
    
    delk = (2*M_PI)/(N*dr);
    halfN = (int) N/2;
    
    
    /*Let us define the plans needed. this works fine for dimensions, but wee need one separate plan per NC*/
    
    for(i=0; i<NC; i++){
        plan1[i].plan = fftw_plan_dft_2d(N, N, comp[i].c, comp[i].c, FFTW_FORWARD, FFTW_ESTIMATE);
        plan2[i].plan = fftw_plan_dft_2d(N, N, g[i].c, g[i].c, FFTW_FORWARD, FFTW_ESTIMATE);
        for(j=0; j<2; j++){
            plan3[i][j].plan = fftw_plan_dft_2d(N, N, chem[i][j].c, chem[i][j].c, FFTW_BACKWARD, FFTW_ESTIMATE);
            if(i<NC-1){
                plan4[i][j].plan = fftw_plan_dft_2d(N, N, y[i][j].c, y[i][j].c, FFTW_FORWARD, FFTW_ESTIMATE);
            }
        }
        if(i<NC-1){
            plan5[i].plan = fftw_plan_dft_2d(N, N, comp[i].c, comp[i].c, FFTW_BACKWARD, FFTW_ESTIMATE);
        }
    }
    
    sprintf(file_name, "output/time%d.txt",temp);
    fw = fopen(file_name, "w");
    if(fw==NULL){
        exit(1);
    }
    
    for(i=0; i<N; i++){
        for(j=0; j<N; j++){
            double max[2];
            max[0]=0;
            max[1]=0;
            for(k=0; k<NC; k++){
                if(max[0]<creal(comp[k].c[j+N*i])){
                    max[0] = creal(comp[k].c[j+N*i]);
                    max[1] = k+1;
                }
            }
            c[j+i*N] = (int)max[1];
            fprintf(fw, "%10d ", c[j+i*N]);
        }
        fprintf(fw,"\n");
    }
    (void) fclose (fw);
    fflush(fw);
    
    
    for(temp=1; temp<time+1; temp++){
        
        /*We copy the values of comp, so we can use them later, and we define the values of g and h, and then transform*/
        for(k=0;k<NC;k++){
            for(i=0; i<N; i++){
                for(j=0; j<N; j++){
                        copy_comp[k].c[j+i*N] = comp[k].c[j+i*N];
                        g[k].c[j+i*N] = 1 + log(creal(comp[k].c[j+i*N]));
                }
            }
        }
        
        /*we execute the first plans, so we are now in the phase space*/
        
        for(i=0; i<NC; i++){
            fftw_execute(plan1[i].plan);
            fftw_execute(plan2[i].plan);
        }
        
        
        /*We define chem takling into account the way kv are saved when doing the transforms.*/
        for(k=0; k<NC; k++){
            for(i=0; i<N; i++){
                kv[0] = ( (i > halfN) ? (i-N) : i )*delk;
                for(j=0; j < N; j++){
                    kv[1] = ( (j > halfN) ? (j-N) : j )*delk;
                    k2 = kv[0]*kv[0] + kv[1]*kv[1];
                    for(int l=0; l<2; l++){
                        chem[k][l].c[j+i*N] = g[k].c[j+i*N];
                        for(int m=0; m<NC; m++){
                            chem[k][l].c[j+i*N] += chi[k][m]*(1-lambda*lambda*k2)*comp[m].c[j+N*i];
                        }
                        chem[k][l].c[j+i*N] = chem[k][l].c[j+i*N] * (_Complex_I) * kv[l];
                    }
                }
            }
        }
        
        /*We execute the following plans, taking into account that we are back to real space, so we have to normalize our functions. Also these functions should be real...*/
        
        for(k=0; k<NC; k++){
            for(j=0; j<2; j++){
                fftw_execute(plan3[k][j].plan);
            }
        }
        
        /*int compt=0;
        printf("hi ha algun valor imaginari no nul per pas de simulació t=%d?\n", temp);*/
        for(k=0; k < NC; k++){
            for(int l=0; l<2; l++){
                for(i=0; i<N; i++){
                    for(j=0; j<N; j++){
                        chem[k][l].c[j + i*N] =  chem[k][l].c[j + i*N]/dimtot;
                        /*if(cimag(chem[k][l].c[j + i*N])!=0){
                            printf(" %5d %5d %5d %5d \n", k, l, i, j);
                            compt++;
                        }*/
                        /*__imag__(chem[k][l].c[j + i*N]) = 0;*/
                    }
                }
            }
        }
        /*printf("There are so many values %5d with non-imaginary part\n", compt);*/
        
        /*We define the following functions y1, y2*/
        
        for(k=0; k < NC-1; k++){
            for(i=0; i<N; i++){
                for(j=0; j<N; j++){
                    for(int l=0; l<2; l++){
                        y[k][l].c[j + i*N] = 0;
                        for(int m=0; m<NC; m++){
                            y[k][l].c[j + i*N] += chem[m][l].c[j+N*i]*(delta(k,m) - copy_comp[m].c[j+N*i]);
                        }
                        y[k][l].c[j + i*N] = y[k][l].c[j + i*N]*copy_comp[k].c[j+N*i];
                    }
                }
            }
        }
        
        
        /*we execute the following plans, going back to Fourier space*/
        for(i=0; i<NC-1; i++){
            for(j=0; j<2; j++){
                fftw_execute(plan4[i][j].plan);
            }
        }
        
        for(k=0; k<NC-1; k++){
            for(i=0; i<N; i++){
                kv[0] = ( (i > halfN) ? (i-N) : i )*delk;
                for(j=0; j < N; ++j){
                    kv[1] = ( (j > halfN) ? (j-N) : j )*delk;
                    k2 = kv[0]*kv[0] + kv[1]*kv[1];
                    k4 = k2*k2;
                    denominator = kappa*dt/(1+A*kappa*lambda*lambda*k4*dt);
                        for(int l=0; l<2; l++){
                            comp[k].c[j + i*N] += denominator*(_Complex_I)*kv[l]*y[k][l].c[j + i*N];
                    }
                }
            }
        }
        
        
        /*last but not least, execute the last plan and normalize, so we are back to normal space*/
        
        for(i=0; i<NC-1; i++){
            fftw_execute(plan5[i].plan);
        }
        
        
        for(i=0; i<N; i++){
            for(j=0; j<N; j++){
                comp[NC-1].c[j + i*N] = 1;
                for(k=0; k<NC-1; k++){
                    comp[k].c[j + i*N]=comp[k].c[j + i*N]/dimtot;
                    __imag__(comp[k].c[j + i*N]) = 0;
                    comp[NC-1].c[j + i*N] -= comp[k].c[j + i*N];
                    __imag__(comp[NC-1].c[j + i*N]) = 0;
                    if(creal(comp[k].c[j + i*N])<=0 || creal(comp[k].c[j + i*N])>=1){
                        over = 1;
                    }else if(creal(comp[NC-1].c[j + i*N])<=0){
                        over = 1;
                    }
                }
            }
        }
        
        if(over!=0){
            printf("We have reached a value over 1 or below 0, there are some numerical problems\n");
            /*printf("The values for the composition in the first 5 points of the grid are:\n");
            for(k=0;k<NC;k++){
                for(i=0; i<NC; i++){
                    printf(" %10.3e + i*%10.3e ", creal(comp[k].c[i]), cimag(comp[k].c[i]));
                }
                printf("\n");
            }
            printf("\n");*/
            break;
        }
        
        
        /*at this point we could save the matrix and it would allow us to see the time evolution as expected*/
        if (temp%1000 == 0){
            
            sprintf(file_name, "output/time%d.txt",temp);
            fw = fopen(file_name, "w");
            if(fw==NULL){
                exit(1);
            }
            for(i=0; i<N; i++){
                for(j=0; j<N; j++){
                    double max[2];
                    max[0]=0;
                    max[1]=0;
                    for(k=0; k<NC; k++){
                        if(max[0]<creal(comp[k].c[j+N*i])){
                            max[0] = creal(comp[k].c[j+N*i]);
                            max[1] = k+1;
                        }
                    }
                    c[j+i*N] = (int)max[1];
                    fprintf(fw, "%10d ", c[j+i*N]);
                }
                fprintf(fw,"\n");
            }
            (void) fclose (fw);
            fflush(fw);
            /*sprintf(file_name, "output/simulation_data.txt");
            fw = fopen(file_name, "a");
            if(fw==NULL){
                exit(1);
            }
            double average[NC];
            for(int i1=0; i1<NC; i1++){
                average[i1]=0;
            }
            for(int i1=0; i1<NC; i1++){
                for(int i2=0; i2<N; i2++){
                    for(int i3=0; i3<N; i3++){
                        average[i1]+=comp[i1].c[i3+N*i2];
                    }
                }
            }
            fprintf(fw,"\n average composition at time %8d is:\n", temp);
            for(int i1=0; i1<NC; i1++){
                fprintf(fw,"%10.3le", average[i1]/(N*N));
            }
            fprintf(fw,"\n");
            (void) fclose (fw);
            fflush(fw);*/
        }
        
    }
    
    
    
    //Let us free the space before we do aNthing else.
    for(i=0; i<NC; i++){
        fftw_free(copy_comp[i].c);
        fftw_free(g[i].c);
        for(j=0; j<2; j++){
            fftw_free(chem[i][j].c);
            if(i<NC-1){
                fftw_free(y[i][j].c);
            }
        }
        
    }
    
    free(c);
    
    for(i=0; i<NC; i++){
        fftw_destroy_plan(plan1[i].plan);
        fftw_destroy_plan(plan2[i].plan);
        for(j=0; j<2; j++){
            fftw_destroy_plan(plan3[i][j].plan);
            if(i<NC-1){
                fftw_destroy_plan(plan4[i][j].plan);
            }
        }
        if(i<NC-1){
            fftw_destroy_plan(plan5[i].plan);
        }
    }
}




int delta(int i, int j){
    if(i==j){
        return 1;
    }else{
        return 0;
    }
}


