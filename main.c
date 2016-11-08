#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef R_ADP
#define R_ADP
#endif

typedef double real;

const int Nx = 100;
const int Nt = 2000;

const real X = 1.;

const real g = 9.8;
FILE* fout3;


real I[3];
real Ip[3];

real dx;
real dt;

real* h;
real* hu;

real* hnext;
real* hunext;


real ic_h(real x)
{
    // if (x * x  < 0.01) return 8.;
    // return 1.0;
    // x = x -0.5;

    //  if (x * x  < 0.1) return 4.;
    // return 3.0;

    return sin(2.*M_PI*x) + 4;
}

real ic_hu(real x)
{
    return 0.;
}

void ffun(real* res, int i)
{
    real h_ = h[i + 1];
    real hu_ = hu[i + 1];
    res[0] = hu_;
    res[1] = hu_*hu_ / h_ + 0.5 * g * h_*h_;
}

void fpfun(real* res, int i) 
{
    real res1[2];
    ffun(res,i);
    ffun(res1,i+1);
    for (int k = 0; k < 2; ++k)
    {
        res[k] = (res1[k] - res[k]) / dx;
    }
}

void fmfun(real* res, int i)
{
    fpfun(res,i-1);
} 

void fjac(real* res, int i)
{
    real h_  = h[i + 1];
    real hu_ = hu[i + 1];
    res[0] = 0.; res[1] = 1.;
    res[2] = -hu_*hu_ /h_/h_ + g*h_; 
    res[3] = 2.*hu_ /h_ ; 
}

void ImM(real* A)
{
    for (int i = 0; i < 4; ++i)
    {
        A[i] *= dt/dx;
        A[i] = i%3 > 0 ? -A[i] : 1.- A[i];
    }
} 
void IpM(real* A)
{
    for (int i = 0; i < 4; ++i)
    {
        A[i] *= dt/dx;
        A[i] = i%3 > 0 ? A[i] : 1.+ A[i];
    }
} 
void Ap(real* res, int i)
{
    real res1[4]; 
    fjac(res,  i);
    fjac(res1, i+1);

    for (int c = 0; c < 4; ++c)
    {
        res[c] = (res[c] + res1[c])/2.;
    }
} 

void Am(real* res, int i)
{
    Ap(res, i-1);
}

void dot(real* res, real M[4], real v[2])
{
    res[0] = M[0] * v[0] + M[1] * v[1];
    res[1] = M[2] * v[0] + M[3] * v[1];
}

void updateBC();

real integrate(real *q)
{
    real res = 0.;
    for(int i = 0; i < Nx-1; ++i)
    {
        res += q[i+1];
    }
    return res;
}

real integrateenegy()
{
    real res = 0.;
    for(int i = 0; i < Nx-1; ++i)
    {
        res +=  hu[i+1]*hu[i+1] /h[i+1] +
                g*h[i+1]*h[i+1];
    }
    return 0.5 * res; 
}
real T = 0.;

void printsol(int tstep)
{
    FILE* fout_sol;
    char ss[100] = "";
    sprintf(ss, "out/h_%4.4d.txt", tstep);
    fout_sol = fopen(ss, "w");
    for (int i = 1; i < Nx+1; ++i)
    {    
        fprintf(fout_sol,"%14.16f\n", h[i]);
    }
    fclose(fout_sol);
}



void advInTime(int tstep)
{
    for (int i = 0; i < Nx - 1; ++i)
    {
        real matl[4];
        real vecl[2];
        real matr[4];
        real vecr[2];
        real resl[2];
        real resr[2];

        real htemp = 0., hutemp = 0.;

        Ap(matl,i); 
        ImM(matl); 
        fpfun(vecl,i); 
        dot(resl, matl, vecl);  

        Am(matr,i);
        IpM(matr);
        fmfun(vecr,i);
        dot(resr, matr, vecr);

         htemp  = dt * (resl[0] + resr[0])/ 2.;
         hutemp = dt * (resl[1] + resr[1])/ 2.;
        
         hnext [i + 1] =  h[i + 1] - htemp;   
         hunext[i + 1] = hu[i + 1] - hutemp;
    }

    real * temp = hnext;
    hnext = h;
    h = temp;

    temp = hunext;
    hunext = hu;
    hu = temp;

    updateBC();

    I[0] = integrate(h);
    I[1] = integrate(hu);
    I[2] = integrateenegy();
    fprintf(fout3, "%14.16f ", T);
    for (int i = 0; i < 3; ++i)
    {
        fprintf(fout3, "%14.16f ",dx* (I[i]- Ip[i]));
    }
    fprintf(fout3, "\n");
    for (int i = 0; i < 3; ++i)
    {
        Ip[i] = I[i];
    }
    T += dt;
    printsol(tstep+1);
}



void adaptTimeStep()
{
    real cfl = 0.;
    real qtt;
    for(int i = 0; i < Nx-1; ++i)
    {

        qtt = abs(hu[i+1] / h[i+1]) + sqrt(g * h[i+1]);
        if (cfl < qtt) cfl = qtt;
    }

    dt = dx / cfl;
}

void updateBC() 
{
   h[0] = h[1];
   h[Nx] = h[Nx-1];
   hu[0] = -hu[1];
   hu[Nx] = -hu[Nx-1];
}

void setUpIC()
{
    for(int i = 0; i < Nx-1; ++i)
    {
        real x = i*dx;
        h [i+1] = ic_h(x);
        hu[i+1] = ic_hu(x);
    }
    updateBC();
}



int main()
{
    h = (real*) malloc(sizeof(real) * (Nx + 1));
    hu = (real*) malloc(sizeof(real) * (Nx + 1));
    hnext = (real*) malloc(sizeof(real) * (Nx + 1));
    hunext = (real*) malloc(sizeof(real) * (Nx + 1));
    dx = X / (Nx-1);
    
    setUpIC();
    updateBC();
    adaptTimeStep();
    fout3 = fopen("cons.txt","w");
    I[0] = integrate(h);
    I[1] = integrate(hu);  
    I[2] = integrateenegy();
    printsol(0);
    for (int i = 0; i < 3; ++i)
    {
        Ip[i] = I[i];
    }
    for (int tstep = 0; tstep < Nt; ++tstep)
    {
        advInTime(tstep);
    #ifdef R_ADP
        adaptTimeStep();
    #endif

    }
    fclose(fout3);

    free(h);
    free(hu);
    free(hnext);
    free(hunext);
    
    FILE* fout1 = fopen("out/x.txt","w");
    
    for (int i = 0; i < Nx; ++i)
    {
        fprintf(fout1, "%f ", dx*i);
    }

    fprintf(fout1, "\n");

    fclose(fout1);
    return 0;
}
