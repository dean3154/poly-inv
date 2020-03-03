#include <stdint.h>
#include "cmsis.h"

extern int jump1521divsteps(int minusdelta, int* M, int* f, int* g);
int polyinv4591761(int16_t* H, int16_t* G);

int16_t mult_reduce(int16_t a, int16_t b){
    return (a*b)%4591;
}

int16_t pow4589(int16_t k){
    int16_t r=k;       //1
    r=mult_reduce(r,r);//10
    r=mult_reduce(r,r);//100
    r=mult_reduce(r,r);//1000
    r=mult_reduce(r,r);//10000
    r=mult_reduce(r,k);//10001
    r=mult_reduce(r,r);//100010
    r=mult_reduce(r,k);//100011
    r=mult_reduce(r,r);//1000110
    r=mult_reduce(r,k);//1000111
    r=mult_reduce(r,r);//10001110
    r=mult_reduce(r,k);//10001111
    r=mult_reduce(r,r);//100011110
    r=mult_reduce(r,r);//1000111100
    r=mult_reduce(r,k);//1000111101
    r=mult_reduce(r,r);//10001111010
    r=mult_reduce(r,k);//10001111011
    r=mult_reduce(r,r);//100011110110
    r=mult_reduce(r,r);//1000111101100
    r=mult_reduce(r,k);//1000111101101
    return r;
}

int polyinv4591761(int16_t* H, int16_t* G){
    int16_t f[768],g[768],M[7680];
    int i,j;
    int k,minusdelta=-1;

    for(i=0;i<768;i++)f[i]=0;
    for(i=761;i<768;i++)g[i]=0;
    f[0]=1;
    f[760]=-1;
    f[761]=-1;
    for(i=0;i<761;i++)g[i]=G[760-i];

    minusdelta = jump1521divsteps(minusdelta,M,f,g);

    k=pow4589(M[0]);
    for(i=0;i<768;i++)H[i]=0;
    for(i=0;i<761;i++){
        H[i] = (M[3847-i]*k)%4591;
    }

    return minusdelta == 0 ? 1 : 0;
}