#include <stdint.h>
#include "cmsis.h"

extern int jump1521divsteps_mod3(int minusdelta, int* M, int* f, int* g);
int polyinv4591761_mod3(uint8_t* H, uint8_t* G);

int polyinv4591761_mod3(uint8_t* H, uint8_t* G){
    uint8_t f[768],g[768],M[7680];
    int i,j;
    int minusdelta=-1;

    for(i=0;i<768;i++)f[i]=0;
    for(i=761;i<768;i++)g[i]=0;
    f[0]=1;
    f[760]=2;
    f[761]=2;
    for(i=0;i<761;i++)g[i]=G[760-i];

    minusdelta = jump1521divsteps_mod3(minusdelta,M,f,g);

    for(i=761;i<768;i++)H[i]=0;
    for(i=0;i<761;i++){
        H[i] = (M[3847-i]*M[0])%3;
    }

    return minusdelta;
}