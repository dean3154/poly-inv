#include <stdint.h>
#include "cmsis.h"
#include "bred_int.h"

extern int jump1521divsteps(int minusdelta, int* M, int* f, int* g);
extern void gf_polymul_768x768(void *h, void *f, void *g);
int polyinv4591761(int16_t* H, int16_t* G);

#define q 4591
#define qR2inv 935519 // round(2^32/q)
#define _2P15 (1 << 15)

static int16_t pow4589(int16_t k){
    int16_t r=k;       //1
    r=bred(r*r,q,-qR2inv);//10
    r=bred(r*r,q,-qR2inv);//100
    r=bred(r*r,q,-qR2inv);//1000
    r=bred(r*r,q,-qR2inv);//10000
    r=bred(r*k,q,-qR2inv);//10001
    r=bred(r*r,q,-qR2inv);//100010
    r=bred(r*k,q,-qR2inv);//100011
    r=bred(r*r,q,-qR2inv);//1000110
    r=bred(r*k,q,-qR2inv);//1000111
    r=bred(r*r,q,-qR2inv);//10001110
    r=bred(r*k,q,-qR2inv);//10001111
    r=bred(r*r,q,-qR2inv);//100011110
    r=bred(r*r,q,-qR2inv);//1000111100
    r=bred(r*k,q,-qR2inv);//1000111101
    r=bred(r*r,q,-qR2inv);//10001111010
    r=bred(r*k,q,-qR2inv);//10001111011
    r=bred(r*r,q,-qR2inv);//100011110110
    r=bred(r*r,q,-qR2inv);//1000111101100
    r=bred(r*k,q,-qR2inv);//1000111101101
    return r;
}

int polyinv4591761(int16_t* H, int16_t* const G){
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
    for(i=0;i<761;i++){
        H[i] = bred(M[3847-i]*k,q,-qR2inv);
    }

    return minusdelta;
}