#include <stdint.h>
#include "cmsis.h"

extern int jump768divsteps(int minusdelta, int *M, int *f, int *g);
extern int jump753divsteps(int minusdelta, int *M, int *f, int *g);
extern void gf_polymul_768x768(void *h, void *f, void *g);

void gf_polymul_768x768_2x2_x_2x2 (int *M, int *M1, int *M2); // M = M2 x M1 
int jump1521divsteps(int minusdelta, int *M, int *f, int *g);

#define q 4591
#define qR2inv 935519 // round(2^32/q)
#define _2P15 (1 << 15)

#if 1
// result range: +- 2295 (note: 3 loads for _2P15 and the longer qR2inv)
static inline int barrett_16x2i(int X) {
  int32_t QL = __SMLAWB(qR2inv,X,_2P15);
  int32_t QH = __SMLAWT(qR2inv,X,_2P15);
  int32_t SL = __SMULBT(q,QL);
  int32_t SH = __SMULBT(q,QH);
  return(__SSUB16(X,__PKHBT(SL,SH,16)));
}

#else 
#define barrett_16x2i(A) (A)
#endif

static int B1521_1[769];
int * BB1521_1 = (int *)((void *)B1521_1 + 2);

void gf_polymul_768x768_2x2_x_2x2 (int *M, int *M1, int *M2){ //only v = g^-1 mod f
    int i, T, *X, *Y;
/*
    gf_polymul_768x768(BB1521_1, M2, M1); // x * u2 * u1 
    gf_polymul_768x768(M, M2+384,M1+768); // v2 * r1
    for (i=768, X=M, Y=B1521_1; i>0; i--) {	// u = x u2 u1 + v2 r1
     T = barrett_16x2i(__SADD16(*X,*(Y++)));
       *(X++) = T;
    }
*/
    gf_polymul_768x768(BB1521_1, M2, M1+384); // x * u2 * v1 
    gf_polymul_768x768(M+768, M2+384,M1+1152); // v2 * s1
    for (i=768, X=M+768, Y=B1521_1; i>0; i--) {	// v = x u2 v1 + v2 s1
     T = barrett_16x2i(__SADD16(*X,*(Y++)));
       *(X++) = T;
    }
/*
    gf_polymul_768x768(BB1521_1, M2+768, M1); // x * r2 * u1 
    gf_polymul_768x768(M+1536, M2+1152,M1+768); // s2 * r1
    for (i=768, Y=B1521_1; i>0; i--) {	// r = x r2 u1 + s2 r1
     T = barrett_16x2i(__SADD16(*X,*(Y++)));
       *(X++) = T;
    }

    gf_polymul_768x768(BB1521_1, M2+768, M1+384); // x * r2 * v1 
    gf_polymul_768x768(M+2304, M2+1152,M1+1152); // s2 * s1
    for (i=768, Y=B1521_1; i>0; i--) {	// s = x r2 v1 + s2 s1
     T = barrett_16x2i(__SADD16(*X,*(Y++)));
       *(X++) = T;
    }
*/
}

int jump1521divsteps(int minusdelta, int *M, int *f, int *g){
    int M1[2304],M2[2304];
    int i;

    minusdelta = jump768divsteps(minusdelta,M1,f,g);
    minusdelta = jump753divsteps(minusdelta,M2,M1,M1+384);

    gf_polymul_768x768_2x2_x_2x2(M+768,M1+768,M2+768);

    for(i=0;i<768;i++)M[i]=M2[i];

    return minusdelta;
}