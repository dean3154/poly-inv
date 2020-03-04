#include <stdint.h>
#include "cmsis.h"

extern void gf_polymul_256x256(void *h, void *f, void *g);
void gf_polymul_256x512(int *h, int *f, int *g);

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

void gf_polymul_256x512(int *h, int *f, int *g){
    int i,*start,*H1,*H2;
    int h1[256],h2[256];
    gf_polymul_256x256(h1,f,g);
    gf_polymul_256x256(h2,f,g+128);
    for(start = h, H1 = h1, H2 = h2, i=0;i<128;i++){
        *(start++) = *(H1++);
    }
    for(i=0;i<128;i++){
        *(start++) = barrett_16x2i(__SADD16(*(H1++),*(H2++)));
    }
    for(i=0;i<128;i++){
        *(start++) = *(H2++);
    }
}