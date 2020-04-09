#include <stdint.h>
#include "red_mod3_int.h"
#include "cmsis.h"

extern void gf_polymul_256x256_mod3(void *h, void *f, void *g);
void gf_polymul_256x512_mod3(int *h, int *f, int *g);

void gf_polymul_256x512_mod3(int *h, int *f, int *g){
    int i,*start,*H1,*H2;
    int h1[128],h2[128];
    gf_polymul_256x256_mod3(h1,f,g);
    gf_polymul_256x256_mod3(h2,f,g+64);
    for(start = h, H1 = h1, H2 = h2, i=0;i<64;i++){
        *(start++) = *(H1++);
    }
    for(i=0;i<64;i++){
        *(start++) = (add_ub3(*(H1++),*(H2++)));
    }
    for(i=0;i<64;i++){
        *(start++) = *(H2++);
    }
}