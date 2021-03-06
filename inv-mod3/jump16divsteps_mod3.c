#include <stdint.h>
#include "red_mod3_int.h"
#include "cmsis.h"

extern int jump8divsteps_mod3(int minusdelta, int* M, int* f, int* g);
extern void gf_polymul_8x8_mod3(int* h, int* f, int* g);
int jump16divsteps_mod3(int minusdelta, int* M, int* f, int* g);
void gf_polymul_8x8_2x2_x2p2_mod3(int* V,int* M,int *fh, int* gh);
void gf_polymul_8x8_2x2_x_2x2_mod3(int* M,int* M1,int* M2);

static int C16_1[5], C16_2[5];
int * CC16_1 = (int *)((void *)C16_1 + 1);
int * CC16_2 = (int *)((void *)C16_2 + 1);

void gf_polymul_8x8_2x2_x2p2_mod3(int *V,int *M,int *fh,int *gh){
  int i, T, *X, *Y, *Z, *W;

  gf_polymul_8x8_mod3(CC16_1, M+4, fh); // x * u * fh
  gf_polymul_8x8_mod3(CC16_2, M+6, gh); // x * v * gh

  for (X=V, Y=C16_1, Z=C16_2, W=M, i=2; i>0; i--) {// x(u fh+v gh)+f1
    *(X++) = add_ub3(add_ub3(*(W++),*(Y++)),*(Z++));
  }
  for (i=2; i>0; i--) {
    *(X++) = add_ub3(*(Y++),*(Z++));
  }

  gf_polymul_8x8_mod3(V+4, M+8, fh); // r * fh
  gf_polymul_8x8_mod3(CC16_1, M+10, gh); // s * gh

  for (Y=CC16_1, i=2; i>0; i--) {// r fh+s gh+g1
    T = add_ub3(add_ub3(*(W++),*(Y++)),*X); *(X++) = T;
  }
  for (i=2; i>0; i--) {
    T = add_ub3(*X, *(Y++)); *(X++) = T;
  }
}

void gf_polymul_8x8_2x2_x_2x2_mod3(int *M, int *M1, int *M2){
  int i, T, *X, *Y;

  gf_polymul_8x8_mod3(CC16_1, M2, M1); // x * u2 * u1
  gf_polymul_8x8_mod3(M, M2+2, M1+4); // v2 * r1
  for (i=4, X=M, Y=C16_1; i>0; i--) { // u = x * u2 * u1 + v2 * r1
    T = add_ub3(*X, *(Y++)); *(X++) = T;
  }
  gf_polymul_8x8_mod3(CC16_1, M2, M1+2); // x * u2 * v1
  gf_polymul_8x8_mod3(M+4, M2+2, M1+6); // v2 * s1
  for (i=4, Y=C16_1; i>0; i--) { // v = x * u2 * v1 + v2 * s1
    T = add_ub3(*X, *(Y++)); *(X++) = T;
  }
  gf_polymul_8x8_mod3(CC16_1, M2+4, M1); // x * r2 * u1
  gf_polymul_8x8_mod3(M+8, M2+6, M1+4); // s2 * r1
  for (i=4, Y=C16_1; i>0; i--) { // r = x * r2 * u1 + s2 * r1
    T = add_ub3(*X, *(Y++)); *(X++) = T;
  }
  gf_polymul_8x8_mod3(CC16_1, M2+4, M1+2); // x * r2 * v1
  gf_polymul_8x8_mod3(M+12, M2+6, M1+6); // s2 * s1
  for (i=4, Y=C16_1; i>0; i--) { // s = x * r2 * v1 + s2 * s1
    T = add_ub3(*X, *(Y++)); *(X++) = T;
  }
}

int jump16divsteps_mod3(int minusdelta, int* M, int* f, int* g){
  int M1[12], M2[12], fg[8];

  minusdelta = jump8divsteps_mod3(minusdelta, M1, f, g);
  gf_polymul_8x8_2x2_x2p2_mod3(fg, M1, f+2, g+2);
  minusdelta = jump8divsteps_mod3(minusdelta, M2, fg, fg+4);
  gf_polymul_8x8_2x2_x2p2_mod3(M, M2, fg+2, fg+6);
  gf_polymul_8x8_2x2_x_2x2_mod3(M+8, M1+4, M2+4);
  return minusdelta;
}
