#include <stdint.h>
#include "cmsis.h"

extern void bs3_mul64(void*, void*, void*);
extern void bs3_mul64_negc(void*, void*, void*);
void bs3_mul64_nega(uint32_t*, uint32_t*, uint32_t*);
void polyf3_add_packed(uint32_t *tritOut, uint32_t *tritIn1, uint32_t *tritIn2);
void polyf3_sub_packed(uint32_t *tritOut, uint32_t *tritIn1, uint32_t *tritIn2);
void fft48_64(int fpad[48][4]);
void unfft48_64(int fpad[48][4]);
void bs3_mul768_schonhage(uint32_t*, uint32_t*,uint32_t*);

void bs3_mul64_nega(uint32_t* h, uint32_t* f, uint32_t* g){
    uint32_t htemp[8];
    bs3_mul64(htemp,f,g);
    polyf3_sub_packed(h,htemp,htemp+4);
}

void polyf3_add_packed(uint32_t *tritOut, uint32_t *tritIn1, uint32_t *tritIn2) {
    uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    uint32_t *readOp1, *readOp2, *writeSum, *bndry;
    readOp1 = tritIn1;
    readOp2 = tritIn2;
    writeSum = tritOut;
    bndry = writeSum + 4;
    while (writeSum != bndry) {
        r0 = *(readOp1++); // a0
        r1 = *(readOp1++); // a1
        r2 = *(readOp2++); // b0
        r3 = *(readOp2++); // b1

        r4 = *(readOp1++); // a0
        r5 = *(readOp1++); // a1
        r6 = *(readOp2++); // b0
        r7 = *(readOp2++); // b1

        r8 = r1 ^ r2;
        r1 = r0 ^ r3;
        r1 = r1 & r8;
        r8 = r8 ^ r3;
        r0 = r0 ^ r2;
        r0 = r0 | r8;

        r8 = r5 ^ r6;
        r5 = r4 ^ r7;
        r5 = r5 & r8;
        r8 = r8 ^ r7;
        r4 = r4 ^ r6;
        r4 = r4 | r8;

        *(writeSum++) = r0;
        *(writeSum++) = r1;
        *(writeSum++) = r4;
        *(writeSum++) = r5;
    }
}

void polyf3_sub_packed(uint32_t *tritOut, uint32_t *tritIn1, uint32_t *tritIn2) {
  uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
  uint32_t *readOp1, *readOp2, *writeSum, *bndry;
  readOp1 = tritIn1;
  readOp2 = tritIn2;
  writeSum = tritOut;
  bndry = writeSum + 4;
  while (writeSum != bndry) {
    r0 = *(readOp1++); // a0
    r1 = *(readOp1++); // a1
    r2 = *(readOp2++); // b0
    r3 = *(readOp2++); // b1

    r4 = *(readOp1++); // a0
    r5 = *(readOp1++); // a1
    r6 = *(readOp2++); // b0
    r7 = *(readOp2++); // b1

    r8 = r0 ^ r2;
    r0 = r1 ^ r3;
    r0 = r0 | r8;
    r8 = r8 ^ r3;
    r1 = r1 ^ r2;
    r1 = r1 & r8;

    r8 = r4 ^ r6;
    r4 = r5 ^ r7;
    r4 = r4 | r8;
    r8 = r8 ^ r7;
    r5 = r5 ^ r6;
    r5 = r5 & r8;

    *(writeSum++) = r0;
    *(writeSum++) = r1;
    *(writeSum++) = r4;
    *(writeSum++) = r5;
  }
}

void fft48_64(int fpad[48][4]){
    int i;
    for(i=0;i<16;i++)polyf3_add_packed(fpad+32+i,fpad+i,fpad+16+i);

    for(i=0;i<16;i++){
        polyf3_butterfly64_CT(fpad+i,fpad+16+i,32);
    }

    for(i=0;i<8;i++){
        polyf3_butterfly64_CT(fpad+i,fpad+8+i,16);
        polyf3_butterfly64_CT(fpad+16+i,fpad+24+i,48);
        polyf3_butterfly64_CT(fpad+32+i,fpad+40+i,0);
    }
    for(i=0;i<4;i++){
        polyf3_butterfly64_CT(fpad+i,fpad+4+i,8);
        polyf3_butterfly64_CT(fpad+8+i,fpad+12+i,40);
        polyf3_butterfly64_CT(fpad+16+i,fpad+20+i,24);
        polyf3_butterfly64_CT(fpad+24+i,fpad+28+i,56);
        polyf3_butterfly64_CT(fpad+32+i,fpad+36+i,0);
        polyf3_butterfly64_CT(fpad+40+i,fpad+44+i,32);
    }
    for(i=0;i<2;i++){
        polyf3_butterfly64_CT(fpad+i,fpad+2+i,4);
        polyf3_butterfly64_CT(fpad+4+i,fpad+6+i,36);
        polyf3_butterfly64_CT(fpad+8+i,fpad+10+i,20);
        polyf3_butterfly64_CT(fpad+12+i,fpad+14+i,52);
        polyf3_butterfly64_CT(fpad+16+i,fpad+18+i,12);
        polyf3_butterfly64_CT(fpad+20+i,fpad+22+i,44);
        polyf3_butterfly64_CT(fpad+24+i,fpad+26+i,28);
        polyf3_butterfly64_CT(fpad+28+i,fpad+30+i,60);
        polyf3_butterfly64_CT(fpad+32+i,fpad+34+i,0);
        polyf3_butterfly64_CT(fpad+36+i,fpad+38+i,32);
        polyf3_butterfly64_CT(fpad+40+i,fpad+42+i,16);
        polyf3_butterfly64_CT(fpad+44+i,fpad+46+i,48);
    }
    polyf3_butterfly64_CT(fpad,fpad+1,2);
    polyf3_butterfly64_CT(fpad+2,fpad+3,34);
    polyf3_butterfly64_CT(fpad+4,fpad+5,18);
    polyf3_butterfly64_CT(fpad+6,fpad+7,50);
    polyf3_butterfly64_CT(fpad+8,fpad+9,10);
    polyf3_butterfly64_CT(fpad+10,fpad+11,42);
    polyf3_butterfly64_CT(fpad+12,fpad+13,26);
    polyf3_butterfly64_CT(fpad+14,fpad+15,58);
    polyf3_butterfly64_CT(fpad+16,fpad+17,6);
    polyf3_butterfly64_CT(fpad+18,fpad+19,38);
    polyf3_butterfly64_CT(fpad+20,fpad+21,22);
    polyf3_butterfly64_CT(fpad+22,fpad+23,54);
    polyf3_butterfly64_CT(fpad+24,fpad+25,14);
    polyf3_butterfly64_CT(fpad+26,fpad+27,46);
    polyf3_butterfly64_CT(fpad+28,fpad+29,30);
    polyf3_butterfly64_CT(fpad+30,fpad+31,62);
    polyf3_butterfly64_CT(fpad+32,fpad+33,0);
    polyf3_butterfly64_CT(fpad+34,fpad+35,32);
    polyf3_butterfly64_CT(fpad+36,fpad+37,16);
    polyf3_butterfly64_CT(fpad+38,fpad+39,48);
    polyf3_butterfly64_CT(fpad+40,fpad+41,8);
    polyf3_butterfly64_CT(fpad+42,fpad+43,40);
    polyf3_butterfly64_CT(fpad+44,fpad+45,24);
    polyf3_butterfly64_CT(fpad+46,fpad+47,56);
    
}

void unfft48_64(int fpad[48][4]){
    int i;

    polyf3_butterfly64_GS(fpad,fpad+1,62);
    polyf3_butterfly64_GS(fpad+2,fpad+3,30);
    polyf3_butterfly64_GS(fpad+4,fpad+5,46);
    polyf3_butterfly64_GS(fpad+6,fpad+7,14);
    polyf3_butterfly64_GS(fpad+8,fpad+9,54);
    polyf3_butterfly64_GS(fpad+10,fpad+11,22);
    polyf3_butterfly64_GS(fpad+12,fpad+13,38);
    polyf3_butterfly64_GS(fpad+14,fpad+15,6);
    polyf3_butterfly64_GS(fpad+16,fpad+17,58);
    polyf3_butterfly64_GS(fpad+18,fpad+19,26);
    polyf3_butterfly64_GS(fpad+20,fpad+21,42);
    polyf3_butterfly64_GS(fpad+22,fpad+23,10);
    polyf3_butterfly64_GS(fpad+24,fpad+25,50);
    polyf3_butterfly64_GS(fpad+26,fpad+27,18);
    polyf3_butterfly64_GS(fpad+28,fpad+29,34);
    polyf3_butterfly64_GS(fpad+30,fpad+31,2);
    polyf3_butterfly64_GS(fpad+32,fpad+33,64);
    polyf3_butterfly64_GS(fpad+34,fpad+35,32);
    polyf3_butterfly64_GS(fpad+36,fpad+37,48);
    polyf3_butterfly64_GS(fpad+38,fpad+39,16);
    polyf3_butterfly64_GS(fpad+40,fpad+41,56);
    polyf3_butterfly64_GS(fpad+42,fpad+43,24);
    polyf3_butterfly64_GS(fpad+44,fpad+45,40);
    polyf3_butterfly64_GS(fpad+46,fpad+47,8);

    for(i=0;i<2;i++){
        polyf3_butterfly64_GS(fpad+i,fpad+2+i,60);
        polyf3_butterfly64_GS(fpad+4+i,fpad+6+i,28);
        polyf3_butterfly64_GS(fpad+8+i,fpad+10+i,44);
        polyf3_butterfly64_GS(fpad+12+i,fpad+14+i,12);
        polyf3_butterfly64_GS(fpad+16+i,fpad+18+i,52);
        polyf3_butterfly64_GS(fpad+20+i,fpad+22+i,20);
        polyf3_butterfly64_GS(fpad+24+i,fpad+26+i,36);
        polyf3_butterfly64_GS(fpad+28+i,fpad+30+i,4);
        polyf3_butterfly64_GS(fpad+32+i,fpad+34+i,64);
        polyf3_butterfly64_GS(fpad+36+i,fpad+38+i,32);
        polyf3_butterfly64_GS(fpad+40+i,fpad+42+i,48);
        polyf3_butterfly64_GS(fpad+44+i,fpad+46+i,16);
    }

    for(i=0;i<4;i++){
        polyf3_butterfly64_GS(fpad+i,fpad+4+i,56);
        polyf3_butterfly64_GS(fpad+8+i,fpad+12+i,24);
        polyf3_butterfly64_GS(fpad+16+i,fpad+20+i,40);
        polyf3_butterfly64_GS(fpad+24+i,fpad+28+i,8);
        polyf3_butterfly64_GS(fpad+32+i,fpad+36+i,64);
        polyf3_butterfly64_GS(fpad+40+i,fpad+44+i,32);
    }

    for(i=0;i<8;i++){
        polyf3_butterfly64_GS(fpad+i,fpad+8+i,48);
        polyf3_butterfly64_GS(fpad+16+i,fpad+24+i,16);
        polyf3_butterfly64_GS(fpad+32+i,fpad+40+i,64);
    }

    for(i=0;i<16;i++){
        polyf3_butterfly64_GS(fpad+i,fpad+16+i,32);
    }

    for(i=0;i<32;i++){
        fpad[i][1] = fpad[i][1] ^ fpad[i][0];
        fpad[i][3] = fpad[i][3] ^ fpad[i][2];
    }

    for(i=0;i<16;i++){
        polyf3_sub_packed(fpad+32+i,fpad+32+i,fpad+16+i);
    }

    for(i=0;i<16;i++){
        polyf3_butterfly64_GS(fpad+i,fpad+32+i,0);
        fpad[i][1] = fpad[i][1] ^ fpad[i][0];
        fpad[i][3] = fpad[i][3] ^ fpad[i][2];
        fpad[i+32][1] = fpad[i+32][1] ^ fpad[i+32][0];
        fpad[i+32][3] = fpad[i+32][3] ^ fpad[i+32][2];
    }
    
}

void bs3_mul768_schonhage(uint32_t h[96], uint32_t f[48],uint32_t g[48]){
    uint32_t fpad[48][4],gpad[48][4];
    int i,j;
    for(i=0;i<24;i++){
		fpad[i][0]=f[2*i];
		fpad[i][1]=f[2*i+1];
		fpad[i][2]=0;
		fpad[i][3]=0;
        gpad[i][0]=g[2*i];
		gpad[i][1]=g[2*i+1];
		gpad[i][2]=0;
		gpad[i][3]=0;
	}
	for(i=24;i<48;i++){
		for(j=0;j<4;j++){
			fpad[i][j] = 0;
            gpad[i][j] = 0;
		}
	}

    fft48_64(fpad);
    fft48_64(gpad);

    for(i=0;i<48;i++){
        bs3_mul64_negc(fpad+i,fpad+i,gpad+i);
    }

    unfft48_64(fpad);
    for(i=0;i<96;i++)h[i]=0;
    for(i=0;i<47;i++){
        polyf3_add_packed(h+2*i,h+2*i,fpad+i);
    }
    uint32_t r0,r1,r2,r3,r8;
    r0 = fpad[47][0];
    r1 = fpad[47][1];
    r2 = h[94];
    r3 = h[95];

    r8 = r1 ^ r2;
    r1 = r0 ^ r3;
    r1 = r1 & r8;
    r8 = r8 ^ r3;
    r0 = r0 ^ r2;
    r0 = r0 | r8;

    h[94]=r0;
    h[95]=r1;
}