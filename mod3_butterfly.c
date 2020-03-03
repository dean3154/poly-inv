#include <stdint.h>
void polyf3_rol32_negc(uint32_t *tritIn, int k);
void polyf3_rol64_negc(uint32_t *tritIn, int k);
void polyf3_butterfly32_CT(uint32_t *tritIn1, uint32_t *tritIn2, int k);
void polyf3_butterfly32_GS(uint32_t *tritIn1, uint32_t *tritIn2, int k);
void polyf3_butterfly64_CT(uint32_t *tritIn1, uint32_t *tritIn2, int k);
void polyf3_butterfly64_GS(uint32_t *tritIn1, uint32_t *tritIn2, int k);

inline void polyf3_rol32_negc(uint32_t *tritIn, int k) {
    uint32_t r0, r1;
    uint32_t *readOp = tritIn, *writeOp = tritIn;
    r0 = *(readOp++);
    r1 = *(readOp++);

    r0 = (r0 << k) | (r0 >> (32 - k));
    r1 = r1 ^ (r0 << (32 - k));
    r1 = (r1 << k) | (r1 >> (32 - k));

    *(writeOp++) = r0;
    *(writeOp++) = r1;
}

inline void polyf3_rol64_negc(uint32_t *tritIn, int k) {
    uint32_t r0, r1, r2, r3, r4, r5, r6;
    uint32_t *readOp = tritIn, *writeOp = tritIn;
    r0 = *(readOp++);
    r1 = *(readOp++);
    r2 = *(readOp++);
    r3 = *(readOp++);
    if(k>32){
        r4 = r2;
        r2 = r0;
        r0 = r4;
        r5 = r3;
        r3 = r1;
        r1 = r5 ^ r0;
        k -= 32;
    }
    r4 = r0;
    r5 = r1;
    r6 = r2 ^ r3;

    r0 = (r0 << k) | (r2 >> (32 - k));
    r2 = (r2 << k) | (r4 >> (32 - k));
    r1 = (r1 << k) | (r6 >> (32 - k));
    r3 = (r3 << k) | (r5 >> (32 - k));

    *(writeOp++) = r0;
    *(writeOp++) = r1;
    *(writeOp++) = r2;
    *(writeOp++) = r3;
}

void polyf3_butterfly32_CT(uint32_t *tritIn1, uint32_t *tritIn2, int k) {
    uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    uint32_t *readOp1 = tritIn1, *readOp2 = tritIn2, *writeOp1 = tritIn1, *writeOp2 = tritIn2;

    polyf3_rol32_negc(tritIn2, k);

    r0 = *(readOp1++);
    r1 = *(readOp1++);
    r2 = *(readOp2++);
    r3 = *(readOp2++);

    r8 = r1 ^ r2; // r8 = a1 ^ b0
    r4 = r8 ^ r3; // r4 = (a1 ^ b0) ^ b1
    r9 = r0 ^ r2; // r9 = a0 ^ b0
    r5 = r0 ^ r3; // r5 = a0 ^ b1
    r5 = r5 & r8; // r5 = (a0 ^ b1) & (a1 ^ b0)
    r7 = r9 ^ r3; // r7 = (a0 ^ b0) ^ b1
    r7 = r7 & r8; // r7 = ((a0 ^ b0) ^ b1) & (a1 ^ b0)
    r4 = r4 | r9; // r4 = ((a1 ^ b0) ^ b1) | (a0 ^ b0)
    r6 = r1 ^ r3; // r6 = a1 ^ b1
    r6 = r6 | r9; // r6 = (a0 ^ b0) | (a1 ^ b1)

    *(writeOp1++) = r4;
    *(writeOp1++) = r5;
    *(writeOp2++) = r6;
    *(writeOp2++) = r7;
}

void polyf3_butterfly32_GS(uint32_t *tritIn1, uint32_t *tritIn2, int k) {
    uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    uint32_t *readOp1 = tritIn1, *readOp2 = tritIn2, *writeOp1 = tritIn1, *writeOp2 = tritIn2;

    r2 = *(readOp1++);
    r3 = *(readOp1++);
    r0 = *(readOp2++);
    r1 = *(readOp2++);

    r8 = r1 ^ r2; // r8 = a1 ^ b0
    r4 = r8 ^ r3; // r4 = (a1 ^ b0) ^ b1
    r9 = r0 ^ r2; // r9 = a0 ^ b0
    r5 = r0 ^ r3; // r5 = a0 ^ b1
    r5 = r5 & r8; // r5 = (a0 ^ b1) & (a1 ^ b0)
    r7 = r9 ^ r3; // r7 = (a0 ^ b0) ^ b1
    r7 = r7 & r8; // r7 = ((a0 ^ b0) ^ b1) & (a1 ^ b0)
    r4 = r4 | r9; // r4 = ((a1 ^ b0) ^ b1) | (a0 ^ b0)
    r6 = r1 ^ r3; // r6 = a1 ^ b1
    r6 = r6 | r9; // r6 = (a0 ^ b0) | (a1 ^ b1)

    *(writeOp1++) = r4;
    *(writeOp1++) = r5;
    *(writeOp2++) = r6;
    *(writeOp2++) = r7;

    polyf3_rol32_negc(tritIn2, k);
}

void polyf3_butterfly64_CT(uint32_t *tritIn1, uint32_t *tritIn2, int k) {
    uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    uint32_t *readOp1 = tritIn1, *readOp2 = tritIn2, *writeOp1 = tritIn1, *writeOp2 = tritIn2;

    polyf3_rol64_negc(tritIn2, k);

    for (int i = 0; i < 2; i++) {
        r0 = *(readOp1++);
        r1 = *(readOp1++);
        r2 = *(readOp2++);
        r3 = *(readOp2++);

        r8 = r1 ^ r2; // r8 = a1 ^ b0
        r4 = r8 ^ r3; // r4 = (a1 ^ b0) ^ b1
        r9 = r0 ^ r2; // r9 = a0 ^ b0
        r5 = r0 ^ r3; // r5 = a0 ^ b1
        r5 = r5 & r8; // r5 = (a0 ^ b1) & (a1 ^ b0)
        r7 = r9 ^ r3; // r7 = (a1 ^ b0) ^ b1
        r7 = r7 & r8; // r7 = ((a0 ^ b0) ^ b1) & (a1 ^ b0)
        r4 = r4 | r9; // r4 = ((a1 ^ b0) ^ b1) | (a0 ^ b0)
        r6 = r1 ^ r3; // r6 = a1 ^ b1
        r6 = r6 | r9; // r6 = (a0 ^ b0) | (a1 ^ b1)

        *(writeOp1++) = r4;
        *(writeOp1++) = r5;
        *(writeOp2++) = r6;
        *(writeOp2++) = r7;
    }
}

void polyf3_butterfly64_GS(uint32_t *tritIn1, uint32_t *tritIn2, int k) {
    uint32_t r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
    uint32_t *readOp1 = tritIn1, *readOp2 = tritIn2, *writeOp1 = tritIn1, *writeOp2 = tritIn2;

    for (int i = 0; i < 2; i++) {
        r2 = *(readOp1++);
        r3 = *(readOp1++);
        r0 = *(readOp2++);
        r1 = *(readOp2++);

        r8 = r1 ^ r2; // r8 = a1 ^ b0
        r4 = r8 ^ r3; // r4 = (a1 ^ b0) ^ b1
        r9 = r0 ^ r2; // r9 = a0 ^ b0
        r5 = r0 ^ r3; // r5 = a0 ^ b1
        r5 = r5 & r8; // r5 = (a0 ^ b1) & (a1 ^ b0)
        r7 = r9 ^ r3; // r7 = (a1 ^ b0) ^ b1
        r7 = r7 & r8; // r7 = ((a0 ^ b0) ^ b1) & (a1 ^ b0)
        r4 = r4 | r9; // r4 = ((a1 ^ b0) ^ b1) | (a0 ^ b0)
        r6 = r1 ^ r3; // r6 = a1 ^ b1
        r6 = r6 | r9; // r6 = (a0 ^ b0) | (a1 ^ b1)

        *(writeOp1++) = r4;
        *(writeOp1++) = r5;
        *(writeOp2++) = r6;
        *(writeOp2++) = r7;
    }

    polyf3_rol64_negc(tritIn2, k);
}
