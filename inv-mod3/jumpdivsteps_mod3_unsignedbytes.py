divsize = 16

with open ("jump"+str(divsize)+"divsteps_mod3.c","w") as output:
    print("#include <stdint.h>",file=output)
    print("#include \"red_mod3_int.h\"",file=output)
    print("#include \"cmsis.h\"\n",file=output)
    print("extern int jump"+str(divsize//2)+"divsteps_mod3(int minusdelta, int* M, int* f, int* g);",file=output)
    print("extern void gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(int* h, int* f, int* g);",file=output)
    print("int jump"+str(divsize)+"divsteps_mod3(int minusdelta, int* M, int* f, int* g);",file=output)
    print("void gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x2p2_mod3(int* V,int* M,int *fh, int* gh);",file=output)
    print("void gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x_2x2_mod3(int* M,int* M1,int* M2);\n",file=output)

    print("static int C"+str(divsize)+"_1["+str(divsize//4+1)+"], C"+str(divsize)+"_2["+str(divsize//4+1)+"];",file=output)
    print("int * CC"+str(divsize)+"_1 = (int *)((void *)C"+str(divsize)+"_1 + 1);",file=output)
    print("int * CC"+str(divsize)+"_2 = (int *)((void *)C"+str(divsize)+"_2 + 1);\n",file=output)

    print("void gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x2p2_mod3(int *V,int *M,int *fh,int *gh){",file=output)
    print("  int i, T, *X, *Y, *Z, *W;\n",file=output)

    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_1, M+"+str(divsize//4)+", fh); // x * u * fh",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_2, M+"+str((divsize//8)*3)+", gh); // x * v * gh\n",file=output)

    print("  for (X=V, Y=C"+str(divsize)+"_1, Z=C"+str(divsize)+"_2, W=M, i="+str(divsize//8)+"; i>0; i--) {// x(u fh+v gh)+f1",file=output)
    print("    *(X++) = add_ub3(add_ub3(*(W++),*(Y++)),*(Z++));",file=output)
    print("  }",file=output)
    print("  for (i="+str(divsize//8)+"; i>0; i--) {",file=output)
    print("    *(X++) = add_ub3(*(Y++),*(Z++));",file=output)
    print("  }\n",file=output)

    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(V+"+str(divsize//4)+", M+"+str(divsize//2)+", fh); // r * fh",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_1, M+"+str((divsize//8)*5)+", gh); // s * gh\n",file=output)

    print("  for (Y=CC"+str(divsize)+"_1, i="+str(divsize//8)+"; i>0; i--) {// r fh+s gh+g1",file=output)
    print("    T = add_ub3(add_ub3(*(W++),*(Y++)),*X); *(X++) = T;",file=output)
    print("  }",file=output)
    print("  for (i="+str(divsize//8)+"; i>0; i--) {",file=output)
    print("    T = add_ub3(*X, *(Y++)); *(X++) = T;",file=output)
    print("  }",file=output)
    print("}\n",file=output)

    print("void gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x_2x2_mod3(int *M, int *M1, int *M2){",file=output)
    print("  int i, T, *X, *Y;\n",file=output)

    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_1, M2, M1); // x * u2 * u1",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(M, M2+"+str(divsize//8)+", M1+"+str(divsize//4)+"); // v2 * r1",file=output)
    print("  for (i="+str(divsize//4)+", X=M, Y=C"+str(divsize)+"_1; i>0; i--) { // u = x * u2 * u1 + v2 * r1",file=output)
    print("    T = add_ub3(*X, *(Y++)); *(X++) = T;",file=output)
    print("  }",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_1, M2, M1+"+str(divsize//8)+"); // x * u2 * v1",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(M+"+str(divsize//4)+", M2+"+str(divsize//8)+", M1+"+str((divsize//8)*3)+"); // v2 * s1",file=output)
    print("  for (i="+str(divsize//4)+", Y=C"+str(divsize)+"_1; i>0; i--) { // v = x * u2 * v1 + v2 * s1",file=output)
    print("    T = add_ub3(*X, *(Y++)); *(X++) = T;",file=output)
    print("  }",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_1, M2+"+str(divsize//4)+", M1); // x * r2 * u1",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(M+"+str(divsize//2)+", M2+"+str((divsize//8)*3)+", M1+"+str(divsize//4)+"); // s2 * r1",file=output)
    print("  for (i="+str(divsize//4)+", Y=C"+str(divsize)+"_1; i>0; i--) { // r = x * r2 * u1 + s2 * r1",file=output)
    print("    T = add_ub3(*X, *(Y++)); *(X++) = T;",file=output)
    print("  }",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(CC"+str(divsize)+"_1, M2+"+str(divsize//4)+", M1+"+str(divsize//8)+"); // x * r2 * v1",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_mod3(M+"+str((divsize//4)*3)+", M2+"+str((divsize//8)*3)+", M1+"+str((divsize//8)*3)+"); // s2 * s1",file=output)
    print("  for (i="+str(divsize//4)+", Y=C"+str(divsize)+"_1; i>0; i--) { // s = x * r2 * v1 + s2 * s1",file=output)
    print("    T = add_ub3(*X, *(Y++)); *(X++) = T;",file=output)
    print("  }",file=output)
    print("}\n",file=output)

    print("int jump"+str(divsize)+"divsteps_mod3(int minusdelta, int* M, int* f, int* g){",file=output)
    print("  int M1["+str((divsize//4)*3)+"], M2["+str((divsize//4)*3)+"], fg["+str(divsize//2)+"];\n",file=output)

    print("  minusdelta = jump"+str(divsize//2)+"divsteps_mod3(minusdelta, M1, f, g);",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x2p2_mod3(fg, M1, f+"+str(divsize//8)+", g+"+str(divsize//8)+");",file=output)
    print("  minusdelta = jump"+str(divsize//2)+"divsteps_mod3(minusdelta, M2, fg, fg+"+str(divsize//4)+");",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x2p2_mod3(M, M2, fg+"+str(divsize//8)+", fg+"+str((divsize//8)*3)+");",file=output)
    print("  gf_polymul_"+str(divsize//2)+"x"+str(divsize//2)+"_2x2_x_2x2_mod3(M+"+str(divsize//2)+", M1+"+str(divsize//4)+", M2+"+str(divsize//4)+");",file=output)
    print("  return minusdelta;",file=output)
    print("}",file=output)