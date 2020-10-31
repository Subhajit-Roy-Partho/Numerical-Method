#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  int cells = 200,kc=(cells/2),steps=100,t0=40,t;
  double E[cells], B[cells],pulse=0,fh=12.0,a=1;


  for(int i=0; i<cells;i++)
    E[cells]= (double)0.0;
    B[cells]= (double)0.0;


  for(int n=0; n < steps; n++)
    pulse = a*exp(-0.5 * (pow((t0 - t) / fh, 2.0)));
    t+=1;
    for(int i=1; i< cells-1; i++)
      E[i]-= 0.5*(B[i-1]-B[i]);
    for(int i=1; i<cells-1;i++)
      B[i] -= 0.5*(B[i]-B[i+1]);



  return 0;
}
