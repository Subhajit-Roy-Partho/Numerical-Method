#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char const *argv[]) {
  int cells = 200,kc=(cells/2),steps=500;
  double Ez[cells], Hx[cells],Hy[cells],pulse=0,fh=12.0,a=1.0,t0=40.0,t=0.0,Eold[2],Bold[2],factor;
  FILE *f;
  f = fopen("Outdata.dat","w");


  for(int i=0; i<cells;i++){
    Ez[i]= 0.0;
    Hx[i]= 0.0;
    Hy[i]= 0.0;
  }


  for(int n=0; n < steps; n++){
    pulse = a*exp(-0.5 * (pow((t0 - t) / fh, 2.0)));
    for(int j=0; j<cells-1;j++){
      for()
    }
  }
  fclose(f);
  f = fopen("gnuE.txt","w");
  fprintf(f, "set term gif animate delay 3 optimize\nset output 'E.gif'\nset xrange[1:%i]\nset yrange[%lf:%lf]\n",cells,-(a+0.2),(a+0.2));
  fprintf(f, "do for [i=0:%i]{\nplot 'Outdata.dat' every ::(((i-1)*200)+1)::(i*200) u 1:2 w l\n}\nexit",steps-1);
  fclose(f);
  system("gnuplot gnuE.txt");
  fflush(f);
  return 0;
}
