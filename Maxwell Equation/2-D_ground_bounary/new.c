#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
  int cells = 200,ic=(cells/2),jc=(cells/2),steps=2000;
  double Ez[cells][cells], Hx[cells][cells],Hy[cells][cells],pulse=0,fh=12.0,a=1.0,t0=40.0,t=0.0,Eold[2],Bold[2],factor;
  FILE *f;
  char file[15]="dat.dat";


  for(int i=0; i<cells;i++){
    for(int j=0; j<cells; j++){
      Ez[i][j]= 0.0;
      Hx[i][j]= 0.0;
      Hy[i][j]= 0.0;
    }
  }


  for(int n=0; n < steps; n++){

    pulse = a*exp(-0.5 * (pow((t0 - t) / fh, 2.0)));

    sprintf(file,"Data/out%d.dat",n+1);
    f = fopen(file,"w");


    for(int j=1; j<cells;j++){
      for(int i=1; i <cells; i++){
        Ez[i][j] += 0.5*(Hy[i][j] - Hy[i-1][j] - Hx[i][j]+Hx[i][j-1]);
      }
    }

    Ez[ic][jc] += pulse;

    for(int j=0; j< cells-1; j++){
      for(int i=0; i< cells-1; i++){
        Hx[i][j] += 0.5*(Ez[i][j]-Ez[i][j+1]);
      }
    }

    for(int j=0; j < cells-1;j++){
      for(int i=0; i< cells-1; i++){
        Hy[i][j] += 0.5*(Ez[i+1][j]-Ez[i][j]);
      }
    }
    for(int j=0; j<cells;j++){
      for(int i=0; i < cells; i++){
        fprintf(f, "%lf ",Ez[i][j]);
      }
      fprintf(f, "\n");
    }
    fflush(f);
    fclose(f);
  }
  f = fopen("gnuE.txt","w");
  fprintf(f, "set hidden3d\nset pm3d\nset isosample 11\nset sample 11\nset term gif large animate delay 3 optimize\nset output 'out.gif'\nset print '-'\ndo for [i=1:%d]{\nsplot 'Data/out'.i.'.dat' matrix w l\n print i\n}\nexit",steps);
  fflush(f);
  fclose(f);
  system("gnuplot gnuE.txt");
  return 0;
}